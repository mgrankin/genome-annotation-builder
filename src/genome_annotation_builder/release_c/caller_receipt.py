"""Deterministic, privacy-minimized receipts for read-backed caller runs.

The caller-run identity is prepared before execution from an exact input gate,
tool identity, resource manifest, and path-free command semantics.  A receipt is
finalized after execution by binding one logical raw artifact (a file, a named
artifact set, or a deterministic derived-evidence payload).  Typed adapter rows
can therefore be checked against both the pre-call boundary and the exact raw
output without embedding local paths or subject identifiers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .portable_model import (
    FORBIDDEN_PORTABLE_TEXT,
    SHA256,
    PortableEvidenceError,
    contains_forbidden_binding,
    sha256_file,
)

RECEIPT_SCHEMA_VERSION = 1
EXECUTION_KINDS = frozenset({"oci_image"})
ARTIFACT_IDENTITY_SCHEMES = frozenset({"file_sha256", "named_sha256_set", "derived_canonical_json"})
GATE_CALLER = re.compile(r"^[a-z0-9][a-z0-9_.-]{0,63}$")
COMPONENT_ROLE = re.compile(r"^[a-z][a-z0-9_.-]{0,95}$")
LOGICAL_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:+*/-]{0,255}$")
STEP_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]{0,127}$")
ENVIRONMENT_KEY = re.compile(r"^[A-Za-z_][A-Za-z0-9_]{0,127}$")
ABSOLUTE_LOCAL_PATH = re.compile(r"(?:^|[\s=])(?:/|[A-Za-z]:[\\/]|~[\\/])")


def canonical_json_bytes(value: Any) -> bytes:
    """Encode a finite JSON value with one deterministic representation."""

    _validate_json_tree(value)
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode()


def canonical_json_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _validate_json_tree(value: Any) -> None:
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PortableEvidenceError("caller receipt contains a non-finite JSON number")
        return
    if isinstance(value, list):
        for item in value:
            _validate_json_tree(item)
        return
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise PortableEvidenceError("caller receipt JSON object keys must be strings")
        for item in value.values():
            _validate_json_tree(item)
        return
    raise PortableEvidenceError("caller receipt contains a non-JSON value")


def _walk_strings(value: Any):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from _walk_strings(key)
            yield from _walk_strings(item)
    elif isinstance(value, (list, tuple, set, frozenset)):
        for item in value:
            yield from _walk_strings(item)


def _assert_portable(value: Any, forbidden_bindings: Sequence[str]) -> None:
    forbidden = tuple(binding for binding in forbidden_bindings if binding)
    for text in _walk_strings(value):
        if FORBIDDEN_PORTABLE_TEXT.search(text) or ABSOLUTE_LOCAL_PATH.search(text):
            raise PortableEvidenceError("caller receipt contains a local path or secret")
        if any(contains_forbidden_binding(text, binding) for binding in forbidden):
            raise PortableEvidenceError("caller receipt contains a forbidden private binding")


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise PortableEvidenceError(f"{label} is not a SHA-256")
    return value


def _regular_file(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_file():
        raise PortableEvidenceError(f"{label} must be a regular file")
    try:
        if path.resolve(strict=True) != path.absolute():
            raise PortableEvidenceError(f"{label} contains a symbolic path component")
    except OSError as exc:
        raise PortableEvidenceError(f"{label} is unreadable") from exc


def _load_json(path: Path, label: str) -> dict[str, Any]:
    _regular_file(path, label)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError(f"{label} is unreadable") from exc
    if not isinstance(value, dict):
        raise PortableEvidenceError(f"{label} must be a JSON object")
    return value


def _validate_input_gate(
    report: Mapping[str, Any], *, expected_caller: str, forbidden_bindings: Sequence[str]
) -> dict[str, str]:
    _assert_portable(report, forbidden_bindings)
    if report.get("schema_version") != 1 or report.get("status") != "PASS":
        raise PortableEvidenceError("caller receipt input gate is not PASS")
    if report.get("caller") != expected_caller:
        raise PortableEvidenceError("caller receipt input-gate caller differs")
    if (
        report.get("logical_subject_count") != 1
        or report.get("subject_identifier_persisted") is not False
        or report.get("private_paths_persisted") is not False
        or report.get("private_manifest_digest_persisted") is not False
    ):
        raise PortableEvidenceError("caller receipt input gate violates its privacy boundary")
    alignment = report.get("selected_alignment")
    reference = report.get("reference")
    if not isinstance(alignment, dict) or not isinstance(reference, dict):
        raise PortableEvidenceError("caller receipt input gate lacks artifact bindings")
    return {
        "gate_run_identity_sha256": _sha(
            report.get("run_identity_sha256"), "input-gate run identity"
        ),
        "source_alignment_sha256": _sha(alignment.get("sha256"), "input-gate alignment digest"),
        "reference_sha256": _sha(reference.get("sha256"), "input-gate reference digest"),
    }


def _validate_tool(tool: Mapping[str, Any]) -> dict[str, Any]:
    expected = {"execution_kind", "logical_name", "identity_sha256", "network_disabled"}
    if set(tool) != expected:
        raise PortableEvidenceError("caller tool identity has the wrong fields")
    execution_kind = tool.get("execution_kind")
    logical_name = tool.get("logical_name")
    if execution_kind not in EXECUTION_KINDS:
        raise PortableEvidenceError("caller tool execution kind is unsupported")
    if not isinstance(logical_name, str) or LOGICAL_NAME.fullmatch(logical_name) is None:
        raise PortableEvidenceError("caller tool logical name is invalid")
    if tool.get("network_disabled") is not True:
        raise PortableEvidenceError("read-backed caller must declare network isolation")
    return {
        "execution_kind": execution_kind,
        "logical_name": logical_name,
        "identity_sha256": _sha(tool.get("identity_sha256"), "caller tool identity"),
        "network_disabled": True,
    }


def _validate_command_semantics(value: Mapping[str, Any]) -> dict[str, Any]:
    if set(value) != {"steps", "environment", "parameters"}:
        raise PortableEvidenceError("caller command semantics has the wrong fields")
    steps = value.get("steps")
    environment = value.get("environment")
    parameters = value.get("parameters")
    if not isinstance(steps, list) or not steps:
        raise PortableEvidenceError("caller command semantics has no steps")
    normalized_steps: list[dict[str, Any]] = []
    seen_names: set[str] = set()
    for step in steps:
        if not isinstance(step, dict) or set(step) != {"name", "entrypoint", "arguments"}:
            raise PortableEvidenceError("caller command step has the wrong fields")
        name = step.get("name")
        entrypoint = step.get("entrypoint")
        arguments = step.get("arguments")
        if not isinstance(name, str) or STEP_NAME.fullmatch(name) is None:
            raise PortableEvidenceError("caller command step name is invalid")
        if name in seen_names:
            raise PortableEvidenceError("caller command step name is duplicated")
        seen_names.add(name)
        if not isinstance(entrypoint, str) or LOGICAL_NAME.fullmatch(entrypoint) is None:
            raise PortableEvidenceError("caller command entrypoint is invalid")
        if not isinstance(arguments, list) or any(
            not isinstance(argument, str) or not argument for argument in arguments
        ):
            raise PortableEvidenceError("caller command arguments are invalid")
        normalized_steps.append(
            {"name": name, "entrypoint": entrypoint, "arguments": list(arguments)}
        )
    if not isinstance(environment, dict) or any(
        not isinstance(key, str)
        or ENVIRONMENT_KEY.fullmatch(key) is None
        or not isinstance(item, str)
        for key, item in environment.items()
    ):
        raise PortableEvidenceError("caller command environment is invalid")
    if not isinstance(parameters, dict):
        raise PortableEvidenceError("caller command parameters must be an object")
    normalized = {
        "steps": normalized_steps,
        "environment": dict(sorted(environment.items())),
        "parameters": parameters,
    }
    _validate_json_tree(normalized)
    return normalized


def validate_caller_tool_identity(value: Mapping[str, Any]) -> dict[str, Any]:
    """Validate one portable immutable caller tool identity."""

    return _validate_tool(value)


def validate_caller_command_semantics(value: Mapping[str, Any]) -> dict[str, Any]:
    """Validate path-free caller semantics before a run boundary is prepared."""

    return _validate_command_semantics(value)


def _validate_boundary(boundary: Mapping[str, Any]) -> dict[str, Any]:
    expected = {
        "schema_version",
        "caller_name",
        "caller_version",
        "gate_caller",
        "input_gate",
        "source_alignment_sha256",
        "reference_sha256",
        "resource_manifest_sha256",
        "resource_manifest_bytes",
        "tool",
        "command_semantics",
        "caller_run_identity_sha256",
    }
    if set(boundary) != expected or boundary.get("schema_version") != RECEIPT_SCHEMA_VERSION:
        raise PortableEvidenceError("caller run boundary has the wrong schema")
    caller_name = boundary.get("caller_name")
    caller_version = boundary.get("caller_version")
    gate_caller = boundary.get("gate_caller")
    if not isinstance(caller_name, str) or not caller_name.strip():
        raise PortableEvidenceError("caller run boundary has an empty caller name")
    if not isinstance(caller_version, str) or not caller_version.strip():
        raise PortableEvidenceError("caller run boundary has an empty caller version")
    if not isinstance(gate_caller, str) or GATE_CALLER.fullmatch(gate_caller) is None:
        raise PortableEvidenceError("caller run boundary gate caller is invalid")
    input_gate = boundary.get("input_gate")
    if not isinstance(input_gate, dict) or set(input_gate) != {
        "report_sha256",
        "run_identity_sha256",
    }:
        raise PortableEvidenceError("caller run boundary input gate is invalid")
    resource_bytes = boundary.get("resource_manifest_bytes")
    if (
        isinstance(resource_bytes, bool)
        or not isinstance(resource_bytes, int)
        or resource_bytes < 1
    ):
        raise PortableEvidenceError("caller resource manifest byte count is invalid")
    tool = boundary.get("tool")
    command = boundary.get("command_semantics")
    if not isinstance(tool, dict) or not isinstance(command, dict):
        raise PortableEvidenceError("caller run boundary lacks tool or command semantics")
    normalized: dict[str, Any] = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "caller_name": caller_name,
        "caller_version": caller_version,
        "gate_caller": gate_caller,
        "input_gate": {
            "report_sha256": _sha(input_gate.get("report_sha256"), "input-gate report digest"),
            "run_identity_sha256": _sha(
                input_gate.get("run_identity_sha256"), "input-gate run identity"
            ),
        },
        "source_alignment_sha256": _sha(
            boundary.get("source_alignment_sha256"), "caller source alignment digest"
        ),
        "reference_sha256": _sha(boundary.get("reference_sha256"), "caller reference digest"),
        "resource_manifest_sha256": _sha(
            boundary.get("resource_manifest_sha256"), "caller resource manifest digest"
        ),
        "resource_manifest_bytes": resource_bytes,
        "tool": _validate_tool(tool),
        "command_semantics": _validate_command_semantics(command),
    }
    expected_identity = canonical_json_sha256(normalized)
    observed_identity = _sha(boundary.get("caller_run_identity_sha256"), "caller run identity")
    if observed_identity != expected_identity:
        raise PortableEvidenceError("caller run identity does not match its boundary")
    normalized["caller_run_identity_sha256"] = observed_identity
    return normalized


def prepare_caller_run_boundary(
    input_gate_report: str | Path,
    resource_manifest: str | Path,
    *,
    caller_name: str,
    caller_version: str,
    gate_caller: str,
    expected_resource_manifest_sha256: str,
    tool: Mapping[str, Any],
    command_semantics: Mapping[str, Any],
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    """Prepare the deterministic identity that an adapter records before execution."""

    gate_path = Path(input_gate_report)
    resource_path = Path(resource_manifest)
    _regular_file(gate_path, "input-gate report")
    _regular_file(resource_path, "caller resource manifest")
    if not isinstance(gate_caller, str) or GATE_CALLER.fullmatch(gate_caller) is None:
        raise PortableEvidenceError("caller gate name is invalid")
    report = _load_json(gate_path, "input-gate report")
    private = (*forbidden_bindings, str(gate_path), str(resource_path))
    gate = _validate_input_gate(report, expected_caller=gate_caller, forbidden_bindings=private)
    expected_resource = _sha(
        expected_resource_manifest_sha256, "expected caller resource manifest digest"
    )
    observed_resource = sha256_file(resource_path)
    if observed_resource != expected_resource:
        raise PortableEvidenceError("caller resource manifest digest changed before execution")
    boundary: dict[str, Any] = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "caller_name": caller_name,
        "caller_version": caller_version,
        "gate_caller": gate_caller,
        "input_gate": {
            "report_sha256": sha256_file(gate_path),
            "run_identity_sha256": gate["gate_run_identity_sha256"],
        },
        "source_alignment_sha256": gate["source_alignment_sha256"],
        "reference_sha256": gate["reference_sha256"],
        "resource_manifest_sha256": observed_resource,
        "resource_manifest_bytes": resource_path.stat().st_size,
        "tool": _validate_tool(tool),
        "command_semantics": _validate_command_semantics(command_semantics),
    }
    boundary["caller_run_identity_sha256"] = canonical_json_sha256(boundary)
    normalized = _validate_boundary(boundary)
    _assert_portable(normalized, private)
    return normalized


def _artifact_identity(
    raw_artifact_id: str,
    identity_scheme: str,
    component_paths: Mapping[str, str | Path],
    derived_payload: Mapping[str, Any] | None,
) -> tuple[dict[str, Any], tuple[str, ...]]:
    if (
        not isinstance(raw_artifact_id, str)
        or not raw_artifact_id
        or "/" in raw_artifact_id
        or "\\" in raw_artifact_id
        or "\n" in raw_artifact_id
    ):
        raise PortableEvidenceError("raw artifact ID must be a path-free logical identifier")
    if identity_scheme not in ARTIFACT_IDENTITY_SCHEMES:
        raise PortableEvidenceError("raw artifact identity scheme is unsupported")
    if not isinstance(component_paths, Mapping):
        raise PortableEvidenceError("raw artifact component mapping is invalid")
    components: list[dict[str, Any]] = []
    private_paths: list[str] = []
    for role, raw_path in sorted(component_paths.items()):
        if not isinstance(role, str) or COMPONENT_ROLE.fullmatch(role) is None:
            raise PortableEvidenceError("raw artifact component role is invalid")
        path = Path(raw_path)
        _regular_file(path, "raw artifact component")
        private_paths.append(str(path))
        components.append(
            {
                "role": role,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    if identity_scheme == "file_sha256":
        if len(components) != 1 or derived_payload is not None:
            raise PortableEvidenceError("file artifact identity requires exactly one component")
        artifact_sha256 = str(components[0]["sha256"])
        portable_derived_payload = None
    elif identity_scheme == "named_sha256_set":
        if not components or derived_payload is not None:
            raise PortableEvidenceError("named artifact-set identity requires components only")
        artifact_sha256 = canonical_json_sha256(
            {component["role"]: component["sha256"] for component in components}
        )
        portable_derived_payload = None
    else:
        if components or not isinstance(derived_payload, Mapping):
            raise PortableEvidenceError("derived artifact identity requires one JSON payload")
        portable_derived_payload = dict(derived_payload)
        artifact_sha256 = canonical_json_sha256(portable_derived_payload)
    artifact = {
        "raw_artifact_id": raw_artifact_id,
        "identity_scheme": identity_scheme,
        "raw_artifact_sha256": artifact_sha256,
        "components": components,
        "derived_payload": portable_derived_payload,
    }
    return artifact, tuple(private_paths)


def _validate_artifact(value: Mapping[str, Any]) -> dict[str, Any]:
    expected = {
        "raw_artifact_id",
        "identity_scheme",
        "raw_artifact_sha256",
        "components",
        "derived_payload",
    }
    if set(value) != expected:
        raise PortableEvidenceError("caller receipt raw artifact has the wrong fields")
    raw_artifact_id = value.get("raw_artifact_id")
    if (
        not isinstance(raw_artifact_id, str)
        or not raw_artifact_id
        or "/" in raw_artifact_id
        or "\\" in raw_artifact_id
        or "\n" in raw_artifact_id
    ):
        raise PortableEvidenceError("caller receipt raw artifact ID is invalid")
    identity_scheme = value.get("identity_scheme")
    if identity_scheme not in ARTIFACT_IDENTITY_SCHEMES:
        raise PortableEvidenceError("caller receipt artifact identity scheme is unsupported")
    raw_components = value.get("components")
    if not isinstance(raw_components, list):
        raise PortableEvidenceError("caller receipt artifact components are invalid")
    components: list[dict[str, Any]] = []
    seen_roles: set[str] = set()
    for component in raw_components:
        if not isinstance(component, dict) or set(component) != {"role", "sha256", "bytes"}:
            raise PortableEvidenceError("caller receipt artifact component is invalid")
        role = component.get("role")
        byte_count = component.get("bytes")
        if not isinstance(role, str) or COMPONENT_ROLE.fullmatch(role) is None:
            raise PortableEvidenceError("caller receipt artifact component role is invalid")
        if role in seen_roles:
            raise PortableEvidenceError("caller receipt artifact component role is duplicated")
        seen_roles.add(role)
        if isinstance(byte_count, bool) or not isinstance(byte_count, int) or byte_count < 0:
            raise PortableEvidenceError("caller receipt artifact byte count is invalid")
        components.append(
            {
                "role": role,
                "sha256": _sha(component.get("sha256"), "raw artifact component digest"),
                "bytes": byte_count,
            }
        )
    if components != sorted(components, key=lambda component: component["role"]):
        raise PortableEvidenceError("caller receipt artifact components are not canonical")
    derived_payload = value.get("derived_payload")
    if identity_scheme == "file_sha256":
        if len(components) != 1 or derived_payload is not None:
            raise PortableEvidenceError("caller receipt file artifact is invalid")
        expected_sha256 = components[0]["sha256"]
    elif identity_scheme == "named_sha256_set":
        if not components or derived_payload is not None:
            raise PortableEvidenceError("caller receipt named artifact set is invalid")
        expected_sha256 = canonical_json_sha256(
            {component["role"]: component["sha256"] for component in components}
        )
    else:
        if components or not isinstance(derived_payload, dict):
            raise PortableEvidenceError("caller receipt derived artifact is invalid")
        expected_sha256 = canonical_json_sha256(derived_payload)
    observed_sha256 = _sha(value.get("raw_artifact_sha256"), "raw artifact digest")
    if observed_sha256 != expected_sha256:
        raise PortableEvidenceError("caller receipt raw artifact identity changed")
    return {
        "raw_artifact_id": raw_artifact_id,
        "identity_scheme": identity_scheme,
        "raw_artifact_sha256": observed_sha256,
        "components": components,
        "derived_payload": derived_payload,
    }


def _validate_execution_binding(
    value: Mapping[str, Any],
    boundary: Mapping[str, Any],
    *,
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    from .oci_runner import validate_oci_execution_report_payload

    if set(value) != {"execution_report_file_sha256", "report"}:
        raise PortableEvidenceError("caller receipt execution binding has the wrong fields")
    report = value.get("report")
    if not isinstance(report, dict):
        raise PortableEvidenceError("caller receipt execution report is invalid")
    normalized = validate_oci_execution_report_payload(
        report, forbidden_bindings=forbidden_bindings
    )
    if normalized != report or report.get("status") != "PASS":
        raise PortableEvidenceError("caller receipt execution report is not canonical PASS")
    if (
        report["caller_run_identity_sha256"] != boundary["caller_run_identity_sha256"]
        or report["source_alignment_sha256"] != boundary["source_alignment_sha256"]
        or report["reference_sha256"] != boundary["reference_sha256"]
        or report["resource_manifest_sha256"] != boundary["resource_manifest_sha256"]
        or report["tool"] != boundary["tool"]
    ):
        raise PortableEvidenceError("caller receipt execution report differs from its run boundary")
    expected_steps = [
        (step["name"], canonical_json_sha256(step))
        for step in boundary["command_semantics"]["steps"]
    ]
    observed_steps = [(step["name"], step["semantic_step_sha256"]) for step in report["steps"]]
    if report["declared_steps"] != len(expected_steps) or observed_steps != expected_steps:
        raise PortableEvidenceError("caller execution steps differ from command semantics")
    observed_file_sha256 = _sha(
        value.get("execution_report_file_sha256"), "caller execution report file digest"
    )
    expected_file_sha256 = hashlib.sha256(
        (json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    ).hexdigest()
    if observed_file_sha256 != expected_file_sha256:
        raise PortableEvidenceError("caller execution report file identity changed")
    return {
        "execution_report_file_sha256": observed_file_sha256,
        "report": normalized,
    }


def _load_execution_binding(
    path: Path,
    boundary: Mapping[str, Any],
    *,
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    from .oci_runner import load_oci_execution_report

    report, file_sha256 = load_oci_execution_report(path, forbidden_bindings=forbidden_bindings)
    return _validate_execution_binding(
        {"execution_report_file_sha256": file_sha256, "report": report},
        boundary,
        forbidden_bindings=forbidden_bindings,
    )


def validate_caller_run_receipt_payload(
    value: Mapping[str, Any], *, forbidden_bindings: Sequence[str] = ()
) -> dict[str, Any]:
    expected = {
        "schema_version",
        "status",
        "run_boundary",
        "execution",
        "raw_artifact",
        "logical_subject_count",
        "subject_identifier_persisted",
        "private_paths_persisted",
        "private_manifest_digest_persisted",
        "receipt_identity_sha256",
    }
    if set(value) != expected or value.get("schema_version") != RECEIPT_SCHEMA_VERSION:
        raise PortableEvidenceError("caller receipt has the wrong schema")
    if value.get("status") != "PASS":
        raise PortableEvidenceError("caller receipt is not PASS")
    if (
        value.get("logical_subject_count") != 1
        or value.get("subject_identifier_persisted") is not False
        or value.get("private_paths_persisted") is not False
        or value.get("private_manifest_digest_persisted") is not False
    ):
        raise PortableEvidenceError("caller receipt violates its privacy/subject boundary")
    boundary = value.get("run_boundary")
    execution = value.get("execution")
    artifact = value.get("raw_artifact")
    if (
        not isinstance(boundary, dict)
        or not isinstance(execution, dict)
        or not isinstance(artifact, dict)
    ):
        raise PortableEvidenceError(
            "caller receipt lacks its run boundary, execution, or raw artifact"
        )
    normalized_boundary = _validate_boundary(boundary)
    normalized: dict[str, Any] = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": "PASS",
        "run_boundary": normalized_boundary,
        "execution": _validate_execution_binding(
            execution,
            normalized_boundary,
            forbidden_bindings=forbidden_bindings,
        ),
        "raw_artifact": _validate_artifact(artifact),
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    expected_identity = canonical_json_sha256(normalized)
    observed_identity = _sha(value.get("receipt_identity_sha256"), "caller receipt identity")
    if observed_identity != expected_identity:
        raise PortableEvidenceError("caller receipt identity changed")
    normalized["receipt_identity_sha256"] = observed_identity
    _assert_portable(normalized, forbidden_bindings)
    return normalized


def load_caller_run_receipt(
    path: str | Path, *, forbidden_bindings: Sequence[str] = ()
) -> tuple[dict[str, Any], str]:
    receipt_path = Path(path)
    value = _load_json(receipt_path, "caller-run receipt")
    normalized = validate_caller_run_receipt_payload(
        value, forbidden_bindings=(*forbidden_bindings, str(receipt_path))
    )
    return normalized, sha256_file(receipt_path)


def _atomic_json_no_replace(path: Path, value: Mapping[str, Any], *, artifact_label: str) -> None:
    if path.exists():
        raise PortableEvidenceError(f"refusing to replace an existing {artifact_label}")
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    descriptor, raw_name = tempfile.mkstemp(prefix=f".{path.name}.partial.", dir=path.parent)
    temporary = Path(raw_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as exc:
            raise PortableEvidenceError(
                f"refusing to replace an existing {artifact_label}"
            ) from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def validate_prepared_caller_run_boundary_payload(
    value: Mapping[str, Any], *, forbidden_bindings: Sequence[str] = ()
) -> dict[str, Any]:
    expected = {
        "schema_version",
        "status",
        "run_boundary",
        "logical_subject_count",
        "subject_identifier_persisted",
        "private_paths_persisted",
        "private_manifest_digest_persisted",
        "boundary_checkpoint_identity_sha256",
    }
    if set(value) != expected or value.get("schema_version") != RECEIPT_SCHEMA_VERSION:
        raise PortableEvidenceError("prepared caller boundary has the wrong schema")
    if value.get("status") != "PREPARED":
        raise PortableEvidenceError("caller boundary checkpoint is not PREPARED")
    if (
        value.get("logical_subject_count") != 1
        or value.get("subject_identifier_persisted") is not False
        or value.get("private_paths_persisted") is not False
        or value.get("private_manifest_digest_persisted") is not False
    ):
        raise PortableEvidenceError("prepared caller boundary violates its privacy boundary")
    boundary = value.get("run_boundary")
    if not isinstance(boundary, dict):
        raise PortableEvidenceError("prepared caller boundary lacks its run boundary")
    normalized: dict[str, Any] = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": "PREPARED",
        "run_boundary": _validate_boundary(boundary),
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    expected_identity = canonical_json_sha256(normalized)
    observed_identity = _sha(
        value.get("boundary_checkpoint_identity_sha256"),
        "prepared caller boundary checkpoint identity",
    )
    if observed_identity != expected_identity:
        raise PortableEvidenceError("prepared caller boundary checkpoint identity changed")
    normalized["boundary_checkpoint_identity_sha256"] = observed_identity
    _assert_portable(normalized, forbidden_bindings)
    return normalized


def write_prepared_caller_run_boundary(
    output_path: str | Path,
    input_gate_report: str | Path,
    resource_manifest: str | Path,
    *,
    caller_name: str,
    caller_version: str,
    gate_caller: str,
    expected_resource_manifest_sha256: str,
    tool: Mapping[str, Any],
    command_semantics: Mapping[str, Any],
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    """Atomically persist a restart-safe, path-free pre-execution checkpoint."""

    output = Path(output_path)
    if output.exists():
        raise PortableEvidenceError("refusing to replace an existing caller boundary")
    boundary = prepare_caller_run_boundary(
        input_gate_report,
        resource_manifest,
        caller_name=caller_name,
        caller_version=caller_version,
        gate_caller=gate_caller,
        expected_resource_manifest_sha256=expected_resource_manifest_sha256,
        tool=tool,
        command_semantics=command_semantics,
        forbidden_bindings=forbidden_bindings,
    )
    checkpoint: dict[str, Any] = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": "PREPARED",
        "run_boundary": boundary,
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    checkpoint["boundary_checkpoint_identity_sha256"] = canonical_json_sha256(checkpoint)
    normalized = validate_prepared_caller_run_boundary_payload(
        checkpoint, forbidden_bindings=(*forbidden_bindings, str(output))
    )
    _atomic_json_no_replace(output, normalized, artifact_label="caller boundary")
    return normalized


def load_prepared_caller_run_boundary(
    path: str | Path, *, forbidden_bindings: Sequence[str] = ()
) -> tuple[dict[str, Any], str]:
    boundary_path = Path(path)
    value = _load_json(boundary_path, "prepared caller boundary")
    normalized = validate_prepared_caller_run_boundary_payload(
        value, forbidden_bindings=(*forbidden_bindings, str(boundary_path))
    )
    return normalized, sha256_file(boundary_path)


def write_caller_run_receipt(
    output_path: str | Path,
    boundary: Mapping[str, Any],
    input_gate_report: str | Path,
    resource_manifest: str | Path,
    *,
    execution_report: str | Path,
    raw_artifact_id: str,
    identity_scheme: str,
    component_paths: Mapping[str, str | Path] | None = None,
    derived_payload: Mapping[str, Any] | None = None,
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    """Re-verify the pre-call boundary and atomically bind its raw output."""

    output = Path(output_path)
    if output.exists():
        raise PortableEvidenceError("refusing to replace an existing caller receipt")
    normalized_boundary = _validate_boundary(boundary)
    gate_path = Path(input_gate_report)
    resource_path = Path(resource_manifest)
    execution_path = Path(execution_report)
    _regular_file(gate_path, "input-gate report")
    _regular_file(resource_path, "caller resource manifest")
    if sha256_file(gate_path) != normalized_boundary["input_gate"]["report_sha256"]:
        raise PortableEvidenceError("input-gate report changed during caller execution")
    if (
        sha256_file(resource_path) != normalized_boundary["resource_manifest_sha256"]
        or resource_path.stat().st_size != normalized_boundary["resource_manifest_bytes"]
    ):
        raise PortableEvidenceError("resource manifest changed during caller execution")
    report = _load_json(gate_path, "input-gate report")
    private = (
        *forbidden_bindings,
        str(output),
        str(gate_path),
        str(resource_path),
        str(execution_path),
    )
    gate = _validate_input_gate(
        report,
        expected_caller=normalized_boundary["gate_caller"],
        forbidden_bindings=private,
    )
    if (
        gate["gate_run_identity_sha256"] != normalized_boundary["input_gate"]["run_identity_sha256"]
        or gate["source_alignment_sha256"] != normalized_boundary["source_alignment_sha256"]
        or gate["reference_sha256"] != normalized_boundary["reference_sha256"]
    ):
        raise PortableEvidenceError("input-gate semantics changed during caller execution")
    execution = _load_execution_binding(
        execution_path,
        normalized_boundary,
        forbidden_bindings=private,
    )
    artifact, component_private_paths = _artifact_identity(
        raw_artifact_id,
        identity_scheme,
        component_paths or {},
        derived_payload,
    )
    private = (*private, *component_private_paths)
    receipt: dict[str, Any] = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": "PASS",
        "run_boundary": normalized_boundary,
        "execution": execution,
        "raw_artifact": artifact,
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    receipt["receipt_identity_sha256"] = canonical_json_sha256(receipt)
    normalized = validate_caller_run_receipt_payload(receipt, forbidden_bindings=private)
    _atomic_json_no_replace(output, normalized, artifact_label="caller receipt")
    return normalized


def _read_forbidden_bindings(path: Path | None) -> tuple[str, ...]:
    if path is None:
        return ()
    _regular_file(path, "forbidden-binding file")
    try:
        return tuple(
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("forbidden-binding file is unreadable") from exc


def _add_forbidden_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--forbidden-bindings-file", type=Path)


def _artifact_spec(path: Path) -> dict[str, Any]:
    value = _load_json(path, "raw artifact specification")
    if (
        set(value)
        != {
            "schema_version",
            "raw_artifact_id",
            "identity_scheme",
            "components",
            "derived_payload",
        }
        or value.get("schema_version") != RECEIPT_SCHEMA_VERSION
    ):
        raise PortableEvidenceError("raw artifact specification has the wrong schema")
    components = value.get("components")
    derived_payload = value.get("derived_payload")
    if not isinstance(components, dict) or any(
        not isinstance(role, str) or not isinstance(component, str)
        for role, component in components.items()
    ):
        raise PortableEvidenceError("raw artifact specification components are invalid")
    if derived_payload is not None and not isinstance(derived_payload, dict):
        raise PortableEvidenceError("raw artifact derived payload is invalid")
    return value


def _print_summary(value: Mapping[str, Any]) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def _cli_forbidden_context(
    bindings: Sequence[str], binding_file: Path | None, *paths: Path
) -> tuple[str, ...]:
    result = (*bindings, *(str(path) for path in paths))
    if binding_file is not None:
        result = (*result, str(binding_file))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--input-gate", type=Path, required=True)
    prepare_parser.add_argument("--resource-manifest", type=Path, required=True)
    prepare_parser.add_argument("--expected-resource-manifest-sha256", required=True)
    prepare_parser.add_argument("--caller-name", required=True)
    prepare_parser.add_argument("--caller-version", required=True)
    prepare_parser.add_argument("--gate-caller", required=True)
    prepare_parser.add_argument("--tool-json", type=Path, required=True)
    prepare_parser.add_argument("--command-json", type=Path, required=True)
    prepare_parser.add_argument("--output", type=Path, required=True)
    _add_forbidden_option(prepare_parser)

    finalize_parser = subparsers.add_parser("finalize")
    finalize_parser.add_argument("--boundary", type=Path, required=True)
    finalize_parser.add_argument("--input-gate", type=Path, required=True)
    finalize_parser.add_argument("--resource-manifest", type=Path, required=True)
    finalize_parser.add_argument("--execution-report", type=Path, required=True)
    finalize_parser.add_argument("--artifact-spec", type=Path, required=True)
    finalize_parser.add_argument("--output", type=Path, required=True)
    _add_forbidden_option(finalize_parser)

    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--receipt", type=Path, required=True)
    _add_forbidden_option(validate_parser)

    boundary_parser = subparsers.add_parser("validate-boundary")
    boundary_parser.add_argument("--boundary", type=Path, required=True)
    _add_forbidden_option(boundary_parser)

    args = parser.parse_args()
    try:
        forbidden = _read_forbidden_bindings(args.forbidden_bindings_file)
        if args.command == "prepare":
            tool = _load_json(args.tool_json, "caller tool specification")
            command_semantics = _load_json(
                args.command_json, "caller command-semantics specification"
            )
            private = _cli_forbidden_context(
                forbidden,
                args.forbidden_bindings_file,
                args.tool_json,
                args.command_json,
            )
            checkpoint = write_prepared_caller_run_boundary(
                args.output,
                args.input_gate,
                args.resource_manifest,
                caller_name=args.caller_name,
                caller_version=args.caller_version,
                gate_caller=args.gate_caller,
                expected_resource_manifest_sha256=args.expected_resource_manifest_sha256,
                tool=tool,
                command_semantics=command_semantics,
                forbidden_bindings=private,
            )
            _print_summary(
                {
                    "schema_version": RECEIPT_SCHEMA_VERSION,
                    "status": "PREPARED",
                    "caller_run_identity_sha256": checkpoint["run_boundary"][
                        "caller_run_identity_sha256"
                    ],
                    "boundary_checkpoint_identity_sha256": checkpoint[
                        "boundary_checkpoint_identity_sha256"
                    ],
                    "boundary_file_sha256": sha256_file(args.output),
                }
            )
            return 0
        if args.command == "finalize":
            checkpoint, _ = load_prepared_caller_run_boundary(
                args.boundary, forbidden_bindings=forbidden
            )
            specification = _artifact_spec(args.artifact_spec)
            private = _cli_forbidden_context(
                forbidden,
                args.forbidden_bindings_file,
                args.boundary,
                args.execution_report,
                args.artifact_spec,
            )
            receipt = write_caller_run_receipt(
                args.output,
                checkpoint["run_boundary"],
                args.input_gate,
                args.resource_manifest,
                execution_report=args.execution_report,
                raw_artifact_id=specification["raw_artifact_id"],
                identity_scheme=specification["identity_scheme"],
                component_paths=specification["components"],
                derived_payload=specification["derived_payload"],
                forbidden_bindings=private,
            )
            _print_summary(
                {
                    "schema_version": RECEIPT_SCHEMA_VERSION,
                    "status": "PASS",
                    "caller_run_identity_sha256": receipt["run_boundary"][
                        "caller_run_identity_sha256"
                    ],
                    "receipt_identity_sha256": receipt["receipt_identity_sha256"],
                    "receipt_file_sha256": sha256_file(args.output),
                }
            )
            return 0
        if args.command == "validate":
            receipt, file_sha256 = load_caller_run_receipt(
                args.receipt, forbidden_bindings=forbidden
            )
            _print_summary(
                {
                    "schema_version": RECEIPT_SCHEMA_VERSION,
                    "status": "PASS",
                    "caller_run_identity_sha256": receipt["run_boundary"][
                        "caller_run_identity_sha256"
                    ],
                    "receipt_identity_sha256": receipt["receipt_identity_sha256"],
                    "receipt_file_sha256": file_sha256,
                }
            )
            return 0
        checkpoint, file_sha256 = load_prepared_caller_run_boundary(
            args.boundary, forbidden_bindings=forbidden
        )
        _print_summary(
            {
                "schema_version": RECEIPT_SCHEMA_VERSION,
                "status": "PREPARED",
                "caller_run_identity_sha256": checkpoint["run_boundary"][
                    "caller_run_identity_sha256"
                ],
                "boundary_checkpoint_identity_sha256": checkpoint[
                    "boundary_checkpoint_identity_sha256"
                ],
                "boundary_file_sha256": file_sha256,
            }
        )
        return 0
    except PortableEvidenceError:
        print(
            "caller receipt operation failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
