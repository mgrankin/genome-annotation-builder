"""Restart-safe orchestration of one checksum-bound Release-C caller run.

The lifecycle directory is private runtime state.  Its final summary is
path-free and contains no dataset or subject identifier, so it may be embedded
in a portable audit report after an independent privacy check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .caller_profile import (
    load_caller_profile,
    load_prepared_caller_profile_run,
    prepare_caller_profile_run,
)
from .caller_receipt import (
    load_caller_run_receipt,
    load_prepared_caller_run_boundary,
    write_caller_run_receipt,
)
from .input_gate import DEFAULT_INPUT_GATE_POLICY, InputGatePolicy, run_input_gate
from .oci_runner import load_oci_execution_report, run_oci_caller
from .portable_model import contains_forbidden_binding, sha256_file
from .profile_catalog import resolve_shipped_profile

LIFECYCLE_SCHEMA_VERSION = 1
CHECKPOINT_NAMES = (
    "input_gate",
    "prepared_run",
    "execution_report",
    "caller_receipt",
    "lifecycle_summary",
)
PRIVATE_VALUE_KEYS = frozenset(
    {
        "dataset_id",
        "dictionary_path",
        "fai_path",
        "fasta_path",
        "host_path",
        "index_path",
        "laboratory",
        "logical_subject",
        "path",
        "sample",
        "tree_manifest_path",
    }
)


class CallerLifecycleError(RuntimeError):
    """A caller lifecycle failed without exposing its private binding."""


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise CallerLifecycleError(f"{label}_unreadable") from exc
    if not isinstance(value, dict):
        raise CallerLifecycleError(f"{label}_invalid")
    return value


def _atomic_json_no_replace(path: Path, value: object) -> None:
    if path.exists():
        raise CallerLifecycleError("lifecycle_summary_already_exists")
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
            raise CallerLifecycleError("lifecycle_summary_already_exists") from exc
        temporary.unlink()
        directory_descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _strings_for_keys(value: Any) -> set[str]:
    result: set[str] = set()
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if key in PRIVATE_VALUE_KEYS and isinstance(nested, str) and nested:
                result.add(nested)
            result.update(_strings_for_keys(nested))
    elif isinstance(value, list):
        for nested in value:
            result.update(_strings_for_keys(nested))
    return result


def _read_forbidden_file(path: Path | None) -> set[str]:
    if path is None:
        return set()
    try:
        return {
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
    except (OSError, UnicodeError) as exc:
        raise CallerLifecycleError("forbidden_binding_file_unreadable") from exc


def _private_bindings(
    manifest_path: Path,
    dataset_id: str,
    bindings_path: Path,
    resource_manifest_path: Path,
    work_directory: Path,
    forbidden_binding_file: Path | None,
) -> tuple[str, ...]:
    manifest = _load_json_object(manifest_path, "private_manifest")
    datasets = manifest.get("datasets")
    if not isinstance(datasets, list):
        raise CallerLifecycleError("private_manifest_datasets_invalid")
    matches = [
        row for row in datasets if isinstance(row, dict) and row.get("dataset_id") == dataset_id
    ]
    if len(matches) != 1:
        raise CallerLifecycleError("dataset_id_not_unique")
    binding_spec = _load_json_object(bindings_path, "caller_bindings")
    private = _strings_for_keys(manifest) | _strings_for_keys(binding_spec)
    private.update(
        {
            dataset_id,
            str(manifest_path),
            str(bindings_path),
            str(resource_manifest_path),
            str(work_directory),
        }
    )
    if forbidden_binding_file is not None:
        private.add(str(forbidden_binding_file))
    private.update(_read_forbidden_file(forbidden_binding_file))
    return tuple(sorted(value for value in private if value))


def _assert_portable(value: object, forbidden_bindings: Sequence[str]) -> None:
    text = json.dumps(value, sort_keys=True, ensure_ascii=True, allow_nan=False)
    if any(contains_forbidden_binding(text, binding) for binding in forbidden_bindings):
        raise CallerLifecycleError("portable_lifecycle_summary_contains_private_binding")


def _load_artifact_spec(path: Path) -> dict[str, Any]:
    value = _load_json_object(path, "artifact_spec")
    expected = {
        "schema_version",
        "raw_artifact_id",
        "identity_scheme",
        "components",
        "derived_payload",
    }
    if set(value) != expected or value.get("schema_version") != 1:
        raise CallerLifecycleError("artifact_spec_invalid")
    components = value.get("components")
    if not isinstance(components, dict) or any(
        not isinstance(role, str) or not isinstance(component, str)
        for role, component in components.items()
    ):
        raise CallerLifecycleError("artifact_spec_components_invalid")
    return value


def _validate_existing_gate(path: Path, expected_caller: str) -> None:
    report = _load_json_object(path, "input_gate")
    if report.get("status") != "PASS" or report.get("caller") != expected_caller:
        raise CallerLifecycleError("existing_input_gate_not_reusable")


def _profile_path(profile_id: str | None, profile_path: Path | None) -> Path:
    if (profile_id is None) == (profile_path is None):
        raise CallerLifecycleError("exactly_one_profile_source_required")
    return resolve_shipped_profile(profile_id) if profile_id is not None else profile_path


def _lifecycle_summary(
    profile: Any,
    prepared: Mapping[str, Any],
    boundary: Mapping[str, Any],
    execution: Mapping[str, Any],
    receipt: Mapping[str, Any],
) -> dict[str, Any]:
    run_boundary = boundary["run_boundary"]
    if (
        prepared["profile_id"] != profile.payload["profile_id"]
        or prepared["profile_identity_sha256"] != profile.identity_sha256
        or prepared["caller_run_identity_sha256"] != run_boundary["caller_run_identity_sha256"]
        or execution["caller_run_identity_sha256"] != run_boundary["caller_run_identity_sha256"]
        or receipt["run_boundary"] != run_boundary
    ):
        raise CallerLifecycleError("caller_lifecycle_checkpoint_identity_mismatch")
    summary: dict[str, Any] = {
        "schema_version": LIFECYCLE_SCHEMA_VERSION,
        "status": "PASS",
        "profile_id": profile.payload["profile_id"],
        "profile_identity_sha256": profile.identity_sha256,
        "caller_name": profile.payload["caller_name"],
        "caller_version": profile.payload["caller_version"],
        "gate_caller": profile.payload["gate_caller"],
        "input_gate_run_identity_sha256": run_boundary["input_gate"]["run_identity_sha256"],
        "caller_run_identity_sha256": run_boundary["caller_run_identity_sha256"],
        "execution_identity_sha256": execution["execution_identity_sha256"],
        "receipt_identity_sha256": receipt["receipt_identity_sha256"],
        "raw_artifact_id": receipt["raw_artifact"]["raw_artifact_id"],
        "raw_artifact_sha256": receipt["raw_artifact"]["raw_artifact_sha256"],
        "checkpoint_names": list(CHECKPOINT_NAMES),
        "checkpoint_count": len(CHECKPOINT_NAMES),
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    summary["lifecycle_identity_sha256"] = _canonical_sha256(summary)
    return summary


def run_caller_lifecycle(
    manifest_path: str | Path,
    dataset_id: str,
    binding_spec_path: str | Path,
    resource_manifest_path: str | Path,
    work_directory: str | Path,
    *,
    profile_id: str | None = None,
    profile_path: str | Path | None = None,
    minimum_free_bytes: int = DEFAULT_INPUT_GATE_POLICY.minimum_free_bytes,
    indexed_probe_bases: int = DEFAULT_INPUT_GATE_POLICY.indexed_probe_bases,
    forbidden_binding_file: str | Path | None = None,
) -> dict[str, Any]:
    """Run or resume one caller without replacing any retained checkpoint."""

    manifest = Path(manifest_path)
    bindings = Path(binding_spec_path)
    resources = Path(resource_manifest_path)
    work = Path(work_directory)
    forbidden_file = Path(forbidden_binding_file) if forbidden_binding_file else None
    selected_profile = _profile_path(
        profile_id, Path(profile_path) if profile_path is not None else None
    )
    private = _private_bindings(manifest, dataset_id, bindings, resources, work, forbidden_file)
    profile = load_caller_profile(selected_profile, forbidden_bindings=private)
    if work.is_symlink():
        raise CallerLifecycleError("lifecycle_directory_must_not_be_symlink")
    work.mkdir(parents=True, exist_ok=True)
    if not work.is_dir():
        raise CallerLifecycleError("lifecycle_directory_invalid")

    gate_path = work / "input-gate.json"
    plan_directory = work / "plan"
    execution_path = work / "execution-report.json"
    log_directory = work / "logs"
    receipt_path = work / "caller-receipt.json"
    summary_path = work / "lifecycle-summary.json"
    created: list[str] = []
    reused: list[str] = []

    if gate_path.exists():
        _validate_existing_gate(gate_path, profile.payload["gate_caller"])
        reused.append("input_gate")
    else:
        run_input_gate(
            manifest,
            dataset_id,
            profile.payload["gate_caller"],
            gate_path,
            policy=InputGatePolicy(
                minimum_free_bytes=minimum_free_bytes,
                indexed_probe_bases=indexed_probe_bases,
            ),
        )
        created.append("input_gate")

    if plan_directory.exists():
        prepared = load_prepared_caller_profile_run(plan_directory, forbidden_bindings=private)
        reused.append("prepared_run")
    else:
        prepared = prepare_caller_profile_run(
            selected_profile,
            bindings,
            gate_path,
            resources,
            plan_directory,
            forbidden_bindings=private,
        )
        created.append("prepared_run")
    boundary, _boundary_sha256 = load_prepared_caller_run_boundary(
        plan_directory / "prepared-boundary.json", forbidden_bindings=private
    )

    if execution_path.exists():
        execution, _execution_file_sha256 = load_oci_execution_report(
            execution_path, forbidden_bindings=private
        )
        if execution.get("status") != "PASS":
            raise CallerLifecycleError("existing_execution_report_not_reusable")
        reused.append("execution_report")
    else:
        execution = run_oci_caller(
            plan_directory / "prepared-boundary.json",
            plan_directory / "execution-spec.json",
            execution_path,
            log_directory,
            forbidden_bindings=private,
        )
        created.append("execution_report")

    if receipt_path.exists():
        receipt, _receipt_file_sha256 = load_caller_run_receipt(
            receipt_path, forbidden_bindings=private
        )
        reused.append("caller_receipt")
    else:
        artifact = _load_artifact_spec(plan_directory / "artifact-spec.json")
        receipt = write_caller_run_receipt(
            receipt_path,
            boundary["run_boundary"],
            gate_path,
            resources,
            execution_report=execution_path,
            raw_artifact_id=artifact["raw_artifact_id"],
            identity_scheme=artifact["identity_scheme"],
            component_paths=artifact["components"],
            derived_payload=artifact["derived_payload"],
            forbidden_bindings=private,
        )
        created.append("caller_receipt")

    summary = _lifecycle_summary(profile, prepared, boundary, execution, receipt)
    _assert_portable(summary, private)
    if summary_path.exists():
        existing = _load_json_object(summary_path, "lifecycle_summary")
        if (
            existing != summary
            or sha256_file(summary_path)
            != hashlib.sha256(
                (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode()
            ).hexdigest()
        ):
            raise CallerLifecycleError("existing_lifecycle_summary_changed")
        reused.append("lifecycle_summary")
    else:
        _atomic_json_no_replace(summary_path, summary)
        created.append("lifecycle_summary")

    return {
        **summary,
        "resume": {
            "created_checkpoints": created,
            "reused_checkpoints": reused,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--dataset-id", required=True)
    profile_source = parser.add_mutually_exclusive_group(required=True)
    profile_source.add_argument("--profile-id")
    profile_source.add_argument("--profile", type=Path)
    parser.add_argument("--bindings", type=Path, required=True)
    parser.add_argument("--resource-manifest", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument(
        "--minimum-free-bytes",
        type=int,
        default=DEFAULT_INPUT_GATE_POLICY.minimum_free_bytes,
    )
    parser.add_argument(
        "--indexed-probe-bases",
        type=int,
        default=DEFAULT_INPUT_GATE_POLICY.indexed_probe_bases,
    )
    parser.add_argument("--forbidden-bindings-file", type=Path)
    args = parser.parse_args()
    try:
        result = run_caller_lifecycle(
            args.manifest,
            args.dataset_id,
            args.bindings,
            args.resource_manifest,
            args.work_directory,
            profile_id=args.profile_id,
            profile_path=args.profile,
            minimum_free_bytes=args.minimum_free_bytes,
            indexed_probe_bases=args.indexed_probe_bases,
            forbidden_binding_file=args.forbidden_bindings_file,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except Exception:
        print(
            "caller lifecycle failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
