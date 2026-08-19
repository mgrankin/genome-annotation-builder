"""Fail-closed OCI execution for prepared read-backed caller boundaries.

This runner turns the portable ``network_disabled`` tool claim into a command
constructed by policy: a locally present digest-pinned image, ``--pull never``,
no network, read-only root, no capabilities, no-new-privileges, a non-root user,
bounded resources, explicit bind mounts, and no shell interpolation.  Private
host paths remain in an ignored execution spec; the retained execution report
contains only digests, counts, generic logical roles, and log hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path, PurePosixPath
from typing import Any

from .caller_receipt import (
    canonical_json_sha256,
    load_prepared_caller_run_boundary,
)
from .portable_model import (
    FORBIDDEN_PORTABLE_TEXT,
    SHA256,
    PortableEvidenceError,
    contains_forbidden_binding,
    sha256_file,
)
from .resource_tree import validate_resource_tree_manifest

EXECUTION_SCHEMA_VERSION = 1
RUNTIMES = frozenset({"docker", "podman"})
MOUNT_ROLE = re.compile(r"^[a-z][a-z0-9_.-]{0,95}$")
STEP_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]{0,127}$")
PLACEHOLDER = re.compile(r"\{[A-Za-z][A-Za-z0-9_.:-]{0,95}\}")
OCI_IMAGE = re.compile(r"^[^\s,@]+@sha256:([0-9a-f]{64})$")
NON_ROOT_USER = re.compile(r"^[1-9][0-9]*:[1-9][0-9]*$")
ABSOLUTE_LOCAL_PATH = re.compile(r"(?:^|[\s=])(?:/|[A-Za-z]:[\\/]|~[\\/])")
BANNED_CONTAINER_PREFIXES = (
    PurePosixPath("/proc"),
    PurePosixPath("/sys"),
    PurePosixPath("/dev"),
    PurePosixPath("/var/run"),
)


class OciExecutionFailure(PortableEvidenceError):
    """A sanitized OCI failure, optionally retained in an execution report."""

    def __init__(self, code: str, *, report_written: bool = False):
        super().__init__(code)
        self.code = code
        self.report_written = report_written


def _fail(code: str) -> None:
    raise OciExecutionFailure(code)


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise OciExecutionFailure(label)
    return value


def _integer(value: Any, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise OciExecutionFailure(label)
    return value


def _signed_integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise OciExecutionFailure(label)
    return value


def _finite_number(value: Any, label: str, *, minimum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise OciExecutionFailure(label)
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < minimum:
        raise OciExecutionFailure(label)
    return parsed


def _regular_file(path: Path, code: str) -> None:
    if path.is_symlink() or not path.is_file():
        _fail(code)
    try:
        if path.resolve(strict=True) != path.absolute():
            _fail(code)
    except OSError as exc:
        raise OciExecutionFailure(code) from exc


def _regular_directory(path: Path, code: str) -> None:
    if path.is_symlink() or not path.is_dir():
        _fail(code)
    try:
        if path.resolve(strict=True) != path.absolute():
            _fail(code)
    except OSError as exc:
        raise OciExecutionFailure(code) from exc


def _load_json(path: Path, code: str) -> dict[str, Any]:
    _regular_file(path, code)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise OciExecutionFailure(code) from exc
    if not isinstance(value, dict):
        _fail(code)
    return value


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
            _fail("portable_execution_report_path_or_secret")
        if any(contains_forbidden_binding(text, binding) for binding in forbidden):
            _fail("portable_execution_report_private_binding")


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    if path.exists():
        _fail("execution_report_already_exists")
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
            raise OciExecutionFailure("execution_report_already_exists") from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _container_path(value: Any, code: str) -> PurePosixPath:
    if not isinstance(value, str) or not value.startswith("/") or "," in value or ":" in value:
        _fail(code)
    path = PurePosixPath(value)
    if str(path) != value or value == "/" or ".." in path.parts:
        _fail(code)
    if any(path == prefix or prefix in path.parents for prefix in BANNED_CONTAINER_PREFIXES):
        _fail(code)
    return path


def _host_path(value: Any, code: str) -> Path:
    if not isinstance(value, str) or not value.startswith("/") or "," in value:
        _fail(code)
    return Path(value)


def _stat_signature(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def _path_is_mounted(path: PurePosixPath, mounts: list[dict[str, Any]]) -> bool:
    for mount in mounts:
        target = mount["container_path"]
        if mount["kind"] == "file" and path == target:
            return True
        if mount["kind"] == "directory" and (path == target or target in path.parents):
            return True
    return False


def _render(value: str, bindings: Mapping[str, str], used: set[str]) -> str:
    rendered = value
    for placeholder in sorted(bindings, key=lambda item: (-len(item), item)):
        if placeholder in rendered:
            rendered = rendered.replace(placeholder, bindings[placeholder])
            used.add(placeholder)
    if PLACEHOLDER.search(rendered):
        _fail("command_placeholder_unbound")
    return rendered


def _validate_execution_spec(
    specification: Mapping[str, Any], boundary: Mapping[str, Any]
) -> dict[str, Any]:
    expected_fields = {
        "schema_version",
        "runtime",
        "image",
        "user",
        "mounts",
        "placeholder_bindings",
        "tmpfs",
        "limits",
        "minimum_free_bytes",
    }
    if set(specification) != expected_fields or specification.get("schema_version") != 1:
        _fail("execution_spec_schema_invalid")
    runtime = specification.get("runtime")
    image = specification.get("image")
    user = specification.get("user")
    if runtime not in RUNTIMES:
        _fail("execution_runtime_unsupported")
    if not isinstance(image, str) or (image_match := OCI_IMAGE.fullmatch(image)) is None:
        _fail("execution_image_not_digest_pinned")
    tool = boundary["tool"]
    if tool["execution_kind"] != "oci_image" or tool["network_disabled"] is not True:
        _fail("caller_boundary_is_not_network_isolated_oci")
    if image_match.group(1) != tool["identity_sha256"]:
        _fail("execution_image_differs_from_caller_boundary")
    if not isinstance(user, str) or NON_ROOT_USER.fullmatch(user) is None:
        _fail("execution_user_must_be_non_root_numeric")
    minimum_free_bytes = _integer(
        specification.get("minimum_free_bytes"), "execution_free_space_policy_invalid"
    )

    raw_mounts = specification.get("mounts")
    if not isinstance(raw_mounts, list) or not raw_mounts:
        _fail("execution_mounts_invalid")
    mounts: list[dict[str, Any]] = []
    evidence: list[dict[str, Any]] = []
    signatures: dict[Path, tuple[int, int, int, int, int]] = {}
    seen_roles: set[str] = set()
    seen_targets: set[PurePosixPath] = set()
    for raw in raw_mounts:
        fields = {
            "role",
            "kind",
            "host_path",
            "container_path",
            "mode",
            "bytes",
            "sha256",
            "tree_manifest_path",
            "tree_manifest_bytes",
            "tree_manifest_sha256",
        }
        if not isinstance(raw, dict) or set(raw) != fields:
            _fail("execution_mount_schema_invalid")
        role = raw.get("role")
        kind = raw.get("kind")
        mode = raw.get("mode")
        if not isinstance(role, str) or MOUNT_ROLE.fullmatch(role) is None:
            _fail("execution_mount_role_invalid")
        if role in seen_roles:
            _fail("execution_mount_role_duplicated")
        seen_roles.add(role)
        if kind not in {"file", "directory"} or mode not in {"ro", "rw"}:
            _fail("execution_mount_kind_or_mode_invalid")
        host = _host_path(raw.get("host_path"), "execution_mount_host_path_invalid")
        target = _container_path(
            raw.get("container_path"), "execution_mount_container_path_invalid"
        )
        if target in seen_targets or any(
            target in existing.parents or existing in target.parents for existing in seen_targets
        ):
            _fail("execution_mount_target_duplicated")
        seen_targets.add(target)
        mount: dict[str, Any] = {
            "role": role,
            "kind": kind,
            "mode": mode,
            "host_path": host,
            "container_path": target,
        }
        if kind == "file":
            if mode != "ro" or any(
                raw.get(key) is not None
                for key in (
                    "tree_manifest_path",
                    "tree_manifest_bytes",
                    "tree_manifest_sha256",
                )
            ):
                _fail("execution_file_mount_must_be_read_only_and_digest_bound")
            _regular_file(host, "execution_file_mount_missing_or_symlinked")
            byte_count = _integer(raw.get("bytes"), "execution_file_mount_binding_invalid")
            digest = _sha(raw.get("sha256"), "execution_file_mount_binding_invalid")
            if host.stat().st_size != byte_count or sha256_file(host) != digest:
                _fail("execution_file_mount_digest_or_size_changed")
            signatures[host] = _stat_signature(host)
            evidence.append(
                {
                    "role": role,
                    "kind": "file",
                    "mode": "ro",
                    "sha256": digest,
                    "bytes": byte_count,
                    "tree_manifest_sha256": None,
                    "tree_identity_sha256": None,
                    "tree_files": None,
                    "tree_directories": None,
                }
            )
        elif mode == "ro":
            if raw.get("bytes") is not None or raw.get("sha256") is not None:
                _fail("execution_directory_mount_binding_invalid")
            manifest_path = _host_path(
                raw.get("tree_manifest_path"), "directory_manifest_path_invalid"
            )
            manifest_bytes = _integer(
                raw.get("tree_manifest_bytes"), "directory_manifest_binding_invalid", minimum=1
            )
            manifest_sha256 = _sha(
                raw.get("tree_manifest_sha256"), "directory_manifest_binding_invalid"
            )
            try:
                validated_tree = validate_resource_tree_manifest(
                    host,
                    manifest_path,
                    expected_manifest_bytes=manifest_bytes,
                    expected_manifest_sha256=manifest_sha256,
                )
            except PortableEvidenceError as exc:
                raise OciExecutionFailure("read_only_directory_manifest_validation_failed") from exc
            signatures.update(validated_tree.signatures)
            mount["tree_manifest_path"] = manifest_path
            evidence.append(
                {
                    "role": role,
                    "kind": "directory",
                    "mode": "ro",
                    "sha256": None,
                    "bytes": validated_tree.logical_bytes,
                    "tree_manifest_sha256": validated_tree.manifest_sha256,
                    "tree_identity_sha256": validated_tree.tree_identity_sha256,
                    "tree_files": validated_tree.files,
                    "tree_directories": validated_tree.directories,
                }
            )
        else:
            if role != "workspace" or any(
                raw.get(key) is not None
                for key in (
                    "bytes",
                    "sha256",
                    "tree_manifest_path",
                    "tree_manifest_bytes",
                    "tree_manifest_sha256",
                )
            ):
                _fail("execution_writable_mount_must_be_unique_workspace")
            _regular_directory(host, "execution_workspace_missing_or_symlinked")
            evidence.append(
                {
                    "role": role,
                    "kind": "directory",
                    "mode": "rw",
                    "sha256": None,
                    "bytes": None,
                    "tree_manifest_sha256": None,
                    "tree_identity_sha256": None,
                    "tree_files": None,
                    "tree_directories": None,
                }
            )
        mounts.append(mount)
    if (
        seen_roles.intersection({"selected_alignment", "reference_fasta", "resource_manifest"})
        != {
            "selected_alignment",
            "reference_fasta",
            "resource_manifest",
        }
        or "workspace" not in seen_roles
    ):
        _fail("execution_required_mount_roles_missing")
    if sum(mount["mode"] == "rw" for mount in mounts) != 1:
        _fail("execution_writable_mount_must_be_unique_workspace")
    evidence_by_role = {item["role"]: item for item in evidence}
    if evidence_by_role["selected_alignment"]["sha256"] != boundary["source_alignment_sha256"]:
        _fail("mounted_alignment_differs_from_caller_boundary")
    if evidence_by_role["reference_fasta"]["sha256"] != boundary["reference_sha256"]:
        _fail("mounted_reference_differs_from_caller_boundary")
    if (
        evidence_by_role["resource_manifest"]["sha256"] != boundary["resource_manifest_sha256"]
        or evidence_by_role["resource_manifest"]["bytes"] != boundary["resource_manifest_bytes"]
    ):
        _fail("mounted_resource_manifest_differs_from_caller_boundary")

    raw_bindings = specification.get("placeholder_bindings")
    if not isinstance(raw_bindings, dict) or not raw_bindings:
        _fail("execution_placeholder_bindings_invalid")
    bindings: dict[str, str] = {}
    for placeholder, raw_target in raw_bindings.items():
        if not isinstance(placeholder, str) or PLACEHOLDER.fullmatch(placeholder) is None:
            _fail("execution_placeholder_name_invalid")
        target = _container_path(raw_target, "execution_placeholder_target_invalid")
        if not _path_is_mounted(target, mounts):
            _fail("execution_placeholder_target_is_not_mounted")
        bindings[placeholder] = str(target)

    raw_tmpfs = specification.get("tmpfs")
    if not isinstance(raw_tmpfs, list) or not raw_tmpfs:
        _fail("execution_tmpfs_policy_invalid")
    tmpfs: list[dict[str, Any]] = []
    seen_tmpfs: set[PurePosixPath] = set()
    for item in raw_tmpfs:
        if not isinstance(item, dict) or set(item) != {
            "container_path",
            "bytes",
            "executable",
        }:
            _fail("execution_tmpfs_policy_invalid")
        target = _container_path(item.get("container_path"), "execution_tmpfs_path_invalid")
        if (
            target in seen_tmpfs
            or any(
                target in existing.parents or existing in target.parents for existing in seen_tmpfs
            )
            or target in seen_targets
            or any(
                target in existing.parents or existing in target.parents
                for existing in seen_targets
            )
        ):
            _fail("execution_tmpfs_path_conflict")
        seen_tmpfs.add(target)
        byte_count = _integer(
            item.get("bytes"), "execution_tmpfs_policy_invalid", minimum=64 * 1024**2
        )
        if item.get("executable") is not False:
            _fail("execution_tmpfs_must_be_noexec")
        tmpfs.append({"container_path": target, "bytes": byte_count, "executable": False})
    if PurePosixPath("/tmp") not in seen_tmpfs:
        _fail("execution_tmpfs_requires_private_tmp")

    raw_limits = specification.get("limits")
    if not isinstance(raw_limits, dict) or set(raw_limits) != {
        "pids",
        "memory_bytes",
        "cpus",
        "shm_bytes",
    }:
        _fail("execution_limits_invalid")
    limits = {
        "pids": _integer(raw_limits.get("pids"), "execution_limits_invalid", minimum=1),
        "memory_bytes": _integer(
            raw_limits.get("memory_bytes"), "execution_limits_invalid", minimum=512 * 1024**2
        ),
        "cpus": _finite_number(raw_limits.get("cpus"), "execution_limits_invalid", minimum=0.1),
        "shm_bytes": _integer(
            raw_limits.get("shm_bytes"), "execution_limits_invalid", minimum=64 * 1024**2
        ),
    }

    used: set[str] = set()
    rendered_steps: list[dict[str, Any]] = []
    for step in boundary["command_semantics"]["steps"]:
        arguments = [_render(argument, bindings, used) for argument in step["arguments"]]
        environment = {
            key: _render(value, bindings, used)
            for key, value in boundary["command_semantics"]["environment"].items()
        }
        rendered_steps.append(
            {
                "name": step["name"],
                "entrypoint": step["entrypoint"],
                "arguments": arguments,
                "environment": environment,
                "semantic_step_sha256": canonical_json_sha256(step),
            }
        )
    workspace = next(mount for mount in mounts if mount["role"] == "workspace")
    if workspace["container_path"] in map(PurePosixPath, bindings.values()):
        used.update(
            placeholder
            for placeholder, target in bindings.items()
            if PurePosixPath(target) == workspace["container_path"]
        )
    if used != set(bindings):
        _fail("execution_placeholder_binding_unused")

    return {
        "runtime": runtime,
        "image": image,
        "image_digest_sha256": image_match.group(1),
        "user": user,
        "mounts": sorted(mounts, key=lambda item: item["role"]),
        "mount_evidence": sorted(evidence, key=lambda item: item["role"]),
        "signatures": signatures,
        "workspace_host": workspace["host_path"],
        "workspace_container": workspace["container_path"],
        "tmpfs": sorted(tmpfs, key=lambda item: str(item["container_path"])),
        "limits": limits,
        "minimum_free_bytes": minimum_free_bytes,
        "rendered_steps": rendered_steps,
    }


def _capture(command: list[str]) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(command, check=False, capture_output=True, text=True)
    except OSError as exc:
        raise OciExecutionFailure("oci_runtime_unavailable") from exc


def _execute(command: list[str], stdout_handle: Any, stderr_handle: Any) -> int:
    try:
        return subprocess.run(
            command,
            check=False,
            stdin=subprocess.DEVNULL,
            stdout=stdout_handle,
            stderr=stderr_handle,
        ).returncode
    except OSError as exc:
        raise OciExecutionFailure("oci_runtime_unavailable") from exc


def _oci_command(specification: Mapping[str, Any], step: Mapping[str, Any]) -> list[str]:
    command = [
        str(specification["runtime"]),
        "run",
        "--rm",
        "--pull",
        "never",
        "--network",
        "none",
        "--hostname",
        "gab-caller",
        "--add-host",
        "gab-caller:127.0.0.1",
        "--read-only",
        "--cap-drop",
        "ALL",
        "--security-opt",
        "no-new-privileges",
        "--user",
        str(specification["user"]),
        "--workdir",
        str(specification["workspace_container"]),
        "--pids-limit",
        str(specification["limits"]["pids"]),
        "--memory",
        str(specification["limits"]["memory_bytes"]),
        "--cpus",
        format(specification["limits"]["cpus"], "g"),
        "--shm-size",
        str(specification["limits"]["shm_bytes"]),
    ]
    for mount in specification["mounts"]:
        option = f"type=bind,src={mount['host_path']},dst={mount['container_path']}" + (
            ",readonly" if mount["mode"] == "ro" else ""
        )
        command.extend(["--mount", option])
    for item in specification["tmpfs"]:
        command.extend(
            [
                "--tmpfs",
                f"{item['container_path']}:rw,nosuid,nodev,noexec,size={item['bytes']}",
            ]
        )
    for key, value in sorted(step["environment"].items()):
        command.extend(["--env", f"{key}={value}"])
    command.extend(
        [
            "--entrypoint",
            str(step["entrypoint"]),
            str(specification["image"]),
            *(str(argument) for argument in step["arguments"]),
        ]
    )
    return command


def _signatures_unchanged(signatures: Mapping[Path, tuple[int, int, int, int, int]]) -> bool:
    try:
        return all(_stat_signature(path) == expected for path, expected in signatures.items())
    except OSError:
        return False


def _report_base(boundary: Mapping[str, Any], specification: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": EXECUTION_SCHEMA_VERSION,
        "caller_run_identity_sha256": boundary["caller_run_identity_sha256"],
        "source_alignment_sha256": boundary["source_alignment_sha256"],
        "reference_sha256": boundary["reference_sha256"],
        "resource_manifest_sha256": boundary["resource_manifest_sha256"],
        "tool": boundary["tool"],
        "runtime": {
            "name": specification["runtime"],
            "version_output_sha256": None,
            "image_digest_sha256": specification["image_digest_sha256"],
            "image_digest_verified": False,
        },
        "security": {
            "pull_never": True,
            "network_none": True,
            "fixed_loopback_hostname": True,
            "read_only_root": True,
            "cap_drop_all": True,
            "no_new_privileges": True,
            "non_root_numeric_user": True,
            "stdin_disabled": True,
            "tmpfs_count": len(specification["tmpfs"]),
            "tmpfs_total_bytes": sum(item["bytes"] for item in specification["tmpfs"]),
            "tmpfs_all_noexec": True,
            "pids_limit": specification["limits"]["pids"],
            "memory_limit_bytes": specification["limits"]["memory_bytes"],
            "cpu_limit": specification["limits"]["cpus"],
            "shm_limit_bytes": specification["limits"]["shm_bytes"],
        },
        "input_verification": {
            "mounts": specification["mount_evidence"],
            "minimum_free_bytes": specification["minimum_free_bytes"],
            "free_space_gate_passed": True,
            "post_execution_stat_unchanged": False,
        },
        "declared_steps": len(specification["rendered_steps"]),
        "steps": [],
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }


def _finalize_report(
    base: dict[str, Any],
    *,
    status: str,
    error_code: str | None,
    output_path: Path,
    forbidden_bindings: Sequence[str],
) -> dict[str, Any]:
    report = {**base, "status": status, "error_code": error_code}
    report["execution_identity_sha256"] = canonical_json_sha256(report)
    normalized = validate_oci_execution_report_payload(
        report, forbidden_bindings=forbidden_bindings
    )
    _atomic_json(output_path, normalized)
    return normalized


def validate_oci_execution_report_payload(
    value: Mapping[str, Any], *, forbidden_bindings: Sequence[str] = ()
) -> dict[str, Any]:
    expected = {
        "schema_version",
        "status",
        "error_code",
        "caller_run_identity_sha256",
        "source_alignment_sha256",
        "reference_sha256",
        "resource_manifest_sha256",
        "tool",
        "runtime",
        "security",
        "input_verification",
        "declared_steps",
        "steps",
        "logical_subject_count",
        "subject_identifier_persisted",
        "private_paths_persisted",
        "private_manifest_digest_persisted",
        "execution_identity_sha256",
    }
    if set(value) != expected or value.get("schema_version") != EXECUTION_SCHEMA_VERSION:
        _fail("execution_report_schema_invalid")
    status = value.get("status")
    if status not in {"PASS", "FAIL_RETAINED"}:
        _fail("execution_report_status_invalid")
    if (
        value.get("logical_subject_count") != 1
        or value.get("subject_identifier_persisted") is not False
        or value.get("private_paths_persisted") is not False
        or value.get("private_manifest_digest_persisted") is not False
    ):
        _fail("execution_report_privacy_boundary_invalid")
    for key in (
        "caller_run_identity_sha256",
        "source_alignment_sha256",
        "reference_sha256",
        "resource_manifest_sha256",
    ):
        _sha(value.get(key), "execution_report_digest_invalid")
    tool = value.get("tool")
    runtime = value.get("runtime")
    security = value.get("security")
    verification = value.get("input_verification")
    steps = value.get("steps")
    declared_steps = _integer(
        value.get("declared_steps"), "execution_report_declared_steps_invalid", minimum=1
    )
    if not isinstance(tool, dict) or set(tool) != {
        "execution_kind",
        "logical_name",
        "identity_sha256",
        "network_disabled",
    }:
        _fail("execution_report_tool_invalid")
    if not isinstance(runtime, dict) or set(runtime) != {
        "name",
        "version_output_sha256",
        "image_digest_sha256",
        "image_digest_verified",
    }:
        _fail("execution_report_runtime_invalid")
    if runtime.get("name") not in RUNTIMES:
        _fail("execution_report_runtime_invalid")
    if runtime.get("version_output_sha256") is not None:
        _sha(runtime.get("version_output_sha256"), "execution_report_runtime_invalid")
    _sha(runtime.get("image_digest_sha256"), "execution_report_runtime_invalid")
    if not isinstance(runtime.get("image_digest_verified"), bool):
        _fail("execution_report_runtime_invalid")
    if (
        tool.get("execution_kind") != "oci_image"
        or tool.get("network_disabled") is not True
        or _sha(tool.get("identity_sha256"), "execution_report_tool_invalid")
        != runtime.get("image_digest_sha256")
    ):
        _fail("execution_report_tool_invalid")
    if not isinstance(security, dict) or set(security) != {
        "pull_never",
        "network_none",
        "fixed_loopback_hostname",
        "read_only_root",
        "cap_drop_all",
        "no_new_privileges",
        "non_root_numeric_user",
        "stdin_disabled",
        "tmpfs_count",
        "tmpfs_total_bytes",
        "tmpfs_all_noexec",
        "pids_limit",
        "memory_limit_bytes",
        "cpu_limit",
        "shm_limit_bytes",
    }:
        _fail("execution_report_security_invalid")
    for claim in (
        "pull_never",
        "network_none",
        "fixed_loopback_hostname",
        "read_only_root",
        "cap_drop_all",
        "no_new_privileges",
        "non_root_numeric_user",
        "stdin_disabled",
        "tmpfs_all_noexec",
    ):
        if security.get(claim) is not True:
            _fail("execution_report_security_invalid")
    _integer(security.get("tmpfs_count"), "execution_report_security_invalid", minimum=1)
    _integer(security.get("tmpfs_total_bytes"), "execution_report_security_invalid", minimum=1)
    _integer(security.get("pids_limit"), "execution_report_security_invalid", minimum=1)
    _integer(security.get("memory_limit_bytes"), "execution_report_security_invalid", minimum=1)
    _finite_number(security.get("cpu_limit"), "execution_report_security_invalid", minimum=0.1)
    _integer(security.get("shm_limit_bytes"), "execution_report_security_invalid", minimum=1)
    if not isinstance(verification, dict) or set(verification) != {
        "mounts",
        "minimum_free_bytes",
        "free_space_gate_passed",
        "post_execution_stat_unchanged",
    }:
        _fail("execution_report_input_verification_invalid")
    mounts = verification.get("mounts")
    if not isinstance(mounts, list) or not mounts:
        _fail("execution_report_input_verification_invalid")
    roles: set[str] = set()
    normalized_mount_order: list[str] = []
    for mount in mounts:
        if not isinstance(mount, dict) or set(mount) != {
            "role",
            "kind",
            "mode",
            "sha256",
            "bytes",
            "tree_manifest_sha256",
            "tree_identity_sha256",
            "tree_files",
            "tree_directories",
        }:
            _fail("execution_report_mount_invalid")
        role = mount.get("role")
        if not isinstance(role, str) or MOUNT_ROLE.fullmatch(role) is None or role in roles:
            _fail("execution_report_mount_invalid")
        roles.add(role)
        if mount.get("kind") not in {"file", "directory"} or mount.get("mode") not in {
            "ro",
            "rw",
        }:
            _fail("execution_report_mount_invalid")
        normalized_mount_order.append(role)
        for digest_key in ("sha256", "tree_manifest_sha256", "tree_identity_sha256"):
            if mount.get(digest_key) is not None:
                _sha(mount.get(digest_key), "execution_report_mount_invalid")
        for integer_key in ("bytes", "tree_files", "tree_directories"):
            if mount.get(integer_key) is not None:
                _integer(mount.get(integer_key), "execution_report_mount_invalid")
        if mount["kind"] == "file":
            if (
                mount["mode"] != "ro"
                or mount["sha256"] is None
                or mount["bytes"] is None
                or any(
                    mount[key] is not None
                    for key in (
                        "tree_manifest_sha256",
                        "tree_identity_sha256",
                        "tree_files",
                        "tree_directories",
                    )
                )
            ):
                _fail("execution_report_mount_invalid")
        elif mount["mode"] == "ro":
            if (
                mount["sha256"] is not None
                or mount["bytes"] is None
                or mount["tree_manifest_sha256"] is None
                or mount["tree_identity_sha256"] is None
                or mount["tree_files"] is None
                or mount["tree_directories"] is None
            ):
                _fail("execution_report_mount_invalid")
        elif role != "workspace" or any(
            mount[key] is not None
            for key in (
                "sha256",
                "bytes",
                "tree_manifest_sha256",
                "tree_identity_sha256",
                "tree_files",
                "tree_directories",
            )
        ):
            _fail("execution_report_mount_invalid")
    if normalized_mount_order != sorted(normalized_mount_order):
        _fail("execution_report_mount_invalid")
    if (
        roles.intersection({"selected_alignment", "reference_fasta", "resource_manifest"})
        != {
            "selected_alignment",
            "reference_fasta",
            "resource_manifest",
        }
        or "workspace" not in roles
        or sum(mount["mode"] == "rw" for mount in mounts) != 1
    ):
        _fail("execution_report_mount_invalid")
    by_role = {mount["role"]: mount for mount in mounts}
    if (
        by_role["selected_alignment"]["sha256"] != value["source_alignment_sha256"]
        or by_role["reference_fasta"]["sha256"] != value["reference_sha256"]
        or by_role["resource_manifest"]["sha256"] != value["resource_manifest_sha256"]
    ):
        _fail("execution_report_mount_invalid")
    _integer(
        verification.get("minimum_free_bytes"),
        "execution_report_input_verification_invalid",
    )
    if verification.get("free_space_gate_passed") is not True or not isinstance(
        verification.get("post_execution_stat_unchanged"), bool
    ):
        _fail("execution_report_input_verification_invalid")
    if not isinstance(steps, list) or len(steps) > declared_steps:
        _fail("execution_report_steps_invalid")
    for ordinal, step in enumerate(steps, start=1):
        if not isinstance(step, dict) or set(step) != {
            "ordinal",
            "name",
            "semantic_step_sha256",
            "exit_code",
            "stdout_sha256",
            "stdout_bytes",
            "stderr_sha256",
            "stderr_bytes",
        }:
            _fail("execution_report_step_invalid")
        if (
            step.get("ordinal") != ordinal
            or not isinstance(step.get("name"), str)
            or STEP_NAME.fullmatch(step["name"]) is None
        ):
            _fail("execution_report_step_invalid")
        _sha(step.get("semantic_step_sha256"), "execution_report_step_invalid")
        _signed_integer(step.get("exit_code"), "execution_report_step_invalid")
        _sha(step.get("stdout_sha256"), "execution_report_step_invalid")
        _sha(step.get("stderr_sha256"), "execution_report_step_invalid")
        _integer(step.get("stdout_bytes"), "execution_report_step_invalid")
        _integer(step.get("stderr_bytes"), "execution_report_step_invalid")
    error_code = value.get("error_code")
    if status == "PASS":
        if (
            error_code is not None
            or len(steps) != declared_steps
            or any(step["exit_code"] != 0 for step in steps)
            or runtime.get("image_digest_verified") is not True
            or verification.get("post_execution_stat_unchanged") is not True
        ):
            _fail("execution_report_pass_claim_invalid")
    elif not isinstance(error_code, str) or not error_code:
        _fail("execution_report_failure_claim_invalid")
    identity = _sha(value.get("execution_identity_sha256"), "execution_report_identity_invalid")
    payload = {key: item for key, item in value.items() if key != "execution_identity_sha256"}
    if identity != canonical_json_sha256(payload):
        _fail("execution_report_identity_changed")
    normalized = dict(value)
    _assert_portable(normalized, forbidden_bindings)
    return normalized


def load_oci_execution_report(
    path: str | Path, *, forbidden_bindings: Sequence[str] = ()
) -> tuple[dict[str, Any], str]:
    report_path = Path(path)
    value = _load_json(report_path, "execution_report_unreadable")
    normalized = validate_oci_execution_report_payload(
        value, forbidden_bindings=(*forbidden_bindings, str(report_path))
    )
    return normalized, sha256_file(report_path)


def run_oci_caller(
    boundary_checkpoint_path: str | Path,
    execution_spec_path: str | Path,
    output_report_path: str | Path,
    private_log_dir: str | Path,
    *,
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    """Run all prepared caller steps and retain one portable execution report."""

    boundary_path = Path(boundary_checkpoint_path)
    spec_path = Path(execution_spec_path)
    output_path = Path(output_report_path)
    log_dir = Path(private_log_dir)
    if output_path.exists():
        _fail("execution_report_already_exists")
    if log_dir.exists():
        _fail("execution_log_directory_already_exists")
    checkpoint, _ = load_prepared_caller_run_boundary(
        boundary_path, forbidden_bindings=forbidden_bindings
    )
    boundary = checkpoint["run_boundary"]
    raw_spec = _load_json(spec_path, "execution_spec_unreadable")
    specification = _validate_execution_spec(raw_spec, boundary)
    private = (
        *forbidden_bindings,
        str(boundary_path),
        str(spec_path),
        str(output_path),
        str(log_dir),
        *(str(mount["host_path"]) for mount in specification["mounts"]),
        *(
            str(mount["tree_manifest_path"])
            for mount in specification["mounts"]
            if "tree_manifest_path" in mount
        ),
    )
    base = _report_base(boundary, specification)
    _assert_portable(base, private)
    log_dir.mkdir(parents=True, exist_ok=False)
    try:
        version = _capture([specification["runtime"], "--version"])
        if version.returncode:
            raise OciExecutionFailure("oci_runtime_version_query_failed")
        base["runtime"]["version_output_sha256"] = hashlib.sha256(
            (version.stdout + version.stderr).encode()
        ).hexdigest()
        inspected = _capture(
            [
                specification["runtime"],
                "image",
                "inspect",
                "--format",
                "{{json .RepoDigests}}",
                specification["image"],
            ]
        )
        if inspected.returncode:
            raise OciExecutionFailure("oci_image_not_present_locally")
        try:
            repo_digests = json.loads(inspected.stdout)
        except json.JSONDecodeError as exc:
            raise OciExecutionFailure("oci_image_inspect_invalid") from exc
        expected_suffix = f"@sha256:{specification['image_digest_sha256']}"
        if not isinstance(repo_digests, list) or not any(
            isinstance(item, str) and item.endswith(expected_suffix) for item in repo_digests
        ):
            raise OciExecutionFailure("oci_image_digest_not_verified")
        base["runtime"]["image_digest_verified"] = True

        for ordinal, step in enumerate(specification["rendered_steps"], start=1):
            if (
                shutil.disk_usage(specification["workspace_host"]).free
                < specification["minimum_free_bytes"]
            ):
                raise OciExecutionFailure("free_space_below_execution_hard_stop")
            stdout_path = log_dir / f"{ordinal:03d}-{step['name']}.stdout.log"
            stderr_path = log_dir / f"{ordinal:03d}-{step['name']}.stderr.log"
            command = _oci_command(specification, step)
            with stdout_path.open("xb") as stdout_handle, stderr_path.open("xb") as stderr_handle:
                exit_code = _execute(command, stdout_handle, stderr_handle)
                stdout_handle.flush()
                stderr_handle.flush()
                os.fsync(stdout_handle.fileno())
                os.fsync(stderr_handle.fileno())
            base["steps"].append(
                {
                    "ordinal": ordinal,
                    "name": step["name"],
                    "semantic_step_sha256": step["semantic_step_sha256"],
                    "exit_code": exit_code,
                    "stdout_sha256": sha256_file(stdout_path),
                    "stdout_bytes": stdout_path.stat().st_size,
                    "stderr_sha256": sha256_file(stderr_path),
                    "stderr_bytes": stderr_path.stat().st_size,
                }
            )
            if exit_code:
                raise OciExecutionFailure("caller_step_failed")
        if not _signatures_unchanged(specification["signatures"]):
            raise OciExecutionFailure("read_only_input_changed_during_execution")
        base["input_verification"]["post_execution_stat_unchanged"] = True
    except OciExecutionFailure as exc:
        _finalize_report(
            base,
            status="FAIL_RETAINED",
            error_code=exc.code,
            output_path=output_path,
            forbidden_bindings=private,
        )
        raise OciExecutionFailure(exc.code, report_written=True) from exc
    except Exception as exc:
        _finalize_report(
            base,
            status="FAIL_RETAINED",
            error_code="unexpected_execution_failure",
            output_path=output_path,
            forbidden_bindings=private,
        )
        raise OciExecutionFailure("unexpected_execution_failure", report_written=True) from exc
    return _finalize_report(
        base,
        status="PASS",
        error_code=None,
        output_path=output_path,
        forbidden_bindings=private,
    )


def _read_forbidden_bindings(path: Path | None) -> tuple[str, ...]:
    if path is None:
        return ()
    _regular_file(path, "forbidden_bindings_unreadable")
    try:
        return tuple(
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    except (OSError, UnicodeError) as exc:
        raise OciExecutionFailure("forbidden_bindings_unreadable") from exc


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--boundary", type=Path, required=True)
    run_parser.add_argument("--execution-spec", type=Path, required=True)
    run_parser.add_argument("--report", type=Path, required=True)
    run_parser.add_argument("--log-dir", type=Path, required=True)
    run_parser.add_argument("--forbidden-bindings-file", type=Path)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--report", type=Path, required=True)
    validate_parser.add_argument("--forbidden-bindings-file", type=Path)
    args = parser.parse_args()
    try:
        forbidden = _read_forbidden_bindings(args.forbidden_bindings_file)
        if args.command == "run":
            report = run_oci_caller(
                args.boundary,
                args.execution_spec,
                args.report,
                args.log_dir,
                forbidden_bindings=(
                    *forbidden,
                    *((str(args.forbidden_bindings_file),) if args.forbidden_bindings_file else ()),
                ),
            )
            file_sha256 = sha256_file(args.report)
        else:
            report, file_sha256 = load_oci_execution_report(
                args.report, forbidden_bindings=forbidden
            )
        print(
            json.dumps(
                {
                    "schema_version": EXECUTION_SCHEMA_VERSION,
                    "status": report["status"],
                    "caller_run_identity_sha256": report["caller_run_identity_sha256"],
                    "execution_identity_sha256": report["execution_identity_sha256"],
                    "execution_report_sha256": file_sha256,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    except (OciExecutionFailure, PortableEvidenceError):
        print(
            "OCI caller execution failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
