"""Compile public, path-free caller profiles into private restart-safe run plans."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

from .caller_receipt import (
    canonical_json_sha256,
    load_prepared_caller_run_boundary,
    prepare_caller_run_boundary,
    validate_caller_command_semantics,
    validate_caller_tool_identity,
    write_prepared_caller_run_boundary,
)
from .portable_model import (
    FORBIDDEN_PORTABLE_TEXT,
    SHA256,
    PortableEvidenceError,
    contains_forbidden_binding,
    sha256_file,
)

PROFILE_SCHEMA_VERSION = 1
PROFILE_ID = re.compile(r"^[a-z0-9][a-z0-9_.-]{0,95}$")
MOUNT_ROLE = re.compile(r"^[a-z][a-z0-9_.-]{0,95}$")
PLACEHOLDER = re.compile(r"\{[A-Za-z][A-Za-z0-9_.:-]{0,95}\}")
COMPONENT_ROLE = re.compile(r"^[a-z][a-z0-9_.-]{0,95}$")
OCI_IMAGE = re.compile(r"^[^\s,@]+@sha256:([0-9a-f]{64})$")
NON_ROOT_USER = re.compile(r"^[1-9][0-9]*:[1-9][0-9]*$")
RUNTIMES = frozenset({"docker", "podman"})
BANNED_CONTAINER_PREFIXES = (
    PurePosixPath("/proc"),
    PurePosixPath("/sys"),
    PurePosixPath("/dev"),
    PurePosixPath("/var/run"),
)
RUN_FILES = (
    "artifact-spec.json",
    "command-semantics.json",
    "execution-spec.json",
    "prepared-boundary.json",
    "tool-spec.json",
)


@dataclass(frozen=True)
class ValidatedCallerProfile:
    payload: dict[str, Any] = field(repr=False)
    identity_sha256: str

    def portable_summary(self) -> dict[str, Any]:
        return {
            "schema_version": PROFILE_SCHEMA_VERSION,
            "status": "PASS",
            "profile_id": self.payload["profile_id"],
            "profile_identity_sha256": self.identity_sha256,
            "caller_name": self.payload["caller_name"],
            "caller_version": self.payload["caller_version"],
            "gate_caller": self.payload["gate_caller"],
            "image_identity_sha256": self.payload["image"].rsplit("sha256:", 1)[1],
            "mounts": len(self.payload["mounts"]),
            "steps": len(self.payload["command_semantics"]["steps"]),
            "private_paths_persisted": False,
        }


def _fail(message: str) -> None:
    raise PortableEvidenceError(message)


def _integer(value: Any, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _fail(label)
    return value


def _finite_number(value: Any, label: str, *, minimum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(label)
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < minimum:
        _fail(label)
    return parsed


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        _fail(label)
    return value


def _load_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        _fail(label)
    try:
        if path.resolve(strict=True) != path.absolute():
            _fail(label)
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError(label) from exc
    if not isinstance(value, dict):
        _fail(label)
    return value


def _walk_strings(value: Any):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, child in value.items():
            yield from _walk_strings(key)
            yield from _walk_strings(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _walk_strings(child)


def _assert_public(value: Any, forbidden_bindings: Sequence[str]) -> None:
    forbidden = tuple(binding for binding in forbidden_bindings if binding)
    for text in _walk_strings(value):
        if FORBIDDEN_PORTABLE_TEXT.search(text):
            _fail("caller profile contains a private path or secret")
        if any(contains_forbidden_binding(text, binding) for binding in forbidden):
            _fail("caller profile contains a forbidden private binding")


def _container_path(value: Any, label: str) -> PurePosixPath:
    if not isinstance(value, str) or not value.startswith("/") or "," in value or ":" in value:
        _fail(label)
    path = PurePosixPath(value)
    if str(path) != value or value == "/" or ".." in path.parts:
        _fail(label)
    if any(path == prefix or prefix in path.parents for prefix in BANNED_CONTAINER_PREFIXES):
        _fail(label)
    return path


def _relative_output(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        _fail(label)
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or str(path) != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        _fail(label)
    return value


def _target_is_mounted(path: PurePosixPath, mounts: Sequence[Mapping[str, Any]]) -> bool:
    for mount in mounts:
        target = PurePosixPath(str(mount["container_path"]))
        if mount["kind"] == "file" and path == target:
            return True
        if mount["kind"] == "directory" and (path == target or target in path.parents):
            return True
    return False


def _validate_mounts(value: Any) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        _fail("caller profile mounts are invalid")
    mounts: list[dict[str, Any]] = []
    roles: set[str] = set()
    targets: set[PurePosixPath] = set()
    for raw in value:
        base_fields = {
            "role",
            "kind",
            "mode",
            "container_path",
        }
        if not isinstance(raw, dict) or frozenset(raw) not in {
            frozenset(base_fields),
            frozenset({*base_fields, "expected_sha256"}),
            frozenset({*base_fields, "expected_tree_manifest_sha256"}),
        }:
            _fail("caller profile mount schema is invalid")
        role = raw.get("role")
        kind = raw.get("kind")
        mode = raw.get("mode")
        if not isinstance(role, str) or MOUNT_ROLE.fullmatch(role) is None or role in roles:
            _fail("caller profile mount role is invalid or duplicated")
        if kind not in {"file", "directory"} or mode not in {"ro", "rw"}:
            _fail("caller profile mount kind or mode is invalid")
        if kind == "file" and mode != "ro":
            _fail("caller profile file mounts must be read-only")
        target = _container_path(raw.get("container_path"), "caller profile mount target invalid")
        if target in targets or any(
            target in existing.parents or existing in target.parents for existing in targets
        ):
            _fail("caller profile mount targets overlap")
        roles.add(role)
        targets.add(target)
        mount = {
            "role": role,
            "kind": kind,
            "mode": mode,
            "container_path": str(target),
        }
        if "expected_sha256" in raw:
            if kind != "file":
                _fail("caller profile digest constraint requires a file mount")
            mount["expected_sha256"] = _sha(
                raw.get("expected_sha256"),
                "caller profile file digest constraint is invalid",
            )
        if "expected_tree_manifest_sha256" in raw:
            if kind != "directory" or mode != "ro":
                _fail("caller profile tree constraint requires a read-only directory mount")
            mount["expected_tree_manifest_sha256"] = _sha(
                raw.get("expected_tree_manifest_sha256"),
                "caller profile tree-manifest digest constraint is invalid",
            )
        mounts.append(mount)
    required = {"selected_alignment", "reference_fasta", "resource_manifest", "workspace"}
    if not required.issubset(roles):
        _fail("caller profile lacks required mount roles")
    by_role = {mount["role"]: mount for mount in mounts}
    for role in ("selected_alignment", "reference_fasta", "resource_manifest"):
        if by_role[role]["kind"] != "file" or by_role[role]["mode"] != "ro":
            _fail("caller profile immutable primary mount is invalid")
    if by_role["workspace"]["kind"] != "directory" or by_role["workspace"]["mode"] != "rw":
        _fail("caller profile workspace mount is invalid")
    if sum(mount["mode"] == "rw" for mount in mounts) != 1:
        _fail("caller profile must have one writable workspace")
    return sorted(mounts, key=lambda mount: mount["role"])


def _validate_placeholders(
    value: Any, mounts: Sequence[Mapping[str, Any]], command: Mapping[str, Any]
) -> dict[str, str]:
    if not isinstance(value, dict) or not value:
        _fail("caller profile placeholder bindings are invalid")
    bindings: dict[str, str] = {}
    for placeholder, target_value in value.items():
        if not isinstance(placeholder, str) or PLACEHOLDER.fullmatch(placeholder) is None:
            _fail("caller profile placeholder name is invalid")
        target = _container_path(target_value, "caller profile placeholder target is invalid")
        if not _target_is_mounted(target, mounts):
            _fail("caller profile placeholder target is not mounted")
        bindings[placeholder] = str(target)
    used: set[str] = set()
    for step in command["steps"]:
        for text in (step["entrypoint"], *step["arguments"]):
            used.update(PLACEHOLDER.findall(text))
    for text in command["environment"].values():
        used.update(PLACEHOLDER.findall(text))
    if used != set(bindings):
        _fail("caller profile placeholder bindings are missing or unused")
    return dict(sorted(bindings.items()))


def _validate_tmpfs(value: Any, mounts: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        _fail("caller profile tmpfs policy is invalid")
    mount_targets = {PurePosixPath(str(mount["container_path"])) for mount in mounts}
    seen: set[PurePosixPath] = set()
    result: list[dict[str, Any]] = []
    for raw in value:
        if not isinstance(raw, dict) or set(raw) != {"container_path", "bytes", "executable"}:
            _fail("caller profile tmpfs policy is invalid")
        target = _container_path(raw.get("container_path"), "caller profile tmpfs target invalid")
        if (
            target in seen
            or any(target in other.parents or other in target.parents for other in seen)
            or target in mount_targets
            or any(target in other.parents or other in target.parents for other in mount_targets)
        ):
            _fail("caller profile tmpfs target overlaps another mount")
        byte_count = _integer(
            raw.get("bytes"), "caller profile tmpfs byte count is invalid", minimum=64 * 1024**2
        )
        if raw.get("executable") is not False:
            _fail("caller profile tmpfs must be noexec")
        seen.add(target)
        result.append({"container_path": str(target), "bytes": byte_count, "executable": False})
    if PurePosixPath("/tmp") not in seen:
        _fail("caller profile requires a private /tmp tmpfs")
    return sorted(result, key=lambda item: item["container_path"])


def _validate_limits(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {
        "pids",
        "memory_bytes",
        "cpus",
        "shm_bytes",
    }:
        _fail("caller profile resource limits are invalid")
    return {
        "pids": _integer(value.get("pids"), "caller profile PID limit is invalid", minimum=1),
        "memory_bytes": _integer(
            value.get("memory_bytes"),
            "caller profile memory limit is invalid",
            minimum=512 * 1024**2,
        ),
        "cpus": _finite_number(
            value.get("cpus"), "caller profile CPU limit is invalid", minimum=0.1
        ),
        "shm_bytes": _integer(
            value.get("shm_bytes"),
            "caller profile shared-memory limit is invalid",
            minimum=64 * 1024**2,
        ),
    }


def _validate_artifact(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {
        "raw_artifact_id",
        "identity_scheme",
        "components",
        "derived_payload",
    }:
        _fail("caller profile artifact contract is invalid")
    raw_artifact_id = value.get("raw_artifact_id")
    if (
        not isinstance(raw_artifact_id, str)
        or not raw_artifact_id
        or "/" in raw_artifact_id
        or "\\" in raw_artifact_id
        or "\n" in raw_artifact_id
    ):
        _fail("caller profile raw artifact ID is invalid")
    scheme = value.get("identity_scheme")
    if scheme not in {"file_sha256", "named_sha256_set"}:
        _fail("caller profile artifact identity scheme is invalid")
    if value.get("derived_payload") is not None:
        _fail("caller profile file artifact cannot contain a derived payload")
    components = value.get("components")
    if not isinstance(components, dict) or not components:
        _fail("caller profile artifact components are invalid")
    normalized_components: dict[str, str] = {}
    paths: set[str] = set()
    for role, relative_value in components.items():
        if not isinstance(role, str) or COMPONENT_ROLE.fullmatch(role) is None:
            _fail("caller profile artifact component role is invalid")
        relative = _relative_output(
            relative_value, "caller profile artifact component path is invalid"
        )
        if relative in paths:
            _fail("caller profile artifact component path is duplicated")
        paths.add(relative)
        normalized_components[role] = relative
    if scheme == "file_sha256" and len(normalized_components) != 1:
        _fail("caller profile file artifact requires exactly one component")
    return {
        "raw_artifact_id": raw_artifact_id,
        "identity_scheme": scheme,
        "components": dict(sorted(normalized_components.items())),
        "derived_payload": None,
    }


def validate_caller_profile_payload(
    value: Mapping[str, Any], *, forbidden_bindings: Sequence[str] = ()
) -> ValidatedCallerProfile:
    """Validate and canonicalize one reusable public caller profile."""

    expected = {
        "schema_version",
        "profile_id",
        "caller_name",
        "caller_version",
        "gate_caller",
        "image",
        "licence",
        "mounts",
        "placeholder_bindings",
        "command_semantics",
        "tmpfs",
        "limits",
        "artifact",
    }
    if set(value) != expected or value.get("schema_version") != PROFILE_SCHEMA_VERSION:
        _fail("caller profile schema is invalid")
    profile_id = value.get("profile_id")
    caller_name = value.get("caller_name")
    caller_version = value.get("caller_version")
    gate_caller = value.get("gate_caller")
    image = value.get("image")
    licence = value.get("licence")
    if not isinstance(profile_id, str) or PROFILE_ID.fullmatch(profile_id) is None:
        _fail("caller profile ID is invalid")
    if not isinstance(caller_name, str) or not caller_name.strip():
        _fail("caller profile name is invalid")
    if not isinstance(caller_version, str) or not caller_version.strip():
        _fail("caller profile version is invalid")
    if not isinstance(gate_caller, str) or PROFILE_ID.fullmatch(gate_caller) is None:
        _fail("caller profile gate caller is invalid")
    if not isinstance(image, str) or (image_match := OCI_IMAGE.fullmatch(image)) is None:
        _fail("caller profile image is not digest pinned")
    if not isinstance(licence, str) or not licence.strip():
        _fail("caller profile licence boundary is invalid")
    command_value = value.get("command_semantics")
    if not isinstance(command_value, Mapping):
        _fail("caller profile command semantics are invalid")
    command = validate_caller_command_semantics(command_value)
    if "profile_identity_sha256" in command["parameters"]:
        _fail("caller profile parameters use a reserved identity field")
    mounts = _validate_mounts(value.get("mounts"))
    placeholders = _validate_placeholders(value.get("placeholder_bindings"), mounts, command)
    normalized = {
        "schema_version": PROFILE_SCHEMA_VERSION,
        "profile_id": profile_id,
        "caller_name": caller_name,
        "caller_version": caller_version,
        "gate_caller": gate_caller,
        "image": image,
        "licence": licence,
        "mounts": mounts,
        "placeholder_bindings": placeholders,
        "command_semantics": command,
        "tmpfs": _validate_tmpfs(value.get("tmpfs"), mounts),
        "limits": _validate_limits(value.get("limits")),
        "artifact": _validate_artifact(value.get("artifact")),
    }
    tool = {
        "execution_kind": "oci_image",
        "logical_name": gate_caller,
        "identity_sha256": image_match.group(1),
        "network_disabled": True,
    }
    validate_caller_tool_identity(tool)
    _assert_public(normalized, forbidden_bindings)
    return ValidatedCallerProfile(
        payload=normalized,
        identity_sha256=canonical_json_sha256(normalized),
    )


def load_caller_profile(
    path: str | Path, *, forbidden_bindings: Sequence[str] = ()
) -> ValidatedCallerProfile:
    profile_path = Path(path)
    value = _load_json(profile_path, "caller profile is unreadable")
    return validate_caller_profile_payload(
        value, forbidden_bindings=(*forbidden_bindings, str(profile_path))
    )


def _host_path(value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value.startswith("/") or "," in value:
        _fail(label)
    return Path(value)


def _regular_file(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_file():
        _fail(label)
    try:
        if path.resolve(strict=True) != path.absolute():
            _fail(label)
    except OSError as exc:
        raise PortableEvidenceError(label) from exc


def _regular_directory(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_dir():
        _fail(label)
    try:
        if path.resolve(strict=True) != path.absolute():
            _fail(label)
    except OSError as exc:
        raise PortableEvidenceError(label) from exc


def _compile_mounts(
    profile_mounts: Sequence[Mapping[str, Any]], raw_bindings: Any
) -> tuple[list[dict[str, Any]], Path]:
    if not isinstance(raw_bindings, dict):
        _fail("caller run mount bindings are invalid")
    expected_roles = {str(mount["role"]) for mount in profile_mounts}
    if set(raw_bindings) != expected_roles:
        _fail("caller run mount bindings are missing or unexpected")
    execution_mounts: list[dict[str, Any]] = []
    workspace: Path | None = None
    read_only_paths: list[Path] = []
    for profile_mount in profile_mounts:
        role = str(profile_mount["role"])
        kind = str(profile_mount["kind"])
        mode = str(profile_mount["mode"])
        binding = raw_bindings[role]
        if not isinstance(binding, dict):
            _fail("caller run mount binding is invalid")
        if kind == "file":
            if set(binding) != {"host_path", "bytes", "sha256"}:
                _fail("caller run file binding schema is invalid")
            host = _host_path(binding.get("host_path"), "caller run file path is invalid")
            _regular_file(host, "caller run file is missing or symlinked")
            byte_count = _integer(binding.get("bytes"), "caller run file byte count is invalid")
            digest = _sha(binding.get("sha256"), "caller run file digest is invalid")
            expected_digest = profile_mount.get("expected_sha256")
            if expected_digest is not None and digest != expected_digest:
                _fail("caller run file differs from its profile digest constraint")
            if host.stat().st_size != byte_count:
                _fail("caller run file byte count changed")
            read_only_paths.append(host)
            execution_mounts.append(
                {
                    "role": role,
                    "kind": "file",
                    "host_path": str(host),
                    "container_path": profile_mount["container_path"],
                    "mode": "ro",
                    "bytes": byte_count,
                    "sha256": digest,
                    "tree_manifest_path": None,
                    "tree_manifest_bytes": None,
                    "tree_manifest_sha256": None,
                }
            )
        elif mode == "ro":
            if set(binding) != {
                "host_path",
                "tree_manifest_path",
                "tree_manifest_bytes",
                "tree_manifest_sha256",
            }:
                _fail("caller run resource-tree binding schema is invalid")
            host = _host_path(binding.get("host_path"), "caller run resource-tree path is invalid")
            manifest = _host_path(
                binding.get("tree_manifest_path"),
                "caller run resource-tree manifest path is invalid",
            )
            _regular_directory(host, "caller run resource tree is missing or symlinked")
            _regular_file(manifest, "caller run resource-tree manifest is missing or symlinked")
            manifest_bytes = _integer(
                binding.get("tree_manifest_bytes"),
                "caller run resource-tree manifest byte count is invalid",
                minimum=1,
            )
            manifest_sha256 = _sha(
                binding.get("tree_manifest_sha256"),
                "caller run resource-tree manifest digest is invalid",
            )
            expected_manifest_sha256 = profile_mount.get("expected_tree_manifest_sha256")
            if expected_manifest_sha256 is not None and manifest_sha256 != expected_manifest_sha256:
                _fail("caller run resource tree differs from its profile digest constraint")
            if (
                manifest.stat().st_size != manifest_bytes
                or sha256_file(manifest) != manifest_sha256
            ):
                _fail("caller run resource-tree manifest binding changed")
            read_only_paths.extend((host, manifest))
            execution_mounts.append(
                {
                    "role": role,
                    "kind": "directory",
                    "host_path": str(host),
                    "container_path": profile_mount["container_path"],
                    "mode": "ro",
                    "bytes": None,
                    "sha256": None,
                    "tree_manifest_path": str(manifest),
                    "tree_manifest_bytes": manifest_bytes,
                    "tree_manifest_sha256": manifest_sha256,
                }
            )
        else:
            if role != "workspace" or set(binding) != {"host_path"}:
                _fail("caller run workspace binding schema is invalid")
            host = _host_path(binding.get("host_path"), "caller run workspace path is invalid")
            _regular_directory(host, "caller run workspace is missing or symlinked")
            if any(host.iterdir()):
                _fail("caller run workspace must be empty")
            workspace = host
            execution_mounts.append(
                {
                    "role": role,
                    "kind": "directory",
                    "host_path": str(host),
                    "container_path": profile_mount["container_path"],
                    "mode": "rw",
                    "bytes": None,
                    "sha256": None,
                    "tree_manifest_path": None,
                    "tree_manifest_bytes": None,
                    "tree_manifest_sha256": None,
                }
            )
    assert workspace is not None
    for path in read_only_paths:
        if path == workspace or workspace in path.parents or path in workspace.parents:
            _fail("caller run workspace overlaps a read-only input")
    return execution_mounts, workspace


def _validate_run_bindings(
    value: Mapping[str, Any], profile: ValidatedCallerProfile
) -> tuple[dict[str, Any], Path]:
    if (
        set(value)
        != {
            "schema_version",
            "runtime",
            "user",
            "minimum_free_bytes",
            "mounts",
        }
        or value.get("schema_version") != PROFILE_SCHEMA_VERSION
    ):
        _fail("caller run binding schema is invalid")
    runtime = value.get("runtime")
    user = value.get("user")
    if runtime not in RUNTIMES:
        _fail("caller run OCI runtime is unsupported")
    if not isinstance(user, str) or NON_ROOT_USER.fullmatch(user) is None:
        _fail("caller run user must be non-root numeric uid:gid")
    mounts, workspace = _compile_mounts(profile.payload["mounts"], value.get("mounts"))
    return (
        {
            "schema_version": PROFILE_SCHEMA_VERSION,
            "runtime": runtime,
            "image": profile.payload["image"],
            "user": user,
            "mounts": mounts,
            "placeholder_bindings": profile.payload["placeholder_bindings"],
            "tmpfs": profile.payload["tmpfs"],
            "limits": profile.payload["limits"],
            "minimum_free_bytes": _integer(
                value.get("minimum_free_bytes"),
                "caller run free-space policy is invalid",
            ),
        },
        workspace,
    )


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    if path.exists():
        _fail("refusing to replace an existing caller run file")
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
            raise PortableEvidenceError("refusing to replace an existing caller run file") from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _fsync_directory(path: Path) -> None:
    """Persist directory entries used by the restart-safe completion boundary."""

    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _assert_output_separate(
    output: Path,
    workspace: Path,
    profile_path: Path,
    bindings_path: Path,
    gate_path: Path,
    resource_path: Path,
) -> None:
    absolute = output.absolute()
    compared = (
        workspace,
        profile_path.absolute(),
        bindings_path.absolute(),
        gate_path.absolute(),
        resource_path.absolute(),
    )
    for path in compared:
        if absolute == path or absolute in path.parents or path in absolute.parents:
            _fail("caller run-plan directory overlaps an input or workspace")


def prepare_caller_profile_run(
    profile_path: str | Path,
    binding_spec_path: str | Path,
    input_gate_report: str | Path,
    resource_manifest: str | Path,
    output_directory: str | Path,
    *,
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    """Create a private run directory and publish its completion marker last."""

    profile_file = Path(profile_path)
    bindings_file = Path(binding_spec_path)
    gate_file = Path(input_gate_report)
    resource_file = Path(resource_manifest)
    output = Path(output_directory)
    if output.exists():
        _fail("refusing to replace an existing caller run-plan directory")
    profile = load_caller_profile(profile_file, forbidden_bindings=forbidden_bindings)
    raw_bindings = _load_json(bindings_file, "caller run bindings are unreadable")
    execution_spec, workspace = _validate_run_bindings(raw_bindings, profile)
    _assert_output_separate(
        output, workspace, profile_file, bindings_file, gate_file, resource_file
    )
    resource_mount = next(
        mount for mount in execution_spec["mounts"] if mount["role"] == "resource_manifest"
    )
    if Path(resource_mount["host_path"]).absolute() != resource_file.absolute():
        _fail("caller run resource-manifest mount differs from the prepared boundary")
    parameters = deepcopy(profile.payload["command_semantics"]["parameters"])
    parameters["profile_id"] = profile.payload["profile_id"]
    parameters["profile_identity_sha256"] = profile.identity_sha256
    command = {
        **profile.payload["command_semantics"],
        "parameters": parameters,
    }
    tool = validate_caller_tool_identity(
        {
            "execution_kind": "oci_image",
            "logical_name": profile.payload["gate_caller"],
            "identity_sha256": profile.payload["image"].rsplit("sha256:", 1)[1],
            "network_disabled": True,
        }
    )
    expected_resource_sha256 = str(resource_mount["sha256"])
    boundary = prepare_caller_run_boundary(
        gate_file,
        resource_file,
        caller_name=profile.payload["caller_name"],
        caller_version=profile.payload["caller_version"],
        gate_caller=profile.payload["gate_caller"],
        expected_resource_manifest_sha256=expected_resource_sha256,
        tool=tool,
        command_semantics=command,
        forbidden_bindings=forbidden_bindings,
    )
    by_role = {mount["role"]: mount for mount in execution_spec["mounts"]}
    if (
        by_role["selected_alignment"]["sha256"] != boundary["source_alignment_sha256"]
        or by_role["reference_fasta"]["sha256"] != boundary["reference_sha256"]
        or by_role["resource_manifest"]["sha256"] != boundary["resource_manifest_sha256"]
        or by_role["resource_manifest"]["bytes"] != boundary["resource_manifest_bytes"]
    ):
        _fail("caller run mounts differ from the prepared input boundary")
    artifact = profile.payload["artifact"]
    artifact_spec = {
        "schema_version": PROFILE_SCHEMA_VERSION,
        "raw_artifact_id": artifact["raw_artifact_id"],
        "identity_scheme": artifact["identity_scheme"],
        "components": {
            role: str(workspace / PurePosixPath(relative))
            for role, relative in artifact["components"].items()
        },
        "derived_payload": None,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        if output.parent.resolve(strict=True) != output.parent.absolute():
            _fail("caller run-plan parent contains a symbolic path component")
    except OSError as exc:
        raise PortableEvidenceError("caller run-plan parent is unreadable") from exc
    output.mkdir(exist_ok=False)
    try:
        _atomic_json(output / "tool-spec.json", tool)
        _atomic_json(output / "command-semantics.json", command)
        _atomic_json(output / "execution-spec.json", execution_spec)
        _atomic_json(output / "artifact-spec.json", artifact_spec)
        checkpoint = write_prepared_caller_run_boundary(
            output / "prepared-boundary.json",
            gate_file,
            resource_file,
            caller_name=profile.payload["caller_name"],
            caller_version=profile.payload["caller_version"],
            gate_caller=profile.payload["gate_caller"],
            expected_resource_manifest_sha256=expected_resource_sha256,
            tool=tool,
            command_semantics=command,
            forbidden_bindings=forbidden_bindings,
        )
        if checkpoint["run_boundary"] != boundary:
            _fail("caller input boundary changed while the run plan was being prepared")
        files = [
            {
                "name": name,
                "bytes": (output / name).stat().st_size,
                "sha256": sha256_file(output / name),
            }
            for name in RUN_FILES
        ]
        marker: dict[str, Any] = {
            "schema_version": PROFILE_SCHEMA_VERSION,
            "status": "PREPARED",
            "scope": "private_runtime",
            "portable": False,
            "profile_id": profile.payload["profile_id"],
            "profile_identity_sha256": profile.identity_sha256,
            "caller_run_identity_sha256": checkpoint["run_boundary"]["caller_run_identity_sha256"],
            "files": files,
        }
        marker["prepared_run_identity_sha256"] = canonical_json_sha256(marker)
        _atomic_json(output / "PREPARED.json", marker)
        _fsync_directory(output)
        _fsync_directory(output.parent)
    except BaseException:
        raise
    return {
        "schema_version": PROFILE_SCHEMA_VERSION,
        "status": "PREPARED",
        "profile_id": profile.payload["profile_id"],
        "profile_identity_sha256": profile.identity_sha256,
        "caller_name": profile.payload["caller_name"],
        "caller_version": profile.payload["caller_version"],
        "caller_run_identity_sha256": boundary["caller_run_identity_sha256"],
        "run_files": len(RUN_FILES),
        "subject_identifier_persisted": False,
        "run_plan_scope": "private_runtime",
        "private_paths_persisted_in_run_plan": True,
        "private_paths_persisted_in_portable_summary": False,
    }


def load_prepared_caller_profile_run(
    run_directory: str | Path, *, forbidden_bindings: Sequence[str] = ()
) -> dict[str, Any]:
    """Verify the final marker and every prepared run file before reuse."""

    run_dir = Path(run_directory)
    _regular_directory(run_dir, "prepared caller run directory is invalid")
    marker = _load_json(run_dir / "PREPARED.json", "prepared caller run marker is unreadable")
    expected = {
        "schema_version",
        "status",
        "scope",
        "portable",
        "profile_id",
        "profile_identity_sha256",
        "caller_run_identity_sha256",
        "files",
        "prepared_run_identity_sha256",
    }
    if (
        set(marker) != expected
        or marker.get("schema_version") != PROFILE_SCHEMA_VERSION
        or marker.get("status") != "PREPARED"
        or marker.get("scope") != "private_runtime"
        or marker.get("portable") is not False
    ):
        _fail("prepared caller run marker schema is invalid")
    if (
        not isinstance(marker.get("profile_id"), str)
        or PROFILE_ID.fullmatch(str(marker["profile_id"])) is None
    ):
        _fail("prepared caller run profile ID is invalid")
    _sha(marker.get("profile_identity_sha256"), "prepared caller profile identity is invalid")
    _sha(marker.get("caller_run_identity_sha256"), "prepared caller run identity is invalid")
    raw_files = marker.get("files")
    if not isinstance(raw_files, list) or len(raw_files) != len(RUN_FILES):
        _fail("prepared caller run file manifest is invalid")
    normalized_files: list[dict[str, Any]] = []
    for item in raw_files:
        if not isinstance(item, dict) or set(item) != {"name", "bytes", "sha256"}:
            _fail("prepared caller run file entry is invalid")
        name = item.get("name")
        if name not in RUN_FILES:
            _fail("prepared caller run file name is invalid")
        path = run_dir / str(name)
        _regular_file(path, "prepared caller run file is invalid")
        byte_count = _integer(item.get("bytes"), "prepared caller run file size is invalid")
        digest = _sha(item.get("sha256"), "prepared caller run file digest is invalid")
        if path.stat().st_size != byte_count or sha256_file(path) != digest:
            _fail("prepared caller run file changed")
        normalized_files.append({"name": name, "bytes": byte_count, "sha256": digest})
    if [item["name"] for item in normalized_files] != list(RUN_FILES):
        _fail("prepared caller run file manifest is not canonical")
    identity = _sha(
        marker.get("prepared_run_identity_sha256"),
        "prepared caller run marker identity is invalid",
    )
    without_identity = {
        key: value for key, value in marker.items() if key != "prepared_run_identity_sha256"
    }
    if identity != canonical_json_sha256(without_identity):
        _fail("prepared caller run marker identity changed")
    checkpoint, _ = load_prepared_caller_run_boundary(
        run_dir / "prepared-boundary.json", forbidden_bindings=forbidden_bindings
    )
    tool = _load_json(run_dir / "tool-spec.json", "prepared caller tool is unreadable")
    command = _load_json(
        run_dir / "command-semantics.json", "prepared caller command is unreadable"
    )
    if (
        validate_caller_tool_identity(tool) != checkpoint["run_boundary"]["tool"]
        or validate_caller_command_semantics(command)
        != checkpoint["run_boundary"]["command_semantics"]
        or marker["caller_run_identity_sha256"]
        != checkpoint["run_boundary"]["caller_run_identity_sha256"]
        or command["parameters"].get("profile_id") != marker["profile_id"]
        or command["parameters"].get("profile_identity_sha256") != marker["profile_identity_sha256"]
    ):
        _fail("prepared caller run boundary differs from its profile marker")
    _assert_public(
        marker,
        (*forbidden_bindings, str(run_dir), *(str(run_dir / name) for name in RUN_FILES)),
    )
    return {
        "schema_version": PROFILE_SCHEMA_VERSION,
        "status": "PREPARED",
        "profile_id": marker["profile_id"],
        "profile_identity_sha256": marker["profile_identity_sha256"],
        "caller_run_identity_sha256": marker["caller_run_identity_sha256"],
        "run_files": len(normalized_files),
        "subject_identifier_persisted": False,
        "run_plan_scope": "private_runtime",
        "private_paths_persisted_in_run_plan": True,
        "private_paths_persisted_in_portable_summary": False,
    }


def _read_forbidden_bindings(path: Path | None) -> tuple[str, ...]:
    if path is None:
        return ()
    if path.is_symlink() or not path.is_file():
        _fail("forbidden-binding file is unreadable")
    try:
        return tuple(
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("forbidden-binding file is unreadable") from exc


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate_profile_parser = subparsers.add_parser("validate-profile")
    validate_profile_source = validate_profile_parser.add_mutually_exclusive_group(required=True)
    validate_profile_source.add_argument("--profile", type=Path)
    validate_profile_source.add_argument("--profile-id")
    validate_profile_parser.add_argument("--forbidden-bindings-file", type=Path)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_profile_source = prepare_parser.add_mutually_exclusive_group(required=True)
    prepare_profile_source.add_argument("--profile", type=Path)
    prepare_profile_source.add_argument("--profile-id")
    prepare_parser.add_argument("--bindings", type=Path, required=True)
    prepare_parser.add_argument("--input-gate", type=Path, required=True)
    prepare_parser.add_argument("--resource-manifest", type=Path, required=True)
    prepare_parser.add_argument("--output-directory", type=Path, required=True)
    prepare_parser.add_argument("--forbidden-bindings-file", type=Path)
    validate_run_parser = subparsers.add_parser("validate-run")
    validate_run_parser.add_argument("--run-directory", type=Path, required=True)
    validate_run_parser.add_argument("--forbidden-bindings-file", type=Path)
    args = parser.parse_args()
    try:
        forbidden = _read_forbidden_bindings(args.forbidden_bindings_file)
        if args.command in {"validate-profile", "prepare"}:
            if args.profile_id is not None:
                from .profile_catalog import resolve_shipped_profile

                profile_path = resolve_shipped_profile(args.profile_id)
            else:
                profile_path = args.profile
        if args.command == "validate-profile":
            result = load_caller_profile(
                profile_path,
                forbidden_bindings=(
                    *forbidden,
                    *((str(args.forbidden_bindings_file),) if args.forbidden_bindings_file else ()),
                ),
            ).portable_summary()
        elif args.command == "prepare":
            result = prepare_caller_profile_run(
                profile_path,
                args.bindings,
                args.input_gate,
                args.resource_manifest,
                args.output_directory,
                forbidden_bindings=(
                    *forbidden,
                    *((str(args.forbidden_bindings_file),) if args.forbidden_bindings_file else ()),
                ),
            )
        else:
            result = load_prepared_caller_profile_run(
                args.run_directory,
                forbidden_bindings=(
                    *forbidden,
                    *((str(args.forbidden_bindings_file),) if args.forbidden_bindings_file else ()),
                ),
            )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except Exception:
        print(
            "caller profile operation failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
