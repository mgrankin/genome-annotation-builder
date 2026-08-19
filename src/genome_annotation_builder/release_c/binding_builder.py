"""Build a private caller binding spec from generic profile mount roles.

Dataset identifiers and host paths are read only from ignored runtime files.
The CLI summary is deliberately path-free and identifier-free.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .caller_profile import ValidatedCallerProfile, load_caller_profile
from .portable_model import PortableEvidenceError, contains_forbidden_binding, sha256_file
from .profile_catalog import resolve_shipped_profile

BINDING_SCHEMA_VERSION = 1
BUILTIN_FILE_ROLES = frozenset(
    {
        "reference_dictionary",
        "reference_fai",
        "reference_fasta",
        "resource_manifest",
        "selected_alignment",
        "selected_alignment_index",
        "source_baf_vcf",
        "source_baf_vcf_index",
    }
)


class BindingBuilderError(PortableEvidenceError):
    """A private binding specification is incomplete or inconsistent."""


def _fail(message: str) -> None:
    raise BindingBuilderError(message)


def _load_object(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        _fail(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BindingBuilderError(f"{label} is unreadable") from exc
    if not isinstance(value, dict):
        _fail(f"{label} must be a JSON object")
    return value


def _atomic_json_no_replace(path: Path, value: object) -> None:
    if path.exists():
        _fail("refusing to replace an existing caller binding specification")
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
            raise BindingBuilderError(
                "refusing to replace an existing caller binding specification"
            ) from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _resolved_regular_file(path: Path, label: str) -> Path:
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        raise BindingBuilderError(f"{label} is not a regular file") from exc
    if resolved.is_symlink() or not resolved.is_file():
        _fail(f"{label} is not a regular file")
    return resolved


def _resolved_regular_directory(path: Path, label: str) -> Path:
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        raise BindingBuilderError(f"{label} is not a regular directory") from exc
    if resolved.is_symlink() or not resolved.is_dir():
        _fail(f"{label} is not a regular directory")
    return resolved


def _manifest_file_binding(spec: Mapping[str, Any], *, path_key: str = "path") -> dict[str, Any]:
    try:
        path = Path(str(spec[path_key]))
        expected_bytes = int(spec["bytes"])
        expected_sha256 = str(spec["sha256"])
    except (KeyError, TypeError, ValueError) as exc:
        raise BindingBuilderError("manifest file binding is incomplete") from exc
    path = _resolved_regular_file(path, "manifest-bound input")
    if path.stat().st_size != expected_bytes:
        _fail("manifest-bound input byte count changed")
    return {
        "host_path": str(path),
        "bytes": expected_bytes,
        "sha256": expected_sha256,
    }


def _manifest_index_binding(spec: Mapping[str, Any]) -> dict[str, Any]:
    try:
        path = Path(str(spec["index_path"]))
        expected_bytes = int(spec["index_bytes"])
        expected_sha256 = str(spec["index_sha256"])
    except (KeyError, TypeError, ValueError) as exc:
        raise BindingBuilderError("manifest index binding is incomplete") from exc
    path = _resolved_regular_file(path, "manifest-bound index")
    if path.stat().st_size != expected_bytes:
        _fail("manifest-bound index byte count changed")
    return {
        "host_path": str(path),
        "bytes": expected_bytes,
        "sha256": expected_sha256,
    }


def _reference_binding(reference: Mapping[str, Any], role: str) -> dict[str, Any]:
    keys = {
        "reference_fasta": ("fasta_path", "fasta_sha256"),
        "reference_fai": ("fai_path", "fai_sha256"),
        "reference_dictionary": ("dictionary_path", "dictionary_sha256"),
    }
    path_key, digest_key = keys[role]
    try:
        path = Path(str(reference[path_key]))
        digest = str(reference[digest_key])
    except KeyError as exc:
        raise BindingBuilderError("reference binding is incomplete") from exc
    path = _resolved_regular_file(path, role)
    return {"host_path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def _resource_file_binding(spec: Mapping[str, Any], role: str) -> dict[str, Any]:
    if set(spec) != {"path"}:
        _fail(f"resource file role {role} has the wrong fields")
    path = _resolved_regular_file(Path(str(spec["path"])), f"resource role {role}")
    return {
        "host_path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _resource_tree_binding(spec: Mapping[str, Any], role: str) -> dict[str, Any]:
    if set(spec) != {"path", "tree_manifest_path"}:
        _fail(f"resource tree role {role} has the wrong fields")
    root = _resolved_regular_directory(Path(str(spec["path"])), f"resource role {role}")
    tree_manifest = _resolved_regular_file(
        Path(str(spec["tree_manifest_path"])), f"resource role {role} tree manifest"
    )
    return {
        "host_path": str(root),
        "tree_manifest_path": str(tree_manifest),
        "tree_manifest_bytes": tree_manifest.stat().st_size,
        "tree_manifest_sha256": sha256_file(tree_manifest),
    }


def _profile_resources(
    config: Mapping[str, Any], profile_id: str, dataset_id: str
) -> dict[str, Any]:
    if config.get("schema_version") != 1:
        _fail("resource-binding configuration schema is unsupported")
    profiles = config.get("profiles", {})
    dataset_profiles = config.get("dataset_profiles", {})
    if not isinstance(profiles, dict) or not isinstance(dataset_profiles, dict):
        _fail("resource-binding configuration maps are invalid")
    shared = profiles.get(profile_id, {})
    per_dataset = dataset_profiles.get(dataset_id, {})
    if not isinstance(shared, dict) or not isinstance(per_dataset, dict):
        _fail("resource-binding profile entry is invalid")
    overlay = per_dataset.get(profile_id, {})
    if not isinstance(overlay, dict):
        _fail("resource-binding dataset/profile entry is invalid")
    shared_roles = shared.get("roles", {})
    overlay_roles = overlay.get("roles", {})
    if not isinstance(shared_roles, dict) or not isinstance(overlay_roles, dict):
        _fail("resource-binding roles must be objects")
    return {**shared_roles, **overlay_roles}


def _builtins(
    manifest: Mapping[str, Any],
    dataset: Mapping[str, Any],
    profile: ValidatedCallerProfile,
    resource_manifest: Path,
    workspace: Path,
) -> dict[str, dict[str, Any]]:
    requested_roles = {mount["role"] for mount in profile.payload["mounts"]}
    gate_caller = profile.payload["gate_caller"]
    policy = manifest.get("caller_alignment_policy")
    if not isinstance(policy, dict) or policy.get(gate_caller) not in {
        "original_markdup_bam",
        "bqsr_cram",
    }:
        _fail("caller has no accepted alignment policy")
    alignment = dataset.get(policy[gate_caller])
    if not isinstance(alignment, dict):
        _fail("selected alignment role is absent from the dataset")
    reference = manifest.get("reference")
    if not isinstance(reference, dict):
        _fail("reference binding is absent from the input manifest")
    result: dict[str, dict[str, Any]] = {}
    if "selected_alignment" in requested_roles:
        result["selected_alignment"] = _manifest_file_binding(alignment)
    if "selected_alignment_index" in requested_roles:
        result["selected_alignment_index"] = _manifest_index_binding(alignment)
    for role in ("reference_fasta", "reference_fai", "reference_dictionary"):
        if role in requested_roles:
            result[role] = _reference_binding(reference, role)
    if "resource_manifest" in requested_roles:
        resource_manifest = _resolved_regular_file(resource_manifest, "resource manifest")
        result["resource_manifest"] = {
            "host_path": str(resource_manifest),
            "bytes": resource_manifest.stat().st_size,
            "sha256": sha256_file(resource_manifest),
        }
    if "workspace" in requested_roles:
        if workspace.is_symlink():
            _fail("workspace must not be a symbolic link")
        workspace.mkdir(parents=True, exist_ok=True)
        workspace = _resolved_regular_directory(workspace, "workspace")
        result["workspace"] = {"host_path": str(workspace)}
    baf = dataset.get("baf_vcf")
    if isinstance(baf, dict):
        if "source_baf_vcf" in requested_roles:
            result["source_baf_vcf"] = _manifest_file_binding(baf)
        if "source_baf_vcf_index" in requested_roles:
            result["source_baf_vcf_index"] = _manifest_index_binding(baf)
    return result


def _check_profile_constraint(
    role: str, mount: Mapping[str, Any], binding: Mapping[str, Any]
) -> None:
    expected_file = mount.get("expected_sha256")
    expected_tree = mount.get("expected_tree_manifest_sha256")
    if expected_file is not None and binding.get("sha256") != expected_file:
        _fail(f"resource file role {role} differs from its profile constraint")
    if expected_tree is not None and binding.get("tree_manifest_sha256") != expected_tree:
        _fail(f"resource tree role {role} differs from its profile constraint")


def _private_tokens(
    manifest: Mapping[str, Any],
    resource_config: Mapping[str, Any],
    paths: Sequence[Path],
) -> tuple[str, ...]:
    tokens = {str(path) for path in paths}
    for dataset in manifest.get("datasets", []):
        if isinstance(dataset, dict):
            for key in ("dataset_id", "logical_subject", "laboratory"):
                value = dataset.get(key)
                if isinstance(value, str) and value:
                    tokens.add(value)

    def visit(value: Any, key: str | None = None) -> None:
        if isinstance(value, dict):
            for nested_key, nested in value.items():
                visit(nested, str(nested_key))
        elif isinstance(value, list):
            for nested in value:
                visit(nested, key)
        elif key in {"path", "tree_manifest_path"} and isinstance(value, str):
            tokens.add(value)

    visit(resource_config)
    return tuple(sorted(token for token in tokens if token))


def _compile_caller_bindings(
    manifest_path: str | Path,
    dataset_id: str,
    profile_path: str | Path,
    resource_binding_config_path: str | Path,
    resource_manifest_path: str | Path,
    workspace_path: str | Path,
    output_path: str | Path,
    *,
    runtime: str = "docker",
    user: str,
    minimum_free_bytes: int,
) -> tuple[dict[str, Any], dict[str, Any], tuple[str, ...]]:
    """Resolve and verify one binding spec without publishing it."""

    manifest_file = Path(manifest_path)
    resource_config_file = Path(resource_binding_config_path)
    resource_manifest = Path(resource_manifest_path)
    workspace = Path(workspace_path)
    output = Path(output_path)
    manifest = _load_object(manifest_file, "private input manifest")
    config = _load_object(resource_config_file, "resource-binding configuration")
    profile = load_caller_profile(profile_path)
    datasets = manifest.get("datasets")
    if not isinstance(datasets, list):
        _fail("input manifest datasets are invalid")
    matches = [
        row for row in datasets if isinstance(row, dict) and row.get("dataset_id") == dataset_id
    ]
    if len(matches) != 1:
        _fail("dataset id must select exactly one input row")
    dataset = matches[0]
    resources = _profile_resources(config, profile.payload["profile_id"], dataset_id)
    builtins = _builtins(manifest, dataset, profile, resource_manifest, workspace)
    mounts: dict[str, dict[str, Any]] = {}
    requested_roles = {mount["role"] for mount in profile.payload["mounts"]}
    static_roles = requested_roles - set(builtins)
    if set(resources) != static_roles:
        _fail("resource-binding roles differ from the profile's unresolved mount roles")
    for mount in profile.payload["mounts"]:
        role = mount["role"]
        if role in builtins:
            binding = builtins[role]
        elif mount["kind"] == "file":
            binding = _resource_file_binding(resources[role], role)
        else:
            binding = _resource_tree_binding(resources[role], role)
        _check_profile_constraint(role, mount, binding)
        mounts[role] = binding
    if set(mounts) != requested_roles or BUILTIN_FILE_ROLES & set(resources):
        _fail("caller binding role reconciliation failed")
    binding_spec = {
        "schema_version": BINDING_SCHEMA_VERSION,
        "runtime": runtime,
        "user": user,
        "minimum_free_bytes": minimum_free_bytes,
        "mounts": mounts,
    }
    summary: dict[str, Any] = {
        "schema_version": BINDING_SCHEMA_VERSION,
        "status": "PASS",
        "profile_id": profile.payload["profile_id"],
        "profile_identity_sha256": profile.identity_sha256,
        "gate_caller": profile.payload["gate_caller"],
        "mount_roles": sorted(mounts),
        "mount_count": len(mounts),
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "binding_spec_scope": "private_runtime",
        "private_paths_persisted_in_binding_spec": True,
    }
    private = _private_tokens(
        manifest,
        config,
        (
            manifest_file,
            resource_config_file,
            resource_manifest,
            workspace,
            output,
        ),
    )
    rendered = json.dumps(summary, sort_keys=True)
    if any(contains_forbidden_binding(rendered, token) for token in private):
        _fail("portable binding summary contains a private binding")
    return binding_spec, summary, private


def build_caller_bindings(
    manifest_path: str | Path,
    dataset_id: str,
    profile_path: str | Path,
    resource_binding_config_path: str | Path,
    resource_manifest_path: str | Path,
    workspace_path: str | Path,
    output_path: str | Path,
    *,
    runtime: str = "docker",
    user: str,
    minimum_free_bytes: int,
) -> dict[str, Any]:
    """Create one exact private binding spec and return a portable summary."""

    binding_spec, summary, _private = _compile_caller_bindings(
        manifest_path,
        dataset_id,
        profile_path,
        resource_binding_config_path,
        resource_manifest_path,
        workspace_path,
        output_path,
        runtime=runtime,
        user=user,
        minimum_free_bytes=minimum_free_bytes,
    )
    _atomic_json_no_replace(Path(output_path), binding_spec)
    return summary


def validate_caller_bindings(
    manifest_path: str | Path,
    dataset_id: str,
    profile_path: str | Path,
    resource_binding_config_path: str | Path,
    resource_manifest_path: str | Path,
    workspace_path: str | Path,
    binding_spec_path: str | Path,
    *,
    runtime: str = "docker",
    user: str,
    minimum_free_bytes: int,
) -> dict[str, Any]:
    """Recompute and validate an immutable existing binding checkpoint."""

    output = Path(binding_spec_path)
    observed = _load_object(output, "caller binding specification")
    expected, summary, _private = _compile_caller_bindings(
        manifest_path,
        dataset_id,
        profile_path,
        resource_binding_config_path,
        resource_manifest_path,
        workspace_path,
        output,
        runtime=runtime,
        user=user,
        minimum_free_bytes=minimum_free_bytes,
    )
    if observed != expected:
        _fail("existing caller binding specification differs from recomputed bindings")
    return {**summary, "existing_binding_revalidated": True}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--dataset-id", required=True)
    profile_source = parser.add_mutually_exclusive_group(required=True)
    profile_source.add_argument("--profile-id")
    profile_source.add_argument("--profile", type=Path)
    parser.add_argument("--resource-bindings", type=Path, required=True)
    parser.add_argument("--resource-manifest", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runtime", choices=("docker", "podman"), default="docker")
    parser.add_argument("--user", required=True)
    parser.add_argument("--minimum-free-bytes", type=int, required=True)
    parser.add_argument("--validate-existing", action="store_true")
    args = parser.parse_args()
    try:
        profile = (
            resolve_shipped_profile(args.profile_id)
            if args.profile_id is not None
            else args.profile
        )
        function = validate_caller_bindings if args.validate_existing else build_caller_bindings
        result = function(
            args.manifest,
            args.dataset_id,
            profile,
            args.resource_bindings,
            args.resource_manifest,
            args.workspace,
            args.output,
            runtime=args.runtime,
            user=args.user,
            minimum_free_bytes=args.minimum_free_bytes,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except Exception:
        print(
            "caller binding construction failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
