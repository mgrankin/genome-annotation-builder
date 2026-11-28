"""Create and validate deterministic, path-free Release-C resource manifests."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

from .portable_model import (
    FORBIDDEN_PORTABLE_TEXT,
    PortableEvidenceError,
    contains_forbidden_binding,
)
from .profile_catalog import CURRENT_CATALOG_VERSION, build_profile_catalog

RESOURCE_STATUS = {"PASS", "PASS_WITH_DOCUMENTED_LIMITATIONS"}
METADATA_FIELDS = {
    "schema_version",
    "status",
    "scope",
    "production_dataset_execution_claimed",
    "tools",
    "resources",
    "runtime_resource_trees",
    "validation_evidence",
    "warnings",
}
MANIFEST_FIELDS = METADATA_FIELDS | {
    "profile_catalog",
    "private_paths_persisted",
    "manifest_identity_sha256",
}


def _canonical_sha256(value: object) -> str:
    import hashlib

    rendered = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode()
    return hashlib.sha256(rendered).hexdigest()


def _load_object(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise PortableEvidenceError(f"{label} is missing or symlinked")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError(f"{label} is unreadable") from exc
    if not isinstance(value, dict):
        raise PortableEvidenceError(f"{label} must be a JSON object")
    return value


def _forbidden_tokens(path: Path | None) -> tuple[str, ...]:
    if path is None:
        return ()
    if path.is_symlink() or not path.is_file():
        raise PortableEvidenceError("forbidden-binding file is missing or symlinked")
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    return tuple(sorted({value for value in values if value}))


def _assert_portable(value: object, forbidden: tuple[str, ...]) -> None:
    rendered = json.dumps(value, sort_keys=True, ensure_ascii=True, allow_nan=False)
    if FORBIDDEN_PORTABLE_TEXT.search(rendered) or any(
        contains_forbidden_binding(rendered, token) for token in forbidden
    ):
        raise PortableEvidenceError("resource manifest contains a private binding")


def _validate_metadata(metadata: dict[str, Any], forbidden: tuple[str, ...]) -> None:
    if set(metadata) != METADATA_FIELDS or metadata.get("schema_version") != 1:
        raise PortableEvidenceError("resource metadata field set or schema is invalid")
    status = metadata.get("status")
    warnings = metadata.get("warnings")
    if (
        status not in RESOURCE_STATUS
        or not isinstance(metadata.get("scope"), str)
        or not metadata["scope"].strip()
        or metadata.get("production_dataset_execution_claimed") is not False
        or any(
            not isinstance(metadata.get(field), dict) or not metadata[field]
            for field in ("tools", "resources", "validation_evidence")
        )
        or not isinstance(metadata.get("runtime_resource_trees"), dict)
        or not isinstance(warnings, list)
        or any(not isinstance(value, str) or not value.strip() for value in warnings)
        or len(warnings) != len(set(warnings))
        or (status == "PASS_WITH_DOCUMENTED_LIMITATIONS" and not warnings)
        or (status == "PASS" and warnings)
    ):
        raise PortableEvidenceError("resource metadata readiness boundary is invalid")
    _assert_portable(metadata, forbidden)


def build_resource_manifest(
    metadata: dict[str, Any],
    *,
    forbidden_bindings: tuple[str, ...] = (),
    catalog_version: int = CURRENT_CATALOG_VERSION,
) -> dict[str, Any]:
    """Bind deployment metadata to the exact caller-profile release."""

    _validate_metadata(metadata, forbidden_bindings)
    manifest = {
        **metadata,
        "profile_catalog": build_profile_catalog(
            forbidden_bindings=forbidden_bindings, catalog_version=catalog_version
        ),
        "private_paths_persisted": False,
    }
    manifest["manifest_identity_sha256"] = _canonical_sha256(manifest)
    _assert_portable(manifest, forbidden_bindings)
    return manifest


def validate_resource_manifest_payload(
    manifest: dict[str, Any], *, forbidden_bindings: tuple[str, ...] = ()
) -> dict[str, Any]:
    """Validate exact field, profile-catalog, privacy, and canonical identities."""

    if set(manifest) != MANIFEST_FIELDS:
        raise PortableEvidenceError("resource manifest field set is invalid")
    metadata = {field: manifest[field] for field in METADATA_FIELDS}
    _validate_metadata(metadata, forbidden_bindings)
    if manifest.get("private_paths_persisted") is not False:
        raise PortableEvidenceError("resource manifest privacy boundary is invalid")
    catalog = manifest.get("profile_catalog")
    if not isinstance(catalog, dict) or catalog != build_profile_catalog(
        forbidden_bindings=forbidden_bindings, catalog_version=catalog.get("schema_version")
    ):
        raise PortableEvidenceError("resource manifest profile catalog changed")
    claimed = manifest.get("manifest_identity_sha256")
    unsigned = dict(manifest)
    unsigned.pop("manifest_identity_sha256", None)
    if claimed != _canonical_sha256(unsigned):
        raise PortableEvidenceError("resource manifest identity changed")
    _assert_portable(manifest, forbidden_bindings)
    return {
        "schema_version": 1,
        "status": manifest["status"],
        "profile_count": catalog["profile_count"],
        "tool_records": len(manifest["tools"]),
        "resource_records": len(manifest["resources"]),
        "runtime_resource_tree_records": len(manifest["runtime_resource_trees"]),
        "validation_evidence_records": len(manifest["validation_evidence"]),
        "warnings": len(manifest["warnings"]),
        "manifest_identity_sha256": claimed,
        "private_paths_persisted": False,
    }


def _atomic_no_replace(path: Path, value: object) -> None:
    if path.exists():
        raise PortableEvidenceError("refusing to replace an existing resource manifest")
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    descriptor, raw_name = tempfile.mkstemp(prefix=f".{path.name}.partial.", dir=path.parent)
    temporary = Path(raw_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as exc:
            raise PortableEvidenceError(
                "refusing to replace an existing resource manifest"
            ) from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="operation", required=True)
    create = subparsers.add_parser("create")
    create.add_argument(
        "--catalog-version",
        type=int,
        choices=(1, CURRENT_CATALOG_VERSION),
        default=CURRENT_CATALOG_VERSION,
    )
    create.add_argument("--metadata", type=Path, required=True)
    create.add_argument("--output", type=Path, required=True)
    create.add_argument("--forbidden-bindings-file", type=Path)
    validate = subparsers.add_parser("validate")
    validate.add_argument("--manifest", type=Path, required=True)
    validate.add_argument("--forbidden-bindings-file", type=Path)
    args = parser.parse_args(argv)
    try:
        forbidden = _forbidden_tokens(args.forbidden_bindings_file)
        if args.operation == "create":
            manifest = build_resource_manifest(
                _load_object(args.metadata, "resource metadata"),
                catalog_version=args.catalog_version,
                forbidden_bindings=forbidden,
            )
            _atomic_no_replace(args.output, manifest)
            result = validate_resource_manifest_payload(
                _load_object(args.output, "resource manifest"),
                forbidden_bindings=forbidden,
            )
            result["resume"] = "created"
        else:
            result = validate_resource_manifest_payload(
                _load_object(args.manifest, "resource manifest"),
                forbidden_bindings=forbidden,
            )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, PortableEvidenceError):
        print("error: release_c_resource_manifest_failed", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
