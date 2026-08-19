"""Build and validate the private, checksum-bound Release C input manifest.

The manifest intentionally contains local paths and subject identifiers. It is
runtime input, not a portable artifact, and must remain outside version control.
Portable input-gate reports contain only cryptographic and structural evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

from .portable_model import PortableEvidenceError

SHA256 = re.compile(r"^[0-9a-f]{64}$")
ALIGNMENT_KEYS = ("original_markdup_bam", "bqsr_cram")
ACCEPTED_TECHNICAL_ROLES = frozenset({"primary", "technical_replicate"})
MANIFEST_STATES = frozenset({"PASS", "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH"})


class ManifestError(PortableEvidenceError):
    """Raised when a private Release C manifest cannot be proven safe to use."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise ManifestError(f"{label} is not a lower-case SHA-256")
    return value


def _string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ManifestError(f"{label} must be a non-empty string")
    return value.strip()


def _load_mapping(path: Path) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ManifestError("private Release C config is unreadable") from exc
    if not isinstance(value, dict):
        raise ManifestError("private Release C config must be a mapping")
    return value


def _bound_regular_file(raw_path: Any, *, base: Path, label: str) -> Path:
    value = _string(raw_path, f"{label} path")
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    probe = path.absolute()
    while True:
        if probe.is_symlink():
            raise ManifestError(f"{label} has a symbolic path component")
        if probe == probe.parent:
            break
        probe = probe.parent
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ManifestError(f"{label} must be a regular non-symbolic file")
    return path


def _index_spec(
    raw: dict[str, Any],
    *,
    base: Path,
    label: str,
) -> tuple[Path, dict[str, Any]]:
    path = _bound_regular_file(raw.get("path"), base=base, label=label)
    index = _bound_regular_file(raw.get("index_path"), base=base, label=f"{label} index")
    index_digest = sha256_file(index)
    supplied_index_digest = raw.get("index_sha256")
    if (
        supplied_index_digest is not None
        and _sha(supplied_index_digest, f"{label} supplied index digest") != index_digest
    ):
        raise ManifestError(f"{label} supplied index digest does not match the file")
    supplied_index_bytes = raw.get("index_bytes")
    if supplied_index_bytes is not None and supplied_index_bytes != index.stat().st_size:
        raise ManifestError(f"{label} supplied index byte count does not match the file")
    return path, {
        "path": str(path),
        "bytes": path.stat().st_size,
        "index_path": str(index),
        "index_bytes": index.stat().st_size,
        "index_sha256": index_digest,
    }


def _alignment_spec(
    raw: Any,
    *,
    base: Path,
    label: str,
    reuse_supplied_digest: bool,
) -> tuple[dict[str, Any], bool]:
    if not isinstance(raw, dict):
        raise ManifestError(f"{label} must be a mapping")
    path, result = _index_spec(raw, base=base, label=label)
    if reuse_supplied_digest:
        result["sha256"] = _sha(raw.get("sha256"), f"{label} supplied digest")
        result["sha256_origin"] = "supplied_upstream_provenance_pending_precall_rehash"
        deferred = True
    else:
        result["sha256"] = sha256_file(path)
        result["sha256_origin"] = "computed_by_manifest_builder"
        deferred = False
        supplied_digest = raw.get("sha256")
        if (
            supplied_digest is not None
            and _sha(supplied_digest, f"{label} supplied digest") != result["sha256"]
        ):
            raise ManifestError(f"{label} supplied digest does not match the file")
    supplied_bytes = raw.get("bytes")
    if supplied_bytes is not None and supplied_bytes != result["bytes"]:
        raise ManifestError(f"{label} supplied byte count does not match the file")
    return result, deferred


def _small_indexed_spec(
    raw: Any,
    *,
    base: Path,
    label: str,
) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise ManifestError(f"{label} must be a mapping")
    path, result = _index_spec(raw, base=base, label=label)
    result["sha256"] = sha256_file(path)
    supplied_digest = raw.get("sha256")
    if (
        supplied_digest is not None
        and _sha(supplied_digest, f"{label} supplied digest") != result["sha256"]
    ):
        raise ManifestError(f"{label} supplied digest does not match the file")
    supplied_bytes = raw.get("bytes")
    if supplied_bytes is not None and supplied_bytes != result["bytes"]:
        raise ManifestError(f"{label} supplied byte count does not match the file")
    return result


def validate_manifest_structure(value: Any) -> dict[str, Any]:
    """Validate schema and internal identity rules without opening bound files."""

    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise ManifestError("Release C manifest schema is invalid")
    if value.get("status") not in MANIFEST_STATES:
        raise ManifestError("Release C manifest status is invalid")
    reference = value.get("reference")
    if not isinstance(reference, dict) or not str(reference.get("assembly", "")).startswith(
        "GRCh38"
    ):
        raise ManifestError("Release C manifest reference is not GRCh38")
    for prefix in ("fasta", "fai", "dictionary"):
        _string(reference.get(f"{prefix}_path"), f"reference {prefix} path")
        _sha(reference.get(f"{prefix}_sha256"), f"reference {prefix} digest")
        size = reference.get(f"{prefix}_bytes")
        if not isinstance(size, int) or isinstance(size, bool) or size < 1:
            raise ManifestError(f"reference {prefix} byte count is invalid")
    policy = value.get("caller_alignment_policy")
    if not isinstance(policy, dict) or not policy:
        raise ManifestError("Release C manifest caller policy is empty")
    for caller, alignment_key in policy.items():
        _string(caller, "caller policy name")
        if alignment_key not in ALIGNMENT_KEYS:
            raise ManifestError("Release C caller policy selects an unsupported alignment role")
    datasets = value.get("datasets")
    if (
        not isinstance(datasets, list)
        or not datasets
        or any(not isinstance(row, dict) for row in datasets)
    ):
        raise ManifestError("Release C manifest datasets are invalid")
    dataset_ids = [_string(row.get("dataset_id"), "dataset id") for row in datasets]
    if len(dataset_ids) != len(set(dataset_ids)):
        raise ManifestError("Release C manifest dataset ids are duplicated")
    for row in datasets:
        _string(row.get("logical_subject"), "logical subject")
        if row.get("technical_role") not in ACCEPTED_TECHNICAL_ROLES:
            raise ManifestError("Release C technical role is invalid")
        present_alignments = [key for key in ALIGNMENT_KEYS if isinstance(row.get(key), dict)]
        if not present_alignments:
            raise ManifestError("Release C dataset has no alignment")
        for key in present_alignments:
            spec = row[key]
            for path_key in ("path", "index_path"):
                _string(spec.get(path_key), f"{key} {path_key}")
            for digest_key in ("sha256", "index_sha256"):
                _sha(spec.get(digest_key), f"{key} {digest_key}")
            for size_key in ("bytes", "index_bytes"):
                size = spec.get(size_key)
                if not isinstance(size, int) or isinstance(size, bool) or size < 1:
                    raise ManifestError(f"{key} {size_key} is invalid")
            if spec.get("sha256_origin") not in {
                "computed_by_manifest_builder",
                "supplied_upstream_provenance_pending_precall_rehash",
            }:
                raise ManifestError(f"{key} digest origin is invalid")
        if "baf_vcf" in row:
            baf = row["baf_vcf"]
            if not isinstance(baf, dict):
                raise ManifestError("BAF VCF binding is invalid")
            _string(baf.get("sample"), "BAF VCF sample")
            for path_key in ("path", "index_path"):
                _string(baf.get(path_key), f"BAF VCF {path_key}")
            for digest_key in ("sha256", "index_sha256"):
                _sha(baf.get(digest_key), f"BAF VCF {digest_key}")
            for size_key in ("bytes", "index_bytes"):
                size = baf.get(size_key)
                if not isinstance(size, int) or isinstance(size, bool) or size < 1:
                    raise ManifestError(f"BAF VCF {size_key} is invalid")
    for caller, alignment_key in policy.items():
        if any(alignment_key not in row for row in datasets):
            raise ManifestError(f"caller {caller} selects an alignment absent from a dataset")
        if caller == "cnvpytor" and any("baf_vcf" not in row for row in datasets):
            raise ManifestError("CNVpytor policy requires a BAF VCF for every dataset")
    deferred = any(
        row[key].get("sha256_origin") == "supplied_upstream_provenance_pending_precall_rehash"
        for row in datasets
        for key in ALIGNMENT_KEYS
        if key in row
    )
    expected_state = (
        "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH" if deferred else "PASS"
    )
    if value.get("status") != expected_state:
        raise ManifestError("Release C manifest status disagrees with alignment digest origins")
    return value


def build_private_manifest(
    config_path: str | Path,
    output_path: str | Path,
    *,
    reuse_supplied_alignment_sha256: bool = False,
) -> dict[str, Any]:
    """Create a no-overwrite private manifest from a generic YAML/JSON config."""

    config_path = Path(config_path).resolve(strict=True)
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError("refusing to replace an existing private manifest")
    config = _load_mapping(config_path)
    if config.get("schema_version") != 1:
        raise ManifestError("private Release C config schema is invalid")
    base = config_path.parent
    raw_reference = config.get("reference")
    if not isinstance(raw_reference, dict):
        raise ManifestError("private Release C config lacks reference bindings")
    assembly = _string(raw_reference.get("assembly"), "reference assembly")
    if not assembly.startswith("GRCh38"):
        raise ManifestError("Release C requires an explicitly declared GRCh38 reference")
    reference_paths = {
        key: _bound_regular_file(
            raw_reference.get(f"{key}_path"), base=base, label=f"reference {key}"
        )
        for key in ("fasta", "fai", "dictionary")
    }
    reference: dict[str, Any] = {"assembly": assembly}
    for key, path in reference_paths.items():
        reference[f"{key}_path"] = str(path)
        reference[f"{key}_bytes"] = path.stat().st_size
        reference[f"{key}_sha256"] = sha256_file(path)
        supplied_digest = raw_reference.get(f"{key}_sha256")
        if (
            supplied_digest is not None
            and _sha(supplied_digest, f"reference {key} supplied digest")
            != reference[f"{key}_sha256"]
        ):
            raise ManifestError(f"reference {key} supplied digest does not match the file")

    raw_policy = config.get("caller_alignment_policy")
    if not isinstance(raw_policy, dict):
        raise ManifestError("private Release C config lacks caller alignment policy")
    policy = {
        _string(caller, "caller name"): _string(role, "alignment role")
        for caller, role in raw_policy.items()
    }
    raw_datasets = config.get("datasets")
    if not isinstance(raw_datasets, list) or not raw_datasets:
        raise ManifestError("private Release C config has no datasets")
    datasets: list[dict[str, Any]] = []
    deferred = False
    for index, raw_row in enumerate(raw_datasets):
        if not isinstance(raw_row, dict):
            raise ManifestError(f"dataset row {index + 1} is not a mapping")
        row: dict[str, Any] = {
            "dataset_id": _string(raw_row.get("dataset_id"), "dataset id"),
            "logical_subject": _string(raw_row.get("logical_subject"), "logical subject"),
            "technical_role": _string(raw_row.get("technical_role"), "technical role"),
        }
        laboratory = raw_row.get("laboratory")
        if laboratory is not None:
            row["laboratory"] = _string(laboratory, "laboratory")
        for alignment_key in ALIGNMENT_KEYS:
            if alignment_key not in raw_row:
                continue
            alignment, was_deferred = _alignment_spec(
                raw_row[alignment_key],
                base=base,
                label=f"{row['dataset_id']} {alignment_key}",
                reuse_supplied_digest=reuse_supplied_alignment_sha256,
            )
            row[alignment_key] = alignment
            deferred = deferred or was_deferred
        if "baf_vcf" in raw_row:
            baf = _small_indexed_spec(
                raw_row["baf_vcf"],
                base=base,
                label=f"{row['dataset_id']} BAF VCF",
            )
            baf["sample"] = _string(raw_row["baf_vcf"].get("sample"), "BAF VCF sample")
            row["baf_vcf"] = baf
        datasets.append(row)

    state = "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH" if deferred else "PASS"
    manifest = {
        "schema_version": 1,
        "status": state,
        "generated_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "privacy": "PRIVATE_RUNTIME_MANIFEST_CONTAINS_IDENTIFIERS_AND_LOCAL_PATHS",
        "large_alignment_hash_policy": (
            "Every selected alignment is fully rehashed by its caller-specific input gate "
            "immediately before execution."
        ),
        "reference": reference,
        "caller_alignment_policy": policy,
        "datasets": datasets,
    }
    validate_manifest_structure(manifest)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.partial.", dir=output_path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, output_path)
        except FileExistsError as exc:
            raise FileExistsError("refusing to replace an existing private manifest") from exc
        temporary.unlink()
        directory_descriptor = os.open(output_path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return {
        "schema_version": 1,
        "status": state,
        "datasets": len(datasets),
        "logical_subjects": len({row["logical_subject"] for row in datasets}),
        "alignment_hashes_deferred_to_precall_gate": deferred,
        "manifest_sha256": hashlib.sha256(payload).hexdigest(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create")
    create.add_argument("--config", type=Path, required=True)
    create.add_argument("--output", type=Path, required=True)
    create.add_argument("--reuse-supplied-alignment-sha256", action="store_true")
    validate = subparsers.add_parser("validate-structure")
    validate.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "create":
            result = build_private_manifest(
                args.config,
                args.output,
                reuse_supplied_alignment_sha256=args.reuse_supplied_alignment_sha256,
            )
        else:
            document = json.loads(args.manifest.read_text(encoding="utf-8"))
            validate_manifest_structure(document)
            result = {
                "schema_version": 1,
                "status": "PASS",
                "datasets": len(document["datasets"]),
            }
    except (ManifestError, FileExistsError, OSError, UnicodeError, json.JSONDecodeError) as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, sort_keys=True))
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
