"""Atomically enrich an immutable .genome/1.1 bundle with Release-C evidence."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from ..constants import PIPELINE_VERSION, SCHEMA_VERSION
from ..query import query_bundle
from ..util import atomic_output_directory, atomic_write_json, directory_digest, sha256_file
from ..validator import validate_bundle
from .extension_stage import MANIFEST_FILENAME, validate_release_c_extension_stage
from .portable_model import DOMAINS, PortableEvidenceError, schema_for
from .readback_capture_v2 import EMPTY_TARGET_SHA256

RELEASE_C_SCHEMA_VERSION = 1
RELEASE_C_PROFILE_VERSION = "1.0.0"
REBUILT_PATHS = frozenset(
    {
        "manifest.json",
        "inventory.json",
        "raw_artifacts.json",
        "transformations.json",
        "omissions.json",
        "relationships.json",
        "field_glossary.json",
        "readmygenome.md",
        "readmygenome/SKILL.md",
        "readmygenome/reference/query_catalog.md",
        "readmygenome/scripts/query.py",
    }
)
DOMAIN_ROLES = {
    "structural_variants": "source-specific read-backed structural-variant observations",
    "copy_number": "source-specific copy-number and allelic-imbalance observations",
    "repeat_expansions": "source-specific short-read repeat-expansion observations",
    "hla_typing": "source-specific HLA typing observations and ambiguity sets",
    "complex_pharmacogenes": "source-specific complex pharmacogene observations",
    "mitochondrial_heteroplasmy": "source-specific mitochondrial small-variant observations",
    "read_backed_evidence": "bounded local fragment evidence for selected candidate alleles",
}
DOMAIN_OMISSIONS = {
    "structural_variants": "read-backed structural-variant calling",
    "copy_number": "read-backed copy-number and loss-of-heterozygosity calling",
    "repeat_expansions": "read-backed repeat-expansion calling",
    "hla_typing": "read-backed HLA typing",
    "complex_pharmacogenes": "licensed complex-pharmacogene calling",
    "mitochondrial_heteroplasmy": "mitochondrial heteroplasmy calling",
    "read_backed_evidence": "bounded local read-back evidence",
}
COMMON_ALLOWED_VALUES: dict[str, list[Any]] = {
    "assembly": ["GRCh38"],
    "technical_source_role": ["primary", "technical_replicate", "combined_evidence"],
    "confidence_tier": ["SUPPORTED", "REVIEW", "UNRESOLVED"],
    "concordance_status": [
        "not_assessed",
        "comparable_concordant",
        "comparable_discordant",
        "not_comparable",
        "source_specific",
    ],
}
DOMAIN_ALLOWED_VALUES: dict[tuple[str, str], list[Any]] = {
    ("structural_variants", "sv_type"): ["DEL", "DUP", "INV", "INS", "BND", "CNV"],
    ("structural_variants", "call_status"): ["non_reference", "reference", "no_call"],
    ("repeat_expansions", "variant_subtype"): [
        "Repeat",
        "RareRepeat",
        "SmallVariant",
        "SMN",
    ],
    ("hla_typing", "call_status"): ["complete", "partial", "no_call"],
    ("mitochondrial_heteroplasmy", "chrom"): ["chrM"],
    ("read_backed_evidence", "variant_type"): [
        "SNV",
        "MNV",
        "insertion",
        "deletion",
        "complex",
    ],
    ("read_backed_evidence", "assessment_status"): [
        "assessed",
        "no_eligible_reads",
        "unsupported_allele_representation",
    ],
    ("read_backed_evidence", "overlapping_mate_policy"): ["collapse_by_query_name"],
}
GENERIC_READ_BACKED_OMISSION = (
    "read-backed SV/CNV, repeat, complex-PGx, mosaicism, and mitochondrial analysis"
)
PORTABLE_CORE_SCOPE_OMISSION = "caller records outside the portable GRCh38 core-contig scope"

README_APPENDIX = """

## Read-backed extensions

If `extensions/release_c/` is present, its typed tables retain one observation
per caller and technical source. Query them with:

```bash
python readmygenome/scripts/query.py read-backed <locus-or-term> --bundle .
```

Concordance annotations never replace or average source rows. Absence of a row
is not a reference, normal, benign, or unaffected result.
"""

SKILL_APPENDIX = """

## Release-C read-backed evidence

- Use `read-backed` for SV, CNV, repeats, HLA, complex pharmacogenes,
  mitochondrial calls, and bounded local read evidence.
- Report `caller_name`, `technical_source_role`, `confidence_tier`,
  `concordance_status`, `filter_reasons`, `limitations`, and
  `confirmation_required` with any finding.
- Never collapse caller or technical-replicate rows into a synthetic consensus.
- Missing rows have no negative meaning; consult `omissions.json` and the
  caller receipts in `extensions/release_c/release_c_extensions.json`.
"""

QUERY_CATALOG_APPENDIX = """

## Read-backed extensions

```bash
python readmygenome/scripts/query.py read-backed chr7:117559593 --bundle .
python readmygenome/scripts/query.py read-backed HLA-B --bundle .
python readmygenome/scripts/query.py read-backed CYP2D6 --bundle .
python readmygenome/scripts/query.py read-backed CAG --bundle .
```

The command is output-bounded and preserves every source observation. It does
not perform diagnosis or infer a negative result from an absent row.
"""


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError(f"{label} is unreadable") from exc
    if not isinstance(value, dict):
        raise PortableEvidenceError(f"{label} must be a JSON object")
    return value


def _canonical_timestamp(value: str | None) -> str:
    if value is None:
        parsed = datetime.now(UTC)
    else:
        try:
            parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        except ValueError as exc:
            raise PortableEvidenceError("generated-at is not ISO-8601") from exc
        if parsed.tzinfo is None:
            raise PortableEvidenceError("generated-at must include a timezone")
        parsed = parsed.astimezone(UTC)
    return parsed.isoformat().replace("+00:00", "Z")


def _regular_tree_files(root: Path) -> list[Path]:
    if root.is_symlink() or not root.is_dir():
        raise PortableEvidenceError(f"not a regular directory: {root}")
    files: list[Path] = []
    for item in sorted(root.rglob("*")):
        if item.is_symlink():
            raise PortableEvidenceError(f"symbolic link is forbidden: {item}")
        if item.is_file():
            files.append(item)
        elif not item.is_dir():
            raise PortableEvidenceError(f"special filesystem entry is forbidden: {item}")
    return files


def _materialize_base(
    base: Path,
    destination: Path,
    *,
    mode: str,
) -> dict[str, int | str]:
    if mode not in {"copy", "hardlink"}:
        raise PortableEvidenceError("base materialization must be copy or hardlink")
    linked = 0
    copied = 0
    bytes_reused = 0
    for directory in sorted(item for item in base.rglob("*") if item.is_dir()):
        if directory.is_symlink():
            raise PortableEvidenceError("base bundle contains a symbolic link")
        (destination / directory.relative_to(base)).mkdir(parents=True, exist_ok=True)
    for source in _regular_tree_files(base):
        relative = source.relative_to(base).as_posix()
        if relative in REBUILT_PATHS:
            continue
        if relative == "extensions/release_c" or relative.startswith("extensions/release_c/"):
            raise PortableEvidenceError("base bundle is already Release-C enriched")
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            raise PortableEvidenceError("base materialization path collision")
        if mode == "hardlink":
            os.link(source, target)
            linked += 1
        else:
            shutil.copy2(source, target)
            copied += 1
        bytes_reused += source.stat().st_size
    return {
        "mode": mode,
        "hardlinked_files": linked,
        "copied_files": copied,
        "reused_bytes": bytes_reused,
    }


def _copy_regular_tree(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    for directory in sorted(item for item in source.rglob("*") if item.is_dir()):
        if directory.is_symlink():
            raise PortableEvidenceError("Release-C stage contains a symbolic link")
        (destination / directory.relative_to(source)).mkdir(parents=True, exist_ok=True)
    for item in _regular_tree_files(source):
        target = destination / item.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        with item.open("rb") as read_handle, target.open("xb") as write_handle:
            shutil.copyfileobj(read_handle, write_handle, length=16 * 1024 * 1024)
            write_handle.flush()
            os.fsync(write_handle.fileno())


def _append_document(base: Path, destination: Path, appendix: str, marker: str) -> None:
    text = base.read_text(encoding="utf-8")
    if marker in text:
        raise PortableEvidenceError("base Agent Profile is already Release-C enriched")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text.rstrip() + appendix.rstrip() + "\n", encoding="utf-8")


def _domain_file_metadata(stage_manifest: dict[str, Any], domain: str) -> dict[str, Any]:
    metadata = stage_manifest["domains"][domain]
    return {
        "bytes": metadata["bytes"],
        "rows": metadata["rows"],
        "sha256": metadata["sha256"],
    }


def _release_c_raw_rows(stage_manifest: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for wrapper in stage_manifest["caller_runs"]:
        receipt = wrapper["receipt"]
        boundary = receipt["run_boundary"]
        artifact = receipt["raw_artifact"]
        rows.append(
            {
                "artifact_type": (
                    "Release-C caller raw artifact: "
                    + str(boundary["caller_name"])
                    + " "
                    + str(boundary["caller_version"])
                ),
                "disposition": "separate",
                "integrity_identifier": (
                    "receipt-sha256:"
                    + str(wrapper["receipt_file_sha256"])
                    + ";receipt-identity-sha256:"
                    + str(receipt["receipt_identity_sha256"])
                    + ";raw-artifact-sha256:"
                    + str(artifact["raw_artifact_sha256"])
                ),
                "reason": (
                    "verbatim caller output remains outside the portable bundle; typed "
                    "source-specific evidence and the complete cryptographic receipt are embedded"
                ),
                "limitations": [
                    "verbatim caller-output recovery requires the separately retained raw artifact",
                    "typed evidence is not a clinical diagnosis and may require orthogonal confirmation",
                ],
                "caller_run_identity_sha256": boundary["caller_run_identity_sha256"],
                "resource_manifest_sha256": boundary["resource_manifest_sha256"],
                "reference_sha256": boundary["reference_sha256"],
                "source_alignment_sha256": boundary["source_alignment_sha256"],
                "raw_artifact_id": artifact["raw_artifact_id"],
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            str(row["artifact_type"]),
            str(row["caller_run_identity_sha256"]),
            str(row["source_alignment_sha256"]),
        ),
    )


def _augment_raw_artifacts(
    document: dict[str, Any], stage_manifest: dict[str, Any]
) -> dict[str, Any]:
    rows = document.get("artifacts")
    if document.get("schema_version") != SCHEMA_VERSION or not isinstance(rows, list):
        raise PortableEvidenceError("base raw_artifacts.json is invalid")
    if any(
        isinstance(row, dict)
        and str(row.get("artifact_type", "")).startswith("Release-C caller raw artifact:")
        for row in rows
    ):
        raise PortableEvidenceError("base raw-artifact provenance already contains Release C")
    return {**document, "artifacts": [*rows, *_release_c_raw_rows(stage_manifest)]}


def _augment_transformations(
    document: dict[str, Any],
    stage_manifest: dict[str, Any],
    generated_at: str,
    base_manifest_sha256: str,
    materialization: dict[str, int | str],
) -> dict[str, Any]:
    stages = document.get("stages")
    if document.get("schema_version") != SCHEMA_VERSION or not isinstance(stages, list):
        raise PortableEvidenceError("base transformations.json is invalid")
    if any(
        isinstance(row, dict) and row.get("stage_id") == "release_c_bundle_extension"
        for row in stages
    ):
        raise PortableEvidenceError("base transformation provenance already contains Release C")
    caller_receipts = sorted(
        wrapper["receipt"]["receipt_identity_sha256"] for wrapper in stage_manifest["caller_runs"]
    )
    stage = {
        "stage_id": "release_c_bundle_extension",
        "software": "genome-annotation-builder",
        "software_version": PIPELINE_VERSION,
        "configuration_digest": stage_manifest["extension_set_identity_sha256"],
        "inputs": [
            "sha256:" + base_manifest_sha256,
            *["sha256:" + digest for digest in caller_receipts],
        ],
        "outputs": [
            "extensions/release_c/" + stage_manifest["domains"][domain]["path"]
            for domain in sorted(stage_manifest["domains"])
        ],
        "completed_at": generated_at,
        "status": "completed",
        "extension_schema_version": RELEASE_C_SCHEMA_VERSION,
        "extension_set_identity_sha256": stage_manifest["extension_set_identity_sha256"],
        "base_bundle_manifest_sha256": base_manifest_sha256,
        "base_materialization": materialization,
        "source_rows_removed": 0,
        "caller_or_replicate_rows_collapsed": 0,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return {**document, "stages": [*stages, stage]}


def _omission_row(domain: str, stage_manifest: dict[str, Any]) -> dict[str, str]:
    scope = DOMAIN_OMISSIONS[domain]
    if domain == "read_backed_evidence" and any(
        row.get("domain") == domain
        and row["summary"].get("target_manifest_sha256") == EMPTY_TARGET_SHA256
        and row["summary"].get("capture_records") == 0
        for run in stage_manifest.get("adapter_runs", [])
        for row in run["adapter_summary"]["parsing"]
    ):
        reason = (
            "a completed bounded read-back capture received a verified target list with zero loci"
        )
        recovery = "select specific loci for an additional bounded read-back assessment when needed"
    elif domain == "complex_pharmacogenes":
        reason = "no licence-accepted, checksum-pinned complex-pharmacogene caller was configured"
        recovery = (
            "accept the applicable caller licence, configure an exact pinned profile, and rebuild"
        )
    else:
        reason = f"no validated {domain} adapter output was supplied to the Release-C stage"
        recovery = "configure a checksum-pinned read-backed caller profile and rebuild"
    return {
        "scope": scope,
        "reason": reason,
        "pipeline_stage": "release_c_read_backed_extensions",
        "interpretation_consequence": (
            "absence of this domain is not a normal, reference, benign, or unaffected result"
        ),
        "recovery": recovery,
    }


def _portable_core_scope_omission(stage_manifest: dict[str, Any]) -> dict[str, Any] | None:
    reconciliation = stage_manifest.get("adapter_reconciliation")
    scope = reconciliation.get("portable_core_scope") if isinstance(reconciliation, dict) else None
    if not isinstance(scope, dict):
        raise PortableEvidenceError("Release-C portable-scope reconciliation is absent")
    omitted = scope.get("omitted_records")
    if not isinstance(omitted, int) or omitted < 0:
        raise PortableEvidenceError("Release-C portable-scope omission count is invalid")
    if omitted == 0:
        return None
    return {
        "scope": PORTABLE_CORE_SCOPE_OMISSION,
        "reason": (
            f"{omitted} raw structural-variant caller records were explicitly omitted because "
            "at least one endpoint was outside chr1-22/X/Y/M; reciprocal BND mates were omitted "
            "as complete pairs"
        ),
        "pipeline_stage": "release_c_strict_adapter",
        "interpretation_consequence": (
            "these aggregate-accounted records are absent from portable typed evidence; their "
            "absence is not a reference, benign, or unaffected result"
        ),
        "recovery": (
            "inspect the separately retained immutable raw caller artifact and its embedded "
            "caller/adapter receipts"
        ),
        "portable_contig_scope": scope["policy"],
        "adapter_runs": scope["adapter_runs"],
        "source_records": scope["source_records"],
        "typed_rows": scope["typed_rows"],
        "omitted_records": omitted,
        "omitted_non_bnd_records": scope["omitted_non_bnd_records"],
        "omitted_bnd_records": scope["omitted_bnd_records"],
        "omitted_bnd_pairs": scope["omitted_bnd_pairs"],
        "core_to_noncanonical_bnd_pairs": scope["core_to_noncanonical_bnd_pairs"],
        "noncanonical_to_noncanonical_bnd_pairs": scope["noncanonical_to_noncanonical_bnd_pairs"],
    }


def _augment_omissions(document: dict[str, Any], stage_manifest: dict[str, Any]) -> dict[str, Any]:
    rows = document.get("omissions")
    if document.get("schema_version") != SCHEMA_VERSION or not isinstance(rows, list):
        raise PortableEvidenceError("base omissions.json is invalid")
    retained = [
        row
        for row in rows
        if not (isinstance(row, dict) and row.get("scope") == GENERIC_READ_BACKED_OMISSION)
    ]
    if len(rows) - len(retained) > 1:
        raise PortableEvidenceError("base bundle duplicates its generic read-backed omission")
    release_c_scopes = set(DOMAIN_OMISSIONS.values()) | {
        "genome-wide mosaic variant calling",
        PORTABLE_CORE_SCOPE_OMISSION,
    }
    if any(isinstance(row, dict) and row.get("scope") in release_c_scopes for row in retained):
        raise PortableEvidenceError("base omissions already contain Release-C domain rows")
    missing = sorted(set(DOMAINS) - set(stage_manifest["domains"]))
    retained.extend(_omission_row(domain, stage_manifest) for domain in missing)
    portable_scope_omission = _portable_core_scope_omission(stage_manifest)
    if portable_scope_omission is not None:
        retained.append(portable_scope_omission)
    retained.append(
        {
            "scope": "genome-wide mosaic variant calling",
            "reason": (
                "bounded local read-back evaluates selected candidates and is not a "
                "genome-wide mosaic discovery caller"
            ),
            "pipeline_stage": "release_c_read_backed_extensions",
            "interpretation_consequence": (
                "absence of a mosaic row is not evidence that mosaic variation is absent"
            ),
            "recovery": (
                "run a separately validated genome-wide mosaic caller with matched quality controls"
            ),
        }
    )
    return {**document, "omissions": retained}


def _augment_relationships(
    document: dict[str, Any], stage_manifest: dict[str, Any]
) -> dict[str, Any]:
    rows = document.get("relationships")
    if document.get("schema_version") != SCHEMA_VERSION or not isinstance(rows, list):
        raise PortableEvidenceError("base relationships.json is invalid")
    existing_ids = {row.get("relationship_id") for row in rows if isinstance(row, dict)}
    additions: list[dict[str, Any]] = []
    for domain in sorted(stage_manifest["domains"]):
        relationship_id = f"release_c_{domain}_to_caller_receipt"
        if relationship_id in existing_ids:
            raise PortableEvidenceError("base relationships already contain Release C")
        additions.append(
            {
                "relationship_id": relationship_id,
                "source": [
                    f"{domain}.caller_name",
                    f"{domain}.caller_version",
                    f"{domain}.caller_run_identity_sha256",
                    f"{domain}.resource_manifest_sha256",
                    f"{domain}.source_alignment_sha256",
                    f"{domain}.reference_sha256",
                    f"{domain}.raw_artifact_id",
                    f"{domain}.raw_artifact_sha256",
                ],
                "target": [
                    "release_c_extensions.caller_runs.receipt.run_boundary.caller_name",
                    "release_c_extensions.caller_runs.receipt.run_boundary.caller_version",
                    (
                        "release_c_extensions.caller_runs.receipt.run_boundary."
                        "caller_run_identity_sha256"
                    ),
                    (
                        "release_c_extensions.caller_runs.receipt.run_boundary."
                        "resource_manifest_sha256"
                    ),
                    (
                        "release_c_extensions.caller_runs.receipt.run_boundary."
                        "source_alignment_sha256"
                    ),
                    "release_c_extensions.caller_runs.receipt.run_boundary.reference_sha256",
                    "release_c_extensions.caller_runs.receipt.raw_artifact.raw_artifact_id",
                    "release_c_extensions.caller_runs.receipt.raw_artifact.raw_artifact_sha256",
                ],
                "cardinality": "many_to_one",
                "completeness": "complete",
                "unmatched_meaning": (
                    "invalid: every typed Release-C row must match exactly one embedded caller receipt"
                ),
            }
        )
    return {**document, "relationships": [*rows, *additions]}


def _field_paths(schema: pa.Schema) -> list[tuple[str, pa.Field]]:
    output: list[tuple[str, pa.Field]] = []

    def visit(field: pa.Field, prefix: str) -> None:
        output.append((prefix, field))
        if pa.types.is_struct(field.type):
            for child in field.type:
                visit(child, f"{prefix}.{child.name}")

    for field in schema:
        visit(field, field.name)
    return output


def _field_entry(domain: str, path: str, field: pa.Field) -> dict[str, Any]:
    top_level = path.split(".", 1)[0]
    derived = top_level in {
        "evidence_id",
        "event_id",
        "confidence_tier",
        "concordance_group_id",
        "concordance_status",
        "confirmation_required",
        "filter_reasons",
        "limitations",
    }
    metadata = top_level in {
        "assembly",
        "reference_sha256",
        "technical_source_role",
        "caller_name",
        "caller_version",
        "caller_run_identity_sha256",
        "resource_manifest_sha256",
        "source_alignment_sha256",
        "raw_artifact_id",
        "raw_artifact_sha256",
        "source_record_digest",
    }
    allowed_values = DOMAIN_ALLOWED_VALUES.get((domain, top_level))
    if allowed_values is None:
        allowed_values = COMMON_ALLOWED_VALUES.get(top_level)
    if not field.nullable:
        null_meaning = "field is required and cannot be null"
    else:
        null_meaning = (
            "the source caller did not emit a reliable value; null is not zero or a negative result"
        )
    return {
        "type": str(field.type),
        "nullable": field.nullable,
        "null_meaning": null_meaning,
        "source": ("embedded Release-C caller receipt and typed source-specific adapter output"),
        "derivation": (
            "lossless typed Release-C adaptation without caller voting, replicate averaging, "
            "or absent-to-reference conversion"
        ),
        "allowed_values": allowed_values,
        "evidence_class": (
            "bundle_metadata" if metadata else "pipeline_derived" if derived else "observed"
        ),
        "observed_or_inferred": "metadata" if metadata else "inferred" if derived else "observed",
        "genome_build": "GRCh38",
        "version_or_snapshot": f"release-c-schema-{RELEASE_C_SCHEMA_VERSION}",
        "confidence_semantics": (
            "confidence_tier is analytical evidence strength only; confirmation_required and "
            "limitations govern downstream use"
        ),
        "absence_meaning": (
            "no matching row has no negative meaning; inspect omissions, callability, caller "
            "receipts, technical-source role, and assessment status"
        ),
        "related_evidence": [
            f"extensions/release_c/{domain}.parquet",
            f"extensions/release_c/{MANIFEST_FILENAME}",
            "raw_artifacts.json",
            "transformations.json",
            "omissions.json",
        ],
    }


def _augment_glossary(document: dict[str, Any], stage_manifest: dict[str, Any]) -> dict[str, Any]:
    fields = document.get("fields")
    if document.get("schema_version") != SCHEMA_VERSION or not isinstance(fields, dict):
        raise PortableEvidenceError("base field_glossary.json is invalid")
    output = dict(fields)
    for domain in sorted(stage_manifest["domains"]):
        for path, field in _field_paths(schema_for(domain)):
            key = f"{domain}.{path}"
            if key in output:
                raise PortableEvidenceError("base field glossary already contains Release C")
            output[key] = _field_entry(domain, path, field)
    return {**document, "fields": output}


def _file_metadata(path: Path, *, rows: int | None = None) -> dict[str, Any]:
    if path.is_dir():
        size = sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
        digest = directory_digest(path)
    else:
        size = path.stat().st_size
        digest = sha256_file(path)
    return {"bytes": size, "rows": rows, "sha256": digest}


def _augment_inventory(
    document: dict[str, Any],
    root: Path,
    stage_manifest: dict[str, Any],
) -> dict[str, Any]:
    rows = document.get("entries")
    if document.get("schema_version") != SCHEMA_VERSION or not isinstance(rows, list):
        raise PortableEvidenceError("base inventory.json is invalid")
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in rows:
        if not isinstance(raw, dict) or not isinstance(raw.get("path"), str):
            raise PortableEvidenceError("base inventory has an invalid entry")
        row = dict(raw)
        relative = row["path"]
        if relative in seen:
            raise PortableEvidenceError("base inventory path is duplicated")
        seen.add(relative)
        target = root / relative
        if relative in {"manifest.json", "inventory.json"}:
            row.update({"bytes": None, "rows": None, "sha256": None})
        elif target.exists():
            metadata = _file_metadata(target, rows=row.get("rows"))
            row.update(metadata)
        output.append(row)
    stage_relative = f"extensions/release_c/{MANIFEST_FILENAME}"
    if stage_relative in seen:
        raise PortableEvidenceError("base inventory already contains Release-C manifest")
    stage_path = root / stage_relative
    output.append(
        {
            "authority": "metadata",
            "bytes": stage_path.stat().st_size,
            "logical_role": "Release-C domain, caller-receipt, and losslessness manifest",
            "media_type": "application/json",
            "path": stage_relative,
            "required": False,
            "rows": None,
            "sha256": sha256_file(stage_path),
        }
    )
    for domain in sorted(stage_manifest["domains"]):
        relative = f"extensions/release_c/{domain}.parquet"
        if relative in seen:
            raise PortableEvidenceError("base inventory already contains a Release-C table")
        metadata = _domain_file_metadata(stage_manifest, domain)
        output.append(
            {
                "authority": "authoritative_data",
                "bytes": metadata["bytes"],
                "logical_role": DOMAIN_ROLES[domain],
                "media_type": "application/vnd.apache.parquet",
                "path": relative,
                "required": False,
                "rows": metadata["rows"],
                "sha256": metadata["sha256"],
            }
        )
    return {**document, "entries": output}


def _augment_manifest(
    document: dict[str, Any],
    root: Path,
    stage_manifest: dict[str, Any],
    generated_at: str,
    base_manifest_sha256: str,
) -> dict[str, Any]:
    if document.get("schema_version") != SCHEMA_VERSION:
        raise PortableEvidenceError("base manifest is not .genome/1.1")
    if "release_c" in document.get("resource_provenance", {}):
        raise PortableEvidenceError("base manifest already contains Release C")
    reference_sha256 = document.get("resource_provenance", {}).get("reference_fasta_sha256")
    if reference_sha256 != stage_manifest["reference_sha256"]:
        raise PortableEvidenceError("base and Release-C reference identities differ")
    files = dict(document.get("files", {}))
    required_manifest_paths = {
        "extensions",
        "field_glossary.json",
        "inventory.json",
        "omissions.json",
        "raw_artifacts.json",
        "readmygenome",
        "readmygenome.md",
        "relationships.json",
        "transformations.json",
    }
    if not required_manifest_paths.issubset(files):
        raise PortableEvidenceError("base manifest lacks a required metadata boundary")
    for relative in required_manifest_paths:
        previous_rows = files[relative].get("rows") if isinstance(files[relative], dict) else None
        files[relative] = _file_metadata(root / relative, rows=previous_rows)

    provenance = dict(document["resource_provenance"])
    provenance["release_c"] = {
        "profile_version": RELEASE_C_PROFILE_VERSION,
        "extension_schema_version": RELEASE_C_SCHEMA_VERSION,
        "extension_set_identity_sha256": stage_manifest["extension_set_identity_sha256"],
        "reference_sha256": stage_manifest["reference_sha256"],
        "base_bundle_manifest_sha256": base_manifest_sha256,
        "logical_subject_count": 1,
        "source_dataset_count": stage_manifest["source_dataset_count"],
        "source_rows_removed": 0,
        "caller_or_replicate_rows_collapsed": 0,
        "input_gate_count": len(stage_manifest["input_gates"]),
        "caller_run_count": len(stage_manifest["caller_runs"]),
        "domains": {
            domain: {
                "rows": metadata["rows"],
                "sha256": metadata["sha256"],
                "canonical_row_multiset": metadata["canonical_row_multiset"],
                "caller_boundaries": metadata["caller_boundaries"],
            }
            for domain, metadata in sorted(stage_manifest["domains"].items())
        },
        "privacy": {
            "subject_identifier_persisted": False,
            "private_paths_persisted": False,
            "private_manifest_digest_persisted": False,
        },
    }
    annotations = dict(document.get("annotation_sources", {}))
    annotations["read_backed_extensions"] = (
        "Release C schema 1; exact caller/resource/alignment receipts; source-specific rows"
    )
    consensus = dict(document.get("consensus_definitions", {}))
    consensus["read_backed_extensions"] = {
        "rule": (
            "retain every caller and technical-source observation; annotate concordance without "
            "voting, averaging, absent-to-reference conversion, or silent genotype replacement"
        ),
        "confidence_tiers": ["SUPPORTED", "REVIEW", "UNRESOLVED"],
        "absence_is_negative": False,
        "clinical_diagnosis_inferred": False,
        "confirmation_required_field": "confirmation_required",
    }
    privacy = dict(document.get("privacy", {}))
    privacy.update(
        {
            "release_c_subject_identifier_embedded": False,
            "release_c_private_paths_embedded": False,
        }
    )
    pipeline_version = str(document.get("pipeline_version", ""))
    if "+release-c." in pipeline_version:
        raise PortableEvidenceError("base pipeline version is already Release-C enriched")
    return {
        **document,
        "generated_at": generated_at,
        "pipeline_version": f"{pipeline_version}+release-c.{RELEASE_C_SCHEMA_VERSION}",
        "annotation_sources": annotations,
        "consensus_definitions": consensus,
        "resource_provenance": provenance,
        "privacy": privacy,
        "files": files,
    }


def _first_query_term(stage_root: Path, stage_manifest: dict[str, Any]) -> str:
    domain = next(iter(sorted(stage_manifest["domains"])))
    row = pq.read_table(stage_root / f"{domain}.parquet").slice(0, 1).to_pylist()[0]
    if domain == "structural_variants":
        return f"{row['chrom1']}:{row['pos1']}"
    if domain in {"copy_number", "repeat_expansions"}:
        return f"{row['chrom']}:{row['start_pos']}"
    if domain == "hla_typing":
        return str(row["locus"])
    if domain == "complex_pharmacogenes":
        return str(row["gene"])
    return f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"


def validate_enriched_bundle(
    bundle: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    run_read_backed_query: bool = True,
) -> dict[str, Any]:
    root = Path(bundle).expanduser().resolve()
    core = validate_bundle(
        root,
        reference_fasta=reference_fasta,
        run_profile_queries=True,
    )
    manifest = _load_json(root / "manifest.json", "enriched manifest")
    expected_reference = manifest.get("resource_provenance", {}).get("reference_fasta_sha256")
    release_c = validate_release_c_extension_stage(
        root / "extensions" / "release_c",
        expected_reference_sha256=expected_reference,
    )
    release_c_manifest = _load_json(
        root / "extensions" / "release_c" / MANIFEST_FILENAME,
        "embedded Release-C manifest",
    )
    claimed = manifest.get("resource_provenance", {}).get("release_c")
    if not isinstance(claimed, dict) or claimed.get(
        "extension_set_identity_sha256"
    ) != release_c.get("extension_set_identity_sha256"):
        raise PortableEvidenceError("bundle manifest and Release-C stage identity differ")
    omissions = _load_json(root / "omissions.json", "enriched omissions")["omissions"]
    scopes = [row.get("scope") for row in omissions if isinstance(row, dict)]
    if GENERIC_READ_BACKED_OMISSION in scopes:
        raise PortableEvidenceError("enriched bundle retains its obsolete generic omission")
    for domain in DOMAINS:
        expected_count = 0 if domain in release_c_manifest["domains"] else 1
        if scopes.count(DOMAIN_OMISSIONS[domain]) != expected_count:
            raise PortableEvidenceError("Release-C domain omission accounting differs")
    expected_scope_omission = _portable_core_scope_omission(release_c_manifest)
    observed_scope_omissions = [
        row
        for row in omissions
        if isinstance(row, dict) and row.get("scope") == PORTABLE_CORE_SCOPE_OMISSION
    ]
    if observed_scope_omissions != (
        [] if expected_scope_omission is None else [expected_scope_omission]
    ):
        raise PortableEvidenceError("Release-C portable-scope omission accounting differs")
    query_report: dict[str, Any] = {"status": "not_run"}
    if run_read_backed_query:
        term = _first_query_term(root / "extensions" / "release_c", release_c_manifest)
        result = query_bundle(root, "read-backed", term=term, limit=10)
        if result.get("count", 0) < 1:
            raise PortableEvidenceError("read-backed query smoke test returned no evidence")
        query_report = {
            "status": "PASS",
            "term": term,
            "rows": result["count"],
            "truncated": result["truncated"],
        }
    return {
        "status": "PASS",
        "core": core,
        "release_c": release_c,
        "query": query_report,
    }


def build_enriched_bundle(
    base_bundle: str | Path,
    release_c_stage: str | Path,
    output_bundle: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    materialization: str = "copy",
    generated_at: str | None = None,
    minimum_free_bytes: int = 0,
    validate_base: bool = True,
) -> dict[str, Any]:
    base = Path(base_bundle).expanduser().resolve()
    stage = Path(release_c_stage).expanduser().resolve()
    output = Path(output_bundle).expanduser().resolve()
    if output.exists():
        raise PortableEvidenceError("refusing to replace an existing enriched bundle")
    if not base.is_dir() or not stage.is_dir():
        raise PortableEvidenceError("base bundle or Release-C stage is missing")
    base_manifest_path = base / "manifest.json"
    base_manifest_sha256 = sha256_file(base_manifest_path)
    base_manifest = _load_json(base_manifest_path, "base manifest")
    reference_sha256 = base_manifest.get("resource_provenance", {}).get("reference_fasta_sha256")
    stage_validation = validate_release_c_extension_stage(
        stage,
        expected_reference_sha256=reference_sha256,
    )
    stage_manifest = _load_json(stage / MANIFEST_FILENAME, "Release-C stage manifest")
    if validate_base:
        validate_bundle(
            base,
            reference_fasta=reference_fasta,
            run_profile_queries=True,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(output.parent).free
    if minimum_free_bytes < 0:
        raise PortableEvidenceError("minimum free bytes cannot be negative")
    estimated_new_bytes = sum(item.stat().st_size for item in _regular_tree_files(stage))
    if materialization == "copy":
        estimated_new_bytes += sum(item.stat().st_size for item in _regular_tree_files(base))
    if free_bytes - estimated_new_bytes < minimum_free_bytes:
        raise PortableEvidenceError("enriched bundle would cross the configured free-space floor")
    timestamp = _canonical_timestamp(generated_at)
    base_documents = {
        name: _load_json(base / name, f"base {name}")
        for name in (
            "inventory.json",
            "raw_artifacts.json",
            "transformations.json",
            "omissions.json",
            "relationships.json",
            "field_glossary.json",
        )
    }
    with atomic_output_directory(output) as staging:
        import_report = _materialize_base(base, staging, mode=materialization)
        release_c_destination = staging / "extensions" / "release_c"
        _copy_regular_tree(stage, release_c_destination)
        if directory_digest(stage) != directory_digest(release_c_destination):
            raise PortableEvidenceError("copied Release-C stage identity changed")

        query_source = Path(__file__).resolve().parents[1] / "query.py"
        shutil.copy2(query_source, staging / "readmygenome" / "scripts" / "query.py")
        _append_document(
            base / "readmygenome.md",
            staging / "readmygenome.md",
            README_APPENDIX,
            "## Read-backed extensions",
        )
        _append_document(
            base / "readmygenome" / "SKILL.md",
            staging / "readmygenome" / "SKILL.md",
            SKILL_APPENDIX,
            "## Release-C read-backed evidence",
        )
        _append_document(
            base / "readmygenome" / "reference" / "query_catalog.md",
            staging / "readmygenome" / "reference" / "query_catalog.md",
            QUERY_CATALOG_APPENDIX,
            "## Read-backed extensions",
        )

        atomic_write_json(
            staging / "raw_artifacts.json",
            _augment_raw_artifacts(base_documents["raw_artifacts.json"], stage_manifest),
        )
        atomic_write_json(
            staging / "transformations.json",
            _augment_transformations(
                base_documents["transformations.json"],
                stage_manifest,
                timestamp,
                base_manifest_sha256,
                import_report,
            ),
        )
        atomic_write_json(
            staging / "omissions.json",
            _augment_omissions(base_documents["omissions.json"], stage_manifest),
        )
        atomic_write_json(
            staging / "relationships.json",
            _augment_relationships(base_documents["relationships.json"], stage_manifest),
        )
        atomic_write_json(
            staging / "field_glossary.json",
            _augment_glossary(base_documents["field_glossary.json"], stage_manifest),
        )
        atomic_write_json(
            staging / "inventory.json",
            _augment_inventory(base_documents["inventory.json"], staging, stage_manifest),
        )
        atomic_write_json(
            staging / "manifest.json",
            _augment_manifest(
                base_manifest,
                staging,
                stage_manifest,
                timestamp,
                base_manifest_sha256,
            ),
        )
        validation = validate_enriched_bundle(
            staging,
            reference_fasta=reference_fasta,
            run_read_backed_query=True,
        )
        if sha256_file(base_manifest_path) != base_manifest_sha256:
            raise PortableEvidenceError("base bundle manifest changed during enrichment")
    return {
        "status": "PASS",
        "output": str(output),
        "base_bundle_manifest_sha256": base_manifest_sha256,
        "extension_set_identity_sha256": stage_validation["extension_set_identity_sha256"],
        "materialization": import_report,
        "validation": validation,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Release-C .genome/1.1 bundle extension")
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build")
    build.add_argument("--base-bundle", type=Path, required=True)
    build.add_argument("--release-c-stage", type=Path, required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--reference-fasta", type=Path)
    build.add_argument("--materialization", choices=["copy", "hardlink"], default="copy")
    build.add_argument("--generated-at")
    build.add_argument("--minimum-free-bytes", type=int, default=0)
    validate = subparsers.add_parser("validate")
    validate.add_argument("--bundle", type=Path, required=True)
    validate.add_argument("--reference-fasta", type=Path)
    validate.add_argument("--skip-query", action="store_true")
    args = parser.parse_args()
    if args.command == "build":
        report = build_enriched_bundle(
            args.base_bundle,
            args.release_c_stage,
            args.output,
            reference_fasta=args.reference_fasta,
            materialization=args.materialization,
            generated_at=args.generated_at,
            minimum_free_bytes=args.minimum_free_bytes,
        )
    else:
        report = validate_enriched_bundle(
            args.bundle,
            reference_fasta=args.reference_fasta,
            run_read_backed_query=not args.skip_query,
        )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
