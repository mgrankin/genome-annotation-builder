"""Select a bounded, deterministic set of small variants for local read-back.

The selector operates on an already validated ``.genome/1.1`` bundle.  It
prioritizes clinical findings and pedigree/caller review signals, never treats
an absent row as reference evidence, and publishes a path-free receipt last.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import tempfile
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import duckdb
import pyarrow.parquet as pq

from ..constants import SCHEMA_VERSION
from ..util import directory_digest, native_pipeline_version, sha256_file
from .portable_model import CHROM_RANK, FORBIDDEN_PORTABLE_TEXT, SHA256
from .readback_capture import TARGET_COLUMNS

TARGET_RECEIPT_SCHEMA_VERSION = 1
CONCRETE_ALLELE = re.compile(r"^[ACGT]+$")
DEFAULT_MAXIMUM_TARGETS = 20_000
MAXIMUM_ALLOWED_TARGETS = 100_000
CATEGORY_PRIORITY = {
    "clinical_finding": 0,
    "de_novo_candidate": 10,
    "parental_mosaicism_candidate": 20,
    "parental_dropout_suspected": 30,
    "child_dropout_suspected": 30,
    "technical_replicate_discordant": 35,
    "loco_source_genotype_discordant": 40,
    "loco_mendelian_conflict": 45,
    "mendelian_conflict": 50,
    "caller_discordant": 60,
}


class ReadbackTargetError(ValueError):
    """A target selection or its completion receipt is invalid."""


def _fail(message: str) -> None:
    raise ReadbackTargetError(message)


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _load_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        _fail(f"{label}_missing")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ReadbackTargetError(f"{label}_unreadable") from exc
    if not isinstance(value, dict):
        _fail(f"{label}_invalid")
    return value


def _atomic_bytes_no_replace(path: Path, payload: bytes, label: str) -> None:
    if path.exists():
        _fail(f"{label}_already_exists")
    path.parent.mkdir(parents=True, exist_ok=True)
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
            raise ReadbackTargetError(f"{label}_already_exists") from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _check_file_identity(path: Path, spec: object, label: str) -> str:
    if path.is_symlink() or not path.is_file() or not isinstance(spec, dict):
        _fail(f"{label}_missing")
    digest = sha256_file(path)
    if path.stat().st_size != spec.get("bytes") or digest != spec.get("sha256"):
        _fail(f"{label}_identity_changed")
    return digest


def _verify_legacy_empty_findings(bundle: Path, manifest: dict[str, Any]) -> None:
    """Native releases through 0.2.7 omitted the table iff there were no decisions."""
    producer = native_pipeline_version(manifest.get("pipeline_version"))
    if not producer or producer >= (0, 2, 8):
        _fail("clinical_findings_missing")

    files = manifest["files"]
    documents = {}
    for name in ("inventory", "transformations"):
        path = bundle / f"{name}.json"
        _check_file_identity(path, files.get(path.name), name)
        documents[name] = _load_json(path, name)

    entries = documents["inventory"].get("entries")
    stages = documents["transformations"].get("stages")
    if not isinstance(entries, list) or not all(
        isinstance(entry, dict) and isinstance(entry.get("path"), str) for entry in entries
    ):
        _fail("clinical_findings_inventory_invalid")
    # The native producer emits one finding and one interpretation for *every*
    # decision, including decisions that did not pass its actionability gates.
    decision_outputs = {
        "clinical_findings.parquet",
        "interpretations.parquet",
        "extensions/acmg_secondary_findings_evidence.parquet",
    }
    inventory = {entry["path"] for entry in entries}
    if (
        not {"variants.parquet", "transformations.json", "extensions/family_evidence.parquet"}
        <= inventory
    ):
        _fail("clinical_findings_inventory_invalid")
    registered = set(files) | inventory
    if decision_outputs & registered or any(os.path.lexists(bundle / p) for p in decision_outputs):
        _fail("clinical_findings_conflicting_inventory")

    if not isinstance(stages, list) or not all(isinstance(stage, dict) for stage in stages):
        _fail("clinical_findings_evaluation_unproven")
    acmg = [s for s in stages if s.get("stage_id") == "acmg_secondary_findings_policy_evaluation"]
    resources = manifest.get("resource_provenance")
    policy = resources.get("acmg_secondary_findings") if isinstance(resources, dict) else None
    if len(acmg) != 1 or not isinstance(policy, dict):
        _fail("clinical_findings_evaluation_unproven")
    for resource, field in (
        ("ensembl_vep", "cache_tree_sha256"),
        ("clinvar", "normalized_vcf_sha256"),
    ):
        source = resources.get(resource)
        identity = source.get(field) if isinstance(source, dict) else None
        if not isinstance(identity, str) or not SHA256.fullmatch(identity):
            _fail("clinical_findings_evaluation_unproven")
    stage = acmg[0]
    digest = policy.get("canonical_sha256")
    sources = manifest.get("annotation_sources")
    if not isinstance(sources, dict) or sources.get("actionability_annotation") != (
        f"ACMG SF v{policy.get('version')}; ClinGen snapshot {policy.get('snapshot')}; "
        f"canonical-sha256:{digest}"
    ):
        _fail("clinical_findings_evaluation_unproven")
    zero_counts = (
        obj.get(field)
        for obj in (stage, policy)
        for field in ("decision_rows", "actionable_candidate_variants")
    )
    if (
        stage.get("status") != "completed"
        or stage.get("software") != "genome-annotation-builder"
        or stage.get("software_version") != ".".join(map(str, producer))
        or not isinstance(digest, str)
        or not SHA256.fullmatch(digest)
        or not isinstance(stage.get("inputs"), list)
        or f"sha256:{digest}" not in stage["inputs"]
        or policy.get("decision_status_counts") != {}
        or any(type(n) is not int or n != 0 for n in zero_counts)
    ):
        _fail("clinical_findings_evaluation_unproven")


def _bundle_boundary(bundle: Path) -> dict[str, Any]:
    if bundle.is_symlink() or not bundle.is_dir():
        _fail("source_bundle_invalid")
    manifest_path = bundle / "manifest.json"
    clinical = bundle / "clinical_findings.parquet"
    family = bundle / "extensions/family_evidence.parquet"
    if not family.is_dir():
        _fail("source_bundle_required_tables_missing")
    manifest = _load_json(manifest_path, "source_bundle_manifest")
    files = manifest.get("files")
    if (
        manifest.get("schema_version") != SCHEMA_VERSION
        or manifest.get("genome_build") != "GRCh38"
        or not isinstance(files, dict)
    ):
        _fail("source_bundle_manifest_boundary_invalid")
    variants_spec = files.get("variants.parquet")
    extensions_spec = files.get("extensions")
    if not all(isinstance(value, dict) for value in (variants_spec, extensions_spec)):
        _fail("source_bundle_manifest_file_boundary_missing")
    if clinical.name in files or os.path.lexists(clinical):
        spec = files.get(clinical.name)
        clinical_sha256 = _check_file_identity(clinical, spec, "clinical_findings")
        rows = spec.get("rows")
        if rows is not None and (
            type(rows) is not int or pq.read_metadata(clinical).num_rows != rows
        ):
            _fail("clinical_findings_row_count_changed")
    else:
        _verify_legacy_empty_findings(bundle, manifest)
        clinical_sha256 = None
    family_sha256 = directory_digest(family)
    return {
        "manifest": manifest,
        "manifest_sha256": sha256_file(manifest_path),
        "clinical_path": clinical if clinical_sha256 is not None else None,
        "clinical_sha256": clinical_sha256,
        "family_path": family,
        "family_sha256": family_sha256,
        "expected_variant_rows": int(variants_spec.get("rows", -1)),
        "extensions_sha256": extensions_spec.get("sha256"),
    }


def _candidate_categories(row: Mapping[str, Any], clinical_ids: set[str]) -> list[str]:
    categories: list[str] = []
    if row["variant_id"] in clinical_ids:
        categories.append("clinical_finding")
    for field, category in (
        ("de_novo_candidate", "de_novo_candidate"),
        ("parental_mosaicism_candidate", "parental_mosaicism_candidate"),
        ("parental_dropout_suspected", "parental_dropout_suspected"),
        ("child_dropout_suspected", "child_dropout_suspected"),
    ):
        if row[field]:
            categories.append(category)
    if row["technical_replicate_concordance"] == "discordant":
        categories.append("technical_replicate_discordant")
    if row["loco_status"] == "source_genotype_discordant":
        categories.append("loco_source_genotype_discordant")
    if row["loco_status"] == "mendelian_conflict":
        categories.append("loco_mendelian_conflict")
    if row["loco_status"] == "de_novo_candidate":
        categories.append("de_novo_candidate")
    if row["mendelian_status"] == "conflict":
        categories.append("mendelian_conflict")
    if row["caller_concordance"] == "discordant":
        categories.append("caller_discordant")
    return sorted(set(categories), key=lambda value: (CATEGORY_PRIORITY[value], value))


def _select_candidates(boundary: Mapping[str, Any], maximum_targets: int) -> dict[str, Any]:
    if maximum_targets < 1 or maximum_targets > MAXIMUM_ALLOWED_TARGETS:
        _fail("maximum_targets_out_of_range")
    family_glob = str(Path(boundary["family_path"]) / "chrom=*" / "part-*.parquet")
    connection = duckdb.connect()
    try:
        clinical_ids: set[str] = set()
        if boundary["clinical_path"] is not None:
            clinical_rows = connection.execute(
                "SELECT DISTINCT variant_id FROM read_parquet(?) WHERE variant_id IS NOT NULL",
                [str(boundary["clinical_path"])],
            ).fetchall()
            clinical_ids = {str(row[0]) for row in clinical_rows}
        family_rows = int(
            connection.execute("SELECT count(*) FROM read_parquet(?)", [family_glob]).fetchone()[0]
        )
        if family_rows != boundary["expected_variant_rows"]:
            _fail("family_evidence_variant_cardinality_mismatch")
        connection.execute("CREATE TEMP TABLE required_clinical(variant_id VARCHAR PRIMARY KEY)")
        if clinical_ids:
            connection.executemany(
                "INSERT INTO required_clinical VALUES (?)",
                [(variant_id,) for variant_id in sorted(clinical_ids)],
            )
        query = """
            SELECT
              f.variant_id, f.chrom, f.pos, f.ref, f.alt,
              f.de_novo_candidate, f.parental_mosaicism_candidate,
              f.parental_dropout_suspected, f.child_dropout_suspected,
              f.mendelian_status, f.loco_status, f.caller_concordance,
              f.technical_replicate_concordance, f.complex_region,
              f.difficult_region_flags, f.phase_confidence,
              (c.variant_id IS NOT NULL) AS clinical_finding
            FROM read_parquet(?) AS f
            LEFT JOIN required_clinical AS c USING (variant_id)
            WHERE c.variant_id IS NOT NULL
               OR f.de_novo_candidate
               OR f.parental_mosaicism_candidate
               OR f.parental_dropout_suspected
               OR f.child_dropout_suspected
               OR f.mendelian_status = 'conflict'
               OR f.loco_status IN (
                    'source_genotype_discordant', 'mendelian_conflict', 'de_novo_candidate'
                  )
               OR f.caller_concordance = 'discordant'
               OR f.technical_replicate_concordance = 'discordant'
        """
        result = connection.execute(query, [family_glob])
        columns = [column[0] for column in result.description]
        rows = [dict(zip(columns, values, strict=True)) for values in result.fetchall()]
    finally:
        connection.close()
    seen_variants: set[str] = set()
    seen_alleles: set[tuple[str, int, str, str]] = set()
    candidates: list[dict[str, Any]] = []
    category_counts: Counter[str] = Counter()
    for row in rows:
        variant_id = str(row["variant_id"])
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        ref = str(row["ref"])
        alt = str(row["alt"])
        if (
            variant_id in seen_variants
            or chrom not in CHROM_RANK
            or pos < 1
            or CONCRETE_ALLELE.fullmatch(ref) is None
            or CONCRETE_ALLELE.fullmatch(alt) is None
            or ref == alt
        ):
            _fail("candidate_variant_boundary_invalid")
        allele = (chrom, pos, ref, alt)
        if allele in seen_alleles:
            _fail("candidate_allele_duplicated")
        seen_variants.add(variant_id)
        seen_alleles.add(allele)
        categories = _candidate_categories(row, clinical_ids)
        if not categories:
            _fail("candidate_has_no_selection_reason")
        category_counts.update(categories)
        difficult = row["difficult_region_flags"] or []
        local_complexity = (
            "complex_or_difficult_region"
            if bool(row["complex_region"]) or bool(difficult)
            else "not_flagged"
        )
        phase_confidence = str(row["phase_confidence"] or "not_phased")
        stable_rank = hashlib.sha256(variant_id.encode()).hexdigest()
        candidates.append(
            {
                "target_id": "rb-"
                + hashlib.sha256(f"{chrom}:{pos}:{ref}:{alt}".encode()).hexdigest()[:24],
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "candidate_category": "+".join(categories),
                "local_complexity_status": local_complexity,
                "local_haplotype_status": phase_confidence,
                "priority": min(CATEGORY_PRIORITY[value] for value in categories),
                "stable_rank": stable_rank,
                "required_clinical": "clinical_finding" in categories,
            }
        )
    observed_clinical = {row["variant_id"] for row in rows if row["clinical_finding"]}
    if observed_clinical != clinical_ids:
        _fail("clinical_findings_missing_from_family_evidence")
    required_clinical = sum(candidate["required_clinical"] for candidate in candidates)
    if required_clinical > maximum_targets:
        _fail("maximum_targets_cannot_retain_all_clinical_findings")
    ranked = sorted(candidates, key=lambda row: (row["priority"], row["stable_rank"]))
    selected = ranked[:maximum_targets]
    if sum(row["required_clinical"] for row in selected) != required_clinical:
        _fail("clinical_finding_would_be_truncated")
    selected.sort(key=lambda row: (CHROM_RANK[row["chrom"]], row["pos"], row["ref"], row["alt"]))
    return {
        "rows": selected,
        "family_rows": family_rows,
        "candidate_rows": len(candidates),
        "selected_rows": len(selected),
        "truncated_rows": len(candidates) - len(selected),
        "required_clinical_rows": required_clinical,
        "category_counts": dict(sorted(category_counts.items())),
    }


def _render_targets(rows: list[dict[str, Any]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=TARGET_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({column: row[column] for column in TARGET_COLUMNS})
    return output.getvalue().encode()


def create_readback_targets(
    bundle_path: str | Path,
    output_path: str | Path,
    receipt_path: str | Path,
    *,
    maximum_targets: int = DEFAULT_MAXIMUM_TARGETS,
) -> dict[str, Any]:
    """Publish the target TSV first and its path-free completion receipt last."""

    bundle = Path(bundle_path).resolve()
    output = Path(output_path)
    receipt_file = Path(receipt_path)
    if output.exists() or receipt_file.exists():
        _fail("target_or_receipt_already_exists")
    boundary = _bundle_boundary(bundle)
    selection = _select_candidates(boundary, maximum_targets)
    payload = _render_targets(selection["rows"])
    target_sha256 = hashlib.sha256(payload).hexdigest()
    receipt: dict[str, Any] = {
        "schema_version": TARGET_RECEIPT_SCHEMA_VERSION,
        "status": "PASS",
        "assembly": "GRCh38",
        "source_bundle_schema_version": boundary["manifest"]["schema_version"],
        "source_pipeline_version": boundary["manifest"].get("pipeline_version"),
        "source_bundle_manifest_sha256": boundary["manifest_sha256"],
        "source_tables": {
            "clinical_findings_sha256": boundary["clinical_sha256"],
            "family_evidence_sha256": boundary["family_sha256"],
            "extensions_manifest_sha256": boundary["extensions_sha256"],
        },
        "selection_policy": {
            "maximum_targets": maximum_targets,
            "maximum_allowed_targets": MAXIMUM_ALLOWED_TARGETS,
            "category_priority": CATEGORY_PRIORITY,
            "within_priority_order": "sha256_of_source_variant_id",
            "output_order": "canonical_genomic",
            "clinical_findings_must_not_be_truncated": True,
        },
        "source_variant_rows": selection["family_rows"],
        "candidate_rows": selection["candidate_rows"],
        "selected_rows": selection["selected_rows"],
        "truncated_rows": selection["truncated_rows"],
        "required_clinical_rows": selection["required_clinical_rows"],
        "category_counts": selection["category_counts"],
        "target_manifest_sha256": target_sha256,
        "target_manifest_bytes": len(payload),
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    if boundary["clinical_sha256"] is None:
        receipt["source_tables"]["clinical_findings_state"] = "legacy_evaluated_empty"
    receipt["receipt_identity_sha256"] = _canonical_sha256(receipt)
    rendered = json.dumps(receipt, sort_keys=True)
    if FORBIDDEN_PORTABLE_TEXT.search(rendered):
        _fail("target_receipt_contains_private_path")
    _atomic_bytes_no_replace(output, payload, "target_manifest")
    _atomic_bytes_no_replace(
        receipt_file,
        (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode(),
        "target_receipt",
    )
    return receipt


def validate_readback_targets(
    target_path: str | Path,
    receipt_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
) -> dict[str, Any]:
    """Validate target bytes, receipt identity, and optionally the source bundle."""

    target = Path(target_path)
    receipt_file = Path(receipt_path)
    receipt = _load_json(receipt_file, "target_receipt")
    expected = {
        "schema_version",
        "status",
        "assembly",
        "source_bundle_schema_version",
        "source_pipeline_version",
        "source_bundle_manifest_sha256",
        "source_tables",
        "selection_policy",
        "source_variant_rows",
        "candidate_rows",
        "selected_rows",
        "truncated_rows",
        "required_clinical_rows",
        "category_counts",
        "target_manifest_sha256",
        "target_manifest_bytes",
        "subject_identifier_persisted",
        "private_paths_persisted",
        "receipt_identity_sha256",
    }
    if (
        set(receipt) != expected
        or receipt.get("schema_version") != TARGET_RECEIPT_SCHEMA_VERSION
        or receipt.get("status") != "PASS"
        or receipt.get("assembly") != "GRCh38"
        or receipt.get("subject_identifier_persisted") is not False
        or receipt.get("private_paths_persisted") is not False
    ):
        _fail("target_receipt_schema_invalid")
    identity = receipt.pop("receipt_identity_sha256")
    if not isinstance(identity, str) or SHA256.fullmatch(identity) is None:
        _fail("target_receipt_identity_invalid")
    if identity != _canonical_sha256(receipt):
        _fail("target_receipt_identity_changed")
    receipt["receipt_identity_sha256"] = identity
    tables = receipt["source_tables"]
    if not isinstance(tables, dict) or "clinical_findings_sha256" not in tables:
        _fail("target_receipt_clinical_source_invalid")
    digest = tables["clinical_findings_sha256"]
    if digest is None:
        if tables.get("clinical_findings_state") != "legacy_evaluated_empty" or (
            receipt["required_clinical_rows"] != 0
        ):
            _fail("target_receipt_clinical_source_invalid")
    elif (
        not isinstance(digest, str)
        or not SHA256.fullmatch(digest)
        or "clinical_findings_state" in tables
    ):
        _fail("target_receipt_clinical_source_invalid")
    if target.is_symlink() or not target.is_file():
        _fail("target_manifest_missing")
    if (
        sha256_file(target) != receipt["target_manifest_sha256"]
        or target.stat().st_size != receipt["target_manifest_bytes"]
    ):
        _fail("target_manifest_identity_changed")
    try:
        with target.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != TARGET_COLUMNS:
                _fail("target_manifest_schema_invalid")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ReadbackTargetError("target_manifest_unreadable") from exc
    if len(rows) != receipt["selected_rows"]:
        _fail("target_manifest_row_count_changed")
    ids = [row["target_id"] for row in rows]
    alleles = [(row["chrom"], row["pos"], row["ref"], row["alt"]) for row in rows]
    if len(ids) != len(set(ids)) or len(alleles) != len(set(alleles)):
        _fail("target_manifest_duplicate_rows")
    if bundle_path is not None:
        boundary = _bundle_boundary(Path(bundle_path).resolve())
        if (
            boundary["manifest_sha256"] != receipt["source_bundle_manifest_sha256"]
            or boundary["clinical_sha256"] != receipt["source_tables"]["clinical_findings_sha256"]
            or boundary["family_sha256"] != receipt["source_tables"]["family_evidence_sha256"]
        ):
            _fail("source_bundle_identity_changed")
    rendered = json.dumps(receipt, sort_keys=True)
    if FORBIDDEN_PORTABLE_TEXT.search(rendered):
        _fail("target_receipt_contains_private_path")
    return {
        "schema_version": TARGET_RECEIPT_SCHEMA_VERSION,
        "status": "PASS",
        "targets": len(rows),
        "target_manifest_sha256": receipt["target_manifest_sha256"],
        "receipt_identity_sha256": identity,
        "source_bundle_revalidated": bundle_path is not None,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create")
    create.add_argument("--bundle", type=Path, required=True)
    create.add_argument("--output", type=Path, required=True)
    create.add_argument("--receipt", type=Path, required=True)
    create.add_argument("--maximum-targets", type=int, default=DEFAULT_MAXIMUM_TARGETS)
    validate = subparsers.add_parser("validate")
    validate.add_argument("--targets", type=Path, required=True)
    validate.add_argument("--receipt", type=Path, required=True)
    validate.add_argument("--bundle", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "create":
            result = create_readback_targets(
                args.bundle,
                args.output,
                args.receipt,
                maximum_targets=args.maximum_targets,
            )
            result = {
                "schema_version": result["schema_version"],
                "status": result["status"],
                "candidate_rows": result["candidate_rows"],
                "selected_rows": result["selected_rows"],
                "truncated_rows": result["truncated_rows"],
                "target_manifest_sha256": result["target_manifest_sha256"],
                "receipt_identity_sha256": result["receipt_identity_sha256"],
                "subject_identifier_persisted": False,
                "private_paths_persisted": False,
            }
        else:
            result = validate_readback_targets(args.targets, args.receipt, bundle_path=args.bundle)
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except Exception:
        print(
            "read-back target operation failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
