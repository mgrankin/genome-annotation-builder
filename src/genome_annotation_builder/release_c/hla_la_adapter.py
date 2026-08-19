"""Strict adapter for HLA*LA ``R1_bestguess_G.txt`` output."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .hla_common import (
    PRIMARY_HLA_LOCI,
    finite_float,
    parse_hla_allele,
    semantic_record_digest,
    two_field_allele,
)
from .portable_model import (
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    sha256_file,
    validate_rows,
    write_extension,
)

EXPECTED_COLUMNS = [
    "Locus",
    "Chromosome",
    "Allele",
    "Q1",
    "Q2",
    "AverageCoverage",
    "CoverageFirstDecile",
    "MinimumCoverage",
    "proportionkMersCovered",
    "LocusAvgColumnError",
    "NColumns_UnaccountedAllele_fGT0.2",
    "perfectG",
]
LOCUS_LEVEL_METRICS = (
    "average_coverage",
    "coverage_first_decile",
    "minimum_coverage",
    "locus_avg_column_error",
    "unaccounted_columns_fraction_gt_0_2",
)


@dataclass(frozen=True)
class HlaLaProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    caller_database_name: str
    caller_database_version: str
    caller_database_identity_sha256: str
    nomenclature_database_name: str
    nomenclature_version: str
    nomenclature_identity_sha256: str
    database_age_limitation: bool = True
    caller_name: str = "HLA*LA"
    caller_version: str = "1.0.4"


def _read_rows(path: str | Path) -> list[dict[str, str]]:
    try:
        with Path(path).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != EXPECTED_COLUMNS:
                raise PortableEvidenceError("HLA*LA best-guess column schema differs")
            rows = list(reader)
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("cannot read HLA*LA best-guess artifact") from exc
    except csv.Error as exc:
        raise PortableEvidenceError("malformed HLA*LA best-guess table") from exc
    if not rows:
        raise PortableEvidenceError("HLA*LA best-guess table is empty")
    if any(None in row for row in rows):
        raise PortableEvidenceError("HLA*LA best-guess row width differs")
    return rows


def _nonnegative(value: str, label: str) -> float:
    parsed = finite_float(value, label)
    if parsed < 0:
        raise PortableEvidenceError(f"negative {label}")
    return parsed


def _parse_perfect_g(value: str) -> bool:
    if value not in {"0", "1"}:
        raise PortableEvidenceError("HLA*LA perfectG must be zero or one")
    return value == "1"


def _parse_source_row(source: dict[str, str]) -> tuple[str, int, str | None, dict[str, Any]]:
    gene = source["Locus"].strip()
    if not gene or not gene.isalnum() or gene.upper() != gene:
        raise PortableEvidenceError("invalid HLA*LA locus")
    if source["Chromosome"] not in {"1", "2"}:
        raise PortableEvidenceError("HLA*LA chromosome-copy slot must be one or two")
    copy = int(source["Chromosome"])
    allele = parse_hla_allele(source["Allele"], expected_gene=gene)
    proportion = finite_float(source["proportionkMersCovered"], "HLA*LA proportionkMersCovered")
    if proportion != -1 and not 0 <= proportion <= 1:
        raise PortableEvidenceError("HLA*LA k-mer coverage proportion is out of bounds")
    metrics: dict[str, Any] = {
        "q1": finite_float(source["Q1"], "HLA*LA Q1"),
        "q2": finite_float(source["Q2"], "HLA*LA Q2"),
        "average_coverage": _nonnegative(source["AverageCoverage"], "HLA*LA average coverage"),
        "coverage_first_decile": _nonnegative(
            source["CoverageFirstDecile"], "HLA*LA coverage first decile"
        ),
        "minimum_coverage": _nonnegative(source["MinimumCoverage"], "HLA*LA minimum coverage"),
        "proportion_kmers_covered": proportion,
        "locus_avg_column_error": _nonnegative(
            source["LocusAvgColumnError"], "HLA*LA locus average column error"
        ),
        "unaccounted_columns_fraction_gt_0_2": _nonnegative(
            source["NColumns_UnaccountedAllele_fGT0.2"],
            "HLA*LA unaccounted-column fraction",
        ),
        "perfect_g": _parse_perfect_g(source["perfectG"]),
    }
    return gene, copy, allele, metrics


def _require_shared_locus_metrics(gene: str, first: dict[str, Any], second: dict[str, Any]) -> None:
    for field in LOCUS_LEVEL_METRICS:
        if first[field] != second[field]:
            raise PortableEvidenceError(f"HLA*LA {gene} repeated locus metrics disagree")


def parse_hla_la_bestguess(
    input_table: str | Path,
    provenance: HlaLaProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not provenance.database_age_limitation:
        raise PortableEvidenceError("pinned HLA*LA graph must retain its database-age limitation")
    source_rows = _read_rows(input_table)
    raw_artifact_sha256 = sha256_file(Path(input_table))
    by_locus: dict[str, dict[int, tuple[str | None, dict[str, Any], dict[str, str]]]] = {}
    for source in source_rows:
        gene, copy, allele, metrics = _parse_source_row(source)
        locus = by_locus.setdefault(gene, {})
        if copy in locus:
            raise PortableEvidenceError("duplicate HLA*LA locus-copy slot")
        locus[copy] = (allele, metrics, source)
    rows: list[dict[str, Any]] = []
    reconciled_slots = 0
    for gene, copies in sorted(by_locus.items()):
        if set(copies) != {1, 2}:
            raise PortableEvidenceError("HLA*LA locus lacks one diploid copy slot")
        allele1, metrics1, source1 = copies[1]
        allele2, metrics2, source2 = copies[2]
        _require_shared_locus_metrics(gene, metrics1, metrics2)
        reconciled_slots += 2
        alleles = (allele1, allele2)
        allele_count = sum(allele is not None for allele in alleles)
        call_status = {0: "no_call", 1: "partial_call", 2: "diploid_call"}[allele_count]
        filter_reasons: list[str] = []
        for copy, (allele, metrics) in enumerate(
            ((allele1, metrics1), (allele2, metrics2)), start=1
        ):
            if allele is None:
                filter_reasons.append(f"ALLELE{copy}:UNRESOLVED")
            if metrics["q1"] < 1:
                filter_reasons.append(f"ALLELE{copy}:Q1_BELOW_ONE")
            if not metrics["perfect_g"]:
                filter_reasons.append(f"ALLELE{copy}:PERFECT_G_FALSE")
            if metrics["proportion_kmers_covered"] == -1:
                filter_reasons.append(f"ALLELE{copy}:KMER_COVERAGE_NOT_ASSESSED")
        if call_status == "no_call":
            confidence = "UNRESOLVED"
        elif call_status == "partial_call" or filter_reasons:
            confidence = "REVIEW"
        else:
            confidence = "SUPPORTED"
        has_g_group = [allele is not None and allele.endswith("G") for allele in alleles]
        if all(has_g_group):
            allele_resolution = "g-group"
        elif any(has_g_group):
            allele_resolution = "mixed-g-group-and-full"
        else:
            allele_resolution = "full-or-expression-status"
        ambiguous = call_status != "diploid_call" or bool(filter_reasons)
        limitations = [
            "short-read HLA typing is research evidence and requires confirmation for medical use",
            "caller copy slots are not a phase-resolved chromosome assignment",
            "the pinned HLA*LA graph is a historical 2017 audit resource",
            "a confident historical graph call cannot override a newer nomenclature call",
        ]
        locus = f"HLA-{gene}"
        row: dict[str, Any] = {
            "locus": locus,
            "call_status": call_status,
            "allele1": allele1,
            "allele2": allele2,
            "allele1_g_group": allele1 if allele1 and allele1.endswith("G") else None,
            "allele2_g_group": allele2 if allele2 and allele2.endswith("G") else None,
            "allele1_two_field": two_field_allele(allele1),
            "allele2_two_field": two_field_allele(allele2),
            "allele_resolution": allele_resolution,
            "phase_resolved": False,
            "ambiguous": ambiguous,
            "primary_locus": gene in PRIMARY_HLA_LOCI,
            "allele1_candidates": [allele1] if allele1 is not None else [],
            "allele2_candidates": [allele2] if allele2 is not None else [],
            "allele1_scored_candidates": [],
            "allele2_scored_candidates": [],
            "allele1_quality_score": None,
            "allele2_quality_score": None,
            "allele1_hlala_metrics": metrics1,
            "allele2_hlala_metrics": metrics2,
            "locus_depth": metrics1["average_coverage"],
            "caller_database_name": provenance.caller_database_name,
            "caller_database_version": provenance.caller_database_version,
            "caller_database_identity_sha256": provenance.caller_database_identity_sha256,
            "nomenclature_database_name": provenance.nomenclature_database_name,
            "nomenclature_version": provenance.nomenclature_version,
            "nomenclature_identity_sha256": provenance.nomenclature_identity_sha256,
            "database_age_limitation": provenance.database_age_limitation,
            "evidence_id": "pending",
            "event_id": f"hla:{locus}",
            "assembly": "GRCh38",
            "reference_sha256": provenance.reference_sha256,
            "technical_source_role": provenance.technical_source_role,
            "caller_name": provenance.caller_name,
            "caller_version": provenance.caller_version,
            "caller_run_identity_sha256": provenance.caller_run_identity_sha256,
            "resource_manifest_sha256": provenance.resource_manifest_sha256,
            "source_alignment_sha256": provenance.source_alignment_sha256,
            "raw_artifact_id": provenance.raw_artifact_id,
            "raw_artifact_sha256": raw_artifact_sha256,
            "source_record_digest": semantic_record_digest(
                {"locus": locus, "copy1": source1, "copy2": source2}
            ),
            "filter_reasons": sorted(filter_reasons),
            "confidence_tier": confidence,
            "limitations": limitations,
            "confirmation_required": True,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
        row["evidence_id"] = make_evidence_id("hla_typing", row)
        rows.append(row)
    validated = validate_rows(
        "hla_typing",
        rows,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    if reconciled_slots != len(source_rows):
        raise PortableEvidenceError("HLA*LA source-to-locus reconciliation failed")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "source_allele_rows": len(source_rows),
        "reconciled_allele_slots": reconciled_slots,
        "typed_loci": len(validated),
        "primary_loci": sum(row["primary_locus"] for row in validated),
        "additional_loci": sum(not row["primary_locus"] for row in validated),
        "supported_loci": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_loci": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "unresolved_loci": sum(row["confidence_tier"] == "UNRESOLVED" for row in validated),
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": raw_artifact_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def materialize_hla_la_bestguess(
    input_table: str | Path,
    output_parquet: str | Path,
    provenance: HlaLaProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_hla_la_bestguess(
        input_table, provenance, forbidden_bindings=forbidden_bindings
    )
    output = write_extension(
        "hla_typing",
        rows,
        output_parquet,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    return {"status": "PASS", "parsing": parsing, "output": output}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-table", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--technical-source-role", required=True)
    parser.add_argument("--caller-run-identity-sha256", required=True)
    parser.add_argument("--resource-manifest-sha256", required=True)
    parser.add_argument("--source-alignment-sha256", required=True)
    parser.add_argument("--raw-artifact-id", required=True)
    parser.add_argument("--caller-database-name", required=True)
    parser.add_argument("--caller-database-version", required=True)
    parser.add_argument("--caller-database-identity-sha256", required=True)
    parser.add_argument("--nomenclature-database-name", required=True)
    parser.add_argument("--nomenclature-version", required=True)
    parser.add_argument("--nomenclature-identity-sha256", required=True)
    args = parser.parse_args()
    provenance = HlaLaProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        caller_database_name=args.caller_database_name,
        caller_database_version=args.caller_database_version,
        caller_database_identity_sha256=args.caller_database_identity_sha256,
        nomenclature_database_name=args.nomenclature_database_name,
        nomenclature_version=args.nomenclature_version,
        nomenclature_identity_sha256=args.nomenclature_identity_sha256,
    )
    result = materialize_hla_la_bestguess(args.input_table, args.output_parquet, provenance)
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
