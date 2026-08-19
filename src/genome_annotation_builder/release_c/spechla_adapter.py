"""Strict, privacy-minimized SpecHLA result-set adapter."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .hla_common import (
    PRIMARY_HLA_LOCI,
    artifact_set_identity,
    finite_float,
    parse_hla_allele,
    require_g_group,
    semantic_record_digest,
    two_field_allele,
)
from .portable_model import (
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    validate_rows,
    write_extension,
)

VERSION_LINE = re.compile(r"^# version: IPD-IMGT/HLA (?P<version>\S+)$")
DETAIL_SLOT = re.compile(r"^HLA_(?P<gene>[A-Z0-9]+)_(?P<copy>[12])$")
RESULT_HEADER = [
    "Sample",
    *(f"HLA_{gene}_{copy}" for gene in PRIMARY_HLA_LOCI for copy in (1, 2)),
]
DETAIL_HEADER = ["Gene", "G_best", "allele", "details:allele;Score;Caucasian;Black;Asian"]


@dataclass(frozen=True)
class SpecHLAProvenance:
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
    database_age_limitation: bool = False
    caller_name: str = "SpecHLA"
    caller_version: str = "1.0.12"


def _read_lines(path: str | Path, label: str) -> list[str]:
    try:
        text = Path(path).read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError(f"cannot read SpecHLA {label}") from exc
    if "\x00" in text:
        raise PortableEvidenceError(f"SpecHLA {label} contains a NUL byte")
    lines = text.splitlines()
    while lines and not lines[-1]:
        lines.pop()
    if any(not line for line in lines):
        raise PortableEvidenceError(f"SpecHLA {label} contains an internal blank line")
    return lines


def _parse_version(line: str, expected: str) -> str:
    match = VERSION_LINE.fullmatch(line)
    if match is None:
        raise PortableEvidenceError("invalid SpecHLA IPD-IMGT/HLA version line")
    observed = match.group("version")
    if observed != expected:
        raise PortableEvidenceError("SpecHLA nomenclature version differs from provenance")
    return observed


def _trim_writer_trailer(cells: list[str], expected: int) -> list[str]:
    if len(cells) == expected + 1 and cells[-1] == "":
        return cells[:-1]
    return cells


def _read_result_table(
    path: str | Path, *, expected_version: str, g_group: bool
) -> tuple[str, dict[str, str]]:
    label = "G-group result" if g_group else "full result"
    lines = _read_lines(path, label)
    if len(lines) != 3:
        raise PortableEvidenceError(f"SpecHLA {label} must contain one sample row")
    _parse_version(lines[0], expected_version)
    header = lines[1].split("\t")
    if header != RESULT_HEADER:
        raise PortableEvidenceError(f"SpecHLA {label} column schema differs")
    cells = _trim_writer_trailer(lines[2].split("\t"), len(header))
    if len(cells) != len(header):
        raise PortableEvidenceError(f"SpecHLA {label} row width differs")
    sample_identifier = cells[0].strip()
    if not sample_identifier:
        raise PortableEvidenceError(f"SpecHLA {label} lacks a sample identifier")
    values = dict(zip(header[1:], cells[1:], strict=True))
    for slot, raw in values.items():
        match = DETAIL_SLOT.fullmatch(slot)
        assert match is not None
        parsed = parse_hla_allele(raw, expected_gene=match.group("gene"))
        if g_group:
            require_g_group(parsed)
    return sample_identifier, values


def _parse_candidates(value: str, gene: str) -> list[str]:
    result: list[str] = []
    for raw in value.split(";"):
        parsed = parse_hla_allele(raw, expected_gene=gene)
        if parsed is not None:
            result.append(parsed)
    if len(result) != len(set(result)):
        raise PortableEvidenceError("SpecHLA detail candidates contain duplicates")
    return result


def _parse_scored_candidate(value: str, gene: str) -> dict[str, Any]:
    fields = value.split(";")
    if len(fields) != 5:
        raise PortableEvidenceError("SpecHLA scored candidate width differs")
    allele = parse_hla_allele(fields[0], expected_gene=gene)
    if allele is None:
        raise PortableEvidenceError("SpecHLA scored candidate has no allele")
    score = finite_float(fields[1], "SpecHLA score")
    frequencies = [
        finite_float(fields[index], "SpecHLA population frequency") for index in range(2, 5)
    ]
    if not 0 <= score <= 100 or any(not 0 <= value <= 1 for value in frequencies):
        raise PortableEvidenceError("SpecHLA score or population frequency is out of bounds")
    return {
        "allele": allele,
        "score": score,
        "caucasian_frequency": frequencies[0],
        "black_frequency": frequencies[1],
        "asian_frequency": frequencies[2],
    }


def _read_details(path: str | Path, *, expected_version: str) -> dict[str, dict[str, Any]]:
    lines = _read_lines(path, "details result")
    if len(lines) < 3:
        raise PortableEvidenceError("SpecHLA details result is incomplete")
    _parse_version(lines[0], expected_version)
    if lines[1].split("\t") != DETAIL_HEADER:
        raise PortableEvidenceError("SpecHLA details column schema differs")
    expected_slots = set(RESULT_HEADER[1:])
    details: dict[str, dict[str, Any]] = {}
    for line in lines[2:]:
        cells = line.split("\t")
        while len(cells) > 4 and cells[-1] == "":
            cells.pop()
        if len(cells) < 4:
            raise PortableEvidenceError("SpecHLA details row is too narrow")
        slot = cells[0]
        match = DETAIL_SLOT.fullmatch(slot)
        if match is None or slot not in expected_slots:
            raise PortableEvidenceError("SpecHLA details contains an unexpected allele slot")
        if slot in details:
            raise PortableEvidenceError("SpecHLA details contains a duplicate allele slot")
        gene = match.group("gene")
        selected = parse_hla_allele(cells[1], expected_gene=gene)
        candidates = _parse_candidates(cells[2], gene)
        scored = [_parse_scored_candidate(value, gene) for value in cells[3:]]
        if selected is not None and selected not in candidates:
            raise PortableEvidenceError("SpecHLA selected allele is absent from detail candidates")
        if selected is not None and not scored:
            raise PortableEvidenceError("SpecHLA selected allele lacks score evidence")
        details[slot] = {
            "selected": selected,
            "candidates": candidates,
            "scored": scored,
        }
    missing = sorted(expected_slots - set(details))
    if missing:
        raise PortableEvidenceError(f"SpecHLA details lacks allele slots: {missing}")
    if len(details) != len(expected_slots):
        raise PortableEvidenceError("SpecHLA detail cardinality differs")
    return details


def _quality_for_selected(selected: str | None, scored: list[dict[str, Any]]) -> float | None:
    if selected is None:
        return None
    selected_two_field = two_field_allele(selected)
    matching = [
        float(candidate["score"])
        for candidate in scored
        if two_field_allele(str(candidate["allele"])) == selected_two_field
    ]
    if not matching:
        raise PortableEvidenceError("SpecHLA selected allele lacks a matching scored group")
    return max(matching)


def parse_spechla_results(
    full_result: str | Path,
    g_group_result: str | Path,
    details_result: str | Path,
    provenance: SpecHLAProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    full_sample, full = _read_result_table(
        full_result, expected_version=provenance.nomenclature_version, g_group=False
    )
    group_sample, groups = _read_result_table(
        g_group_result, expected_version=provenance.nomenclature_version, g_group=True
    )
    if full_sample != group_sample:
        raise PortableEvidenceError("SpecHLA full and G-group sample identifiers differ")
    details = _read_details(details_result, expected_version=provenance.nomenclature_version)
    raw_artifact_sha256, component_sha256 = artifact_set_identity(
        {
            "full_result": full_result,
            "g_group_result": g_group_result,
            "details_result": details_result,
        }
    )
    rows: list[dict[str, Any]] = []
    slot_reconciliation = 0
    for gene in PRIMARY_HLA_LOCI:
        locus = f"HLA-{gene}"
        alleles: list[str | None] = []
        g_groups: list[str | None] = []
        candidates_by_slot: list[list[str]] = []
        scored_by_slot: list[list[dict[str, Any]]] = []
        qualities: list[float | None] = []
        filter_reasons: list[str] = []
        for copy in (1, 2):
            slot = f"HLA_{gene}_{copy}"
            allele = parse_hla_allele(full[slot], expected_gene=gene)
            g_group = require_g_group(parse_hla_allele(groups[slot], expected_gene=gene))
            detail = details[slot]
            if allele != detail["selected"]:
                raise PortableEvidenceError("SpecHLA full result and detail selection disagree")
            if (allele is None) != (g_group is None):
                filter_reasons.append(f"ALLELE{copy}:FULL_G_GROUP_CARDINALITY_MISMATCH")
            quality = _quality_for_selected(allele, detail["scored"])
            if quality is not None and quality < 95:
                filter_reasons.append(f"ALLELE{copy}:LOW_CALLER_SCORE")
            if allele is None:
                filter_reasons.append(f"ALLELE{copy}:UNRESOLVED")
            alleles.append(allele)
            g_groups.append(g_group)
            candidates_by_slot.append(detail["candidates"])
            scored_by_slot.append(detail["scored"])
            qualities.append(quality)
            slot_reconciliation += 1
        allele_count = sum(allele is not None for allele in alleles)
        call_status = {0: "no_call", 1: "partial_call", 2: "diploid_call"}[allele_count]
        ambiguous = call_status != "diploid_call" or any(
            len(candidates) > 1 for candidates in candidates_by_slot
        )
        if call_status == "no_call":
            confidence = "UNRESOLVED"
        elif call_status == "partial_call" or filter_reasons:
            confidence = "REVIEW"
        else:
            confidence = "SUPPORTED"
        limitations = [
            "short-read HLA typing is research evidence and requires confirmation for medical use",
            "caller allele slots are not a phase-resolved chromosome assignment",
            "mapping reference and allele nomenclature are separately digest-bound resources",
        ]
        if ambiguous:
            limitations.append("one or more higher-field allele assignments remain ambiguous")
        if provenance.database_age_limitation:
            limitations.append("the caller database has an explicit age limitation")
        semantic_source = {
            "locus": locus,
            "full": alleles,
            "g_group": g_groups,
            "details": [
                {
                    "candidates": candidates_by_slot[index],
                    "scored": scored_by_slot[index],
                }
                for index in range(2)
            ],
        }
        row: dict[str, Any] = {
            "locus": locus,
            "call_status": call_status,
            "allele1": alleles[0],
            "allele2": alleles[1],
            "allele1_g_group": g_groups[0],
            "allele2_g_group": g_groups[1],
            "allele1_two_field": two_field_allele(alleles[0]),
            "allele2_two_field": two_field_allele(alleles[1]),
            "allele_resolution": "full-sequence-with-g-group-projection",
            "phase_resolved": False,
            "ambiguous": ambiguous,
            "primary_locus": True,
            "allele1_candidates": candidates_by_slot[0],
            "allele2_candidates": candidates_by_slot[1],
            "allele1_scored_candidates": scored_by_slot[0],
            "allele2_scored_candidates": scored_by_slot[1],
            "allele1_quality_score": qualities[0],
            "allele2_quality_score": qualities[1],
            "allele1_hlala_metrics": None,
            "allele2_hlala_metrics": None,
            "locus_depth": None,
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
            "source_record_digest": semantic_record_digest(semantic_source),
            "filter_reasons": sorted(filter_reasons),
            "confidence_tier": confidence,
            "limitations": limitations,
            "confirmation_required": True,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
        row["evidence_id"] = make_evidence_id("hla_typing", row)
        rows.append(row)
    private_bindings = (*forbidden_bindings, full_sample)
    validated = validate_rows(
        "hla_typing",
        rows,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=private_bindings,
    )
    expected_slots = len(PRIMARY_HLA_LOCI) * 2
    if slot_reconciliation != expected_slots:
        raise PortableEvidenceError("SpecHLA allele-slot reconciliation failed")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "source_artifacts": 3,
        "full_result_slots": len(full),
        "g_group_result_slots": len(groups),
        "detail_result_slots": len(details),
        "reconciled_slots": slot_reconciliation,
        "typed_loci": len(validated),
        "supported_loci": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_loci": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "unresolved_loci": sum(row["confidence_tier"] == "UNRESOLVED" for row in validated),
        "ambiguous_loci": sum(row["ambiguous"] for row in validated),
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": raw_artifact_sha256,
        "component_sha256": component_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    rendered = json.dumps(summary, sort_keys=True)
    if full_sample in rendered:
        raise PortableEvidenceError("SpecHLA summary persisted the sample identifier")
    return validated, summary


def materialize_spechla_results(
    full_result: str | Path,
    g_group_result: str | Path,
    details_result: str | Path,
    output_parquet: str | Path,
    provenance: SpecHLAProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_spechla_results(
        full_result,
        g_group_result,
        details_result,
        provenance,
        forbidden_bindings=forbidden_bindings,
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
    parser.add_argument("--full-result", type=Path, required=True)
    parser.add_argument("--g-group-result", type=Path, required=True)
    parser.add_argument("--details-result", type=Path, required=True)
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
    parser.add_argument("--database-age-limitation", action="store_true")
    args = parser.parse_args()
    provenance = SpecHLAProvenance(
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
        database_age_limitation=args.database_age_limitation,
    )
    result = materialize_spechla_results(
        args.full_result,
        args.g_group_result,
        args.details_result,
        args.output_parquet,
        provenance,
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
