"""Output-only adapter for a separately installed, licence-gated Aldy caller."""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .hla_common import artifact_set_identity, semantic_record_digest
from .portable_model import (
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    validate_rows,
    write_extension,
)

OUTPUT_COLUMNS = (
    "Sample",
    "Gene",
    "SolutionID",
    "Major",
    "Minor",
    "Copy",
    "Allele",
    "Location",
    "Type",
    "Coverage",
    "Effect",
    "dbSNP",
    "Code",
    "Status",
)
OUTPUT_HEADER = "#" + "\t".join(OUTPUT_COLUMNS)
SOLUTION_PREFIX = "#Solution "
STAR_ALLELE = re.compile(r"\*[0-9A-Za-z.]+(?:x[0-9]+)?")
GENE_NAME = re.compile(r"^[A-Z][A-Z0-9-]*$")


@dataclass(frozen=True)
class AldyProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    definition_version: str
    definition_identity_sha256: str
    licence_gated_caller: bool = True
    caller_name: str = "Aldy"
    caller_version: str = "4.8.3"


def _plain_text(value: str, label: str, *, optional: bool = False) -> str | None:
    if not value:
        if optional:
            return None
        raise PortableEvidenceError(f"Aldy {label} is empty")
    if value != value.strip() or "\x00" in value:
        raise PortableEvidenceError(f"invalid Aldy {label}")
    return value


def _nonnegative_int(value: str, label: str) -> int:
    if not value.isdigit():
        raise PortableEvidenceError(f"invalid Aldy {label}")
    return int(value)


def _nonnegative_float(value: str, label: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise PortableEvidenceError(f"invalid Aldy {label}") from exc
    if not math.isfinite(number) or number < 0:
        raise PortableEvidenceError(f"invalid Aldy {label}")
    return number


def _parse_solution_header(line: str) -> dict[str, Any]:
    if not line.startswith(SOLUTION_PREFIX):
        raise PortableEvidenceError("unexpected Aldy comment line")
    identifier_text, separator, payload = line[len(SOLUTION_PREFIX) :].partition(": ")
    if not separator or not identifier_text.isdigit() or int(identifier_text) <= 0:
        raise PortableEvidenceError("invalid Aldy solution header")
    expression, separator, metadata = payload.rpartition("; gene=")
    if not separator or not expression or not metadata:
        raise PortableEvidenceError("incomplete Aldy solution header")
    segments = metadata.split("; ")
    definition = _plain_text(segments[0], "definition")
    annotations: dict[str, str] = {}
    for segment in segments[1:]:
        key, equals, value = segment.partition("=")
        if not equals or key not in {"cpic", "cpic_score"} or key in annotations:
            raise PortableEvidenceError("unknown or duplicate Aldy solution annotation")
        annotations[key] = str(_plain_text(value, f"solution annotation {key}"))
    if "cpic_score" in annotations and "cpic" not in annotations:
        raise PortableEvidenceError("Aldy CPIC score lacks a phenotype")
    activity_score = (
        _nonnegative_float(annotations["cpic_score"], "CPIC score")
        if "cpic_score" in annotations
        else None
    )
    return {
        "solution_id": int(identifier_text),
        "solution_expression": expression,
        "definition": definition,
        "phenotype": annotations.get("cpic"),
        "activity_score": activity_score,
        "rows": [],
        "major_diplotype": None,
        "minor_alleles": None,
    }


def _parse_data_line(line: str) -> dict[str, Any]:
    fields = line.split("\t")
    if len(fields) == len(OUTPUT_COLUMNS) - 1:
        # Aldy 4.8.3 omits the trailing Status cell on mutation-bearing rows.
        fields.append("")
    if len(fields) != len(OUTPUT_COLUMNS):
        raise PortableEvidenceError("Aldy data-column count differs")
    raw = dict(zip(OUTPUT_COLUMNS, fields, strict=True))
    sample = _plain_text(raw["Sample"], "sample identifier")
    gene = _plain_text(raw["Gene"], "gene")
    if GENE_NAME.fullmatch(str(gene)) is None:
        raise PortableEvidenceError("invalid Aldy gene")
    solution_id = _nonnegative_int(raw["SolutionID"], "solution identifier")
    if solution_id <= 0:
        raise PortableEvidenceError("invalid Aldy solution identifier")
    major = _plain_text(raw["Major"], "major diplotype")
    minor = _plain_text(raw["Minor"], "minor alleles")
    copy_index = _nonnegative_int(raw["Copy"], "copy index")
    allele = _plain_text(raw["Allele"], "allele")
    location = _nonnegative_int(raw["Location"], "zero-based location") if raw["Location"] else None
    if location is None:
        if any(raw[column] for column in ("Type", "Coverage", "Effect", "dbSNP", "Code", "Status")):
            raise PortableEvidenceError("Aldy empty-location row contains mutation fields")
        nucleotide_change = coverage = effect = dbsnp = code = status = None
    else:
        nucleotide_change = _plain_text(raw["Type"], "mutation type")
        coverage = _nonnegative_float(raw["Coverage"], "mutation coverage")
        effect = _plain_text(raw["Effect"], "mutation effect", optional=True)
        dbsnp = _plain_text(raw["dbSNP"], "dbSNP", optional=True)
        code = _plain_text(raw["Code"], "code", optional=True)
        status = _plain_text(raw["Status"], "status", optional=True)
    return {
        "sample": sample,
        "gene": gene,
        "solution_id": solution_id,
        "major_diplotype": major,
        "minor_text": minor,
        "copy_index": copy_index,
        "allele": allele,
        "position_zero_based": location,
        "nucleotide_change": nucleotide_change,
        "coverage": coverage,
        "effect": effect,
        "dbsnp": dbsnp,
        "code": code,
        "status": status,
    }


def _load_output(path: str | Path) -> tuple[str, str, list[dict[str, Any]]]:
    try:
        text = Path(path).read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("invalid Aldy decomposition artifact") from exc
    lines = text.splitlines()
    if not lines or lines[0] != OUTPUT_HEADER or any(not line for line in lines):
        raise PortableEvidenceError("Aldy decomposition header or line structure differs")
    solutions: list[dict[str, Any]] = []
    active: dict[str, Any] | None = None
    sample_identifier: str | None = None
    gene: str | None = None
    for line in lines[1:]:
        if line.startswith("#"):
            if active is not None and not active["rows"]:
                raise PortableEvidenceError("Aldy solution has no decomposition rows")
            active = _parse_solution_header(line)
            solutions.append(active)
            continue
        if active is None:
            raise PortableEvidenceError("Aldy data row precedes its solution header")
        parsed = _parse_data_line(line)
        if parsed["solution_id"] != active["solution_id"]:
            raise PortableEvidenceError("Aldy row references the wrong solution")
        if sample_identifier is None:
            sample_identifier = str(parsed["sample"])
            gene = str(parsed["gene"])
        elif parsed["sample"] != sample_identifier or parsed["gene"] != gene:
            raise PortableEvidenceError("Aldy rows contain multiple samples or genes")
        if active["major_diplotype"] is None:
            active["major_diplotype"] = parsed["major_diplotype"]
            active["minor_alleles"] = str(parsed["minor_text"]).split(";")
            if any(not allele for allele in active["minor_alleles"]):
                raise PortableEvidenceError("Aldy minor-allele list is malformed")
        elif (
            parsed["major_diplotype"] != active["major_diplotype"]
            or str(parsed["minor_text"]).split(";") != active["minor_alleles"]
        ):
            raise PortableEvidenceError("Aldy solution rows disagree on diplotype")
        active["rows"].append(parsed)
    if active is not None and not active["rows"]:
        raise PortableEvidenceError("Aldy solution has no decomposition rows")
    if sample_identifier is None or gene is None or not solutions:
        raise PortableEvidenceError("Aldy decomposition contains no solutions")
    expected_ids = list(range(1, len(solutions) + 1))
    if [solution["solution_id"] for solution in solutions] != expected_ids:
        raise PortableEvidenceError("Aldy solution identifiers are not contiguous")
    definitions = {solution["definition"] for solution in solutions}
    if len(definitions) != 1:
        raise PortableEvidenceError("Aldy solutions use different definitions")
    definition = str(next(iter(definitions)))
    if not definition.startswith(f"{gene}-"):
        raise PortableEvidenceError("Aldy definition does not match its gene")
    for solution in solutions:
        minor_alleles = solution["minor_alleles"]
        copies = sorted({row["copy_index"] for row in solution["rows"]})
        if copies != list(range(len(minor_alleles))):
            raise PortableEvidenceError("Aldy solution copy indexes are incomplete")
        for row in solution["rows"]:
            if row["allele"] != minor_alleles[row["copy_index"]]:
                raise PortableEvidenceError("Aldy copy allele differs from the minor list")
    return sample_identifier, gene, solutions


def parse_aldy_results(
    input_aldy: str | Path,
    provenance: AldyProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not provenance.licence_gated_caller:
        raise PortableEvidenceError("Aldy must retain its external licence-gated boundary")
    sample_identifier, gene, solutions = _load_output(input_aldy)
    definitions = {str(solution["definition"]) for solution in solutions}
    if definitions != {provenance.definition_version}:
        raise PortableEvidenceError("Aldy output definition differs from frozen provenance")
    artifact_sha256, component_sha256 = artifact_set_identity({"aldy_decomposition": input_aldy})
    candidate_diplotypes: list[str] = []
    star_alleles: list[str] = []
    candidate_details: list[dict[str, Any]] = []
    decomposition: list[dict[str, Any]] = []
    copy_counts: set[int] = set()
    variants_called: set[str] = set()
    for solution in solutions:
        major = str(solution["major_diplotype"])
        if major not in candidate_diplotypes:
            candidate_diplotypes.append(major)
        for allele in STAR_ALLELE.findall(major):
            if allele not in star_alleles:
                star_alleles.append(allele)
        copy_count = len(solution["minor_alleles"])
        copy_counts.add(copy_count)
        candidate_details.append(
            {
                "solution_id": solution["solution_id"],
                "solution_expression": solution["solution_expression"],
                "major_diplotype": major,
                "minor_alleles": list(solution["minor_alleles"]),
                "copy_count": copy_count,
                "activity_score": solution["activity_score"],
                "phenotype": solution["phenotype"],
            }
        )
        for source in solution["rows"]:
            item = {
                "solution_id": source["solution_id"],
                "copy_index": source["copy_index"],
                "allele": source["allele"],
                "position_zero_based": source["position_zero_based"],
                "nucleotide_change": source["nucleotide_change"],
                "coverage": source["coverage"],
                "effect": source["effect"],
                "dbsnp": source["dbsnp"],
                "code": source["code"],
                "status": source["status"],
            }
            decomposition.append(item)
            if source["position_zero_based"] is not None:
                marker = source["dbsnp"] or (
                    f"position0:{source['position_zero_based']}:{source['nucleotide_change']}"
                )
                variants_called.add(str(marker))
    unique_solution = len(solutions) == 1
    selected = candidate_diplotypes[0] if unique_solution else None
    selected_solution = solutions[0] if unique_solution else None
    call_status = "called" if unique_solution else "ambiguous"
    caller_filter = "UNIQUE_OPTIMAL_SOLUTION" if unique_solution else "MULTIPLE_OPTIMAL_SOLUTIONS"
    semantic_source = {
        "gene": gene,
        "definition": provenance.definition_version,
        "candidate_solution_details": candidate_details,
        "allele_decomposition": decomposition,
    }
    row: dict[str, Any] = {
        "gene": gene,
        "call_status": call_status,
        "diplotype": selected,
        "candidate_diplotypes": candidate_diplotypes,
        "star_alleles": star_alleles,
        "candidate_solution_details": candidate_details,
        "allele_decomposition": decomposition,
        "gene_copy_number": float(next(iter(copy_counts))) if len(copy_counts) == 1 else None,
        "locus_total_copy_number": None,
        "paralog_marker_copy_number": None,
        "raw_normalized_locus_depth": None,
        "raw_normalized_paralog_depth": None,
        "structural_configurations": [
            f"ALDY_MAJOR:{candidate}" for candidate in candidate_diplotypes
        ],
        "activity_score": selected_solution["activity_score"] if unique_solution else None,
        "phenotype": selected_solution["phenotype"] if unique_solution else None,
        "coverage_mad": None,
        "median_depth": None,
        "variants_called": sorted(variants_called),
        "raw_star_allele": selected_solution["solution_expression"] if unique_solution else None,
        "caller_filter": caller_filter,
        "call_info": f"optimal_solution_count={len(solutions)}",
        "exon_specific_copy_number": None,
        "copy_number_consensus": [],
        "differentiating_site_copy_numbers": [],
        "differentiating_site_raw_depths": [],
        "raw_variant_counts": [],
        "definition_version": provenance.definition_version,
        "definition_identity_sha256": provenance.definition_identity_sha256,
        "licence_gated_caller": True,
        "ambiguous": not unique_solution,
        "outside_call_eligible": unique_solution,
        "evidence_id": "pending",
        "event_id": f"complex-pgx:{gene}",
        "assembly": "GRCh38",
        "reference_sha256": provenance.reference_sha256,
        "technical_source_role": provenance.technical_source_role,
        "caller_name": provenance.caller_name,
        "caller_version": provenance.caller_version,
        "caller_run_identity_sha256": provenance.caller_run_identity_sha256,
        "resource_manifest_sha256": provenance.resource_manifest_sha256,
        "source_alignment_sha256": provenance.source_alignment_sha256,
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": artifact_sha256,
        "source_record_digest": semantic_record_digest(semantic_source),
        "filter_reasons": [] if unique_solution else ["CALLER:MULTIPLE_OPTIMAL_SOLUTIONS"],
        "confidence_tier": "SUPPORTED" if unique_solution else "REVIEW",
        "limitations": [
            "complex pharmacogene short-read evidence requires medical-grade confirmation",
            "Aldy is a separately installed licence-gated caller and is not redistributed",
            "caller diplotype-side grouping is retained but is not independent read-backed phasing",
            "multiple optimal solutions are retained without automatic selection",
            "caller-derived activity is not a clinical prescribing recommendation",
        ],
        "confirmation_required": True,
        "concordance_group_id": None,
        "concordance_status": "not_assessed",
    }
    row["evidence_id"] = make_evidence_id("complex_pharmacogenes", row)
    private_bindings = (*forbidden_bindings, sample_identifier)
    validated = validate_rows(
        "complex_pharmacogenes",
        [row],
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=private_bindings,
    )
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "gene": gene,
        "solution_count": len(solutions),
        "candidate_diplotype_count": len(candidate_diplotypes),
        "allele_decomposition_rows": len(decomposition),
        "call_status": call_status,
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": artifact_sha256,
        "component_sha256": component_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
        "caller_executed_by_adapter": False,
        "licence_gated_caller": True,
    }
    if sample_identifier in json.dumps(summary, sort_keys=True):
        raise PortableEvidenceError("Aldy summary persisted the sample identifier")
    return validated, summary


def materialize_aldy_results(
    input_aldy: str | Path,
    output_parquet: str | Path,
    provenance: AldyProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_aldy_results(
        input_aldy,
        provenance,
        forbidden_bindings=forbidden_bindings,
    )
    output = write_extension(
        "complex_pharmacogenes",
        rows,
        output_parquet,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    return {"status": "PASS", "parsing": parsing, "output": output}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-aldy", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--technical-source-role", required=True)
    parser.add_argument("--caller-run-identity-sha256", required=True)
    parser.add_argument("--resource-manifest-sha256", required=True)
    parser.add_argument("--source-alignment-sha256", required=True)
    parser.add_argument("--raw-artifact-id", required=True)
    parser.add_argument("--definition-version", required=True)
    parser.add_argument("--definition-identity-sha256", required=True)
    args = parser.parse_args()
    provenance = AldyProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        definition_version=args.definition_version,
        definition_identity_sha256=args.definition_identity_sha256,
    )
    result = materialize_aldy_results(
        args.input_aldy,
        args.output_parquet,
        provenance,
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
