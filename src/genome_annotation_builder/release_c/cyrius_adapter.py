"""Output-only adapter for separately installed, licence-gated Cyrius."""

from __future__ import annotations

import argparse
import csv
import json
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

CYRIUS_FIELDS = (
    "Coverage_MAD",
    "Median_depth",
    "Total_CN",
    "Spacer_CN",
    "Total_CN_raw",
    "Spacer_CN_raw",
    "Variants_called",
    "CNV_group",
    "Genotype",
    "Filter",
    "Raw_star_allele",
    "Call_info",
    "Exon9_CN",
    "CNV_consensus",
    "d67_snp_call",
    "d67_snp_raw",
    "Variant_raw_count",
)
TSV_COLUMNS = ["Sample", "Genotype", "Filter"]
RAW_COUNT = re.compile(
    r"^(?P<alt>[0-9]+)(?:\((?P<forward>[0-9]+):(?P<reverse>[0-9]+)\))?,(?P<ref>[0-9]+)$"
)
STAR_ALLELE = re.compile(r"\*[0-9A-Za-z.]+(?:x[0-9]+)?")
FILTERS = {
    "PASS",
    "More_than_one_possible_genotype",
    "Not_assigned_to_haplotypes",
    "LowQ_high_CN",
}


@dataclass(frozen=True)
class CyriusProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    definition_version: str
    definition_identity_sha256: str
    licence_gated_caller: bool = True
    caller_name: str = "Cyrius"
    caller_version: str = "1.1.1"


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PortableEvidenceError(f"duplicate Cyrius JSON key: {key}")
        result[key] = value
    return result


def _load_json(path: str | Path) -> tuple[str, dict[str, Any]]:
    try:
        with Path(path).open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_reject_duplicate_keys)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError("invalid Cyrius JSON artifact") from exc
    if not isinstance(payload, dict) or len(payload) != 1:
        raise PortableEvidenceError("Cyrius JSON must contain exactly one sample")
    sample_identifier, result = next(iter(payload.items()))
    if not isinstance(sample_identifier, str) or not sample_identifier:
        raise PortableEvidenceError("Cyrius JSON sample identifier is empty")
    if not isinstance(result, dict) or set(result) != set(CYRIUS_FIELDS):
        raise PortableEvidenceError("Cyrius JSON result schema differs")
    return sample_identifier, result


def _load_tsv(path: str | Path) -> tuple[str, str | None, str | None]:
    try:
        with Path(path).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != TSV_COLUMNS:
                raise PortableEvidenceError("Cyrius TSV column schema differs")
            rows = list(reader)
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("invalid Cyrius TSV artifact") from exc
    except csv.Error as exc:
        raise PortableEvidenceError("malformed Cyrius TSV artifact") from exc
    if len(rows) != 1 or None in rows[0]:
        raise PortableEvidenceError("Cyrius TSV must contain exactly one complete sample row")
    row = rows[0]
    sample_identifier = row["Sample"].strip()
    if not sample_identifier:
        raise PortableEvidenceError("Cyrius TSV sample identifier is empty")
    genotype = None if row["Genotype"] == "None" else row["Genotype"]
    caller_filter = None if row["Filter"] == "None" else row["Filter"]
    return sample_identifier, genotype, caller_filter


def _optional_float(value: Any, label: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise PortableEvidenceError(f"invalid Cyrius {label}")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid Cyrius {label}") from exc
    if result != result or result in {float("inf"), float("-inf")} or result < 0:
        raise PortableEvidenceError(f"invalid Cyrius {label}")
    return result


def _string_list(value: Any, label: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
        raise PortableEvidenceError(f"invalid Cyrius {label}")
    if len(value) != len(set(value)):
        raise PortableEvidenceError(f"duplicate Cyrius {label}")
    return list(value)


def _numeric_list(value: Any, label: str) -> list[float]:
    if value is None:
        return []
    if isinstance(value, str):
        if value == "":
            return []
        raw_values = value.split(",")
    elif isinstance(value, list):
        raw_values = value
    else:
        raise PortableEvidenceError(f"invalid Cyrius {label}")
    parsed = [_optional_float(item, label) for item in raw_values]
    if any(item is None for item in parsed):
        raise PortableEvidenceError(f"invalid Cyrius {label}")
    return [float(item) for item in parsed if item is not None]


def _raw_variant_counts(value: Any) -> list[dict[str, Any]]:
    if value is None:
        return []
    if not isinstance(value, dict):
        raise PortableEvidenceError("invalid Cyrius Variant_raw_count")
    rows: list[dict[str, Any]] = []
    for variant_id, raw in value.items():
        if not isinstance(variant_id, str) or not variant_id or not isinstance(raw, str):
            raise PortableEvidenceError("invalid Cyrius raw variant-count entry")
        match = RAW_COUNT.fullmatch(raw)
        if match is None:
            raise PortableEvidenceError("invalid Cyrius raw variant-count value")
        alt = int(match.group("alt"))
        ref = int(match.group("ref"))
        forward = int(match.group("forward")) if match.group("forward") else None
        reverse = int(match.group("reverse")) if match.group("reverse") else None
        if forward is not None and forward + int(reverse) != alt:
            raise PortableEvidenceError("Cyrius raw ALT strand counts disagree")
        rows.append(
            {
                "variant_id": variant_id,
                "alt_reads": alt,
                "ref_reads": ref,
                "forward_alt_reads": forward,
                "reverse_alt_reads": reverse,
            }
        )
    return sorted(rows, key=lambda row: row["variant_id"])


def _call_state(
    genotype: str | None, caller_filter: str | None
) -> tuple[str, str | None, list[str], bool, bool]:
    if genotype is None and caller_filter is None:
        return "no_call", None, [], False, False
    if genotype is None or caller_filter not in FILTERS:
        raise PortableEvidenceError("Cyrius genotype/filter state is inconsistent")
    if caller_filter == "PASS":
        if ";" in genotype or "/" not in genotype:
            raise PortableEvidenceError("Cyrius PASS genotype is not one assigned diplotype")
        return "called", genotype, [genotype], False, True
    if caller_filter == "More_than_one_possible_genotype":
        candidates = genotype.split(";")
        if len(candidates) < 2 or any("/" not in item for item in candidates):
            raise PortableEvidenceError("Cyrius multi-genotype ambiguity is malformed")
        return "ambiguous", None, candidates, True, False
    if caller_filter == "Not_assigned_to_haplotypes":
        if "/" in genotype or "_" not in genotype:
            raise PortableEvidenceError("Cyrius unassigned-haplotype genotype is malformed")
        return "ambiguous", None, [genotype], True, False
    if "/" not in genotype:
        raise PortableEvidenceError("Cyrius high-copy-number genotype is malformed")
    return "ambiguous", genotype, [genotype], True, False


def parse_cyrius_results(
    input_json: str | Path,
    input_tsv: str | Path,
    provenance: CyriusProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not provenance.licence_gated_caller:
        raise PortableEvidenceError("Cyrius must retain its external licence-gated boundary")
    json_sample, result = _load_json(input_json)
    tsv_sample, tsv_genotype, tsv_filter = _load_tsv(input_tsv)
    if json_sample != tsv_sample:
        raise PortableEvidenceError("Cyrius JSON and TSV sample identifiers differ")
    genotype = result["Genotype"]
    caller_filter = result["Filter"]
    if genotype is not None and not isinstance(genotype, str):
        raise PortableEvidenceError("invalid Cyrius Genotype")
    if caller_filter is not None and not isinstance(caller_filter, str):
        raise PortableEvidenceError("invalid Cyrius Filter")
    if genotype != tsv_genotype or caller_filter != tsv_filter:
        raise PortableEvidenceError("Cyrius JSON and TSV genotype/filter differ")
    call_status, diplotype, candidates, ambiguous, outside_eligible = _call_state(
        genotype, caller_filter
    )
    variants_called = _string_list(result["Variants_called"], "Variants_called")
    raw_counts = _raw_variant_counts(result["Variant_raw_count"])
    star_alleles: list[str] = []
    for candidate in candidates:
        for allele in STAR_ALLELE.findall(candidate):
            if allele not in star_alleles:
                star_alleles.append(allele)
    cnv_group = result["CNV_group"]
    if cnv_group is not None and not isinstance(cnv_group, str):
        raise PortableEvidenceError("invalid Cyrius CNV_group")
    raw_star = result["Raw_star_allele"]
    call_info = result["Call_info"]
    for value, label in ((raw_star, "Raw_star_allele"), (call_info, "Call_info")):
        if value is not None and not isinstance(value, str):
            raise PortableEvidenceError(f"invalid Cyrius {label}")
    artifact_set_sha256, component_sha256 = artifact_set_identity(
        {"json_result": input_json, "tsv_result": input_tsv}
    )
    structural_configurations = []
    if cnv_group:
        structural_configurations.append(f"CNV_GROUP:{cnv_group}")
    if genotype and ("+" in genotype or "x" in genotype):
        structural_configurations.append(f"GENOTYPE_STRUCTURE:{genotype}")
    filter_reasons = [] if caller_filter == "PASS" else [f"CALLER:{caller_filter or 'NO_CALL'}"]
    confidence = (
        "SUPPORTED"
        if call_status == "called"
        else "UNRESOLVED"
        if call_status == "no_call"
        else "REVIEW"
    )
    semantic_source = {field: result[field] for field in CYRIUS_FIELDS}
    row: dict[str, Any] = {
        "gene": "CYP2D6",
        "call_status": call_status,
        "diplotype": diplotype,
        "candidate_diplotypes": candidates,
        "star_alleles": star_alleles,
        "candidate_solution_details": [],
        "allele_decomposition": [],
        "gene_copy_number": None,
        "locus_total_copy_number": _optional_float(result["Total_CN"], "Total_CN"),
        "paralog_marker_copy_number": _optional_float(result["Spacer_CN"], "Spacer_CN"),
        "raw_normalized_locus_depth": _optional_float(result["Total_CN_raw"], "Total_CN_raw"),
        "raw_normalized_paralog_depth": _optional_float(result["Spacer_CN_raw"], "Spacer_CN_raw"),
        "structural_configurations": structural_configurations,
        "activity_score": None,
        "phenotype": None,
        "coverage_mad": _optional_float(result["Coverage_MAD"], "Coverage_MAD"),
        "median_depth": _optional_float(result["Median_depth"], "Median_depth"),
        "variants_called": variants_called,
        "raw_star_allele": raw_star,
        "caller_filter": caller_filter or "NO_CALL",
        "call_info": call_info,
        "exon_specific_copy_number": _optional_float(result["Exon9_CN"], "Exon9_CN"),
        "copy_number_consensus": _numeric_list(result["CNV_consensus"], "CNV_consensus"),
        "differentiating_site_copy_numbers": _numeric_list(result["d67_snp_call"], "d67_snp_call"),
        "differentiating_site_raw_depths": _numeric_list(result["d67_snp_raw"], "d67_snp_raw"),
        "raw_variant_counts": raw_counts,
        "definition_version": provenance.definition_version,
        "definition_identity_sha256": provenance.definition_identity_sha256,
        "licence_gated_caller": provenance.licence_gated_caller,
        "ambiguous": ambiguous,
        "outside_call_eligible": outside_eligible,
        "evidence_id": "pending",
        "event_id": "complex-pgx:CYP2D6",
        "assembly": "GRCh38",
        "reference_sha256": provenance.reference_sha256,
        "technical_source_role": provenance.technical_source_role,
        "caller_name": provenance.caller_name,
        "caller_version": provenance.caller_version,
        "caller_run_identity_sha256": provenance.caller_run_identity_sha256,
        "resource_manifest_sha256": provenance.resource_manifest_sha256,
        "source_alignment_sha256": provenance.source_alignment_sha256,
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": artifact_set_sha256,
        "source_record_digest": semantic_record_digest(semantic_source),
        "filter_reasons": filter_reasons,
        "confidence_tier": confidence,
        "limitations": [
            "complex CYP2D6 short-read evidence requires medical-grade confirmation",
            "Cyrius is a separately installed licence-gated caller and is not redistributed",
            "activity score and phenotype are not inferred by this output adapter",
            "Total_CN represents the combined CYP2D6+CYP2D7 locus, not CYP2D6 gene copy number",
            "Spacer_CN is retained only as a paralog-region marker copy number",
        ],
        "confirmation_required": True,
        "concordance_group_id": None,
        "concordance_status": "not_assessed",
    }
    row["evidence_id"] = make_evidence_id("complex_pharmacogenes", row)
    private_bindings = (*forbidden_bindings, json_sample)
    validated = validate_rows(
        "complex_pharmacogenes",
        [row],
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=private_bindings,
    )
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "json_samples": 1,
        "tsv_samples": 1,
        "typed_rows": 1,
        "call_status": call_status,
        "candidate_diplotypes": len(candidates),
        "raw_variant_count_rows": len(raw_counts),
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": artifact_set_sha256,
        "component_sha256": component_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
        "caller_executed_by_adapter": False,
        "licence_gated_caller": True,
    }
    if json_sample in json.dumps(summary, sort_keys=True):
        raise PortableEvidenceError("Cyrius summary persisted the sample identifier")
    return validated, summary


def materialize_cyrius_results(
    input_json: str | Path,
    input_tsv: str | Path,
    output_parquet: str | Path,
    provenance: CyriusProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_cyrius_results(
        input_json,
        input_tsv,
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
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--input-tsv", type=Path, required=True)
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
    provenance = CyriusProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        definition_version=args.definition_version,
        definition_identity_sha256=args.definition_identity_sha256,
    )
    result = materialize_cyrius_results(
        args.input_json,
        args.input_tsv,
        args.output_parquet,
        provenance,
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
