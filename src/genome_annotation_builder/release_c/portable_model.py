"""Typed, privacy-minimized portable evidence model for read-backed extensions.

This is a generic Release-C candidate implementation.  It deliberately models
one caller/source observation per row: cross-caller or replicate agreement is an
annotation, never a vote that destroys the source observation.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA_VERSION = 1
DOMAINS = (
    "structural_variants",
    "copy_number",
    "repeat_expansions",
    "hla_typing",
    "complex_pharmacogenes",
    "mitochondrial_heteroplasmy",
    "read_backed_evidence",
)
STANDARD_CHROMS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY", "chrM"])
CHROM_RANK = {chrom: index for index, chrom in enumerate(STANDARD_CHROMS)}
SHA256 = re.compile(r"^[0-9a-f]{64}$")
CONCRETE_ALLELE = re.compile(r"^[ACGT]+$")
IUPAC_DNA_MOTIF = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
FORBIDDEN_PORTABLE_TEXT = re.compile(
    r"(?:/(?:raid|home|mnt|media|Users)/|file://|s3://|gs://|"
    r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|"
    r"gh[pousr]_[A-Za-z0-9_]{20,}|AKIA[0-9A-Z]{16})"
)
PRIVATE_TOKEN_CHAR = re.compile(r"[A-Za-z0-9]")
SOURCE_ROLES = frozenset({"primary", "technical_replicate", "combined_evidence"})
CONFIDENCE_TIERS = frozenset({"SUPPORTED", "REVIEW", "UNRESOLVED"})
CONCORDANCE_STATES = frozenset(
    {
        "not_assessed",
        "comparable_concordant",
        "comparable_discordant",
        "not_comparable",
        "source_specific",
    }
)


class PortableEvidenceError(ValueError):
    """The typed evidence would violate a Release-C invariant."""


def contains_forbidden_binding(text: str, binding: str) -> bool:
    """Detect private bindings without treating a short ID as any substring."""

    if not binding:
        return False
    if len(binding) >= 16 or "/" in binding or "\\" in binding or SHA256.fullmatch(binding):
        return binding in text
    start = 0
    while True:
        index = text.find(binding, start)
        if index < 0:
            return False
        before = text[index - 1] if index else ""
        after_index = index + len(binding)
        after = text[after_index] if after_index < len(text) else ""
        if (not before or PRIVATE_TOKEN_CHAR.fullmatch(before) is None) and (
            not after or PRIVATE_TOKEN_CHAR.fullmatch(after) is None
        ):
            return True
        start = index + 1


def _field(name: str, data_type: pa.DataType, *, nullable: bool = True) -> pa.Field:
    return pa.field(name, data_type, nullable=nullable)


def _common_fields() -> list[pa.Field]:
    return [
        _field("evidence_id", pa.string(), nullable=False),
        _field("event_id", pa.string(), nullable=False),
        _field("assembly", pa.string(), nullable=False),
        _field("reference_sha256", pa.string(), nullable=False),
        _field("technical_source_role", pa.string(), nullable=False),
        _field("caller_name", pa.string(), nullable=False),
        _field("caller_version", pa.string(), nullable=False),
        _field("caller_run_identity_sha256", pa.string(), nullable=False),
        _field("resource_manifest_sha256", pa.string(), nullable=False),
        _field("source_alignment_sha256", pa.string(), nullable=False),
        _field("raw_artifact_id", pa.string(), nullable=False),
        _field("raw_artifact_sha256", pa.string(), nullable=False),
        _field("source_record_digest", pa.string(), nullable=False),
        _field("filter_reasons", pa.list_(pa.string()), nullable=False),
        _field("confidence_tier", pa.string(), nullable=False),
        _field("limitations", pa.list_(pa.string()), nullable=False),
        _field("confirmation_required", pa.bool_(), nullable=False),
        _field("concordance_group_id", pa.string()),
        _field("concordance_status", pa.string(), nullable=False),
    ]


def _schema(fields: list[pa.Field]) -> pa.Schema:
    return pa.schema(
        fields + _common_fields(),
        metadata={
            b"gab.extension.schema_version": str(SCHEMA_VERSION).encode(),
            b"gab.extension.subject_scope": b"bundle_subject_only",
        },
    )


def structural_variant_schema() -> pa.Schema:
    return _schema(
        [
            _field("caller_record_id", pa.string(), nullable=False),
            _field("chrom1", pa.string(), nullable=False),
            _field("pos1", pa.int32(), nullable=False),
            _field("chrom2", pa.string(), nullable=False),
            _field("pos2", pa.int32(), nullable=False),
            _field("sv_type", pa.string(), nullable=False),
            _field("orientation", pa.string()),
            _field("length", pa.int64()),
            _field("genotype", pa.list_(pa.int8())),
            _field("call_status", pa.string(), nullable=False),
            _field("genotype_quality", pa.float64()),
            _field("genotype_likelihoods", pa.list_(pa.int32()), nullable=False),
            _field("copy_number", pa.float64()),
            _field("quality", pa.float64()),
            _field("imprecise", pa.bool_(), nullable=False),
            _field("cipos", pa.list_(pa.int32()), nullable=False),
            _field("ciend", pa.list_(pa.int32()), nullable=False),
            _field("paired_end_support", pa.int32()),
            _field("caller_site_paired_end_support", pa.int32()),
            _field("paired_end_reference_reads", pa.int32()),
            _field("split_read_support", pa.int32()),
            _field("caller_site_split_read_support", pa.int32()),
            _field("split_read_reference_reads", pa.int32()),
            _field("read_depth_support", pa.int32()),
            _field("regional_depth", pa.int32()),
            _field("mapping_quality", pa.float64()),
            _field("split_read_mapping_quality", pa.float64()),
            _field("insertion_length", pa.int32()),
            _field("homology_length", pa.int32()),
            _field("bnd_mate_id", pa.string()),
            _field("mate_link_status", pa.string(), nullable=False),
        ]
    )


def copy_number_schema() -> pa.Schema:
    model_type = pa.struct(
        [
            _field("rank", pa.int8(), nullable=False),
            _field("total_copy_number", pa.int16(), nullable=False),
            _field("haplotype1_copy_number", pa.int16(), nullable=False),
            _field("haplotype2_copy_number", pa.int16(), nullable=False),
            _field("likelihood", pa.float64(), nullable=False),
            _field("cell_fraction", pa.float64(), nullable=False),
            _field("reported_in_tsv", pa.bool_(), nullable=False),
            _field("padding", pa.bool_(), nullable=False),
        ]
    )
    return _schema(
        [
            _field("call_type", pa.string(), nullable=False),
            _field("caller_model", pa.string(), nullable=False),
            _field("chrom", pa.string(), nullable=False),
            _field("start_pos", pa.int32(), nullable=False),
            _field("end_pos", pa.int32(), nullable=False),
            _field("copy_number", pa.float64()),
            _field("copy_number_basis", pa.string(), nullable=False),
            _field("read_depth_copy_number", pa.float64()),
            _field("major_copy_number", pa.float64()),
            _field("minor_copy_number", pa.float64()),
            _field("normalized_read_depth", pa.float64()),
            _field("log2_ratio", pa.float64()),
            _field("baf", pa.float64()),
            _field("baf_deviation_from_half", pa.float64()),
            _field("combined_p_value", pa.float64()),
            _field("caller_reported_rd_p_value", pa.float64()),
            _field("baf_p_value", pa.float64()),
            _field("caller_reported_rd_p_value_duplicates_combined", pa.bool_()),
            _field("likelihood_deletion", pa.float64()),
            _field("likelihood_copy_neutral_loh", pa.float64()),
            _field("likelihood_duplication", pa.float64()),
            _field("e_value_1", pa.float64()),
            _field("e_value_2", pa.float64()),
            _field("e_value_3", pa.float64()),
            _field("e_value_4", pa.float64()),
            _field("q0_fraction", pa.float64()),
            _field("reference_gap_fraction", pa.float64()),
            _field("segment_gap_fraction", pa.float64()),
            _field("strict_mask_fraction", pa.float64()),
            _field("distance_to_large_gap", pa.int64()),
            _field("rd_mask_applied", pa.bool_(), nullable=False),
            _field("baf_mask_applied", pa.bool_(), nullable=False),
            _field("bin_size", pa.int32()),
            _field("bins_observed", pa.int32()),
            _field("heterozygous_sites", pa.int32()),
            _field("homozygous_sites", pa.int32()),
            _field("copy_number_models", pa.list_(model_type), nullable=False),
            _field("source_baf_vcf_sha256", pa.string()),
            _field("sex_ploidy", pa.int8()),
        ]
    )


def repeat_expansion_schema() -> pa.Schema:
    return _schema(
        [
            _field("locus_id", pa.string(), nullable=False),
            _field("variant_id", pa.string(), nullable=False),
            _field("variant_subtype", pa.string(), nullable=False),
            _field("chrom", pa.string(), nullable=False),
            _field("start_pos", pa.int32(), nullable=False),
            _field("end_pos", pa.int32(), nullable=False),
            _field("motif", pa.string(), nullable=False),
            _field("ploidy", pa.int8(), nullable=False),
            _field("reference_repeat_count", pa.int32()),
            _field("allele1_repeat_count", pa.int32()),
            _field("allele1_ci_lower", pa.int32()),
            _field("allele1_ci_upper", pa.int32()),
            _field("allele1_count_is_lower_bound", pa.bool_(), nullable=False),
            _field("allele1_support_type", pa.string()),
            _field("allele1_spanning_reads", pa.int32()),
            _field("allele1_flanking_reads", pa.int32()),
            _field("allele1_in_repeat_reads", pa.int32()),
            _field("allele2_repeat_count", pa.int32()),
            _field("allele2_ci_lower", pa.int32()),
            _field("allele2_ci_upper", pa.int32()),
            _field("allele2_count_is_lower_bound", pa.bool_(), nullable=False),
            _field("allele2_support_type", pa.string()),
            _field("allele2_spanning_reads", pa.int32()),
            _field("allele2_flanking_reads", pa.int32()),
            _field("allele2_in_repeat_reads", pa.int32()),
            _field("locus_coverage", pa.float64()),
            _field("spanning_reads", pa.int32()),
            _field("flanking_reads", pa.int32()),
            _field("in_repeat_reads", pa.int32()),
            _field("off_target_evidence", pa.bool_()),
            _field("short_read_resolvable", pa.bool_(), nullable=False),
        ]
    )


def _hla_scored_candidate_type() -> pa.StructType:
    return pa.struct(
        [
            _field("allele", pa.string(), nullable=False),
            _field("score", pa.float64(), nullable=False),
            _field("caucasian_frequency", pa.float64()),
            _field("black_frequency", pa.float64()),
            _field("asian_frequency", pa.float64()),
        ]
    )


def _hla_la_allele_metrics_type() -> pa.StructType:
    return pa.struct(
        [
            _field("q1", pa.float64(), nullable=False),
            _field("q2", pa.float64(), nullable=False),
            _field("average_coverage", pa.float64(), nullable=False),
            _field("coverage_first_decile", pa.float64(), nullable=False),
            _field("minimum_coverage", pa.float64(), nullable=False),
            _field("proportion_kmers_covered", pa.float64(), nullable=False),
            _field("locus_avg_column_error", pa.float64(), nullable=False),
            _field("unaccounted_columns_fraction_gt_0_2", pa.float64(), nullable=False),
            _field("perfect_g", pa.bool_(), nullable=False),
        ]
    )


def hla_typing_schema() -> pa.Schema:
    return _schema(
        [
            _field("locus", pa.string(), nullable=False),
            _field("call_status", pa.string(), nullable=False),
            _field("allele1", pa.string()),
            _field("allele2", pa.string()),
            _field("allele1_g_group", pa.string()),
            _field("allele2_g_group", pa.string()),
            _field("allele1_two_field", pa.string()),
            _field("allele2_two_field", pa.string()),
            _field("allele_resolution", pa.string(), nullable=False),
            _field("phase_resolved", pa.bool_(), nullable=False),
            _field("ambiguous", pa.bool_(), nullable=False),
            _field("primary_locus", pa.bool_(), nullable=False),
            _field("allele1_candidates", pa.list_(pa.string()), nullable=False),
            _field("allele2_candidates", pa.list_(pa.string()), nullable=False),
            _field(
                "allele1_scored_candidates",
                pa.list_(_hla_scored_candidate_type()),
                nullable=False,
            ),
            _field(
                "allele2_scored_candidates",
                pa.list_(_hla_scored_candidate_type()),
                nullable=False,
            ),
            _field("allele1_quality_score", pa.float64()),
            _field("allele2_quality_score", pa.float64()),
            _field("allele1_hlala_metrics", _hla_la_allele_metrics_type()),
            _field("allele2_hlala_metrics", _hla_la_allele_metrics_type()),
            _field("locus_depth", pa.float64()),
            _field("caller_database_name", pa.string(), nullable=False),
            _field("caller_database_version", pa.string(), nullable=False),
            _field("caller_database_identity_sha256", pa.string(), nullable=False),
            _field("nomenclature_database_name", pa.string(), nullable=False),
            _field("nomenclature_version", pa.string(), nullable=False),
            _field("nomenclature_identity_sha256", pa.string(), nullable=False),
            _field("database_age_limitation", pa.bool_(), nullable=False),
        ]
    )


def _pgx_candidate_solution_type() -> pa.StructType:
    return pa.struct(
        [
            _field("solution_id", pa.int32(), nullable=False),
            _field("solution_expression", pa.string(), nullable=False),
            _field("major_diplotype", pa.string(), nullable=False),
            _field("minor_alleles", pa.list_(pa.string()), nullable=False),
            _field("copy_count", pa.int32(), nullable=False),
            _field("activity_score", pa.float64()),
            _field("phenotype", pa.string()),
        ]
    )


def _pgx_allele_decomposition_type() -> pa.StructType:
    return pa.struct(
        [
            _field("solution_id", pa.int32(), nullable=False),
            _field("copy_index", pa.int32(), nullable=False),
            _field("allele", pa.string(), nullable=False),
            _field("position_zero_based", pa.int64()),
            _field("nucleotide_change", pa.string()),
            _field("coverage", pa.float64()),
            _field("effect", pa.string()),
            _field("dbsnp", pa.string()),
            _field("code", pa.string()),
            _field("status", pa.string()),
        ]
    )


def complex_pharmacogene_schema() -> pa.Schema:
    return _schema(
        [
            _field("gene", pa.string(), nullable=False),
            _field("call_status", pa.string(), nullable=False),
            _field("diplotype", pa.string()),
            _field("candidate_diplotypes", pa.list_(pa.string()), nullable=False),
            _field("star_alleles", pa.list_(pa.string()), nullable=False),
            _field(
                "candidate_solution_details",
                pa.list_(_pgx_candidate_solution_type()),
                nullable=False,
            ),
            _field(
                "allele_decomposition",
                pa.list_(_pgx_allele_decomposition_type()),
                nullable=False,
            ),
            _field("gene_copy_number", pa.float64()),
            _field("locus_total_copy_number", pa.float64()),
            _field("paralog_marker_copy_number", pa.float64()),
            _field("raw_normalized_locus_depth", pa.float64()),
            _field("raw_normalized_paralog_depth", pa.float64()),
            _field("structural_configurations", pa.list_(pa.string()), nullable=False),
            _field("activity_score", pa.float64()),
            _field("phenotype", pa.string()),
            _field("coverage_mad", pa.float64()),
            _field("median_depth", pa.float64()),
            _field("variants_called", pa.list_(pa.string()), nullable=False),
            _field("raw_star_allele", pa.string()),
            _field("caller_filter", pa.string(), nullable=False),
            _field("call_info", pa.string()),
            _field("exon_specific_copy_number", pa.float64()),
            _field("copy_number_consensus", pa.list_(pa.float64()), nullable=False),
            _field(
                "differentiating_site_copy_numbers",
                pa.list_(pa.float64()),
                nullable=False,
            ),
            _field(
                "differentiating_site_raw_depths",
                pa.list_(pa.float64()),
                nullable=False,
            ),
            _field(
                "raw_variant_counts",
                pa.list_(
                    pa.struct(
                        [
                            _field("variant_id", pa.string(), nullable=False),
                            _field("alt_reads", pa.int32(), nullable=False),
                            _field("ref_reads", pa.int32(), nullable=False),
                            _field("forward_alt_reads", pa.int32()),
                            _field("reverse_alt_reads", pa.int32()),
                        ]
                    )
                ),
                nullable=False,
            ),
            _field("definition_version", pa.string(), nullable=False),
            _field("definition_identity_sha256", pa.string(), nullable=False),
            _field("licence_gated_caller", pa.bool_(), nullable=False),
            _field("ambiguous", pa.bool_(), nullable=False),
            _field("outside_call_eligible", pa.bool_(), nullable=False),
        ]
    )


def mitochondrial_heteroplasmy_schema() -> pa.Schema:
    return _schema(
        [
            _field("chrom", pa.string(), nullable=False),
            _field("pos", pa.int32(), nullable=False),
            _field("ref", pa.string(), nullable=False),
            _field("alt", pa.string(), nullable=False),
            _field("allele_fraction", pa.float64(), nullable=False),
            _field("depth", pa.int32(), nullable=False),
            _field("ref_depth", pa.int32()),
            _field("alt_depth", pa.int32()),
            _field("forward_alt", pa.int32()),
            _field("reverse_alt", pa.int32()),
            _field("f1r2_alt", pa.int32()),
            _field("f2r1_alt", pa.int32()),
            _field("median_base_quality", pa.float64()),
            _field("median_fragment_length", pa.float64()),
            _field("median_mapping_quality", pa.float64()),
            _field("median_read_position", pa.float64()),
            _field("tumor_lod", pa.float64()),
            _field("population_af_neg_log10", pa.float64()),
            _field("unique_alt_read_count", pa.int32()),
            _field("original_alignment_mismatch_reads", pa.int32()),
            _field("allele_filter_status", pa.list_(pa.string()), nullable=False),
            _field("short_tandem_repeat", pa.bool_(), nullable=False),
            _field("repeat_unit", pa.string()),
            _field("repeat_count_ref", pa.int32()),
            _field("repeat_count_alt", pa.int32()),
            _field("numt_or_paralogous_context", pa.bool_(), nullable=False),
            _field("circular_boundary_audited", pa.bool_(), nullable=False),
        ]
    )


def read_backed_evidence_schema() -> pa.Schema:
    return _schema(
        [
            _field("target_id", pa.string(), nullable=False),
            _field("target_manifest_sha256", pa.string(), nullable=False),
            _field("chrom", pa.string(), nullable=False),
            _field("pos", pa.int32(), nullable=False),
            _field("ref", pa.string(), nullable=False),
            _field("alt", pa.string(), nullable=False),
            _field("variant_type", pa.string(), nullable=False),
            _field("assessment_status", pa.string(), nullable=False),
            _field("eligible_fragments", pa.int32(), nullable=False),
            _field("ref_reads", pa.int32(), nullable=False),
            _field("alt_reads", pa.int32(), nullable=False),
            _field("other_reads", pa.int32(), nullable=False),
            _field("forward_ref_reads", pa.int32(), nullable=False),
            _field("reverse_ref_reads", pa.int32(), nullable=False),
            _field("forward_alt_reads", pa.int32(), nullable=False),
            _field("reverse_alt_reads", pa.int32(), nullable=False),
            _field("f1r2_ref_reads", pa.int32(), nullable=False),
            _field("f2r1_ref_reads", pa.int32(), nullable=False),
            _field("unclassified_orientation_ref_reads", pa.int32(), nullable=False),
            _field("f1r2_alt_reads", pa.int32(), nullable=False),
            _field("f2r1_alt_reads", pa.int32(), nullable=False),
            _field("unclassified_orientation_alt_reads", pa.int32(), nullable=False),
            _field("mean_ref_mapping_quality", pa.float64()),
            _field("mean_alt_mapping_quality", pa.float64()),
            _field("mean_ref_base_quality", pa.float64()),
            _field("mean_alt_base_quality", pa.float64()),
            _field("mean_ref_read_position", pa.float64()),
            _field("mean_alt_read_position", pa.float64()),
            _field("soft_clipped_ref_reads", pa.int32(), nullable=False),
            _field("soft_clipped_alt_reads", pa.int32(), nullable=False),
            _field("nearby_non_target_indel_reads", pa.int32(), nullable=False),
            _field("nearby_indel_window", pa.int32(), nullable=False),
            _field("excluded_low_mapping_quality_reads", pa.int32(), nullable=False),
            _field("excluded_low_base_quality_reads", pa.int32(), nullable=False),
            _field("excluded_duplicate_reads", pa.int32(), nullable=False),
            _field("excluded_nonprimary_reads", pa.int32(), nullable=False),
            _field("excluded_qc_fail_reads", pa.int32(), nullable=False),
            _field("collapsed_overlapping_mates", pa.int32(), nullable=False),
            _field("discordant_overlapping_mates", pa.int32(), nullable=False),
            _field("minimum_mapping_quality", pa.int32(), nullable=False),
            _field("minimum_base_quality", pa.int32(), nullable=False),
            _field("duplicate_reads_included", pa.bool_(), nullable=False),
            _field("overlapping_mate_policy", pa.string(), nullable=False),
            _field("local_complexity_status", pa.string(), nullable=False),
            _field("local_haplotype_status", pa.string(), nullable=False),
            _field("candidate_category", pa.string(), nullable=False),
        ]
    )


SCHEMAS = {
    "structural_variants": structural_variant_schema,
    "copy_number": copy_number_schema,
    "repeat_expansions": repeat_expansion_schema,
    "hla_typing": hla_typing_schema,
    "complex_pharmacogenes": complex_pharmacogene_schema,
    "mitochondrial_heteroplasmy": mitochondrial_heteroplasmy_schema,
    "read_backed_evidence": read_backed_evidence_schema,
}

SORT_FIELDS = {
    "structural_variants": ("chrom1", "pos1", "chrom2", "pos2"),
    "copy_number": ("chrom", "start_pos", "end_pos"),
    "repeat_expansions": ("chrom", "start_pos", "end_pos", "locus_id"),
    "hla_typing": ("locus",),
    "complex_pharmacogenes": ("gene",),
    "mitochondrial_heteroplasmy": ("chrom", "pos", "ref", "alt"),
    "read_backed_evidence": ("chrom", "pos", "ref", "alt", "target_id"),
}


def schema_for(domain: str) -> pa.Schema:
    try:
        return SCHEMAS[domain]()
    except KeyError as exc:
        raise PortableEvidenceError(f"unknown read-backed domain: {domain}") from exc


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, float) and (value != value or value in {float("inf"), float("-inf")}):
        raise PortableEvidenceError("non-finite numeric evidence is forbidden")
    return value


def make_evidence_id(domain: str, row: Mapping[str, Any]) -> str:
    if domain not in SCHEMAS:
        raise PortableEvidenceError(f"unknown read-backed domain: {domain}")
    payload = {key: value for key, value in row.items() if key != "evidence_id"}
    encoded = json.dumps(
        _jsonable(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode()
    return f"{domain}:{hashlib.sha256(encoded).hexdigest()}"


def _walk_strings(value: Any) -> Iterable[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _walk_strings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _walk_strings(item)


def _require_sha256(row: Mapping[str, Any], field: str) -> None:
    value = row.get(field)
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise PortableEvidenceError(f"{field} is not a lowercase SHA-256")


def _require_nonnegative(row: Mapping[str, Any], fields: Sequence[str]) -> None:
    for field in fields:
        value = row.get(field)
        if value is not None and (not isinstance(value, (int, float)) or value < 0):
            raise PortableEvidenceError(f"{field} must be non-negative when present")


def _require_interval(
    row: Mapping[str, Any], chrom_field: str, start_field: str, end_field: str
) -> None:
    chrom = row.get(chrom_field)
    start = row.get(start_field)
    end = row.get(end_field)
    if chrom not in CHROM_RANK:
        raise PortableEvidenceError(f"non-canonical chromosome in {chrom_field}")
    if not isinstance(start, int) or not isinstance(end, int) or start < 1 or end < start:
        raise PortableEvidenceError(f"invalid interval {chrom_field}:{start_field}-{end_field}")


def _domain_checks(domain: str, row: Mapping[str, Any]) -> None:
    if domain == "structural_variants":
        for chrom_field, pos_field in (("chrom1", "pos1"), ("chrom2", "pos2")):
            if row.get(chrom_field) not in CHROM_RANK or not isinstance(row.get(pos_field), int):
                raise PortableEvidenceError("invalid structural-variant breakpoint")
            if int(row[pos_field]) < 1:
                raise PortableEvidenceError("invalid structural-variant breakpoint")
        if row.get("sv_type") not in {"DEL", "DUP", "INV", "INS", "BND", "CNV"}:
            raise PortableEvidenceError("unsupported structural-variant class")
        if not row.get("caller_record_id"):
            raise PortableEvidenceError("structural-variant caller record ID is missing")
        if row.get("call_status") not in {"non_reference", "reference", "no_call"}:
            raise PortableEvidenceError("invalid structural-variant call status")
        genotype = row.get("genotype")
        if genotype is not None and (
            not isinstance(genotype, list)
            or not genotype
            or any(allele not in {None, 0, 1} for allele in genotype)
        ):
            raise PortableEvidenceError("invalid structural-variant genotype")
        if row.get("call_status") == "non_reference" and (
            not genotype or not any(allele == 1 for allele in genotype)
        ):
            raise PortableEvidenceError("non-reference SV status lacks an ALT genotype")
        if row.get("call_status") == "reference" and (
            not genotype or any(allele != 0 for allele in genotype)
        ):
            raise PortableEvidenceError("reference SV status lacks a 0/0 genotype")
        if row.get("call_status") == "no_call" and genotype is not None:
            raise PortableEvidenceError("no-call SV status contains a genotype")
        if len(row.get("cipos") or ()) != 2 or len(row.get("ciend") or ()) != 2:
            raise PortableEvidenceError("SV breakpoint confidence intervals require two bounds")
        if row.get("mate_link_status") not in {
            "not_applicable",
            "single_record_partner_coordinates",
            "paired_mate_id",
            "unresolved",
        }:
            raise PortableEvidenceError("invalid BND mate-link status")
        _require_nonnegative(
            row,
            (
                "quality",
                "genotype_quality",
                "copy_number",
                "paired_end_support",
                "caller_site_paired_end_support",
                "paired_end_reference_reads",
                "split_read_support",
                "caller_site_split_read_support",
                "split_read_reference_reads",
                "read_depth_support",
                "regional_depth",
                "mapping_quality",
                "split_read_mapping_quality",
                "insertion_length",
                "homology_length",
            ),
        )
        likelihoods = row.get("genotype_likelihoods")
        if not isinstance(likelihoods, list) or any(
            not isinstance(value, int) or value < 0 for value in likelihoods
        ):
            raise PortableEvidenceError("invalid structural-variant genotype likelihoods")
    elif domain == "copy_number":
        _require_interval(row, "chrom", "start_pos", "end_pos")
        if row.get("call_type") not in {"deletion", "duplication", "cnnloh"}:
            raise PortableEvidenceError("invalid copy-number call type")
        if not row.get("caller_model") or not row.get("copy_number_basis"):
            raise PortableEvidenceError("copy-number model and basis are required")
        _require_nonnegative(
            row,
            (
                "copy_number",
                "read_depth_copy_number",
                "major_copy_number",
                "minor_copy_number",
                "bin_size",
                "bins_observed",
                "heterozygous_sites",
                "homozygous_sites",
                "sex_ploidy",
                "normalized_read_depth",
                "baf_deviation_from_half",
                "combined_p_value",
                "caller_reported_rd_p_value",
                "baf_p_value",
                "likelihood_deletion",
                "likelihood_copy_neutral_loh",
                "likelihood_duplication",
                "e_value_1",
                "e_value_2",
                "e_value_3",
                "e_value_4",
                "q0_fraction",
                "reference_gap_fraction",
                "segment_gap_fraction",
                "strict_mask_fraction",
                "distance_to_large_gap",
            ),
        )
        for field in (
            "baf",
            "baf_deviation_from_half",
            "q0_fraction",
            "reference_gap_fraction",
            "segment_gap_fraction",
            "strict_mask_fraction",
        ):
            value = row.get(field)
            if value is not None and not 0 <= value <= 1:
                raise PortableEvidenceError(f"{field} must be within [0,1]")
        baf_source = row.get("source_baf_vcf_sha256")
        if baf_source is not None and (
            not isinstance(baf_source, str) or SHA256.fullmatch(baf_source) is None
        ):
            raise PortableEvidenceError("invalid copy-number BAF source digest")
        models = row.get("copy_number_models")
        if not isinstance(models, list):
            raise PortableEvidenceError("copy-number models must be an explicit list")
        for expected_rank, model in enumerate(models, start=1):
            if not isinstance(model, dict) or model.get("rank") != expected_rank:
                raise PortableEvidenceError("copy-number model ranks are not contiguous")
            values = (
                model.get("total_copy_number"),
                model.get("haplotype1_copy_number"),
                model.get("haplotype2_copy_number"),
            )
            if any(not isinstance(value, int) or value < 0 for value in values):
                raise PortableEvidenceError("invalid copy-number model genotype")
            if values[1] + values[2] != values[0]:
                raise PortableEvidenceError("copy-number model haplotypes do not sum to total")
            for field in ("likelihood", "cell_fraction"):
                value = model.get(field)
                if not isinstance(value, (int, float)) or value < 0:
                    raise PortableEvidenceError("invalid copy-number model score")
            if model["cell_fraction"] > 1:
                raise PortableEvidenceError("copy-number model cell fraction exceeds one")
            if not isinstance(model.get("reported_in_tsv"), bool) or not isinstance(
                model.get("padding"), bool
            ):
                raise PortableEvidenceError("copy-number model provenance is incomplete")
        if row.get("caller_model") == "combined_rd_baf":
            required = (
                "copy_number",
                "read_depth_copy_number",
                "major_copy_number",
                "minor_copy_number",
                "baf_deviation_from_half",
                "combined_p_value",
                "caller_reported_rd_p_value",
                "baf_p_value",
                "heterozygous_sites",
                "homozygous_sites",
                "source_baf_vcf_sha256",
            )
            if any(row.get(field) is None for field in required) or len(models) != 10:
                raise PortableEvidenceError("combined RD+BAF copy-number evidence is incomplete")
            if row.get("caller_reported_rd_p_value_duplicates_combined") is not True:
                raise PortableEvidenceError(
                    "combined RD+BAF caller p-value duplication is undeclared"
                )
            if row.get("baf_mask_applied") is not True:
                raise PortableEvidenceError("combined RD+BAF evidence lacks its strict SNP mask")
            if not math.isclose(float(row["copy_number"]), float(models[0]["total_copy_number"])):
                raise PortableEvidenceError(
                    "selected copy number differs from the top combined model"
                )
    elif domain == "repeat_expansions":
        _require_interval(row, "chrom", "start_pos", "end_pos")
        if not row.get("locus_id") or not row.get("variant_id"):
            raise PortableEvidenceError("repeat locus and variant identity are required")
        motif = row.get("motif")
        if not isinstance(motif, str) or IUPAC_DNA_MOTIF.fullmatch(motif.upper()) is None:
            raise PortableEvidenceError("repeat motif must contain only valid IUPAC DNA symbols")
        if row.get("ploidy") not in {1, 2}:
            raise PortableEvidenceError("repeat ploidy must be one or two")
        for prefix in ("allele1", "allele2"):
            count = row.get(f"{prefix}_repeat_count")
            lower = row.get(f"{prefix}_ci_lower")
            upper = row.get(f"{prefix}_ci_upper")
            lower_bound = row.get(f"{prefix}_count_is_lower_bound")
            if not isinstance(lower_bound, bool):
                raise PortableEvidenceError("repeat lower-bound status must be explicit")
            if count is None:
                if lower is not None or upper is not None or lower_bound:
                    raise PortableEvidenceError("repeat count and confidence interval are atomic")
            else:
                if not isinstance(count, int) or count < 0:
                    raise PortableEvidenceError("repeat count must be non-negative")
                if not isinstance(lower, int) or lower < 0:
                    raise PortableEvidenceError("repeat confidence lower bound is required")
                if upper is not None and (not isinstance(upper, int) or upper < lower):
                    raise PortableEvidenceError("repeat confidence upper bound is invalid")
                if upper is None and not lower_bound:
                    raise PortableEvidenceError("finite repeat call lacks an upper bound")
                if count < lower or (upper is not None and count > upper):
                    raise PortableEvidenceError("repeat count lies outside its confidence interval")
        if row.get("ploidy") == 1 and row.get("allele2_repeat_count") is not None:
            raise PortableEvidenceError("haploid repeat call contains a second allele")
        _require_nonnegative(
            row,
            (
                "reference_repeat_count",
                "locus_coverage",
                "spanning_reads",
                "flanking_reads",
                "in_repeat_reads",
                "allele1_spanning_reads",
                "allele1_flanking_reads",
                "allele1_in_repeat_reads",
                "allele2_spanning_reads",
                "allele2_flanking_reads",
                "allele2_in_repeat_reads",
            ),
        )
    elif domain == "hla_typing":
        locus = row.get("locus")
        if not isinstance(locus, str) or re.fullmatch(r"HLA-[A-Z0-9]+", locus) is None:
            raise PortableEvidenceError("HLA locus must use the HLA- prefix")
        if row.get("call_status") not in {"diploid_call", "partial_call", "no_call"}:
            raise PortableEvidenceError("invalid HLA call status")
        if row.get("allele1") is None and row.get("allele2") is not None:
            raise PortableEvidenceError("HLA allele2 cannot exist without allele1")
        allele_count = sum(row.get(field) is not None for field in ("allele1", "allele2"))
        expected_status = {0: "no_call", 1: "partial_call", 2: "diploid_call"}[allele_count]
        if row.get("call_status") != expected_status:
            raise PortableEvidenceError("HLA call status disagrees with allele cardinality")
        if expected_status != "diploid_call" and row.get("confidence_tier") == "SUPPORTED":
            raise PortableEvidenceError("incomplete HLA call cannot be supported")
        gene = locus.removeprefix("HLA-")
        for slot in ("allele1", "allele2"):
            allele = row.get(slot)
            if allele is not None and not str(allele).startswith(f"{gene}*"):
                raise PortableEvidenceError("HLA allele and locus disagree")
            two_field = row.get(f"{slot}_two_field")
            if allele is None and two_field is not None:
                raise PortableEvidenceError("missing HLA allele has a two-field projection")
            if two_field is not None and not str(two_field).startswith(f"{gene}*"):
                raise PortableEvidenceError("HLA two-field allele and locus disagree")
            g_group = row.get(f"{slot}_g_group")
            if g_group is not None and (
                not str(g_group).startswith(f"{gene}*") or not str(g_group).endswith("G")
            ):
                raise PortableEvidenceError("invalid HLA G-group allele")
            candidates = row.get(f"{slot}_candidates")
            scored = row.get(f"{slot}_scored_candidates")
            if not isinstance(candidates, list) or not isinstance(scored, list):
                raise PortableEvidenceError("HLA candidate lists must be explicit")
            if allele is not None and allele not in candidates:
                raise PortableEvidenceError("selected HLA allele is absent from its candidates")
            for candidate in candidates:
                if not isinstance(candidate, str) or not candidate.startswith(f"{gene}*"):
                    raise PortableEvidenceError("invalid HLA candidate allele")
            for metric in scored:
                if not isinstance(metric, dict) or not metric.get("allele"):
                    raise PortableEvidenceError("invalid scored HLA candidate")
                score = metric.get("score")
                if not isinstance(score, (int, float)) or not 0 <= score <= 100:
                    raise PortableEvidenceError("HLA candidate score must be within [0,100]")
                for frequency_field in (
                    "caucasian_frequency",
                    "black_frequency",
                    "asian_frequency",
                ):
                    frequency = metric.get(frequency_field)
                    if frequency is not None and not 0 <= frequency <= 1:
                        raise PortableEvidenceError("HLA population frequency must be within [0,1]")
            quality = row.get(f"{slot}_quality_score")
            if quality is not None and not 0 <= quality <= 100:
                raise PortableEvidenceError("HLA quality score must be within [0,100]")
            metrics = row.get(f"{slot}_hlala_metrics")
            if metrics is not None:
                if not isinstance(metrics, dict):
                    raise PortableEvidenceError("HLA*LA metrics must be a typed object")
                for coverage_field in (
                    "average_coverage",
                    "coverage_first_decile",
                    "minimum_coverage",
                    "locus_avg_column_error",
                    "unaccounted_columns_fraction_gt_0_2",
                ):
                    value = metrics.get(coverage_field)
                    if not isinstance(value, (int, float)) or value < 0:
                        raise PortableEvidenceError("HLA*LA non-negative metric is invalid")
                proportion = metrics.get("proportion_kmers_covered")
                if not isinstance(proportion, (int, float)) or not (
                    proportion == -1 or 0 <= proportion <= 1
                ):
                    raise PortableEvidenceError("HLA*LA k-mer coverage proportion is invalid")
                for required in ("q1", "q2"):
                    if not isinstance(metrics.get(required), (int, float)):
                        raise PortableEvidenceError("HLA*LA quality metric is invalid")
                if not isinstance(metrics.get("perfect_g"), bool):
                    raise PortableEvidenceError("HLA*LA perfectG status is invalid")
        for field in (
            "caller_database_name",
            "caller_database_version",
            "nomenclature_database_name",
            "nomenclature_version",
            "allele_resolution",
        ):
            if not isinstance(row.get(field), str) or not row[field].strip():
                raise PortableEvidenceError(f"HLA {field} is missing")
        _require_sha256(row, "caller_database_identity_sha256")
        _require_sha256(row, "nomenclature_identity_sha256")
        _require_nonnegative(row, ("locus_depth",))
    elif domain == "complex_pharmacogenes":
        if not row.get("gene") or not row.get("definition_version"):
            raise PortableEvidenceError("complex pharmacogene identity is incomplete")
        if row.get("call_status") not in {"called", "ambiguous", "no_call"}:
            raise PortableEvidenceError("invalid complex pharmacogene call status")
        _require_sha256(row, "definition_identity_sha256")
        _require_nonnegative(
            row,
            (
                "gene_copy_number",
                "locus_total_copy_number",
                "paralog_marker_copy_number",
                "raw_normalized_locus_depth",
                "raw_normalized_paralog_depth",
                "activity_score",
                "coverage_mad",
                "median_depth",
                "exon_specific_copy_number",
            ),
        )
        if not isinstance(row.get("candidate_diplotypes"), list):
            raise PortableEvidenceError("complex pharmacogene candidates must be explicit")
        if row.get("diplotype") is not None and row["diplotype"] not in row["candidate_diplotypes"]:
            raise PortableEvidenceError("selected diplotype is absent from candidate diplotypes")
        if row.get("call_status") == "no_call" and row.get("diplotype") is not None:
            raise PortableEvidenceError("complex pharmacogene no-call contains a diplotype")
        if row.get("call_status") == "called" and row.get("diplotype") is None:
            raise PortableEvidenceError("called complex pharmacogene lacks a diplotype")
        if row.get("call_status") == "ambiguous" and not row.get("ambiguous"):
            raise PortableEvidenceError("ambiguous complex call lacks its ambiguity flag")
        candidate_details = row.get("candidate_solution_details")
        if not isinstance(candidate_details, list):
            raise PortableEvidenceError("complex pharmacogene solution details must be explicit")
        solution_ids: set[int] = set()
        for candidate in candidate_details:
            if not isinstance(candidate, dict):
                raise PortableEvidenceError("invalid complex pharmacogene candidate solution")
            solution_id = candidate.get("solution_id")
            if not isinstance(solution_id, int) or solution_id <= 0 or solution_id in solution_ids:
                raise PortableEvidenceError("invalid complex pharmacogene solution identifier")
            solution_ids.add(solution_id)
            for field in ("solution_expression", "major_diplotype"):
                if not isinstance(candidate.get(field), str) or not candidate[field]:
                    raise PortableEvidenceError("incomplete complex pharmacogene candidate")
            minor_alleles = candidate.get("minor_alleles")
            copy_count = candidate.get("copy_count")
            if (
                not isinstance(minor_alleles, list)
                or not all(isinstance(allele, str) and allele for allele in minor_alleles)
                or not isinstance(copy_count, int)
                or copy_count < 0
                or len(minor_alleles) != copy_count
            ):
                raise PortableEvidenceError("invalid complex pharmacogene candidate copies")
            score = candidate.get("activity_score")
            if score is not None and (not isinstance(score, (int, float)) or score < 0):
                raise PortableEvidenceError("invalid complex pharmacogene activity score")
        decomposition = row.get("allele_decomposition")
        if not isinstance(decomposition, list):
            raise PortableEvidenceError("complex pharmacogene decomposition must be explicit")
        decomposition_keys: set[tuple[Any, ...]] = set()
        for item in decomposition:
            if not isinstance(item, dict):
                raise PortableEvidenceError("invalid complex pharmacogene decomposition row")
            solution_id = item.get("solution_id")
            copy_index = item.get("copy_index")
            if (
                not isinstance(solution_id, int)
                or solution_id not in solution_ids
                or not isinstance(copy_index, int)
                or copy_index < 0
                or not isinstance(item.get("allele"), str)
                or not item["allele"]
            ):
                raise PortableEvidenceError("invalid complex pharmacogene decomposition identity")
            position = item.get("position_zero_based")
            coverage = item.get("coverage")
            if position is not None and (not isinstance(position, int) or position < 0):
                raise PortableEvidenceError("invalid complex pharmacogene decomposition position")
            if coverage is not None and (not isinstance(coverage, (int, float)) or coverage < 0):
                raise PortableEvidenceError("invalid complex pharmacogene decomposition coverage")
            key = (
                solution_id,
                copy_index,
                item["allele"],
                position,
                item.get("nucleotide_change"),
            )
            if key in decomposition_keys:
                raise PortableEvidenceError("duplicate complex pharmacogene decomposition row")
            decomposition_keys.add(key)
        for count in row.get("raw_variant_counts") or ():
            if not isinstance(count, dict) or not count.get("variant_id"):
                raise PortableEvidenceError("invalid complex pharmacogene raw variant count")
            for field in ("alt_reads", "ref_reads", "forward_alt_reads", "reverse_alt_reads"):
                value = count.get(field)
                if value is not None and (not isinstance(value, int) or value < 0):
                    raise PortableEvidenceError("negative complex pharmacogene read count")
            forward = count.get("forward_alt_reads")
            reverse = count.get("reverse_alt_reads")
            if (forward is None) != (reverse is None):
                raise PortableEvidenceError("partial complex pharmacogene strand accounting")
            if forward is not None and forward + reverse != count["alt_reads"]:
                raise PortableEvidenceError("complex pharmacogene ALT strand accounting differs")
        if row.get("ambiguous") and row.get("outside_call_eligible"):
            raise PortableEvidenceError("ambiguous complex call cannot become an outside call")
        if row.get("call_status") != "called" and row.get("outside_call_eligible"):
            raise PortableEvidenceError("non-called complex evidence cannot become an outside call")
    elif domain == "mitochondrial_heteroplasmy":
        if row.get("chrom") != "chrM" or not isinstance(row.get("pos"), int):
            raise PortableEvidenceError("mitochondrial evidence must use chrM coordinates")
        if not 1 <= int(row["pos"]) <= 16_569:
            raise PortableEvidenceError("mitochondrial coordinate is out of bounds")
        for field in ("ref", "alt"):
            allele = row.get(field)
            if not isinstance(allele, str) or CONCRETE_ALLELE.fullmatch(allele) is None:
                raise PortableEvidenceError("mitochondrial alleles must be concrete uppercase DNA")
        fraction = row.get("allele_fraction")
        if not isinstance(fraction, (int, float)) or not 0 <= fraction <= 1:
            raise PortableEvidenceError("mitochondrial allele fraction must be within [0,1]")
        _require_nonnegative(
            row,
            (
                "depth",
                "ref_depth",
                "alt_depth",
                "forward_alt",
                "reverse_alt",
                "f1r2_alt",
                "f2r1_alt",
                "median_base_quality",
                "median_fragment_length",
                "median_mapping_quality",
                "median_read_position",
                "population_af_neg_log10",
                "unique_alt_read_count",
                "original_alignment_mismatch_reads",
                "repeat_count_ref",
                "repeat_count_alt",
            ),
        )
        if int(row["depth"]) == 0:
            raise PortableEvidenceError("mitochondrial evidence has zero depth")
        if (
            row.get("ref_depth") is not None
            and row.get("alt_depth") is not None
            and int(row["ref_depth"]) + int(row["alt_depth"]) > int(row["depth"])
        ):
            raise PortableEvidenceError("mitochondrial allele depths exceed total depth")
        if not isinstance(row.get("allele_filter_status"), list) or not row["allele_filter_status"]:
            raise PortableEvidenceError("mitochondrial allele filter status is missing")
        if row.get("short_tandem_repeat") and not row.get("repeat_unit"):
            raise PortableEvidenceError("mitochondrial STR evidence lacks its repeat unit")
    elif domain == "read_backed_evidence":
        if row.get("chrom") not in CHROM_RANK or not isinstance(row.get("pos"), int):
            raise PortableEvidenceError("invalid local read-back coordinate")
        if int(row["pos"]) < 1:
            raise PortableEvidenceError("invalid local read-back coordinate")
        for field in ("ref", "alt"):
            allele = row.get(field)
            if not isinstance(allele, str) or CONCRETE_ALLELE.fullmatch(allele) is None:
                raise PortableEvidenceError("local read-back alleles must be concrete DNA")
        if row.get("ref") == row.get("alt"):
            raise PortableEvidenceError("local read-back REF and ALT are identical")
        if row.get("variant_type") not in {
            "SNV",
            "MNV",
            "insertion",
            "deletion",
            "complex",
        }:
            raise PortableEvidenceError("invalid local read-back variant type")
        if row.get("assessment_status") not in {
            "assessed",
            "no_eligible_reads",
            "unsupported_allele_representation",
        }:
            raise PortableEvidenceError("invalid local read-back assessment status")
        _require_nonnegative(
            row,
            (
                "eligible_fragments",
                "ref_reads",
                "alt_reads",
                "other_reads",
                "forward_ref_reads",
                "reverse_ref_reads",
                "forward_alt_reads",
                "reverse_alt_reads",
                "f1r2_ref_reads",
                "f2r1_ref_reads",
                "unclassified_orientation_ref_reads",
                "f1r2_alt_reads",
                "f2r1_alt_reads",
                "unclassified_orientation_alt_reads",
                "mean_ref_mapping_quality",
                "mean_alt_mapping_quality",
                "mean_ref_base_quality",
                "mean_alt_base_quality",
                "mean_ref_read_position",
                "mean_alt_read_position",
                "soft_clipped_ref_reads",
                "soft_clipped_alt_reads",
                "nearby_non_target_indel_reads",
                "nearby_indel_window",
                "excluded_low_mapping_quality_reads",
                "excluded_low_base_quality_reads",
                "excluded_duplicate_reads",
                "excluded_nonprimary_reads",
                "excluded_qc_fail_reads",
                "collapsed_overlapping_mates",
                "discordant_overlapping_mates",
                "minimum_mapping_quality",
                "minimum_base_quality",
            ),
        )
        _require_sha256(row, "target_manifest_sha256")
        if row["ref_reads"] + row["alt_reads"] + row["other_reads"] != row["eligible_fragments"]:
            raise PortableEvidenceError("local read-back fragment accounting differs")
        if row["forward_ref_reads"] + row["reverse_ref_reads"] != row["ref_reads"]:
            raise PortableEvidenceError("local read-back REF strand accounting differs")
        if row["forward_alt_reads"] + row["reverse_alt_reads"] != row["alt_reads"]:
            raise PortableEvidenceError("local read-back ALT strand accounting differs")
        if (
            row["f1r2_ref_reads"]
            + row["f2r1_ref_reads"]
            + row["unclassified_orientation_ref_reads"]
            != row["ref_reads"]
        ):
            raise PortableEvidenceError("local read-back REF orientation accounting differs")
        if (
            row["f1r2_alt_reads"]
            + row["f2r1_alt_reads"]
            + row["unclassified_orientation_alt_reads"]
            != row["alt_reads"]
        ):
            raise PortableEvidenceError("local read-back ALT orientation accounting differs")
        if row["assessment_status"] == "no_eligible_reads" and row["eligible_fragments"] != 0:
            raise PortableEvidenceError("no-read assessment contains eligible fragments")
        if row["assessment_status"] == "assessed" and row["eligible_fragments"] == 0:
            raise PortableEvidenceError("assessed local target has no eligible fragments")
        if row["overlapping_mate_policy"] not in {"collapse_by_query_name"}:
            raise PortableEvidenceError("invalid overlapping-mate policy")
        if row["duplicate_reads_included"] and row["excluded_duplicate_reads"] != 0:
            raise PortableEvidenceError("duplicate inclusion policy contradicts exclusion count")
        if not row.get("candidate_category") or not row.get("local_haplotype_status"):
            raise PortableEvidenceError("local read-back target annotations are incomplete")


def _sort_value(field: str, value: Any) -> Any:
    if field.startswith("chrom") and value in CHROM_RANK:
        return CHROM_RANK[value]
    return value


def row_sort_key(domain: str, row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        *(_sort_value(field, row.get(field)) for field in SORT_FIELDS[domain]),
        row["event_id"],
        row["technical_source_role"],
        row["caller_name"],
        row["evidence_id"],
    )


def validate_rows(
    domain: str,
    rows: Sequence[Mapping[str, Any]],
    *,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: Iterable[str] = (),
) -> list[dict[str, Any]]:
    schema = schema_for(domain)
    if not rows:
        raise PortableEvidenceError("empty optional read-backed extensions are forbidden")
    expected_fields = set(schema.names)
    forbidden = {value for value in forbidden_bindings if value}
    forbidden.update(hashlib.sha256(value.encode()).hexdigest() for value in tuple(forbidden))
    validated: list[dict[str, Any]] = []
    seen: set[str] = set()
    for ordinal, source in enumerate(rows, start=1):
        row = dict(source)
        if set(row) != expected_fields:
            missing = sorted(expected_fields - set(row))
            extra = sorted(set(row) - expected_fields)
            raise PortableEvidenceError(
                f"row {ordinal} physical fields differ; missing={missing}, extra={extra}"
            )
        for field in schema:
            if not field.nullable and row[field.name] is None:
                raise PortableEvidenceError(f"row {ordinal} required field is null: {field.name}")
        for field in (
            "event_id",
            "caller_name",
            "caller_version",
            "raw_artifact_id",
        ):
            if not isinstance(row[field], str) or not row[field].strip():
                raise PortableEvidenceError(f"row {ordinal} has an empty {field}")
        if row.get("assembly") != "GRCh38":
            raise PortableEvidenceError("read-backed evidence is not explicitly GRCh38")
        for field in (
            "reference_sha256",
            "caller_run_identity_sha256",
            "resource_manifest_sha256",
            "source_alignment_sha256",
            "raw_artifact_sha256",
            "source_record_digest",
        ):
            _require_sha256(row, field)
        if expected_reference_sha256 and row["reference_sha256"] != expected_reference_sha256:
            raise PortableEvidenceError("read-backed row uses the wrong reference digest")
        if row.get("technical_source_role") not in SOURCE_ROLES:
            raise PortableEvidenceError("invalid technical source role")
        if row.get("confidence_tier") not in CONFIDENCE_TIERS:
            raise PortableEvidenceError("invalid read-backed confidence tier")
        if row.get("concordance_status") not in CONCORDANCE_STATES:
            raise PortableEvidenceError("invalid concordance state")
        if row["concordance_status"] == "not_assessed" and row["concordance_group_id"] is not None:
            raise PortableEvidenceError("unassessed evidence cannot claim a concordance group")
        if row["concordance_status"] != "not_assessed" and not row["concordance_group_id"]:
            raise PortableEvidenceError("assessed evidence lacks a concordance group")
        if not isinstance(row.get("filter_reasons"), list):
            raise PortableEvidenceError("filter reasons must be an explicit list")
        if not isinstance(row.get("limitations"), list) or not row["limitations"]:
            raise PortableEvidenceError("at least one explicit limitation is required")
        if row["confidence_tier"] != "SUPPORTED" and row.get("confirmation_required") is not True:
            raise PortableEvidenceError("review/unresolved evidence requires confirmation")
        artifact_id = row.get("raw_artifact_id")
        if (
            not isinstance(artifact_id, str)
            or not artifact_id
            or "/" in artifact_id
            or "\\" in artifact_id
        ):
            raise PortableEvidenceError("raw artifact ID must be a path-free logical identifier")
        for text in _walk_strings(row):
            if FORBIDDEN_PORTABLE_TEXT.search(text):
                raise PortableEvidenceError("portable evidence contains a private path or secret")
            if any(contains_forbidden_binding(text, binding) for binding in forbidden):
                raise PortableEvidenceError(
                    "portable evidence contains a forbidden private binding"
                )
        _domain_checks(domain, row)
        expected_evidence_id = make_evidence_id(domain, row)
        if row.get("evidence_id") != expected_evidence_id:
            raise PortableEvidenceError("evidence ID is not the canonical row digest")
        if expected_evidence_id in seen:
            raise PortableEvidenceError("duplicate read-backed evidence row")
        seen.add(expected_evidence_id)
        validated.append(row)
    return sorted(validated, key=lambda row: row_sort_key(domain, row))


def canonical_row_multiset(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    digests = [
        hashlib.sha256(
            json.dumps(
                _jsonable(dict(row)), sort_keys=True, separators=(",", ":"), ensure_ascii=True
            ).encode()
        ).digest()
        for row in rows
    ]
    xor = bytearray(32)
    total = 0
    modulus = 1 << 256
    for digest in digests:
        total = (total + int.from_bytes(digest, "big")) % modulus
        for index, value in enumerate(digest):
            xor[index] ^= value
    return {
        "count": len(digests),
        "sha256_xor": bytes(xor).hex(),
        "sha256_sum_mod_2_256": total.to_bytes(32, "big").hex(),
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_extension(
    domain: str,
    rows: Sequence[Mapping[str, Any]],
    destination: str | Path,
    *,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: Iterable[str] = (),
) -> dict[str, Any]:
    destination = Path(destination)
    if destination.exists():
        raise PortableEvidenceError("refusing to replace an existing extension")
    destination.parent.mkdir(parents=True, exist_ok=True)
    validated = validate_rows(
        domain,
        rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    table = pa.Table.from_pylist(validated, schema=schema_for(domain))
    descriptor, name = tempfile.mkstemp(
        prefix=f".{destination.name}.partial.", suffix=".parquet", dir=destination.parent
    )
    os.close(descriptor)
    temporary = Path(name)
    try:
        pq.write_table(
            table,
            temporary,
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
            write_statistics=True,
            row_group_size=65_536,
            data_page_version="2.0",
        )
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    observed = pq.read_table(destination)
    if observed.schema != schema_for(domain):
        raise PortableEvidenceError("written extension physical schema changed")
    observed_rows = observed.to_pylist()
    if observed_rows != validated:
        raise PortableEvidenceError("written extension row content/order changed")
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "PASS",
        "domain": domain,
        "rows": len(validated),
        "parquet_sha256": sha256_file(destination),
        "canonical_row_multiset": canonical_row_multiset(validated),
        "private_paths_persisted": False,
        "subject_scope": "bundle_subject_only",
    }
