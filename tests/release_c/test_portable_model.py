from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.portable_model import (
    DOMAINS,
    PortableEvidenceError,
    contains_forbidden_binding,
    make_evidence_id,
    schema_for,
    validate_rows,
    write_extension,
)

SHA = "a" * 64
REFERENCE = "1" * 64


def common(event_id: str) -> dict[str, object]:
    return {
        "evidence_id": "pending",
        "event_id": event_id,
        "assembly": "GRCh38",
        "reference_sha256": REFERENCE,
        "technical_source_role": "primary",
        "caller_name": "synthetic-caller",
        "caller_version": "1.0.0",
        "caller_run_identity_sha256": "2" * 64,
        "resource_manifest_sha256": "3" * 64,
        "source_alignment_sha256": "4" * 64,
        "raw_artifact_id": "caller-output-primary",
        "raw_artifact_sha256": "5" * 64,
        "source_record_digest": "6" * 64,
        "filter_reasons": [],
        "confidence_tier": "SUPPORTED",
        "limitations": ["synthetic short-read evidence"],
        "confirmation_required": True,
        "concordance_group_id": None,
        "concordance_status": "not_assessed",
    }


def finalize(domain: str, specific: dict[str, object]) -> dict[str, object]:
    row = {**specific, **common(str(specific.get("event_id", "synthetic-event")))}
    row["event_id"] = str(specific.get("event_id", "synthetic-event"))
    row["evidence_id"] = make_evidence_id(domain, row)
    return row


def rows() -> dict[str, dict[str, object]]:
    return {
        "structural_variants": finalize(
            "structural_variants",
            {
                "event_id": "sv:chr1:100-200:DEL",
                "caller_record_id": "synthetic-sv-1",
                "chrom1": "chr1",
                "pos1": 100,
                "chrom2": "chr1",
                "pos2": 200,
                "sv_type": "DEL",
                "orientation": "+-",
                "length": -100,
                "genotype": [0, 1],
                "call_status": "non_reference",
                "genotype_quality": 60.0,
                "genotype_likelihoods": [60, 0, 100],
                "copy_number": 1.0,
                "quality": 60.0,
                "imprecise": False,
                "cipos": [0, 0],
                "ciend": [0, 0],
                "paired_end_support": 12,
                "caller_site_paired_end_support": 12,
                "paired_end_reference_reads": 20,
                "split_read_support": 8,
                "caller_site_split_read_support": 8,
                "split_read_reference_reads": 15,
                "read_depth_support": 20,
                "regional_depth": 30,
                "mapping_quality": 60.0,
                "split_read_mapping_quality": 55.0,
                "insertion_length": 0,
                "homology_length": 2,
                "bnd_mate_id": None,
                "mate_link_status": "not_applicable",
            },
        ),
        "copy_number": finalize(
            "copy_number",
            {
                "event_id": "cnv:chr2:100-1000",
                "call_type": "duplication",
                "caller_model": "synthetic-depth-baf",
                "chrom": "chr2",
                "start_pos": 100,
                "end_pos": 1000,
                "copy_number": 3.0,
                "copy_number_basis": "caller_model",
                "read_depth_copy_number": 3.0,
                "major_copy_number": 2.0,
                "minor_copy_number": 1.0,
                "normalized_read_depth": 1.5,
                "log2_ratio": 0.58,
                "baf": 0.33,
                "baf_deviation_from_half": None,
                "combined_p_value": None,
                "caller_reported_rd_p_value": None,
                "baf_p_value": None,
                "caller_reported_rd_p_value_duplicates_combined": None,
                "likelihood_deletion": None,
                "likelihood_copy_neutral_loh": None,
                "likelihood_duplication": None,
                "e_value_1": 1e-8,
                "e_value_2": 1e-7,
                "e_value_3": 1e-6,
                "e_value_4": 1e-5,
                "q0_fraction": 0.1,
                "reference_gap_fraction": 0.0,
                "segment_gap_fraction": 0.0,
                "strict_mask_fraction": 0.9,
                "distance_to_large_gap": 100000,
                "rd_mask_applied": True,
                "baf_mask_applied": True,
                "bin_size": 100,
                "bins_observed": 9,
                "heterozygous_sites": 4,
                "homozygous_sites": None,
                "copy_number_models": [],
                "source_baf_vcf_sha256": None,
                "sex_ploidy": 2,
            },
        ),
        "repeat_expansions": finalize(
            "repeat_expansions",
            {
                "event_id": "repeat:chr3:100-120:CAG",
                "locus_id": "SYNTHETIC_REPEAT",
                "variant_id": "SYNTHETIC_REPEAT_VARIANT",
                "variant_subtype": "Repeat",
                "chrom": "chr3",
                "start_pos": 100,
                "end_pos": 120,
                "motif": "CAG",
                "ploidy": 2,
                "reference_repeat_count": 7,
                "allele1_repeat_count": 10,
                "allele1_ci_lower": 9,
                "allele1_ci_upper": 11,
                "allele1_count_is_lower_bound": False,
                "allele1_support_type": "SPANNING",
                "allele1_spanning_reads": 7,
                "allele1_flanking_reads": 5,
                "allele1_in_repeat_reads": 1,
                "allele2_repeat_count": 20,
                "allele2_ci_lower": 18,
                "allele2_ci_upper": 22,
                "allele2_count_is_lower_bound": False,
                "allele2_support_type": "SPANNING",
                "allele2_spanning_reads": 5,
                "allele2_flanking_reads": 3,
                "allele2_in_repeat_reads": 1,
                "locus_coverage": 30.0,
                "spanning_reads": 12,
                "flanking_reads": 8,
                "in_repeat_reads": 2,
                "off_target_evidence": False,
                "short_read_resolvable": True,
            },
        ),
        "hla_typing": finalize(
            "hla_typing",
            {
                "event_id": "hla:HLA-A",
                "locus": "HLA-A",
                "call_status": "diploid_call",
                "allele1": "A*01:01",
                "allele2": "A*02:01",
                "allele1_g_group": "A*01:01:01G",
                "allele2_g_group": "A*02:01:01G",
                "allele1_two_field": "A*01:01",
                "allele2_two_field": "A*02:01",
                "allele_resolution": "full-and-g-group",
                "phase_resolved": False,
                "ambiguous": False,
                "primary_locus": True,
                "allele1_candidates": ["A*01:01"],
                "allele2_candidates": ["A*02:01"],
                "allele1_scored_candidates": [],
                "allele2_scored_candidates": [],
                "allele1_quality_score": 100.0,
                "allele2_quality_score": 100.0,
                "allele1_hlala_metrics": None,
                "allele2_hlala_metrics": None,
                "locus_depth": 28.0,
                "caller_database_name": "synthetic HLA mapping graph",
                "caller_database_version": "synthetic-1",
                "caller_database_identity_sha256": "7" * 64,
                "nomenclature_database_name": "synthetic HLA nomenclature",
                "nomenclature_version": "synthetic-1",
                "nomenclature_identity_sha256": "8" * 64,
                "database_age_limitation": False,
            },
        ),
        "complex_pharmacogenes": finalize(
            "complex_pharmacogenes",
            {
                "event_id": "pgx:SYNTH1",
                "gene": "SYNTH1",
                "call_status": "called",
                "diplotype": "*1/*2",
                "candidate_diplotypes": ["*1/*2"],
                "star_alleles": ["*1", "*2"],
                "candidate_solution_details": [
                    {
                        "solution_id": 1,
                        "solution_expression": "synthetic solution",
                        "major_diplotype": "*1/*2",
                        "minor_alleles": ["1.001", "2.001"],
                        "copy_count": 2,
                        "activity_score": 1.5,
                        "phenotype": "synthetic",
                    }
                ],
                "allele_decomposition": [
                    {
                        "solution_id": 1,
                        "copy_index": 0,
                        "allele": "1.001",
                        "position_zero_based": 99,
                        "nucleotide_change": "A>G",
                        "coverage": 30.0,
                        "effect": "synthetic",
                        "dbsnp": "synthetic-variant",
                        "code": None,
                        "status": None,
                    }
                ],
                "gene_copy_number": 2.0,
                "locus_total_copy_number": 2.0,
                "paralog_marker_copy_number": 0.0,
                "raw_normalized_locus_depth": 2.0,
                "raw_normalized_paralog_depth": 0.0,
                "structural_configurations": [],
                "activity_score": 1.5,
                "phenotype": "synthetic",
                "coverage_mad": 0.1,
                "median_depth": 30.0,
                "variants_called": ["synthetic-variant"],
                "raw_star_allele": "*1_*2",
                "caller_filter": "PASS",
                "call_info": "synthetic",
                "exon_specific_copy_number": 2.0,
                "copy_number_consensus": [2.0],
                "differentiating_site_copy_numbers": [2.0],
                "differentiating_site_raw_depths": [2.0],
                "raw_variant_counts": [
                    {
                        "variant_id": "synthetic-variant",
                        "alt_reads": 10,
                        "ref_reads": 20,
                        "forward_alt_reads": 5,
                        "reverse_alt_reads": 5,
                    }
                ],
                "definition_version": "synthetic-1",
                "definition_identity_sha256": "8" * 64,
                "licence_gated_caller": False,
                "ambiguous": False,
                "outside_call_eligible": True,
            },
        ),
        "mitochondrial_heteroplasmy": finalize(
            "mitochondrial_heteroplasmy",
            {
                "event_id": "mt:chrM:100:A:G",
                "chrom": "chrM",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "allele_fraction": 0.2,
                "depth": 100,
                "ref_depth": 80,
                "alt_depth": 20,
                "forward_alt": 10,
                "reverse_alt": 10,
                "f1r2_alt": 10,
                "f2r1_alt": 10,
                "median_base_quality": 35.0,
                "median_fragment_length": 350.0,
                "median_mapping_quality": 60.0,
                "median_read_position": 75.0,
                "tumor_lod": 50.0,
                "population_af_neg_log10": 2.4,
                "unique_alt_read_count": 20,
                "original_alignment_mismatch_reads": 0,
                "allele_filter_status": ["PASS"],
                "short_tandem_repeat": False,
                "repeat_unit": None,
                "repeat_count_ref": None,
                "repeat_count_alt": None,
                "numt_or_paralogous_context": False,
                "circular_boundary_audited": True,
            },
        ),
        "read_backed_evidence": finalize(
            "read_backed_evidence",
            {
                "event_id": "readback:chr4:100:A:T",
                "target_id": "priority-1",
                "target_manifest_sha256": "9" * 64,
                "chrom": "chr4",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "variant_type": "SNV",
                "assessment_status": "assessed",
                "eligible_fragments": 30,
                "ref_reads": 20,
                "alt_reads": 10,
                "other_reads": 0,
                "forward_ref_reads": 10,
                "reverse_ref_reads": 10,
                "forward_alt_reads": 5,
                "reverse_alt_reads": 5,
                "f1r2_ref_reads": 10,
                "f2r1_ref_reads": 10,
                "unclassified_orientation_ref_reads": 0,
                "f1r2_alt_reads": 5,
                "f2r1_alt_reads": 5,
                "unclassified_orientation_alt_reads": 0,
                "mean_ref_mapping_quality": 60.0,
                "mean_alt_mapping_quality": 60.0,
                "mean_ref_base_quality": 35.0,
                "mean_alt_base_quality": 35.0,
                "mean_ref_read_position": 75.0,
                "mean_alt_read_position": 75.0,
                "soft_clipped_ref_reads": 0,
                "soft_clipped_alt_reads": 0,
                "nearby_non_target_indel_reads": 0,
                "nearby_indel_window": 10,
                "excluded_low_mapping_quality_reads": 0,
                "excluded_low_base_quality_reads": 0,
                "excluded_duplicate_reads": 0,
                "excluded_nonprimary_reads": 0,
                "excluded_qc_fail_reads": 0,
                "collapsed_overlapping_mates": 0,
                "discordant_overlapping_mates": 0,
                "minimum_mapping_quality": 20,
                "minimum_base_quality": 20,
                "duplicate_reads_included": False,
                "overlapping_mate_policy": "collapse_by_query_name",
                "local_complexity_status": "not_assessed",
                "local_haplotype_status": "supported",
                "candidate_category": "MOSAIC_CANDIDATE",
            },
        ),
    }


def test_all_domains_have_exact_typed_schema_and_validate() -> None:
    examples = rows()
    assert set(examples) == set(DOMAINS)
    for domain, row in examples.items():
        assert set(row) == set(schema_for(domain).names)
        assert validate_rows(domain, [row], expected_reference_sha256=REFERENCE) == [row]


def test_deterministic_atomic_parquet_and_roundtrip(tmp_path: Path) -> None:
    row = rows()["mitochondrial_heteroplasmy"]
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    first_report = write_extension(
        "mitochondrial_heteroplasmy", [row], first, expected_reference_sha256=REFERENCE
    )
    second_report = write_extension(
        "mitochondrial_heteroplasmy", [row], second, expected_reference_sha256=REFERENCE
    )
    assert first.read_bytes() == second.read_bytes()
    assert first_report["parquet_sha256"] == second_report["parquet_sha256"]
    assert first_report["canonical_row_multiset"] == second_report["canonical_row_multiset"]
    assert pq.read_table(first).schema == schema_for("mitochondrial_heteroplasmy")
    with pytest.raises(PortableEvidenceError, match="replace"):
        write_extension("mitochondrial_heteroplasmy", [row], first)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"raw_artifact_id": "/" + "raid/private/output.vcf"}, "path-free"),
        ({"limitations": ["inspect /" + "home/person/sample.bam"]}, "private path"),
        ({"confidence_tier": "REVIEW", "confirmation_required": False}, "requires confirmation"),
        ({"reference_sha256": SHA}, "wrong reference"),
    ],
)
def test_fail_closed_common_invariants(mutation: dict[str, object], message: str) -> None:
    row = deepcopy(rows()["structural_variants"])
    row.update(mutation)
    row["evidence_id"] = make_evidence_id("structural_variants", row)
    with pytest.raises(PortableEvidenceError, match=message):
        validate_rows("structural_variants", [row], expected_reference_sha256=REFERENCE)


def test_private_identifier_and_unsalted_hash_are_both_rejected() -> None:
    private_identifier = "PRIVATE_SUBJECT"
    for leaked in (
        private_identifier,
        __import__("hashlib").sha256(private_identifier.encode()).hexdigest(),
    ):
        row = deepcopy(rows()["hla_typing"])
        row["limitations"] = [f"source {leaked}"]
        row["evidence_id"] = make_evidence_id("hla_typing", row)
        with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
            validate_rows("hla_typing", [row], forbidden_bindings=[private_identifier])


def test_short_private_identifier_uses_token_boundaries() -> None:
    assert contains_forbidden_binding("PASS", "A") is False
    assert contains_forbidden_binding("sampleS1suffix", "S1") is False
    assert contains_forbidden_binding("artifact-A", "A") is True
    assert contains_forbidden_binding("artifact_S1", "S1") is True

    row = deepcopy(rows()["structural_variants"])
    assert validate_rows("structural_variants", [row], forbidden_bindings=["A"]) == [row]
    row["raw_artifact_id"] = "caller-output-A"
    row["evidence_id"] = make_evidence_id("structural_variants", row)
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        validate_rows("structural_variants", [row], forbidden_bindings=["A"])


def test_empty_duplicate_and_noncanonical_id_are_rejected() -> None:
    with pytest.raises(PortableEvidenceError, match="empty optional"):
        validate_rows("copy_number", [])
    row = rows()["copy_number"]
    with pytest.raises(PortableEvidenceError, match="duplicate"):
        validate_rows("copy_number", [row, row])
    altered = deepcopy(row)
    altered["copy_number"] = 4.0
    with pytest.raises(PortableEvidenceError, match="canonical row digest"):
        validate_rows("copy_number", [altered])

    required_null = deepcopy(row)
    required_null["caller_name"] = None
    with pytest.raises(PortableEvidenceError, match="required field is null"):
        validate_rows("copy_number", [required_null])

    extra_field = deepcopy(row)
    extra_field["private_subject"] = "forbidden"
    with pytest.raises(PortableEvidenceError, match="physical fields differ"):
        validate_rows("copy_number", [extra_field])


def test_domain_specific_negative_controls() -> None:
    repeat = deepcopy(rows()["repeat_expansions"])
    repeat["allele1_repeat_count"] = 12
    repeat["evidence_id"] = make_evidence_id("repeat_expansions", repeat)
    with pytest.raises(PortableEvidenceError, match="outside"):
        validate_rows("repeat_expansions", [repeat])

    mt = deepcopy(rows()["mitochondrial_heteroplasmy"])
    mt["pos"] = 16_570
    mt["evidence_id"] = make_evidence_id("mitochondrial_heteroplasmy", mt)
    with pytest.raises(PortableEvidenceError, match="out of bounds"):
        validate_rows("mitochondrial_heteroplasmy", [mt])

    pgx = deepcopy(rows()["complex_pharmacogenes"])
    pgx["ambiguous"] = True
    pgx["evidence_id"] = make_evidence_id("complex_pharmacogenes", pgx)
    with pytest.raises(PortableEvidenceError, match="outside call"):
        validate_rows("complex_pharmacogenes", [pgx])
