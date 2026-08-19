"""Explicit Arrow schemas; no inference from nullable production data."""

from __future__ import annotations

import pyarrow as pa


def _dictionary_string() -> pa.DataType:
    return pa.dictionary(pa.int32(), pa.string())


def transcript_annotation_type() -> pa.DataType:
    return pa.struct(
        [
            pa.field("allele", pa.string()),
            pa.field("consequence_terms", pa.list_(pa.string()), nullable=False),
            pa.field("impact", pa.string()),
            pa.field("symbol", pa.string()),
            pa.field("gene_id", pa.string()),
            pa.field("feature_type", pa.string()),
            pa.field("transcript_id", pa.string()),
            pa.field("biotype", pa.string()),
            pa.field("exon", pa.string()),
            pa.field("intron", pa.string()),
            pa.field("hgvsg", pa.string()),
            pa.field("hgvsc", pa.string()),
            pa.field("hgvsp", pa.string()),
            pa.field("cdna_position", pa.string()),
            pa.field("cds_position", pa.string()),
            pa.field("protein_position", pa.string()),
            pa.field("amino_acids", pa.string()),
            pa.field("codons", pa.string()),
            pa.field("existing_variation", pa.list_(pa.string())),
            pa.field("distance", pa.int32()),
            pa.field("strand", pa.int8()),
            pa.field("canonical", pa.bool_()),
            pa.field("mane_select", pa.string()),
            pa.field("mane_plus_clinical", pa.string()),
            pa.field("appris", pa.string()),
            pa.field("transcript_support_level", pa.int8()),
            pa.field("transcript_source", pa.string()),
            pa.field("refseq_transcript_id", pa.string()),
            pa.field("refseq_match", pa.string()),
            pa.field("ccds", pa.string()),
            pa.field("protein_id", pa.string()),
            pa.field("hgnc_id", pa.string()),
            pa.field("symbol_source", pa.string()),
            pa.field("sift_prediction", pa.string()),
            pa.field("sift_score", pa.float64()),
            pa.field("polyphen_prediction", pa.string()),
            pa.field("polyphen_score", pa.float64()),
            pa.field("domains", pa.list_(pa.string())),
            pa.field("lof", pa.string()),
            pa.field("lof_filter", pa.string()),
            pa.field("lof_flags", pa.string()),
            pa.field("lof_info", pa.string()),
            pa.field("flags", pa.list_(pa.string())),
            pa.field("variant_class", pa.string()),
        ]
    )


def variant_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("ref", pa.string(), nullable=False),
            pa.field("alt", pa.string(), nullable=False),
            pa.field("variant_type", _dictionary_string(), nullable=False),
            pa.field("rsid", pa.string()),
            pa.field(
                "quality",
                pa.struct(
                    [
                        pa.field("call_confidence", _dictionary_string(), nullable=False),
                        pa.field("origin", _dictionary_string(), nullable=False),
                        pa.field("imputation_r2", pa.float32()),
                        pa.field("imputation_panel", _dictionary_string()),
                        pa.field("dp", pa.int16()),
                        pa.field("gq", pa.int8()),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "genotype",
                pa.struct(
                    [
                        pa.field("gt", pa.list_(pa.int8()), nullable=False),
                        pa.field("phased", pa.bool_(), nullable=False),
                        pa.field("zygosity", _dictionary_string(), nullable=False),
                        pa.field("allele_balance", pa.float32()),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "gene",
                pa.struct(
                    [
                        pa.field("symbol", _dictionary_string()),
                        pa.field("ensembl_id", pa.string()),
                        pa.field("hgnc_id", pa.string()),
                        pa.field("biotype", _dictionary_string()),
                        pa.field("is_canonical", pa.bool_()),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "consequence",
                pa.struct(
                    [
                        pa.field("most_severe", _dictionary_string()),
                        pa.field("impact", _dictionary_string()),
                        pa.field("hgvsc", pa.string()),
                        pa.field("hgvsp", pa.string()),
                        pa.field("all_annotations", pa.list_(transcript_annotation_type())),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "population_frequency",
                pa.struct(
                    [
                        pa.field("gnomad_af", pa.float32()),
                        pa.field("gnomad_af_popmax", pa.float32()),
                        pa.field("gnomad_af_nfe", pa.float32()),
                        pa.field("gnomad_af_eas", pa.float32()),
                        pa.field("gnomad_af_afr", pa.float32()),
                        pa.field("gnomad_af_sas", pa.float32()),
                        pa.field("gnomad_af_amr", pa.float32()),
                        pa.field("is_rare", pa.bool_(), nullable=False),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "pathogenicity",
                pa.struct(
                    [
                        pa.field("clinvar_significance", _dictionary_string()),
                        pa.field("clinvar_has_conflicts", pa.bool_(), nullable=False),
                        pa.field("clinvar_conflict_summary", _dictionary_string()),
                        pa.field("clinvar_review_stars", pa.int8()),
                        pa.field("clinvar_submitters_count", pa.int16()),
                        pa.field("clinvar_id", pa.string()),
                        pa.field("is_actionable", pa.bool_(), nullable=False),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "predictions",
                pa.struct(
                    [
                        pa.field("cadd_phred", pa.float32()),
                        pa.field("revel", pa.float32()),
                        pa.field("alphamissense", pa.float32()),
                        pa.field("spliceai_max", pa.float32()),
                        pa.field("polyphen", pa.float32()),
                        pa.field("sift", pa.float32()),
                        pa.field("consensus_damaging", pa.bool_(), nullable=False),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "pharmacogenomics",
                pa.struct(
                    [
                        pa.field("is_pgx", pa.bool_(), nullable=False),
                        pa.field("contributes_to_gene", _dictionary_string()),
                    ]
                ),
                nullable=False,
            ),
            pa.field(
                "trait_associations",
                pa.struct(
                    [
                        pa.field("is_gwas_hit", pa.bool_(), nullable=False),
                        pa.field("traits", pa.list_(pa.string())),
                        pa.field("effect_size", pa.float32()),
                        pa.field("effect_type", _dictionary_string()),
                        pa.field("effect_allele", pa.string()),
                        pa.field("ci_lower", pa.float32()),
                        pa.field("ci_upper", pa.float32()),
                        pa.field("p_value", pa.float32()),
                        pa.field("sample_size", pa.int32()),
                        pa.field("study_pmids", pa.list_(pa.string())),
                    ]
                ),
                nullable=False,
            ),
            pa.field("clinical_grade", pa.bool_(), nullable=False),
        ]
    )


def reference_evidence_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("evidence_id", pa.string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("start_pos", pa.int32(), nullable=False),
            pa.field("end_pos", pa.int32(), nullable=False),
            pa.field("state", pa.string(), nullable=False),
            pa.field("reference_observed", pa.bool_()),
            pa.field("callable", pa.bool_()),
            pa.field("call_confidence", pa.string()),
            pa.field("ploidy", pa.int8()),
            pa.field("gt", pa.list_(pa.int16())),
            pa.field("phased", pa.bool_()),
            pa.field("dp", pa.int32()),
            pa.field("minimum_dp", pa.int32()),
            pa.field("gq", pa.int32()),
            pa.field("ad", pa.list_(pa.int32())),
            pa.field("pl", pa.list_(pa.int32())),
            pa.field("record_filters", pa.list_(pa.string())),
            pa.field("genotype_filters", pa.list_(pa.string())),
            pa.field("source_info_json", pa.string()),
            pa.field("other_format_json", pa.string()),
            pa.field("source_record_digest", pa.string(), nullable=False),
            pa.field("source_ordinal", pa.int64(), nullable=False),
            pa.field("source_alts", pa.list_(pa.string())),
            pa.field("reason", pa.string()),
        ]
    )


def non_core_variant_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("evidence_id", pa.string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("ref", pa.string(), nullable=False),
            pa.field("alt", pa.string(), nullable=False),
            pa.field("category", pa.string(), nullable=False),
            pa.field("source_alt_index", pa.int16(), nullable=False),
            pa.field("source_gt", pa.list_(pa.int16()), nullable=False),
            pa.field("phased", pa.bool_(), nullable=False),
            pa.field("dp", pa.int32()),
            pa.field("gq", pa.int32()),
            pa.field("ad", pa.list_(pa.int32())),
            pa.field("pl", pa.list_(pa.int32())),
            pa.field("source_info_json", pa.string()),
            pa.field("other_format_json", pa.string()),
            pa.field("source_record_digest", pa.string(), nullable=False),
            pa.field("source_ordinal", pa.int64(), nullable=False),
        ]
    )


def variant_provenance_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("source_id", pa.string(), nullable=False),
            pa.field("source_record_digest", pa.string(), nullable=False),
            pa.field("source_ordinal", pa.int64(), nullable=False),
            pa.field("source_alt_index", pa.int16(), nullable=False),
            pa.field("source_chrom", pa.string(), nullable=False),
            pa.field("source_pos", pa.int64(), nullable=False),
            pa.field("source_ref", pa.string(), nullable=False),
            pa.field("source_alt", pa.string(), nullable=False),
            pa.field("source_gt", pa.string(), nullable=False),
            pa.field("source_dp", pa.int32()),
            pa.field("source_gq", pa.int32()),
            pa.field("source_ad", pa.list_(pa.int32())),
            pa.field("source_pl", pa.list_(pa.int32())),
            pa.field("source_info_json", pa.string()),
            pa.field("source_other_format_json", pa.string()),
            pa.field("projected_gt", pa.list_(pa.int8()), nullable=False),
            pa.field("phased", pa.bool_(), nullable=False),
            pa.field("normalization_changed", pa.bool_(), nullable=False),
            pa.field("source_filters", pa.list_(pa.string())),
        ]
    )


def callable_region_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("region_id", pa.string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("start_pos", pa.int32(), nullable=False),
            pa.field("end_pos", pa.int32(), nullable=False),
            pa.field("callable", pa.bool_(), nullable=False),
            pa.field("minimum_depth", pa.int32()),
            pa.field("assay_scope", pa.string(), nullable=False),
            pa.field("derivation", pa.string(), nullable=False),
        ]
    )


def callability_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("callability_id", pa.string(), nullable=False),
            pa.field("variant_id", pa.string()),
            pa.field("gene_symbol", pa.string()),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("reference_observed", pa.bool_()),
            pa.field("callable", pa.bool_()),
            pa.field("call_confidence", _dictionary_string()),
            pa.field("depth", pa.int32()),
            pa.field("genotype_quality", pa.float32()),
            pa.field("assay_scope", pa.string(), nullable=False),
            pa.field("evidence_scope", _dictionary_string(), nullable=False),
        ]
    )


def target_catalog_schema() -> pa.Schema:
    """Internal normalized representation of a frozen target-catalogue row."""

    return pa.schema(
        [
            pa.field("callability_id", pa.string(), nullable=False),
            pa.field("target_id", pa.string(), nullable=False),
            pa.field("target_source", pa.string(), nullable=False),
            pa.field("target_source_version", pa.string(), nullable=False),
            pa.field("catalog_sha256", pa.string(), nullable=False),
            pa.field("catalog_row_digest", pa.string(), nullable=False),
            pa.field("catalog_row_number", pa.int64(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("chrom_rank", pa.int8(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("end_pos", pa.int32(), nullable=False),
            pa.field("ref", pa.string()),
            pa.field("alt", pa.string()),
            pa.field("variant_id", pa.string()),
            pa.field("gene_symbol", pa.string()),
            pa.field("assay_scope", pa.string(), nullable=False),
            pa.field("evidence_scope", _dictionary_string(), nullable=False),
        ]
    )


def target_callability_evidence_schema() -> pa.Schema:
    """Lossless target-to-source proof behind normative callability rows."""

    return pa.schema(
        [
            pa.field("callability_id", pa.string(), nullable=False),
            pa.field("target_id", pa.string(), nullable=False),
            pa.field("target_source", pa.string(), nullable=False),
            pa.field("target_source_version", pa.string(), nullable=False),
            pa.field("catalog_sha256", pa.string(), nullable=False),
            pa.field("catalog_row_digest", pa.string(), nullable=False),
            pa.field("catalog_row_number", pa.int64(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("end_pos", pa.int32(), nullable=False),
            pa.field("ref", pa.string()),
            pa.field("alt", pa.string()),
            pa.field("variant_id", pa.string()),
            pa.field("gene_symbol", pa.string()),
            pa.field("resolution", _dictionary_string(), nullable=False),
            pa.field("reference_observed", pa.bool_()),
            pa.field("callable", pa.bool_()),
            pa.field("call_confidence", _dictionary_string()),
            pa.field("depth", pa.int32()),
            pa.field("genotype_quality", pa.float32()),
            pa.field("observed_ploidy", pa.int8()),
            pa.field("assay_scope", pa.string(), nullable=False),
            pa.field("evidence_scope", _dictionary_string(), nullable=False),
            pa.field("observed_variant_ids", pa.list_(pa.string()), nullable=False),
            pa.field("observed_variant_source_digests", pa.list_(pa.string()), nullable=False),
            pa.field("observed_non_core_evidence_ids", pa.list_(pa.string()), nullable=False),
            pa.field("overlapping_no_call_evidence_ids", pa.list_(pa.string()), nullable=False),
            pa.field("selected_reference_evidence_ids", pa.list_(pa.string()), nullable=False),
            pa.field("callable_region_ids", pa.list_(pa.string()), nullable=False),
            pa.field("callable_min_dp", pa.int32(), nullable=False),
            pa.field("callable_min_gq", pa.int32(), nullable=False),
            pa.field("precedence_rule", pa.string(), nullable=False),
            pa.field("reason", pa.string(), nullable=False),
        ]
    )


def gene_index_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("gene_symbol", _dictionary_string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("start_pos", pa.int32(), nullable=False),
            pa.field("end_pos", pa.int32(), nullable=False),
            pa.field("variant_count", pa.int32(), nullable=False),
            pa.field("actionable_count", pa.int32(), nullable=False),
            pa.field("pgx_count", pa.int32(), nullable=False),
        ]
    )


def clinical_finding_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("finding_id", pa.string(), nullable=False),
            pa.field("variant_id", pa.string()),
            pa.field("gene_symbol", pa.string()),
            pa.field("condition", pa.string(), nullable=False),
            pa.field("claim_type", _dictionary_string(), nullable=False),
            pa.field("classification", pa.string()),
            pa.field("clinical_grade", pa.bool_(), nullable=False),
            pa.field("evidence_ids", pa.list_(pa.string()), nullable=False),
        ]
    )


def clinical_evidence_schema() -> pa.Schema:
    """Normative .genome/1.1 clinical-evidence retrieval table."""

    return pa.schema(
        [
            pa.field("evidence_id", pa.string(), nullable=False),
            pa.field("finding_id", pa.string()),
            pa.field("variant_id", pa.string()),
            pa.field("source", pa.string(), nullable=False),
            pa.field("source_record_id", pa.string()),
            pa.field("source_version", pa.string(), nullable=False),
            pa.field("assertion", pa.string(), nullable=False),
            pa.field("review_status", pa.string()),
            pa.field("retrieved_at", pa.timestamp("us", tz="UTC"), nullable=False),
        ]
    )


def interpretation_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("interpretation_id", pa.string(), nullable=False),
            pa.field("finding_id", pa.string()),
            pa.field("variant_id", pa.string()),
            pa.field("gene_symbol", pa.string()),
            pa.field("interpretation_status", _dictionary_string(), nullable=False),
            pa.field("claim_type", _dictionary_string(), nullable=False),
            pa.field("condition", pa.string()),
            pa.field("inheritance_model", pa.string()),
            pa.field("required_genotype", pa.string()),
            pa.field("observed_genotype", pa.string()),
            pa.field("genotype_matches_model", pa.bool_()),
            pa.field("penetrance", pa.string()),
            pa.field("effect_size", pa.float64()),
            pa.field("analytical_confidence", _dictionary_string()),
            pa.field("classification_review_status", pa.string()),
            pa.field("gene_disease_validity", pa.string()),
            pa.field("actionability_basis", pa.string()),
            pa.field("confirmation_required", pa.bool_(), nullable=False),
            pa.field("reviewed_by", pa.string()),
            pa.field("reviewed_at", pa.timestamp("us", tz="UTC")),
            pa.field("evidence_snapshot", pa.string(), nullable=False),
            pa.field("evidence_ids", pa.list_(pa.string()), nullable=False),
        ]
    )


def clinvar_assertion_schema() -> pa.Schema:
    """Lossless ClinVar fields retained in a namespaced extension."""

    return pa.schema(
        [
            pa.field("evidence_id", pa.string(), nullable=False),
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("source", pa.string(), nullable=False),
            pa.field("snapshot", pa.string(), nullable=False),
            pa.field("clinvar_variation_id", pa.string(), nullable=False),
            pa.field("clinvar_accession", pa.string(), nullable=False),
            pa.field("allele_id", pa.int64()),
            pa.field("raw_significance", pa.list_(pa.string()), nullable=False),
            pa.field("aggregate_significance", pa.string()),
            pa.field("review_status", pa.list_(pa.string()), nullable=False),
            pa.field("review_stars", pa.int8(), nullable=False),
            pa.field("has_conflicts", pa.bool_(), nullable=False),
            pa.field("conflict_classifications", pa.list_(pa.string()), nullable=False),
            pa.field("condition_names", pa.list_(pa.string()), nullable=False),
            pa.field("condition_identifiers", pa.list_(pa.string()), nullable=False),
            pa.field("submitter_accessions", pa.list_(pa.string()), nullable=False),
            pa.field("gene_info", pa.list_(pa.string()), nullable=False),
            pa.field("rsids", pa.list_(pa.string()), nullable=False),
            pa.field("origin_codes", pa.list_(pa.int32()), nullable=False),
            pa.field("hgvs", pa.list_(pa.string()), nullable=False),
            pa.field("variation_type", pa.string()),
            pa.field("source_record_digest", pa.string(), nullable=False),
        ]
    )


def acmg_secondary_findings_evidence_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("decision_id", pa.string(), nullable=False),
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("gene_symbol", pa.string(), nullable=False),
            pa.field("policy_key", pa.string(), nullable=False),
            pa.field("disease", pa.string(), nullable=False),
            pa.field("disease_identifiers", pa.list_(pa.string()), nullable=False),
            pa.field("inheritance", pa.string(), nullable=False),
            pa.field("phenotype_category", pa.string()),
            pa.field("variants_to_report", pa.string(), nullable=False),
            pa.field("decision_status", pa.string(), nullable=False),
            pa.field("is_actionable", pa.bool_(), nullable=False),
            pa.field("clinical_grade", pa.bool_(), nullable=False),
            pa.field("confirmation_required", pa.bool_(), nullable=False),
            pa.field("analytical_confidence", pa.string(), nullable=False),
            pa.field("genotype", pa.list_(pa.int8()), nullable=False),
            pa.field("phased", pa.bool_(), nullable=False),
            pa.field("clinvar_significance", pa.string(), nullable=False),
            pa.field("clinvar_snapshot", pa.string(), nullable=False),
            pa.field("clinvar_review_stars", pa.int8(), nullable=False),
            pa.field("clinvar_has_conflicts", pa.bool_(), nullable=False),
            pa.field("clinvar_evidence_ids", pa.list_(pa.string()), nullable=False),
            pa.field("clinvar_accessions", pa.list_(pa.string()), nullable=False),
            pa.field("matched_condition_identifiers", pa.list_(pa.string()), nullable=False),
            pa.field("policy_version", pa.string(), nullable=False),
            pa.field("policy_snapshot", pa.string(), nullable=False),
            pa.field("policy_gene_set_sha256", pa.string(), nullable=False),
            pa.field("guidance", pa.string()),
            pa.field("reasons", pa.list_(pa.string()), nullable=False),
        ]
    )


def prs_schema() -> pa.Schema:
    """Normative .genome/1.1 polygenic-score table."""

    return pa.schema(
        [
            pa.field("trait", pa.string(), nullable=False),
            pa.field("score_value", pa.float32(), nullable=False),
            pa.field("percentile", pa.float32(), nullable=False),
            pa.field("reference_population", _dictionary_string(), nullable=False),
            pa.field("training_source", pa.string(), nullable=False),
            pa.field("training_date", pa.date32(), nullable=False),
            pa.field("contributing_variants", pa.list_(pa.string())),
            pa.field("model_type", _dictionary_string(), nullable=False),
        ]
    )


def pgs_locus_resource_schema() -> pa.Schema:
    """Internal, checksum-pinned PGS scoring-locus representation."""

    return pa.schema(
        [
            pa.field("pgs_id", pa.string(), nullable=False),
            pa.field("scoring_row_number", pa.int32(), nullable=False),
            pa.field("scoring_row_digest", pa.string(), nullable=False),
            pa.field("scoring_file_sha256", pa.string(), nullable=False),
            pa.field("supported", pa.bool_(), nullable=False),
            pa.field("unsupported_reason", pa.string()),
            pa.field("target_id", pa.string()),
            pa.field("chrom", _dictionary_string()),
            pa.field("pos", pa.int32()),
            pa.field("ref", pa.string()),
            pa.field("alt", pa.string()),
            pa.field("variant_id", pa.string()),
            pa.field("effect_allele", pa.string(), nullable=False),
            pa.field("other_allele", pa.string()),
            pa.field("source_effect_allele", pa.string(), nullable=False),
            pa.field("source_other_allele", pa.string()),
            pa.field("strand_flipped", pa.bool_()),
            pa.field("effect_weight", pa.float64(), nullable=False),
            pa.field("effect_is_reference", pa.bool_()),
            pa.field("source_rsid", pa.string()),
            pa.field("source_locus_name", pa.string()),
            pa.field("harmonization_source", pa.string()),
            pa.field("harmonization_rsid", pa.string()),
        ]
    )


def pgs_locus_audit_schema() -> pa.Schema:
    """Per-score-locus evidence and dosage decision without silent imputation."""

    return pa.schema(
        [
            pa.field("pgs_id", pa.string(), nullable=False),
            pa.field("scoring_row_number", pa.int32(), nullable=False),
            pa.field("scoring_row_digest", pa.string(), nullable=False),
            pa.field("target_id", pa.string()),
            pa.field("callability_id", pa.string()),
            pa.field("target_catalog_sha256", pa.string()),
            pa.field("target_catalog_row_digest", pa.string()),
            pa.field("chrom", _dictionary_string()),
            pa.field("pos", pa.int32()),
            pa.field("ref", pa.string()),
            pa.field("alt", pa.string()),
            pa.field("variant_id", pa.string()),
            pa.field("effect_allele", pa.string(), nullable=False),
            pa.field("other_allele", pa.string()),
            pa.field("source_effect_allele", pa.string(), nullable=False),
            pa.field("source_other_allele", pa.string()),
            pa.field("strand_flipped", pa.bool_()),
            pa.field("effect_weight", pa.float64(), nullable=False),
            pa.field("effect_is_reference", pa.bool_()),
            pa.field("resolution", _dictionary_string()),
            pa.field("status", _dictionary_string(), nullable=False),
            pa.field("dosage", pa.int8()),
            pa.field("observed_ploidy", pa.int8()),
            pa.field("contribution", pa.float64()),
            pa.field("observed_variant_ids", pa.list_(pa.string()), nullable=False),
            pa.field("depth", pa.int32()),
            pa.field("genotype_quality", pa.float32()),
            pa.field("reason", pa.string(), nullable=False),
        ]
    )


def pgs_model_audit_schema() -> pa.Schema:
    """Per-model completeness, calibration, and publication decision."""

    return pa.schema(
        [
            pa.field("pgs_id", pa.string(), nullable=False),
            pa.field("trait", pa.string(), nullable=False),
            pa.field("scoring_file_sha256", pa.string(), nullable=False),
            pa.field("scoring_format_version", pa.string(), nullable=False),
            pa.field("source_version", pa.string(), nullable=False),
            pa.field("training_date", pa.date32(), nullable=False),
            pa.field("model_type", _dictionary_string(), nullable=False),
            pa.field("license", pa.string(), nullable=False),
            pa.field("redistribution_allowed", pa.bool_(), nullable=False),
            pa.field("expected_loci", pa.int32(), nullable=False),
            pa.field("supported_loci", pa.int32(), nullable=False),
            pa.field("strand_flipped_loci", pa.int32(), nullable=False),
            pa.field("alternate_called_loci", pa.int32(), nullable=False),
            pa.field("reference_called_loci", pa.int32(), nullable=False),
            pa.field("missing_loci", pa.int32(), nullable=False),
            pa.field("ambiguous_loci", pa.int32(), nullable=False),
            pa.field("low_confidence_loci", pa.int32(), nullable=False),
            pa.field("unsupported_loci", pa.int32(), nullable=False),
            pa.field("expected_absolute_weight", pa.float64(), nullable=False),
            pa.field("used_absolute_weight", pa.float64(), nullable=False),
            pa.field("missing_absolute_weight", pa.float64(), nullable=False),
            pa.field("weight_completeness", pa.float64(), nullable=False),
            pa.field("max_missing_weight_fraction", pa.float64(), nullable=False),
            pa.field("partial_raw_score", pa.float64(), nullable=False),
            pa.field("score_valid", pa.bool_(), nullable=False),
            pa.field("calibrated", pa.bool_(), nullable=False),
            pa.field("percentile", pa.float64()),
            pa.field("reference_population", pa.string()),
            pa.field("calibration_source", pa.string()),
            pa.field("calibration_source_version", pa.string()),
            pa.field("calibration_date", pa.date32()),
            pa.field("calibration_sha256", pa.string()),
            pa.field("calibration_max_missing_weight_fraction", pa.float64()),
            pa.field("publication_status", _dictionary_string(), nullable=False),
            pa.field("reason", pa.string(), nullable=False),
            pa.field("contributing_variant_count", pa.int64(), nullable=False),
            pa.field("contributing_variant_ids_complete", pa.bool_(), nullable=False),
            pa.field("contributing_variant_ids", pa.list_(pa.string())),
            pa.field("target_catalog_sha256", pa.string(), nullable=False),
        ]
    )


def pharmacogenomics_schema() -> pa.Schema:
    """Normative .genome/1.1 gene/diplotype table."""

    return pa.schema(
        [
            pa.field("gene_symbol", _dictionary_string(), nullable=False),
            pa.field("diplotype", pa.string(), nullable=False),
            pa.field("phenotype", _dictionary_string(), nullable=False),
            pa.field("activity_score", pa.float32()),
            pa.field("copy_number", pa.int8()),
            pa.field("cpic_level", _dictionary_string(), nullable=False),
            pa.field("affected_drugs", pa.list_(pa.string()), nullable=False),
            pa.field("guideline_url", pa.string()),
            pa.field("caller", _dictionary_string(), nullable=False),
            pa.field("caller_version", pa.string(), nullable=False),
        ]
    )


def pgx_locus_resource_schema() -> pa.Schema:
    """One allele from the immutable PharmCAT positions resource."""

    return pa.schema(
        [
            pa.field("locus_id", pa.string(), nullable=False),
            pa.field("record_number", pa.int32(), nullable=False),
            pa.field("allele_index", pa.int16(), nullable=False),
            pa.field("source_record_digest", pa.string(), nullable=False),
            pa.field("positions_vcf_sha256", pa.string(), nullable=False),
            pa.field("target_id", pa.string(), nullable=False),
            pa.field("target_variant_id", pa.string(), nullable=False),
            pa.field("normalization_changed", pa.bool_(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("source_id", pa.string()),
            pa.field("ref", pa.string(), nullable=False),
            pa.field("alt", pa.string(), nullable=False),
            pa.field("genes", pa.list_(pa.string()), nullable=False),
            pa.field("position_of_interest", pa.bool_(), nullable=False),
            pa.field("structurally_complex", pa.bool_(), nullable=False),
        ]
    )


def cpic_pair_resource_schema() -> pa.Schema:
    """Pinned CPIC pair_view projection used for normative PGx evidence levels."""

    return pa.schema(
        [
            pa.field("pair_id", pa.int64(), nullable=False),
            pa.field("drug_id", pa.string()),
            pa.field("drug_name", pa.string(), nullable=False),
            pa.field("gene_symbol", pa.string(), nullable=False),
            pa.field("guideline_name", pa.string()),
            pa.field("guideline_url", pa.string()),
            pa.field("cpic_level", pa.string()),
            pa.field("clinpgx_level", pa.string()),
            pa.field("pgx_testing", pa.string()),
            pa.field("pmids", pa.list_(pa.string()), nullable=False),
            # CPIC pair_view uses Yes/No for released rows and n/a for
            # provisional rows.  Keep both a nullable semantic projection and
            # the exact source status so n/a is never silently coerced to No.
            pa.field("used_for_recommendation", pa.bool_()),
            pa.field("used_for_recommendation_status", pa.string(), nullable=False),
            pa.field("provisional", pa.bool_(), nullable=False),
        ]
    )


def pgx_locus_audit_schema() -> pa.Schema:
    """Evidence decision used to materialize one PharmCAT input locus."""

    return pa.schema(
        [
            pa.field("locus_id", pa.string(), nullable=False),
            pa.field("record_number", pa.int32(), nullable=False),
            pa.field("source_record_digest", pa.string(), nullable=False),
            pa.field("positions_vcf_sha256", pa.string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("source_id", pa.string()),
            pa.field("ref", pa.string(), nullable=False),
            pa.field("alts", pa.list_(pa.string()), nullable=False),
            pa.field("genes", pa.list_(pa.string()), nullable=False),
            pa.field("target_ids", pa.list_(pa.string()), nullable=False),
            pa.field("target_variant_ids", pa.list_(pa.string()), nullable=False),
            pa.field("normalization_changed", pa.list_(pa.bool_()), nullable=False),
            pa.field("callability_ids", pa.list_(pa.string()), nullable=False),
            pa.field("evidence_resolutions", pa.list_(pa.string()), nullable=False),
            pa.field("status", _dictionary_string(), nullable=False),
            pa.field("genotype", pa.list_(pa.int8())),
            pa.field("phased", pa.bool_()),
            pa.field("phase_set", pa.string()),
            pa.field("observed_ploidy", pa.int8()),
            pa.field("depth", pa.int32()),
            pa.field("genotype_quality", pa.float32()),
            pa.field("observed_variant_ids", pa.list_(pa.string()), nullable=False),
            pa.field("structurally_complex", pa.bool_(), nullable=False),
            pa.field("eligible_for_matcher", pa.bool_(), nullable=False),
            pa.field("reason", pa.string(), nullable=False),
            pa.field("materialized_record_digest", pa.string(), nullable=False),
        ]
    )


def pgx_gene_call_audit_schema() -> pa.Schema:
    """Lossless-enough PharmCAT call decision preceding normative projection."""

    return pa.schema(
        [
            pa.field("gene_symbol", pa.string(), nullable=False),
            pa.field("source_diplotypes", pa.list_(pa.string()), nullable=False),
            pa.field("selected_diplotype", pa.string()),
            pa.field("phenotypes", pa.list_(pa.string()), nullable=False),
            pa.field("selected_phenotype", pa.string()),
            pa.field("activity_score_raw", pa.string()),
            pa.field("activity_score", pa.float64()),
            pa.field("outside_call", pa.bool_(), nullable=False),
            pa.field("match_score", pa.string()),
            pa.field("missing_positions", pa.list_(pa.string()), nullable=False),
            pa.field("structurally_complex", pa.bool_(), nullable=False),
            pa.field("call_status", _dictionary_string(), nullable=False),
            pa.field("cpic_levels", pa.list_(pa.string()), nullable=False),
            pa.field("affected_drugs", pa.list_(pa.string()), nullable=False),
            pa.field("guideline_urls", pa.list_(pa.string()), nullable=False),
            pa.field("caller", pa.string(), nullable=False),
            pa.field("caller_version", pa.string(), nullable=False),
            pa.field("data_version", pa.string(), nullable=False),
            pa.field("calls_tsv_sha256", pa.string(), nullable=False),
            pa.field("report_json_sha256", pa.string(), nullable=False),
            pa.field("reason", pa.string(), nullable=False),
        ]
    )


def pgx_guidance_audit_schema() -> pa.Schema:
    """Matched PharmCAT guidance retained separately from the compact core table."""

    return pa.schema(
        [
            pa.field("guidance_id", pa.string(), nullable=False),
            pa.field("source", _dictionary_string(), nullable=False),
            pa.field("gene_symbol", pa.string(), nullable=False),
            pa.field("drug", pa.string(), nullable=False),
            pa.field("guideline_id", pa.string()),
            pa.field("guideline_url", pa.string()),
            pa.field("cpic_level", pa.string()),
            pa.field("recommendation", pa.string()),
            pa.field("implications", pa.list_(pa.string()), nullable=False),
            pa.field("classification", pa.string()),
            pa.field("population", pa.string()),
            pa.field("matched", pa.bool_(), nullable=False),
            pa.field("report_json_sha256", pa.string(), nullable=False),
        ]
    )


def family_evidence_schema() -> pa.Schema:
    """Privacy-minimized pedigree evidence for the bundle subject only."""

    return pa.schema(
        [
            pa.field("evidence_id", pa.string(), nullable=False),
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("ref", pa.string(), nullable=False),
            pa.field("alt", pa.string(), nullable=False),
            pa.field("source_version", pa.string(), nullable=False),
            pa.field("snapshot_date", pa.date32(), nullable=False),
            pa.field("joint_vcf_sha256", pa.string(), nullable=False),
            pa.field("pedigree_sha256", pa.string(), nullable=False),
            pa.field("loco_review_sha256", pa.string()),
            pa.field("loco_summary_sha256", pa.string()),
            pa.field("subject_gt", pa.list_(pa.int8()), nullable=False),
            pa.field("subject_phased", pa.bool_(), nullable=False),
            pa.field("phase_set", pa.string()),
            pa.field("phase_confidence", _dictionary_string(), nullable=False),
            pa.field("family_confidence", _dictionary_string()),
            pa.field("family_reasons", pa.list_(pa.string()), nullable=False),
            pa.field("inheritance_state", _dictionary_string(), nullable=False),
            pa.field("parental_origin", _dictionary_string(), nullable=False),
            pa.field("mendelian_status", _dictionary_string(), nullable=False),
            pa.field("loco_status", _dictionary_string(), nullable=False),
            pa.field("loco_source_status", _dictionary_string()),
            pa.field("loco_subject_gt_concordant", pa.bool_()),
            pa.field("loco_exact_allele_in_model", pa.bool_()),
            pa.field("de_novo_score", pa.float32()),
            pa.field("de_novo_candidate", pa.bool_(), nullable=False),
            pa.field("parental_mosaicism_candidate", pa.bool_(), nullable=False),
            pa.field("parental_dropout_suspected", pa.bool_(), nullable=False),
            pa.field("child_dropout_suspected", pa.bool_(), nullable=False),
            pa.field("caller_concordance", _dictionary_string(), nullable=False),
            pa.field("technical_replicate_concordance", _dictionary_string(), nullable=False),
            pa.field("complex_region", pa.bool_(), nullable=False),
            pa.field("difficult_region_flags", pa.list_(pa.string()), nullable=False),
            pa.field("read_evidence_status", _dictionary_string(), nullable=False),
            pa.field("subject_evidence_digest", pa.string(), nullable=False),
        ]
    )


def clingen_gene_disease_schema() -> pa.Schema:
    """Expert-curated gene-disease validity context, never a variant assertion."""

    return pa.schema(
        [
            pa.field("evidence_id", pa.string(), nullable=False),
            pa.field("assertion_id", pa.string(), nullable=False),
            pa.field("gene_symbol", pa.string(), nullable=False),
            pa.field("hgnc_id", pa.string(), nullable=False),
            pa.field("disease_label", pa.string(), nullable=False),
            pa.field("mondo_id", pa.string(), nullable=False),
            pa.field("mode_of_inheritance", pa.string(), nullable=False),
            pa.field("sop", pa.string(), nullable=False),
            pa.field("classification", pa.string(), nullable=False),
            pa.field("online_report", pa.string(), nullable=False),
            pa.field("classification_date", pa.date32(), nullable=False),
            pa.field("classification_datetime_source", pa.string(), nullable=False),
            pa.field("gcep", pa.string(), nullable=False),
            pa.field("snapshot", pa.string(), nullable=False),
            pa.field("source_row_number", pa.int32(), nullable=False),
            pa.field("source_row_sha256", pa.string(), nullable=False),
        ]
    )


def clingen_dosage_sensitivity_schema() -> pa.Schema:
    """GRCh38 gene dosage sensitivity context with lossless PMID source cells."""

    return pa.schema(
        [
            pa.field("evidence_id", pa.string(), nullable=False),
            pa.field("gene_symbol", pa.string(), nullable=False),
            pa.field("ncbi_gene_id", pa.int64(), nullable=False),
            pa.field("cytoband", pa.string(), nullable=False),
            pa.field("chrom", pa.string(), nullable=False),
            pa.field("start", pa.int64(), nullable=False),
            pa.field("end", pa.int64(), nullable=False),
            pa.field("haploinsufficiency_score", pa.int8()),
            pa.field("haploinsufficiency_source_score", pa.string(), nullable=False),
            pa.field("haploinsufficiency_description", pa.string(), nullable=False),
            pa.field("haploinsufficiency_pmids", pa.list_(pa.int64()), nullable=False),
            pa.field(
                "haploinsufficiency_pmid_source_values",
                pa.list_(pa.string()),
                nullable=False,
            ),
            pa.field("haploinsufficiency_disease_id", pa.string()),
            pa.field("triplosensitivity_score", pa.int8()),
            pa.field("triplosensitivity_source_score", pa.string(), nullable=False),
            pa.field("triplosensitivity_description", pa.string(), nullable=False),
            pa.field("triplosensitivity_pmids", pa.list_(pa.int64()), nullable=False),
            pa.field(
                "triplosensitivity_pmid_source_values",
                pa.list_(pa.string()),
                nullable=False,
            ),
            pa.field("triplosensitivity_disease_id", pa.string()),
            pa.field("last_evaluated", pa.date32(), nullable=False),
            pa.field("snapshot", pa.string(), nullable=False),
            pa.field("source_row_number", pa.int32(), nullable=False),
            pa.field("source_row_sha256", pa.string(), nullable=False),
        ]
    )


def gwas_catalog_resource_schema() -> pa.Schema:
    """Normalized immutable NHGRI-EBI GWAS Catalog resource rows."""

    return pa.schema(
        [
            pa.field("chrom", pa.string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("source_row_number", pa.int64(), nullable=False),
            pa.field("risk_rsid", pa.string(), nullable=False),
            pa.field("rsid_aliases", pa.list_(pa.string()), nullable=False),
            pa.field("effect_allele", pa.string(), nullable=False),
            pa.field("reported_trait", pa.string(), nullable=False),
            pa.field("mapped_traits", pa.list_(pa.string()), nullable=False),
            pa.field("mapped_trait_uris", pa.list_(pa.string()), nullable=False),
            pa.field("p_value", pa.float64(), nullable=False),
            pa.field("effect_size", pa.float64()),
            pa.field("effect_type", pa.string()),
            pa.field("effect_type_basis", pa.string()),
            pa.field("ci_lower", pa.float64()),
            pa.field("ci_upper", pa.float64()),
            pa.field("sample_size", pa.int64()),
            pa.field("pmid", pa.string(), nullable=False),
            pa.field("study_accession", pa.string(), nullable=False),
            pa.field("raw_effect_text", pa.string()),
            pa.field("raw_ci_text", pa.string()),
            pa.field("source_record_digest", pa.string(), nullable=False),
        ]
    )


def gwas_association_schema() -> pa.Schema:
    """Normative ``gwas_associations.parquet`` schema from genome-spec 1.1."""

    return pa.schema(
        [
            pa.field("association_id", pa.string(), nullable=False),
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("trait", pa.string(), nullable=False),
            pa.field("effect_allele", pa.string(), nullable=False),
            pa.field("effect_size", pa.float64(), nullable=False),
            pa.field("effect_type", _dictionary_string(), nullable=False),
            pa.field("p_value", pa.float64()),
            pa.field("sample_size", pa.int64()),
            pa.field("study_pmids", pa.list_(pa.string()), nullable=False),
            pa.field("source_version", pa.string(), nullable=False),
        ]
    )


def gwas_catalog_evidence_schema() -> pa.Schema:
    """Lossless exact-match audit, including rows unusable in the compact table."""

    return pa.schema(
        [
            pa.field("association_id", pa.string(), nullable=False),
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("risk_rsid", pa.string(), nullable=False),
            pa.field("rsid_aliases", pa.list_(pa.string()), nullable=False),
            pa.field("effect_allele", pa.string(), nullable=False),
            pa.field("effect_allele_dosage", pa.int8(), nullable=False),
            pa.field("reported_trait", pa.string(), nullable=False),
            pa.field("mapped_traits", pa.list_(pa.string()), nullable=False),
            pa.field("mapped_trait_uris", pa.list_(pa.string()), nullable=False),
            pa.field("p_value", pa.float64(), nullable=False),
            pa.field("effect_size", pa.float64()),
            pa.field("effect_type", _dictionary_string()),
            pa.field("effect_type_basis", pa.string()),
            pa.field("ci_lower", pa.float64()),
            pa.field("ci_upper", pa.float64()),
            pa.field("sample_size", pa.int64()),
            pa.field("study_pmids", pa.list_(pa.string()), nullable=False),
            pa.field("study_accession", pa.string(), nullable=False),
            pa.field("source_version", pa.string(), nullable=False),
            pa.field("source_row_number", pa.int64(), nullable=False),
            pa.field("source_record_digest", pa.string(), nullable=False),
            pa.field("raw_effect_text", pa.string()),
            pa.field("raw_ci_text", pa.string()),
            pa.field("normative_row_emitted", pa.bool_(), nullable=False),
            pa.field("omission_reason", pa.string()),
        ]
    )
