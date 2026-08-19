from __future__ import annotations

import pyarrow as pa

from genome_annotation_builder.arrow_schemas import (
    callability_schema,
    callable_region_schema,
    clingen_dosage_sensitivity_schema,
    clingen_gene_disease_schema,
    clinical_evidence_schema,
    clinical_finding_schema,
    family_evidence_schema,
    gene_index_schema,
    gwas_association_schema,
    gwas_catalog_evidence_schema,
    interpretation_schema,
    pgs_locus_audit_schema,
    pgs_model_audit_schema,
    pgx_gene_call_audit_schema,
    pgx_guidance_audit_schema,
    pgx_locus_audit_schema,
    pharmacogenomics_schema,
    prs_schema,
    target_callability_evidence_schema,
    variant_schema,
)

DICTIONARY_STRING = pa.dictionary(pa.int32(), pa.string())


def _struct_field(schema: pa.Schema, root: str, child: str) -> pa.Field:
    root_field = schema.field(root)
    assert pa.types.is_struct(root_field.type)
    return root_field.type.field(child)


def test_normative_variant_physical_widths_and_dictionary_fields() -> None:
    schema = variant_schema()
    assert schema.field("chrom").type == DICTIONARY_STRING
    assert schema.field("pos").type == pa.int32()
    assert schema.field("variant_type").type == DICTIONARY_STRING
    assert _struct_field(schema, "quality", "call_confidence").type == DICTIONARY_STRING
    assert _struct_field(schema, "quality", "origin").type == DICTIONARY_STRING
    assert _struct_field(schema, "quality", "imputation_r2").type == pa.float32()
    assert _struct_field(schema, "quality", "imputation_panel").type == DICTIONARY_STRING
    assert _struct_field(schema, "quality", "dp").type == pa.int16()
    assert _struct_field(schema, "quality", "gq").type == pa.int8()
    assert _struct_field(schema, "genotype", "zygosity").type == DICTIONARY_STRING
    assert _struct_field(schema, "genotype", "allele_balance").type == pa.float32()
    assert _struct_field(schema, "gene", "symbol").type == DICTIONARY_STRING
    assert _struct_field(schema, "gene", "biotype").type == DICTIONARY_STRING
    assert _struct_field(schema, "consequence", "most_severe").type == DICTIONARY_STRING
    assert _struct_field(schema, "consequence", "impact").type == DICTIONARY_STRING
    for field in (
        "gnomad_af",
        "gnomad_af_popmax",
        "gnomad_af_nfe",
        "gnomad_af_eas",
        "gnomad_af_afr",
        "gnomad_af_sas",
        "gnomad_af_amr",
    ):
        assert _struct_field(schema, "population_frequency", field).type == pa.float32()
    assert _struct_field(schema, "pathogenicity", "clinvar_significance").type == (
        DICTIONARY_STRING
    )
    assert _struct_field(schema, "pathogenicity", "clinvar_submitters_count").type == pa.int16()
    for field in ("cadd_phred", "revel", "alphamissense", "spliceai_max", "polyphen", "sift"):
        assert _struct_field(schema, "predictions", field).type == pa.float32()
    assert _struct_field(schema, "trait_associations", "sample_size").type == pa.int32()


def test_normative_index_callability_and_clinical_schemas() -> None:
    index = gene_index_schema()
    assert index.names == [
        "gene_symbol",
        "chrom",
        "start_pos",
        "end_pos",
        "variant_count",
        "actionable_count",
        "pgx_count",
    ]
    assert index.field("gene_symbol").type == DICTIONARY_STRING
    assert index.field("chrom").type == DICTIONARY_STRING
    for name in index.names[2:]:
        assert index.field(name).type == pa.int32()

    callable_regions = callable_region_schema()
    assert callable_regions.field("chrom").type == DICTIONARY_STRING
    assert callable_regions.field("start_pos").type == pa.int32()
    assert callable_regions.field("end_pos").type == pa.int32()

    callability = callability_schema()
    assert callability.field("chrom").type == DICTIONARY_STRING
    assert callability.field("pos").type == pa.int32()
    assert callability.field("call_confidence").type == DICTIONARY_STRING
    assert callability.field("genotype_quality").type == pa.float32()
    assert callability.field("evidence_scope").type == DICTIONARY_STRING

    target_evidence = target_callability_evidence_schema()
    assert target_evidence.field("chrom").type == DICTIONARY_STRING
    assert target_evidence.field("resolution").type == DICTIONARY_STRING
    assert target_evidence.field("observed_ploidy").type == pa.int8()
    assert target_evidence.field("catalog_row_number").type == pa.int64()

    findings = clinical_finding_schema()
    assert findings.field("claim_type").type == DICTIONARY_STRING
    evidence = clinical_evidence_schema()
    assert evidence.field("retrieved_at").type == pa.timestamp("us", tz="UTC")
    interpretations = interpretation_schema()
    assert interpretations.field("interpretation_status").type == DICTIONARY_STRING
    assert interpretations.field("claim_type").type == DICTIONARY_STRING
    assert interpretations.field("analytical_confidence").type == DICTIONARY_STRING

    prs = prs_schema()
    assert prs.field("score_value").type == pa.float32()
    assert prs.field("percentile").type == pa.float32()
    assert prs.field("reference_population").type == DICTIONARY_STRING
    assert prs.field("training_date").type == pa.date32()
    assert prs.field("model_type").type == DICTIONARY_STRING

    model_audit = pgs_model_audit_schema()
    assert model_audit.field("publication_status").type == DICTIONARY_STRING
    assert model_audit.field("training_date").type == pa.date32()
    assert model_audit.field("contributing_variant_count").type == pa.int64()
    assert model_audit.field("contributing_variant_ids").nullable is True
    locus_audit = pgs_locus_audit_schema()
    assert locus_audit.field("status").type == DICTIONARY_STRING
    assert locus_audit.field("dosage").type == pa.int8()

    pharmacogenomics = pharmacogenomics_schema()
    assert pharmacogenomics.field("gene_symbol").type == DICTIONARY_STRING
    assert pharmacogenomics.field("activity_score").type == pa.float32()
    assert pharmacogenomics.field("copy_number").type == pa.int8()
    assert pharmacogenomics.field("cpic_level").type == DICTIONARY_STRING

    pgx_locus = pgx_locus_audit_schema()
    assert pgx_locus.field("status").type == DICTIONARY_STRING
    assert pgx_locus.field("genotype").type == pa.list_(pa.int8())
    assert pgx_locus.field("observed_ploidy").type == pa.int8()
    pgx_gene = pgx_gene_call_audit_schema()
    assert pgx_gene.field("call_status").type == DICTIONARY_STRING
    assert pgx_gene.field("activity_score").type == pa.float64()
    pgx_guidance = pgx_guidance_audit_schema()
    assert pgx_guidance.field("source").type == DICTIONARY_STRING

    family = family_evidence_schema()
    assert family.field("chrom").type == DICTIONARY_STRING
    assert family.field("snapshot_date").type == pa.date32()
    assert family.field("subject_gt").type == pa.list_(pa.int8())
    assert family.field("phase_confidence").type == DICTIONARY_STRING
    assert family.field("mendelian_status").type == DICTIONARY_STRING
    assert {"father_gt", "mother_gt", "relative_identifier"}.isdisjoint(family.names)

    gwas = gwas_association_schema()
    assert gwas.field("effect_size").type == pa.float64()
    assert gwas.field("effect_type").type == DICTIONARY_STRING
    assert gwas.field("sample_size").type == pa.int64()
    gwas_audit = gwas_catalog_evidence_schema()
    assert gwas_audit.field("effect_allele_dosage").type == pa.int8()
    assert gwas_audit.field("effect_type").type == DICTIONARY_STRING

    clingen_gene = clingen_gene_disease_schema()
    assert clingen_gene.field("classification_date").type == pa.date32()
    assert clingen_gene.field("source_row_number").type == pa.int32()
    assert clingen_gene.field("gene_symbol").nullable is False
    clingen_dosage = clingen_dosage_sensitivity_schema()
    assert clingen_dosage.field("start").type == pa.int64()
    assert clingen_dosage.field("haploinsufficiency_score").type == pa.int8()
    assert clingen_dosage.field("haploinsufficiency_score").nullable is True
    assert clingen_dosage.field("haploinsufficiency_pmids").type == pa.list_(pa.int64())
