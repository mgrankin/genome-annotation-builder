"""Atomic assembly of a self-describing .genome/1.1 bundle."""

from __future__ import annotations

import shutil
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .config import BuildConfig
from .constants import (
    ACMG_SF_GENE_SET_SHA256,
    GENOME_SPEC_COMMIT,
    GENOME_SPEC_REPOSITORY,
    PHARMCAT_VERSION,
    PINNED_UPSTREAM_FILES,
    PIPELINE_VERSION,
    SCHEMA_VERSION,
    STRUCTURALLY_COMPLEX_PGX_GENES,
)
from .profile_assets import QUERY_CATALOG, README, SKILL, SYNONYMS, canonical_questions
from .resources import resource_path
from .util import (
    atomic_output_directory,
    atomic_write_bytes,
    atomic_write_json,
    directory_digest,
    sha256_file,
    sha256_json,
)

EMPTY_CANONICAL_LIST_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


def _copy_tree(source: Path, destination: Path) -> None:
    for item in source.rglob("*"):
        if item.is_symlink():
            raise ValueError(f"refusing to copy symbolic link into bundle: {item}")
    shutil.copytree(source, destination)


def _directory_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _parquet_rows(path: Path) -> int:
    files = [path] if path.is_file() else sorted(path.rglob("*.parquet"))
    return sum(pq.ParquetFile(item).metadata.num_rows for item in files)


def _first_parquet_schema(path: Path) -> pa.Schema:
    files = [path] if path.is_file() else sorted(path.rglob("*.parquet"))
    if not files:
        raise ValueError(f"Parquet dataset has no files: {path}")
    return pq.read_schema(files[0])


def _predictor_semantic_control_stage(
    predictor_resource: dict[str, Any], configuration_digest: str
) -> dict[str, Any]:
    """Build the transformation row from the same immutable control boundary."""

    semantic_controls = predictor_resource["semantic_controls"]
    return {
        "stage_id": "offline_predictor_semantic_controls",
        "software": "genome-annotation-builder",
        "software_version": PIPELINE_VERSION,
        "configuration_digest": configuration_digest,
        "inputs": [
            "sha256:" + predictor_resource["manifest_sha256"],
            "sha256:" + semantic_controls["control_definition_sha256"],
        ],
        "outputs": ["validated offline predictor runtime boundary"],
        "control_set": semantic_controls["control_set"],
        "control_definition_sha256": semantic_controls["control_definition_sha256"],
        "doctor_run_key": semantic_controls["doctor_run_key"],
        "report_sha256": semantic_controls["report_sha256"],
        "enabled_plugins": semantic_controls["enabled_plugins"],
        "record_count": semantic_controls["record_count"],
        "record_identity_sha256": semantic_controls["record_identity_sha256"],
        "reference_fasta_sha256": semantic_controls["reference_fasta_sha256"],
        "vep_cache_manifest_sha256": semantic_controls["vep_cache_manifest_sha256"],
        "control_results": semantic_controls["control_results"],
        "semantic_boundary": semantic_controls["semantic_boundary"],
        "status": "completed",
    }


ALLOWED_VALUES: dict[str, list[Any]] = {
    "variants.variant_type": ["SNV", "INS", "DEL", "MNV", "SV"],
    "variants.quality.call_confidence": ["high", "medium", "low"],
    "variants.quality.origin": ["observed", "imputed"],
    "variants.genotype.zygosity": ["hom_ref", "het", "hom_alt"],
    "variants.consequence.impact": ["HIGH", "MODERATE", "LOW", "MODIFIER"],
    "variants.pathogenicity.clinvar_significance": [
        "Pathogenic",
        "Likely_pathogenic",
        "Uncertain_significance",
        "Likely_benign",
        "Benign",
        "Conflicting_interpretations",
    ],
    "variants.trait_associations.effect_type": ["beta", "odds_ratio", "hazard_ratio"],
    "reference_confidence_regions.state": [
        "reference_block",
        "explicit_reference_site",
        "no_call",
    ],
    "clinical_findings.claim_type": [
        "variant_classification",
        "carrier_status",
        "disease_risk",
        "pharmacogenomic",
        "trait_association",
        "computational_prediction",
    ],
    "interpretations.interpretation_status": [
        "algorithmic",
        "draft",
        "expert_reviewed",
        "guideline_based",
    ],
    "interpretations.claim_type": [
        "variant_classification",
        "carrier_status",
        "disease_risk",
        "pharmacogenomic",
        "trait_association",
        "computational_prediction",
    ],
    "interpretations.analytical_confidence": ["high", "medium", "low"],
    "callability.call_confidence": ["high", "medium", "low"],
    "callability.evidence_scope": ["site", "target_panel", "derived_interval"],
    "target_callability_evidence.resolution": [
        "non_reference_observed",
        "overlapping_variant_observed",
        "no_call",
        "reference_observed",
        "low_confidence_reference",
        "callable_without_reference",
        "uncallable_without_reference",
        "unknown",
    ],
    "clingen_gene_disease.classification": [
        "Definitive",
        "Strong",
        "Moderate",
        "Limited",
        "Disputed",
        "Refuted",
        "No Known Disease Relationship",
    ],
    "clingen_dosage_sensitivity.haploinsufficiency_score": [0, 1, 2, 3, 30, 40],
    "clingen_dosage_sensitivity.triplosensitivity_score": [0, 1, 2, 3, 30, 40],
    "target_callability_evidence.call_confidence": ["high", "medium", "low"],
    "target_callability_evidence.evidence_scope": [
        "site",
        "target_panel",
        "derived_interval",
    ],
    "pgs_locus_audit.status": [
        "alternate_called",
        "reference_called",
        "missing_no_call",
        "missing_interval_only",
        "missing_unknown",
        "ambiguous_other_variant",
        "ambiguous_genotype",
        "low_confidence",
        "unsupported",
    ],
    "pgs_model_audit.model_type": ["linear"],
    "pgs_model_audit.publication_status": [
        "published",
        "valid_uncalibrated",
        "calibration_inapplicable",
        "incomplete",
    ],
    "prs.model_type": ["linear"],
    "pharmacogenomics.cpic_level": ["A", "B", "C", "D"],
    "pgx_locus_audit.status": [
        "variant_call",
        "reference_call",
        "no_call_structurally_complex",
        "no_call_other_variant",
        "no_call_low_confidence_variant",
        "no_call_inconsistent_genotype",
        "no_call_low_confidence_reference",
        "no_call_source_missing",
        "no_call_interval_only",
        "no_call_inconsistent_reference",
        "no_call_unknown",
    ],
    "pgx_gene_call_audit.call_status": [
        "published",
        "complex_requires_validated_outside_call",
        "ambiguous_diplotype",
        "incomplete_positions",
        "phenotype_unavailable",
        "no_normative_cpic_level",
    ],
    "family_evidence.phase_confidence": [
        "phased_with_phase_set",
        "phased_without_phase_set",
        "unphased",
    ],
    "family_evidence.mendelian_status": [
        "consistent",
        "conflict",
        "parent_or_subject_no_call",
        "not_assessed",
        "not_applicable_founder",
    ],
    "family_evidence.caller_concordance": [
        "exact_genotype_concordant",
        "discordant",
        "not_assessed",
    ],
    "family_evidence.technical_replicate_concordance": [
        "both_sources_exact_genotype_concordant",
        "single_source_support",
        "discordant",
        "not_assessed",
    ],
}


def _glossary_entry(
    table: str, field_path: str, field: pa.Field, config: BuildConfig
) -> dict[str, Any]:
    key = f"{table}.{field_path}"
    annotation_prefixes = (
        "gene",
        "consequence",
        "population_frequency",
        "pathogenicity",
        "predictions",
        "pharmacogenomics",
        "trait_associations",
        "clinical_grade",
    )
    first = field_path.split(".", maxsplit=1)[0]
    vep_field = first in {"gene", "consequence", "population_frequency"} or field_path in {
        "predictions.polyphen",
        "predictions.sift",
    }
    if table == "family_evidence":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "digest-pinned family joint VCF, PED, and optional child-independent LOCO audit"
        derivation = (
            "exact subject genotype reconciliation plus privacy-minimized inheritance, phasing, "
            "concordance, and review statuses; no relative genotype is persisted"
        )
    elif table == "pharmacogenomics":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "digest-pinned PharmCAT and frozen CPIC pair snapshot"
        derivation = (
            "complete unambiguous non-complex PharmCAT diplotype projected only when an exact "
            "non-provisional CPIC A-D recommendation mapping exists"
        )
    elif table == "pgx_locus_audit":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "official PharmCAT positions joined to target-level VCF evidence"
        derivation = (
            "exact carried genotype or complete explicit-reference evidence; no-call, interval, "
            "other-ALT, low-confidence, and structurally complex loci remain non-reference claims"
        )
    elif table == "pgx_gene_call_audit":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "offline digest-pinned PharmCAT calls report and frozen CPIC pair snapshot"
        derivation = (
            "gene-call publication gate retaining ambiguity, missingness, and exclusion reason"
        )
    elif table == "pgx_guidance_audit":
        evidence_class = "external_annotation"
        observed_or_inferred = "inferred"
        source = "offline digest-pinned PharmCAT report and frozen CPIC guideline mapping"
        derivation = "lossless structured guidance projection linked to the gated gene-call audit"
    elif table == "prs":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "checksum-pinned PGS Catalog model and explicit calibration resource"
        derivation = (
            "allele-oriented additive scoring from explicit carried/reference genotypes; "
            "missing and ambiguous loci are never converted to reference"
        )
    elif table == "pgs_model_audit":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "checksum-pinned PGS Catalog model, locus audit, and optional calibration"
        derivation = (
            "model-level expected/used weight, missingness, calibration, and publication gate"
        )
    elif table == "pgs_locus_audit":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "PGS scoring row joined to target-level reference-confidence proof"
        derivation = (
            "exact GRCh38 allele orientation and variant/reference/no-call/unknown dosage decision"
        )
    elif table in {"clinical_findings", "interpretations"}:
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "checksum-pinned evidence and policy resources declared in manifest.json"
        derivation = (
            "algorithmic retrieval/interpretation record; clinical_grade is false and "
            "confirmation remains required unless separately reviewed"
        )
    elif table == "acmg_secondary_findings_evidence":
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "checksum-pinned ACMG SF and ClinGen policy plus exact-allele ClinVar evidence"
        derivation = (
            "conservative disease, review, conflict, genotype, and gene-specific policy gates; "
            "clinical confirmation remains required"
        )
    elif table == "clinical_evidence":
        evidence_class = "external_annotation"
        observed_or_inferred = "observed"
        source = "checksum-pinned evidence snapshots declared in manifest.json"
        derivation = (
            "normative retrieval row retaining assertion, version, review, and snapshot time"
        )
    elif table == "clinvar_assertions":
        evidence_class = "external_annotation"
        observed_or_inferred = "observed"
        source = "checksum-pinned exact-allele NCBI ClinVar snapshot declared in manifest.json"
        derivation = "lossless typed extraction without collapsing conflicts or review status"
    elif table in {"clingen_gene_disease", "clingen_dosage_sensitivity"}:
        evidence_class = "external_annotation"
        observed_or_inferred = "observed"
        source = "checksum-pinned CC0 ClinGen gene-level curation snapshot"
        derivation = (
            "exact HGNC gene-symbol projection; gene-level context only and never a "
            "variant-pathogenicity or subject-diagnosis assertion"
        )
    elif table in {"gwas_associations", "gwas_catalog_evidence"}:
        evidence_class = "external_annotation"
        observed_or_inferred = "observed"
        source = "checksum-pinned dated NHGRI-EBI GWAS Catalog release"
        derivation = (
            "GRCh38 coordinate plus exact rsID alias and exact effect allele matching; "
            "ambiguous legacy OR-or-BETA rows remain explicitly untyped"
        )
    elif table == "variants" and config.clinvar_enabled and first == "pathogenicity":
        evidence_class = "external_annotation"
        observed_or_inferred = "inferred"
        source = "checksum-pinned exact-allele NCBI ClinVar snapshot declared in manifest.json"
        derivation = (
            "documented aggregation backed by clinical_evidence.parquet and "
            "extensions/clinvar_assertions.parquet"
        )
    elif table in {"variants", "observed_partial_variants"} and config.vep_enabled and vep_field:
        evidence_class = "external_annotation"
        observed_or_inferred = "inferred"
        source = "checksum-pinned offline Ensembl VEP/cache declared in manifest.json"
        derivation = "typed parsing of the exact-allele VEP CSQ annotation"
    elif table == "variants" and config.gwas.enabled and first == "trait_associations":
        evidence_class = "external_annotation"
        observed_or_inferred = "inferred"
        source = "checksum-pinned dated NHGRI-EBI GWAS Catalog release"
        derivation = (
            "compact deterministic summary of exact matches backed by "
            "gwas_associations.parquet and extensions/gwas_catalog_evidence.parquet"
        )
    elif table == "variants" and first in annotation_prefixes:
        evidence_class = "external_annotation"
        observed_or_inferred = "inferred"
        source = "annotation adapters; absent in Release-A core"
        derivation = "null/false unless a checksum-pinned adapter populated the field"
    elif table in {"callability", "target_callability_evidence"}:
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "frozen target catalogue plus retained locus-specific VCF evidence"
        derivation = (
            "carried-variant > no-call > complete-reference > interval-only > unknown "
            "precedence with explicit DP/GQ and ploidy provenance"
        )
    elif table in {"callable_regions"}:
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "reference-confidence VCF record metrics"
        derivation = "configured FILTER/FT/DP/GQ thresholds"
    elif table == "observed_partial_variants" and (
        first.startswith("source_") or first in {"gt", "phased"}
    ):
        evidence_class = "observed"
        observed_or_inferred = "observed"
        source = "immutable incomplete source genotype"
        derivation = "known selected ALT with original missing slot and phasing preserved"
    elif table in {"variant_provenance", "reference_confidence_regions", "non_core_variants"}:
        evidence_class = "observed"
        observed_or_inferred = "observed"
        source = "immutable source VCF record"
        derivation = "lossless extraction or identifier-preserving normalization"
    elif table == "variants" and first in {
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_type",
    }:
        evidence_class = "pipeline_derived"
        observed_or_inferred = "inferred"
        source = "carried concrete ALT from the source VCF"
        derivation = "biallelic projection and bcftools normalization against GRCh38"
    else:
        evidence_class = "observed"
        observed_or_inferred = "observed"
        source = "selected sample fields in the immutable source VCF"
        derivation = "typed extraction without missing-to-reference coercion"
    return {
        "type": str(field.type),
        "nullable": field.nullable,
        "null_meaning": (
            "source value was missing or exceeded the normative physical width; the exact "
            "value remains in extensions/variant_provenance.parquet"
            if key in {"variants.quality.dp", "variants.quality.gq"}
            else (
                "value unavailable, not applicable, or not produced; never a benign/reference claim"
                if field.nullable
                else "field is required and cannot be null"
            )
        ),
        "source": source,
        "derivation": derivation,
        "allowed_values": ALLOWED_VALUES.get(key),
        "evidence_class": evidence_class,
        "observed_or_inferred": observed_or_inferred,
        "genome_build": "GRCh38" if table != "variant_provenance" or first != "source_id" else None,
        "version_or_snapshot": PIPELINE_VERSION,
        "confidence_semantics": (
            "missing AF is conservatively retained; true does not prove observed rarity"
            if key == "variants.population_frequency.is_rare"
            else (
                "technical thresholds declared in manifest; not clinical confidence"
                if "confidence" in field_path or field_path in {"callable", "quality"}
                else "no additional confidence semantics beyond the field definition"
            )
        ),
        "absence_meaning": (
            "no matching row is only a retrieval outcome; inspect callability and omissions"
        ),
        "related_evidence": [
            "manifest.json",
            "omissions.json",
            "extensions/reconciliation.json",
        ],
    }


def _schema_fields(table: str, schema: pa.Schema, config: BuildConfig) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}

    def visit(field: pa.Field, prefix: str) -> None:
        output[f"{table}.{prefix}"] = _glossary_entry(table, prefix, field, config)
        if pa.types.is_struct(field.type):
            for child in field.type:
                visit(child, f"{prefix}.{child.name}")

    for root_field in schema:
        visit(root_field, root_field.name)
    return output


def _configuration_digest(config: BuildConfig) -> str:
    return sha256_json(
        {
            "genome_build": config.genome_build,
            "genome_spec_commit": config.genome_spec_commit,
            "normalization": {"check_ref": config.check_ref},
            "quality": vars(config.quality),
            "annotation": {
                "vep_enabled": config.vep_enabled,
                "vep_cache_version": config.vep_cache_version,
                "predictors": {
                    name: {
                        "enabled": getattr(config.predictors, name).enabled,
                        "license": getattr(config.predictors, name).license,
                        "resource_redistribution_allowed": getattr(
                            config.predictors, name
                        ).redistribution_allowed,
                    }
                    for name in (
                        "loftee",
                        "cadd",
                        "revel",
                        "alphamissense",
                        "spliceai",
                    )
                },
                "clinvar_enabled": config.clinvar_enabled,
                "clinvar_snapshot": config.clinvar_snapshot,
                "acmg_sf_enabled": config.acmg_sf_enabled,
                "acmg_sf_version": config.acmg_sf_version,
                "acmg_sf_snapshot": config.acmg_sf_snapshot,
                "acmg_sf_min_clinvar_review_stars": (config.acmg_sf_min_clinvar_review_stars),
                "pgs_enabled": config.pgs_enabled,
                "pgs_models": [
                    {
                        "pgs_id": model.pgs_id,
                        "scoring_file_sha256": model.expected_sha256,
                        "source_version": model.source_version,
                        "training_date": model.training_date.isoformat(),
                        "trait": model.trait,
                        "model_type": model.model_type,
                        "license": model.license,
                        "redistribution_allowed": model.redistribution_allowed,
                        "max_missing_weight_fraction": model.max_missing_weight_fraction,
                        "calibration_sha256": (
                            model.calibration.expected_sha256 if model.calibration else None
                        ),
                    }
                    for model in config.pgs_models
                ],
                "pgx": {
                    "enabled": config.pgx.enabled,
                    "pharmcat_version": config.pgx.version,
                    "pharmcat_release_commit": config.pgx.release_commit,
                    "positions_vcf_sha256": config.pgx.positions_expected_sha256,
                    "container_image": config.pgx.container_image,
                    "cpic_pairs_sha256": config.pgx.cpic_pairs_expected_sha256,
                    "cpic_source_version": config.pgx.cpic_source_version,
                    "cpic_snapshot_date": config.pgx.cpic_snapshot_date.isoformat(),
                },
                "family_evidence": {
                    "enabled": config.family_evidence.enabled,
                    "joint_vcf_sha256": config.family_evidence.joint_vcf_expected_sha256,
                    "pedigree_sha256": config.family_evidence.pedigree_expected_sha256,
                    "source_version": config.family_evidence.source_version,
                    "snapshot_date": config.family_evidence.snapshot_date.isoformat(),
                    "tags": vars(config.family_evidence.tags),
                    "loco_review_sha256": (
                        config.family_evidence.loco.review_expected_sha256
                        if config.family_evidence.loco is not None
                        else None
                    ),
                    "loco_summary_sha256": (
                        config.family_evidence.loco.summary_expected_sha256
                        if config.family_evidence.loco is not None
                        else None
                    ),
                },
                "gwas_catalog": {
                    "enabled": config.gwas.enabled,
                    "snapshot": config.gwas.snapshot,
                    "archive_sha256": config.gwas.archive_expected_sha256,
                    "license": config.gwas.license,
                    "redistribution_allowed": config.gwas.redistribution_allowed,
                },
                "target_catalog_enabled": config.target_catalog_enabled,
                "target_catalog_expected_sha256": config.target_catalog_expected_sha256,
                "target_catalog_assay_scope": config.target_catalog_assay_scope,
            },
            "output": {
                "agent_exploration_profile": config.agent_exploration_profile,
                "reference_confidence_extension": config.reference_confidence_extension,
            },
        }
    )


def _write_profile(staging: Path) -> None:
    atomic_write_bytes(staging / "readmygenome.md", README.encode("utf-8"))
    profile = staging / "readmygenome"
    atomic_write_bytes(profile / "SKILL.md", SKILL.encode("utf-8"))
    atomic_write_bytes(profile / "reference" / "query_catalog.md", QUERY_CATALOG.encode("utf-8"))
    atomic_write_json(profile / "reference" / "synonyms.json", SYNONYMS)
    atomic_write_json(profile / "evals" / "canonical_questions.json", canonical_questions())
    query_source = Path(__file__).with_name("query.py")
    atomic_write_bytes(profile / "scripts" / "query.py", query_source.read_bytes(), mode=0o755)


def _inventory_entry(
    root: Path,
    relative: str,
    media_type: str,
    logical_role: str,
    required: bool,
    authority: str,
    rows: int | None = None,
) -> dict[str, Any]:
    path = root / relative
    entry: dict[str, Any] = {
        "path": relative,
        "media_type": media_type,
        "logical_role": logical_role,
        "required": required,
        "authority": authority,
        "rows": rows,
    }
    if path.exists():
        entry["bytes"] = _directory_bytes(path) if path.is_dir() else path.stat().st_size
        entry["sha256"] = directory_digest(path) if path.is_dir() else sha256_file(path)
    else:
        entry["bytes"] = None
        entry["sha256"] = None
    return entry


def assemble_bundle(
    config: BuildConfig,
    *,
    run_key: str,
    preflight: dict[str, Any],
    ingest_dir: Path,
    ingest_summary: dict[str, Any],
    normalization_summary: dict[str, Any],
    vep_summary: dict[str, Any] | None,
    variant_dir: Path,
    variant_summary: dict[str, Any],
    target_callability_dir: Path | None,
    target_callability_summary: dict[str, Any] | None,
    pgs_score_dir: Path | None,
    pgs_score_summary: dict[str, Any] | None,
    pgx_materialization_dir: Path | None,
    pgx_materialization_summary: dict[str, Any] | None,
    pharmcat_dir: Path | None,
    pharmcat_summary: dict[str, Any] | None,
    pgx_projection_dir: Path | None,
    pgx_projection_summary: dict[str, Any] | None,
    family_evidence_dir: Path | None,
    family_evidence_summary: dict[str, Any] | None,
) -> dict[str, Any]:
    pgx_parts = (
        pgx_materialization_dir,
        pgx_materialization_summary,
        pharmcat_dir,
        pharmcat_summary,
        pgx_projection_dir,
        pgx_projection_summary,
    )
    if any(value is not None for value in pgx_parts) and not all(
        value is not None for value in pgx_parts
    ):
        raise ValueError("PGx bundle inputs must be supplied as one complete provenance chain")
    if (family_evidence_dir is None) != (family_evidence_summary is None):
        raise ValueError("family evidence bundle inputs must be supplied together")
    generated_at = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    with atomic_output_directory(config.bundle_path, overwrite=config.overwrite) as staging:
        schema_source = resource_path(config, "schema.json")
        if sha256_file(schema_source) != PINNED_UPSTREAM_FILES["schema.json"]:
            raise ValueError("normative schema digest changed before bundle assembly")
        atomic_write_bytes(staging / "schema.json", schema_source.read_bytes())

        _copy_tree(variant_dir / "variants.parquet", staging / "variants.parquet")
        gene_index_source = variant_dir / "gene_index.parquet"
        if gene_index_source.is_file():
            atomic_write_bytes(staging / "gene_index.parquet", gene_index_source.read_bytes())
        clinical_evidence_source = variant_dir / "clinical_evidence.parquet"
        if clinical_evidence_source.is_dir():
            _copy_tree(clinical_evidence_source, staging / "clinical_evidence.parquet")
        clinical_findings_source = variant_dir / "clinical_findings.parquet"
        if clinical_findings_source.is_file():
            atomic_write_bytes(
                staging / "clinical_findings.parquet", clinical_findings_source.read_bytes()
            )
        interpretations_source = variant_dir / "interpretations.parquet"
        if interpretations_source.is_file():
            atomic_write_bytes(
                staging / "interpretations.parquet", interpretations_source.read_bytes()
            )
        gwas_associations_source = variant_dir / "gwas_associations.parquet"
        if gwas_associations_source.is_dir():
            _copy_tree(gwas_associations_source, staging / "gwas_associations.parquet")
        gwas_evidence_source = variant_dir / "gwas_catalog_evidence.parquet"
        clinvar_assertions_source = variant_dir / "clinvar_assertions.parquet"
        acmg_evidence_source = variant_dir / "acmg_secondary_findings_evidence.parquet"
        clingen_gene_disease_source = variant_dir / "clingen_gene_disease.parquet"
        clingen_dosage_source = variant_dir / "clingen_dosage_sensitivity.parquet"
        callability_source = (
            target_callability_dir / "callability.parquet"
            if target_callability_dir is not None
            else None
        )
        target_callability_evidence_source = (
            target_callability_dir / "target_callability_evidence.parquet"
            if target_callability_dir is not None
            else None
        )
        if callability_source is not None and callability_source.is_file():
            atomic_write_bytes(staging / "callability.parquet", callability_source.read_bytes())
        prs_source = pgs_score_dir / "prs.parquet" if pgs_score_dir is not None else None
        pgs_model_audit_source = (
            pgs_score_dir / "pgs_model_audit.parquet" if pgs_score_dir is not None else None
        )
        pgs_locus_audit_source = (
            pgs_score_dir / "pgs_locus_audit.parquet" if pgs_score_dir is not None else None
        )
        if prs_source is not None and prs_source.is_file():
            atomic_write_bytes(staging / "prs.parquet", prs_source.read_bytes())
        pharmacogenomics_source = (
            pgx_projection_dir / "pharmacogenomics.parquet"
            if pgx_projection_dir is not None
            else None
        )
        pgx_locus_audit_source = (
            pgx_materialization_dir / "pgx_locus_audit.parquet"
            if pgx_materialization_dir is not None
            else None
        )
        pgx_gene_call_audit_source = (
            pgx_projection_dir / "pgx_gene_call_audit.parquet"
            if pgx_projection_dir is not None
            else None
        )
        pgx_guidance_audit_source = (
            pgx_projection_dir / "pgx_guidance_audit.parquet"
            if pgx_projection_dir is not None
            else None
        )
        family_evidence_source = (
            family_evidence_dir / "family_evidence.parquet"
            if family_evidence_dir is not None
            else None
        )
        if pharmacogenomics_source is not None and pharmacogenomics_source.is_file():
            atomic_write_bytes(
                staging / "pharmacogenomics.parquet", pharmacogenomics_source.read_bytes()
            )
        _copy_tree(ingest_dir / "callable_regions.parquet", staging / "callable_regions.parquet")
        extensions = staging / "extensions"
        extensions.mkdir(parents=True, exist_ok=False)
        if clinvar_assertions_source.is_dir():
            _copy_tree(
                clinvar_assertions_source,
                extensions / "clinvar_assertions.parquet",
            )
        if acmg_evidence_source.is_file():
            atomic_write_bytes(
                extensions / "acmg_secondary_findings_evidence.parquet",
                acmg_evidence_source.read_bytes(),
            )
        if (
            target_callability_evidence_source is not None
            and target_callability_evidence_source.is_file()
        ):
            atomic_write_bytes(
                extensions / "target_callability_evidence.parquet",
                target_callability_evidence_source.read_bytes(),
            )
        if pgs_model_audit_source is not None and pgs_model_audit_source.is_file():
            atomic_write_bytes(
                extensions / "pgs_model_audit.parquet", pgs_model_audit_source.read_bytes()
            )
        if pgs_locus_audit_source is not None and pgs_locus_audit_source.is_file():
            atomic_write_bytes(
                extensions / "pgs_locus_audit.parquet", pgs_locus_audit_source.read_bytes()
            )
        if pgx_locus_audit_source is not None and pgx_locus_audit_source.is_file():
            atomic_write_bytes(
                extensions / "pgx_locus_audit.parquet", pgx_locus_audit_source.read_bytes()
            )
        if pgx_gene_call_audit_source is not None and pgx_gene_call_audit_source.is_file():
            atomic_write_bytes(
                extensions / "pgx_gene_call_audit.parquet",
                pgx_gene_call_audit_source.read_bytes(),
            )
        if pgx_guidance_audit_source is not None and pgx_guidance_audit_source.is_file():
            atomic_write_bytes(
                extensions / "pgx_guidance_audit.parquet",
                pgx_guidance_audit_source.read_bytes(),
            )
        if family_evidence_source is not None and family_evidence_source.is_dir():
            _copy_tree(family_evidence_source, extensions / "family_evidence.parquet")
        if gwas_evidence_source.is_dir():
            _copy_tree(
                gwas_evidence_source,
                extensions / "gwas_catalog_evidence.parquet",
            )
        if clingen_gene_disease_source.is_file():
            atomic_write_bytes(
                extensions / "clingen_gene_disease.parquet",
                clingen_gene_disease_source.read_bytes(),
            )
        if clingen_dosage_source.is_file():
            atomic_write_bytes(
                extensions / "clingen_dosage_sensitivity.parquet",
                clingen_dosage_source.read_bytes(),
            )
        _copy_tree(
            variant_dir / "variant_provenance.parquet",
            extensions / "variant_provenance.parquet",
        )
        _copy_tree(
            ingest_dir / "reference_evidence.parquet",
            extensions / "reference_confidence_regions.parquet",
        )
        partial_name = "observed_partial_variants.parquet"
        partial_rows = variant_summary.get("observed_partial_variant_rows", 0)
        if partial_rows:
            _copy_tree(variant_dir / partial_name, extensions / partial_name)
        non_core_rows = _parquet_rows(ingest_dir / "non_core_variants.parquet")
        if non_core_rows:
            _copy_tree(
                ingest_dir / "non_core_variants.parquet",
                extensions / "non_core_variants.parquet",
            )
        reconciliation = {
            "status": "PASS",
            "input_sha256": preflight["input"]["sha256"],
            "accepted_source_alleles": ingest_summary["accepted_evidence_multiset"],
            "normalized_output_alleles": variant_summary["source_allele_multiset"],
            "unique_source_alleles": variant_summary["unique_source_alleles"],
            "unique_normalized_variants": variant_summary["unique_normalized_variants"],
            "non_core_variant_rows": non_core_rows,
            "observed_partial_variant_rows": partial_rows,
            "source_record_counts": ingest_summary["counters"],
            "missing_to_reference_conversions": 0,
        }
        atomic_write_json(extensions / "reconciliation.json", reconciliation)

        if config.agent_exploration_profile:
            _write_profile(staging)

        table_paths = {
            "variants": staging / "variants.parquet",
            "callable_regions": staging / "callable_regions.parquet",
            "variant_provenance": extensions / "variant_provenance.parquet",
            "reference_confidence_regions": extensions / "reference_confidence_regions.parquet",
        }
        if callability_source is not None and callability_source.is_file():
            table_paths["callability"] = staging / "callability.parquet"
        if (
            target_callability_evidence_source is not None
            and target_callability_evidence_source.is_file()
        ):
            table_paths["target_callability_evidence"] = (
                extensions / "target_callability_evidence.parquet"
            )
        if prs_source is not None and prs_source.is_file():
            table_paths["prs"] = staging / "prs.parquet"
        if pgs_model_audit_source is not None and pgs_model_audit_source.is_file():
            table_paths["pgs_model_audit"] = extensions / "pgs_model_audit.parquet"
        if pgs_locus_audit_source is not None and pgs_locus_audit_source.is_file():
            table_paths["pgs_locus_audit"] = extensions / "pgs_locus_audit.parquet"
        if pharmacogenomics_source is not None and pharmacogenomics_source.is_file():
            table_paths["pharmacogenomics"] = staging / "pharmacogenomics.parquet"
        if pgx_locus_audit_source is not None and pgx_locus_audit_source.is_file():
            table_paths["pgx_locus_audit"] = extensions / "pgx_locus_audit.parquet"
        if pgx_gene_call_audit_source is not None and pgx_gene_call_audit_source.is_file():
            table_paths["pgx_gene_call_audit"] = extensions / "pgx_gene_call_audit.parquet"
        if pgx_guidance_audit_source is not None and pgx_guidance_audit_source.is_file():
            table_paths["pgx_guidance_audit"] = extensions / "pgx_guidance_audit.parquet"
        if family_evidence_source is not None and family_evidence_source.is_dir():
            table_paths["family_evidence"] = extensions / "family_evidence.parquet"
        if gwas_associations_source.is_dir():
            table_paths["gwas_associations"] = staging / "gwas_associations.parquet"
        if gwas_evidence_source.is_dir():
            table_paths["gwas_catalog_evidence"] = extensions / "gwas_catalog_evidence.parquet"
        if gene_index_source.is_file():
            table_paths["gene_index"] = staging / "gene_index.parquet"
        if clinical_evidence_source.is_dir():
            table_paths["clinical_evidence"] = staging / "clinical_evidence.parquet"
        if clinical_findings_source.is_file():
            table_paths["clinical_findings"] = staging / "clinical_findings.parquet"
        if interpretations_source.is_file():
            table_paths["interpretations"] = staging / "interpretations.parquet"
        if clinvar_assertions_source.is_dir():
            table_paths["clinvar_assertions"] = extensions / "clinvar_assertions.parquet"
        if acmg_evidence_source.is_file():
            table_paths["acmg_secondary_findings_evidence"] = (
                extensions / "acmg_secondary_findings_evidence.parquet"
            )
        if clingen_gene_disease_source.is_file():
            table_paths["clingen_gene_disease"] = extensions / "clingen_gene_disease.parquet"
        if clingen_dosage_source.is_file():
            table_paths["clingen_dosage_sensitivity"] = (
                extensions / "clingen_dosage_sensitivity.parquet"
            )
        if partial_rows:
            table_paths["observed_partial_variants"] = extensions / partial_name
        if non_core_rows:
            table_paths["non_core_variants"] = extensions / "non_core_variants.parquet"
        glossary_fields: dict[str, dict[str, Any]] = {}
        for table, path in table_paths.items():
            glossary_fields.update(_schema_fields(table, _first_parquet_schema(path), config))
        atomic_write_json(
            staging / "field_glossary.json",
            {"schema_version": SCHEMA_VERSION, "fields": glossary_fields},
        )

        configuration_digest = _configuration_digest(config)
        raw_artifacts = {
            "schema_version": SCHEMA_VERSION,
            "artifacts": [
                {
                    "artifact_type": "reference-confidence VCF/gVCF",
                    "disposition": "separate",
                    "reason": "immutable sovereignty record is retained outside the derived bundle",
                    "integrity_identifier": f"sha256:{preflight['input']['sha256']}",
                    "limitations": [
                        "the source VCF is required to recover omitted source fields or re-run conversion"
                    ],
                },
                {
                    "artifact_type": "GRCh38 reference FASTA",
                    "disposition": "referenced",
                    "reason": "large public reference is identified by digest rather than embedded",
                    "integrity_identifier": f"sha256:{preflight['reference']['sha256']}",
                    "limitations": ["independent REF revalidation requires the exact FASTA"],
                },
                {
                    "artifact_type": "aligned sequencing reads",
                    "disposition": "never_received",
                    "reason": "the generic small-variant core build accepts a VCF/gVCF input",
                    "integrity_identifier": None,
                    "limitations": [
                        "read-backed SV/CNV, repeat, mosaicism, and local alignment review are unavailable"
                    ],
                },
            ],
        }
        if family_evidence_summary is not None:
            raw_artifacts["artifacts"].extend(
                [
                    {
                        "artifact_type": "multisample family joint VCF",
                        "disposition": "referenced",
                        "reason": (
                            "relative genotypes remain outside the portable subject bundle; "
                            "only privacy-minimized derived statuses are embedded"
                        ),
                        "integrity_identifier": (
                            "sha256:" + family_evidence_summary["joint_vcf_sha256"]
                        ),
                        "limitations": [
                            "re-derivation requires authorized access to the byte-identical family VCF"
                        ],
                    },
                    {
                        "artifact_type": "pedigree relationship file",
                        "disposition": "referenced",
                        "reason": "relative identifiers are not embedded in the portable bundle",
                        "integrity_identifier": (
                            "sha256:" + family_evidence_summary["pedigree_sha256"]
                        ),
                        "limitations": [
                            "relationship revalidation requires the byte-identical PED"
                        ],
                    },
                ]
            )
            if family_evidence_summary.get("loco") is not None:
                raw_artifacts["artifacts"].append(
                    {
                        "artifact_type": "child-independent leave-one-out audit",
                        "disposition": "referenced",
                        "reason": (
                            "the large exception ledger and private audit summary remain separate"
                        ),
                        "integrity_identifier": (
                            "sha256:"
                            + family_evidence_summary["loco"]["review_sha256"]
                            + ";summary-sha256:"
                            + family_evidence_summary["loco"]["summary_sha256"]
                        ),
                        "limitations": [
                            "absence from the exception ledger is not exposed as per-site consistency"
                        ],
                    }
                )
        if target_callability_summary is not None:
            raw_artifacts["artifacts"].append(
                {
                    "artifact_type": "frozen target catalogue",
                    "disposition": "referenced",
                    "reason": (
                        "the catalogue is identified by digest; normalized target identity and "
                        "row-level provenance are retained in the callability evidence extension"
                    ),
                    "integrity_identifier": (
                        "sha256:" + target_callability_summary["catalog"]["sha256"]
                    ),
                    "limitations": [
                        "re-materialization requires the byte-identical target catalogue"
                    ],
                }
            )
        if config.gwas.enabled:
            gwas_resource = preflight["resources"]["gwas_catalog"]
            raw_artifacts["artifacts"].extend(
                [
                    {
                        "artifact_type": "dated NHGRI-EBI GWAS Catalog source archive",
                        "disposition": "referenced",
                        "reason": (
                            "the full public source archive is never embedded; exact source "
                            "identity, terms, and normalization boundary are retained"
                        ),
                        "integrity_identifier": "sha256:" + gwas_resource["archive_sha256"],
                        "source_url": gwas_resource["source_url"],
                        "license": gwas_resource["license"],
                        "redistribution_allowed": False,
                        "limitations": [
                            "re-normalization requires the byte-identical dated source archive",
                            "downstream redistribution requires an independent review of the "
                            "recorded EMBL-EBI terms",
                        ],
                    },
                    {
                        "artifact_type": "normalized GWAS Catalog association resource",
                        "disposition": "referenced",
                        "reason": (
                            "the complete normalized resource remains outside the sample bundle; "
                            "only exact matches and their lossless audit rows are embedded"
                        ),
                        "integrity_identifier": (
                            "sha256:" + gwas_resource["normalized_parquet_sha256"]
                        ),
                        "license": gwas_resource["license"],
                        "redistribution_allowed": config.gwas.redistribution_allowed,
                        "limitations": [
                            "re-projection requires the byte-identical normalized resource",
                            "an absent exact match is not a negative trait or clinical result",
                        ],
                    },
                ]
            )
        if config.clingen.enabled:
            clingen_resource = preflight["resources"]["clingen"]
            source_digests = ",".join(
                f"{name}=sha256:{metadata['sha256']}"
                for name, metadata in sorted(clingen_resource["sources"].items())
            )
            raw_artifacts["artifacts"].append(
                {
                    "artifact_type": "ClinGen gene-level curation snapshot",
                    "disposition": "referenced",
                    "reason": (
                        "the complete CC0 snapshot remains in the versioned builder resource; "
                        "only exact gene-symbol matches are embedded in this subject bundle"
                    ),
                    "integrity_identifier": (
                        "canonical-sha256:"
                        + clingen_resource["canonical_sha256"]
                        + ";manifest-sha256:"
                        + clingen_resource["manifest_sha256"]
                        + ";sources:"
                        + source_digests
                    ),
                    "license": clingen_resource["license"],
                    "redistribution_allowed": True,
                    "attribution": clingen_resource["attribution"],
                    "coverage": {
                        "evidence_level": "gene",
                        "domains": ["gene-disease validity", "dosage sensitivity"],
                        "assembly": "GRCh38 for dosage coordinates",
                    },
                    "limitations": [
                        clingen_resource["semantic_boundary"],
                        "an absent gene-symbol match is unavailable context, not a negative result",
                    ],
                }
            )
        if config.predictors.enabled:
            predictor_resources = preflight["resources"]["predictors"]
            for plugin in config.predictors.enabled_plugins:
                resource = predictor_resources["plugins"][plugin]
                artifact_digests = resource["artifact_sha256"]
                limitations = [
                    "the complete predictor resource remains outside the portable bundle",
                    "the score is computational evidence, not a diagnosis or clinical classification",
                    "re-annotation requires every byte-identical artifact named by digest",
                ]
                if plugin == "spliceai":
                    limitations.append(
                        "the pinned Ensembl MANE v1.4 resource covers SNVs only; indel "
                        "scores are explicitly unavailable"
                    )
                raw_artifacts["artifacts"].append(
                    {
                        "artifact_type": f"offline VEP predictor resource: {plugin}",
                        "disposition": "referenced",
                        "reason": (
                            "large and independently licensed source files are identified by "
                            "a private checksum manifest rather than embedded"
                        ),
                        "integrity_identifier": (
                            "manifest-sha256:"
                            + predictor_resources["manifest_sha256"]
                            + ";artifacts:"
                            + ",".join(
                                f"{key}=sha256:{digest}"
                                for key, digest in sorted(artifact_digests.items())
                            )
                        ),
                        "license": resource["license"],
                        "redistribution_allowed": False,
                        "coverage": resource["coverage"],
                        "limitations": limitations,
                    }
                )
            semantic_controls = predictor_resources["semantic_controls"]
            raw_artifacts["artifacts"].append(
                {
                    "artifact_type": "offline VEP predictor semantic control report",
                    "disposition": "referenced",
                    "reason": (
                        "the fixed public positive-control run remains outside the subject "
                        "bundle; its deterministic results and report digest are retained"
                    ),
                    "integrity_identifier": "sha256:" + semantic_controls["report_sha256"],
                    "control_set": semantic_controls["control_set"],
                    "control_definition_sha256": semantic_controls["control_definition_sha256"],
                    "doctor_run_key": semantic_controls["doctor_run_key"],
                    "enabled_plugins": semantic_controls["enabled_plugins"],
                    "record_identity_sha256": semantic_controls["record_identity_sha256"],
                    "reference_fasta_sha256": semantic_controls["reference_fasta_sha256"],
                    "vep_cache_manifest_sha256": semantic_controls["vep_cache_manifest_sha256"],
                    "limitations": [semantic_controls["semantic_boundary"]],
                }
            )
        if pgs_score_summary is not None:
            for model in config.pgs_models:
                raw_artifacts["artifacts"].append(
                    {
                        "artifact_type": f"PGS Catalog scoring file {model.pgs_id}",
                        "disposition": "referenced",
                        "reason": (
                            "the licensed scoring resource remains separate; its exact digest and "
                            "derived locus decisions are retained"
                        ),
                        "integrity_identifier": f"sha256:{model.expected_sha256}",
                        "limitations": [
                            "re-scoring requires the byte-identical scoring file and GRCh38 FASTA"
                        ],
                    }
                )
                if model.calibration is not None:
                    raw_artifacts["artifacts"].append(
                        {
                            "artifact_type": f"PGS calibration table {model.pgs_id}",
                            "disposition": "referenced",
                            "reason": (
                                "population calibration remains separate and is bound by digest"
                            ),
                            "integrity_identifier": (f"sha256:{model.calibration.expected_sha256}"),
                            "limitations": [
                                "the percentile is applicable only to the declared reference population"
                            ],
                        }
                    )
        if pgx_projection_summary is not None:
            pgx_resource = preflight["resources"]["pgx"]
            raw_artifacts["artifacts"].extend(
                [
                    {
                        "artifact_type": "official PharmCAT positions VCF",
                        "disposition": "referenced",
                        "reason": (
                            "the immutable matcher boundary remains outside the bundle; every "
                            "derived locus retains its source and materialization digest"
                        ),
                        "integrity_identifier": (
                            "sha256:" + pgx_resource["pharmcat"]["positions_vcf_sha256"]
                        ),
                        "limitations": [
                            "re-materialization requires the byte-identical PharmCAT positions asset"
                        ],
                    },
                    {
                        "artifact_type": "frozen CPIC pair-view snapshot",
                        "disposition": "referenced",
                        "reason": (
                            "the licensed source snapshot remains separate and its exact digest, "
                            "version, date, and normative row boundary are retained"
                        ),
                        "integrity_identifier": "sha256:" + pgx_resource["cpic"]["sha256"],
                        "limitations": [
                            "re-projection requires the byte-identical frozen CPIC snapshot"
                        ],
                    },
                    {
                        "artifact_type": "materialized PharmCAT input VCF",
                        "disposition": "separate",
                        "reason": (
                            "sample-derived matcher input remains in the restart-safe work area; "
                            "the complete per-locus decision ledger is embedded"
                        ),
                        "integrity_identifier": (
                            "sha256:" + pgx_materialization_summary["outputs"]["pgx-input.vcf.gz"]
                        ),
                        "limitations": [
                            "re-running PharmCAT requires the retained input VCF and index"
                        ],
                    },
                    {
                        "artifact_type": "raw and canonical PharmCAT reports",
                        "disposition": "separate",
                        "reason": (
                            "full sample-derived reports remain outside the portable bundle; "
                            "typed gene and guidance audits retain their content digests"
                        ),
                        "integrity_identifier": (
                            "sha256:"
                            + pharmcat_summary["canonical_report_sha256"]
                            + ";raw-sha256:"
                            + pharmcat_summary["raw_report_sha256"]
                        ),
                        "limitations": [
                            "verbatim report recovery requires the retained restart-safe work area"
                        ],
                    },
                ]
            )
        atomic_write_json(staging / "raw_artifacts.json", raw_artifacts)

        transformation_stages = [
            {
                "stage_id": "preflight",
                "software": "genome-annotation-builder",
                "software_version": PIPELINE_VERSION,
                "configuration_digest": configuration_digest,
                "inputs": [f"sha256:{preflight['input']['sha256']}"],
                "outputs": ["preflight:PASS"],
                "status": "completed",
            },
            {
                "stage_id": "lossless_vcf_decomposition",
                "software": "genome-annotation-builder/pysam",
                "software_version": PIPELINE_VERSION,
                "configuration_digest": configuration_digest,
                "inputs": [f"sha256:{preflight['input']['sha256']}"],
                "outputs": [
                    "projected carried alleles",
                    "reference-confidence evidence",
                    "non-core evidence",
                ],
                "status": "completed",
            },
            {
                "stage_id": "grch38_normalization",
                "software": "bcftools",
                "software_version": preflight["software"]["bcftools"],
                "configuration_digest": configuration_digest,
                "inputs": ["projected carried alleles", preflight["reference"]["sha256"]],
                "outputs": ["normalized carried alleles"],
                "started_at": normalization_summary.get("started_at"),
                "completed_at": normalization_summary.get("completed_at"),
                "status": "completed",
            },
        ]
        if family_evidence_summary is not None:
            transformation_stages.append(
                {
                    "stage_id": "privacy_minimized_family_evidence",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "normalized carried alleles",
                        "sha256:" + family_evidence_summary["joint_vcf_sha256"],
                        "sha256:" + family_evidence_summary["pedigree_sha256"],
                    ],
                    "outputs": ["extensions/family_evidence.parquet"],
                    "subject_genotype_refinements": family_evidence_summary["genotype_refinements"],
                    "relative_identifiers_persisted": family_evidence_summary[
                        "relative_identifiers_persisted"
                    ],
                    "relative_genotypes_persisted": family_evidence_summary[
                        "relative_genotypes_persisted"
                    ],
                    "normalization_projection": family_evidence_summary["normalization_projection"],
                    "status": "completed",
                }
            )
        if vep_summary is not None:
            predictor_inputs = []
            if config.predictors.enabled:
                predictor_inputs.append("sha256:" + vep_summary["predictors"]["manifest_sha256"])
            transformation_stages.append(
                {
                    "stage_id": "ensembl_vep_annotation",
                    "software": "Ensembl VEP",
                    "software_version": vep_summary["vep_version"],
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "normalized carried alleles",
                        f"sha256:{vep_summary['cache_tree_sha256']}",
                        *predictor_inputs,
                    ],
                    "outputs": [
                        "typed transcript and population annotations",
                        *(
                            ["typed offline predictor annotations"]
                            if config.predictors.enabled
                            else []
                        ),
                    ],
                    "enabled_predictors": list(config.predictors.enabled_plugins),
                    "predictor_coverage": {
                        plugin: vep_summary["predictors"]["plugins"][plugin]["coverage"]
                        for plugin in config.predictors.enabled_plugins
                    },
                    "started_at": vep_summary.get("started_at"),
                    "completed_at": vep_summary.get("completed_at"),
                    "container_image": vep_summary.get("container_image"),
                    "command_semantics": vep_summary.get("command_semantics"),
                    "status": "completed",
                }
            )
        if config.predictors.enabled:
            transformation_stages.append(
                _predictor_semantic_control_stage(
                    preflight["resources"]["predictors"], configuration_digest
                )
            )
        if config.clinvar_enabled:
            transformation_stages.append(
                {
                    "stage_id": "clinvar_exact_allele_join",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "normalized carried alleles",
                        "sha256:" + preflight["resources"]["clinvar"]["normalized_vcf_sha256"],
                    ],
                    "outputs": [
                        "normative clinical_evidence retrieval rows",
                        "lossless ClinVar assertion extension",
                        "typed variant pathogenicity fields",
                    ],
                    "matched_variants": variant_summary["counters"].get(
                        "clinvar_matched_variants", 0
                    ),
                    "status": "completed",
                }
            )
        if config.acmg_sf_enabled:
            acmg_summary = variant_summary["annotation_resources"]["acmg_sf"]
            transformation_stages.append(
                {
                    "stage_id": "acmg_secondary_findings_policy_evaluation",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "typed transcript annotations",
                        "exact-allele ClinVar evidence",
                        "sha256:" + preflight["resources"]["acmg_sf"]["canonical_sha256"],
                    ],
                    "outputs": [
                        "conservative ACMG reportable-candidate decisions",
                        "variant-level is_actionable policy flags",
                        "clinical_findings and algorithmic interpretations",
                    ],
                    "decision_rows": acmg_summary["decision_rows"],
                    "actionable_candidate_variants": acmg_summary["actionable_candidate_variants"],
                    "status": "completed",
                }
            )
        if config.clingen.enabled:
            clingen_summary = variant_summary["annotation_resources"]["clingen"]
            transformation_stages.append(
                {
                    "stage_id": "clingen_gene_context_projection",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "exact HGNC symbols from typed VEP transcript annotations",
                        "sha256:" + clingen_summary["canonical_sha256"],
                        "sha256:" + clingen_summary["manifest_sha256"],
                    ],
                    "outputs": sorted("extensions/" + name for name in clingen_summary["outputs"]),
                    "matching_policy": "exact case-sensitive HGNC gene symbol",
                    "input_gene_symbols": clingen_summary["input_gene_symbols"],
                    "matched_gene_symbols": clingen_summary["matched_gene_symbols"],
                    "unmatched_gene_symbols": clingen_summary["unmatched_gene_symbols"],
                    "gene_disease_rows": clingen_summary["gene_disease_rows"],
                    "dosage_rows": clingen_summary["dosage_rows"],
                    "semantic_boundary": clingen_summary["semantic_boundary"],
                    "variant_pathogenicity_claims_created": 0,
                    "subject_diagnosis_claims_created": 0,
                    "status": "completed",
                }
            )
        if config.gwas.enabled:
            gwas_resource = preflight["resources"]["gwas_catalog"]
            transformation_stages.append(
                {
                    "stage_id": "gwas_catalog_exact_association_projection",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "normalized carried alleles",
                        "sha256:" + gwas_resource["normalized_parquet_sha256"],
                    ],
                    "outputs": [
                        "gwas_associations.parquet",
                        "extensions/gwas_catalog_evidence.parquet",
                        "typed variant trait_associations summary",
                    ],
                    "matching_policy": gwas_resource["matching_policy"],
                    "effect_type_policy": gwas_resource["effect_type_policy"],
                    "matched_variants": variant_summary["counters"].get("gwas_matched_variants", 0),
                    "exact_association_rows": variant_summary["counters"].get(
                        "gwas_exact_association_rows", 0
                    ),
                    "normative_association_rows": variant_summary["counters"].get(
                        "gwas_normative_association_rows", 0
                    ),
                    "untyped_exact_rows": variant_summary["counters"].get(
                        "gwas_untyped_exact_rows", 0
                    ),
                    "clinical_grade_claims": 0,
                    "status": "completed",
                }
            )
        if target_callability_summary is not None:
            transformation_stages.append(
                {
                    "stage_id": "target_level_callability_materialization",
                    "software": "genome-annotation-builder/DuckDB",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "reference-confidence evidence",
                        "normalized carried alleles",
                        "sha256:" + target_callability_summary["catalog"]["sha256"],
                    ],
                    "outputs": [
                        "callability.parquet",
                        "extensions/target_callability_evidence.parquet",
                    ],
                    "target_rows": target_callability_summary["catalog"]["target_rows"],
                    "resolution_counts": target_callability_summary["resolution_counts"],
                    "precedence_rule": target_callability_summary["precedence_rule"],
                    "status": "completed",
                }
            )
        if pgs_score_summary is not None:
            transformation_stages.append(
                {
                    "stage_id": "polygenic_score_materialization",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "checksum-pinned PGS Catalog scoring rows",
                        "extensions/target_callability_evidence.parquet",
                        "normalized carried alleles",
                    ],
                    "outputs": [
                        "extensions/pgs_locus_audit.parquet",
                        "extensions/pgs_model_audit.parquet",
                        *(["prs.parquet"] if pgs_score_summary["published_model_count"] else []),
                    ],
                    "model_count": pgs_score_summary["model_count"],
                    "published_model_count": pgs_score_summary["published_model_count"],
                    "publication_status_counts": pgs_score_summary["publication_status_counts"],
                    "missing_genotypes_imputed": 0,
                    "interval_only_reference_conversions": 0,
                    "status": "completed",
                }
            )
        if pgx_materialization_summary is not None:
            transformation_stages.append(
                {
                    "stage_id": "pgx_explicit_genotype_materialization",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "official digest-pinned PharmCAT positions",
                        "extensions/target_callability_evidence.parquet",
                        "normalized carried alleles",
                    ],
                    "outputs": [
                        "separate materialized PharmCAT input VCF",
                        "extensions/pgx_locus_audit.parquet",
                    ],
                    "record_count": pgx_materialization_summary["record_count"],
                    "allele_count": pgx_materialization_summary["allele_count"],
                    "status_counts": pgx_materialization_summary["status_counts"],
                    "missing_to_reference_conversions": 0,
                    "interval_only_reference_conversions": 0,
                    "structurally_complex_matcher_enabled": False,
                    "status": "completed",
                }
            )
        if pharmcat_summary is not None:
            transformation_stages.append(
                {
                    "stage_id": "offline_pharmcat_matching",
                    "software": "PharmCAT",
                    "software_version": pharmcat_summary["version"],
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "separate materialized PharmCAT input VCF",
                        "official PharmCAT allele definitions",
                    ],
                    "outputs": ["separate raw and canonical PharmCAT reports"],
                    "container_image": pharmcat_summary["container_image"],
                    "command_semantics": pharmcat_summary["command_semantics"],
                    "network": pharmcat_summary["network"],
                    "read_only_root": pharmcat_summary["read_only_root"],
                    "dangerous_missing_to_reference_options_used": False,
                    "research_mode_used": False,
                    "excluded_complex_genes": pharmcat_summary["excluded_complex_genes"],
                    "status": "completed",
                }
            )
        if pgx_projection_summary is not None:
            transformation_stages.append(
                {
                    "stage_id": "cpic_gated_pgx_projection",
                    "software": "genome-annotation-builder",
                    "software_version": PIPELINE_VERSION,
                    "configuration_digest": configuration_digest,
                    "inputs": [
                        "canonical PharmCAT report",
                        "PharmCAT calls-only report",
                        "frozen CPIC pair-view snapshot",
                    ],
                    "outputs": [
                        "extensions/pgx_gene_call_audit.parquet",
                        "optional pharmacogenomics.parquet",
                        "optional extensions/pgx_guidance_audit.parquet",
                    ],
                    "gene_rows": pgx_projection_summary["gene_rows"],
                    "published_gene_rows": pgx_projection_summary["published_gene_rows"],
                    "call_status_counts": pgx_projection_summary["call_status_counts"],
                    "complex_gene_core_rows": 0,
                    "status": "completed",
                }
            )
        transformation_stages.append(
            {
                "stage_id": "genome_bundle_assembly",
                "software": "genome-annotation-builder",
                "software_version": PIPELINE_VERSION,
                "configuration_digest": configuration_digest,
                "inputs": [
                    (
                        "VEP-annotated normalized carried alleles"
                        if vep_summary is not None
                        else "normalized carried alleles"
                    ),
                    "reference-confidence evidence",
                ],
                "outputs": [
                    "variants.parquet",
                    "optional normative clinical tables",
                    "extensions",
                    "metadata",
                ],
                "completed_at": generated_at,
                "status": "completed",
            }
        )
        transformations = {
            "schema_version": SCHEMA_VERSION,
            "stages": transformation_stages,
        }
        atomic_write_json(staging / "transformations.json", transformations)

        predictor_labels = {
            "loftee": "LOFTEE",
            "cadd": "CADD",
            "revel": "REVEL",
            "alphamissense": "AlphaMissense",
            "spliceai": "SpliceAI",
        }
        missing_predictors = [
            predictor_labels[name]
            for name in predictor_labels
            if name not in config.predictors.enabled_plugins
        ]
        if config.vep_enabled and config.clinvar_enabled:
            missing_scope = (
                "; disabled predictor adapters: " + ", ".join(missing_predictors)
                if missing_predictors
                else ""
            )
            annotation_omission = {
                "scope": "external annotation/evidence sources outside the pinned adapter set"
                + missing_scope,
                "reason": (
                    "only explicitly enabled, checksum-pinned resources are projected into the bundle"
                ),
                "pipeline_stage": "annotation",
                "interpretation_consequence": (
                    "null/false fields do not assert benignity or lack of actionability"
                ),
                "recovery": "add a licensed, checksum-pinned adapter and rebuild",
            }
        elif config.vep_enabled:
            missing_scope = (
                "; disabled predictors: " + ", ".join(missing_predictors)
                if missing_predictors
                else ""
            )
            annotation_omission = {
                "scope": "ClinVar classifications and external clinical evidence" + missing_scope,
                "reason": (
                    "functional VEP annotation was enabled, but ClinVar and any listed predictor "
                    "adapters were not part of this build"
                ),
                "pipeline_stage": "annotation",
                "interpretation_consequence": (
                    "null/false fields do not assert benignity or lack of actionability"
                ),
                "recovery": "enable the corresponding frozen resources and rebuild",
            }
        elif config.clinvar_enabled:
            annotation_omission = {
                "scope": "functional, population, and prediction annotations",
                "reason": "ClinVar was enabled, but no functional/population adapter was enabled",
                "pipeline_stage": "annotation",
                "interpretation_consequence": (
                    "null fields do not assert benignity, commonness, or lack of effect"
                ),
                "recovery": "re-run with frozen VEP/gnomAD/predictor resources",
            }
        else:
            annotation_omission = {
                "scope": "functional, population, clinical, and prediction annotations",
                "reason": "no checksum-pinned annotation adapter was enabled for this build",
                "pipeline_stage": "annotation",
                "interpretation_consequence": (
                    "null/false fields do not assert benignity, commonness, or lack of actionability"
                ),
                "recovery": "re-run with frozen VEP/ClinVar/gnomAD/predictor resources",
            }
        omission_rows = [annotation_omission]
        for plugin_name, label in predictor_labels.items():
            if plugin_name in config.predictors.enabled_plugins:
                continue
            omission_rows.append(
                {
                    "scope": f"{label} predictor annotation",
                    "reason": f"the checksum-pinned {label} adapter was not enabled",
                    "pipeline_stage": "functional_prediction",
                    "interpretation_consequence": (
                        f"a null {label} field is unavailable evidence, not a benign prediction"
                    ),
                    "recovery": (
                        f"review the {label} resource terms, enable the pinned adapter, and rebuild"
                    ),
                }
            )
        if config.predictors.spliceai.enabled:
            omission_rows.append(
                {
                    "scope": "SpliceAI indel predictions",
                    "reason": (
                        "the pinned Ensembl MANE v1.4 release contains SNV scores only; an "
                        "explicit zero-record sentinel satisfies the VEP plugin interface"
                    ),
                    "pipeline_stage": "functional_prediction",
                    "interpretation_consequence": (
                        "null SpliceAI scores for indels are unavailable, never zero-effect calls"
                    ),
                    "recovery": (
                        "supply a separately licensed, official, checksum-pinned GRCh38 indel "
                        "resource through a future validated adapter"
                    ),
                }
            )
        if config.acmg_sf_enabled:
            omission_rows.append(
                {
                    "scope": "clinically reviewed ACMG secondary findings",
                    "reason": (
                        "the automated ACMG policy layer emits reportable candidates, not "
                        "professional case-level interpretations"
                    ),
                    "pipeline_stage": "clinical_interpretation",
                    "interpretation_consequence": (
                        "is_actionable marks only a reproducible policy candidate; clinical_grade "
                        "remains false and confirmation is required"
                    ),
                    "recovery": (
                        "perform variant-level clinical review, orthogonal confirmation where "
                        "appropriate, and publish a separately provenance-tracked interpretation"
                    ),
                }
            )
        else:
            omission_rows.append(
                {
                    "scope": "ACMG secondary-findings policy evaluation",
                    "reason": "the checksum-pinned ACMG SF adapter was not enabled",
                    "pipeline_stage": "clinical_interpretation",
                    "interpretation_consequence": (
                        "is_actionable remains false and cannot be read as a negative result"
                    ),
                    "recovery": "enable VEP, ClinVar, and the frozen ACMG SF policy adapter",
                }
            )
        if not gene_index_source.is_file():
            omission_rows.append(
                {
                    "scope": "gene_index.parquet",
                    "reason": "a derived gene index has not yet been materialized",
                    "pipeline_stage": "bundle_assembly",
                    "interpretation_consequence": (
                        "gene-scoped retrieval cannot resolve unannotated rows"
                    ),
                    "recovery": "re-annotate variants and rebuild the derived index",
                }
            )
        if target_callability_summary is None:
            omission_rows.append(
                {
                    "scope": "site-targeted callability.parquet",
                    "reason": "no explicit target catalogue was configured",
                    "pipeline_stage": "callability",
                    "interpretation_consequence": (
                        "exact queries use the retained reference-confidence ledger; "
                        "empty lookup is not negative"
                    ),
                    "recovery": (
                        "materialize target sites from a frozen catalogue against the retained ledger"
                    ),
                }
            )
        if family_evidence_summary is None:
            omission_rows.append(
                {
                    "scope": "pedigree-aware family evidence",
                    "reason": "no digest-pinned family joint VCF and PED were configured",
                    "pipeline_stage": "family_evidence",
                    "interpretation_consequence": (
                        "inheritance, parental origin, Mendelian, phasing, LOCO, and replicate "
                        "statuses are unavailable; absence is not family confirmation"
                    ),
                    "recovery": (
                        "configure authorized immutable family inputs and rebuild without "
                        "embedding relative identifiers or genotypes"
                    ),
                }
            )
        if not config.clingen.enabled:
            omission_rows.append(
                {
                    "scope": "ClinGen gene-disease validity and dosage-sensitivity context",
                    "reason": "the checksum-pinned ClinGen gene-level adapter was not enabled",
                    "pipeline_stage": "clinical_gene_context",
                    "interpretation_consequence": (
                        "absence of ClinGen context is unavailable evidence, not evidence that a "
                        "gene lacks a disease relationship or dosage mechanism"
                    ),
                    "recovery": "enable and verify the vendored CC0 ClinGen snapshot, then rebuild",
                }
            )
        else:
            for source, scope in (
                (clingen_gene_disease_source, "ClinGen gene-disease projection table"),
                (clingen_dosage_source, "ClinGen dosage-sensitivity projection table"),
            ):
                if source.is_file():
                    continue
                omission_rows.append(
                    {
                        "scope": scope,
                        "reason": (
                            "the adapter completed, but no exact HGNC gene-symbol match existed "
                            "for this bundle"
                        ),
                        "pipeline_stage": "clinical_gene_context",
                        "interpretation_consequence": (
                            "the absent optional table is not a negative clinical or genetic result"
                        ),
                        "recovery": (
                            "inspect VEP gene-symbol coverage and the immutable ClinGen snapshot; "
                            "do not infer a match through aliases without a versioned mapping"
                        ),
                    }
                )
        if not config.gwas.enabled:
            omission_rows.append(
                {
                    "scope": "GWAS Catalog association evidence",
                    "reason": "no checksum-pinned dated GWAS Catalog release was configured",
                    "pipeline_stage": "population_association_annotation",
                    "interpretation_consequence": (
                        "absence of GWAS rows is not evidence for or against any trait"
                    ),
                    "recovery": (
                        "freeze and verify a dated NHGRI-EBI GWAS Catalog release, then rebuild"
                    ),
                }
            )
        else:
            gwas_resource = preflight["resources"]["gwas_catalog"]
            gwas_counts = gwas_resource["counts"]
            excluded_counts = {
                key: int(gwas_counts.get(key, 0))
                for key in (
                    "omitted_non_single_core_chromosome",
                    "omitted_unknown_or_complex_risk_allele",
                    "omitted_missing_or_invalid_p_value",
                    "omitted_missing_reported_trait",
                )
                if int(gwas_counts.get(key, 0))
            }
            if excluded_counts:
                omission_rows.append(
                    {
                        "scope": "GWAS Catalog rows outside the exact projection boundary",
                        "reason": "; ".join(
                            f"{key}={value}" for key, value in sorted(excluded_counts.items())
                        ),
                        "pipeline_stage": "gwas_resource_normalization",
                        "interpretation_consequence": (
                            "multi-locus, non-core, invalid, or ambiguous-effect-allele source "
                            "rows were not guessed into a sample association"
                        ),
                        "recovery": (
                            "inspect the immutable source release and add a separately validated "
                            "representation for a currently unsupported row class"
                        ),
                    }
                )
            untyped_rows = int(variant_summary["counters"].get("gwas_untyped_exact_rows", 0))
            if untyped_rows:
                omission_rows.append(
                    {
                        "scope": "normative effect estimates for exact GWAS matches",
                        "reason": f"ambiguous_or_missing_effect_type_or_size={untyped_rows}",
                        "pipeline_stage": "population_association_annotation",
                        "interpretation_consequence": (
                            "the exact source rows remain in the lossless audit extension but "
                            "are withheld from gwas_associations.parquet"
                        ),
                        "recovery": (
                            "use a future source field or independently validated mapping that "
                            "states effect type without inferring it from the legacy OR-or-BETA column"
                        ),
                    }
                )
            omission_rows.append(
                {
                    "scope": "personal causal or clinical interpretation of GWAS associations",
                    "reason": (
                        "the source contains literature-curated population associations, not "
                        "case-level diagnostic evidence"
                    ),
                    "pipeline_stage": "population_association_annotation",
                    "interpretation_consequence": (
                        "a matched effect allele is not a diagnosis, causal finding, penetrance "
                        "estimate, or clinical recommendation"
                    ),
                    "recovery": (
                        "perform a separately governed evidence review with ancestry, study design, "
                        "replication, effect scale, and clinical context"
                    ),
                }
            )
        if pgs_score_summary is None:
            omission_rows.append(
                {
                    "scope": "polygenic scores and prs.parquet",
                    "reason": "no checksum-pinned PGS Catalog model was configured",
                    "pipeline_stage": "polygenic_scoring",
                    "interpretation_consequence": (
                        "absence of a score is not a low-risk or average-risk result"
                    ),
                    "recovery": (
                        "configure a GRCh38-harmonized scoring file, exact digest, licence, "
                        "missingness policy, and population calibration"
                    ),
                }
            )
        else:
            unpublished = [
                row
                for row in pgs_score_summary["models"]
                if row["publication_status"] != "published"
            ]
            if unpublished:
                omission_rows.append(
                    {
                        "scope": "unpublished configured polygenic scores",
                        "reason": "; ".join(
                            f"{row['pgs_id']}={row['publication_status']}" for row in unpublished
                        ),
                        "pipeline_stage": "polygenic_scoring",
                        "interpretation_consequence": (
                            "raw partial scores remain in the audit extension and are not exposed "
                            "as normative percentile-bearing PRS rows"
                        ),
                        "recovery": (
                            "resolve explicit missing/ambiguous loci or supply an applicable "
                            "digest-bound population calibration"
                        ),
                    }
                )
        if pgx_projection_summary is None:
            omission_rows.append(
                {
                    "scope": "pharmacogenomics and pharmacogenomics.parquet",
                    "reason": "no digest-pinned PharmCAT/CPIC adapter was configured",
                    "pipeline_stage": "pharmacogenomics",
                    "interpretation_consequence": (
                        "absence of a PGx row is not a normal-function or standard-dose result"
                    ),
                    "recovery": (
                        "configure the exact PharmCAT image and positions plus a frozen CPIC "
                        "pair snapshot, then rebuild from explicit genotype/reference evidence"
                    ),
                }
            )
        else:
            withheld_statuses = {
                status: count
                for status, count in pgx_projection_summary["call_status_counts"].items()
                if status != "published" and count
            }
            if withheld_statuses:
                omission_rows.append(
                    {
                        "scope": "withheld PharmCAT gene calls",
                        "reason": "; ".join(
                            f"{status}={count}"
                            for status, count in sorted(withheld_statuses.items())
                        ),
                        "pipeline_stage": "pharmacogenomics",
                        "interpretation_consequence": (
                            "only complete unambiguous non-complex calls with a normative exact "
                            "CPIC A-D mapping appear in pharmacogenomics.parquet; audit rows are "
                            "not negative findings"
                        ),
                        "recovery": (
                            "inspect the PGx gene/locus audits and use validated read-backed or "
                            "outside-call methods where the recorded reason requires them"
                        ),
                    }
                )
        omission_rows.extend(
            [
                {
                    "scope": "alternate/unplaced/decoy and symbolic calls in core variants.parquet",
                    "reason": (
                        "the .genome/1.1 core schema permits only concrete alleles on primary contigs"
                    ),
                    "pipeline_stage": "vcf_decomposition",
                    "interpretation_consequence": (
                        "these calls are not queryable as core small variants"
                    ),
                    "recovery": (
                        "inspect extensions/non_core_variants.parquet and the retained source VCF"
                    ),
                },
                {
                    "scope": (
                        "read-backed SV/CNV, repeat, complex-PGx, mosaicism, and mitochondrial analysis"
                    ),
                    "reason": "aligned reads were not an input to the small-variant build",
                    "pipeline_stage": "read_backed_extensions",
                    "interpretation_consequence": "the bundle makes no claims in those domains",
                    "recovery": "run the read-backed expansion with retained BAM/CRAM inputs",
                },
            ]
        )
        omissions = {
            "schema_version": SCHEMA_VERSION,
            "omissions": omission_rows,
        }
        atomic_write_json(staging / "omissions.json", omissions)

        relationships = {
            "schema_version": SCHEMA_VERSION,
            "relationships": [
                {
                    "relationship_id": "variant_to_source_provenance",
                    "source": ["variants.variant_id"],
                    "target": ["variant_provenance.variant_id"],
                    "cardinality": "one_to_one",
                    "completeness": "complete",
                    "unmatched_meaning": "an unmatched row is a reconciliation failure",
                },
                *(
                    [
                        {
                            "relationship_id": "variant_to_family_evidence",
                            "source": ["variants.variant_id"],
                            "target": ["family_evidence.variant_id"],
                            "cardinality": "one_to_one",
                            "completeness": "complete",
                            "unmatched_meaning": (
                                "an unmatched row means the family source no longer reconciles "
                                "to the accepted subject genotype"
                            ),
                        }
                    ]
                    if family_evidence_summary is not None
                    else []
                ),
                {
                    "relationship_id": "callable_region_to_reference_evidence",
                    "source": [
                        "callable_regions.chrom",
                        "callable_regions.start_pos",
                        "callable_regions.end_pos",
                    ],
                    "target": [
                        "reference_confidence_regions.chrom",
                        "reference_confidence_regions.start_pos",
                        "reference_confidence_regions.end_pos",
                    ],
                    "cardinality": "many_to_many",
                    "completeness": "partial",
                    "unmatched_meaning": (
                        "callability may derive from a called variant record rather than reference evidence"
                    ),
                },
                *(
                    [
                        {
                            "relationship_id": "callability_to_target_proof",
                            "source": ["callability.callability_id"],
                            "target": ["target_callability_evidence.callability_id"],
                            "cardinality": "one_to_one",
                            "completeness": "complete",
                            "unmatched_meaning": "an unmatched row is a validation failure",
                        },
                        {
                            "relationship_id": "target_callability_to_variant_evidence",
                            "source": ["target_callability_evidence.observed_variant_ids[]"],
                            "target": ["variants.variant_id"],
                            "cardinality": "many_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "the target was resolved by reference/no-call/interval evidence, "
                                "not by a carried core variant"
                            ),
                        },
                        {
                            "relationship_id": "target_callability_to_reference_evidence",
                            "source": [
                                "target_callability_evidence.selected_reference_evidence_ids[]"
                            ],
                            "target": ["reference_confidence_regions.evidence_id"],
                            "cardinality": "many_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "the target was resolved by variant, no-call, interval, or no evidence"
                            ),
                        },
                        {
                            "relationship_id": "target_callability_to_no_call_evidence",
                            "source": [
                                "target_callability_evidence.overlapping_no_call_evidence_ids[]"
                            ],
                            "target": ["reference_confidence_regions.evidence_id"],
                            "cardinality": "many_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": "the target had no overlapping no-call record",
                        },
                        {
                            "relationship_id": "target_callability_to_callable_region",
                            "source": ["target_callability_evidence.callable_region_ids[]"],
                            "target": ["callable_regions.region_id"],
                            "cardinality": "many_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": "no covering interval was needed or available",
                        },
                        {
                            "relationship_id": "target_callability_to_non_core_evidence",
                            "source": [
                                "target_callability_evidence.observed_non_core_evidence_ids[]"
                            ],
                            "target": ["non_core_variants.evidence_id"],
                            "cardinality": "many_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": "no symbolic or non-core carried call overlapped",
                        },
                    ]
                    if target_callability_summary is not None
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "pgx_locus_to_target_proof",
                            "source": ["pgx_locus_audit.target_ids[]"],
                            "target": ["target_callability_evidence.target_id"],
                            "cardinality": "many_to_many",
                            "completeness": "complete",
                            "unmatched_meaning": "an unmatched official PharmCAT allele is a validation failure",
                        },
                        {
                            "relationship_id": "pgx_locus_to_carried_variant",
                            "source": ["pgx_locus_audit.observed_variant_ids[]"],
                            "target": ["variants.variant_id"],
                            "cardinality": "many_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "the locus resolved as reference/no-call, or an overlapping "
                                "non-target variant is retained only as an audit identifier"
                            ),
                        },
                        *(
                            [
                                {
                                    "relationship_id": "pharmacogenomics_to_gene_call_audit",
                                    "source": ["pharmacogenomics.gene_symbol"],
                                    "target": ["pgx_gene_call_audit.gene_symbol"],
                                    "cardinality": "one_to_one",
                                    "completeness": "complete",
                                    "unmatched_meaning": "an unmatched normative PGx row is a validation failure",
                                }
                            ]
                            if pharmacogenomics_source is not None
                            and pharmacogenomics_source.is_file()
                            else []
                        ),
                        *(
                            [
                                {
                                    "relationship_id": "pgx_guidance_to_gene_call_audit",
                                    "source": ["pgx_guidance_audit.gene_symbol"],
                                    "target": ["pgx_gene_call_audit.gene_symbol"],
                                    "cardinality": "many_to_one",
                                    "completeness": "partial",
                                    "unmatched_meaning": (
                                        "guidance may describe a CPIC gene not emitted by the "
                                        "configured PharmCAT matcher"
                                    ),
                                }
                            ]
                            if pgx_guidance_audit_source is not None
                            and pgx_guidance_audit_source.is_file()
                            else []
                        ),
                    ]
                    if pgx_projection_summary is not None
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "pgs_model_to_locus_audit",
                            "source": ["pgs_model_audit.pgs_id"],
                            "target": ["pgs_locus_audit.pgs_id"],
                            "cardinality": "one_to_many",
                            "completeness": "complete",
                            "unmatched_meaning": "an unmatched configured model is a validation failure",
                        },
                        {
                            "relationship_id": "pgs_locus_to_target_proof",
                            "source": ["pgs_locus_audit.target_id"],
                            "target": ["target_callability_evidence.target_id"],
                            "cardinality": "one_to_one",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "the scoring row was unsupported and therefore has no materialized target"
                            ),
                        },
                        {
                            "relationship_id": "pgs_locus_to_carried_variant",
                            "source": ["pgs_locus_audit.variant_id"],
                            "target": ["variants.variant_id"],
                            "cardinality": "many_to_one",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "the locus was reference, missing, ambiguous, or unsupported"
                            ),
                        },
                        *(
                            [
                                {
                                    "relationship_id": "prs_to_pgs_model_audit",
                                    "source": ["prs.training_source"],
                                    "target": ["pgs_model_audit.pgs_id"],
                                    "cardinality": "one_to_one",
                                    "completeness": "complete",
                                    "unmatched_meaning": "an unmatched published score is a validation failure",
                                }
                            ]
                            if prs_source is not None and prs_source.is_file()
                            else []
                        ),
                    ]
                    if pgs_score_summary is not None
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "variant_to_clinical_finding",
                            "source": ["variants.variant_id"],
                            "target": ["clinical_findings.variant_id"],
                            "cardinality": "one_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "no pipeline-selected finding was materialized; this is not a "
                                "negative clinical result"
                            ),
                        },
                        {
                            "relationship_id": "clinical_finding_to_evidence",
                            "source": ["clinical_findings.evidence_ids[]"],
                            "target": ["clinical_evidence.evidence_id"],
                            "cardinality": "many_to_many",
                            "completeness": "complete",
                            "unmatched_meaning": "an unmatched evidence identifier is a validation failure",
                        },
                    ]
                    if clinical_findings_source.is_file()
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "clinical_finding_to_interpretation",
                            "source": ["clinical_findings.finding_id"],
                            "target": ["interpretations.finding_id"],
                            "cardinality": "one_to_many",
                            "completeness": "complete",
                            "unmatched_meaning": "an unmatched generated finding is a validation failure",
                        }
                    ]
                    if clinical_findings_source.is_file() and interpretations_source.is_file()
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "clinical_evidence_to_clinvar_assertion",
                            "source": ["clinical_evidence.evidence_id"],
                            "target": ["clinvar_assertions.evidence_id"],
                            "cardinality": "one_to_one",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "the evidence row came from another declared source, such as ACMG SF"
                            ),
                        }
                    ]
                    if clinvar_assertions_source.is_dir()
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "variant_to_acmg_secondary_finding_decision",
                            "source": ["variants.variant_id"],
                            "target": ["acmg_secondary_findings_evidence.variant_id"],
                            "cardinality": "one_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "no disease-specific ACMG policy candidate was materialized; "
                                "this is not a negative clinical result"
                            ),
                        }
                    ]
                    if acmg_evidence_source.is_file()
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "variant_to_gwas_association",
                            "source": ["variants.variant_id"],
                            "target": ["gwas_associations.variant_id"],
                            "cardinality": "one_to_many",
                            "completeness": "partial",
                            "unmatched_meaning": (
                                "no exact coordinate+rsID+effect-allele GWAS match; this is not "
                                "evidence that the variant has no trait association"
                            ),
                        },
                        {
                            "relationship_id": "gwas_association_to_lossless_evidence",
                            "source": ["gwas_associations.association_id"],
                            "target": ["gwas_catalog_evidence.association_id"],
                            "cardinality": "one_to_one",
                            "completeness": "complete",
                            "unmatched_meaning": (
                                "invalid: every compact typed row must have a lossless audit row"
                            ),
                        },
                    ]
                    if gwas_associations_source.is_dir()
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "clingen_gene_disease_to_gene_index",
                            "source": ["clingen_gene_disease.gene_symbol"],
                            "target": ["gene_index.gene_symbol"],
                            "cardinality": "many_to_many",
                            "completeness": "complete",
                            "unmatched_meaning": (
                                "invalid: every projected ClinGen curation must match an exact "
                                "gene symbol present in this bundle"
                            ),
                        }
                    ]
                    if clingen_gene_disease_source.is_file()
                    else []
                ),
                *(
                    [
                        {
                            "relationship_id": "clingen_dosage_to_gene_index",
                            "source": ["clingen_dosage_sensitivity.gene_symbol"],
                            "target": ["gene_index.gene_symbol"],
                            "cardinality": "many_to_many",
                            "completeness": "complete",
                            "unmatched_meaning": (
                                "invalid: every projected ClinGen dosage row must match an exact "
                                "gene symbol present in this bundle"
                            ),
                        }
                    ]
                    if clingen_dosage_source.is_file()
                    else []
                ),
            ],
        }
        atomic_write_json(staging / "relationships.json", relationships)

        asset_rows = {
            "variants.parquet": _parquet_rows(staging / "variants.parquet"),
            "callable_regions.parquet": _parquet_rows(staging / "callable_regions.parquet"),
            "extensions/variant_provenance.parquet": _parquet_rows(
                extensions / "variant_provenance.parquet"
            ),
            "extensions/reference_confidence_regions.parquet": _parquet_rows(
                extensions / "reference_confidence_regions.parquet"
            ),
        }
        if family_evidence_source is not None and family_evidence_source.is_dir():
            asset_rows["extensions/family_evidence.parquet"] = _parquet_rows(
                extensions / "family_evidence.parquet"
            )
        if gwas_associations_source.is_dir():
            asset_rows["gwas_associations.parquet"] = _parquet_rows(
                staging / "gwas_associations.parquet"
            )
            asset_rows["extensions/gwas_catalog_evidence.parquet"] = _parquet_rows(
                extensions / "gwas_catalog_evidence.parquet"
            )
        if callability_source is not None and callability_source.is_file():
            asset_rows["callability.parquet"] = _parquet_rows(staging / "callability.parquet")
        if (
            target_callability_evidence_source is not None
            and target_callability_evidence_source.is_file()
        ):
            asset_rows["extensions/target_callability_evidence.parquet"] = _parquet_rows(
                extensions / "target_callability_evidence.parquet"
            )
        if prs_source is not None and prs_source.is_file():
            asset_rows["prs.parquet"] = _parquet_rows(staging / "prs.parquet")
        if pgs_model_audit_source is not None and pgs_model_audit_source.is_file():
            asset_rows["extensions/pgs_model_audit.parquet"] = _parquet_rows(
                extensions / "pgs_model_audit.parquet"
            )
        if pgs_locus_audit_source is not None and pgs_locus_audit_source.is_file():
            asset_rows["extensions/pgs_locus_audit.parquet"] = _parquet_rows(
                extensions / "pgs_locus_audit.parquet"
            )
        if pharmacogenomics_source is not None and pharmacogenomics_source.is_file():
            asset_rows["pharmacogenomics.parquet"] = _parquet_rows(
                staging / "pharmacogenomics.parquet"
            )
        if pgx_locus_audit_source is not None and pgx_locus_audit_source.is_file():
            asset_rows["extensions/pgx_locus_audit.parquet"] = _parquet_rows(
                extensions / "pgx_locus_audit.parquet"
            )
        if pgx_gene_call_audit_source is not None and pgx_gene_call_audit_source.is_file():
            asset_rows["extensions/pgx_gene_call_audit.parquet"] = _parquet_rows(
                extensions / "pgx_gene_call_audit.parquet"
            )
        if pgx_guidance_audit_source is not None and pgx_guidance_audit_source.is_file():
            asset_rows["extensions/pgx_guidance_audit.parquet"] = _parquet_rows(
                extensions / "pgx_guidance_audit.parquet"
            )
        if gene_index_source.is_file():
            asset_rows["gene_index.parquet"] = _parquet_rows(staging / "gene_index.parquet")
        if clinical_evidence_source.is_dir():
            asset_rows["clinical_evidence.parquet"] = _parquet_rows(
                staging / "clinical_evidence.parquet"
            )
        if clinical_findings_source.is_file():
            asset_rows["clinical_findings.parquet"] = _parquet_rows(
                staging / "clinical_findings.parquet"
            )
        if interpretations_source.is_file():
            asset_rows["interpretations.parquet"] = _parquet_rows(
                staging / "interpretations.parquet"
            )
        if clinvar_assertions_source.is_dir():
            asset_rows["extensions/clinvar_assertions.parquet"] = _parquet_rows(
                extensions / "clinvar_assertions.parquet"
            )
        if acmg_evidence_source.is_file():
            asset_rows["extensions/acmg_secondary_findings_evidence.parquet"] = _parquet_rows(
                extensions / "acmg_secondary_findings_evidence.parquet"
            )
        if clingen_gene_disease_source.is_file():
            asset_rows["extensions/clingen_gene_disease.parquet"] = _parquet_rows(
                extensions / "clingen_gene_disease.parquet"
            )
        if clingen_dosage_source.is_file():
            asset_rows["extensions/clingen_dosage_sensitivity.parquet"] = _parquet_rows(
                extensions / "clingen_dosage_sensitivity.parquet"
            )
        if partial_rows:
            asset_rows["extensions/" + partial_name] = partial_rows
        if non_core_rows:
            asset_rows["extensions/non_core_variants.parquet"] = non_core_rows
        inventory_specs = [
            ("manifest.json", "application/json", "reproducibility anchor", True, "metadata"),
            ("schema.json", "application/schema+json", "normative schema", True, "metadata"),
            (
                "variants.parquet",
                "application/vnd.apache.parquet",
                "normalized carried alleles",
                True,
                "authoritative_data",
            ),
            (
                "callable_regions.parquet",
                "application/vnd.apache.parquet",
                "record-derived callability intervals",
                False,
                "derived_data",
            ),
            (
                "extensions/variant_provenance.parquet",
                "application/vnd.apache.parquet",
                "one-to-one source allele provenance",
                False,
                "authoritative_data",
            ),
            (
                "extensions/reference_confidence_regions.parquet",
                "application/vnd.apache.parquet",
                "reference and no-call evidence",
                False,
                "authoritative_data",
            ),
            (
                "extensions/reconciliation.json",
                "application/json",
                "source/output accounting",
                False,
                "metadata",
            ),
            ("field_glossary.json", "application/json", "field-level semantics", False, "metadata"),
            ("inventory.json", "application/json", "logical bundle inventory", False, "metadata"),
            (
                "raw_artifacts.json",
                "application/json",
                "source artifact disposition",
                False,
                "metadata",
            ),
            (
                "transformations.json",
                "application/json",
                "ordered transformation provenance",
                False,
                "metadata",
            ),
            (
                "omissions.json",
                "application/json",
                "explicit unsupported and absent domains",
                False,
                "metadata",
            ),
            ("relationships.json", "application/json", "declared joins", False, "metadata"),
        ]
        if family_evidence_source is not None and family_evidence_source.is_dir():
            inventory_specs.append(
                (
                    "extensions/family_evidence.parquet",
                    "application/vnd.apache.parquet",
                    "privacy-minimized pedigree, phase, LOCO, and concordance evidence",
                    False,
                    "derived_data",
                )
            )
        if gwas_associations_source.is_dir():
            inventory_specs.extend(
                [
                    (
                        "gwas_associations.parquet",
                        "application/vnd.apache.parquet",
                        "typed exact-allele GWAS associations with unambiguous effect type",
                        False,
                        "authoritative_data",
                    ),
                    (
                        "extensions/gwas_catalog_evidence.parquet",
                        "application/vnd.apache.parquet",
                        "lossless exact-match GWAS audit including untyped legacy rows",
                        False,
                        "authoritative_data",
                    ),
                ]
            )
        if callability_source is not None and callability_source.is_file():
            inventory_specs.extend(
                [
                    (
                        "callability.parquet",
                        "application/vnd.apache.parquet",
                        "target-level reference, variant, no-call, and unknown state",
                        False,
                        "derived_data",
                    ),
                    (
                        "extensions/target_callability_evidence.parquet",
                        "application/vnd.apache.parquet",
                        "target-catalogue and source-record proof chain",
                        False,
                        "authoritative_data",
                    ),
                ]
            )
        if prs_source is not None and prs_source.is_file():
            inventory_specs.append(
                (
                    "prs.parquet",
                    "application/vnd.apache.parquet",
                    "calibrated polygenic scores passing explicit completeness gates",
                    False,
                    "derived_data",
                )
            )
        if pgs_model_audit_source is not None and pgs_model_audit_source.is_file():
            inventory_specs.append(
                (
                    "extensions/pgs_model_audit.parquet",
                    "application/vnd.apache.parquet",
                    "polygenic model completeness and publication audit",
                    False,
                    "derived_data",
                )
            )
        if pgs_locus_audit_source is not None and pgs_locus_audit_source.is_file():
            inventory_specs.append(
                (
                    "extensions/pgs_locus_audit.parquet",
                    "application/vnd.apache.parquet",
                    "per-locus allele, genotype, missingness, and contribution audit",
                    False,
                    "derived_data",
                )
            )
        if pharmacogenomics_source is not None and pharmacogenomics_source.is_file():
            inventory_specs.append(
                (
                    "pharmacogenomics.parquet",
                    "application/vnd.apache.parquet",
                    "gated PharmCAT diplotype/phenotype calls with exact CPIC evidence level",
                    False,
                    "derived_data",
                )
            )
        if pgx_locus_audit_source is not None and pgx_locus_audit_source.is_file():
            inventory_specs.append(
                (
                    "extensions/pgx_locus_audit.parquet",
                    "application/vnd.apache.parquet",
                    "per-locus explicit variant/reference/no-call PharmCAT input proof",
                    False,
                    "derived_data",
                )
            )
        if pgx_gene_call_audit_source is not None and pgx_gene_call_audit_source.is_file():
            inventory_specs.append(
                (
                    "extensions/pgx_gene_call_audit.parquet",
                    "application/vnd.apache.parquet",
                    "PharmCAT gene-call publication decisions and withheld-call reasons",
                    False,
                    "derived_data",
                )
            )
        if pgx_guidance_audit_source is not None and pgx_guidance_audit_source.is_file():
            inventory_specs.append(
                (
                    "extensions/pgx_guidance_audit.parquet",
                    "application/vnd.apache.parquet",
                    "structured PharmCAT drug guidance linked to gated gene calls",
                    False,
                    "derived_data",
                )
            )
        if gene_index_source.is_file():
            inventory_specs.append(
                (
                    "gene_index.parquet",
                    "application/vnd.apache.parquet",
                    "gene-scoped query index",
                    False,
                    "derived_data",
                )
            )
        if clinical_evidence_source.is_dir():
            inventory_specs.append(
                (
                    "clinical_evidence.parquet",
                    "application/vnd.apache.parquet",
                    "normative versioned clinical evidence retrieval rows",
                    False,
                    "derived_data",
                )
            )
        if clinical_findings_source.is_file():
            inventory_specs.append(
                (
                    "clinical_findings.parquet",
                    "application/vnd.apache.parquet",
                    "pipeline-selected clinical finding retrieval index",
                    False,
                    "derived_data",
                )
            )
        if interpretations_source.is_file():
            inventory_specs.append(
                (
                    "interpretations.parquet",
                    "application/vnd.apache.parquet",
                    "structured algorithmic interpretations with explicit review status",
                    False,
                    "derived_data",
                )
            )
        if clinvar_assertions_source.is_dir():
            inventory_specs.append(
                (
                    "extensions/clinvar_assertions.parquet",
                    "application/vnd.apache.parquet",
                    "lossless exact-allele ClinVar assertions",
                    False,
                    "authoritative_data",
                )
            )
        if acmg_evidence_source.is_file():
            inventory_specs.append(
                (
                    "extensions/acmg_secondary_findings_evidence.parquet",
                    "application/vnd.apache.parquet",
                    "ACMG SF disease-aware policy candidate evidence",
                    False,
                    "derived_data",
                )
            )
        if clingen_gene_disease_source.is_file():
            inventory_specs.append(
                (
                    "extensions/clingen_gene_disease.parquet",
                    "application/vnd.apache.parquet",
                    "ClinGen expert-curated gene-disease validity context",
                    False,
                    "authoritative_data",
                )
            )
        if clingen_dosage_source.is_file():
            inventory_specs.append(
                (
                    "extensions/clingen_dosage_sensitivity.parquet",
                    "application/vnd.apache.parquet",
                    "ClinGen gene-level dosage-sensitivity context",
                    False,
                    "authoritative_data",
                )
            )
        if partial_rows:
            inventory_specs.append(
                (
                    "extensions/" + partial_name,
                    "application/vnd.apache.parquet",
                    "observed concrete alleles with incomplete genotypes and transcript consequences",
                    False,
                    "derived_data",
                )
            )
        if non_core_rows:
            inventory_specs.append(
                (
                    "extensions/non_core_variants.parquet",
                    "application/vnd.apache.parquet",
                    "non-core carried calls",
                    False,
                    "authoritative_data",
                )
            )
        if config.agent_exploration_profile:
            inventory_specs.extend(
                [
                    (
                        "readmygenome.md",
                        "text/markdown",
                        "agent entry point",
                        False,
                        "operational_guidance",
                    ),
                    (
                        "readmygenome/SKILL.md",
                        "text/markdown",
                        "agent operating guide",
                        False,
                        "operational_guidance",
                    ),
                    (
                        "readmygenome/scripts/query.py",
                        "text/x-python",
                        "bounded read-only query CLI",
                        False,
                        "operational_guidance",
                    ),
                    (
                        "readmygenome/reference/query_catalog.md",
                        "text/markdown",
                        "tested query catalogue",
                        False,
                        "operational_guidance",
                    ),
                    (
                        "readmygenome/reference/synonyms.json",
                        "application/json",
                        "retrieval synonym map",
                        False,
                        "operational_guidance",
                    ),
                    (
                        "readmygenome/evals/canonical_questions.json",
                        "application/json",
                        "expected-answer evaluation suite",
                        False,
                        "evaluation",
                    ),
                ]
            )
        inventory = {
            "schema_version": SCHEMA_VERSION,
            "entries": [
                _inventory_entry(
                    staging,
                    relative,
                    media,
                    role,
                    required,
                    authority,
                    asset_rows.get(relative),
                )
                for relative, media, role, required, authority in inventory_specs
            ],
        }
        atomic_write_json(staging / "inventory.json", inventory)

        files: dict[str, dict[str, Any]] = {}
        for item in sorted(staging.iterdir()):
            if item.name == "manifest.json":
                continue
            files[item.name] = {
                "sha256": directory_digest(item) if item.is_dir() else sha256_file(item),
                "rows": asset_rows.get(item.name),
                "bytes": _directory_bytes(item) if item.is_dir() else item.stat().st_size,
            }

        if vep_summary is not None:
            functional_annotation = (
                f"Ensembl VEP {vep_summary['vep_version']}; GRCh38 cache "
                f"{vep_summary['cache_version']}; tree-sha256:{vep_summary['cache_tree_sha256']}"
            )
            population_annotation = (
                "Ensembl VEP cached gnomAD frequencies; exact cache identity is recorded in "
                "resource_provenance"
            )
        else:
            functional_annotation = "not_run"
            population_annotation = "not_run"
        clinical_annotation = (
            "NCBI ClinVar "
            + preflight["resources"]["clinvar"]["snapshot"]
            + "; normalized-vcf-sha256:"
            + preflight["resources"]["clinvar"]["normalized_vcf_sha256"]
            if config.clinvar_enabled
            else "not_run"
        )
        actionability_annotation = (
            "ACMG SF v"
            + preflight["resources"]["acmg_sf"]["version"]
            + "; ClinGen snapshot "
            + preflight["resources"]["acmg_sf"]["snapshot"]
            + "; canonical-sha256:"
            + preflight["resources"]["acmg_sf"]["canonical_sha256"]
            if config.acmg_sf_enabled
            else "not_run"
        )
        resource_provenance: dict[str, Any] = {
            "build_run_key": run_key,
            "builder_implementation_sha256": preflight["software"]["builder_implementation_sha256"],
            "genome_spec_commit": GENOME_SPEC_COMMIT,
            "schema_sha256": PINNED_UPSTREAM_FILES["schema.json"],
            "reference_fasta_sha256": preflight["reference"]["sha256"],
            "reference_fai_sha256": preflight["reference_fai_sha256"],
        }
        if vep_summary is not None:
            resource_provenance["ensembl_vep"] = {
                "software_version": vep_summary["vep_version"],
                "cache_version": vep_summary["cache_version"],
                "cache_tree_sha256": vep_summary["cache_tree_sha256"],
                "cache_metadata": vep_summary["cache_metadata"],
                "container_image": vep_summary.get("container_image"),
            }
        if config.predictors.enabled:
            predictor_resource = preflight["resources"]["predictors"]
            semantic_controls = predictor_resource["semantic_controls"]
            resource_provenance["offline_predictors"] = {
                "manifest_sha256": predictor_resource["manifest_sha256"],
                "enabled_plugins": list(config.predictors.enabled_plugins),
                "semantic_controls": {
                    "control_set": semantic_controls["control_set"],
                    "control_definition_sha256": semantic_controls["control_definition_sha256"],
                    "doctor_run_key": semantic_controls["doctor_run_key"],
                    "report_sha256": semantic_controls["report_sha256"],
                    "record_count": semantic_controls["record_count"],
                    "record_identity_sha256": semantic_controls["record_identity_sha256"],
                    "reference_fasta_sha256": semantic_controls["reference_fasta_sha256"],
                    "vep_cache_manifest_sha256": semantic_controls["vep_cache_manifest_sha256"],
                    "control_results": semantic_controls["control_results"],
                    "semantic_boundary": semantic_controls["semantic_boundary"],
                },
            }
        if config.clinvar_enabled:
            resource_provenance["clinvar"] = {
                "snapshot": preflight["resources"]["clinvar"]["snapshot"],
                "archive_sha256": preflight["resources"]["clinvar"]["archive_sha256"],
                "normalized_vcf_sha256": preflight["resources"]["clinvar"]["normalized_vcf_sha256"],
                "normalized_records": preflight["resources"]["clinvar"]["normalized_records"],
            }
        if config.acmg_sf_enabled:
            acmg_summary = variant_summary["annotation_resources"]["acmg_sf"]
            resource_provenance["acmg_secondary_findings"] = {
                "version": preflight["resources"]["acmg_sf"]["version"],
                "snapshot": preflight["resources"]["acmg_sf"]["snapshot"],
                "canonical_sha256": preflight["resources"]["acmg_sf"]["canonical_sha256"],
                "source_workbook_sha256": preflight["resources"]["acmg_sf"][
                    "source_workbook_sha256"
                ],
                "policy_pmid": preflight["resources"]["acmg_sf"]["policy_pmid"],
                "policy_doi": preflight["resources"]["acmg_sf"]["policy_doi"],
                "counts": preflight["resources"]["acmg_sf"]["counts"],
                "minimum_clinvar_review_stars": acmg_summary["minimum_clinvar_review_stars"],
                "decision_rows": acmg_summary["decision_rows"],
                "actionable_candidate_variants": acmg_summary["actionable_candidate_variants"],
                "decision_status_counts": acmg_summary["status_counts"],
            }
        if config.clingen.enabled:
            clingen_resource = preflight["resources"]["clingen"]
            clingen_summary = variant_summary["annotation_resources"]["clingen"]
            resource_provenance["clingen_gene_level"] = {
                "snapshot": clingen_resource["snapshot"],
                "canonical_sha256": clingen_resource["canonical_sha256"],
                "manifest_sha256": clingen_resource["manifest_sha256"],
                "license": clingen_resource["license"],
                "redistribution_allowed": clingen_resource["redistribution_allowed"],
                "attribution": clingen_resource["attribution"],
                "sources": clingen_resource["sources"],
                "resource_outputs": clingen_resource["outputs"],
                "bundle_projection": {
                    "input_gene_symbols": clingen_summary["input_gene_symbols"],
                    "matched_gene_symbols": clingen_summary["matched_gene_symbols"],
                    "unmatched_gene_symbols": clingen_summary["unmatched_gene_symbols"],
                    "gene_disease_rows": clingen_summary["gene_disease_rows"],
                    "dosage_rows": clingen_summary["dosage_rows"],
                    "outputs": clingen_summary["outputs"],
                },
                "semantic_boundary": clingen_resource["semantic_boundary"],
                "variant_pathogenicity_claims": 0,
                "diagnosis_claims": 0,
            }
        if config.gwas.enabled:
            gwas_resource = preflight["resources"]["gwas_catalog"]
            resource_provenance["gwas_catalog"] = {
                "snapshot": gwas_resource["snapshot"],
                "source_url": gwas_resource["source_url"],
                "archive_sha256": gwas_resource["archive_sha256"],
                "archive_member_crc32": gwas_resource["archive_member_crc32"],
                "normalized_parquet_sha256": gwas_resource["normalized_parquet_sha256"],
                "license": gwas_resource["license"],
                "attribution": gwas_resource["attribution"],
                "source_archive_redistributed": False,
                "derived_rows_redistribution_allowed": config.gwas.redistribution_allowed,
                "matching_policy": gwas_resource["matching_policy"],
                "effect_type_policy": gwas_resource["effect_type_policy"],
                "resource_counts": gwas_resource["counts"],
                "bundle_match_counts": {
                    key: value
                    for key, value in variant_summary["counters"].items()
                    if key.startswith("gwas_")
                },
            }
        if target_callability_summary is not None:
            resource_provenance["target_catalog"] = {
                "sha256": target_callability_summary["catalog"]["sha256"],
                "bytes": target_callability_summary["catalog"]["bytes"],
                "target_rows": target_callability_summary["catalog"]["target_rows"],
                "allele_rows": target_callability_summary["catalog"]["allele_rows"],
                "normalization_checked": target_callability_summary["catalog"][
                    "normalization_checked"
                ],
                "assay_scope": target_callability_summary["catalog"]["assay_scope"],
                "precedence_rule": target_callability_summary["precedence_rule"],
                "resolution_counts": target_callability_summary["resolution_counts"],
            }
        if family_evidence_summary is not None:
            resource_provenance["family_evidence"] = {
                "source_version": family_evidence_summary["source_version"],
                "snapshot_date": family_evidence_summary["snapshot_date"],
                "joint_vcf_sha256": family_evidence_summary["joint_vcf_sha256"],
                "pedigree_sha256": family_evidence_summary["pedigree_sha256"],
                "subject_role": family_evidence_summary["subject_role"],
                "family_evidence_rows": family_evidence_summary["counters"].get(
                    "family_evidence_rows", 0
                ),
                "relative_identifiers_persisted": family_evidence_summary[
                    "relative_identifiers_persisted"
                ],
                "relative_genotypes_persisted": family_evidence_summary[
                    "relative_genotypes_persisted"
                ],
                "genotype_refinements": family_evidence_summary["genotype_refinements"],
                "normalization_projection": family_evidence_summary["normalization_projection"],
                "loco": family_evidence_summary.get("loco"),
            }
        if pgs_score_summary is not None:
            resource_provenance["polygenic_scores"] = {
                "source": "PGS Catalog",
                "scoring_format_version": "2.0",
                "target_catalog_sha256": pgs_score_summary["target_catalog_sha256"],
                "model_count": pgs_score_summary["model_count"],
                "published_model_count": pgs_score_summary["published_model_count"],
                "publication_status_counts": pgs_score_summary["publication_status_counts"],
                "models": [
                    {
                        "pgs_id": model.pgs_id,
                        "scoring_file_sha256": model.expected_sha256,
                        "source_version": model.source_version,
                        "training_date": model.training_date.isoformat(),
                        "trait": model.trait,
                        "model_type": model.model_type,
                        "license": model.license,
                        "redistribution_allowed": model.redistribution_allowed,
                        "max_missing_weight_fraction": model.max_missing_weight_fraction,
                        "calibration_sha256": (
                            model.calibration.expected_sha256 if model.calibration else None
                        ),
                    }
                    for model in config.pgs_models
                ],
            }
        if pgx_projection_summary is not None:
            pgx_resource = preflight["resources"]["pgx"]
            resource_provenance["pharmacogenomics"] = {
                "pharmcat_version": PHARMCAT_VERSION,
                "pharmcat_release_commit": config.pgx.release_commit,
                "pharmcat_positions_vcf_sha256": pgx_resource["pharmcat"]["positions_vcf_sha256"],
                "pharmcat_positions_records": pgx_resource["pharmcat"]["record_count"],
                "pharmcat_positions_alleles": pgx_resource["pharmcat"]["allele_count"],
                "pharmcat_normalization_changed_alleles": pgx_resource["pharmcat"][
                    "normalization_changed_alleles"
                ],
                "pharmcat_container_image": config.pgx.container_image,
                "pharmcat_license": config.pgx.license,
                "pharmcat_redistribution_allowed": config.pgx.redistribution_allowed,
                "pharmcat_data_version": pharmcat_summary["data_version"],
                "canonical_report_sha256": pharmcat_summary["canonical_report_sha256"],
                "calls_tsv_sha256": pharmcat_summary["calls_tsv_sha256"],
                "cpic_source": pgx_resource["cpic"]["source"],
                "cpic_source_version": pgx_resource["cpic"]["source_version"],
                "cpic_snapshot_date": pgx_resource["cpic"]["snapshot_date"],
                "cpic_snapshot_sha256": pgx_resource["cpic"]["sha256"],
                "cpic_rowset_canonical_sha256": pgx_resource["cpic"]["rowset_canonical_sha256"],
                "cpic_source_release": pgx_resource["cpic"]["source_release"],
                "cpic_api_identity": pgx_resource["cpic"]["api_identity"],
                "cpic_license": config.pgx.cpic_license,
                "cpic_redistribution_allowed": config.pgx.cpic_redistribution_allowed,
                "cpic_normative_recommendation_rows": pgx_resource["cpic"][
                    "normative_recommendation_rows"
                ],
                "materialization_status_counts": pgx_materialization_summary["status_counts"],
                "published_gene_rows": pgx_projection_summary["published_gene_rows"],
                "gene_call_status_counts": pgx_projection_summary["call_status_counts"],
                "complex_matcher_genes_excluded": sorted(STRUCTURALLY_COMPLEX_PGX_GENES),
                "missing_to_reference_conversion": False,
                "network_during_matching": "none",
            }

        manifest: dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "pipeline_version": f"genome-annotation-builder-{PIPELINE_VERSION}",
            "source_bundle_sha256": preflight["input"]["sha256"],
            "genome_build": "GRCh38",
            "generated_at": generated_at,
            "annotation_sources": {
                "genome_spec": f"{GENOME_SPEC_REPOSITORY}@{GENOME_SPEC_COMMIT}",
                "functional_annotation": functional_annotation,
                "population_annotation": population_annotation,
                "clinical_annotation": clinical_annotation,
                "actionability_annotation": actionability_annotation,
                "polygenic_scores": (
                    "checksum-pinned PGS Catalog models with explicit locus and calibration audits"
                    if pgs_score_summary is not None
                    else "not_run"
                ),
                "pharmacogenomics": (
                    f"PharmCAT {PHARMCAT_VERSION} offline with frozen CPIC "
                    f"{config.pgx.cpic_source_version} pair snapshot"
                    if pgx_projection_summary is not None
                    else "not_run"
                ),
                "family_evidence": (
                    "digest-pinned joint VCF/PED with privacy-minimized subject-only projection"
                    if family_evidence_summary is not None
                    else "not_run"
                ),
                "gwas_catalog": (
                    f"NHGRI-EBI GWAS Catalog {config.gwas.snapshot}; exact allele projection"
                    if config.gwas.enabled
                    else "not_run"
                ),
                "clingen_gene_level": (
                    f"ClinGen gene-disease validity and dosage sensitivity "
                    f"{config.clingen.snapshot}; exact-symbol gene context only"
                    if config.clingen.enabled
                    else "not_run"
                ),
            },
            "consensus_definitions": {
                "is_rare": {
                    "rule": (
                        "true when gnomad_af_popmax is missing or < 0.01; missing is "
                        "conservatively retained and is not evidence of observed rarity"
                    ),
                    "threshold": 0.01,
                    "population": "gnomad_af_popmax",
                    "null_treatment": "rare",
                },
                "is_actionable": (
                    {
                        "rule": (
                            "true only for an automated ACMG SF reportable candidate with an "
                            "exact disease-mapped ClinVar Pathogenic/Likely_pathogenic assertion, "
                            "no conflict, sufficient review stars, a high-confidence genotype, "
                            "and all implemented inheritance/variant-class gates; gene membership "
                            "alone never qualifies"
                        ),
                        "gene_list": f"ACMG_SF_v{config.acmg_sf_version}",
                        "gene_list_sha256": ACMG_SF_GENE_SET_SHA256,
                        "policy_snapshot": config.acmg_sf_snapshot,
                        "policy_canonical_sha256": preflight["resources"]["acmg_sf"][
                            "canonical_sha256"
                        ],
                        "significance_set": ["Pathogenic", "Likely_pathogenic"],
                        "min_review_stars": config.acmg_sf_min_clinvar_review_stars,
                        "conflicts_allowed": False,
                        "analytical_confidence": "high",
                        "compound_heterozygous_rule": (
                            "not actionable unless trans phase is independently established; "
                            "the current evaluator emits PHASE_REVIEW_REQUIRED"
                        ),
                        "clinical_confirmation_required": True,
                    }
                    if config.acmg_sf_enabled
                    else {
                        "rule": (
                            "false because no versioned actionability policy and gene list were applied"
                        ),
                        "gene_list": "EMPTY_NO_ACTIONABILITY_POLICY",
                        "gene_list_sha256": EMPTY_CANONICAL_LIST_SHA256,
                        "significance_set": ["Pathogenic", "Likely_pathogenic"],
                        "min_review_stars": 2,
                    }
                ),
                "consensus_damaging": {
                    "rule": "at least two populated predictor scores exceed thresholds; missing scores cast no vote",
                    "applies_to": ["missense_variant"],
                    "tools": [
                        {"name": "REVEL", "field": "predictions.revel", "threshold": 0.75},
                        {
                            "name": "AlphaMissense",
                            "field": "predictions.alphamissense",
                            "threshold": 0.564,
                        },
                        {"name": "CADD", "field": "predictions.cadd_phred", "threshold": 20.0},
                    ],
                    "aggregation": "at_least_n",
                    "n": 2,
                },
                "clinical_grade": {
                    "rule": (
                        "false because no pipeline-reviewed interpretation policy was applied; "
                        "raw ClinVar evidence alone is not clinical-grade"
                    ),
                    "qualifying_sources": [],
                },
                "target_reference_observed": {
                    "rule": (
                        "true only when one or more complete explicit-reference source records "
                        "cover the full normalized target span, no carried core/non-core variant "
                        "overlaps it, and no overlapping no-call record supersedes reference evidence"
                    ),
                    "precedence": (
                        target_callability_summary["precedence_rule"]
                        if target_callability_summary is not None
                        else "not_run_no_target_catalog"
                    ),
                    "callable_min_dp": config.quality.callable_min_dp,
                    "callable_min_gq": config.quality.callable_min_gq,
                    "interval_callability_is_reference": False,
                },
                "polygenic_score": {
                    "rule": (
                        "additive effect-allele dosage is derived only from an exact passing "
                        "carried genotype or complete passing explicit-reference evidence; "
                        "no-call, interval-only, unknown, low-confidence, unsupported, and "
                        "other-ALT loci contribute no invented dosage"
                    ),
                    "missing_genotypes_imputed": False,
                    "silent_reference_fill": False,
                    "percentile_requires_digest_bound_population_calibration": True,
                    "normative_prs_requires_score_valid_and_calibrated": True,
                },
                "pharmacogenomic_call": {
                    "rule": (
                        "publish only one complete unambiguous non-complex PharmCAT diplotype "
                        "with phenotype and at least one non-provisional CPIC recommendation at "
                        "an exact A-D level; compact core drugs come only from the strongest level"
                    ),
                    "locus_reference_requires_explicit_reference_evidence": True,
                    "missing_to_reference_conversion": False,
                    "interval_callability_is_reference": False,
                    "other_alt_is_reference": False,
                    "structurally_complex_vcf_only_calls_published": False,
                    "cpic_levels": ["A", "B", "C", "D"],
                    "clinical_confirmation_required": True,
                },
                "family_evidence": {
                    "rule": (
                        "family evidence may add a subject-only confidence or review status but "
                        "never changes CHROM/POS/REF/ALT/GT, ploidy, or phase"
                    ),
                    "relative_identifiers_embedded": False,
                    "relative_genotypes_embedded": False,
                    "missing_loco_exception_means_consistent": False,
                    "genotype_refinements": 0,
                },
                "call_confidence": {
                    "rule": "technical confidence from source PASS/FT plus DP and GQ thresholds",
                    "high_min_dp": config.quality.high_min_dp,
                    "high_min_gq": config.quality.high_min_gq,
                    "medium_min_dp": config.quality.medium_min_dp,
                    "medium_min_gq": config.quality.medium_min_gq,
                    "core_physical_projection": (
                        "DP > 32767 or GQ > 127 becomes null only in the normative core field; "
                        "the exact source value remains in extensions/variant_provenance.parquet"
                    ),
                },
                "callable": {
                    "rule": "complete GT, passing record/genotype filters, effective DP and GQ thresholds",
                    "minimum_dp": config.quality.callable_min_dp,
                    "minimum_gq": config.quality.callable_min_gq,
                    "reference_block_depth_field_preference": ["MIN_DP", "DP"],
                },
            },
            "normalization": (
                "bcftools norm --no-version -f GRCh38.fa --atomize --check-ref "
                f"{'e' if config.check_ref == 'error' else 'w'}; bcftools sort"
            ),
            "profiles": (
                {"agent_exploration": "1.0.0"} if config.agent_exploration_profile else {}
            ),
            "files": files,
            "resource_provenance": resource_provenance,
            "privacy": {
                "input_paths_embedded": config.include_input_paths,
                "sample_identifier_embedded": config.include_sample_identifier,
            },
        }
        atomic_write_json(staging / "manifest.json", manifest)

        from .validator import validate_bundle

        validation = validate_bundle(
            staging,
            reference_fasta=config.reference_fasta,
            expected_schema_path=schema_source,
            run_profile_queries=True,
        )
        if not validation["valid"]:
            raise ValueError("staged bundle validation did not pass")
    return validation
