"""Convert normalized projected calls into typed core and provenance datasets."""

from __future__ import annotations

import base64
import re
from collections import Counter
from contextlib import ExitStack
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from . import native_observation
from .acmg_evaluator import (
    AcmgSecondaryFindingsEvaluator,
    materialize_acmg_clinical_artifacts,
)
from .acmg_resources import load_acmg_sf
from .arrow_schemas import (
    acmg_secondary_findings_evidence_schema,
    clinical_evidence_schema,
    clinical_finding_schema,
    clinvar_assertion_schema,
    gene_index_schema,
    gwas_association_schema,
    gwas_catalog_evidence_schema,
    interpretation_schema,
    observed_partial_variant_schema,
    variant_provenance_schema,
    variant_schema,
)
from .clingen_project import project_clingen_gene_context
from .clinvar_parser import SortedClinVarLookup, aggregate_clinvar
from .clinvar_resources import verify_clinvar
from .config import BuildConfig
from .errors import InputValidationError, ReconciliationError
from .gwas_parser import SortedGwasCatalogLookup
from .gwas_resources import verify_gwas_catalog
from .parquet_io import PartitionedParquetWriter
from .util import (
    MultisetDigest,
    atomic_output_directory,
    atomic_write_json,
    directory_digest,
    sha256_file,
)
from .vcf_ingest import canonical_core_chrom, source_id_from_digest
from .vep_parser import VepParser, apply_vep_annotation

SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _snapshot_datetime(snapshot: str) -> datetime:
    try:
        return datetime.strptime(snapshot, "%Y%m%d").replace(tzinfo=UTC)
    except ValueError as exc:
        raise InputValidationError(f"invalid resource snapshot date: {snapshot!r}") from exc


def _standard_clinvar_evidence(row: dict[str, Any]) -> dict[str, Any]:
    assertion = row.get("aggregate_significance") or "; ".join(row["raw_significance"])
    if not assertion:
        assertion = "classification unavailable"
    return {
        "evidence_id": row["evidence_id"],
        "finding_id": None,
        "variant_id": row["variant_id"],
        "source": "NCBI ClinVar",
        "source_record_id": row["clinvar_accession"],
        "source_version": row["snapshot"],
        "assertion": assertion,
        "review_status": "; ".join(row["review_status"]) or None,
        "retrieved_at": _snapshot_datetime(str(row["snapshot"])),
    }


def _evaluate_acmg_candidates(config: BuildConfig, normalized_vcf: Path) -> dict[str, Any] | None:
    if not config.acmg_sf_enabled:
        return None
    evaluator = AcmgSecondaryFindingsEvaluator(
        load_acmg_sf(config),
        min_review_stars=config.acmg_sf_min_clinvar_review_stars,
    )
    source = pysam.VariantFile(str(normalized_vcf))
    try:
        parser = VepParser.from_header(source.header)
        if not parser.enabled:
            raise InputValidationError("ACMG SF evaluation requires typed VEP annotations")
        parser.require_pipeline_fields()
        with SortedClinVarLookup(config.clinvar_vcf, config.clinvar_snapshot) as lookup:
            for record in source:
                if not record.alts or len(record.alts) != 1:
                    raise InputValidationError("ACMG input record is not normalized biallelic")
                chrom = canonical_core_chrom(record.contig)
                if chrom is None:
                    raise InputValidationError("non-core contig reached ACMG evaluator")
                ref = record.ref.upper()
                alt = record.alts[0].upper()
                annotation = parser.parse(record)
                if annotation is None:
                    raise InputValidationError(
                        f"ACMG input lacks VEP allele annotation at {record.contig}:{record.pos}"
                    )
                transcripts = list(annotation["consequence"]["all_annotations"])
                symbols = {
                    str(item["symbol"]).upper()
                    for item in transcripts
                    if item.get("symbol") is not None
                }
                selected_symbol = annotation["gene"].get("symbol")
                if selected_symbol:
                    symbols.add(str(selected_symbol).upper())
                if not symbols & evaluator.genes:
                    continue
                sample = record.samples["SAMPLE"]
                raw_gt = sample.get("GT")
                gt = tuple(raw_gt) if raw_gt is not None else ()
                if None in gt and (
                    1 in gt or (gt == (None,) and native_observation.PROJECTION_FLAG in record.info)
                ):
                    # Allele annotations are available separately for an incomplete genotype.
                    continue
                if len(gt) not in {1, 2} or any(allele not in {0, 1} for allele in gt):
                    raise InputValidationError("invalid genotype reached ACMG evaluator")
                dp_value = sample.get("DP")
                gq_value = sample.get("GQ")
                dp = int(dp_value) if dp_value is not None else None
                gq = int(gq_value) if gq_value is not None else None
                filters = _source_filters(record)
                passed = all(item == "PASS" for item in filters)
                confidence = config.quality.confidence(dp, gq) if passed else "low"
                records = lookup.lookup(chrom, record.pos, ref, alt)
                if not records:
                    continue
                _pathogenicity, evidence_rows, _rsid = aggregate_clinvar(
                    records, f"{chrom}:{record.pos}:{ref}:{alt}", config.clinvar_snapshot
                )
                evaluator.observe(
                    variant_id=f"{chrom}:{record.pos}:{ref}:{alt}",
                    gene_symbols=symbols,
                    gt=tuple(int(allele) for allele in gt),
                    phased=bool(sample.phased),
                    call_confidence=confidence,
                    transcripts=transcripts,
                    clinvar_rows=evidence_rows,
                )
    finally:
        source.close()
    return evaluator.finalize()


def _info_scalar(record: pysam.VariantRecord, key: str) -> Any:
    value = record.info.get(key)
    if isinstance(value, tuple):
        if len(value) != 1:
            raise InputValidationError(f"normalized INFO/{key} is not scalar")
        return value[0]
    return value


def _info_integer_list(record: pysam.VariantRecord, key: str) -> list[int] | None:
    value = record.info.get(key)
    if value is None:
        return None
    raw = value if isinstance(value, tuple) else (value,)
    return [int(item) for item in raw]


def _decode_text(value: Any) -> str | None:
    if value is None:
        return None
    encoded = str(value)
    padding = "=" * (-len(encoded) % 4)
    try:
        return base64.urlsafe_b64decode(encoded + padding).decode("utf-8")
    except (ValueError, UnicodeDecodeError) as exc:
        raise InputValidationError("invalid base64url provenance metadata") from exc


def _variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt):
        return "INS"
    if len(ref) > len(alt):
        return "DEL"
    return "MNV"


def _zygosity(gt: tuple[int, ...]) -> str:
    if all(allele == 0 for allele in gt):
        return "hom_ref"
    if all(allele == 1 for allele in gt):
        return "hom_alt"
    return "het"


def _source_filters(record: pysam.VariantRecord) -> list[str]:
    value = _info_scalar(record, "GAB_SOURCE_FILTERS")
    if value in {None, "", ".", "PASS"}:
        return ["PASS"]
    return [item for item in str(value).split("|") if item]


def _blank_variant_row(
    *,
    variant_id: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    rsid: str | None,
    call_confidence: str,
    dp: int | None,
    gq: int | None,
    gt: tuple[int, ...],
    phased: bool,
    allele_balance: float | None,
) -> dict[str, Any]:
    return {
        "variant_id": variant_id,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "variant_type": _variant_type(ref, alt),
        "rsid": rsid,
        "quality": {
            "call_confidence": call_confidence,
            "origin": "observed",
            "imputation_r2": None,
            "imputation_panel": None,
            "dp": dp,
            "gq": gq,
        },
        "genotype": {
            "gt": list(gt),
            "phased": phased,
            "zygosity": _zygosity(gt),
            "allele_balance": allele_balance,
        },
        "gene": {
            "symbol": None,
            "ensembl_id": None,
            "hgnc_id": None,
            "biotype": None,
            "is_canonical": None,
        },
        "consequence": {
            "most_severe": None,
            "impact": None,
            "hgvsc": None,
            "hgvsp": None,
            "all_annotations": None,
        },
        "population_frequency": {
            "gnomad_af": None,
            "gnomad_af_popmax": None,
            "gnomad_af_nfe": None,
            "gnomad_af_eas": None,
            "gnomad_af_afr": None,
            "gnomad_af_sas": None,
            "gnomad_af_amr": None,
            # The normative schema requires a boolean. Missing population data
            # is retained under the conservative discovery rule, never treated
            # as evidence that an allele is common.
            "is_rare": True,
        },
        "pathogenicity": {
            "clinvar_significance": None,
            "clinvar_has_conflicts": False,
            "clinvar_conflict_summary": None,
            "clinvar_review_stars": None,
            "clinvar_submitters_count": None,
            "clinvar_id": None,
            "is_actionable": False,
        },
        "predictions": {
            "cadd_phred": None,
            "revel": None,
            "alphamissense": None,
            "spliceai_max": None,
            "polyphen": None,
            "sift": None,
            "consensus_damaging": False,
        },
        "pharmacogenomics": {"is_pgx": False, "contributes_to_gene": None},
        "trait_associations": {
            "is_gwas_hit": False,
            "traits": None,
            "effect_size": None,
            "effect_type": None,
            "effect_allele": None,
            "ci_lower": None,
            "ci_upper": None,
            "p_value": None,
            "sample_size": None,
            "study_pmids": None,
        },
        "clinical_grade": False,
    }


LOFTEE_APPLICABLE_TERMS = frozenset(
    {
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
    }
)


def _count_predictor_completeness(
    config: BuildConfig,
    row: dict[str, Any],
    vep_annotation: dict[str, Any],
    counters: Counter[str],
) -> None:
    transcripts = vep_annotation["consequence"]["all_annotations"]
    terms = {term for transcript in transcripts for term in transcript["consequence_terms"] or []}
    predictions = row["predictions"]
    variant_type = row["variant_type"]

    applicability = {
        "loftee": bool(terms & LOFTEE_APPLICABLE_TERMS),
        "cadd": variant_type in {"SNV", "INS", "DEL"},
        "revel": "missense_variant" in terms and variant_type == "SNV",
        "alphamissense": "missense_variant" in terms and variant_type == "SNV",
        "spliceai": variant_type == "SNV"
        and any(transcript["mane_select"] for transcript in transcripts),
    }
    observed = {
        "loftee": any(transcript["lof"] is not None for transcript in transcripts),
        "cadd": predictions["cadd_phred"] is not None,
        "revel": predictions["revel"] is not None,
        "alphamissense": predictions["alphamissense"] is not None,
        "spliceai": predictions["spliceai_max"] is not None,
    }
    for plugin in config.predictors.enabled_plugins:
        if not applicability[plugin]:
            continue
        counters[f"predictor_{plugin}_eligible_variants"] += 1
        status = "scored" if observed[plugin] else "missing"
        counters[f"predictor_{plugin}_{status}_variants"] += 1
    if config.predictors.spliceai.enabled and variant_type in {"INS", "DEL", "MNV"}:
        counters["predictor_spliceai_indel_or_mnv_unavailable_variants"] += 1
    if config.predictors.cadd.enabled and variant_type == "MNV":
        counters["predictor_cadd_mnv_outside_resource_scope_variants"] += 1


def build_variant_datasets(
    config: BuildConfig,
    normalized_vcf: str | Path,
    input_sha256: str,
    expected_multiset: dict[str, Any],
    destination: str | Path,
) -> dict[str, Any]:
    normalized_vcf = Path(normalized_vcf)
    destination = Path(destination)
    with atomic_output_directory(destination) as staging:
        variants_path = staging / "variants.parquet"
        provenance_path = staging / "variant_provenance.parquet"
        partial_path = staging / "observed_partial_variants.parquet"
        counters: Counter[str] = Counter()
        gene_stats: dict[tuple[str, str], dict[str, Any]] = {}
        clinvar_manifest = verify_clinvar(config) if config.clinvar_enabled else None
        gwas_manifest = verify_gwas_catalog(config) if config.gwas.enabled else None
        acmg_evaluation = _evaluate_acmg_candidates(config, normalized_vcf)
        clinical_artifacts = (
            materialize_acmg_clinical_artifacts(acmg_evaluation)
            if acmg_evaluation is not None
            else {"clinical_findings": [], "policy_evidence": [], "interpretations": []}
        )
        actionable_variant_ids = (
            set(acmg_evaluation["actionable_variant_ids"]) if acmg_evaluation is not None else set()
        )
        actionable_variant_gene_pairs = (
            {
                (str(row["variant_id"]), str(row["gene_symbol"]))
                for row in acmg_evaluation["actionable_variant_gene_pairs"]
            }
            if acmg_evaluation is not None
            else set()
        )

        source = pysam.VariantFile(str(normalized_vcf))
        if list(source.header.samples) != ["SAMPLE"]:
            source.close()
            raise InputValidationError("normalized projected VCF must contain exactly SAMPLE")
        vep_parser = VepParser.from_header(source.header)
        if config.vep_enabled and not vep_parser.enabled:
            source.close()
            raise InputValidationError("VEP was enabled but normalized input has no CSQ annotation")
        if config.vep_enabled:
            vep_parser.require_pipeline_fields()
        fasta = pysam.FastaFile(str(config.reference_fasta))

        previous_locus: tuple[str, int] | None = None
        locus_ids: set[str] = set()
        clinical_evidence_path = staging / "clinical_evidence.parquet"
        clinvar_assertions_path = staging / "clinvar_assertions.parquet"
        clinical_findings_path = staging / "clinical_findings.parquet"
        interpretations_path = staging / "interpretations.parquet"
        gwas_associations_path = staging / "gwas_associations.parquet"
        gwas_evidence_path = staging / "gwas_catalog_evidence.parquet"
        try:
            with ExitStack() as stack:
                variant_writer = stack.enter_context(
                    PartitionedParquetWriter(variants_path, variant_schema())
                )
                partial_writer = stack.enter_context(
                    PartitionedParquetWriter(partial_path, observed_partial_variant_schema())
                )
                provenance_writer = stack.enter_context(
                    PartitionedParquetWriter(
                        provenance_path, variant_provenance_schema(), partition_field=None
                    )
                )
                clinvar_lookup = (
                    stack.enter_context(
                        SortedClinVarLookup(config.clinvar_vcf, config.clinvar_snapshot)
                    )
                    if config.clinvar_enabled
                    else None
                )
                clinical_evidence_writer = (
                    stack.enter_context(
                        PartitionedParquetWriter(
                            clinical_evidence_path,
                            clinical_evidence_schema(),
                            partition_field=None,
                        )
                    )
                    if config.clinvar_enabled
                    else None
                )
                clinvar_assertions_writer = (
                    stack.enter_context(
                        PartitionedParquetWriter(
                            clinvar_assertions_path,
                            clinvar_assertion_schema(),
                            partition_field=None,
                        )
                    )
                    if config.clinvar_enabled
                    else None
                )
                gwas_lookup = (
                    stack.enter_context(
                        SortedGwasCatalogLookup(
                            config.gwas.normalized_parquet,
                            config.gwas.snapshot,
                        )
                    )
                    if config.gwas.enabled
                    else None
                )
                gwas_associations_writer = (
                    stack.enter_context(
                        PartitionedParquetWriter(
                            gwas_associations_path,
                            gwas_association_schema(),
                            partition_field=None,
                        )
                    )
                    if config.gwas.enabled
                    else None
                )
                gwas_evidence_writer = (
                    stack.enter_context(
                        PartitionedParquetWriter(
                            gwas_evidence_path,
                            gwas_catalog_evidence_schema(),
                            partition_field=None,
                        )
                    )
                    if config.gwas.enabled
                    else None
                )
                if clinical_evidence_writer is not None:
                    for evidence in clinical_artifacts["policy_evidence"]:
                        clinical_evidence_writer.write(evidence)
                    counters["acmg_policy_evidence_rows"] = len(
                        clinical_artifacts["policy_evidence"]
                    )
                for record in source:
                    counters["normalized_records"] += 1
                    if not record.alts or len(record.alts) != 1:
                        raise InputValidationError("normalized projected record is not biallelic")
                    ref = record.ref.upper()
                    alt = record.alts[0].upper()
                    if not re.fullmatch(r"[ACGT]+", ref) or not re.fullmatch(r"[ACGT]+", alt):
                        raise InputValidationError(
                            "symbolic or non-ACGT allele reached core variant builder"
                        )
                    chrom = canonical_core_chrom(record.contig)
                    if chrom is None:
                        raise InputValidationError("non-core contig reached core variant builder")
                    if record.contig not in fasta.references:
                        raise InputValidationError(
                            f"normalized contig absent from FASTA: {record.contig}"
                        )
                    observed_ref = fasta.fetch(
                        record.contig, record.pos - 1, record.pos - 1 + len(ref)
                    ).upper()
                    if observed_ref != ref:
                        raise InputValidationError(
                            f"normalized REF mismatch at {record.contig}:{record.pos}: "
                            f"{ref} != {observed_ref}"
                        )

                    sample = record.samples["SAMPLE"]
                    raw_gt = sample.get("GT")
                    gt = tuple(raw_gt) if raw_gt is not None else ()
                    if len(gt) not in {1, 2} or any(allele not in {None, 0, 1} for allele in gt):
                        raise InputValidationError("invalid projected genotype after normalization")
                    native_no_call = native_observation.PROJECTION_FLAG in record.info
                    if native_no_call and gt != (None,):
                        raise InputValidationError("native no-call observation acquired a genotype")
                    if 1 not in gt and not native_no_call:
                        raise InputValidationError(
                            "homozygous-reference call reached variants.parquet"
                        )

                    source_id = str(_info_scalar(record, "GAB_SOURCE_ID"))
                    source_digest = str(_info_scalar(record, "GAB_SOURCE_DIGEST"))
                    source_ordinal = int(_info_scalar(record, "GAB_SOURCE_ORDINAL"))
                    source_alt_index = int(_info_scalar(record, "GAB_ALT_INDEX"))
                    if not SHA256_PATTERN.fullmatch(source_id) or not SHA256_PATTERN.fullmatch(
                        source_digest
                    ):
                        raise InputValidationError(
                            "invalid source digest metadata after normalization"
                        )
                    expected_source_id = source_id_from_digest(
                        input_sha256, source_ordinal, source_digest
                    )
                    if source_id != expected_source_id:
                        raise ReconciliationError(
                            "normalized record source identifier is not reproducible"
                        )

                    variant_id = f"{chrom}:{record.pos}:{ref}:{alt}"
                    vep_annotation = vep_parser.parse(record) if config.vep_enabled else None
                    if config.vep_enabled and vep_annotation is None:
                        raise InputValidationError(
                            f"VEP output record has no allele annotation at {record.contig}:{record.pos}"
                        )
                    if None in gt:
                        # Keep allele-level VEP results outside the normative complete-GT
                        # variants and all genotype-dependent clinical/trait consumers.
                        partial_writer.write(
                            {
                                "variant_id": variant_id,
                                "chrom": chrom,
                                "pos": record.pos,
                                "ref": ref,
                                "alt": alt,
                                "rsid": _info_scalar(record, "GAB_RSID")
                                or (vep_annotation or {}).get("rsid"),
                                "gt": list(gt),
                                "phased": False if native_no_call else bool(sample.phased),
                                "source_id": source_id,
                                "source_record_digest": source_digest,
                                "source_ordinal": source_ordinal,
                                "source_alt_index": source_alt_index,
                                "source_chrom": _info_scalar(record, "GAB_ORIG_CONTIG"),
                                "source_pos": _info_scalar(record, "GAB_ORIG_POS"),
                                "source_ref": _info_scalar(record, "GAB_ORIG_REF"),
                                "source_alt": _info_scalar(record, "GAB_ORIG_ALT"),
                                "source_gt": _info_scalar(record, "GAB_ORIG_GT"),
                                "gene": (vep_annotation or {}).get("gene", {}),
                                "consequence": (vep_annotation or {}).get("consequence", {}),
                            }
                        )
                        counters["observed_partial_variant_rows"] += 1
                        continue
                    locus = (chrom, record.pos)
                    if locus != previous_locus:
                        locus_ids.clear()
                        previous_locus = locus
                    if variant_id in locus_ids:
                        raise ReconciliationError(f"duplicate normalized variant: {variant_id}")
                    locus_ids.add(variant_id)

                    dp_value = sample.get("DP")
                    gq_value = sample.get("GQ")
                    dp = int(dp_value) if dp_value is not None else None
                    gq = int(gq_value) if gq_value is not None else None
                    if dp is not None and dp < 0:
                        raise InputValidationError("negative DP reached core variant builder")
                    if gq is not None and gq < 0:
                        raise InputValidationError("negative GQ reached core variant builder")
                    core_dp = dp if dp is None or dp <= 32_767 else None
                    core_gq = gq if gq is None or gq <= 127 else None
                    if dp is not None and core_dp is None:
                        counters["core_dp_out_of_int16_range_preserved_in_provenance"] += 1
                    if gq is not None and core_gq is None:
                        counters["core_gq_out_of_int8_range_preserved_in_provenance"] += 1
                    filters = _source_filters(record)
                    passed = all(item == "PASS" for item in filters)
                    confidence = config.quality.confidence(dp, gq) if passed else "low"
                    alt_ad_value = _info_scalar(record, "GAB_ALT_AD")
                    total_ad_value = _info_scalar(record, "GAB_TOTAL_AD")
                    allele_balance = None
                    if alt_ad_value is not None and total_ad_value is not None:
                        alt_ad = int(alt_ad_value)
                        total_ad = int(total_ad_value)
                        if total_ad > 0 and 0 <= alt_ad <= total_ad:
                            allele_balance = alt_ad / total_ad

                    rsid_value = _info_scalar(record, "GAB_RSID")
                    row = _blank_variant_row(
                        variant_id=variant_id,
                        chrom=chrom,
                        pos=record.pos,
                        ref=ref,
                        alt=alt,
                        rsid=str(rsid_value) if rsid_value else None,
                        call_confidence=confidence,
                        dp=core_dp,
                        gq=core_gq,
                        gt=tuple(int(allele) for allele in gt),
                        phased=bool(sample.phased),
                        allele_balance=allele_balance,
                    )
                    apply_vep_annotation(row, vep_annotation)
                    if clinvar_lookup is not None:
                        clinvar_records = clinvar_lookup.lookup(chrom, record.pos, ref, alt)
                        if clinvar_records:
                            pathogenicity, evidence_rows, clinvar_rsid = aggregate_clinvar(
                                clinvar_records, variant_id, config.clinvar_snapshot
                            )
                            row["pathogenicity"].update(pathogenicity)
                            if row["rsid"] is None and clinvar_rsid is not None:
                                row["rsid"] = clinvar_rsid
                            assert clinical_evidence_writer is not None
                            assert clinvar_assertions_writer is not None
                            for evidence in evidence_rows:
                                clinvar_assertions_writer.write(evidence)
                                clinical_evidence_writer.write(_standard_clinvar_evidence(evidence))
                            counters["clinvar_matched_variants"] += 1
                            counters["clinvar_evidence_rows"] += len(evidence_rows)
                            if pathogenicity["clinvar_has_conflicts"]:
                                counters["clinvar_conflicting_variants"] += 1
                            if pathogenicity["clinvar_significance"]:
                                counters[
                                    "clinvar_significance_" + pathogenicity["clinvar_significance"]
                                ] += 1
                        else:
                            counters["clinvar_unmatched_variants"] += 1
                    if gwas_lookup is not None:
                        gwas_core, gwas_evidence, gwas_normative = gwas_lookup.lookup(
                            chrom=chrom,
                            pos=record.pos,
                            ref=ref,
                            alt=alt,
                            rsid=row["rsid"],
                            gt=tuple(int(allele) for allele in gt),
                        )
                        row["trait_associations"].update(gwas_core)
                        assert gwas_associations_writer is not None
                        assert gwas_evidence_writer is not None
                        for evidence in gwas_evidence:
                            gwas_evidence_writer.write(evidence)
                        for association in gwas_normative:
                            gwas_associations_writer.write(association)
                        if gwas_evidence:
                            counters["gwas_matched_variants"] += 1
                            counters["gwas_exact_association_rows"] += len(gwas_evidence)
                            counters["gwas_normative_association_rows"] += len(gwas_normative)
                            counters["gwas_untyped_exact_rows"] += len(gwas_evidence) - len(
                                gwas_normative
                            )
                        else:
                            counters["gwas_unmatched_variants"] += 1
                    if variant_id in actionable_variant_ids:
                        row["pathogenicity"]["is_actionable"] = True
                        counters["acmg_actionable_candidate_variants"] += 1
                    if vep_annotation is not None:
                        counters["vep_annotated_variants"] += 1
                        counters["vep_annotation_rows"] += len(
                            vep_annotation["consequence"]["all_annotations"]
                        )
                        _count_predictor_completeness(config, row, vep_annotation, counters)
                        symbols = {
                            annotation["symbol"]
                            for annotation in vep_annotation["consequence"]["all_annotations"]
                            if annotation["symbol"]
                        }
                        if row["gene"]["symbol"]:
                            symbols.add(row["gene"]["symbol"])
                        for symbol in symbols:
                            key = (symbol, chrom)
                            values = gene_stats.setdefault(
                                key,
                                {
                                    "gene_symbol": symbol,
                                    "chrom": chrom,
                                    "start_pos": record.pos,
                                    "end_pos": record.pos,
                                    "variant_count": 0,
                                    "actionable_count": 0,
                                    "pgx_count": 0,
                                },
                            )
                            values["start_pos"] = min(values["start_pos"], record.pos)
                            values["end_pos"] = max(values["end_pos"], record.pos)
                            values["variant_count"] += 1
                            values["actionable_count"] += int(
                                (variant_id, symbol) in actionable_variant_gene_pairs
                            )
                            values["pgx_count"] += int(
                                row["pharmacogenomics"]["is_pgx"]
                                and row["pharmacogenomics"]["contributes_to_gene"] == symbol
                            )
                    variant_writer.write(row)

                    source_chrom = str(_info_scalar(record, "GAB_ORIG_CONTIG"))
                    source_pos = int(_info_scalar(record, "GAB_ORIG_POS"))
                    source_ref = str(_info_scalar(record, "GAB_ORIG_REF"))
                    source_alt = str(_info_scalar(record, "GAB_ORIG_ALT"))
                    source_gt = str(_info_scalar(record, "GAB_ORIG_GT"))
                    source_ad = _info_integer_list(record, "GAB_ORIG_AD")
                    source_pl = _info_integer_list(record, "GAB_ORIG_PL")
                    source_info_json = _decode_text(_info_scalar(record, "GAB_SOURCE_INFO_B64"))
                    source_other_format_json = _decode_text(
                        _info_scalar(record, "GAB_OTHER_FORMAT_B64")
                    )
                    provenance_writer.write(
                        {
                            "variant_id": variant_id,
                            "source_id": source_id,
                            "source_record_digest": source_digest,
                            "source_ordinal": source_ordinal,
                            "source_alt_index": source_alt_index,
                            "source_chrom": source_chrom,
                            "source_pos": source_pos,
                            "source_ref": source_ref,
                            "source_alt": source_alt,
                            "source_gt": source_gt,
                            "source_dp": dp,
                            "source_gq": gq,
                            "source_ad": source_ad,
                            "source_pl": source_pl,
                            "source_info_json": source_info_json,
                            "source_other_format_json": source_other_format_json,
                            "projected_gt": [int(allele) for allele in gt],
                            "phased": bool(sample.phased),
                            "normalization_changed": (
                                source_chrom != record.contig
                                or source_pos != record.pos
                                or source_ref != ref
                                or source_alt != alt
                            ),
                            "source_filters": filters,
                        }
                    )
                    counters[f"variant_type_{row['variant_type']}"] += 1
                    counters[f"confidence_{confidence}"] += 1
        finally:
            source.close()
            fasta.close()

        gene_index_path = staging / "gene_index.parquet"
        if config.vep_enabled:
            gene_table = pa.Table.from_pylist(
                [gene_stats[key] for key in sorted(gene_stats)], schema=gene_index_schema()
            )
            pq.write_table(
                gene_table,
                gene_index_path,
                compression="zstd",
                compression_level=9,
                use_dictionary=True,
                write_statistics=True,
                version="2.6",
            )
            counters["gene_index_rows"] = len(gene_stats)

        clingen_projection = project_clingen_gene_context(
            config,
            {symbol for symbol, _chrom in gene_stats},
            staging,
        )
        if config.clingen.enabled:
            counters["clingen_gene_disease_rows"] = clingen_projection["gene_disease_rows"]
            counters["clingen_dosage_rows"] = clingen_projection["dosage_rows"]
            counters["clingen_matched_gene_symbols"] = clingen_projection["matched_gene_symbols"]
            counters["clingen_unmatched_gene_symbols"] = clingen_projection[
                "unmatched_gene_symbols"
            ]

        acmg_evidence_path = staging / "acmg_secondary_findings_evidence.parquet"
        if acmg_evaluation is not None:
            for status, count in acmg_evaluation["status_counts"].items():
                counters[f"acmg_status_{status}"] = count
            counters["acmg_decision_rows"] = len(acmg_evaluation["rows"])
            if acmg_evaluation["rows"]:
                acmg_table = pa.Table.from_pylist(
                    acmg_evaluation["rows"],
                    schema=acmg_secondary_findings_evidence_schema(),
                )
                pq.write_table(
                    acmg_table,
                    acmg_evidence_path,
                    compression="zstd",
                    compression_level=9,
                    use_dictionary=True,
                    write_statistics=True,
                    version="2.6",
                )
            findings_table = pa.Table.from_pylist(
                clinical_artifacts["clinical_findings"],
                schema=clinical_finding_schema(),
            )
            pq.write_table(
                findings_table,
                clinical_findings_path,
                compression="zstd",
                compression_level=9,
                use_dictionary=True,
                write_statistics=True,
                version="2.6",
            )
            if clinical_artifacts["interpretations"]:
                interpretations_table = pa.Table.from_pylist(
                    clinical_artifacts["interpretations"],
                    schema=interpretation_schema(),
                )
                pq.write_table(
                    interpretations_table,
                    interpretations_path,
                    compression="zstd",
                    compression_level=9,
                    use_dictionary=True,
                    write_statistics=True,
                    version="2.6",
                )
            counters["clinical_finding_rows"] = len(clinical_artifacts["clinical_findings"])
            counters["interpretation_rows"] = len(clinical_artifacts["interpretations"])

        connection = duckdb.connect()
        try:
            provenance_glob = str(provenance_path / "**" / "*.parquet").replace("'", "''")
            variants_glob = str(variants_path / "**" / "*.parquet").replace("'", "''")
            partial_glob = str(partial_path / "**" / "*.parquet").replace("'", "''")
            total_provenance = connection.execute(
                f"SELECT count(*) FROM read_parquet('{provenance_glob}')"
            ).fetchone()[0]
            partial_rows, distinct_partial = connection.execute(
                "SELECT count(*), count(DISTINCT (source_id, source_alt_index, variant_id)) "
                f"FROM read_parquet('{partial_glob}')"
            ).fetchone()
            if partial_rows != distinct_partial:
                raise ReconciliationError(
                    "duplicate normalized partial allele from the same source"
                )

            total_variants, distinct_variants = connection.execute(
                f"SELECT count(*), count(DISTINCT variant_id) FROM read_parquet('{variants_glob}')"
            ).fetchone()
            source_cursor = connection.execute(
                "SELECT DISTINCT source_id, source_alt_index "
                f"FROM read_parquet('{provenance_glob}') "
                f"UNION SELECT source_id, source_alt_index FROM read_parquet('{partial_glob}') "
                "ORDER BY source_id, source_alt_index"
            )
            actual_sources = MultisetDigest()
            while source_rows := source_cursor.fetchmany(100_000):
                for source_id, source_alt_index in source_rows:
                    actual_sources.add(f"{source_id}:a{source_alt_index}")
        finally:
            connection.close()
        actual_multiset = actual_sources.as_dict()
        if actual_multiset != expected_multiset:
            raise ReconciliationError(
                "source/output allele set mismatch after atomization: "
                f"expected {expected_multiset}, got {actual_multiset}"
            )
        if total_provenance != total_variants:
            raise ReconciliationError("every normalized variant must have one provenance row")
        if total_variants != distinct_variants:
            raise ReconciliationError("duplicate variant_id exists in variants.parquet")

        predictor_completeness = {
            plugin: {
                "eligible_variants": counters.get(f"predictor_{plugin}_eligible_variants", 0),
                "scored_variants": counters.get(f"predictor_{plugin}_scored_variants", 0),
                "missing_variants": counters.get(f"predictor_{plugin}_missing_variants", 0),
                "scored_fraction_of_eligible": (
                    counters.get(f"predictor_{plugin}_scored_variants", 0)
                    / counters[f"predictor_{plugin}_eligible_variants"]
                    if counters.get(f"predictor_{plugin}_eligible_variants", 0)
                    else None
                ),
            }
            for plugin in config.predictors.enabled_plugins
        }

        summary = {
            **config.checkpoint_identity,
            "status": "PASS",
            "counters": dict(sorted(counters.items())),
            "source_allele_multiset": actual_multiset,
            "unique_source_alleles": actual_multiset["count"],
            "observed_partial_variant_rows": partial_rows,
            "unique_normalized_variants": distinct_variants,
            "normalized_variant_rows": total_variants,
            "outputs": {
                "variants.parquet": directory_digest(variants_path),
                "variant_provenance.parquet": directory_digest(provenance_path),
                "observed_partial_variants.parquet": directory_digest(partial_path),
                **(
                    {"gene_index.parquet": sha256_file(gene_index_path)}
                    if config.vep_enabled
                    else {}
                ),
                **(
                    {"clinical_evidence.parquet": directory_digest(clinical_evidence_path)}
                    if config.clinvar_enabled
                    else {}
                ),
                **(
                    {"clinvar_assertions.parquet": directory_digest(clinvar_assertions_path)}
                    if config.clinvar_enabled
                    else {}
                ),
                **(
                    {"clinical_findings.parquet": sha256_file(clinical_findings_path)}
                    if clinical_findings_path.is_file()
                    else {}
                ),
                **(
                    {"interpretations.parquet": sha256_file(interpretations_path)}
                    if interpretations_path.is_file()
                    else {}
                ),
                **(
                    {"acmg_secondary_findings_evidence.parquet": sha256_file(acmg_evidence_path)}
                    if acmg_evidence_path.is_file()
                    else {}
                ),
                **(
                    {"gwas_associations.parquet": directory_digest(gwas_associations_path)}
                    if config.gwas.enabled
                    else {}
                ),
                **(
                    {"gwas_catalog_evidence.parquet": directory_digest(gwas_evidence_path)}
                    if config.gwas.enabled
                    else {}
                ),
                **{
                    name: metadata["sha256"]
                    for name, metadata in clingen_projection.get("outputs", {}).items()
                },
            },
            "annotation_resources": {
                "predictors": {
                    "enabled_plugins": list(config.predictors.enabled_plugins),
                    "completeness": predictor_completeness,
                    "missing_score_semantics": (
                        "unavailable evidence; never interpreted as benign or zero effect"
                    ),
                },
                "clinvar": (
                    {
                        "snapshot": config.clinvar_snapshot,
                        "normalized_vcf_sha256": clinvar_manifest["normalized_vcf_sha256"],
                    }
                    if clinvar_manifest is not None
                    else None
                ),
                "acmg_sf": (
                    {
                        "version": acmg_evaluation["policy_version"],
                        "snapshot": acmg_evaluation["policy_snapshot"],
                        "minimum_clinvar_review_stars": acmg_evaluation[
                            "minimum_clinvar_review_stars"
                        ],
                        "decision_rows": len(acmg_evaluation["rows"]),
                        "actionable_candidate_variants": len(
                            acmg_evaluation["actionable_variant_ids"]
                        ),
                        "status_counts": acmg_evaluation["status_counts"],
                    }
                    if acmg_evaluation is not None
                    else None
                ),
                "gwas_catalog": (
                    {
                        "snapshot": config.gwas.snapshot,
                        "archive_sha256": gwas_manifest["archive_sha256"],
                        "normalized_parquet_sha256": gwas_manifest["normalized_parquet_sha256"],
                        "license": gwas_manifest["license"],
                        "attribution": gwas_manifest["attribution"],
                        "counts": gwas_manifest["counts"],
                    }
                    if gwas_manifest is not None
                    else None
                ),
                "clingen": clingen_projection if config.clingen.enabled else None,
            },
        }
        atomic_write_json(staging / "variant_build_summary.json", summary)
    return summary
