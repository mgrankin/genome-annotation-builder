from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import replace
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
import pytest

import genome_annotation_builder.pgs as pgs_module
import genome_annotation_builder.pipeline as pipeline_module
import genome_annotation_builder.preflight as preflight_module
import genome_annotation_builder.variant_build as variant_build_module
from genome_annotation_builder.arrow_schemas import gwas_catalog_resource_schema
from genome_annotation_builder.clingen_resources import sync_clingen_resources
from genome_annotation_builder.clinvar_parser import SortedClinVarLookup
from genome_annotation_builder.config import load_config
from genome_annotation_builder.errors import (
    BundleValidationError,
    InputValidationError,
    ReconciliationError,
)
from genome_annotation_builder.family_evidence import _load_loco_summary
from genome_annotation_builder.normalization import normalize_projected_vcf
from genome_annotation_builder.pgs import _resolve_locus
from genome_annotation_builder.pipeline import build_bundle
from genome_annotation_builder.query import query_bundle
from genome_annotation_builder.resources import sync_pinned_resources
from genome_annotation_builder.target_callability import (
    _resolve,
    _validate_normalized_alleles,
    _write_catalog_table,
)
from genome_annotation_builder.validator import validate_bundle
from genome_annotation_builder.variant_build import build_variant_datasets
from genome_annotation_builder.vcf_ingest import (
    _CallableRegionMerger,
    canonical_core_chrom,
    ingest_reference_confidence_vcf,
)
from genome_annotation_builder.vep_parser import PIPELINE_REQUIRED_FIELDS

FIXTURES = Path(__file__).parent / "fixtures"


def _synthetic_family_config(
    tmp_path: Path,
    *,
    genotype_mismatch: bool = False,
    privacy_leak: bool = False,
    haploid_projection: bool = False,
    projection_phase_mismatch: bool = False,
) -> str:
    family_vcf = tmp_path / "synthetic-family.vcf"
    family_vcf_text = """##fileformat=VCFv4.2
##contig=<ID=1,length=200>
##INFO=<ID=DIFFICULT,Number=0,Type=Flag,Description="Synthetic difficult region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=FC,Number=1,Type=String,Description="Family confidence">
##FORMAT=<ID=FR,Number=.,Type=String,Description="Family reasons">
##FORMAT=<ID=GS,Number=1,Type=Integer,Description="Independent caller support">
##FORMAT=<ID=XS,Number=1,Type=Integer,Description="Technical replicate support">
##FORMAT=<ID=DNM,Number=1,Type=Float,Description="De-novo score">
##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFATHER\tMOTHER\tSAMPLE_ID
1\t10\t.\tC\tT\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t0/1:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t0/1:HIGH_CONFIDENCE:NO_MENDELIAN_CONFLICT:1:2:0:.
1\t11\t.\tG\tA\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t0/0:.:.:.:.:.:.\t0/1:.:.:.:.:.:.\t1/0:SUPPORTED:NO_MENDELIAN_CONFLICT:1:1:0:.
1\t11\t.\tG\tT\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t0/1:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t0/1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.
1\t22\t.\tC\tG\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t0/1:.:.:.:.:.:.\t0/1:.:.:.:.:.:.\t1/1:SUPPORTED:NO_MENDELIAN_CONFLICT:0:.:0:.
1\t23\t.\tG\tGA\t60\tPASS\tDIFFICULT\tGT:FC:FR:GS:XS:DNM:PS\t0/0:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t0/1:DE_NOVO_CANDIDATE:PARENTAL_MOSAICISM_CANDIDATE:1:0:50:.
1\t24\t.\tTA\tT\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t./.:.:.:.:.:.:.\t0/1:.:.:.:.:.:.\t0/1:REVIEW:PARENTAL_DROPOUT_SUSPECTED:.:.:0:.
1\t25\t.\tA\tGG\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t0/1:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t0/1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.
1\t26\t.\tC\tT\t20\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t0/0:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t1/.:REVIEW:UNCALLABLE_MEMBER:.:.:.:.
1\t30\t.\tC\tA\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t0/1:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t1:SUPPORTED:CHILD_DROPOUT_SUSPECTED:1:.:0:.
"""
    if genotype_mismatch:
        family_vcf_text = family_vcf_text.replace(
            "0/1:HIGH_CONFIDENCE:NO_MENDELIAN_CONFLICT:1:2:0:.",
            "1/1:HIGH_CONFIDENCE:NO_MENDELIAN_CONFLICT:1:2:0:.",
            1,
        )
    if haploid_projection:
        diploid = (
            "1\t25\t.\tA\tGG\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t"
            "0/1:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t"
            "0/1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.\n"
        )
        haploid = diploid.replace(
            "\t0/1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.\n",
            "\t1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.\n",
        )
        if family_vcf_text.count(diploid) != 1:
            raise AssertionError("synthetic haploid projection boundary changed")
        family_vcf_text = family_vcf_text.replace(diploid, haploid, 1)
    if projection_phase_mismatch:
        if haploid_projection:
            raise AssertionError("synthetic projection modes are mutually exclusive")
        diploid = (
            "1\t25\t.\tA\tGG\t60\tPASS\t.\tGT:FC:FR:GS:XS:DNM:PS\t"
            "0/1:.:.:.:.:.:.\t0/0:.:.:.:.:.:.\t"
            "0/1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.\n"
        )
        phased = diploid.replace(
            "\t0/1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.\n",
            "\t0|1:SUPPORTED:NO_MENDELIAN_CONFLICT:1:.:0:.\n",
        )
        if family_vcf_text.count(diploid) != 1:
            raise AssertionError("synthetic diploid projection boundary changed")
        family_vcf_text = family_vcf_text.replace(diploid, phased, 1)
    if privacy_leak:
        family_vcf_text = family_vcf_text.replace("NO_MENDELIAN_CONFLICT", "FATHER_SUPPORT", 1)
    family_vcf.write_text(family_vcf_text, encoding="utf-8")
    joint_vcf = Path(pysam.tabix_index(str(family_vcf), preset="vcf", force=True))
    pedigree = tmp_path / "synthetic-family.ped"
    pedigree.write_text(
        """FAM001 FATHER 0 0 1 0
FAM001 MOTHER 0 0 2 0
FAM001 SAMPLE_ID FATHER MOTHER 0 0
""",
        encoding="utf-8",
    )
    review = tmp_path / "synthetic-loco-review.tsv"
    review.write_text(
        """chrom\tpos\tref\talt\tstatus\tchild_gt\tfather_gt\tmother_gt\texact_allele_in_loco_model
1\t25\tA\tGG\tuninformative_parent_not_callable\t1/1\t0/1\t0/0\tTrue
1\t23\tG\tGA\tde_novo_candidate\t0/1\t0/0\t0/0\tFalse
1\t24\tTA\tT\tuninformative_parent_not_callable\t0/1\t./.\t0/1\tTrue
""",
        encoding="utf-8",
    )
    summary = tmp_path / "synthetic-loco-summary.json"
    summary.write_text(
        json.dumps(
            {
                "heldout_child": "SAMPLE_ID",
                "validated_unit": "SAMPLE_ID",
                "counters": {
                    "heldout_carried_alleles": 8,
                    "mendelian_consistent": 5,
                    "de_novo_candidate": 1,
                    "mendelian_conflict": 0,
                    "uninformative_parent_not_callable": 2,
                },
                "review_rows": 3,
                "informative_alleles": 6,
                "mendelian_consistency_fraction": 5 / 6,
                "parental_haplotype_crosscheck": {"model_is_child_independent": True},
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    return f"""
  family_evidence:
    enabled: true
    joint_vcf: {joint_vcf}
    joint_vcf_expected_sha256: {digest(joint_vcf)}
    pedigree: {pedigree}
    pedigree_expected_sha256: {digest(pedigree)}
    source_version: synthetic-family-v1
    snapshot_date: "2026-08-15"
    tags:
      confidence_format: FC
      reasons_format: FR
      caller_concordance_format: GS
      technical_replicate_concordance_format: XS
      de_novo_score_format: DNM
      phase_set_format: PS
      difficult_info_flags: [DIFFICULT]
      de_novo_confidence_values: [DE_NOVO_CANDIDATE]
      conflict_confidence_values: [MENDELIAN_CONFLICT]
      parental_mosaicism_reason_tokens: [PARENTAL_MOSAICISM_CANDIDATE]
      parental_dropout_reason_tokens: [PARENTAL_DROPOUT_SUSPECTED]
      child_dropout_reason_tokens: [CHILD_DROPOUT_SUSPECTED]
    loco:
      review_tsv: {review}
      review_expected_sha256: {digest(review)}
      summary_json: {summary}
      summary_expected_sha256: {digest(summary)}
"""


def _environment(
    tmp_path: Path,
    *,
    target_catalog: bool = False,
    pgs: bool = False,
    pgs_calibration: bool = True,
    family: bool = False,
    family_genotype_mismatch: bool = False,
    family_privacy_leak: bool = False,
    family_haploid_projection: bool = False,
    family_projection_phase_mismatch: bool = False,
    gwas: bool = False,
    clingen: bool = False,
):
    vep_cache_version = '"116"' if clingen else "null"
    reference = tmp_path / "synthetic.fa"
    shutil.copyfile(FIXTURES / "synthetic-grch38.fa", reference)
    shutil.copyfile(FIXTURES / "synthetic-grch38.fa.fai", Path(f"{reference}.fai"))
    plain_vcf = tmp_path / "synthetic.vcf"
    shutil.copyfile(FIXTURES / "synthetic-reference-confidence.vcf", plain_vcf)
    if family:
        text = plain_vcf.read_text(encoding="utf-8")
        marker = "1\t30\t.\tC\tA\t50\tPASS"
        additions = (
            "1\t25\t.\tA\tGG\t50\tPASS\t.\tGT:DP:GQ:AD:FT\t"
            f"{'1' if family_haploid_projection else '0/1'}:24:50:12,12:PASS\n"
            "1\t26\t.\tC\tT\t20\tPASS\t.\tGT:DP:GQ:AD:FT\t"
            "1/.:9:15:0,5:PASS\n"
        )
        if text.count(marker) != 1:
            raise AssertionError("synthetic family fixture insertion boundary changed")
        plain_vcf.write_text(text.replace(marker, additions + marker), encoding="utf-8")
    compressed_vcf = Path(pysam.tabix_index(str(plain_vcf), preset="vcf", force=True))
    callability_config = ""
    if target_catalog:
        targets = tmp_path / "synthetic-targets.tsv"
        shutil.copyfile(FIXTURES / "synthetic-targets.tsv", targets)
        targets_sha256 = hashlib.sha256(targets.read_bytes()).hexdigest()
        callability_config = f"""
callability:
  target_catalog:
    enabled: true
    path: {targets}
    expected_sha256: {targets_sha256}
    assay_scope: synthetic_reference_confidence_vcf
"""
    pgs_config = ""
    if pgs:
        scoring_file = tmp_path / "synthetic-pgs.txt"
        shutil.copyfile(FIXTURES / "synthetic-pgs.txt", scoring_file)
        scoring_sha256 = hashlib.sha256(scoring_file.read_bytes()).hexdigest()
        calibration_config = ""
        if pgs_calibration:
            calibration_file = tmp_path / "synthetic-pgs-calibration.tsv"
            calibration_file.write_text(
                f"""#pgs_id=PGS000001
#scoring_file_sha256={scoring_sha256}
#reference_population=Synthetic reference cohort
#calibration_source=Synthetic empirical CDF
#calibration_source_version=1.0
#calibration_date=2026-08-15
#max_missing_weight_fraction=0.1
score\tpercentile
0\t10
2\t50
3\t90
""",
                encoding="utf-8",
            )
            calibration_sha256 = hashlib.sha256(calibration_file.read_bytes()).hexdigest()
            calibration_config = f"""
        calibration:
          file: {calibration_file}
          expected_sha256: {calibration_sha256}"""
        pgs_config = f"""
  pgs:
    enabled: true
    models:
      - pgs_id: PGS000001
        scoring_file: {scoring_file}
        expected_sha256: {scoring_sha256}
        source_version: "2026-08-15"
        training_date: "2015-02-18"
        trait: Synthetic trait
        model_type: linear
        license: CC0-1.0
        redistribution_allowed: true
        max_missing_weight_fraction: 0.1{calibration_config}
"""
    family_config = (
        _synthetic_family_config(
            tmp_path,
            genotype_mismatch=family_genotype_mismatch,
            privacy_leak=family_privacy_leak,
            haploid_projection=family_haploid_projection,
            projection_phase_mismatch=family_projection_phase_mismatch,
        )
        if family
        else ""
    )
    gwas_config = ""
    if gwas:
        gwas_resource = tmp_path / "synthetic-gwas.parquet"
        pq.write_table(
            pa.Table.from_pylist(
                [
                    {
                        "chrom": "chr1",
                        "pos": 10,
                        "source_row_number": 2,
                        "risk_rsid": "rs100001",
                        "rsid_aliases": ["rs100001"],
                        "effect_allele": "T",
                        "reported_trait": "Synthetic GWAS trait",
                        "mapped_traits": ["synthetic trait measurement"],
                        "mapped_trait_uris": ["http://example.org/EFO_1"],
                        "p_value": 1e-9,
                        "effect_size": 0.5,
                        "effect_type": "beta",
                        "effect_type_basis": "synthetic explicit increase marker",
                        "ci_lower": 0.3,
                        "ci_upper": 0.7,
                        "sample_size": 1000,
                        "pmid": "12345678",
                        "study_accession": "GCST000001",
                        "raw_effect_text": "0.5",
                        "raw_ci_text": "[0.3-0.7] unit increase",
                        "source_record_digest": "a" * 64,
                    }
                ],
                schema=gwas_catalog_resource_schema(),
            ),
            gwas_resource,
            compression="zstd",
            version="2.6",
        )
        gwas_config = f"""
  gwas:
    enabled: true
    snapshot: "2026-08-03"
    archive: {tmp_path / "synthetic-gwas.zip"}
    archive_expected_sha256: 3ed9a1335926181d36334fe66fbf6754faf5120c780bc385241c5ad3fafcfe47
    normalized_parquet: {gwas_resource}
    manifest: {tmp_path / "synthetic-gwas-manifest.json"}
    license: EMBL-EBI Terms of Use
    redistribution_allowed: false
"""
    clingen_config = ""
    if clingen:
        clingen_config = """
  clingen:
    enabled: true
    snapshot: "2026-08-15"
    license: CC0-1.0
    redistribution_allowed: true
"""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        f"""
project:
  output_dir: {tmp_path / "output"}
  work_dir: {tmp_path / "work"}
  bundle_name: synthetic.genome
input:
  vcf: {compressed_vcf}
  sample: SAMPLE_ID
  genome_build: GRCh38
reference:
  fasta: {reference}
  expected_sha256: null
  contig_style: bare
resources:
  cache_dir: {tmp_path / "resources"}
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
normalization:
  bcftools: bcftools
  check_ref: error
quality:
  callable_min_dp: 10
  callable_min_gq: 20
  medium_min_dp: 10
  medium_min_gq: 20
  high_min_dp: 20
  high_min_gq: 30
annotation:
  vep:
    enabled: {str(clingen).lower()}
    executable: vep
    cache_dir: {tmp_path / "resources" / "vep"}
    cache_version: {vep_cache_version}
{clingen_config}
{pgs_config}
{family_config}
{gwas_config}
{callability_config}
output:
  agent_exploration_profile: true
  reference_confidence_extension: true
  partition_variants_by_chromosome: true
  overwrite: false
privacy:
  include_input_paths: false
  include_sample_identifier: false
""",
        encoding="utf-8",
    )
    return load_config(config_path), compressed_vcf


def test_canonical_core_chrom() -> None:
    assert canonical_core_chrom("1") == "chr1"
    assert canonical_core_chrom("chrX") == "chrX"
    assert canonical_core_chrom("MT") == "chrM"
    assert canonical_core_chrom("chrMT") == "chrM"
    assert canonical_core_chrom("GL000220.1") is None


def test_callable_region_sweep_is_disjoint_and_deduplicated() -> None:
    class RowCollector:
        def __init__(self) -> None:
            self.rows: list[dict[str, object]] = []

        def write(self, row: dict[str, object]) -> None:
            self.rows.append(row)

    writer = RowCollector()
    merger = _CallableRegionMerger(writer)  # type: ignore[arg-type]
    merger.add(chrom="1", start_pos=1, end_pos=10, callable_value=True, minimum_depth=10)
    merger.add(chrom="1", start_pos=5, end_pos=5, callable_value=False, minimum_depth=20)
    merger.add(chrom="1", start_pos=5, end_pos=5, callable_value=False, minimum_depth=20)
    merger.add(chrom="1", start_pos=11, end_pos=12, callable_value=True, minimum_depth=10)
    merger.add(chrom="2", start_pos=2, end_pos=3, callable_value=True, minimum_depth=10)
    merger.flush()

    assert [
        (row["chrom"], row["start_pos"], row["end_pos"], row["callable"], row["minimum_depth"])
        for row in writer.rows
    ] == [
        ("chr1", 1, 4, True, 10),
        ("chr1", 5, 5, False, 20),
        ("chr1", 6, 12, True, 10),
        ("chr2", 2, 3, True, 10),
    ]
    assert len({row["region_id"] for row in writer.rows}) == len(writer.rows)


def test_ingest_normalize_and_reconcile(tmp_path: Path) -> None:
    config, input_vcf = _environment(tmp_path)

    digest = hashlib.sha256(input_vcf.read_bytes()).hexdigest()
    ingest_dir = tmp_path / "ingest"
    ingest = ingest_reference_confidence_vcf(config, digest, ingest_dir)
    assert ingest["counters"]["records_total"] == 11
    assert ingest["counters"]["called_concrete_alleles_projected"] == 7
    assert ingest["counters"]["called_alleles_symbolic_or_non_acgt_allele"] == 1
    assert ingest["counters"]["records_no_call"] == 1
    assert ingest["counters"]["records_reference_block"] == 2
    assert ingest["counters"]["records_explicit_reference_site"] == 1

    normalized_dir = tmp_path / "normalized"
    normalization = normalize_projected_vcf(config, ingest_dir / "projected.vcf.gz", normalized_dir)
    assert normalization["record_count"] == 7

    variant_dir = tmp_path / "variant-build"
    built = build_variant_datasets(
        config,
        normalized_dir / "normalized.vcf.gz",
        digest,
        ingest["accepted_evidence_multiset"],
        variant_dir,
    )
    assert built["unique_normalized_variants"] == 7
    assert built["unique_source_alleles"] == 7

    connection = duckdb.connect()
    try:
        variants = connection.execute(
            "SELECT variant_id, genotype.gt, genotype.zygosity, quality.call_confidence "
            f"FROM read_parquet('{variant_dir / 'variants.parquet'}/**/*.parquet') "
            "ORDER BY chrom, pos, alt"
        ).fetchall()
        evidence = connection.execute(
            "SELECT state, count(*) FROM "
            f"read_parquet('{ingest_dir / 'reference_evidence.parquet'}/**/*.parquet') "
            "GROUP BY state ORDER BY state"
        ).fetchall()
    finally:
        connection.close()
    assert len(variants) == 7
    assert ("chr1:11:G:A", [1, 0], "het", "high") in variants
    assert ("chr1:11:G:T", [0, 1], "het", "high") in variants
    assert ("chr1:22:C:G", [1, 1], "hom_alt", "low") in variants
    assert ("chr1:30:C:A", [1], "hom_alt", "medium") in variants
    assert evidence == [
        ("explicit_reference_site", 1),
        ("no_call", 1),
        ("reference_block", 2),
    ]


def test_normalization_reorders_lexical_reference_header_for_streaming_joins(
    tmp_path: Path,
) -> None:
    config, _input_vcf = _environment(tmp_path)
    reference = tmp_path / "lexical-reference.fa"
    reference.write_text(
        ">1\n" + "A" * 100 + "\n>10\n" + "C" * 100 + "\n>2\n" + "G" * 100 + "\n",
        encoding="ascii",
    )
    pysam.faidx(str(reference))
    source = tmp_path / "lexical-source.vcf"
    source.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=1,length=100>
##contig=<ID=10,length=100>
##contig=<ID=2,length=100>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID
1\t10\t.\tA\tC\t60\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:50:15,15
10\t10\t.\tC\tA\t60\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:50:15,15
2\t10\t.\tG\tA\t60\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:50:15,15
""",
        encoding="ascii",
    )
    compressed = Path(pysam.tabix_index(str(source), preset="vcf", force=True))
    config = replace(config, input_vcf=compressed, reference_fasta=reference)
    digest = hashlib.sha256(compressed.read_bytes()).hexdigest()

    ingest_dir = tmp_path / "lexical-ingest"
    ingest_reference_confidence_vcf(config, digest, ingest_dir)
    projected = pysam.VariantFile(str(ingest_dir / "projected.vcf.gz"))
    try:
        assert list(projected.header.contigs) == ["1", "10", "2"]
    finally:
        projected.close()

    normalized_dir = tmp_path / "lexical-normalized"
    summary = normalize_projected_vcf(config, ingest_dir / "projected.vcf.gz", normalized_dir)
    assert summary["record_count"] == 3
    assert summary["canonical_order"] == {
        "status": "PASS",
        "canonical_core_header_order": ["chr1", "chr2", "chr10"],
        "header_contigs": 3,
        "non_core_header_contigs": 0,
        "records": 3,
        "first_locus": "chr1:10",
        "last_locus": "chr10:10",
        "duplicate_query_alleles": 0,
    }
    normalized_path = normalized_dir / "normalized.vcf.gz"
    normalized = pysam.VariantFile(str(normalized_path))
    try:
        assert list(normalized.header.contigs) == ["1", "2", "10"]
        assert [(record.contig, record.pos) for record in normalized] == [
            ("1", 10),
            ("2", 10),
            ("10", 10),
        ]
    finally:
        normalized.close()

    clinvar = tmp_path / "natural-clinvar.vcf"
    clinvar.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=1,length=100>
##contig=<ID=2,length=100>
##contig=<ID=10,length=100>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10\t1\tA\tC\t.\tPASS\t.
2\t10\t2\tG\tA\t.\tPASS\t.
10\t10\t3\tC\tA\t.\tPASS\t.
""",
        encoding="ascii",
    )
    normalized = pysam.VariantFile(str(normalized_path))
    try:
        with SortedClinVarLookup(clinvar, "20260808") as lookup:
            matched = [
                len(
                    lookup.lookup(
                        canonical_core_chrom(record.contig) or "",
                        record.pos,
                        record.ref,
                        record.alts[0],
                    )
                )
                for record in normalized
            ]
    finally:
        normalized.close()
    assert matched == [1, 1, 1]


def test_complex_allele_atomization_restores_biallelic_gt_and_phase(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path)
    source = tmp_path / "complex-allele.vcf"
    source.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=1,length=200>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID
1\t27\t.\tG\tATAA\t50\tPASS\t.\tGT:DP:GQ:AD\t1|0:22:45:12,10
""",
        encoding="utf-8",
    )
    compressed = Path(pysam.tabix_index(str(source), preset="vcf", force=True))
    config = replace(config, input_vcf=compressed)
    digest = hashlib.sha256(compressed.read_bytes()).hexdigest()

    ingest_dir = tmp_path / "complex-ingest"
    ingest = ingest_reference_confidence_vcf(config, digest, ingest_dir)
    assert ingest["counters"]["called_concrete_alleles_projected"] == 1

    normalized_dir = tmp_path / "complex-normalized"
    normalization = normalize_projected_vcf(config, ingest_dir / "projected.vcf.gz", normalized_dir)
    assert normalization["record_count"] == 2
    assert normalization["projection"] == {
        "atomized_records": 2,
        "biallelic_records": 2,
        "spanning_star_companion_records_removed": 2,
    }

    normalized = pysam.VariantFile(str(normalized_dir / "normalized.vcf.gz"))
    try:
        rows = [
            (
                record.ref,
                record.alts,
                record.samples["SAMPLE"]["GT"],
                record.samples["SAMPLE"].phased,
            )
            for record in normalized
        ]
    finally:
        normalized.close()
    assert rows == [
        ("G", ("A",), (1, 0), True),
        ("G", ("GTAA",), (1, 0), True),
    ]

    variant_dir = tmp_path / "complex-variant-build"
    built = build_variant_datasets(
        config,
        normalized_dir / "normalized.vcf.gz",
        digest,
        ingest["accepted_evidence_multiset"],
        variant_dir,
    )
    assert built["unique_source_alleles"] == 1
    assert built["unique_normalized_variants"] == 2
    assert built["normalized_variant_rows"] == 2
    connection = duckdb.connect()
    try:
        variants = connection.execute(
            "SELECT variant_id, genotype.gt, genotype.phased "
            f"FROM read_parquet('{variant_dir / 'variants.parquet'}/**/*.parquet') "
            "ORDER BY variant_id"
        ).fetchall()
    finally:
        connection.close()
    assert variants == [
        ("chr1:27:G:A", [1, 0], True),
        ("chr1:27:G:GTAA", [1, 0], True),
    ]


def test_full_release_a_bundle(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path)
    resources = sync_pinned_resources(config)
    assert len(resources["files"]) == 6
    result = build_bundle(config)
    assert result["status"] == "PASS"
    assert result["validation"]["valid"] is True
    bundle = config.bundle_path

    overview = query_bundle(bundle, "overview", limit=10)
    assert overview["rows"][0]["variant_count"] == 7
    assert query_bundle(bundle, "callability", term="chr1:10:C:T")["conclusion"] == (
        "non_reference_observed"
    )
    assert query_bundle(bundle, "callability", term="chr1:2")["conclusion"] == (
        "reference_observed"
    )
    assert query_bundle(bundle, "callability", term="chr1:13")["conclusion"] == ("unknown_no_call")
    assert query_bundle(bundle, "callability", term="chr1:16")["conclusion"] == (
        "low_confidence_reference_observation"
    )
    assert query_bundle(bundle, "callability", term="chr1:100")["conclusion"] == "unknown"

    text_payload = "\n".join(
        path.read_text(encoding="utf-8")
        for path in bundle.rglob("*")
        if path.is_file() and path.suffix in {".json", ".md", ".py"}
    )
    assert "SAMPLE_ID" not in text_payload
    assert str(tmp_path) not in text_payload


def test_gwas_catalog_exact_projection_is_audited_and_queryable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config, _input_vcf = _environment(tmp_path, gwas=True)
    normalized_sha256 = hashlib.sha256(config.gwas.normalized_parquet.read_bytes()).hexdigest()
    fake_manifest = {
        "status": "verified",
        "snapshot": "2026-08-03",
        "source_url": "https://example.invalid/synthetic-gwas",
        "archive_sha256": config.gwas.archive_expected_sha256,
        "archive_member_crc32": "17484714",
        "normalized_parquet_sha256": normalized_sha256,
        "license": "EMBL-EBI Terms of Use",
        "attribution": "NHGRI-EBI GWAS Catalog; synthetic test row",
        "matching_policy": (
            "exact GRCh38 coordinate, rsID alias, and effect allele equal to REF or ALT"
        ),
        "effect_type_policy": ("emit a normative row only when effect type and value are explicit"),
        "counts": {
            "source_rows": 1,
            "normalized_rows": 1,
            "typed_effect_rows": 1,
            "untyped_effect_rows": 0,
        },
    }

    def verify_synthetic_gwas(_config: object) -> dict[str, object]:
        return fake_manifest

    monkeypatch.setattr(preflight_module, "verify_gwas_catalog", verify_synthetic_gwas)
    monkeypatch.setattr(variant_build_module, "verify_gwas_catalog", verify_synthetic_gwas)

    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["status"] == "PASS"
    assert result["validation"]["valid"] is True
    gwas_metrics = result["validation"]["metrics"]["gwas_catalog"]
    assert gwas_metrics == {
        "status": "validated",
        "association_rows": 1,
        "evidence_rows": 1,
        "exact_join_failures": 0,
        "dosage_failures": 0,
        "variant_summary_failures": 0,
    }

    bundle = config.bundle_path
    connection = duckdb.connect()
    try:
        compact = connection.execute(
            "SELECT trait_associations.is_gwas_hit, trait_associations.effect_size, "
            "trait_associations.effect_type, trait_associations.effect_allele, "
            "trait_associations.p_value FROM "
            f"read_parquet('{bundle / 'variants.parquet'}/**/*.parquet') "
            "WHERE variant_id='chr1:10:C:T'"
        ).fetchone()
        evidence = connection.execute(
            "SELECT effect_allele, effect_allele_dosage, normative_row_emitted, "
            "raw_effect_text, raw_ci_text FROM "
            f"read_parquet('{bundle / 'extensions' / 'gwas_catalog_evidence.parquet'}"
            "/**/*.parquet')"
        ).fetchone()
    finally:
        connection.close()
    assert compact == (True, 0.5, "beta", "T", pytest.approx(1e-9))
    assert evidence == ("T", 1, True, "0.5", "[0.3-0.7] unit increase")

    query = query_bundle(bundle, "gwas", term="rs100001")
    assert query["count"] == 1
    assert query["rows"][0]["variant_id"] == "chr1:10:C:T"
    assert query["rows"][0]["effect_allele_dosage"] == 1
    assert query["typed_associations"][0]["effect_type"] == "beta"
    assert "no row is a personal diagnosis" in query["interpretation"]

    audit = query_bundle(bundle, "evidence", term="chr1:10:C:T", limit=20)
    assert "extensions/gwas_catalog_evidence" in audit["tables_searched"]
    assert any(row.get("_table") == "extensions/gwas_catalog_evidence" for row in audit["rows"])

    raw_artifacts = json.loads((bundle / "raw_artifacts.json").read_text(encoding="utf-8"))
    gwas_artifacts = [
        row for row in raw_artifacts["artifacts"] if "GWAS Catalog" in row["artifact_type"]
    ]
    assert len(gwas_artifacts) == 2
    assert all(row["disposition"] == "referenced" for row in gwas_artifacts)
    assert all("integrity_identifier" in row for row in gwas_artifacts)

    transformations = json.loads((bundle / "transformations.json").read_text(encoding="utf-8"))
    gwas_stage = next(
        row
        for row in transformations["stages"]
        if row["stage_id"] == "gwas_catalog_exact_association_projection"
    )
    assert gwas_stage["exact_association_rows"] == 1
    assert gwas_stage["normative_association_rows"] == 1
    assert gwas_stage["clinical_grade_claims"] == 0

    omissions = json.loads((bundle / "omissions.json").read_text(encoding="utf-8"))
    assert any(
        row["scope"] == "personal causal or clinical interpretation of GWAS associations"
        for row in omissions["omissions"]
    )


def test_clingen_gene_context_is_projected_queryable_and_semantically_gated(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config, _input_vcf = _environment(tmp_path, clingen=True)
    fake_cache = {
        "cache_version": "116",
        "tree_sha256": "a" * 64,
        "archive_sha256": "b" * 64,
    }
    monkeypatch.setattr(preflight_module, "verify_vep_cache", lambda _config: fake_cache)
    monkeypatch.setattr(
        preflight_module,
        "verify_vep_runtime",
        lambda _config: ("116.0", ["synthetic-vep", "--help"]),
    )

    def synthetic_vep(_config: object, input_vcf: Path, output_dir: Path) -> dict[str, object]:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_vcf = output_dir / "annotated.vcf.gz"
        fields = tuple(sorted(PIPELINE_REQUIRED_FIELDS))
        with pysam.VariantFile(str(input_vcf)) as source:
            header = source.header.copy()
            header.info.add(
                "CSQ",
                ".",
                "String",
                "Synthetic VEP annotations. Format: " + "|".join(fields),
            )
            with pysam.VariantFile(str(output_vcf), "wz", header=header) as destination:
                record_count = 0
                for record in source:
                    record.translate(header)
                    values = {field: "" for field in fields}
                    values.update(
                        {
                            "Allele": record.alts[0],
                            "Consequence": "missense_variant",
                            "IMPACT": "MODERATE",
                            "SYMBOL": "BRCA1",
                            "Gene": "ENSG00000012048",
                            "Feature_type": "Transcript",
                            "Feature": "ENST00000357654",
                            "BIOTYPE": "protein_coding",
                            "HGVSc": "ENST00000357654:c.1A>G",
                            "HGVSp": "ENSP00000350283:p.Met1Val",
                            "CANONICAL": "YES",
                            "MANE_SELECT": "NM_007294.4",
                            "APPRIS": "principal1",
                            "TSL": "1",
                            "HGNC_ID": "1100",
                            "SIFT": "tolerated(0.2)",
                            "PolyPhen": "benign(0.1)",
                            "ALLELE_NUM": "1",
                            "gnomADg_AF": "0.0001",
                            "gnomADe_AF": "0.0001",
                        }
                    )
                    record.info["CSQ"] = ("|".join(values[field] for field in fields),)
                    destination.write(record)
                    record_count += 1
        pysam.tabix_index(str(output_vcf), preset="vcf", force=True)
        index = Path(str(output_vcf) + ".tbi")
        summary: dict[str, object] = {
            "status": "PASS",
            "record_count": record_count,
            "vep_version": "116.0",
            "cache_version": "116",
            "cache_tree_sha256": fake_cache["tree_sha256"],
            "cache_metadata": fake_cache,
            "container_image": None,
            "completed_at": "2026-08-15T00:00:00Z",
            "command_semantics": {"runtime": "synthetic-test-double", "network": "none"},
            "outputs": {
                output_vcf.name: hashlib.sha256(output_vcf.read_bytes()).hexdigest(),
                index.name: hashlib.sha256(index.read_bytes()).hexdigest(),
            },
        }
        (output_dir / "vep_summary.json").write_text(
            json.dumps(summary, sort_keys=True) + "\n", encoding="utf-8"
        )
        return summary

    monkeypatch.setattr(pipeline_module, "run_vep_annotation", synthetic_vep)
    sync_pinned_resources(config)
    sync_clingen_resources(config)
    result = build_bundle(config)
    assert result["status"] == "PASS"
    metrics = result["validation"]["metrics"]["clingen_gene_context"]
    assert metrics["status"] == "validated"
    assert metrics["gene_disease_rows"] > 0
    assert metrics["dosage_rows"] == 1
    assert metrics["matched_gene_symbols"] == 1
    assert metrics["variant_pathogenicity_claims"] == 0
    assert metrics["diagnosis_claims"] == 0

    gene_query = query_bundle(config.bundle_path, "gene", term="BRCA1", limit=20)
    context = gene_query["clingen_gene_context"]
    assert {row["_table"] for row in context} == {
        "extensions/clingen_gene_disease",
        "extensions/clingen_dosage_sensitivity",
    }
    assert "do not diagnose" in gene_query["clingen_semantic_boundary"]
    clinical_query = query_bundle(config.bundle_path, "clinical", term="BRCA1", limit=20)
    assert any(row["_table"].startswith("extensions/clingen_") for row in clinical_query["rows"])
    assert "gene-level context only" in clinical_query["interpretation"]

    manifest_path = config.bundle_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["resource_provenance"]["clingen_gene_level"]["semantic_boundary"] = (
        "a gene match proves a diagnosis"
    )
    manifest_path.write_text(json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(
        BundleValidationError,
        match="ClinGen provenance boundary mismatch: semantic_boundary",
    ):
        validate_bundle(config.bundle_path, reference_fasta=config.reference_fasta)


def test_family_evidence_is_subject_only_and_exactly_reconciled(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path, family=True)
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["status"] == "PASS"
    family_metrics = result["validation"]["metrics"]["family_evidence"]
    assert family_metrics["rows"] == 9
    assert family_metrics["orphan_family_rows"] == 0
    assert family_metrics["variants_without_family_rows"] == 0
    assert family_metrics["subject_ledger_mismatches"] == 0

    bundle = config.bundle_path
    family_path = bundle / "extensions" / "family_evidence.parquet"
    connection = duckdb.connect()
    try:
        rows = {
            row[0]: row[1:]
            for row in connection.execute(
                "SELECT variant_id, inheritance_state, parental_origin, mendelian_status, "
                "loco_status, de_novo_candidate, parental_mosaicism_candidate, "
                "parental_dropout_suspected, child_dropout_suspected, caller_concordance, "
                "technical_replicate_concordance, complex_region "
                f"FROM read_parquet('{family_path}/**/*.parquet') ORDER BY variant_id"
            ).fetchall()
        }
        columns = connection.execute(
            f"DESCRIBE SELECT * FROM read_parquet('{family_path}/**/*.parquet')"
        ).fetchall()
        projected_loco = connection.execute(
            "SELECT loco_source_status, loco_subject_gt_concordant FROM "
            f"read_parquet('{family_path}/**/*.parquet') WHERE pos=25 ORDER BY variant_id"
        ).fetchall()
    finally:
        connection.close()
    assert rows["chr1:10:C:T"] == (
        "inherited_paternal",
        "paternal_observed",
        "consistent",
        "not_individually_persisted",
        False,
        False,
        False,
        False,
        "exact_genotype_concordant",
        "both_sources_exact_genotype_concordant",
        False,
    )
    assert rows["chr1:23:G:GA"] == (
        "de_novo_candidate",
        "neither_parent_observed",
        "conflict",
        "de_novo_candidate",
        True,
        True,
        False,
        False,
        "exact_genotype_concordant",
        "discordant",
        True,
    )
    assert rows["chr1:24:TA:T"][0:4] == (
        "unresolved_parent_no_call",
        "unknown",
        "parent_or_subject_no_call",
        "uninformative_parent_not_callable",
    )
    for variant_id in ("chr1:25:A:G", "chr1:25:A:AG"):
        assert rows[variant_id][0:4] == (
            "inherited_paternal",
            "paternal_observed",
            "consistent",
            "source_genotype_discordant",
        )
    assert projected_loco == [
        ("uninformative_parent_not_callable", False),
        ("uninformative_parent_not_callable", False),
    ]
    assert rows["chr1:30:C:A"][0:4] == (
        "not_assessed",
        "unknown",
        "not_assessed",
        "not_individually_persisted",
    )
    column_names = {row[0] for row in columns}
    assert {"father_gt", "mother_gt", "relative_identifier"}.isdisjoint(column_names)

    evidence = query_bundle(bundle, "evidence", term="chr1:23:G:GA", limit=20)
    assert "extensions/family_evidence" in evidence["tables_searched"]
    family_rows = [
        row for row in evidence["rows"] if row.get("_table") == "extensions/family_evidence"
    ]
    assert len(family_rows) == 1
    assert family_rows[0]["de_novo_candidate"] is True

    text_payload = "\n".join(
        path.read_text(encoding="utf-8")
        for path in bundle.rglob("*")
        if path.is_file() and path.suffix in {".json", ".md", ".py"}
    )
    for private_identifier in ("FATHER", "MOTHER", "SAMPLE_ID", "FAM001"):
        assert private_identifier not in text_payload

    family_summary = json.loads(
        (
            config.work_dir
            / f"run-{result['run_key']}"
            / "05-family-evidence"
            / "family_evidence_summary.json"
        ).read_text(encoding="utf-8")
    )
    projection = family_summary["normalization_projection"]
    assert projection["source_alleles"] == 1
    assert projection["output_alleles"] == 2
    assert projection["extra_output_alleles"] == 1
    assert projection["family_source_records_projected"] == 1
    assert projection["incomplete_subject_genotype_records_omitted"] == 1
    assert projection["incomplete_subject_genotype_records_with_alt_evidence"] == 1
    assert family_summary["loco"]["materialization"]["review_rows"] == 3
    assert family_summary["loco"]["materialization"]["projected_review_rows"] == 4
    assert family_summary["loco"]["materialization"]["normalization_projected_source_rows"] == 1
    assert (
        family_summary["loco"]["materialization"]["normalization_projection_source_gt_discordant"]
        == 1
    )
    assert family_summary["loco"]["materialization"]["matched_genotype_discordant"] == 2
    assert family_summary["loco"]["materialization"]["source_order_inversions_normalized"] == 1


def test_family_loco_unscorable_source_projection_and_sparse_counters(tmp_path: Path) -> None:
    config, _ = _environment(tmp_path, family=True)
    family = config.family_evidence
    loco = family.loco
    assert loco is not None
    # Two carried ALTs from 1|2 and one from 1/. have no complete target
    # genotype in the held-out source. The separate family master has its own calls.
    review = (
        "chrom\tpos\tref\talt\tstatus\tchild_gt\texact_allele_in_loco_model"
        "\tchild_source_gt\tchild_projection_status\n"
        "1\t11\tG\tT\tuninformative_child_genotype\t.\tTrue\t1|2\tother_known_alt\n"
        "1\t11\tG\tA\tuninformative_child_genotype\t.\tFalse\t1|2\tother_known_alt\n"
        "1\t26\tC\tT\tuninformative_child_genotype\t.\tFalse\t1/.\tmissing_gt\n"
    )
    loco.review_tsv.write_text(review, encoding="utf-8")
    review_digest = hashlib.sha256(loco.review_tsv.read_bytes()).hexdigest()
    summary = json.loads(loco.summary_json.read_bytes())
    summary.update(
        counters={
            "source_records": 2,
            "heldout_carried_alleles": 3,
            "uninformative_child_genotype": 3,
        },
        informative_alleles=0,
        mendelian_consistency_fraction=None,
        artifacts_sha256={"review.tsv": review_digest},
    )
    loco.summary_json.write_text(json.dumps(summary), encoding="utf-8")
    loco = replace(
        loco,
        review_expected_sha256=review_digest,
        summary_expected_sha256=hashlib.sha256(loco.summary_json.read_bytes()).hexdigest(),
    )
    config = replace(config, family_evidence=replace(family, loco=loco))
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["status"] == "PASS"
    family_stage = config.work_dir / f"run-{result['run_key']}" / "05-family-evidence"
    metrics = json.loads((family_stage / "family_evidence_summary.json").read_bytes())["loco"][
        "materialization"
    ]
    assert metrics["review_rows"] == 3
    assert metrics["matched_genotype_unassessable"] == 2
    assert metrics["projected_review_rows_outside_subject_core"] == 1
    assert metrics["matched_genotype_discordant"] == metrics["matched_genotype_concordant"] == 0
    assert metrics["status_counts"] == {"uninformative_child_genotype": 3}
    family_rows = pq.read_table(
        config.bundle_path / "extensions" / "family_evidence.parquet"
    ).to_pylist()
    matching = sorted((row for row in family_rows if row["pos"] == 11), key=lambda row: row["alt"])
    assert [row["loco_exact_allele_in_model"] for row in matching] == [False, True]
    for row in matching:
        assert row["loco_status"] == row["loco_source_status"] == "uninformative_child_genotype"
        assert row["loco_subject_gt_concordant"] is None
        assert row["de_novo_candidate"] is False
    for counters, total, fraction in [
        ({"source_records": 2, "heldout_carried_alleles": 2, "mendelian_consistent": 2}, 2, 1.0),
        ({"source_records": 2}, 0, None),
    ]:
        summary.update(
            counters=counters,
            informative_alleles=total,
            review_rows=0,
            mendelian_consistency_fraction=fraction,
        )
        loco.summary_json.write_text(json.dumps(summary), encoding="utf-8")
        verified = _load_loco_summary(config)
        assert verified is not None
        assert verified["counters"]["heldout_carried_alleles"] == total
        assert verified["counters"]["uninformative_child_genotype"] == 0
    summary.update(
        validated_unit="TECHNICAL_UNIT",
        excluded_person_units=["SAMPLE_ID", "TECHNICAL_UNIT"],
        model_samples=["FATHER", "MOTHER"],
        model_header_excludes_all_person_units=True,
    )
    loco.summary_json.write_text(json.dumps(summary), encoding="utf-8")
    assert _load_loco_summary(config) is not None
    summary["model_samples"].append("TECHNICAL_UNIT")
    loco.summary_json.write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(InputValidationError, match="validated unit"):
        _load_loco_summary(config)
    summary["model_samples"].pop()
    summary["artifacts_sha256"]["review.tsv"] = "0" * 64
    loco.summary_json.write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(InputValidationError, match="identify the configured review TSV"):
        _load_loco_summary(config)


def test_family_evidence_accepts_separatorless_haploid_normalization_projection(
    tmp_path: Path,
) -> None:
    """A haploid GT has no phase separator; htslib's phase bit is representational."""

    config, _input_vcf = _environment(
        tmp_path,
        family=True,
        family_haploid_projection=True,
    )
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["status"] == "PASS"

    connection = duckdb.connect()
    try:
        rows = connection.execute(
            "SELECT variant_id, subject_gt, subject_phased "
            f"FROM read_parquet('{config.bundle_path / 'extensions' / 'family_evidence.parquet'}/**/*.parquet') "
            "WHERE pos=25 ORDER BY variant_id"
        ).fetchall()
    finally:
        connection.close()
    assert rows == [
        ("chr1:25:A:AG", [1], True),
        ("chr1:25:A:G", [1], True),
    ]


def test_family_evidence_still_rejects_diploid_projection_phase_change(
    tmp_path: Path,
) -> None:
    config, _input_vcf = _environment(
        tmp_path,
        family=True,
        family_projection_phase_mismatch=True,
    )
    sync_pinned_resources(config)
    with pytest.raises(ReconciliationError, match="genotype/phasing differs"):
        build_bundle(config)


def test_family_evidence_rejects_subject_genotype_refinement(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path, family=True, family_genotype_mismatch=True)
    sync_pinned_resources(config)
    with pytest.raises(ReconciliationError, match="genotype/phasing differs"):
        build_bundle(config)


def test_family_evidence_rejects_relative_identifier_in_tags(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path, family=True, family_privacy_leak=True)
    sync_pinned_resources(config)
    with pytest.raises(InputValidationError, match="privacy leak"):
        build_bundle(config)


def test_target_callability_truth_table_and_provenance(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path, target_catalog=True)
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["status"] == "PASS"
    assert result["validation"]["valid"] is True
    bundle = config.bundle_path

    expected = {
        "chr1:10:C:T": "non_reference_observed",
        "chr1:10:C:G": "overlapping_variant_observed",
        "chr1:2:C:G": "reference_observed",
        "chr1:12:T:G": "reference_observed",
        "chr1:13:A:G": "unknown_no_call",
        "chr1:16:T:C": "low_confidence_reference_observation",
        "chr1:21:A:G": "overlapping_variant_observed",
        "chr1:100:T:C": "unknown",
    }
    for locus, conclusion in expected.items():
        assert query_bundle(bundle, "callability", term=locus)["conclusion"] == conclusion
    assert query_bundle(bundle, "callability", term="rs100001")["conclusion"] == (
        "non_reference_observed"
    )

    connection = duckdb.connect()
    try:
        core = connection.execute(
            f"SELECT count(*), count(DISTINCT callability_id) FROM "
            f"read_parquet('{bundle / 'callability.parquet'}')"
        ).fetchone()
        resolutions = dict(
            connection.execute(
                "SELECT resolution, count(*) FROM "
                f"read_parquet('{bundle / 'extensions' / 'target_callability_evidence.parquet'}') "
                "GROUP BY resolution ORDER BY resolution"
            ).fetchall()
        )
        reference_row = connection.execute(
            "SELECT reference_observed, callable, selected_reference_evidence_ids, "
            "observed_ploidy FROM "
            f"read_parquet('{bundle / 'extensions' / 'target_callability_evidence.parquet'}') "
            "WHERE target_id='explicit_reference_block'"
        ).fetchone()
        no_call_row = connection.execute(
            "SELECT reference_observed, callable, overlapping_no_call_evidence_ids FROM "
            f"read_parquet('{bundle / 'extensions' / 'target_callability_evidence.parquet'}') "
            "WHERE target_id='no_call'"
        ).fetchone()
    finally:
        connection.close()
    assert core == (12, 12)
    assert resolutions == {
        "low_confidence_reference": 1,
        "no_call": 1,
        "non_reference_observed": 5,
        "overlapping_variant_observed": 2,
        "reference_observed": 2,
        "unknown": 1,
    }
    assert reference_row[0:2] == (True, True)
    assert len(reference_row[2]) == 1
    assert reference_row[3] == 2
    assert no_call_row[0:2] == (None, False)
    assert len(no_call_row[2]) == 1


def test_target_callability_precedence_is_conservative(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path)
    base = {
        "callability_id": "callability:" + "a" * 64,
        "target_id": "target",
        "target_source": "Synthetic Catalogue",
        "target_source_version": "1.0",
        "catalog_sha256": "b" * 64,
        "catalog_row_digest": "c" * 64,
        "catalog_row_number": 2,
        "chrom": "chr1",
        "pos": 10,
        "end_pos": 10,
        "ref": "C",
        "alt": "T",
        "variant_id": "chr1:10:C:T",
        "gene_symbol": "GENE1",
        "assay_scope": "test",
        "evidence_scope": "site",
    }

    variant_over_nocall = {
        **base,
        "variant_hit_count": 1,
        "exact_target_observed": True,
        "variants_callable": True,
        "variant_confidence": "high",
        "variant_depth": 30,
        "variant_gq": 50,
        "variant_ploidy": 2,
        "observed_variant_ids": ["chr1:10:C:T"],
        "observed_variant_source_digests": ["d" * 64],
        "no_call_count": 1,
        "no_call_evidence_ids": ["no-call"],
    }
    core, evidence = _resolve(variant_over_nocall, config)
    assert evidence["resolution"] == "non_reference_observed"
    assert core["reference_observed"] is False

    no_call_over_reference = {
        **base,
        "no_call_count": 1,
        "no_call_evidence_ids": ["no-call"],
        "reference_count": 1,
        "reference_evidence_ids": ["reference"],
        "references_callable": True,
        "region_count": 1,
        "regions_callable": True,
        "region_ids": ["region"],
    }
    core, evidence = _resolve(no_call_over_reference, config)
    assert evidence["resolution"] == "no_call"
    assert core["reference_observed"] is None
    assert core["callable"] is False

    callable_interval_only = {
        **base,
        "region_count": 1,
        "regions_callable": True,
        "region_ids": ["region"],
    }
    core, evidence = _resolve(callable_interval_only, config)
    assert evidence["resolution"] == "callable_without_reference"
    assert core["reference_observed"] is None
    assert core["callable"] is True


def test_target_catalog_rejects_non_normalized_allele(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path)
    catalog = tmp_path / "non-normalized-targets.tsv"
    catalog.write_text(
        "\t".join(
            [
                "target_id",
                "chrom",
                "pos",
                "ref",
                "alt",
                "gene_symbol",
                "target_source",
                "target_source_version",
                "evidence_scope",
            ]
        )
        + "\nnonminimal\t1\t22\tCG\tCA\tGENE1\tSynthetic Catalogue\t1.0\tsite\n",
        encoding="utf-8",
    )
    from dataclasses import replace

    configured = replace(config, target_catalog_path=catalog)
    raw_vcf = tmp_path / "targets.raw.vcf"
    targets = tmp_path / "targets.parquet"
    rows, allele_rows = _write_catalog_table(
        configured, hashlib.sha256(catalog.read_bytes()).hexdigest(), targets, raw_vcf
    )
    assert (rows, allele_rows) == (1, 1)
    with pytest.raises(InputValidationError, match="non-normalized allele"):
        _validate_normalized_alleles(
            configured,
            raw_vcf,
            tmp_path / "targets.normalized.vcf",
            allele_rows,
        )


def test_pgs_scoring_uses_only_explicit_genotypes_and_calibration(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path, pgs=True)
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["status"] == "PASS"
    assert result["validation"]["valid"] is True
    bundle = config.bundle_path

    query = query_bundle(bundle, "prs", term="synthetic")
    assert query["count"] == 1
    assert query["rows"][0]["training_source"] == "PGS000001"
    assert query["rows"][0]["score_value"] == pytest.approx(2.45)
    assert query["rows"][0]["percentile"] == pytest.approx(68.0)

    connection = duckdb.connect()
    try:
        model = connection.execute(
            "SELECT expected_loci, supported_loci, strand_flipped_loci, "
            "alternate_called_loci, "
            "reference_called_loci, missing_loci, ambiguous_loci, "
            "low_confidence_loci, unsupported_loci, partial_raw_score, "
            "weight_completeness, score_valid, calibrated, publication_status, "
            "contributing_variant_count, contributing_variant_ids_complete, "
            "contributing_variant_ids FROM "
            f"read_parquet('{bundle / 'extensions' / 'pgs_model_audit.parquet'}')"
        ).fetchone()
        loci = connection.execute(
            "SELECT scoring_row_number, effect_is_reference, strand_flipped, status, "
            "dosage, contribution "
            "FROM "
            f"read_parquet('{bundle / 'extensions' / 'pgs_locus_audit.parquet'}') "
            "ORDER BY scoring_row_number"
        ).fetchall()
    finally:
        connection.close()
    assert model[:9] == (6, 5, 1, 2, 3, 0, 0, 0, 1)
    assert model[9] == pytest.approx(2.45)
    assert model[10] == pytest.approx(2.1 / 2.3)
    assert model[11:14] == (True, True, "published")
    assert model[14:16] == (2, True)
    assert model[16] == ["chr1:10:C:T", "chr1:30:C:A"]
    assert loci == [
        (1, False, False, "alternate_called", 1, 1.0),
        (2, True, False, "reference_called", 2, 1.0),
        (3, False, False, "alternate_called", 1, -0.25),
        (4, True, False, "reference_called", 2, 0.5),
        (5, True, True, "reference_called", 2, 0.2),
        (6, None, None, "unsupported", None, None),
    ]

    replay = build_bundle(config)
    assert replay["reused"] is True
    assert replay["run_key"] == result["run_key"]


def test_pgs_caps_embedded_ids_but_preserves_complete_locus_audit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(pgs_module, "MAX_EMBEDDED_CONTRIBUTING_VARIANT_IDS", 1)
    config, _input_vcf = _environment(tmp_path, pgs=True)
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["validation"]["valid"] is True

    connection = duckdb.connect()
    try:
        model = connection.execute(
            "SELECT contributing_variant_count, contributing_variant_ids_complete, "
            "contributing_variant_ids FROM "
            f"read_parquet('{config.bundle_path / 'extensions' / 'pgs_model_audit.parquet'}')"
        ).fetchone()
        normative = connection.execute(
            "SELECT contributing_variants FROM "
            f"read_parquet('{config.bundle_path / 'prs.parquet'}')"
        ).fetchone()
        contributing_loci = connection.execute(
            "SELECT count(DISTINCT variant_id) FROM "
            f"read_parquet('{config.bundle_path / 'extensions' / 'pgs_locus_audit.parquet'}') "
            "WHERE status='alternate_called' AND contribution IS NOT NULL"
        ).fetchone()[0]
    finally:
        connection.close()
    assert model == (2, False, None)
    assert normative == (None,)
    assert contributing_loci == 2


def test_pgs_without_calibration_retains_audit_but_omits_prs(tmp_path: Path) -> None:
    config, _input_vcf = _environment(tmp_path, pgs=True, pgs_calibration=False)
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["validation"]["valid"] is True
    assert not (config.bundle_path / "prs.parquet").exists()
    assert (config.bundle_path / "extensions" / "pgs_model_audit.parquet").is_file()

    query = query_bundle(config.bundle_path, "prs", term="PGS000001")
    assert query["count"] == 1
    assert query["rows"][0]["publication_status"] == "valid_uncalibrated"
    assert "no matching calibrated score" in query["unavailable"]


@pytest.mark.parametrize(
    ("resolution", "expected_status"),
    [
        ("no_call", "missing_no_call"),
        ("callable_without_reference", "missing_interval_only"),
        ("uncallable_without_reference", "missing_interval_only"),
        ("unknown", "missing_unknown"),
        ("overlapping_variant_observed", "ambiguous_other_variant"),
        ("low_confidence_reference", "low_confidence"),
    ],
)
def test_pgs_never_scores_unresolved_or_other_alt_loci(
    resolution: str, expected_status: str
) -> None:
    row = {
        "pgs_id": "PGS000001",
        "scoring_row_number": 1,
        "scoring_row_digest": "a" * 64,
        "supported": True,
        "unsupported_reason": None,
        "target_id": "pgs:PGS000001:row:000000001",
        "callability_id": "callability:" + "b" * 64,
        "target_catalog_sha256": "c" * 64,
        "target_catalog_row_digest": "d" * 64,
        "chrom": "chr1",
        "pos": 13,
        "ref": "A",
        "alt": "G",
        "variant_id": "chr1:13:A:G",
        "effect_allele": "G",
        "other_allele": "A",
        "source_effect_allele": "G",
        "source_other_allele": "A",
        "strand_flipped": False,
        "effect_weight": 2.0,
        "effect_is_reference": False,
        "resolution": resolution,
        "reference_observed": None,
        "callable": resolution == "callable_without_reference",
        "call_confidence": None,
        "depth": None,
        "genotype_quality": None,
        "observed_ploidy": None,
        "observed_variant_ids": (
            ["chr1:13:A:T"] if resolution == "overlapping_variant_observed" else []
        ),
        "variant_gt": None,
    }
    result = _resolve_locus(row)
    assert result["status"] == expected_status
    assert result["dosage"] is None
    assert result["contribution"] is None


def _partial_fixture(tmp_path: Path):
    config, input_vcf = _environment(tmp_path)
    plain = tmp_path / "partial.vcf"
    header = (FIXTURES / "synthetic-reference-confidence.vcf").read_text().split("#CHROM")[0]
    plain.write_text(
        header + '##FORMAT=<ID=NDOS,Number=1,Type=Integer,Description="Native dosage">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID\n"
        "1\t9\t.\tAC\tTG\t60\tPASS\t.\tGT:DP:GQ:NDOS\t1|.:30:60:1\n"
        "1\t15\t.\tG\tC,A\t60\tPASS\t.\tGT:DP:GQ:NDOS\t./2:30:60:1\n"
        "1\t19\t.\tG\tA\t60\tPASS\t.\tGT:DP:GQ:NDOS\t1/.:30:60:1\n"
        "1\t22\t.\tC\tA\t60\tPASS\t.\tGT:DP:GQ:NDOS\t.|1:30:60:1\n"
        "1\t25\t.\tA\tG\t60\tPASS\t.\tGT:DP:GQ:NDOS\t0/.:30:60:1\n"
        "1\t30\t.\tC\tT\t60\tPASS\t.\tGT:DP:GQ:NDOS\t./.:30:60:1\n"
        "1\t32\t.\tT\tG\t60\tPASS\t.\tGT:DP:GQ:NDOS\t.:30:60:1\n"
    )
    pysam.tabix_compress(str(plain), str(input_vcf), force=True)
    pysam.tabix_index(str(input_vcf), preset="vcf", force=True)
    return config, input_vcf


def test_partial_only_bundle_reconciles_source_and_atoms(tmp_path: Path) -> None:
    config, _ = _partial_fixture(tmp_path)
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["validation"]["valid"]
    bundle = config.bundle_path
    partial = bundle / "extensions/observed_partial_variants.parquet"
    with duckdb.connect() as conn:
        assert conn.execute(
            "SELECT count(*) FROM read_parquet(?)", [str(bundle / "variants.parquet/**/*.parquet")]
        ).fetchone() == (0,)
        data = conn.execute(
            "SELECT pos, source_gt, gt, phased, source_alt_index FROM read_parquet(?) ORDER BY pos",
            [str(partial / "**/*.parquet")],
        ).fetchall()
        assert data == [
            (9, "1|.", [1, None], True, 1),
            (10, "1|.", [1, None], True, 1),
            (15, "./2", [None, 1], False, 2),
            (19, "1/.", [1, None], False, 1),
            (22, ".|1", [None, 1], True, 1),
        ]
        assert conn.execute(
            "SELECT count(*) FROM read_parquet(?) WHERE state='no_call' AND NOT callable AND reference_observed IS NULL",
            [str(bundle / "extensions/reference_confidence_regions.parquet/**/*.parquet")],
        ).fetchone() == (7,)
    reconcile = json.loads((bundle / "extensions/reconciliation.json").read_bytes())
    assert reconcile["observed_partial_variant_rows"] == 5
    assert reconcile["unique_source_alleles"] == 4
    assert reconcile["unique_normalized_variants"] == 0


def test_partial_and_complete_sources_can_annotate_the_same_allele(tmp_path: Path) -> None:
    config, vcf = _partial_fixture(tmp_path)
    with pysam.VariantFile(str(vcf)) as f:
        header = str(f.header)
    plain = tmp_path / "overlap.vcf"
    plain.write_text(
        header
        + "".join(
            f"1\t19\t.\tG\tA\t60\tPASS\t.\tGT:DP:GQ:NDOS\t{gt}:30:60:1\n"
            for gt in ("1/.", "./1", "0/1")
        )
    )
    pysam.tabix_compress(str(plain), str(vcf), force=True)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    sync_pinned_resources(config)
    assert build_bundle(config)["validation"]["valid"]
    root = config.bundle_path
    with duckdb.connect() as db:
        assert db.execute(
            "SELECT genotype.gt FROM read_parquet(?)", [str(root / "variants.parquet/**/*.parquet")]
        ).fetchall() == [([0, 1],)]
        assert db.execute(
            "SELECT count(*), count(DISTINCT source_id), count(DISTINCT variant_id) FROM read_parquet(?)",
            [str(root / "extensions/observed_partial_variants.parquet/**/*.parquet")],
        ).fetchone() == (2, 2, 1)


def test_partial_allele_annotations_consume_normal_vep_csq(tmp_path: Path) -> None:
    from genome_annotation_builder.vep_parser import PIPELINE_REQUIRED_FIELDS

    config, input_vcf = _partial_fixture(tmp_path)
    digest = hashlib.sha256(input_vcf.read_bytes()).hexdigest()
    ingest = ingest_reference_confidence_vcf(config, digest, tmp_path / "ingest")
    normalize_projected_vcf(config, tmp_path / "ingest/projected.vcf.gz", tmp_path / "norm")
    annotated = tmp_path / "annotated.vcf.gz"
    fields = sorted(PIPELINE_REQUIRED_FIELDS)
    with pysam.VariantFile(str(tmp_path / "norm/normalized.vcf.gz")) as src:
        hdr = src.header.copy()
        hdr.info.add("CSQ", ".", "String", "Consequence annotations. Format: " + "|".join(fields))
        with pysam.VariantFile(str(annotated), "wz", header=hdr) as out:
            for r in src:
                r.translate(hdr)
                csq = {
                    "ALLELE_NUM": "1",
                    "Allele": r.alts[0],
                    "SYMBOL": "CONTROL",
                    "Consequence": "missense_variant",
                    "IMPACT": "MODERATE",
                    "Feature_type": "Transcript",
                    "Feature": "ENST000001",
                }
                r.info["CSQ"] = ("|".join(csq.get(k, "") for k in fields),)
                out.write(r)
    built = build_variant_datasets(
        replace(config, vep_enabled=True),
        annotated,
        digest,
        ingest["accepted_evidence_multiset"],
        tmp_path / "build",
    )
    assert built["observed_partial_variant_rows"] == 5
    path = tmp_path / "build/observed_partial_variants.parquet/**/*.parquet"
    with duckdb.connect() as conn:
        assert conn.execute(
            "SELECT DISTINCT gene.symbol, consequence.most_severe, len(consequence.all_annotations) FROM read_parquet(?)",
            [str(path)],
        ).fetchall() == [("CONTROL", "missense_variant", 1)]
