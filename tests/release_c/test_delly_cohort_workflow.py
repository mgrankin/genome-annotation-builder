from pathlib import Path

import pytest

from genome_annotation_builder.release_c.delly_cohort_workflow import (
    DellyCohortPolicy,
    DellyTechnicalInput,
    build_delly_cohort_command_plan,
    command_plan_semantic_summary,
    validate_delly_candidate_union,
    validate_force_genotype_outputs,
    validate_sample_free_projection,
    validate_tag_only_germline_filter,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError

HEADER = """##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##ALT=<ID=DEL,Description="Deletion">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=FailGermline,Description="Failed germline filter">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=CT,Number=1,Type=String,Description="Orientation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CI">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="CI">
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="Allele ID">
##INFO=<ID=NALLELE,Number=1,Type=Integer,Description="Alleles">
##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Support">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Reference pairs">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Variant pairs">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference splits">
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Variant splits">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Depth">
"""


def write_discovery(path: Path, *, sample_free: bool, moved: bool = False) -> Path:
    sample = "" if sample_free else "\tFORMAT\tSYNTHETIC"
    start = 103 if moved else 100
    path.write_text(
        HEADER
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{sample}\n"
        + f"1\t{start}\tdiscovery-1\tA\t<DEL>\t500\tPASS\t"
        + "SVTYPE=DEL;END=300;SVLEN=-200;CT=3to5;PRECISE;CIPOS=0,0;CIEND=0,0"
        + ("\n" if sample_free else "\tGT:GQ:PL:DR:DV:RR:RV:RC\t0/1:60:60,0,90:20:8:14:6:32\n")
    )
    return path


def write_sites(path: Path, *, end: int = 300) -> Path:
    path.write_text(
        HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + f"1\t100\tsite-1\tA\t<DEL>\t500\tPASS\tSVTYPE=DEL;END={end};SVLEN=-200;"
        + "CT=3to5;PRECISE;CIPOS=0,3;CIEND=0,0;ALLELEID=0;NALLELE=1;SUPP=1\n"
    )
    return path


def write_genotyped(
    path: Path, *, sample: str, filter_value: str = "PASS", gt: str = "0/1"
) -> Path:
    path.write_text(
        HEADER
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + f"1\t100\tsite-1\tA\t<DEL>\t500\t{filter_value}\t"
        + "SVTYPE=DEL;END=300;SVLEN=-200;CT=3to5;PRECISE;CIPOS=0,3;CIEND=0,0;"
        + "ALLELEID=0;NALLELE=1;SUPP=1\tGT:GQ:PL:DR:DV:RR:RV:RC\t"
        + f"{gt}:60:60,0,90:20:8:14:6:32\n"
    )
    return path


def write_multisample(path: Path, *, filter_value: str = "PASS", changed_gt: bool = False) -> Path:
    second_gt = "1/1" if changed_gt else "0/0"
    path.write_text(
        HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tONE\tTWO\n"
        + f"1\t100\tsite-1\tA\t<DEL>\t500\t{filter_value}\t"
        + "SVTYPE=DEL;END=300;SVLEN=-200;CT=3to5;PRECISE;CIPOS=0,3;CIEND=0,0;"
        + "ALLELEID=0;NALLELE=1;SUPP=1\tGT:GQ:PL:DR:DV:RR:RV:RC\t"
        + f"0/1:60:60,0,90:20:8:14:6:32\t{second_gt}:55:0,55,90:28:0:20:0:31\n"
    )
    return path


def test_command_plan_is_loss_minimized_and_summary_is_private_free(tmp_path: Path) -> None:
    inputs = [
        DellyTechnicalInput(
            discovery_vcf=tmp_path / "source.bcf",
            sample_free_vcf=tmp_path / "source.sites.bcf",
            alignment=tmp_path / "source.bam",
            sample_name="SYNTHETIC",
            genotyped_vcf=tmp_path / "source.geno.bcf",
            tagged_single_sample_vcf=tmp_path / "source.tagged.bcf",
        )
    ]
    plan = build_delly_cohort_command_plan(
        delly="delly",
        bcftools="bcftools",
        reference_fasta=tmp_path / "reference.fa",
        exclude_regions=tmp_path / "exclude.tsv",
        technical_inputs=inputs,
        candidate_sites_vcf=tmp_path / "sites.bcf",
        raw_multisample_vcf=tmp_path / "raw.bcf",
        tagged_multisample_vcf=tmp_path / "tagged.bcf",
    )
    merge = next(item["argv"] for item in plan if item["stage"] == "merge_candidate_sites")
    tagged = next(
        item["argv"]
        for item in plan
        if item["stage"] == "tag_germline_filters_without_refinement_or_collapse"
    )
    assert merge[merge.index("-y") + 1] == "0"
    assert merge[merge.index("-a") + 1] == "0.0"
    assert merge[merge.index("-v") + 1] == "0"
    assert merge[merge.index("--rep-min-af") + 1] == "0.0"
    assert {"-t", "--no-collapse", "--no-refine"}.issubset(tagged)
    summary = command_plan_semantic_summary(plan)
    assert summary["status"] == "PASS"
    assert summary["argv_persisted"] is False
    assert "SYNTHETIC" not in str(summary)


def test_policy_rejects_lossy_candidate_merge() -> None:
    with pytest.raises(PortableEvidenceError):
        DellyCohortPolicy(merge_min_quality=200)


def test_projection_and_candidate_union_reconcile_every_source_record(tmp_path: Path) -> None:
    discovery = write_discovery(tmp_path / "discovery.vcf", sample_free=False)
    projected = write_discovery(tmp_path / "projected.vcf", sample_free=True)
    sites = write_sites(tmp_path / "sites.vcf")
    projection = validate_sample_free_projection(discovery, projected)
    union = validate_delly_candidate_union([projected], sites)
    assert projection["records"] == 1
    assert union["projected_discovery_records"] == union["candidate_sites"] == 1
    assert union["unmatched_discovery_records"] == 0
    assert union["sample_identifier_persisted"] is False


def test_projection_change_and_unmatched_candidate_fail_closed(tmp_path: Path) -> None:
    discovery = write_discovery(tmp_path / "discovery.vcf", sample_free=False)
    changed = write_discovery(tmp_path / "changed.vcf", sample_free=True, moved=True)
    with pytest.raises(PortableEvidenceError):
        validate_sample_free_projection(discovery, changed)
    sites = write_sites(tmp_path / "sites.vcf", end=900)
    with pytest.raises(PortableEvidenceError):
        validate_delly_candidate_union([changed], sites)


def test_force_genotype_requires_exact_site_universe_and_unique_samples(tmp_path: Path) -> None:
    sites = write_sites(tmp_path / "sites.vcf")
    first = write_genotyped(tmp_path / "first.vcf", sample="ONE")
    second = write_genotyped(tmp_path / "second.vcf", sample="TWO", gt="0/0")
    summary = validate_force_genotype_outputs(sites, [first, second])
    assert summary["expected_genotype_rows"] == summary["validated_genotype_rows"] == 2
    duplicate = write_genotyped(tmp_path / "duplicate.vcf", sample="ONE", gt="0/0")
    with pytest.raises(PortableEvidenceError):
        validate_force_genotype_outputs(sites, [first, duplicate])


def test_tag_only_filter_preserves_biology_and_detects_gt_change(tmp_path: Path) -> None:
    raw = write_multisample(tmp_path / "raw.vcf")
    tagged = write_multisample(tmp_path / "tagged.vcf", filter_value="FailGermline")
    summary = validate_tag_only_germline_filter(raw, tagged)
    assert summary["removed_sites"] == 0
    assert summary["changed_genotypes_or_evidence"] == 0
    assert summary["tagged_filter_counts"] == {"FailGermline": 1}
    assert summary["biological_payload_sha256_before"] == summary["biological_payload_sha256_after"]
    changed = write_multisample(tmp_path / "changed.vcf", changed_gt=True)
    with pytest.raises(PortableEvidenceError):
        validate_tag_only_germline_filter(raw, changed)


def test_small_delly_cohort_smoke_is_reconciled(tmp_path: Path) -> None:
    sites = write_sites(tmp_path / "sites.vcf")
    genotyped = write_genotyped(tmp_path / "genotyped.vcf", sample="SYNTHETIC")
    summary = validate_force_genotype_outputs(sites, [genotyped])
    assert summary["candidate_sites"] == 1
    assert summary["missing_or_extra_rows"] == 0
