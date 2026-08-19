from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.mt_adapter import (
    MitochondrialProvenance,
    materialize_filtered_mutect_vcf,
    parse_filtered_mutect_vcf,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

REFERENCE = "1" * 64
PROVENANCE = MitochondrialProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="filtered-mitochondrial-vcf",
)


def write_vcf(path: Path, records: str, *, samples: str = "SAMPLE") -> Path:
    path.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=possible_numt,Description="Possible NUMT">
##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description="Allele filter">
##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Strand table">
##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description="Unique reads">
##INFO=<ID=MBQ,Number=R,Type=Integer,Description="Median base quality">
##INFO=<ID=MFRL,Number=R,Type=Integer,Description="Median fragment length">
##INFO=<ID=MMQ,Number=R,Type=Integer,Description="Median mapping quality">
##INFO=<ID=MPOS,Number=A,Type=Integer,Description="Median read position">
##INFO=<ID=OCM,Number=1,Type=Integer,Description="Original alignment mismatch">
##INFO=<ID=POPAF,Number=A,Type=Float,Description="Population AF -log10">
##INFO=<ID=TLOD,Number=A,Type=Float,Description="Tumor LOD">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="F1R2">
##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="F2R1">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Strand bias">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"""
        + samples
        + "\n"
        + records,
        encoding="ascii",
    )
    return path


def test_multiallelic_filter_and_metric_projection(tmp_path: Path) -> None:
    vcf = write_vcf(
        tmp_path / "multi.vcf",
        "chrM\t100\t.\tA\tC,G\t.\tPASS\t"
        "AS_FilterStatus=PASS,possible_numt;AS_SB_TABLE=35,35|5,5|10,10;"
        "AS_UNIQ_ALT_READ_COUNT=9,19;MBQ=30,35,36;MFRL=340,350,360;"
        "MMQ=60,59,58;MPOS=40,41;OCM=0;POPAF=2.4,2.5;TLOD=20,30\t"
        "GT:AD:AF:DP:F1R2:F2R1:SB\t0/1:70,10,20:0.1,0.2:100:70,10,20:0,0,0:0,0,5,5\n",
    )
    rows, summary = parse_filtered_mutect_vcf(vcf, PROVENANCE)
    assert summary["source_records"] == 1
    assert summary["allele_rows"] == 2
    assert summary["supported_rows"] == 1
    assert summary["review_rows"] == 1
    assert rows[0]["alt"] == "C"
    assert rows[0]["allele_fraction"] == pytest.approx(0.1)
    assert rows[0]["forward_alt"] == 5
    assert rows[0]["reverse_alt"] == 5
    assert rows[0]["median_mapping_quality"] == 59
    assert rows[0]["tumor_lod"] == pytest.approx(20)
    assert rows[1]["alt"] == "G"
    assert rows[1]["allele_filter_status"] == ["possible_numt"]
    assert rows[1]["filter_reasons"] == ["ALLELE:possible_numt"]
    assert rows[1]["numt_or_paralogous_context"] is True
    assert rows[1]["f1r2_alt"] == 20


def test_gatk_pipe_delimited_allele_filter_groups_are_lossless(tmp_path: Path) -> None:
    vcf = write_vcf(
        tmp_path / "gatk-filter-groups.vcf",
        "chrM\t200\t.\tA\tC,G,T\t.\tPASS\t"
        "AS_FilterStatus=base_qual,strand_bias|PASS|weak_evidence,possible_numt;"
        "AS_SB_TABLE=20,20|2,2|3,3|4,4;AS_UNIQ_ALT_READ_COUNT=4,6,8;"
        "MBQ=30,31,32,33;MFRL=340,341,342,343;MMQ=60,59,58,57;"
        "MPOS=40,41,42;OCM=0;POPAF=2.1,2.2,2.3;TLOD=10,20,30\t"
        "GT:AD:AF:DP:F1R2:F2R1:SB\t0/1/2/3:40,4,6,10:0.067,0.1,0.167:60:"
        "40,4,6,10:0,0,0,0:0,0,2,2\n"
        "chrM\t300\t.\tG\tA\t.\tPASS\t"
        "AS_FilterStatus=base_qual,strand_bias;AS_SB_TABLE=10,10|5,5;"
        "AS_UNIQ_ALT_READ_COUNT=9;MBQ=30,35;MFRL=340,350;MMQ=60,59;"
        "MPOS=40;OCM=0;POPAF=2.4;TLOD=20\t"
        "GT:AD:AF:DP:F1R2:F2R1:SB\t0/1:20,10:0.333:30:20,10:0,0:0,0,5,5\n",
    )
    rows, summary = parse_filtered_mutect_vcf(vcf, PROVENANCE)
    assert summary["source_records"] == 2
    assert summary["allele_rows"] == 4
    assert [row["allele_filter_status"] for row in rows] == [
        ["base_qual", "strand_bias"],
        ["PASS"],
        ["weak_evidence", "possible_numt"],
        ["base_qual", "strand_bias"],
    ]
    assert rows[1]["filter_reasons"] == []
    assert rows[2]["numt_or_paralogous_context"] is True


def test_real_end_to_end_circular_truth_control_materializes(tmp_path: Path) -> None:
    records = "".join(
        f"chrM\t{position}\t.\tT\tA\t.\tPASS\t"
        "AS_FilterStatus=PASS;AS_SB_TABLE=0,0|50,50;MBQ=0,40;MFRL=0,350;"
        "MMQ=60,60;MPOS=37;POPAF=2.4;TLOD=445\t"
        "GT:AD:AF:DP:F1R2:F2R1:SB\t0/1:0,100:0.99:100:0,100:0,0:0,0,50,50\n"
        for position in (20, 1000, 16550)
    )
    source = write_vcf(tmp_path / "circular-control.vcf", records, samples="CONTROL")
    output = tmp_path / "mitochondrial_heteroplasmy.parquet"
    result = materialize_filtered_mutect_vcf(source, output, PROVENANCE)
    assert result["parsing"]["source_records"] == 3
    assert result["parsing"]["supported_rows"] == 3
    assert result["parsing"]["review_rows"] == 0
    assert result["parsing"]["numt_context_rows"] == 0
    table = pq.read_table(output)
    assert table.schema == schema_for("mitochondrial_heteroplasmy")
    rows = table.to_pylist()
    assert [row["pos"] for row in rows] == [20, 1000, 16550]
    assert all(row["allele_fraction"] == pytest.approx(0.99) for row in rows)
    assert all(row["depth"] == 100 for row in rows)
    assert all(row["forward_alt"] == row["reverse_alt"] == 50 for row in rows)
    assert all(row["circular_boundary_audited"] for row in rows)


@pytest.mark.parametrize(
    ("records", "message"),
    [
        (
            "chrM\t100\t.\tA\tC\t.\tPASS\tAS_FilterStatus=PASS;"
            "AS_SB_TABLE=5,5|5,5;MBQ=30,35;MFRL=340,350;MMQ=60,59;"
            "MPOS=40;POPAF=2.4;TLOD=20\tGT:AD:DP:F1R2:F2R1:SB\t"
            "0/1:10,10:20:10,10:0,0:0,0,5,5\n",
            "FORMAT/AF",
        ),
        (
            "chrM\t16570\t.\tA\tC\t.\tPASS\tAS_FilterStatus=PASS;"
            "AS_SB_TABLE=5,5|5,5;MBQ=30,35;MFRL=340,350;MMQ=60,59;"
            "MPOS=40;POPAF=2.4;TLOD=20\tGT:AD:AF:DP:F1R2:F2R1:SB\t"
            "0/1:10,10:0.5:20:10,10:0,0:0,0,5,5\n",
            "out of bounds",
        ),
    ],
)
def test_missing_metric_and_coordinate_fail_closed(
    tmp_path: Path, records: str, message: str
) -> None:
    vcf = write_vcf(tmp_path / "bad.vcf", records)
    with pytest.raises(PortableEvidenceError, match=message):
        parse_filtered_mutect_vcf(vcf, PROVENANCE)


def test_multisample_and_private_binding_fail_closed(tmp_path: Path) -> None:
    records = (
        "chrM\t100\t.\tA\tC\t.\tPASS\tAS_FilterStatus=PASS;AS_SB_TABLE=5,5|5,5;"
        "MBQ=30,35;MFRL=340,350;MMQ=60,59;MPOS=40;POPAF=2.4;TLOD=20\t"
        "GT:AD:AF:DP:F1R2:F2R1:SB\t0/1:10,10:0.5:20:10,10:0,0:0,0,5,5"
    )
    multi = write_vcf(
        tmp_path / "multi-sample.vcf",
        records + "\t" + records.split("\t")[-1] + "\n",
        samples="ONE\tTWO",
    )
    with pytest.raises(PortableEvidenceError, match="exactly one sample"):
        parse_filtered_mutect_vcf(multi, PROVENANCE)

    private_identifier = "PRIVATE_SUBJECT"
    leaking = MitochondrialProvenance(
        **{
            **PROVENANCE.__dict__,
            "raw_artifact_id": f"artifact-{private_identifier}",
        }
    )
    single = write_vcf(tmp_path / "single.vcf", records + "\n")
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        parse_filtered_mutect_vcf(single, leaking, forbidden_bindings=(private_identifier,))
