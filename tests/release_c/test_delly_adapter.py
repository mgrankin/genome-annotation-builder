from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pysam
import pytest

from genome_annotation_builder.release_c.delly_adapter import (
    DellyProvenance,
    materialize_delly_calls,
    parse_delly_calls,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

REFERENCE = "1" * 64
PROVENANCE = DellyProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="delly-single-sample-bcf",
)


def write_vcf(path: Path, *, broken_bnd: bool = False) -> Path:
    partner = "chr2:500" if not broken_bnd else "chr2:501"
    path.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="CIEND">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CIPOS">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Partner chromosome">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Partner position">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=PE,Number=1,Type=Integer,Description="PE">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="MAPQ">
##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description="SR MAPQ">
##INFO=<ID=SR,Number=1,Type=Integer,Description="SR">
##INFO=<ID=CT,Number=1,Type=String,Description="Orientation">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Insertion length">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Homology length">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Filter">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=RDCN,Number=1,Type=Integer,Description="Copy number">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Reference pairs">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Variant pairs">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference splits">
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Variant splits">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t100\tDEL1\tA\t<DEL>\t80\tPASS\tIMPRECISE;SVTYPE=DEL;END=300;SVLEN=-200;PE=9;SR=7;MAPQ=60;SRMAPQ=55;CT=3to5;CIPOS=-5,5;CIEND=-10,10;INSLEN=0;HOMLEN=2\tGT:PL:GQ:FT:RC:RDCN:DR:DV:RR:RV\t0/1:80,0,120:60:PASS:30:1:20:12:15:8
chr1\t400\tBND1\tG\t["""
        + partner
        + """[G\t70\tPASS\tPRECISE;SVTYPE=BND;END=400;CHR2=chr2;POS2=500;PE=0;SR=6;MAPQ=0;SRMAPQ=60;CT=5to5;CIPOS=0,0;CIEND=0,0;INSLEN=0;HOMLEN=0\tGT:PL:GQ:FT:RC:RDCN:DR:DV:RR:RV\t0/1:70,0,100:50:PASS:25:2:0:0:12:6
chr1\t700\tDEL2\tC\t<DEL>\t20\tLowQual\tIMPRECISE;SVTYPE=DEL;END=900;SVLEN=-200;PE=2;SR=0;MAPQ=40;CT=3to5;CIPOS=-20,20;CIEND=-20,20;INSLEN=0;HOMLEN=0\tGT:PL:GQ:FT:RC:RDCN:DR:DV:RR:RV\t0/0:0,20,100:20:LowQual:30:2:30:0:0:0
""",
        encoding="ascii",
    )
    return path


def test_delly_projection_preserves_call_status_and_bnd(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "calls.vcf")
    rows, summary = parse_delly_calls(source, PROVENANCE)
    assert summary["source_records"] == summary["typed_rows"] == 3
    assert summary["record_reconciliation_delta"] == 0
    assert summary["call_status_counts"] == {
        "non_reference": 2,
        "reference": 1,
        "no_call": 0,
    }
    deletion = next(row for row in rows if row["caller_record_id"] == "DEL1")
    assert deletion["genotype"] == [0, 1]
    assert deletion["genotype_likelihoods"] == [80, 0, 120]
    assert deletion["paired_end_support"] == 12
    assert deletion["split_read_support"] == 8
    assert deletion["caller_site_paired_end_support"] == 9
    assert deletion["caller_site_split_read_support"] == 7
    assert deletion["cipos"] == [-5, 5]
    bnd = next(row for row in rows if row["caller_record_id"] == "BND1")
    assert (bnd["chrom2"], bnd["pos2"]) == ("chr2", 500)
    assert bnd["mate_link_status"] == "single_record_partner_coordinates"
    reference = next(row for row in rows if row["caller_record_id"] == "DEL2")
    assert reference["call_status"] == "reference"
    assert reference["confidence_tier"] == "REVIEW"
    assert reference["filter_reasons"] == ["GENOTYPE:LowQual", "SITE:LowQual"]


def test_materialization_is_deterministic_and_source_sample_free(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "synthetic-delly.vcf")
    with pysam.VariantFile(str(source)) as handle:
        raw_sample = next(iter(handle.header.samples))
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    result_one = materialize_delly_calls(
        source, first, PROVENANCE, forbidden_bindings=(raw_sample,)
    )
    result_two = materialize_delly_calls(
        source, second, PROVENANCE, forbidden_bindings=(raw_sample,)
    )
    assert result_one["parsing"]["source_records"] == 3
    assert result_one["parsing"]["record_reconciliation_delta"] == 0
    assert result_one["parsing"]["call_status_counts"] == {
        "non_reference": 2,
        "reference": 1,
        "no_call": 0,
    }
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(first)
    assert table.schema == schema_for("structural_variants")
    rendered = json.dumps(table.to_pylist(), sort_keys=True)
    assert raw_sample not in rendered
    bnd = next(row for row in table.to_pylist() if row["sv_type"] == "BND")
    assert (bnd["chrom1"], bnd["pos1"], bnd["chrom2"], bnd["pos2"]) == (
        "chr1",
        400,
        "chr2",
        500,
    )


def test_bnd_partner_disagreement_fails_closed(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "broken.vcf", broken_bnd=True)
    with pytest.raises(PortableEvidenceError, match="ALT and CHR2/POS2 disagree"):
        parse_delly_calls(source, PROVENANCE)


def test_delly_sex_chromosome_genotype_is_forced_to_review(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "autosomal.vcf")
    sex_chromosome = tmp_path / "sex-chromosome.vcf"
    payload = source.read_text().replace("chr1", "chrX")
    payload = payload.replace(":25:2:0:0:12:6", ":25:-1:0:0:12:6")
    sex_chromosome.write_text(payload, encoding="ascii")
    rows, summary = parse_delly_calls(sex_chromosome, PROVENANCE)
    assert rows
    for row in rows:
        assert row["confidence_tier"] == "REVIEW"
        assert "MODEL:DIPLOID_ONLY_SEX_CHROMOSOME" in row["filter_reasons"]
        assert any("assumes diploidy" in value for value in row["limitations"])
    bnd = next(row for row in rows if row["caller_record_id"] == "BND1")
    assert bnd["copy_number"] is None
    assert any("RDCN=-1" in value for value in bnd["limitations"])
    assert summary["rdcn_missing_sentinel_rows"] == 1


def test_delly_rejects_unknown_negative_copy_number(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "negative-copy-number.vcf")
    payload = source.read_text().replace(":25:2:0:0:12:6", ":25:-2:0:0:12:6")
    source.write_text(payload, encoding="ascii")
    with pytest.raises(PortableEvidenceError, match="invalid negative DELLY FORMAT/RDCN"):
        parse_delly_calls(source, PROVENANCE)
