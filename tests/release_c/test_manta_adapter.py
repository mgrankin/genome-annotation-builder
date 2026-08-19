from __future__ import annotations

from pathlib import Path

import pytest

from genome_annotation_builder.release_c.manta_adapter import (
    MantaProvenance,
    parse_manta_germline_calls,
    validate_manta_bnd_pairing,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError

REFERENCE = "1" * 64
PROVENANCE = MantaProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="manta-diploid-sv-vcf",
)


def write_vcf(path: Path, *, reciprocal: bool = True) -> Path:
    mate_for_second = "BND1" if reciprocal else "MISSING"
    path.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CIPOS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="CIEND">
##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Homology">
##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description="Insertion length">
##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=PAIR_COUNT,Number=1,Type=Integer,Description="Pair count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Pair support">
##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split support">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t100\tDEL1\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=300;SVLEN=-200;CIPOS=-2,2;CIEND=-3,3;HOMLEN=1;PAIR_COUNT=8\tGT:FT:GQ:PL:PR:SR\t0/1:PASS:50:50,0,90:20,8:15,5
chr1\t500\tBND1\tG\tG]chr2:700]\t70\tPASS\tSVTYPE=BND;MATEID=BND2;CIPOS=0,0;HOMLEN=0;BND_DEPTH=30;PAIR_COUNT=4\tGT:FT:GQ:PL:PR:SR\t0/1:PASS:45:45,0,80:20,4:15,3
chr2\t700\tBND2\tT\tT]chr1:500]\t70\tPASS\tSVTYPE=BND;MATEID="""
        + mate_for_second
        + """;CIPOS=0,0;HOMLEN=0;BND_DEPTH=28;PAIR_COUNT=4\tGT:FT:GQ:PL:PR:SR\t0/1:PASS:45:45,0,80:20,4:15,3
""",
        encoding="ascii",
    )
    return path


def write_three_bnd_pairs(path: Path) -> Path:
    records = []
    for index, offset in enumerate((100, 300, 500), start=1):
        left = f"PAIR{index}A"
        right = f"PAIR{index}B"
        records.extend(
            [
                f"chr1\t{offset}\t{left}\tA\tA]chr2:{offset + 50}]\t.\tPASS\t"
                f"SVTYPE=BND;MATEID={right}",
                f"chr2\t{offset + 50}\t{right}\tC\tC]chr1:{offset}]\t.\tPASS\t"
                f"SVTYPE=BND;MATEID={left}",
            ]
        )
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "##contig=<ID=chr2,length=1000000>\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(records) + "\n",
        encoding="ascii",
    )
    return path


def write_noncanonical_vcf(path: Path, *, reciprocal: bool = True) -> Path:
    noncanonical_mate = "CORE_ALT" if reciprocal else "MISSING"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "##contig=<ID=chrUn_KI270442v1,length=392061>\n"
        "##contig=<ID=chr2_KI270774v1_alt,length=223625>\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t100\tCORE_DEL\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200\tGT\t0/1\n"
        "chrUn_KI270442v1\t300\tNONCORE_INS\tA\t<INS>\t.\tPASS\t"
        "SVTYPE=INS;END=300\tGT\t0/1\n"
        "chr1\t500\tCORE_ALT\tA\tA]chrUn_KI270442v1:700]\t.\tPASS\t"
        "SVTYPE=BND;MATEID=NONCORE_ALT\tGT\t0/1\n"
        "chrUn_KI270442v1\t700\tNONCORE_ALT\tC\tC]chr1:500]\t.\tPASS\t"
        f"SVTYPE=BND;MATEID={noncanonical_mate}\tGT\t0/1\n"
        "chrUn_KI270442v1\t900\tNONCORE_A\tG\tG]chr2_KI270774v1_alt:1100]"
        "\t.\tPASS\tSVTYPE=BND;MATEID=NONCORE_B\tGT\t0/1\n"
        "chr2_KI270774v1_alt\t1100\tNONCORE_B\tT\tT]chrUn_KI270442v1:900]"
        "\t.\tPASS\tSVTYPE=BND;MATEID=NONCORE_A\tGT\t0/1\n",
        encoding="ascii",
    )
    return path


def test_typed_manta_rows_and_reciprocal_event_identity(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "manta.vcf")
    rows, summary = parse_manta_germline_calls(source, PROVENANCE)
    assert summary["source_records"] == summary["typed_rows"] == 3
    assert summary["record_reconciliation_delta"] == 0
    assert summary["bnd_records"] == 2
    assert summary["reciprocal_bnd_pairs"] == 1
    deletion = next(row for row in rows if row["caller_record_id"] == "DEL1")
    assert deletion["paired_end_support"] == 8
    assert deletion["paired_end_reference_reads"] == 20
    assert deletion["split_read_support"] == 5
    bnds = [row for row in rows if row["sv_type"] == "BND"]
    assert bnds[0]["event_id"] == bnds[1]["event_id"]
    assert {row["bnd_mate_id"] for row in bnds} == {"BND1", "BND2"}
    assert all(row["mate_link_status"] == "paired_mate_id" for row in bnds)


def test_three_pair_bnd_mechanics_control(tmp_path: Path) -> None:
    source = write_three_bnd_pairs(tmp_path / "three-pairs.vcf")
    summary = validate_manta_bnd_pairing(source)
    assert summary["source_records"] == summary["bnd_records"] == 6
    assert summary["reciprocal_bnd_pairs"] == 3
    assert summary["sample_identifiers_persisted"] is False


def test_nonreciprocal_mate_fails_closed(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "broken.vcf", reciprocal=False)
    with pytest.raises(PortableEvidenceError, match="not reciprocal"):
        parse_manta_germline_calls(source, PROVENANCE)


def test_noncanonical_records_are_omitted_as_complete_pairs(tmp_path: Path) -> None:
    source = write_noncanonical_vcf(tmp_path / "noncanonical.vcf")
    rows, summary = parse_manta_germline_calls(source, PROVENANCE)
    assert [row["caller_record_id"] for row in rows] == ["CORE_DEL"]
    assert summary["source_records"] == 6
    assert summary["typed_rows"] == 1
    assert summary["record_reconciliation_delta"] == 5
    assert summary["omitted_noncanonical_records"] == 5
    assert summary["omitted_noncanonical_non_bnd_records"] == 1
    assert summary["omitted_noncanonical_bnd_records"] == 4
    assert summary["omitted_noncanonical_bnd_pairs"] == 2
    assert summary["core_to_noncanonical_bnd_pairs"] == 1
    assert summary["noncanonical_to_noncanonical_bnd_pairs"] == 1
    assert summary["source_bnd_records"] == 4
    assert summary["source_reciprocal_bnd_pairs"] == 2
    assert summary["bnd_records"] == 0
    assert summary["reciprocal_bnd_pairs"] == 0


def test_nonreciprocal_noncanonical_pair_still_fails_closed(tmp_path: Path) -> None:
    source = write_noncanonical_vcf(tmp_path / "broken-noncanonical.vcf", reciprocal=False)
    with pytest.raises(PortableEvidenceError, match="not reciprocal"):
        parse_manta_germline_calls(source, PROVENANCE)
