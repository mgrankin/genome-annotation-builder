from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.expansionhunter_adapter import (
    ExpansionHunterProvenance,
    materialize_expansionhunter_calls,
    parse_expansionhunter_calls,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

REFERENCE = "1" * 64
PROVENANCE = ExpansionHunterProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="expansionhunter-call-pair",
)


def write_pair(root: Path, *, genotype: str = "2/10") -> tuple[Path, Path]:
    payload = {
        "LocusResults": {
            "SYNTHETIC_LOCUS": {
                "AlleleCount": 2,
                "Coverage": 30.0,
                "LocusId": "SYNTHETIC_LOCUS",
                "Variants": {
                    "small": {
                        "Genotype": "0/1",
                        "ReferenceRegion": "chr1:100-101",
                        "VariantId": "small",
                        "VariantSubtype": "Swap",
                        "VariantType": "SmallVariant",
                    },
                    "repeat": {
                        "CountsOfFlankingReads": "(2, 3), (10, 4)",
                        "CountsOfInrepeatReads": "()",
                        "CountsOfSpanningReads": "(2, 20), (10, 15)",
                        "Genotype": genotype,
                        "GenotypeConfidenceInterval": "2-2/10-10",
                        "ReferenceRegion": "chr1:200-203",
                        "RepeatUnit": "CAG",
                        "VariantId": "repeat",
                        "VariantSubtype": "Repeat",
                        "VariantType": "Repeat",
                    },
                },
            }
        },
        "SampleParameters": {"SampleId": "SYNTHETIC_SAMPLE", "Sex": "Female"},
    }
    json_path = root / "calls.json"
    json_path.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    vcf_path = root / "calls.vcf"
    vcf_path.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeats">
##INFO=<ID=REPID,Number=1,Type=String,Description="Repeat ID">
##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
##INFO=<ID=VARID,Number=1,Type=String,Description="Variant ID">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depth">
##FORMAT=<ID=LC,Number=1,Type=Float,Description="Coverage">
##FORMAT=<ID=SO,Number=1,Type=String,Description="Support">
##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Repeat count">
##FORMAT=<ID=REPCI,Number=1,Type=String,Description="Repeat CI">
##FORMAT=<ID=ADSP,Number=1,Type=String,Description="Spanning">
##FORMAT=<ID=ADFL,Number=1,Type=String,Description="Flanking">
##FORMAT=<ID=ADIR,Number=1,Type=String,Description="In repeat">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t101\t.\tA\tT\t.\tPASS\tVARID=small\tGT:AD:LC\t0/1:20,10:30
chr1\t200\t.\tC\t<STR2>,<STR10>\t.\tPASS\tEND=203;REF=1;RU=CAG;VARID=repeat;REPID=repeat\tGT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC\t1/2:SPANNING/SPANNING:2/10:2-2/10-10:20/15:3/4:0/0:30
""",
        encoding="ascii",
    )
    return json_path, vcf_path


def test_json_vcf_reconciliation_and_typed_projection(tmp_path: Path) -> None:
    source_json, source_vcf = write_pair(tmp_path)
    rows, summary = parse_expansionhunter_calls(source_json, source_vcf, PROVENANCE)
    assert summary["json_variants"] == summary["vcf_records"] == 2
    assert summary["repeat_rows"] == summary["supported_rows"] == 1
    assert summary["omission_reason_counts"] == {"OUTSIDE_REPEAT_EXTENSION:SmallVariant": 1}
    row = rows[0]
    assert row["variant_id"] == "repeat"
    assert row["start_pos"] == 200
    assert row["allele1_repeat_count"] == 2
    assert row["allele2_repeat_count"] == 10
    assert row["allele1_spanning_reads"] == 20
    assert row["allele2_flanking_reads"] == 4
    assert row["spanning_reads"] == 35
    assert row["flanking_reads"] == 7
    assert row["short_read_resolvable"] is True


def test_deterministic_control_ignores_source_sample_binding(tmp_path: Path) -> None:
    source_json, source_vcf = write_pair(tmp_path)
    raw_sample = json.loads(source_json.read_text(encoding="utf-8"))["SampleParameters"]["SampleId"]
    output_one = tmp_path / "one.parquet"
    output_two = tmp_path / "two.parquet"
    result_one = materialize_expansionhunter_calls(
        source_json,
        source_vcf,
        output_one,
        PROVENANCE,
        forbidden_bindings=(raw_sample,),
    )
    result_two = materialize_expansionhunter_calls(
        source_json,
        source_vcf,
        output_two,
        PROVENANCE,
        forbidden_bindings=(raw_sample,),
    )
    assert result_one["parsing"]["repeat_rows"] == 1
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(output_one)
    assert table.schema == schema_for("repeat_expansions")
    row = table.to_pylist()[0]
    assert [row["allele1_repeat_count"], row["allele2_repeat_count"]] == [2, 10]
    assert [row["allele1_spanning_reads"], row["allele2_spanning_reads"]] == [20, 15]
    assert row["locus_coverage"] == pytest.approx(30.0)
    assert raw_sample not in json.dumps(row, sort_keys=True)


def test_mismatch_and_no_call_fail_or_omit_closed(tmp_path: Path) -> None:
    source_json, source_vcf = write_pair(tmp_path)
    payload = json.loads(source_json.read_text(encoding="utf-8"))
    repeat = payload["LocusResults"]["SYNTHETIC_LOCUS"]["Variants"]["repeat"]
    repeat["Genotype"] = "3/10"
    repeat["GenotypeConfidenceInterval"] = "3-3/10-10"
    source_json.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="genotype mismatch"):
        parse_expansionhunter_calls(source_json, source_vcf, PROVENANCE)

    payload = json.loads(source_json.read_text(encoding="utf-8"))
    payload["LocusResults"]["SYNTHETIC_LOCUS"]["Variants"]["repeat"]["Genotype"] = None
    payload["LocusResults"]["SYNTHETIC_LOCUS"]["Variants"]["repeat"][
        "GenotypeConfidenceInterval"
    ] = None
    source_json.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    rows, summary = parse_expansionhunter_calls(source_json, source_vcf, PROVENANCE)
    assert rows == []
    assert summary["extension_required"] is False
    assert summary["omission_reason_counts"] == {
        "NO_REPEAT_CALL": 1,
        "OUTSIDE_REPEAT_EXTENSION:SmallVariant": 1,
    }
    output = tmp_path / "must-not-exist.parquet"
    result = materialize_expansionhunter_calls(source_json, source_vcf, output, PROVENANCE)
    assert result["output"] is None
    assert not output.exists()


def test_open_ended_repeat_is_review_not_exact(tmp_path: Path) -> None:
    source_json, source_vcf = write_pair(tmp_path)
    payload = json.loads(source_json.read_text(encoding="utf-8"))
    repeat = payload["LocusResults"]["SYNTHETIC_LOCUS"]["Variants"]["repeat"]
    repeat["Genotype"] = "2/10+"
    repeat["GenotypeConfidenceInterval"] = "2-2/10+"
    source_json.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    rendered = source_vcf.read_text(encoding="ascii")
    rendered = rendered.replace(":2/10:2-2/10-10:", ":2/10+:2-2/10+:")
    source_vcf.write_text(rendered, encoding="ascii")
    rows, summary = parse_expansionhunter_calls(source_json, source_vcf, PROVENANCE)
    assert summary["review_rows"] == 1
    assert rows[0]["allele2_count_is_lower_bound"] is True
    assert rows[0]["allele2_ci_upper"] is None
    assert rows[0]["short_read_resolvable"] is False


@pytest.mark.parametrize("motif", ["GCN", "AARRG"])
def test_iupac_repeat_motif_is_preserved_and_reviewed(tmp_path: Path, motif: str) -> None:
    source_json, source_vcf = write_pair(tmp_path)
    payload = json.loads(source_json.read_text(encoding="utf-8"))
    payload["LocusResults"]["SYNTHETIC_LOCUS"]["Variants"]["repeat"]["RepeatUnit"] = motif
    source_json.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    source_vcf.write_text(
        source_vcf.read_text(encoding="ascii").replace("RU=CAG", f"RU={motif}"),
        encoding="ascii",
    )

    rows, summary = parse_expansionhunter_calls(source_json, source_vcf, PROVENANCE)

    assert summary["repeat_rows"] == 1
    assert summary["review_rows"] == 1
    assert summary["motif_class_counts"] == {"concrete": 0, "iupac_ambiguous": 1}
    assert rows[0]["motif"] == motif
    assert "CATALOG:IUPAC_AMBIGUOUS_MOTIF" in rows[0]["filter_reasons"]
    assert any("IUPAC ambiguity codes" in item for item in rows[0]["limitations"])


def test_non_iupac_repeat_motif_fails_closed(tmp_path: Path) -> None:
    source_json, source_vcf = write_pair(tmp_path)
    payload = json.loads(source_json.read_text(encoding="utf-8"))
    payload["LocusResults"]["SYNTHETIC_LOCUS"]["Variants"]["repeat"]["RepeatUnit"] = "CAG?"
    source_json.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    source_vcf.write_text(
        source_vcf.read_text(encoding="ascii").replace("RU=CAG", "RU=CAG?"),
        encoding="ascii",
    )

    with pytest.raises(PortableEvidenceError, match="invalid IUPAC DNA symbols"):
        parse_expansionhunter_calls(source_json, source_vcf, PROVENANCE)
