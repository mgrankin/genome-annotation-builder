from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.cyrius_adapter import (
    CyriusProvenance,
    materialize_cyrius_results,
    parse_cyrius_results,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

SAMPLE = "SYNTHETIC_SUBJECT"
REFERENCE = "1" * 64
PROVENANCE = CyriusProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="cyrius-result-set",
    definition_version="Cyrius 1.1.1 separately installed definition snapshot",
    definition_identity_sha256="5" * 64,
)


def base_result(
    *, genotype: str | None = "*1/*2", caller_filter: str | None = "PASS"
) -> dict[str, object]:
    return {
        "Coverage_MAD": 0.12,
        "Median_depth": 31.5,
        "Total_CN": 4,
        "Spacer_CN": 2,
        "Total_CN_raw": 4.08,
        "Spacer_CN_raw": 1.96,
        "Variants_called": ["g.100A>G", "g.200C>T"],
        "CNV_group": "cn2",
        "Genotype": genotype,
        "Filter": caller_filter,
        "Raw_star_allele": "*1_*2",
        "Call_info": "matched",
        "Exon9_CN": 2,
        "CNV_consensus": "2,2,2",
        "d67_snp_call": "2,2,2",
        "d67_snp_raw": "2.01,1.98,2.00",
        "Variant_raw_count": {
            "g.100A>G": "10(6:4),20",
            "g.200C>T": "8,22",
        },
    }


def write_results(
    root: Path,
    result: dict[str, object],
    *,
    tsv_genotype: str | None = None,
    tsv_filter: str | None = None,
) -> tuple[Path, Path]:
    json_path = root / "result.json"
    tsv_path = root / "result.tsv"
    json_path.write_text(json.dumps({SAMPLE: result}, sort_keys=True), encoding="utf-8")
    genotype = result["Genotype"] if tsv_genotype is None else tsv_genotype
    caller_filter = result["Filter"] if tsv_filter is None else tsv_filter
    tsv_path.write_text(
        "Sample\tGenotype\tFilter\n"
        f"{SAMPLE}\t{genotype if genotype is not None else 'None'}\t"
        f"{caller_filter if caller_filter is not None else 'None'}\n",
        encoding="utf-8",
    )
    return json_path, tsv_path


def test_pass_call_is_lossless_deterministic_and_sample_private(tmp_path: Path) -> None:
    json_path, tsv_path = write_results(tmp_path, base_result())
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    result_one = materialize_cyrius_results(json_path, tsv_path, first, PROVENANCE)
    result_two = materialize_cyrius_results(json_path, tsv_path, second, PROVENANCE)
    assert result_one["parsing"]["call_status"] == "called"
    assert result_one["parsing"]["raw_variant_count_rows"] == 2
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(first)
    assert table.schema == schema_for("complex_pharmacogenes")
    row = table.to_pylist()[0]
    assert row["gene"] == "CYP2D6"
    assert row["diplotype"] == "*1/*2"
    assert row["candidate_diplotypes"] == ["*1/*2"]
    assert row["star_alleles"] == ["*1", "*2"]
    assert row["gene_copy_number"] is None
    assert row["locus_total_copy_number"] == pytest.approx(4)
    assert row["paralog_marker_copy_number"] == pytest.approx(2)
    assert row["outside_call_eligible"] is True
    assert row["licence_gated_caller"] is True
    assert row["activity_score"] is row["phenotype"] is None
    assert row["raw_variant_counts"][0] == {
        "variant_id": "g.100A>G",
        "alt_reads": 10,
        "ref_reads": 20,
        "forward_alt_reads": 6,
        "reverse_alt_reads": 4,
    }
    rendered = json.dumps({"row": row, "summary": result_one["parsing"]}, sort_keys=True)
    assert SAMPLE not in rendered
    assert "/raid/" not in rendered


@pytest.mark.parametrize(
    ("genotype", "caller_filter", "candidate_count", "selected"),
    [
        ("*1/*46;*43/*45", "More_than_one_possible_genotype", 2, None),
        ("*1_*2_*68", "Not_assigned_to_haplotypes", 1, None),
        ("*1/*2x5", "LowQ_high_CN", 1, "*1/*2x5"),
    ],
)
def test_ambiguous_filter_states_never_become_outside_calls(
    tmp_path: Path,
    genotype: str,
    caller_filter: str,
    candidate_count: int,
    selected: str | None,
) -> None:
    json_path, tsv_path = write_results(
        tmp_path, base_result(genotype=genotype, caller_filter=caller_filter)
    )
    rows, summary = parse_cyrius_results(json_path, tsv_path, PROVENANCE)
    row = rows[0]
    assert summary["call_status"] == "ambiguous"
    assert len(row["candidate_diplotypes"]) == candidate_count
    assert row["diplotype"] == selected
    assert row["ambiguous"] is True
    assert row["outside_call_eligible"] is False
    assert row["confidence_tier"] == "REVIEW"


def test_no_call_is_explicit_and_not_promoted(tmp_path: Path) -> None:
    result = base_result(genotype=None, caller_filter=None)
    for field in (
        "Total_CN",
        "Spacer_CN",
        "Total_CN_raw",
        "Spacer_CN_raw",
        "Variants_called",
        "CNV_group",
        "Raw_star_allele",
        "Call_info",
        "Exon9_CN",
        "CNV_consensus",
        "d67_snp_call",
        "d67_snp_raw",
        "Variant_raw_count",
    ):
        result[field] = None
    json_path, tsv_path = write_results(tmp_path, result)
    rows, _ = parse_cyrius_results(json_path, tsv_path, PROVENANCE)
    row = rows[0]
    assert row["call_status"] == "no_call"
    assert row["diplotype"] is None
    assert row["candidate_diplotypes"] == []
    assert row["star_alleles"] == []
    assert row["confidence_tier"] == "UNRESOLVED"
    assert row["outside_call_eligible"] is False


def test_schema_cross_artifact_counts_licence_and_privacy_fail_closed(
    tmp_path: Path,
) -> None:
    result = base_result()
    json_path, tsv_path = write_results(tmp_path, result, tsv_genotype="*1/*4")
    with pytest.raises(PortableEvidenceError, match="genotype/filter differ"):
        parse_cyrius_results(json_path, tsv_path, PROVENANCE)

    result["Unexpected"] = "value"
    json_path, tsv_path = write_results(tmp_path, result)
    with pytest.raises(PortableEvidenceError, match="schema differs"):
        parse_cyrius_results(json_path, tsv_path, PROVENANCE)
    result.pop("Unexpected")

    result["Variant_raw_count"] = {"g.100A>G": "10(9:9),20"}
    json_path, tsv_path = write_results(tmp_path, result)
    with pytest.raises(PortableEvidenceError, match="strand counts disagree"):
        parse_cyrius_results(json_path, tsv_path, PROVENANCE)

    result = base_result()
    json_path, tsv_path = write_results(tmp_path, result)
    no_licence_boundary = CyriusProvenance(**{**PROVENANCE.__dict__, "licence_gated_caller": False})
    with pytest.raises(PortableEvidenceError, match="licence-gated boundary"):
        parse_cyrius_results(json_path, tsv_path, no_licence_boundary)

    leaking = CyriusProvenance(**{**PROVENANCE.__dict__, "raw_artifact_id": f"cyrius-{SAMPLE}"})
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        parse_cyrius_results(json_path, tsv_path, leaking)
