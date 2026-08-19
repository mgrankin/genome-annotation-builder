from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.aldy_adapter import (
    AldyProvenance,
    materialize_aldy_results,
    parse_aldy_results,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

SAMPLE = "SYNTHETIC_SUBJECT"
REFERENCE = "1" * 64
DEFINITION = "CYP2D6-synthetic-definition-1"
PROVENANCE = AldyProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="aldy-result-set",
    definition_version=DEFINITION,
    definition_identity_sha256="5" * 64,
)
HEADER = (
    "#Sample\tGene\tSolutionID\tMajor\tMinor\tCopy\tAllele\tLocation\tType"
    "\tCoverage\tEffect\tdbSNP\tCode\tStatus"
)


def mutation_row(
    solution: int,
    major: str,
    minor: str,
    copy: int,
    allele: str,
    position: int,
    change: str,
    coverage: float,
    effect: str,
    dbsnp: str,
    *,
    sample: str = SAMPLE,
) -> str:
    # Aldy 4.8.3 mutation rows have 13 cells: the final Status cell is omitted.
    fields = [
        sample,
        "CYP2D6",
        str(solution),
        major,
        minor,
        str(copy),
        allele,
        str(position),
        change,
        str(coverage),
        effect,
        dbsnp,
        "",
    ]
    assert len(fields) == 13
    return "\t".join(fields)


def solution_block(
    solution: int,
    major: str,
    minor: str,
    expression: str,
    *,
    phenotype: str = "normal",
    score: float = 2.0,
    sample: str = SAMPLE,
) -> list[str]:
    alleles = minor.split(";")
    rows = [
        f"#Solution {solution}: {expression}; gene={DEFINITION}; "
        f"cpic={phenotype}; cpic_score={score}"
    ]
    for copy, allele in enumerate(alleles):
        rows.append(
            mutation_row(
                solution,
                major,
                minor,
                copy,
                allele,
                100 + copy,
                "A>G" if copy == 0 else "C>T",
                31 + copy,
                "synthetic_effect",
                f"rsSYNTHETIC{solution}{copy}",
                sample=sample,
            )
        )
    return rows


def write_output(path: Path, blocks: list[list[str]]) -> Path:
    path.write_text(
        "\n".join([HEADER, *(line for block in blocks for line in block)]) + "\n",
        encoding="utf-8",
    )
    return path


def test_unique_solution_is_lossless_deterministic_and_private(tmp_path: Path) -> None:
    source = write_output(
        tmp_path / "result.aldy",
        [solution_block(1, "*1/*2", "1.001;2.001", "*1.001, *2.001")],
    )
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    one = materialize_aldy_results(source, first, PROVENANCE)
    two = materialize_aldy_results(source, second, PROVENANCE)
    assert one["output"]["parquet_sha256"] == two["output"]["parquet_sha256"]
    assert one["parsing"]["solution_count"] == 1
    assert one["parsing"]["allele_decomposition_rows"] == 2
    table = pq.read_table(first)
    assert table.schema == schema_for("complex_pharmacogenes")
    row = table.to_pylist()[0]
    assert row["call_status"] == "called"
    assert row["diplotype"] == "*1/*2"
    assert row["candidate_diplotypes"] == ["*1/*2"]
    assert row["star_alleles"] == ["*1", "*2"]
    assert row["gene_copy_number"] == pytest.approx(2)
    assert row["phenotype"] == "normal"
    assert row["activity_score"] == pytest.approx(2)
    assert row["outside_call_eligible"] is True
    assert row["candidate_solution_details"][0]["minor_alleles"] == ["1.001", "2.001"]
    assert row["allele_decomposition"][0]["position_zero_based"] == 100
    assert row["allele_decomposition"][0]["status"] is None
    rendered = json.dumps({"row": row, "summary": one["parsing"]}, sort_keys=True)
    assert SAMPLE not in rendered
    assert "/raid/" not in rendered


def test_multiple_optimal_solutions_remain_ambiguous(tmp_path: Path) -> None:
    source = write_output(
        tmp_path / "result.aldy",
        [
            solution_block(1, "*1/*2", "1.001;2.001", "*1.001, *2.001"),
            solution_block(
                2,
                "*1/*4+*4X",
                "1.001;4.001;4.002X",
                "*1.001, *4.001, *4.002X",
                phenotype="intermediate",
                score=1.0,
            ),
        ],
    )
    rows, summary = parse_aldy_results(source, PROVENANCE)
    row = rows[0]
    assert summary["solution_count"] == 2
    assert row["call_status"] == "ambiguous"
    assert row["diplotype"] is None
    assert row["candidate_diplotypes"] == ["*1/*2", "*1/*4+*4X"]
    assert row["gene_copy_number"] is None
    assert row["phenotype"] is row["activity_score"] is None
    assert row["confidence_tier"] == "REVIEW"
    assert row["outside_call_eligible"] is False
    assert len(row["candidate_solution_details"]) == 2
    assert len(row["allele_decomposition"]) == 5


def test_definition_solution_and_sample_reconciliation_fail_closed(tmp_path: Path) -> None:
    source = write_output(
        tmp_path / "result.aldy",
        [solution_block(1, "*1/*2", "1.001;2.001", "*1.001, *2.001")],
    )
    wrong_definition = AldyProvenance(
        **{**PROVENANCE.__dict__, "definition_version": "CYP2D6-other-definition"}
    )
    with pytest.raises(PortableEvidenceError, match="definition differs"):
        parse_aldy_results(source, wrong_definition)

    source = write_output(
        tmp_path / "mixed.aldy",
        [
            [
                *solution_block(1, "*1/*2", "1.001;2.001", "*1.001, *2.001")[:-1],
                mutation_row(
                    1,
                    "*1/*2",
                    "1.001;2.001",
                    1,
                    "2.001",
                    101,
                    "C>T",
                    32,
                    "synthetic_effect",
                    "rsSYNTHETIC11",
                    sample="SECOND_SYNTHETIC_SUBJECT",
                ),
            ]
        ],
    )
    with pytest.raises(PortableEvidenceError, match="multiple samples"):
        parse_aldy_results(source, PROVENANCE)

    source = write_output(
        tmp_path / "gap.aldy",
        [solution_block(2, "*1/*2", "1.001;2.001", "*1.001, *2.001")],
    )
    with pytest.raises(PortableEvidenceError, match="not contiguous"):
        parse_aldy_results(source, PROVENANCE)


def test_schema_copy_and_licence_boundaries_fail_closed(tmp_path: Path) -> None:
    source = write_output(
        tmp_path / "result.aldy",
        [solution_block(1, "*1/*2", "1.001;2.001", "*1.001, *2.001")],
    )
    text = source.read_text(encoding="utf-8")
    source.write_text(text.replace("; cpic=", "; unexpected=value; cpic="), encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="unknown or duplicate"):
        parse_aldy_results(source, PROVENANCE)

    source = write_output(
        tmp_path / "copy.aldy",
        [solution_block(1, "*1/*2", "1.001;2.001", "*1.001, *2.001")],
    )
    source.write_text(
        source.read_text(encoding="utf-8").replace("\t1\t2.001\t101\t", "\t2\t2.001\t101\t"),
        encoding="utf-8",
    )
    with pytest.raises(PortableEvidenceError, match="copy indexes are incomplete"):
        parse_aldy_results(source, PROVENANCE)

    source = write_output(
        tmp_path / "licence.aldy",
        [solution_block(1, "*1/*2", "1.001;2.001", "*1.001, *2.001")],
    )
    no_boundary = AldyProvenance(**{**PROVENANCE.__dict__, "licence_gated_caller": False})
    with pytest.raises(PortableEvidenceError, match="licence-gated boundary"):
        parse_aldy_results(source, no_boundary)

    leaking = AldyProvenance(**{**PROVENANCE.__dict__, "raw_artifact_id": f"aldy-{SAMPLE}"})
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        parse_aldy_results(source, leaking)
