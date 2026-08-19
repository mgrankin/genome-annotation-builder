from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from genome_annotation_builder.arrow_schemas import gwas_catalog_resource_schema
from genome_annotation_builder.gwas_parser import SortedGwasCatalogLookup
from genome_annotation_builder.gwas_resources import EXPECTED_HEADER, _normalized_row


def _source_row(**updates: str) -> dict[str, str]:
    row = {field: "" for field in EXPECTED_HEADER}
    row.update(
        {
            "PUBMEDID": "12345678",
            "DISEASE/TRAIT": "Synthetic quantitative trait",
            "CHR_ID": "1",
            "CHR_POS": "100",
            "STRONGEST SNP-RISK ALLELE": "rs123-C",
            "SNPS": "rs123",
            "P-VALUE": "1E-8",
            "OR or BETA": "0.5",
            "95% CI (TEXT)": "[0.3-0.7] unit decrease",
            "MAPPED_TRAIT": "synthetic trait measurement",
            "MAPPED_TRAIT_URI": "http://example.org/EFO_1",
            "STUDY ACCESSION": "GCST000001",
        }
    )
    row.update(updates)
    return row


def test_gwas_beta_direction_and_ci_are_oriented_without_losing_raw_text() -> None:
    row, reason = _normalized_row(_source_row(), 2)
    assert reason is None
    assert row is not None
    assert row["effect_type"] == "beta"
    assert row["effect_size"] == -0.5
    assert row["ci_lower"] == -0.7
    assert row["ci_upper"] == -0.3
    assert row["raw_effect_text"] == "0.5"
    assert row["raw_ci_text"] == "[0.3-0.7] unit decrease"

    signed, reason = _normalized_row(
        _source_row(**{"95% CI (TEXT)": "[-0.7--0.3] unit decrease"}), 3
    )
    assert reason is None
    assert signed is not None
    assert signed["effect_size"] == -0.5
    assert signed["ci_lower"] == -0.7
    assert signed["ci_upper"] == -0.3


def test_gwas_unknown_or_complex_effect_allele_is_not_guessed() -> None:
    row, reason = _normalized_row(_source_row(**{"STRONGEST SNP-RISK ALLELE": "rs123-?"}), 2)
    assert row is None
    assert reason == "unknown_or_complex_risk_allele"


def test_gwas_plain_conjunction_or_is_not_misread_as_odds_ratio() -> None:
    row, reason = _normalized_row(
        _source_row(
            **{
                "95% CI (TEXT)": "[0.3-0.7] adjusted for age or sex",
            }
        ),
        2,
    )
    assert reason is None
    assert row is not None
    assert row["effect_size"] == 0.5
    assert row["effect_type"] is None
    assert row["effect_type_basis"] is None


def _resource_row(
    *,
    effect_allele: str,
    p_value: float,
    effect_size: float | None,
    effect_type: str | None,
    source_row: int,
) -> dict[str, object]:
    return {
        "chrom": "chr1",
        "pos": 100,
        "source_row_number": source_row,
        "risk_rsid": "rs123",
        "rsid_aliases": ["rs123"],
        "effect_allele": effect_allele,
        "reported_trait": f"Trait {source_row}",
        "mapped_traits": [f"mapped trait {source_row}"],
        "mapped_trait_uris": [f"http://example.org/EFO_{source_row}"],
        "p_value": p_value,
        "effect_size": effect_size,
        "effect_type": effect_type,
        "effect_type_basis": "synthetic explicit marker" if effect_type else None,
        "ci_lower": None,
        "ci_upper": None,
        "sample_size": 1000,
        "pmid": str(10000000 + source_row),
        "study_accession": f"GCST{source_row:06d}",
        "raw_effect_text": str(effect_size) if effect_size is not None else None,
        "raw_ci_text": "unit increase" if effect_type else None,
        "source_record_digest": f"{source_row:064x}",
    }


def test_gwas_merge_join_requires_coordinate_rsid_and_exact_effect_allele(
    tmp_path: Path,
) -> None:
    resource = tmp_path / "gwas.parquet"
    rows = [
        _resource_row(
            effect_allele="A",
            p_value=1e-9,
            effect_size=0.25,
            effect_type="beta",
            source_row=2,
        ),
        _resource_row(
            effect_allele="C",
            p_value=1e-12,
            effect_size=1.2,
            effect_type=None,
            source_row=3,
        ),
        _resource_row(
            effect_allele="G",
            p_value=1e-20,
            effect_size=0.9,
            effect_type="beta",
            source_row=4,
        ),
    ]
    pq.write_table(
        pa.Table.from_pylist(rows, schema=gwas_catalog_resource_schema()),
        resource,
        compression="zstd",
        version="2.6",
    )

    with SortedGwasCatalogLookup(resource, "2026-08-03") as lookup:
        core, evidence, normative = lookup.lookup(
            chrom="chr1",
            pos=100,
            ref="A",
            alt="C",
            rsid="rs123",
            gt=(0, 1),
        )
    assert core["is_gwas_hit"] is True
    assert core["effect_allele"] == "A"
    assert core["effect_size"] == 0.25
    assert core["p_value"] == 1e-9
    assert len(evidence) == 2
    assert len(normative) == 1
    assert {row["effect_allele_dosage"] for row in evidence} == {1}
    assert {row["effect_allele"] for row in evidence} == {"A", "C"}
    assert any(row["omission_reason"] for row in evidence)

    with SortedGwasCatalogLookup(resource, "2026-08-03") as lookup:
        core, evidence, normative = lookup.lookup(
            chrom="chr1",
            pos=100,
            ref="A",
            alt="C",
            rsid="rs999",
            gt=(1, 1),
        )
    assert core["is_gwas_hit"] is False
    assert evidence == []
    assert normative == []


def test_gwas_semantic_duplicate_source_rows_keep_unique_stable_ids(tmp_path: Path) -> None:
    resource = tmp_path / "gwas-duplicates.parquet"
    first = _resource_row(
        effect_allele="C",
        p_value=1e-12,
        effect_size=1.2,
        effect_type="odds_ratio",
        source_row=10,
    )
    second = {**first, "source_row_number": 11}
    assert first["source_record_digest"] == second["source_record_digest"]
    pq.write_table(
        pa.Table.from_pylist([first, second], schema=gwas_catalog_resource_schema()),
        resource,
        compression="zstd",
        version="2.6",
    )

    def project() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
        with SortedGwasCatalogLookup(resource, "2026-08-03") as lookup:
            _core, evidence, normative = lookup.lookup(
                chrom="chr1",
                pos=100,
                ref="A",
                alt="C",
                rsid="rs123",
                gt=(0, 1),
            )
        return evidence, normative

    evidence, normative = project()
    assert len(evidence) == 2
    assert len(normative) == 2
    assert {row["source_row_number"] for row in evidence} == {10, 11}
    assert len({row["association_id"] for row in evidence}) == 2
    assert {row["association_id"] for row in normative} == {
        row["association_id"] for row in evidence
    }
    replay_evidence, replay_normative = project()
    assert replay_evidence == evidence
    assert replay_normative == normative
