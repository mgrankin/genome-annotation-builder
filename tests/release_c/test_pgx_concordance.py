from __future__ import annotations

from copy import deepcopy

import pytest

from genome_annotation_builder.release_c.pgx_concordance import annotate_pgx_concordance
from genome_annotation_builder.release_c.portable_model import (
    PortableEvidenceError,
    make_evidence_id,
)

from .test_portable_model import rows as example_rows


def pgx_row(
    caller: str,
    diplotypes: list[str],
    *,
    source_digit: str,
    selected: str | None,
    definition_digit: str,
) -> dict[str, object]:
    row = deepcopy(example_rows()["complex_pharmacogenes"])
    ambiguous = selected is None
    row.update(
        {
            "caller_name": caller,
            "call_status": "ambiguous" if ambiguous else "called",
            "diplotype": selected,
            "candidate_diplotypes": diplotypes,
            "star_alleles": [],
            "candidate_solution_details": [],
            "allele_decomposition": [],
            "caller_filter": "AMBIGUOUS" if ambiguous else "PASS",
            "ambiguous": ambiguous,
            "outside_call_eligible": not ambiguous,
            "confidence_tier": "REVIEW" if ambiguous else "SUPPORTED",
            "source_alignment_sha256": source_digit * 64,
            "source_record_digest": source_digit * 64,
            "definition_identity_sha256": definition_digit * 64,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
    )
    row["evidence_id"] = make_evidence_id("complex_pharmacogenes", row)
    return row


def test_suballele_projection_and_unordered_sides_are_concordant() -> None:
    cyrius = pgx_row(
        "Cyrius",
        ["*1/*2"],
        source_digit="a",
        selected="*1/*2",
        definition_digit="c",
    )
    aldy = pgx_row(
        "Aldy",
        ["*2.001/*1.018"],
        source_digit="b",
        selected="*2.001/*1.018",
        definition_digit="d",
    )
    annotated, summary = annotate_pgx_concordance(
        [cyrius, aldy], expected_reference_sha256="1" * 64
    )
    assert summary["row_status_counts"] == {"comparable_concordant": 2}
    assert summary["definition_boundary_mismatch_groups"] == 1
    assert summary["majority_vote_applied"] is False
    assert summary["source_rows_removed"] == 0
    assert summary["biological_fields_changed"] is False
    assert all(row["confidence_tier"] == "SUPPORTED" for row in annotated)


def test_tandem_side_disagreement_is_review_without_overwrite() -> None:
    first = pgx_row(
        "Cyrius",
        ["*36+*10/*1"],
        source_digit="a",
        selected="*36+*10/*1",
        definition_digit="c",
    )
    second = pgx_row(
        "Aldy",
        ["*36/*10+*1"],
        source_digit="b",
        selected="*36/*10+*1",
        definition_digit="d",
    )
    annotated, summary = annotate_pgx_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["row_status_counts"] == {"comparable_discordant": 2}
    assert {row["diplotype"] for row in annotated} == {"*36+*10/*1", "*36/*10+*1"}
    assert all(row["confidence_tier"] == "REVIEW" for row in annotated)
    assert all(row["outside_call_eligible"] is False for row in annotated)
    assert all("CROSS_CALLER_COMPLEX_PGX_DISCORDANT" in row["filter_reasons"] for row in annotated)


def test_candidate_overlap_does_not_resolve_ambiguity() -> None:
    called = pgx_row(
        "Cyrius",
        ["*1/*2"],
        source_digit="a",
        selected="*1/*2",
        definition_digit="c",
    )
    ambiguous = pgx_row(
        "Aldy",
        ["*1.018/*2.001", "*1/*4"],
        source_digit="b",
        selected=None,
        definition_digit="d",
    )
    annotated, summary = annotate_pgx_concordance(
        [called, ambiguous], expected_reference_sha256="1" * 64
    )
    by_caller = {row["caller_name"]: row for row in annotated}
    assert summary["row_status_counts"] == {"not_comparable": 2}
    assert summary["candidate_overlap_groups_retained_ambiguous"] == 1
    assert by_caller["Cyrius"]["confidence_tier"] == "SUPPORTED"
    assert by_caller["Aldy"]["confidence_tier"] == "REVIEW"
    assert by_caller["Aldy"]["diplotype"] is None
    assert by_caller["Aldy"]["candidate_diplotypes"] == ["*1.018/*2.001", "*1/*4"]


def test_opaque_structural_suffix_and_duplicate_source_are_fail_safe() -> None:
    first = pgx_row(
        "Cyrius",
        ["*1/*4"],
        source_digit="a",
        selected="*1/*4",
        definition_digit="c",
    )
    opaque = pgx_row(
        "Aldy",
        ["*1/*4NX"],
        source_digit="b",
        selected="*1/*4NX",
        definition_digit="d",
    )
    annotated, summary = annotate_pgx_concordance(
        [first, opaque], expected_reference_sha256="1" * 64
    )
    assert summary["row_status_counts"] == {"not_comparable": 2}
    assert all(row["concordance_status"] == "not_comparable" for row in annotated)

    duplicate = deepcopy(first)
    duplicate["diplotype"] = "*1/*3"
    duplicate["candidate_diplotypes"] = ["*1/*3"]
    duplicate["source_record_digest"] = "e" * 64
    duplicate["evidence_id"] = make_evidence_id("complex_pharmacogenes", duplicate)
    with pytest.raises(PortableEvidenceError, match="duplicate gene rows"):
        annotate_pgx_concordance([first, duplicate, opaque], expected_reference_sha256="1" * 64)
