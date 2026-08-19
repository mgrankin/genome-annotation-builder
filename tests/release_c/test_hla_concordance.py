from __future__ import annotations

from copy import deepcopy

from genome_annotation_builder.release_c.hla_concordance import annotate_hla_concordance
from genome_annotation_builder.release_c.portable_model import make_evidence_id

from .test_portable_model import rows as example_rows


def hla_row(
    caller: str,
    *,
    full: tuple[str, str],
    g_group: tuple[str | None, str | None],
    two_field: tuple[str, str],
    source_digit: str,
    database_age_limitation: bool,
) -> dict[str, object]:
    row = deepcopy(example_rows()["hla_typing"])
    gene = full[0].split("*", maxsplit=1)[0]
    locus = f"HLA-{gene}"
    row.update(
        {
            "locus": locus,
            "event_id": f"hla:{locus}",
            "caller_name": caller,
            "allele1": full[0],
            "allele2": full[1],
            "allele1_g_group": g_group[0],
            "allele2_g_group": g_group[1],
            "allele1_two_field": two_field[0],
            "allele2_two_field": two_field[1],
            "allele1_candidates": [full[0]],
            "allele2_candidates": [full[1]],
            "source_alignment_sha256": source_digit * 64,
            "source_record_digest": source_digit * 64,
            "database_age_limitation": database_age_limitation,
            "caller_database_identity_sha256": source_digit * 64,
            "nomenclature_identity_sha256": source_digit * 64,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
    )
    row["evidence_id"] = make_evidence_id("hla_typing", row)
    return row


def test_g_group_comparison_handles_different_full_two_field_representatives() -> None:
    current = hla_row(
        "SpecHLA",
        full=("DQA1*05:05:17:01", "DQA1*05:01:08:02"),
        g_group=("DQA1*05:01:01G", "DQA1*05:01:01G"),
        two_field=("DQA1*05:05", "DQA1*05:01"),
        source_digit="a",
        database_age_limitation=False,
    )
    historical = hla_row(
        "HLA*LA",
        full=("DQA1*05:01:01G", "DQA1*05:01:01G"),
        g_group=("DQA1*05:01:01G", "DQA1*05:01:01G"),
        two_field=("DQA1*05:01", "DQA1*05:01"),
        source_digit="b",
        database_age_limitation=True,
    )
    annotated, summary = annotate_hla_concordance(
        [current, historical], expected_reference_sha256="1" * 64
    )
    assert summary["row_status_counts"] == {"comparable_concordant": 2}
    assert summary["row_boundary_counts"] == {"g-group": 2}
    assert summary["biological_fields_changed"] is False
    assert all(row["concordance_status"] == "comparable_concordant" for row in annotated)
    assert all(row["confidence_tier"] == "SUPPORTED" for row in annotated)


def test_historical_disagreement_cannot_downgrade_current_call() -> None:
    current = hla_row(
        "SpecHLA",
        full=("A*01:01:01:01", "A*02:01:01:01"),
        g_group=("A*01:01:01G", "A*02:01:01G"),
        two_field=("A*01:01", "A*02:01"),
        source_digit="a",
        database_age_limitation=False,
    )
    historical = hla_row(
        "HLA*LA",
        full=("A*03:01:01G", "A*24:02:01G"),
        g_group=("A*03:01:01G", "A*24:02:01G"),
        two_field=("A*03:01", "A*24:02"),
        source_digit="b",
        database_age_limitation=True,
    )
    annotated, summary = annotate_hla_concordance(
        [current, historical], expected_reference_sha256="1" * 64
    )
    by_caller = {row["caller_name"]: row for row in annotated}
    assert summary["row_status_counts"] == {"comparable_discordant": 2}
    assert summary["current_database_rows_protected_from_historical_override"] == 1
    assert by_caller["SpecHLA"]["confidence_tier"] == "SUPPORTED"
    assert by_caller["HLA*LA"]["confidence_tier"] == "REVIEW"
    assert by_caller["SpecHLA"]["allele1"] == "A*01:01:01:01"
    assert by_caller["HLA*LA"]["allele1"] == "A*03:01:01G"


def test_current_source_disagreement_is_review_and_incomplete_is_not_comparable() -> None:
    first = hla_row(
        "CallerOne",
        full=("A*01:01", "A*02:01"),
        g_group=(None, None),
        two_field=("A*01:01", "A*02:01"),
        source_digit="a",
        database_age_limitation=False,
    )
    second = hla_row(
        "CallerTwo",
        full=("A*01:01", "A*03:01"),
        g_group=(None, None),
        two_field=("A*01:01", "A*03:01"),
        source_digit="b",
        database_age_limitation=False,
    )
    annotated, summary = annotate_hla_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["row_boundary_counts"] == {"two-field": 2}
    assert all(row["confidence_tier"] == "REVIEW" for row in annotated)

    incomplete = deepcopy(second)
    incomplete.update(
        {
            "allele2": None,
            "allele2_two_field": None,
            "allele2_candidates": [],
            "call_status": "partial_call",
            "confidence_tier": "REVIEW",
        }
    )
    incomplete["evidence_id"] = make_evidence_id("hla_typing", incomplete)
    annotated, summary = annotate_hla_concordance(
        [first, incomplete], expected_reference_sha256="1" * 64
    )
    assert summary["row_status_counts"] == {"not_comparable": 2}
    assert all(row["concordance_status"] == "not_comparable" for row in annotated)
