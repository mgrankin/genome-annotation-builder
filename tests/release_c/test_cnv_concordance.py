from __future__ import annotations

from copy import deepcopy

from genome_annotation_builder.release_c.cnv_concordance import annotate_cnv_concordance
from genome_annotation_builder.release_c.portable_model import make_evidence_id

from .test_portable_model import rows as example_rows


def cnv_row(
    model: str,
    *,
    caller: str = "CNVpytor",
    source_digit: str,
    run_digit: str,
    bin_size: int,
    interval: tuple[int, int],
    call_type: str,
    copy_number: float,
) -> dict[str, object]:
    row = deepcopy(example_rows()["copy_number"])
    start, end = interval
    row.update(
        {
            "event_id": f"cnv:chr2:{start}:{end}:{model}:{bin_size}",
            "caller_name": caller,
            "caller_model": model,
            "caller_run_identity_sha256": run_digit * 64,
            "source_alignment_sha256": source_digit * 64,
            "source_record_digest": run_digit * 64,
            "start_pos": start,
            "end_pos": end,
            "call_type": call_type,
            "copy_number": copy_number,
            "read_depth_copy_number": copy_number,
            "bin_size": bin_size,
            "filter_reasons": [],
            "confidence_tier": "SUPPORTED",
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
    )
    row["evidence_id"] = make_evidence_id("copy_number", row)
    return row


def test_same_caller_cross_configuration_agreement_is_not_independent_replication() -> None:
    first = cnv_row(
        "mean_shift_read_depth",
        source_digit="a",
        run_digit="b",
        bin_size=1000,
        interval=(100_000, 200_000),
        call_type="deletion",
        copy_number=1.1,
    )
    second = cnv_row(
        "combined_rd_baf_control",
        source_digit="a",
        run_digit="c",
        bin_size=10_000,
        interval=(105_000, 205_000),
        call_type="deletion",
        copy_number=1.0,
    )
    annotated, summary = annotate_cnv_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["callers"] == 1
    assert summary["configurations"] == 2
    assert summary["same_caller_cross_configuration_concordant_components"] == 1
    assert summary["cross_caller_concordant_components"] == 0
    assert summary["row_status_counts"] == {"comparable_concordant": 2}
    assert summary["majority_vote_applied"] is False
    assert summary["source_rows_removed"] == 0
    assert summary["biological_fields_changed"] is False
    assert all(row["confidence_tier"] == "SUPPORTED" for row in annotated)
    assert all(
        any("not independent caller replication" in value for value in row["limitations"])
        for row in annotated
    )


def test_cross_caller_agreement_is_counted_separately() -> None:
    first = cnv_row(
        "depth",
        caller="CallerOne",
        source_digit="a",
        run_digit="b",
        bin_size=1000,
        interval=(100_000, 200_000),
        call_type="duplication",
        copy_number=3.0,
    )
    second = cnv_row(
        "breakpoint-depth",
        caller="CallerTwo",
        source_digit="a",
        run_digit="c",
        bin_size=1000,
        interval=(100_100, 199_900),
        call_type="duplication",
        copy_number=3.2,
    )
    _, summary = annotate_cnv_concordance([first, second], expected_reference_sha256="1" * 64)
    assert summary["callers"] == 2
    assert summary["cross_caller_concordant_components"] == 1


def test_overlapping_type_or_copy_number_disagreement_is_review() -> None:
    deletion = cnv_row(
        "depth",
        source_digit="a",
        run_digit="b",
        bin_size=1000,
        interval=(100_000, 200_000),
        call_type="deletion",
        copy_number=1.0,
    )
    duplication = cnv_row(
        "combined",
        source_digit="a",
        run_digit="c",
        bin_size=10_000,
        interval=(100_000, 200_000),
        call_type="duplication",
        copy_number=3.0,
    )
    annotated, summary = annotate_cnv_concordance(
        [deletion, duplication], expected_reference_sha256="1" * 64
    )
    assert summary["row_status_counts"] == {"comparable_discordant": 2}
    assert all(row["confidence_tier"] == "REVIEW" for row in annotated)
    assert all("CROSS_CONFIGURATION_CNV_DISCORDANT" in row["filter_reasons"] for row in annotated)


def test_nonoverlapping_calls_remain_source_specific_not_refuted() -> None:
    first = cnv_row(
        "depth",
        source_digit="a",
        run_digit="b",
        bin_size=1000,
        interval=(100_000, 200_000),
        call_type="deletion",
        copy_number=1.0,
    )
    second = cnv_row(
        "combined",
        source_digit="a",
        run_digit="c",
        bin_size=10_000,
        interval=(500_000, 600_000),
        call_type="deletion",
        copy_number=1.0,
    )
    annotated, summary = annotate_cnv_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["row_status_counts"] == {"source_specific": 2}
    assert all(row["confidence_tier"] == "SUPPORTED" for row in annotated)
    assert all(
        any("absence is not refutation" in value for value in row["limitations"])
        for row in annotated
    )
