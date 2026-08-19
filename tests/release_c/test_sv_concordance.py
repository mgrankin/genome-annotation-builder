from __future__ import annotations

from copy import deepcopy

import pytest

from genome_annotation_builder.release_c.portable_model import (
    PortableEvidenceError,
    make_evidence_id,
)
from genome_annotation_builder.release_c.sv_concordance import (
    ConcordancePolicy,
    _biological_digest,
    _build_events,
    _events_match,
    _matching_event_pairs,
    _source,
    annotate_sv_concordance,
)

from .test_portable_model import rows as example_rows


def sv_row(
    caller: str,
    record_id: str,
    *,
    start: int,
    end: int,
    genotype: list[int] | None = None,
    orientation: str = "3to5",
    sv_type: str = "DEL",
    source_digit: str,
) -> dict[str, object]:
    row = deepcopy(example_rows()["structural_variants"])
    row.update(
        {
            "caller_name": caller,
            "caller_record_id": record_id,
            "chrom1": "chr1",
            "pos1": start,
            "chrom2": "chr1",
            "pos2": end,
            "sv_type": sv_type,
            "orientation": orientation,
            "length": end - start if sv_type != "DEL" else start - end,
            "genotype": genotype if genotype is not None else [0, 1],
            "call_status": "non_reference",
            "source_alignment_sha256": source_digit * 64,
            "source_record_digest": record_id.encode().hex().ljust(64, "0")[:64],
            "event_id": f"sv:{sv_type}:chr1:{start}:chr1:{end}:{orientation}",
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
    )
    row["evidence_id"] = make_evidence_id("structural_variants", row)
    return row


def test_type_aware_matching_preserves_rows_and_biological_fields() -> None:
    rows = [
        sv_row("DELLY", "d1", start=100, end=300, source_digit="a"),
        sv_row("Manta", "m1", start=102, end=302, source_digit="b"),
        sv_row("DELLY", "d2", start=1000, end=1200, source_digit="a"),
        sv_row("Manta", "m2", start=5000, end=5200, source_digit="b"),
    ]
    annotated, summary = annotate_sv_concordance(
        rows,
        policy=ConcordancePolicy(max_breakpoint_distance=10),
        expected_reference_sha256="1" * 64,
    )
    assert summary["input_rows"] == summary["output_rows"] == 4
    assert summary["candidate_edges"] == 1
    assert summary["biological_fields_changed"] is False
    by_id = {row["caller_record_id"]: row for row in annotated}
    assert by_id["d1"]["concordance_status"] == "comparable_concordant"
    assert by_id["m1"]["concordance_group_id"] == by_id["d1"]["concordance_group_id"]
    assert by_id["d2"]["concordance_status"] == "source_specific"
    assert by_id["m2"]["concordance_status"] == "source_specific"
    assert [by_id[key]["genotype"] for key in ("d1", "m1", "d2", "m2")] == [
        [0, 1],
        [0, 1],
        [0, 1],
        [0, 1],
    ]


def test_genotype_disagreement_is_review_not_intersection_drop() -> None:
    first = sv_row("DELLY", "d1", start=100, end=300, source_digit="a")
    second = sv_row("Manta", "m1", start=100, end=300, genotype=[1, 1], source_digit="b")
    annotated, summary = annotate_sv_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["output_rows"] == 2
    assert summary["row_status_counts"] == {"comparable_discordant": 2}
    assert all(row["confidence_tier"] == "REVIEW" for row in annotated)
    assert {tuple(row["genotype"]) for row in annotated} == {(0, 1), (1, 1)}


def test_delly_and_manta_bnd_orientation_notation_is_comparable() -> None:
    first = sv_row(
        "DELLY",
        "d-bnd",
        start=100,
        end=500,
        orientation="5to5",
        sv_type="BND",
        source_digit="a",
    )
    second = sv_row(
        "Manta",
        "m-bnd",
        start=102,
        end=498,
        orientation="--",
        sv_type="BND",
        source_digit="b",
    )
    annotated, summary = annotate_sv_concordance(
        [first, second],
        policy=ConcordancePolicy(max_breakpoint_distance=10),
        expected_reference_sha256="1" * 64,
    )
    assert summary["candidate_edges"] == 1
    assert all(row["concordance_status"] == "comparable_concordant" for row in annotated)


def test_delly_and_manta_insertion_orientation_notation_is_comparable() -> None:
    first = sv_row(
        "DELLY",
        "d-ins",
        start=100,
        end=140,
        orientation="NtoN",
        sv_type="INS",
        source_digit="a",
    )
    second = sv_row(
        "Manta",
        "m-ins",
        start=102,
        end=142,
        orientation="not_applicable",
        sv_type="INS",
        source_digit="b",
    )
    annotated, summary = annotate_sv_concordance(
        [first, second],
        policy=ConcordancePolicy(max_breakpoint_distance=10),
        expected_reference_sha256="1" * 64,
    )
    assert summary["candidate_edges"] == 1
    assert all(row["concordance_status"] == "comparable_concordant" for row in annotated)


def test_unknown_orientation_remains_fail_closed() -> None:
    rows = [
        sv_row(
            "DELLY",
            "d-unknown",
            start=100,
            end=140,
            orientation="unknown",
            sv_type="INS",
            source_digit="a",
        ),
        sv_row(
            "Manta",
            "m-ins",
            start=102,
            end=142,
            orientation="not_applicable",
            sv_type="INS",
            source_digit="b",
        ),
    ]
    with pytest.raises(PortableEvidenceError, match="unsupported SV orientation"):
        annotate_sv_concordance(rows, expected_reference_sha256="1" * 64)


def test_sweep_matches_the_reference_pairwise_definition() -> None:
    rows = [
        sv_row("DELLY", "d1", start=100, end=300, source_digit="a"),
        sv_row("Manta", "m1", start=102, end=302, source_digit="b"),
        sv_row("DELLY", "d2", start=700, end=900, source_digit="a"),
        sv_row("Manta", "m2", start=1190, end=1390, source_digit="b"),
        sv_row(
            "DELLY",
            "d-ins",
            start=2000,
            end=2040,
            orientation="NtoN",
            sv_type="INS",
            source_digit="a",
        ),
        sv_row(
            "Manta",
            "m-ins",
            start=2004,
            end=2044,
            orientation="not_applicable",
            sv_type="INS",
            source_digit="b",
        ),
    ]
    rows[2]["cipos"] = [-10, 20]
    rows[2]["ciend"] = [-5, 10]
    rows[2]["evidence_id"] = make_evidence_id("structural_variants", rows[2])
    policy = ConcordancePolicy(max_breakpoint_distance=500)
    events = _build_events(rows)
    expected = {
        (first_index, second_index)
        for first_index, first in enumerate(events)
        for second_index, second in enumerate(events[first_index + 1 :], first_index + 1)
        if _source(first) != _source(second) and _events_match(first, second, policy)
    }
    assert set(_matching_event_pairs(events, policy)) == expected


def test_reversed_breakpoint_interval_is_rejected() -> None:
    first = sv_row("DELLY", "d1", start=100, end=300, source_digit="a")
    first["cipos"] = [10, -10]
    first["evidence_id"] = make_evidence_id("structural_variants", first)
    second = sv_row("Manta", "m1", start=100, end=300, source_digit="b")
    with pytest.raises(PortableEvidenceError, match="interval bounds are reversed"):
        annotate_sv_concordance([first, second], expected_reference_sha256="1" * 64)


def test_colliding_event_id_with_overlapping_ci_is_retained_conservatively() -> None:
    first = sv_row(
        "Manta",
        "m-ins-1",
        start=100,
        end=100,
        orientation="not_applicable",
        sv_type="INS",
        source_digit="a",
    )
    second = sv_row(
        "Manta",
        "m-ins-2",
        start=100,
        end=100,
        orientation="not_applicable",
        sv_type="INS",
        source_digit="a",
    )
    second["cipos"] = [0, 1]
    second["event_id"] = first["event_id"]
    second["evidence_id"] = make_evidence_id("structural_variants", second)
    audit = sv_row(
        "DELLY",
        "d-ins",
        start=100,
        end=100,
        orientation="NtoN",
        sv_type="INS",
        source_digit="b",
    )

    annotated, summary = annotate_sv_concordance(
        [first, second, audit], expected_reference_sha256="1" * 64
    )

    assert summary["input_rows"] == summary["output_rows"] == 3
    assert summary["events"] == 3
    assert summary["candidate_edges"] == 2
    assert summary["event_id_collision_groups_partitioned"] == 1
    assert summary["event_id_collision_rows_partitioned"] == 2
    assert summary["row_status_counts"] == {"comparable_discordant": 3}
    by_id = {row["caller_record_id"]: row for row in annotated}
    for record_id in ("m-ins-1", "m-ins-2"):
        assert by_id[record_id]["confidence_tier"] == "REVIEW"
        assert any(
            "retained as distinct source records" in limitation
            for limitation in by_id[record_id]["limitations"]
        )


def test_colliding_event_id_with_nonoverlapping_ci_remains_fail_closed() -> None:
    first = sv_row(
        "Manta",
        "m-ins-1",
        start=100,
        end=100,
        orientation="not_applicable",
        sv_type="INS",
        source_digit="a",
    )
    second = sv_row(
        "Manta",
        "m-ins-2",
        start=100,
        end=100,
        orientation="not_applicable",
        sv_type="INS",
        source_digit="a",
    )
    second["cipos"] = [2, 3]
    second["event_id"] = first["event_id"]
    second["evidence_id"] = make_evidence_id("structural_variants", second)
    audit = sv_row(
        "DELLY",
        "d-ins",
        start=100,
        end=100,
        orientation="NtoN",
        sv_type="INS",
        source_digit="b",
    )

    with pytest.raises(PortableEvidenceError, match="inconsistent breakend records"):
        annotate_sv_concordance([first, second, audit], expected_reference_sha256="1" * 64)


def test_biological_digest_is_order_independent_and_content_sensitive() -> None:
    first = sv_row("DELLY", "d1", start=100, end=300, source_digit="a")
    second = sv_row("DELLY", "d2", start=100, end=300, source_digit="a")
    assert _biological_digest([first, second]) == _biological_digest([second, first])

    changed = deepcopy(second)
    changed["genotype"] = [1, 1]
    assert _biological_digest([first, second]) != _biological_digest([first, changed])
