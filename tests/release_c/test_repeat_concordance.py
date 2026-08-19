from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from genome_annotation_builder.release_c.expansionhunter_adapter import (
    ExpansionHunterProvenance,
    parse_expansionhunter_calls,
)
from genome_annotation_builder.release_c.portable_model import (
    PortableEvidenceError,
    make_evidence_id,
)
from genome_annotation_builder.release_c.repeat_concordance import (
    annotate_repeat_replicate_concordance,
)

from .test_expansionhunter_adapter import write_pair
from .test_portable_model import rows as example_rows


def repeat_row(
    source_digit: str,
    *,
    first_count: int = 10,
    first_ci: tuple[int, int | None] = (10, 10),
    second_count: int = 20,
    second_ci: tuple[int, int | None] = (20, 20),
    variant_id: str = "SYNTHETIC_REPEAT_VARIANT",
) -> dict[str, object]:
    row = deepcopy(example_rows()["repeat_expansions"])
    row.update(
        {
            "caller_name": "ExpansionHunter",
            "caller_version": "5.0.0",
            "variant_id": variant_id,
            "event_id": f"repeat:chr3:100-120:{variant_id}",
            "source_alignment_sha256": source_digit * 64,
            "caller_run_identity_sha256": source_digit * 64,
            "source_record_digest": source_digit * 64,
            "raw_artifact_id": f"repeat-artifact-{source_digit}",
            "allele1_repeat_count": first_count,
            "allele1_ci_lower": first_ci[0],
            "allele1_ci_upper": first_ci[1],
            "allele1_count_is_lower_bound": first_ci[1] is None,
            "allele2_repeat_count": second_count,
            "allele2_ci_lower": second_ci[0],
            "allele2_ci_upper": second_ci[1],
            "allele2_count_is_lower_bound": second_ci[1] is None,
            "short_read_resolvable": first_ci == (first_count, first_count)
            and second_ci == (second_count, second_count),
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
    )
    row["evidence_id"] = make_evidence_id("repeat_expansions", row)
    return row


def test_exact_replicate_agreement_preserves_rows_and_biology() -> None:
    annotated, summary = annotate_repeat_replicate_concordance(
        [repeat_row("a"), repeat_row("b")], expected_reference_sha256="1" * 64
    )
    assert summary["input_rows"] == summary["output_rows"] == 2
    assert summary["group_relation_counts"] == {"exact": 1}
    assert summary["row_status_counts"] == {"comparable_concordant": 2}
    assert summary["biological_fields_changed"] is False
    assert summary["source_rows_removed"] == 0
    assert summary["majority_vote_used"] is False
    assert all(row["confidence_tier"] == "SUPPORTED" for row in annotated)


def test_unordered_interval_overlap_is_concordant_without_exact_claim() -> None:
    first = repeat_row("a", first_count=10, first_ci=(9, 12), second_count=20, second_ci=(18, 22))
    second = repeat_row("b", first_count=21, first_ci=(19, 23), second_count=11, second_ci=(10, 13))
    annotated, summary = annotate_repeat_replicate_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["group_relation_counts"] == {"interval_compatible": 1}
    assert all(row["concordance_status"] == "comparable_concordant" for row in annotated)
    assert all(
        any("overlapping repeat confidence intervals" in value for value in row["limitations"])
        for row in annotated
    )


def test_disjoint_repeat_intervals_are_review_without_voting() -> None:
    first = repeat_row("a")
    second = repeat_row("b", first_count=30, first_ci=(29, 31))
    annotated, summary = annotate_repeat_replicate_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["group_relation_counts"] == {"discordant": 1}
    assert all(row["concordance_status"] == "comparable_discordant" for row in annotated)
    assert all(row["confidence_tier"] == "REVIEW" for row in annotated)
    assert all("REPLICATE:REPEAT_GENOTYPE_DISCORDANT" in row["filter_reasons"] for row in annotated)
    assert {row["allele1_repeat_count"] for row in annotated} == {10, 30}


def test_open_ended_and_source_specific_are_not_false_negatives() -> None:
    first = repeat_row("a", second_count=40, second_ci=(40, None))
    second = repeat_row("b", second_count=45, second_ci=(45, None))
    unique = repeat_row("a", variant_id="UNIQUE_REPEAT")
    annotated, summary = annotate_repeat_replicate_concordance(
        [first, second, unique], expected_reference_sha256="1" * 64
    )
    assert summary["group_relation_counts"] == {
        "open_ended_not_comparable": 1,
        "single_source": 1,
    }
    by_variant: dict[str, list[dict[str, object]]] = {}
    for row in annotated:
        by_variant.setdefault(str(row["variant_id"]), []).append(row)
    assert all(
        row["concordance_status"] == "not_comparable" and row["confidence_tier"] == "REVIEW"
        for row in by_variant["SYNTHETIC_REPEAT_VARIANT"]
    )
    assert by_variant["UNIQUE_REPEAT"][0]["concordance_status"] == "source_specific"
    assert any(
        "absence is not refutation" in value
        for value in by_variant["UNIQUE_REPEAT"][0]["limitations"]
    )


def test_same_alignment_is_not_an_independent_replicate() -> None:
    first = repeat_row("a")
    replay = repeat_row("a")
    replay["caller_run_identity_sha256"] = "b" * 64
    replay["source_record_digest"] = "b" * 64
    replay["evidence_id"] = make_evidence_id("repeat_expansions", replay)
    with pytest.raises(PortableEvidenceError, match="two technical alignments"):
        annotate_repeat_replicate_concordance([first, replay], expected_reference_sha256="1" * 64)


def test_two_adapter_outputs_flow_into_replicate_concordance(tmp_path: Path) -> None:
    parsed_rows = []
    for ordinal, digit in enumerate(("a", "b"), start=1):
        source_root = tmp_path / f"source-{ordinal}"
        source_root.mkdir()
        source_json, source_vcf = write_pair(source_root)
        provenance = ExpansionHunterProvenance(
            reference_sha256="1" * 64,
            technical_source_role="primary" if ordinal == 1 else "technical_replicate",
            caller_run_identity_sha256=digit * 64,
            resource_manifest_sha256="3" * 64,
            source_alignment_sha256=digit * 64,
            raw_artifact_id=f"repeat-source-{ordinal}",
        )
        rows, _ = parse_expansionhunter_calls(source_json, source_vcf, provenance)
        parsed_rows.extend(rows)
    annotated, summary = annotate_repeat_replicate_concordance(
        parsed_rows, expected_reference_sha256="1" * 64
    )
    assert len(annotated) == 2
    assert summary["group_relation_counts"] == {"exact": 1}
    assert summary["row_status_counts"] == {"comparable_concordant": 2}
