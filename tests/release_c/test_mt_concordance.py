from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from genome_annotation_builder.release_c.mt_adapter import (
    MitochondrialProvenance,
    parse_filtered_mutect_vcf,
)
from genome_annotation_builder.release_c.mt_concordance import (
    MitochondrialReplicatePolicy,
    annotate_mitochondrial_replicate_concordance,
)
from genome_annotation_builder.release_c.portable_model import (
    PortableEvidenceError,
    make_evidence_id,
)

from .test_mt_adapter import write_vcf
from .test_portable_model import rows as example_rows


def mt_row(
    source_digit: str,
    *,
    allele_fraction: float = 0.2,
    ref_depth: int = 80,
    alt_depth: int = 20,
    depth: int | None = None,
    pos: int = 100,
    resource_digit: str = "3",
    circular_boundary_audited: bool = True,
) -> dict[str, object]:
    row = deepcopy(example_rows()["mitochondrial_heteroplasmy"])
    row.update(
        {
            "event_id": f"mt:chrM:{pos}:A:G",
            "pos": pos,
            "caller_name": "GATK Mutect2 mitochondrial mode",
            "caller_version": "4.6.2.0",
            "allele_fraction": allele_fraction,
            "depth": depth if depth is not None else ref_depth + alt_depth,
            "ref_depth": ref_depth,
            "alt_depth": alt_depth,
            "source_alignment_sha256": source_digit * 64,
            "caller_run_identity_sha256": source_digit * 64,
            "resource_manifest_sha256": resource_digit * 64,
            "source_record_digest": source_digit * 64,
            "raw_artifact_id": f"mt-artifact-{source_digit}",
            "circular_boundary_audited": circular_boundary_audited,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
    )
    row["evidence_id"] = make_evidence_id("mitochondrial_heteroplasmy", row)
    return row


def test_exact_fraction_replicates_preserve_rows_and_biology() -> None:
    annotated, summary = annotate_mitochondrial_replicate_concordance(
        [mt_row("a"), mt_row("b", ref_depth=160, alt_depth=40, depth=200)],
        expected_reference_sha256="1" * 64,
    )
    assert summary["group_relation_counts"] == {"exact_fraction": 1}
    assert summary["row_status_counts"] == {"comparable_concordant": 2}
    assert summary["source_rows_removed"] == 0
    assert summary["allele_fraction_averaging_used"] is False
    assert summary["majority_vote_used"] is False
    assert summary["biological_fields_changed"] is False
    assert {row["allele_fraction"] for row in annotated} == {0.2}
    assert {row["alt_depth"] for row in annotated} == {20, 40}


def test_overlapping_wilson_intervals_are_compatible_without_averaging() -> None:
    annotated, summary = annotate_mitochondrial_replicate_concordance(
        [mt_row("a"), mt_row("b", allele_fraction=0.25, ref_depth=75, alt_depth=25)],
        expected_reference_sha256="1" * 64,
    )
    assert summary["group_relation_counts"] == {"sampling_interval_compatible": 1}
    assert all(row["concordance_status"] == "comparable_concordant" for row in annotated)
    assert {row["allele_fraction"] for row in annotated} == {0.2, 0.25}
    assert all(
        any("Wilson sampling intervals" in value for value in row["limitations"])
        for row in annotated
    )


def test_nonoverlapping_fraction_intervals_are_review_without_call_replacement() -> None:
    first = mt_row("a", allele_fraction=0.05, ref_depth=95, alt_depth=5)
    second = mt_row("b", allele_fraction=0.4, ref_depth=60, alt_depth=40)
    annotated, summary = annotate_mitochondrial_replicate_concordance(
        [first, second], expected_reference_sha256="1" * 64
    )
    assert summary["group_relation_counts"] == {"allele_fraction_discordant": 1}
    assert all(row["concordance_status"] == "comparable_discordant" for row in annotated)
    assert all(row["confidence_tier"] == "REVIEW" for row in annotated)
    assert all(
        "REPLICATE:MT_ALLELE_FRACTION_DISCORDANT" in row["filter_reasons"] for row in annotated
    )
    assert {row["allele_fraction"] for row in annotated} == {0.05, 0.4}


def test_low_depth_and_source_specific_rows_are_not_false_negative_evidence() -> None:
    low_first = mt_row("a", allele_fraction=0.2, ref_depth=8, alt_depth=2, pos=100)
    low_second = mt_row("b", allele_fraction=0.2, ref_depth=8, alt_depth=2, pos=100)
    unique = mt_row("a", pos=200)
    annotated, summary = annotate_mitochondrial_replicate_concordance(
        [low_first, low_second, unique], expected_reference_sha256="1" * 64
    )
    assert summary["group_relation_counts"] == {
        "insufficient_allelic_depth": 1,
        "single_source": 1,
    }
    assert summary["variant_only_absence_used_as_reference"] is False
    by_pos: dict[int, list[dict[str, object]]] = {}
    for row in annotated:
        by_pos.setdefault(int(row["pos"]), []).append(row)
    assert all(
        row["concordance_status"] == "not_comparable" and row["confidence_tier"] == "REVIEW"
        for row in by_pos[100]
    )
    assert by_pos[200][0]["concordance_status"] == "source_specific"
    assert any("not reference evidence" in value for value in by_pos[200][0]["limitations"])


@pytest.mark.parametrize(
    "rows",
    [
        [mt_row("a"), mt_row("b", resource_digit="4")],
        [mt_row("a"), mt_row("b", circular_boundary_audited=False)],
    ],
)
def test_resource_or_circular_boundary_mismatch_is_not_comparable(
    rows: list[dict[str, object]],
) -> None:
    annotated, summary = annotate_mitochondrial_replicate_concordance(
        rows, expected_reference_sha256="1" * 64
    )
    assert summary["group_relation_counts"] == {"caller_resource_or_circular_boundary_mismatch": 1}
    assert all(row["concordance_status"] == "not_comparable" for row in annotated)
    assert all("REPLICATE:MT_NOT_COMPARABLE" in row["filter_reasons"] for row in annotated)


def test_same_alignment_is_not_an_independent_replicate() -> None:
    first = mt_row("a")
    replay = mt_row("a")
    replay["caller_run_identity_sha256"] = "b" * 64
    replay["source_record_digest"] = "b" * 64
    replay["evidence_id"] = make_evidence_id("mitochondrial_heteroplasmy", replay)
    with pytest.raises(PortableEvidenceError, match="two technical alignments"):
        annotate_mitochondrial_replicate_concordance(
            [first, replay], expected_reference_sha256="1" * 64
        )


def test_invalid_policy_and_reference_boundary_fail_closed() -> None:
    with pytest.raises(PortableEvidenceError, match="invalid mitochondrial replicate policy"):
        annotate_mitochondrial_replicate_concordance(
            [mt_row("a"), mt_row("b")],
            policy=MitochondrialReplicatePolicy(minimum_allelic_depth=0),
            expected_reference_sha256="1" * 64,
        )
    second_reference = mt_row("b")
    second_reference["reference_sha256"] = "2" * 64
    second_reference["evidence_id"] = make_evidence_id(
        "mitochondrial_heteroplasmy", second_reference
    )
    with pytest.raises(PortableEvidenceError, match="different references"):
        annotate_mitochondrial_replicate_concordance([mt_row("a"), second_reference])


def test_two_adapter_outputs_flow_into_mitochondrial_concordance(tmp_path: Path) -> None:
    record = (
        "chrM\t100\t.\tA\tG\t.\tPASS\tAS_FilterStatus=PASS;"
        "AS_SB_TABLE=40,40|10,10;AS_UNIQ_ALT_READ_COUNT=20;"
        "MBQ=30,35;MFRL=340,350;MMQ=60,59;MPOS=40;OCM=0;POPAF=2.4;TLOD=20\t"
        "GT:AD:AF:DP:F1R2:F2R1:SB\t0/1:80,20:0.2:100:80,20:0,0:0,0,10,10\n"
    )
    parsed_rows = []
    for ordinal, digit in enumerate(("a", "b"), start=1):
        source = write_vcf(tmp_path / f"source-{ordinal}.vcf", record)
        provenance = MitochondrialProvenance(
            reference_sha256="1" * 64,
            technical_source_role="primary" if ordinal == 1 else "technical_replicate",
            caller_run_identity_sha256=digit * 64,
            resource_manifest_sha256="3" * 64,
            source_alignment_sha256=digit * 64,
            raw_artifact_id=f"mt-source-{ordinal}",
        )
        rows, _ = parse_filtered_mutect_vcf(source, provenance)
        parsed_rows.extend(rows)
    annotated, summary = annotate_mitochondrial_replicate_concordance(
        parsed_rows, expected_reference_sha256="1" * 64
    )
    assert len(annotated) == 2
    assert summary["group_relation_counts"] == {"exact_fraction": 1}
    assert summary["row_status_counts"] == {"comparable_concordant": 2}
