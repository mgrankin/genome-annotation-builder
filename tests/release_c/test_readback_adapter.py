from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq
import pysam
import pytest

from genome_annotation_builder.release_c import readback_capture
from genome_annotation_builder.release_c.portable_model import (
    PortableEvidenceError,
    make_evidence_id,
    schema_for,
    sha256_file,
)
from genome_annotation_builder.release_c.readback_adapter import (
    ReadbackPolicy,
    ReadbackProvenance,
    materialize_readback_capture,
    materialize_readback_targets,
    parse_readback_capture,
    parse_readback_targets,
)
from genome_annotation_builder.release_c.readback_capture import (
    build_capture,
    write_capture,
)
from genome_annotation_builder.release_c.readback_concordance import (
    ReadbackReplicatePolicy,
    annotate_readback_replicate_concordance,
)

SYNTHETIC_SAMPLE = "SYNTHETIC_SUBJECT"


def _read(
    name: str,
    start: int,
    *,
    sequence: str,
    cigar: str = "50M",
    reverse: bool = False,
    read2: bool = False,
    mapping_quality: int = 60,
    duplicate: bool = False,
) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment()
    record.query_name = name
    record.query_sequence = sequence
    record.flag = 128 if read2 else 64
    if reverse:
        record.flag |= 16
    if duplicate:
        record.flag |= 1024
    record.reference_id = 0
    record.reference_start = start
    record.mapping_quality = mapping_quality
    record.cigarstring = cigar
    record.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    record.set_tag("RG", "synthetic-rg")
    return record


def _substitution_sequence(base: str, query_position: int = 20) -> str:
    sequence = list("A" * 50)
    sequence[query_position] = base
    return "".join(sequence)


def test_unpaired_read_orientation_has_a_matching_metrics_bucket() -> None:
    record = pysam.AlignedSegment()
    record.flag = 0
    orientation = readback_capture._orientation(record)
    metrics = readback_capture._empty_metrics(readback_capture.DEFAULT_CAPTURE_POLICY)
    assert orientation == "unclassified_orientation"
    assert f"{orientation}_ref_reads" in metrics


def build_inputs(tmp_path: Path) -> tuple[Path, Path, Path, ReadbackProvenance]:
    reference = tmp_path / "synthetic.fa"
    reference.write_text(">chr1\n" + "A" * 500 + "\n", encoding="ascii")
    pysam.faidx(str(reference))

    alignment = tmp_path / "synthetic.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 500}],
        "RG": [{"ID": "synthetic-rg", "SM": SYNTHETIC_SAMPLE}],
    }
    records: list[pysam.AlignedSegment] = []
    for index in range(4):
        records.append(_read(f"snv-ref-{index}", 80, sequence=_substitution_sequence("A")))
    for index in range(4):
        records.append(
            _read(
                f"snv-alt-{index}",
                80,
                sequence=_substitution_sequence("T"),
                reverse=index >= 2,
            )
        )
    records.append(_read("snv-other", 80, sequence=_substitution_sequence("C")))
    records.append(
        _read(
            "snv-low-mapq",
            80,
            sequence=_substitution_sequence("T"),
            mapping_quality=5,
        )
    )
    records.append(
        _read(
            "snv-duplicate",
            80,
            sequence=_substitution_sequence("T"),
            duplicate=True,
        )
    )
    records.extend(
        [
            _read("snv-overlap-match", 80, sequence=_substitution_sequence("T")),
            _read(
                "snv-overlap-match",
                80,
                sequence=_substitution_sequence("T"),
                read2=True,
                reverse=True,
            ),
            _read("snv-overlap-discordant", 80, sequence=_substitution_sequence("T")),
            _read(
                "snv-overlap-discordant",
                80,
                sequence=_substitution_sequence("A"),
                read2=True,
                reverse=True,
            ),
        ]
    )

    insertion_alt = list("A" * 50)
    insertion_alt[21:23] = list("TG")
    for index in range(4):
        records.append(
            _read(
                f"ins-alt-{index}",
                130,
                sequence="".join(insertion_alt),
                cigar="21M2I27M",
                reverse=index >= 2,
            )
        )
        records.append(
            _read(
                f"ins-ref-{index}",
                130,
                sequence="A" * 50,
                reverse=index >= 2,
            )
        )
    for index in range(4):
        records.append(
            _read(
                f"del-alt-{index}",
                180,
                sequence="A" * 50,
                cigar="21M2D29M",
                reverse=index >= 2,
            )
        )
        records.append(
            _read(
                f"del-ref-{index}",
                180,
                sequence="A" * 50,
                reverse=index >= 2,
            )
        )
    with pysam.AlignmentFile(str(alignment), "wb", header=header) as output:
        for record in sorted(records, key=lambda item: (item.reference_start, item.query_name)):
            output.write(record)
    pysam.index(str(alignment))

    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "target_id\tchrom\tpos\tref\talt\tcandidate_category\t"
        "local_complexity_status\tlocal_haplotype_status\n"
        "snv-1\tchr1\t101\tA\tT\tMOSAIC_CANDIDATE\tstandard\tnot_assessed\n"
        "ins-1\tchr1\t151\tA\tATG\tPRIORITY_VARIANT\tstandard\tnot_assessed\n"
        "del-1\tchr1\t201\tAAA\tA\tPRIORITY_VARIANT\tstandard\tnot_assessed\n"
        "complex-1\tchr1\t251\tAA\tATG\tPRIORITY_VARIANT\tstandard\tnot_assessed\n"
        "empty-1\tchr1\t301\tA\tT\tPRIORITY_VARIANT\tstandard\tnot_assessed\n",
        encoding="ascii",
    )
    provenance = ReadbackProvenance(
        reference_sha256=sha256_file(reference),
        technical_source_role="primary",
        caller_run_identity_sha256="2" * 64,
        resource_manifest_sha256="3" * 64,
        source_alignment_sha256=sha256_file(alignment),
        raw_artifact_id="bounded-readback-observations",
    )
    return alignment, reference, targets, provenance


def _replicate_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    replicated: list[dict[str, object]] = []
    for source in rows:
        row = dict(source)
        row["technical_source_role"] = "technical_replicate"
        row["caller_run_identity_sha256"] = "7" * 64
        row["source_alignment_sha256"] = "8" * 64
        row["raw_artifact_sha256"] = "9" * 64
        row["source_record_digest"] = hashlib.sha256(
            f"replicate:{row['target_id']}".encode()
        ).hexdigest()
        row["evidence_id"] = make_evidence_id("read_backed_evidence", row)
        replicated.append(row)
    return replicated


def _set_fragment_counts(row: dict[str, object], *, ref: int, alt: int) -> None:
    row["assessment_status"] = "assessed"
    row["eligible_fragments"] = ref + alt
    row["ref_reads"] = ref
    row["alt_reads"] = alt
    row["other_reads"] = 0
    row["forward_ref_reads"] = (ref + 1) // 2
    row["reverse_ref_reads"] = ref // 2
    row["forward_alt_reads"] = (alt + 1) // 2
    row["reverse_alt_reads"] = alt // 2
    row["f1r2_ref_reads"] = (ref + 1) // 2
    row["f2r1_ref_reads"] = ref // 2
    row["unclassified_orientation_ref_reads"] = 0
    row["f1r2_alt_reads"] = (alt + 1) // 2
    row["f2r1_alt_reads"] = alt // 2
    row["unclassified_orientation_alt_reads"] = 0
    row["mean_ref_mapping_quality"] = 60.0 if ref else None
    row["mean_alt_mapping_quality"] = 60.0 if alt else None
    row["mean_ref_base_quality"] = 40.0 if ref else None
    row["mean_alt_base_quality"] = 40.0 if alt else None
    row["mean_ref_read_position"] = 20.0 if ref else None
    row["mean_alt_read_position"] = 20.0 if alt else None
    row["filter_reasons"] = [] if alt >= 3 else ["NO_ALT_SUPPORT"]
    row["confidence_tier"] = "SUPPORTED" if alt >= 3 else "REVIEW"
    row["source_record_digest"] = hashlib.sha256(
        f"counts:{row['target_id']}:{row['source_alignment_sha256']}:{ref}:{alt}".encode()
    ).hexdigest()
    row["evidence_id"] = make_evidence_id("read_backed_evidence", row)


def test_snv_indel_overlap_and_no_evidence_accounting_is_deterministic(
    tmp_path: Path,
) -> None:
    alignment, reference, targets, provenance = build_inputs(tmp_path)
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    result_one = materialize_readback_targets(alignment, reference, targets, first, provenance)
    result_two = materialize_readback_targets(alignment, reference, targets, second, provenance)
    parsing = result_one["parsing"]
    assert parsing["targets"] == parsing["typed_rows"] == 5
    assert parsing["record_reconciliation_delta"] == 0
    assert parsing["assessment_status_counts"] == {
        "assessed": 3,
        "no_eligible_reads": 1,
        "unsupported_allele_representation": 1,
    }
    assert parsing["supported_rows"] == 2
    assert parsing["review_rows"] == 1
    assert parsing["unresolved_rows"] == 2
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(first)
    assert table.schema == schema_for("read_backed_evidence")
    rows = {row["target_id"]: row for row in table.to_pylist()}
    snv = rows["snv-1"]
    assert snv["eligible_fragments"] == 11
    assert (snv["ref_reads"], snv["alt_reads"], snv["other_reads"]) == (4, 5, 2)
    assert (snv["forward_alt_reads"], snv["reverse_alt_reads"]) == (3, 2)
    assert snv["excluded_low_mapping_quality_reads"] == 1
    assert snv["excluded_duplicate_reads"] == 1
    assert snv["collapsed_overlapping_mates"] == 2
    assert snv["discordant_overlapping_mates"] == 1
    assert snv["confidence_tier"] == "REVIEW"
    assert rows["ins-1"]["variant_type"] == "insertion"
    assert (rows["ins-1"]["ref_reads"], rows["ins-1"]["alt_reads"]) == (4, 4)
    assert rows["ins-1"]["nearby_non_target_indel_reads"] == 0
    assert rows["del-1"]["variant_type"] == "deletion"
    assert (rows["del-1"]["ref_reads"], rows["del-1"]["alt_reads"]) == (4, 4)
    assert rows["complex-1"]["assessment_status"] == "unsupported_allele_representation"
    assert rows["empty-1"]["assessment_status"] == "no_eligible_reads"
    rendered = json.dumps({"rows": rows, "summary": parsing}, sort_keys=True)
    assert SYNTHETIC_SAMPLE not in rendered
    assert "/raid/" not in rendered


def test_digest_ref_duplicate_target_and_private_binding_fail_closed(tmp_path: Path) -> None:
    alignment, reference, targets, provenance = build_inputs(tmp_path)
    wrong_digest = ReadbackProvenance(
        **{**provenance.__dict__, "source_alignment_sha256": "f" * 64}
    )
    with pytest.raises(PortableEvidenceError, match="alignment SHA-256 differs"):
        parse_readback_targets(alignment, reference, targets, wrong_digest)

    bad_ref = tmp_path / "bad-ref.tsv"
    bad_ref.write_text(
        targets.read_text(encoding="ascii").replace(
            "snv-1\tchr1\t101\tA\tT", "snv-1\tchr1\t101\tG\tT"
        ),
        encoding="ascii",
    )
    with pytest.raises(PortableEvidenceError, match="REF differs"):
        parse_readback_targets(alignment, reference, bad_ref, provenance)

    duplicate = tmp_path / "duplicate.tsv"
    lines = targets.read_text(encoding="ascii").splitlines()
    duplicate.write_text("\n".join([*lines, lines[1]]) + "\n", encoding="ascii")
    with pytest.raises(PortableEvidenceError, match="duplicate read-back target ID"):
        parse_readback_targets(alignment, reference, duplicate, provenance)

    leaking = ReadbackProvenance(
        **{
            **provenance.__dict__,
            "raw_artifact_id": f"readback-{SYNTHETIC_SAMPLE}",
        }
    )
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        parse_readback_targets(alignment, reference, targets, leaking)


def test_duplicate_inclusion_policy_is_explicit(tmp_path: Path) -> None:
    alignment, reference, targets, provenance = build_inputs(tmp_path)
    rows, _ = parse_readback_targets(
        alignment,
        reference,
        targets,
        provenance,
        policy=ReadbackPolicy(include_duplicates=True),
    )
    snv = next(row for row in rows if row["target_id"] == "snv-1")
    assert snv["duplicate_reads_included"] is True
    assert snv["excluded_duplicate_reads"] == 0
    assert snv["alt_reads"] == 6


def test_standalone_raw_capture_matches_direct_evidence_without_identifiers(
    tmp_path: Path,
) -> None:
    alignment, reference, targets, provenance = build_inputs(tmp_path)
    direct_rows, _ = parse_readback_targets(alignment, reference, targets, provenance)
    payload = build_capture(alignment, reference, targets)
    assert payload["status"] == "PASS"
    assert payload["record_count"] == len(direct_rows) == 5
    assert payload["sample_identifier_persisted"] is False
    assert payload["read_names_persisted"] is False
    assert payload["read_sequences_persisted"] is False

    by_target = {row["target_id"]: row for row in direct_rows}
    for captured in payload["records"]:
        direct = by_target[captured["target"]["target_id"]]
        assert captured["variant_type"] == direct["variant_type"]
        assert captured["assessment_status"] == direct["assessment_status"]
        assert captured["filter_reasons"] == direct["filter_reasons"]
        assert all(direct[key] == value for key, value in captured["metrics"].items())

    first = tmp_path / "first-capture.json"
    second = tmp_path / "second-capture.json"
    first_summary = write_capture(first, payload)
    second_summary = write_capture(second, build_capture(alignment, reference, targets))
    assert first.read_bytes() == second.read_bytes()
    assert first_summary["raw_artifact_sha256"] == second_summary["raw_artifact_sha256"]
    rendered = first.read_text(encoding="utf-8")
    assert SYNTHETIC_SAMPLE not in rendered
    assert str(tmp_path) not in rendered
    assert "snv-alt-0" not in rendered


def test_receipt_bound_capture_materializes_typed_rows_and_rejects_tamper(
    tmp_path: Path,
) -> None:
    alignment, reference, targets, provenance = build_inputs(tmp_path)
    capture = tmp_path / "capture.json"
    write_capture(capture, build_capture(alignment, reference, targets))

    rows, summary = parse_readback_capture(capture, provenance)
    assert summary["capture_records"] == summary["typed_rows"] == len(rows) == 5
    assert summary["record_reconciliation_delta"] == 0
    assert summary["raw_artifact_sha256"] == sha256_file(capture)
    assert all(row["raw_artifact_sha256"] == sha256_file(capture) for row in rows)
    assert all(row["raw_artifact_id"] == provenance.raw_artifact_id for row in rows)

    output = tmp_path / "capture.parquet"
    result = materialize_readback_capture(capture, output, provenance)
    assert result["status"] == "PASS"
    assert pq.read_table(output).schema == schema_for("read_backed_evidence")

    tampered_payload = json.loads(capture.read_text(encoding="utf-8"))
    tampered_payload["records"][0]["filter_reasons"] = ["NO_ALT_SUPPORT"]
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(tampered_payload), encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="assessment differs"):
        parse_readback_capture(tampered, provenance)


def test_readback_replicate_concordance_retains_every_source_row(tmp_path: Path) -> None:
    alignment, reference, targets, provenance = build_inputs(tmp_path)
    primary, _ = parse_readback_targets(alignment, reference, targets, provenance)
    replicate = _replicate_rows(primary)
    annotated, summary = annotate_readback_replicate_concordance(
        [*primary, *replicate],
        policy=ReadbackReplicatePolicy(minimum_ref_alt_fragments=5),
        expected_reference_sha256=provenance.reference_sha256,
    )
    assert len(annotated) == len(primary) + len(replicate) == 10
    assert summary["source_rows_removed"] == 0
    assert summary["allele_fraction_averaging_used"] is False
    assert summary["majority_vote_used"] is False
    assert summary["absence_used_as_reference"] is False
    assert summary["row_status_counts"] == {
        "comparable_concordant": 6,
        "not_comparable": 4,
    }
    assert summary["biological_fields_changed"] is False


def test_readback_replicate_discordance_is_review_not_a_vote(tmp_path: Path) -> None:
    alignment, reference, targets, provenance = build_inputs(tmp_path)
    primary_rows, _ = parse_readback_targets(alignment, reference, targets, provenance)
    primary = dict(next(row for row in primary_rows if row["target_id"] == "snv-1"))
    replicate = _replicate_rows([primary])[0]
    _set_fragment_counts(primary, ref=15, alt=15)
    _set_fragment_counts(replicate, ref=30, alt=0)
    annotated, summary = annotate_readback_replicate_concordance(
        [primary, replicate],
        expected_reference_sha256=provenance.reference_sha256,
    )
    assert summary["row_status_counts"] == {"comparable_discordant": 2}
    assert summary["group_relation_counts"] == {"alt_fraction_interval_discordant": 1}
    assert all(row["confidence_tier"] == "REVIEW" for row in annotated)
    assert all(
        "REPLICATE:READBACK_ALT_FRACTION_DISCORDANT" in row["filter_reasons"] for row in annotated
    )
    assert {row["alt_reads"] for row in annotated} == {0, 15}
