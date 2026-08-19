from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for
from genome_annotation_builder.release_c.spechla_adapter import (
    SpecHLAProvenance,
    materialize_spechla_results,
    parse_spechla_results,
)

REFERENCE = "1" * 64
PROVENANCE = SpecHLAProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="spechla-result-set",
    caller_database_name="SpecHLA mapping and assembly scaffold",
    caller_database_version="1.0.12 immutable container database",
    caller_database_identity_sha256="5" * 64,
    nomenclature_database_name="IPD-IMGT/HLA annotation catalogue",
    nomenclature_version="3.65.0",
    nomenclature_identity_sha256="6" * 64,
)


PRIMARY_LOCI = ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")


def control_files(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Create three internally consistent, fully synthetic SpecHLA outputs."""

    sample = "SYNTHETIC_SUBJECT"
    header = ["Sample", *(f"HLA_{gene}_{copy}" for gene in PRIMARY_LOCI for copy in (1, 2))]
    full_calls: dict[str, str] = {}
    group_calls: dict[str, str] = {}
    for locus_index, gene in enumerate(PRIMARY_LOCI, start=1):
        for copy in (1, 2):
            first = locus_index + copy
            full_calls[f"HLA_{gene}_{copy}"] = f"{gene}*{first:02d}:01:01:01"
            group_calls[f"HLA_{gene}_{copy}"] = f"{gene}*{first:02d}:01:01G"
    full_calls["HLA_DQA1_1"] = "DQA1*05:05:17:01"
    group_calls["HLA_DQA1_1"] = "DQA1*05:01:01G"
    full_calls["HLA_B_1"] = "B*18:01:02"
    group_calls["HLA_B_1"] = "B*18:01:01G"

    version = "# version: IPD-IMGT/HLA 3.65.0"
    full = tmp_path / "synthetic-spechla-full.txt"
    group = tmp_path / "synthetic-spechla-g-group.txt"
    details = tmp_path / "synthetic-spechla-details.txt"
    full.write_text(
        "\n".join(
            (
                version,
                "\t".join(header),
                "\t".join([sample, *(full_calls[slot] for slot in header[1:])]),
            )
        )
        + "\n",
        encoding="utf-8",
    )
    group.write_text(
        "\n".join(
            (
                version,
                "\t".join(header),
                "\t".join([sample, *(group_calls[slot] for slot in header[1:])]),
            )
        )
        + "\n",
        encoding="utf-8",
    )

    detail_rows = [
        version,
        "Gene\tG_best\tallele\tdetails:allele;Score;Caucasian;Black;Asian",
    ]
    for slot in header[1:]:
        selected = full_calls[slot]
        candidates = [selected]
        scored = [f"{selected};100.000;0.10;0.20;0.30"]
        if slot == "HLA_B_1":
            candidates.append("B*18:40")
            scored = [
                "B*18:01;99.965;0.10;0.20;0.30",
                "B*18:40;99.965;0.01;0.02;0.03",
            ]
        detail_rows.append(
            "\t".join(
                (
                    slot,
                    selected,
                    ";".join(candidates) + ";",
                    *scored,
                )
            )
        )
    details.write_text("\n".join(detail_rows) + "\n", encoding="utf-8")
    return full, group, details


def test_synthetic_control_reconciles_all_three_artifacts_without_sample_leak(
    tmp_path: Path,
) -> None:
    full, group, details = control_files(tmp_path)
    raw_sample = full.read_text(encoding="utf-8").splitlines()[2].split("\t")[0]
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    result_one = materialize_spechla_results(
        full,
        group,
        details,
        first,
        PROVENANCE,
        forbidden_bindings=(raw_sample,),
    )
    result_two = materialize_spechla_results(
        full,
        group,
        details,
        second,
        PROVENANCE,
        forbidden_bindings=(raw_sample,),
    )
    parsing = result_one["parsing"]
    assert parsing["full_result_slots"] == 16
    assert parsing["g_group_result_slots"] == 16
    assert parsing["detail_result_slots"] == 16
    assert parsing["reconciled_slots"] == 16
    assert parsing["typed_loci"] == 8
    assert parsing["supported_loci"] == 8
    assert parsing["review_loci"] == parsing["unresolved_loci"] == 0
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(first)
    assert table.schema == schema_for("hla_typing")
    rows = table.to_pylist()
    assert len(rows) == 8
    assert all(row["call_status"] == "diploid_call" for row in rows)
    assert all(row["phase_resolved"] is False for row in rows)
    assert all(row["nomenclature_version"] == "3.65.0" for row in rows)
    assert all(row["caller_database_identity_sha256"] == "5" * 64 for row in rows)
    assert all(row["nomenclature_identity_sha256"] == "6" * 64 for row in rows)
    dqa1 = next(row for row in rows if row["locus"] == "HLA-DQA1")
    assert dqa1["allele1"] == "DQA1*05:05:17:01"
    assert dqa1["allele1_g_group"] == "DQA1*05:01:01G"
    assert dqa1["allele1_two_field"] == "DQA1*05:05"
    assert dqa1["filter_reasons"] == []
    b_locus = next(row for row in rows if row["locus"] == "HLA-B")
    assert b_locus["allele1_quality_score"] == pytest.approx(99.965)
    assert len(b_locus["allele1_scored_candidates"]) == 2
    assert raw_sample not in json.dumps(rows, sort_keys=True)
    assert raw_sample not in json.dumps(parsing, sort_keys=True)


def test_version_and_selected_detail_mismatches_fail_closed(tmp_path: Path) -> None:
    full, group, details = control_files(tmp_path)
    bad_version = tmp_path / "bad-version.txt"
    bad_version.write_text(
        group.read_text(encoding="utf-8").replace("3.65.0", "3.64.0", 1),
        encoding="utf-8",
    )
    with pytest.raises(PortableEvidenceError, match="version differs"):
        parse_spechla_results(full, bad_version, details, PROVENANCE)

    bad_details = tmp_path / "bad-details.txt"
    bad_details.write_text(
        details.read_text(encoding="utf-8").replace(
            "HLA_A_1\tA*02:01:01:01", "HLA_A_1\tA*31:01:02:01", 1
        ),
        encoding="utf-8",
    )
    with pytest.raises(PortableEvidenceError, match="selected allele"):
        parse_spechla_results(full, group, bad_details, PROVENANCE)


def test_missing_slot_and_dynamic_private_binding_fail_closed(tmp_path: Path) -> None:
    full, group, details = control_files(tmp_path)
    details_lines = details.read_text(encoding="utf-8").splitlines()
    missing = tmp_path / "missing-details.txt"
    missing.write_text("\n".join(details_lines[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="lacks allele slots"):
        parse_spechla_results(full, group, missing, PROVENANCE)

    raw_sample = full.read_text(encoding="utf-8").splitlines()[2].split("\t")[0]
    leaking = SpecHLAProvenance(
        **{
            **PROVENANCE.__dict__,
            "raw_artifact_id": f"spechla-{raw_sample}",
        }
    )
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        parse_spechla_results(full, group, details, leaking)
