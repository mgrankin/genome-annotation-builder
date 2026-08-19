from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from genome_annotation_builder.release_c.extension_stage import (
    validate_release_c_extension_stage,
)
from genome_annotation_builder.release_c.portable_model import (
    sha256_file,
)
from genome_annotation_builder.release_c.subject_stage import (
    SubjectStageError,
    _apply_concordance,
    build_subject_stage,
)

from .test_extension_stage import (
    ALIGNMENT_A,
    ALIGNMENT_B,
    PRIVATE_SUBJECT,
    caller_receipts,
    extension_inputs,
    gate_report,
    private_manifest,
    write_adapter_summary,
    write_json,
)
from .test_sv_concordance import sv_row


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def make_completed_matrix(tmp_path: Path, manifest: Path) -> Path:
    gate_manta = gate_report(tmp_path, "manta", ALIGNMENT_A)
    gate_mt = gate_report(tmp_path, "gatk_mitochondria", ALIGNMENT_B)
    receipts, provenance = caller_receipts(tmp_path, [gate_manta, gate_mt])
    manta_receipt = next(path for path in receipts if "manta" in path.name)
    typed_inputs = extension_inputs(tmp_path, provenance)

    matrix = tmp_path / "matrix"
    attempt = matrix / "job-0001/attempt-0001"
    lifecycle = attempt / "lifecycle"
    typed = attempt / "typed"
    lifecycle.mkdir(parents=True)
    typed.mkdir()
    (attempt / "bindings.private.json").write_text('{"schema_version":1}\n', encoding="utf-8")
    (lifecycle / "input-gate.json").write_bytes(gate_manta.read_bytes())
    (lifecycle / "caller-receipt.json").write_bytes(manta_receipt.read_bytes())
    structural = typed / "manta.parquet"
    structural.write_bytes(typed_inputs["structural_variants"].read_bytes())
    write_adapter_summary(
        typed / "adapter-summary.json",
        manta_receipt,
        profile_id="manta-germline-1.6.0",
        technical_source_role="primary",
        domain_files={"structural_variants": structural},
        manta_omitted_non_bnd_records=1,
    )
    plan = {
        "schema_version": 1,
        "status": "READY",
        "input_identities": {"private_manifest_sha256": sha256_file(manifest)},
        "job_count": 1,
        "jobs": [
            {
                "ordinal": 1,
                "dataset_id": "synthetic-source-a",
                "logical_subject": PRIVATE_SUBJECT,
                "technical_source_role": "primary",
                "expected_karyotype": "XX",
                "profile_id": "manta-germline-1.6.0",
            }
        ],
    }
    plan["plan_identity_sha256"] = canonical_sha256(plan)
    write_json(matrix / "matrix-plan.private.json", plan)
    write_json(
        matrix / "matrix-runtime.private.json",
        {
            "schema_version": 1,
            "state": "complete",
            "plan_identity_sha256": plan["plan_identity_sha256"],
            "completed_jobs": 1,
            "total_jobs": 1,
            "jobs": [
                {
                    "ordinal": 1,
                    "dataset_id": "synthetic-source-a",
                    "profile_id": "manta-germline-1.6.0",
                    "state": "complete",
                    "attempt": 1,
                    "last_completed_stage": "adapter",
                }
            ],
        },
    )
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "job_count": 1,
        "completed_jobs": 1,
        "dataset_identifiers_persisted": False,
        "logical_subject_identifiers_persisted": False,
        "private_paths_persisted": False,
    }
    summary["matrix_summary_identity_sha256"] = canonical_sha256(summary)
    write_json(matrix / "matrix-summary.json", summary)
    return matrix


def test_completed_matrix_builds_atomic_lossless_subject_stage_and_replays(
    tmp_path: Path,
) -> None:
    manifest = private_manifest(tmp_path)
    matrix = make_completed_matrix(tmp_path, manifest)
    work = tmp_path / "subject-work"
    output = tmp_path / "portable-stage"
    first = build_subject_stage(
        manifest,
        matrix,
        ["synthetic-source-a"],
        work,
        output,
        minimum_free_bytes=1,
    )
    assert first == {
        "schema_version": 1,
        "status": "PASS",
        "resume": "created",
        "domains": 1,
        "rows": 1,
        "source_dataset_count": 1,
        "source_rows_removed": 0,
        "private_paths_persisted": False,
        "subject_identifier_persisted": False,
    }
    validation = validate_release_c_extension_stage(output)
    assert validation["status"] == "PASS"
    assert validation["source_rows_removed"] == 0
    assembly_summary = json.loads((work / "subject-stage-summary.json").read_text(encoding="utf-8"))
    assert assembly_summary["input_rows"] == assembly_summary["output_rows"] == 1
    assert assembly_summary["concordance"]["structural_variants"]["status"] == (
        "NOT_APPLICABLE_SINGLE_SOURCE"
    )
    payload = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in output.rglob("*")
        if path.is_file()
    )
    assert PRIVATE_SUBJECT not in payload
    assert str(tmp_path) not in payload

    replay = build_subject_stage(
        manifest,
        matrix,
        ["synthetic-source-a"],
        work,
        output,
        minimum_free_bytes=1,
    )
    assert replay["resume"] == "reused"
    assert replay["rows"] == first["rows"]


def test_subject_concordance_dispatch_preserves_independent_caller_rows() -> None:
    rows = [
        sv_row("DELLY", "d1", start=100, end=300, source_digit="a"),
        sv_row("Manta", "m1", start=102, end=302, source_digit="b"),
    ]
    annotated, summary = _apply_concordance(
        "structural_variants",
        rows,
        reference_sha256="1" * 64,
        forbidden_bindings=(),
    )
    assert len(annotated) == len(rows) == 2
    assert summary["source_rows_removed"] == 0
    assert summary["row_status_counts"] == {"comparable_concordant": 2}


def test_subject_stage_refuses_an_incomplete_matrix_before_creating_work(
    tmp_path: Path,
) -> None:
    manifest = private_manifest(tmp_path)
    matrix = make_completed_matrix(tmp_path, manifest)
    runtime_path = matrix / "matrix-runtime.private.json"
    runtime = json.loads(runtime_path.read_text(encoding="utf-8"))
    runtime["state"] = "running"
    write_json(runtime_path, runtime)
    work = tmp_path / "must-not-exist"
    with pytest.raises(SubjectStageError, match="matrix_not_complete"):
        build_subject_stage(
            manifest,
            matrix,
            ["synthetic-source-a"],
            work,
            tmp_path / "must-not-publish",
            minimum_free_bytes=1,
        )
    assert not work.exists()
