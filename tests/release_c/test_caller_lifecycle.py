from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pysam
import pytest

from genome_annotation_builder.release_c import oci_runner
from genome_annotation_builder.release_c.caller_lifecycle import (
    CallerLifecycleError,
    run_caller_lifecycle,
)
from genome_annotation_builder.release_c.caller_profile import load_caller_profile

from .test_caller_profile import make_fixture, write_json
from .test_manta_adapter import write_vcf

PRIVATE_DATASET = "PRIVATE_DATASET_001"
PRIVATE_SUBJECT = "PRIVATE_SUBJECT_001"
PRIVATE_LAB = "PRIVATE_LAB_001"


def private_manifest(path: Path, fixture: dict[str, Path]) -> Path:
    manifest = path / "private-manifest.json"
    write_json(
        manifest,
        {
            "schema_version": 1,
            "status": "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH",
            "caller_alignment_policy": {"manta": "original_markdup_bam"},
            "reference": {
                "fasta_path": str(fixture["reference"]),
                "fai_path": str(fixture["reference_fai"]),
            },
            "datasets": [
                {
                    "dataset_id": PRIVATE_DATASET,
                    "logical_subject": PRIVATE_SUBJECT,
                    "laboratory": PRIVATE_LAB,
                    "original_markdup_bam": {
                        "path": str(fixture["alignment"]),
                        "index_path": str(fixture["alignment_index"]),
                    },
                }
            ],
        },
    )
    return manifest


def install_manta_runtime(
    monkeypatch: pytest.MonkeyPatch, fixture: dict[str, Path], commands: list[list[str]]
) -> None:
    profile = load_caller_profile(fixture["profile"])

    def capture(command: list[str]) -> subprocess.CompletedProcess[str]:
        if command[1:] == ["--version"]:
            return subprocess.CompletedProcess(command, 0, "Docker synthetic\n", "")
        return subprocess.CompletedProcess(
            command, 0, json.dumps([profile.payload["image"]]) + "\n", ""
        )

    def execute(command: list[str], stdout_handle, stderr_handle) -> int:
        commands.append(command)
        stdout_handle.write(b"synthetic caller output\n")
        stderr_handle.write(b"")
        if len(commands) == 3:
            artifact = fixture["workspace"] / "manta/results/variants/diploidSV.vcf.gz"
            artifact.parent.mkdir(parents=True)
            source = artifact.with_suffix("")
            write_vcf(source)
            pysam.tabix_compress(str(source), str(artifact), force=True)
        return 0

    monkeypatch.setattr(oci_runner, "_capture", capture)
    monkeypatch.setattr(oci_runner, "_execute", execute)


def test_lifecycle_completes_then_reuses_every_immutable_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_fixture(tmp_path / "inputs")
    manifest = private_manifest(tmp_path, fixture)
    lifecycle = tmp_path / "lifecycle"
    lifecycle.mkdir()
    (lifecycle / "input-gate.json").write_bytes(fixture["gate"].read_bytes())
    commands: list[list[str]] = []
    install_manta_runtime(monkeypatch, fixture, commands)

    first = run_caller_lifecycle(
        manifest,
        PRIVATE_DATASET,
        fixture["bindings"],
        fixture["resources"],
        lifecycle,
        profile_path=fixture["profile"],
        minimum_free_bytes=1,
    )
    assert first["status"] == "PASS"
    assert first["resume"] == {
        "created_checkpoints": [
            "prepared_run",
            "execution_report",
            "caller_receipt",
            "lifecycle_summary",
        ],
        "reused_checkpoints": ["input_gate"],
    }
    assert len(commands) == 3
    portable = (lifecycle / "lifecycle-summary.json").read_text(encoding="utf-8")
    for private in (
        PRIVATE_DATASET,
        PRIVATE_SUBJECT,
        PRIVATE_LAB,
        str(manifest),
        str(fixture["alignment"]),
        str(fixture["workspace"]),
    ):
        assert private not in portable

    second = run_caller_lifecycle(
        manifest,
        PRIVATE_DATASET,
        fixture["bindings"],
        fixture["resources"],
        lifecycle,
        profile_path=fixture["profile"],
        minimum_free_bytes=1,
    )
    assert second["status"] == "PASS"
    assert second["resume"] == {
        "created_checkpoints": [],
        "reused_checkpoints": [
            "input_gate",
            "prepared_run",
            "execution_report",
            "caller_receipt",
            "lifecycle_summary",
        ],
    }
    assert len(commands) == 3


def test_lifecycle_rejects_changed_final_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_fixture(tmp_path / "inputs")
    manifest = private_manifest(tmp_path, fixture)
    lifecycle = tmp_path / "lifecycle"
    lifecycle.mkdir()
    (lifecycle / "input-gate.json").write_bytes(fixture["gate"].read_bytes())
    commands: list[list[str]] = []
    install_manta_runtime(monkeypatch, fixture, commands)
    run_caller_lifecycle(
        manifest,
        PRIVATE_DATASET,
        fixture["bindings"],
        fixture["resources"],
        lifecycle,
        profile_path=fixture["profile"],
        minimum_free_bytes=1,
    )
    summary_path = lifecycle / "lifecycle-summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["status"] = "CHANGED"
    write_json(summary_path, summary)

    with pytest.raises(CallerLifecycleError, match="existing_lifecycle_summary_changed"):
        run_caller_lifecycle(
            manifest,
            PRIVATE_DATASET,
            fixture["bindings"],
            fixture["resources"],
            lifecycle,
            profile_path=fixture["profile"],
            minimum_free_bytes=1,
        )
