from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome_annotation_builder.cli import main as root_cli
from genome_annotation_builder.release_c import __main__ as release_c_cli


def test_dispatcher_help_lists_generic_commands(capsys: pytest.CaptureFixture[str]) -> None:
    assert release_c_cli.main([]) == 0
    output = capsys.readouterr().out
    assert "manifest" in output
    assert "input-gate" in output
    assert "binding-builder" in output
    assert "caller-lifecycle" in output
    assert "profile-catalog" in output
    assert "profile-adapter" in output
    assert "resource-manifest" in output
    assert "readback-capture" in output
    assert "readback-targets" in output
    assert "readback-concordance" in output
    assert "extension-stage" in output
    assert "bundle-extension" in output
    assert "matrix-runner" in output
    assert "subject-stage" in output
    assert "PRIVATE_SUBJECT_TOKEN" not in output


def test_dispatcher_runs_manifest_structure_validation(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    manifest = {
        "schema_version": 1,
        "status": "PASS",
        "reference": {
            "assembly": "GRCh38 synthetic control",
            "fasta_path": "/synthetic/reference.fa",
            "fasta_bytes": 1,
            "fasta_sha256": "1" * 64,
            "fai_path": "/synthetic/reference.fa.fai",
            "fai_bytes": 1,
            "fai_sha256": "2" * 64,
            "dictionary_path": "/synthetic/reference.dict",
            "dictionary_bytes": 1,
            "dictionary_sha256": "3" * 64,
        },
        "caller_alignment_policy": {"manta": "original_markdup_bam"},
        "datasets": [
            {
                "dataset_id": "DATASET_01",
                "logical_subject": "SUBJECT_01",
                "technical_role": "primary",
                "original_markdup_bam": {
                    "path": "/synthetic/alignment.bam",
                    "bytes": 1,
                    "sha256": "4" * 64,
                    "sha256_origin": "computed_by_manifest_builder",
                    "index_path": "/synthetic/alignment.bam.bai",
                    "index_bytes": 1,
                    "index_sha256": "5" * 64,
                },
            }
        ],
    }
    path = tmp_path / "private-manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    assert release_c_cli.main(["manifest", "validate-structure", "--manifest", str(path)]) == 0
    assert json.loads(capsys.readouterr().out) == {
        "datasets": 1,
        "schema_version": 1,
        "status": "PASS",
    }


def test_root_cli_dispatches_release_c_help(capsys: pytest.CaptureFixture[str]) -> None:
    assert root_cli(["release-c"]) == 0
    output = capsys.readouterr().out
    assert "checksum-bound" in output
    assert "caller-profile" in output
