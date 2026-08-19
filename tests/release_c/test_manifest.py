from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from genome_annotation_builder.release_c.manifest import (
    ManifestError,
    build_private_manifest,
    validate_manifest_structure,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    files = {
        name: tmp_path / name
        for name in (
            "reference.fa",
            "reference.fa.fai",
            "reference.dict",
            "alignment.bam",
            "alignment.bam.bai",
            "baf.vcf.gz",
            "baf.vcf.gz.tbi",
        )
    }
    for name, path in files.items():
        path.write_bytes(f"synthetic {name}\n".encode())
    config = {
        "schema_version": 1,
        "reference": {
            "assembly": "GRCh38 synthetic control",
            "fasta_path": files["reference.fa"].name,
            "fai_path": files["reference.fa.fai"].name,
            "dictionary_path": files["reference.dict"].name,
        },
        "caller_alignment_policy": {
            "manta": "original_markdup_bam",
            "cnvpytor": "original_markdup_bam",
        },
        "datasets": [
            {
                "dataset_id": "TECHNICAL_DATASET_01",
                "logical_subject": "SUBJECT_01",
                "technical_role": "primary",
                "laboratory": "LAB_01",
                "original_markdup_bam": {
                    "path": files["alignment.bam"].name,
                    "index_path": files["alignment.bam.bai"].name,
                    "sha256": digest(files["alignment.bam"]),
                },
                "baf_vcf": {
                    "path": files["baf.vcf.gz"].name,
                    "index_path": files["baf.vcf.gz.tbi"].name,
                    "sample": "VCF_SAMPLE_01",
                },
            }
        ],
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=True), encoding="utf-8")
    return path, files


def test_manifest_builder_hashes_every_bound_file_and_resolves_paths(tmp_path: Path) -> None:
    config, files = fixture(tmp_path)
    output = tmp_path / "private-manifest.json"
    summary = build_private_manifest(config, output)
    assert summary["status"] == "PASS"
    assert summary["alignment_hashes_deferred_to_precall_gate"] is False
    document = json.loads(output.read_text(encoding="utf-8"))
    assert validate_manifest_structure(document) == document
    alignment = document["datasets"][0]["original_markdup_bam"]
    assert alignment["path"] == str(files["alignment.bam"].resolve())
    assert alignment["sha256"] == digest(files["alignment.bam"])
    assert alignment["index_sha256"] == digest(files["alignment.bam.bai"])
    assert alignment["sha256_origin"] == "computed_by_manifest_builder"
    assert document["reference"]["fasta_sha256"] == digest(files["reference.fa"])
    assert document["datasets"][0]["baf_vcf"]["sha256"] == digest(files["baf.vcf.gz"])
    with pytest.raises(FileExistsError):
        build_private_manifest(config, output)


def test_precomputed_alignment_digest_is_explicitly_deferred_to_input_gate(
    tmp_path: Path,
) -> None:
    config, files = fixture(tmp_path)
    output = tmp_path / "private-manifest.json"
    summary = build_private_manifest(
        config,
        output,
        reuse_supplied_alignment_sha256=True,
    )
    assert summary["status"] == "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH"
    document = json.loads(output.read_text(encoding="utf-8"))
    alignment = document["datasets"][0]["original_markdup_bam"]
    assert alignment["sha256"] == digest(files["alignment.bam"])
    assert alignment["sha256_origin"].endswith("pending_precall_rehash")


def test_manifest_rejects_wrong_reference_policy_and_duplicate_dataset(tmp_path: Path) -> None:
    config, _files = fixture(tmp_path)
    value = yaml.safe_load(config.read_text(encoding="utf-8"))
    value["reference"]["assembly"] = "hg19"
    config.write_text(yaml.safe_dump(value), encoding="utf-8")
    with pytest.raises(ManifestError, match="GRCh38"):
        build_private_manifest(config, tmp_path / "bad-reference.json")

    config, _files = fixture(tmp_path)
    value = yaml.safe_load(config.read_text(encoding="utf-8"))
    value["datasets"].append(dict(value["datasets"][0]))
    config.write_text(yaml.safe_dump(value), encoding="utf-8")
    with pytest.raises(ManifestError, match="duplicated"):
        build_private_manifest(config, tmp_path / "duplicate.json")


def test_manifest_rejects_symlink_and_wrong_supplied_digest(tmp_path: Path) -> None:
    config, files = fixture(tmp_path)
    value = yaml.safe_load(config.read_text(encoding="utf-8"))
    value["datasets"][0]["original_markdup_bam"]["sha256"] = "0" * 64
    config.write_text(yaml.safe_dump(value), encoding="utf-8")
    with pytest.raises(ManifestError, match="digest does not match"):
        build_private_manifest(config, tmp_path / "wrong-digest.json")

    config, files = fixture(tmp_path)
    alias = tmp_path / "alignment-alias.bam"
    alias.symlink_to(files["alignment.bam"])
    value = yaml.safe_load(config.read_text(encoding="utf-8"))
    value["datasets"][0]["original_markdup_bam"]["path"] = alias.name
    config.write_text(yaml.safe_dump(value), encoding="utf-8")
    with pytest.raises(ManifestError, match="symbolic"):
        build_private_manifest(config, tmp_path / "symlink.json")


def test_manifest_cli_emits_only_summary(tmp_path: Path) -> None:
    config, _files = fixture(tmp_path)
    output = tmp_path / "PRIVATE_SUBJECT_MANIFEST.json"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "genome_annotation_builder.release_c.manifest",
            "create",
            "--config",
            str(config),
            "--output",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    summary = json.loads(completed.stdout)
    assert summary["status"] == "PASS"
    assert "SUBJECT_01" not in completed.stdout
    assert str(output) not in completed.stdout
