from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import readback_capture
from genome_annotation_builder.release_c.binding_builder import (
    BindingBuilderError,
    build_caller_bindings,
    validate_caller_bindings,
)
from genome_annotation_builder.release_c.caller_profile import prepare_caller_profile_run

from .test_caller_profile import (
    READBACK_PROFILE_PATH,
    make_fixture,
    sha256_file,
    write_json,
)

PRIVATE_DATASET = "PRIVATE_DATASET_002"
PRIVATE_SUBJECT = "PRIVATE_SUBJECT_002"


def manifest_file(path: Path) -> dict[str, object]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def manifest_index(path: Path) -> dict[str, object]:
    return {
        "index_path": str(path),
        "index_bytes": path.stat().st_size,
        "index_sha256": sha256_file(path),
    }


def make_manifest(tmp_path: Path, fixture: dict[str, Path], policies: dict[str, str]) -> Path:
    manifest = tmp_path / "private-input-manifest.json"
    alignment = {
        **manifest_file(fixture["alignment"]),
        **manifest_index(fixture["alignment_index"]),
    }
    write_json(
        manifest,
        {
            "schema_version": 1,
            "status": "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH",
            "caller_alignment_policy": policies,
            "reference": {
                "fasta_path": str(fixture["reference"]),
                "fasta_sha256": sha256_file(fixture["reference"]),
                "fai_path": str(fixture["reference_fai"]),
                "fai_sha256": sha256_file(fixture["reference_fai"]),
            },
            "datasets": [
                {
                    "dataset_id": PRIVATE_DATASET,
                    "logical_subject": PRIVATE_SUBJECT,
                    "technical_role": "primary",
                    "original_markdup_bam": alignment,
                }
            ],
        },
    )
    return manifest


def test_builder_compiles_manta_bindings_that_prepare_exact_profile(
    tmp_path: Path,
) -> None:
    fixture = make_fixture(tmp_path / "inputs")
    manifest = make_manifest(tmp_path, fixture, {"manta": "original_markdup_bam"})
    config = tmp_path / "resource-bindings.json"
    write_json(
        config,
        {
            "schema_version": 1,
            "profiles": {"manta-germline-1.6.0": {"roles": {}}},
            "dataset_profiles": {},
        },
    )
    workspace = tmp_path / "workspace"
    output = tmp_path / "bindings.json"
    summary = build_caller_bindings(
        manifest,
        PRIVATE_DATASET,
        fixture["profile"],
        config,
        fixture["resources"],
        workspace,
        output,
        user="1000:1000",
        minimum_free_bytes=1,
    )
    assert summary["status"] == "PASS"
    assert summary["mount_count"] == 6
    assert summary["private_paths_persisted"] is False
    assert PRIVATE_DATASET not in json.dumps(summary)
    assert PRIVATE_SUBJECT not in json.dumps(summary)
    validated = validate_caller_bindings(
        manifest,
        PRIVATE_DATASET,
        fixture["profile"],
        config,
        fixture["resources"],
        workspace,
        output,
        user="1000:1000",
        minimum_free_bytes=1,
    )
    assert validated["existing_binding_revalidated"] is True
    prepared = prepare_caller_profile_run(
        fixture["profile"],
        output,
        fixture["gate"],
        fixture["resources"],
        tmp_path / "plan",
    )
    assert prepared["status"] == "PREPARED"

    with pytest.raises(BindingBuilderError, match="refusing to replace"):
        build_caller_bindings(
            manifest,
            PRIVATE_DATASET,
            fixture["profile"],
            config,
            fixture["resources"],
            workspace,
            output,
            user="1000:1000",
            minimum_free_bytes=1,
        )

    changed = json.loads(output.read_text(encoding="utf-8"))
    changed["minimum_free_bytes"] = 2
    write_json(output, changed)
    with pytest.raises(BindingBuilderError, match="differs from recomputed"):
        validate_caller_bindings(
            manifest,
            PRIVATE_DATASET,
            fixture["profile"],
            config,
            fixture["resources"],
            workspace,
            output,
            user="1000:1000",
            minimum_free_bytes=1,
        )


def test_builder_merges_shared_and_dataset_specific_readback_resources(
    tmp_path: Path,
) -> None:
    fixture = make_fixture(tmp_path / "inputs")
    manifest = make_manifest(tmp_path, fixture, {"readback": "original_markdup_bam"})
    targets = tmp_path / "private-targets.tsv"
    targets.write_text(
        "target_id\tchrom\tpos\tref\talt\tcandidate_category\t"
        "local_complexity_status\tlocal_haplotype_status\n"
        "candidate-1\tchr1\t1\tA\tC\treview\tclear\tconsistent\n",
        encoding="utf-8",
    )
    capture = Path(readback_capture.__file__)
    profile_payload = json.loads(READBACK_PROFILE_PATH.read_text(encoding="utf-8"))
    reference_digests = {
        "reference_fasta": sha256_file(fixture["reference"]),
        "reference_fai": sha256_file(fixture["reference_fai"]),
    }
    for mount in profile_payload["mounts"]:
        if mount["role"] in reference_digests:
            mount["expected_sha256"] = reference_digests[mount["role"]]
    readback_profile = tmp_path / "readback-test-profile.json"
    write_json(readback_profile, profile_payload)
    config = tmp_path / "resource-bindings.json"
    write_json(
        config,
        {
            "schema_version": 1,
            "profiles": {
                "readback-bounded-pysam-0.24.0": {
                    "roles": {"readback_capture": {"path": str(capture)}}
                }
            },
            "dataset_profiles": {
                PRIVATE_DATASET: {
                    "readback-bounded-pysam-0.24.0": {
                        "roles": {"target_manifest": {"path": str(targets)}}
                    }
                }
            },
        },
    )
    output = tmp_path / "readback-bindings.json"
    summary = build_caller_bindings(
        manifest,
        PRIVATE_DATASET,
        readback_profile,
        config,
        fixture["resources"],
        tmp_path / "readback-workspace",
        output,
        user="1000:1000",
        minimum_free_bytes=1,
    )
    assert summary["status"] == "PASS"
    assert set(json.loads(output.read_text(encoding="utf-8"))["mounts"]) == {
        "readback_capture",
        "reference_fai",
        "reference_fasta",
        "resource_manifest",
        "selected_alignment",
        "selected_alignment_index",
        "target_manifest",
        "workspace",
    }


def test_builder_rejects_extra_or_missing_resource_roles(tmp_path: Path) -> None:
    fixture = make_fixture(tmp_path / "inputs")
    manifest = make_manifest(tmp_path, fixture, {"manta": "original_markdup_bam"})
    config = tmp_path / "resource-bindings.json"
    write_json(
        config,
        {
            "schema_version": 1,
            "profiles": {
                "manta-germline-1.6.0": {
                    "roles": {"unexpected_role": {"path": str(fixture["reference"])}}
                }
            },
            "dataset_profiles": {},
        },
    )
    with pytest.raises(BindingBuilderError, match="unresolved mount roles"):
        build_caller_bindings(
            manifest,
            PRIVATE_DATASET,
            fixture["profile"],
            config,
            fixture["resources"],
            tmp_path / "workspace",
            tmp_path / "bindings.json",
            user="1000:1000",
            minimum_free_bytes=1,
        )
