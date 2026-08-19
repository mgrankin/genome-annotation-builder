from __future__ import annotations

import json
from pathlib import Path

import yaml

from genome_annotation_builder.release_c.matrix_runner import (
    DEFAULT_FIXED_PROFILES,
    DEFAULT_KARYOTYPE_PROFILES,
)
from genome_annotation_builder.release_c.profile_catalog import build_profile_catalog


def test_release_c_example_config_covers_every_shipped_gate_caller() -> None:
    root = Path(__file__).resolve().parents[2]
    example = yaml.safe_load((root / "config/release-c.example.yaml").read_text(encoding="utf-8"))
    catalog = build_profile_catalog()
    expected_callers = {row["gate_caller"] for row in catalog["profiles"]}
    assert example["schema_version"] == 1
    assert set(example["caller_alignment_policy"]) == expected_callers
    assert len(example["datasets"]) == 1
    dataset = example["datasets"][0]
    assert dataset["dataset_id"] == "TECHNICAL_DATASET_01"
    assert dataset["logical_subject"] == "SUBJECT_01"
    rendered = json.dumps(example, sort_keys=True)
    assert "/raid/" not in rendered
    assert "/home/" not in rendered


def test_release_c_matrix_example_matches_the_default_public_matrix() -> None:
    root = Path(__file__).resolve().parents[2]
    example = json.loads(
        (root / "config/release-c-matrix.example.json").read_text(encoding="utf-8")
    )
    assert example == {
        "schema_version": 1,
        "fixed_profiles": list(DEFAULT_FIXED_PROFILES),
        "karyotype_profiles": DEFAULT_KARYOTYPE_PROFILES,
    }
