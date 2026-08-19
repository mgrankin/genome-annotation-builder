from __future__ import annotations

import json
from pathlib import Path

from genome_annotation_builder.release_c import binding_builder, caller_profile
from genome_annotation_builder.release_c.profile_adapter import HLA_OPTION_FIELDS
from genome_annotation_builder.release_c.profile_catalog import EXPECTED_PROFILE_IDENTITIES

EXAMPLE = Path(__file__).resolve().parents[2] / "config" / "release-c-resources.example.json"
ADAPTER_EXAMPLE = Path(__file__).resolve().parents[2] / "config" / "release-c-adapters.example.json"
PROFILE_ROOT = Path(caller_profile.__file__).with_name("profiles")
BUILTIN_ROLES = binding_builder.BUILTIN_FILE_ROLES | {"workspace"}


def test_resource_example_resolves_every_shipped_profile_mount_role() -> None:
    config = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    assert config["schema_version"] == 1
    assert set(config["profiles"]) == {
        path.removesuffix(".json") for path in EXPECTED_PROFILE_IDENTITIES
    }
    dataset = config["dataset_profiles"]["TECHNICAL_DATASET_01"]
    for profile_file in EXPECTED_PROFILE_IDENTITIES:
        payload = json.loads((PROFILE_ROOT / profile_file).read_text(encoding="utf-8"))
        profile_id = payload["profile_id"]
        expected = {mount["role"] for mount in payload["mounts"]} - BUILTIN_ROLES
        shared = config["profiles"][profile_id]["roles"]
        overlay = dataset.get(profile_id, {}).get("roles", {})
        assert set(shared) | set(overlay) == expected


def test_resource_example_has_only_generic_placeholders() -> None:
    text = EXAMPLE.read_text(encoding="utf-8")
    assert "/raid/" not in text
    assert "/home/" not in text
    assert "TECHNICAL_DATASET_01" in text
    assert "SUBJECT_" not in text


def test_adapter_example_is_generic_and_covers_both_hla_profiles() -> None:
    config = json.loads(ADAPTER_EXAMPLE.read_text(encoding="utf-8"))
    assert config["schema_version"] == 1
    assert set(config["profiles"]) == {
        "hla-la-short-read-1.0.4",
        "spechla-short-read-bam-1.0.12",
    }
    for options in config["profiles"].values():
        assert set(options) == HLA_OPTION_FIELDS
    text = ADAPTER_EXAMPLE.read_text(encoding="utf-8")
    assert "/raid/" not in text
    assert "/home/" not in text
    assert "SUBJECT_" not in text
