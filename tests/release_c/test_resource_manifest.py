from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import resource_manifest
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError
from genome_annotation_builder.release_c.resource_manifest import (
    build_resource_manifest,
    validate_resource_manifest_payload,
)

EXAMPLE = Path(__file__).resolve().parents[2] / "config/release-c-resource-metadata.example.json"


def metadata() -> dict[str, object]:
    return json.loads(EXAMPLE.read_text(encoding="utf-8"))


def test_resource_manifest_is_deterministic_exact_and_path_free() -> None:
    first = build_resource_manifest(metadata(), forbidden_bindings=("PRIVATE_SUBJECT",))
    second = build_resource_manifest(metadata(), forbidden_bindings=("PRIVATE_SUBJECT",))
    assert first == second
    summary = validate_resource_manifest_payload(first, forbidden_bindings=("PRIVATE_SUBJECT",))
    assert summary["status"] == "PASS_WITH_DOCUMENTED_LIMITATIONS"
    assert summary["profile_count"] == 9
    assert summary["private_paths_persisted"] is False
    rendered = json.dumps(first, sort_keys=True)
    assert "/raid/" not in rendered
    assert "/home/" not in rendered
    assert "PRIVATE_SUBJECT" not in rendered


def test_resource_manifest_rejects_catalog_identity_and_privacy_drift() -> None:
    valid = build_resource_manifest(metadata())
    changed = deepcopy(valid)
    changed["profile_catalog"]["profiles"][0]["profile_identity_sha256"] = "0" * 64
    with pytest.raises(PortableEvidenceError, match="profile catalog changed"):
        validate_resource_manifest_payload(changed)

    private = metadata()
    private["scope"] = "/" + "raid" + "/private/resource"
    with pytest.raises(PortableEvidenceError, match="private binding"):
        build_resource_manifest(private)


def test_resource_manifest_cli_creates_once_and_validates(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    output = tmp_path / "resource-manifest.json"
    assert (
        resource_manifest.main(["create", "--metadata", str(EXAMPLE), "--output", str(output)]) == 0
    )
    created = json.loads(capsys.readouterr().out)
    assert created["resume"] == "created"
    assert resource_manifest.main(["validate", "--manifest", str(output)]) == 0
    assert (
        json.loads(capsys.readouterr().out)["manifest_identity_sha256"]
        == created["manifest_identity_sha256"]
    )
    original = output.read_bytes()
    assert (
        resource_manifest.main(["create", "--metadata", str(EXAMPLE), "--output", str(output)]) == 2
    )
    assert output.read_bytes() == original
    assert "Traceback" not in capsys.readouterr().err
