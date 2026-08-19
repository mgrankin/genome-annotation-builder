from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import caller_profile, profile_catalog
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError
from genome_annotation_builder.release_c.profile_catalog import (
    EXPECTED_PROFILE_IDENTITIES,
    build_profile_catalog,
)


def test_shipped_profile_catalog_is_exact_path_free_and_deterministic() -> None:
    first = build_profile_catalog(forbidden_bindings=("PRIVATE_SUBJECT",))
    second = build_profile_catalog(forbidden_bindings=("PRIVATE_SUBJECT",))
    assert first == second
    assert first["status"] == "PASS"
    assert first["profile_count"] == len(EXPECTED_PROFILE_IDENTITIES) == 9
    assert {row["profile_file"] for row in first["profiles"]} == set(EXPECTED_PROFILE_IDENTITIES)
    assert all(row["status"] == "PASS" for row in first["profiles"])
    rendered = json.dumps(first, sort_keys=True)
    assert "PRIVATE_SUBJECT" not in rendered
    assert "/raid/" not in rendered
    assert "/home/" not in rendered


def test_catalog_rejects_missing_extra_and_changed_profiles(tmp_path: Path) -> None:
    source = Path(profile_catalog.__file__).with_name("profiles")
    copied = tmp_path / "profiles"
    shutil.copytree(source, copied)

    missing = copied / next(iter(EXPECTED_PROFILE_IDENTITIES))
    missing.unlink()
    with pytest.raises(PortableEvidenceError, match="set differs"):
        build_profile_catalog(copied)

    shutil.copytree(source, copied, dirs_exist_ok=True)
    (copied / "unexpected.json").write_text("{}\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="set differs"):
        build_profile_catalog(copied)

    (copied / "unexpected.json").unlink()
    changed = copied / "manta-germline-1.6.0.json"
    payload = json.loads(changed.read_text(encoding="utf-8"))
    payload["licence"] += " changed"
    changed.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="identity differs"):
        build_profile_catalog(copied)


def test_profile_catalog_cli_emits_only_portable_json(capsys: pytest.CaptureFixture[str]) -> None:
    assert profile_catalog.main([]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "PASS"
    assert payload["profile_count"] == 9


def test_caller_profile_cli_resolves_a_validated_shipped_profile_id(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "caller-profile",
            "validate-profile",
            "--profile-id",
            "manta-germline-1.6.0",
        ],
    )
    assert caller_profile.main() == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "PASS"
    assert payload["profile_id"] == "manta-germline-1.6.0"

    monkeypatch.setattr(
        "sys.argv",
        ["caller-profile", "validate-profile", "--profile-id", "missing-profile"],
    )
    assert caller_profile.main() == 2
    assert "failed" in capsys.readouterr().err
