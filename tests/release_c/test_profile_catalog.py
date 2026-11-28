from __future__ import annotations

import hashlib
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


@pytest.mark.parametrize(
    ("version", "profile_id", "profile_hash", "catalog_hash"),
    [
        (
            1,
            "readback-bounded-pysam-0.24.0",
            "4a1a8d2e478d1158edebbcc496e3a7dd98b706d2f53967e433cd8c4ec1f79080",
            "5500f314a26eb7148cb5f4fe945460a11e7ddcd20d56ff6274b9e2eacff8d5af",
        ),
        (
            2,
            "readback-bounded-pysam-0.24.0-v2",
            "b163b2afcfc402a019dc6c79a18a8b28e27725b9b61abc01b2f51e9b00d92747",
            "72bcbafb808168ec1bc3176847fca0ac9806cf2ced71d8abad54f1a85ac9a441",
        ),
    ],
)
def test_explicit_readback_profile_and_catalog_replay(
    version: int,
    profile_id: str,
    profile_hash: str,
    catalog_hash: str,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert profile_catalog.main(["--catalog-version", str(version)]) == 0
    catalog = json.loads(capsys.readouterr().out)
    assert catalog["profile_count"] == 9 and catalog["catalog_identity_sha256"] == catalog_hash
    assert [p["profile_id"] for p in catalog["profiles"] if p["gate_caller"] == "readback"] == [
        profile_id
    ]
    monkeypatch.setattr(
        "sys.argv", ["caller-profile", "validate-profile", "--profile-id", profile_id]
    )
    assert caller_profile.main() == 0
    actual = json.loads(capsys.readouterr().out)
    assert actual["profile_id"] == profile_id and actual["profile_identity_sha256"] == profile_hash
    path = profile_catalog.resolve_shipped_profile(profile_id)
    doc = json.loads(path.read_bytes())
    assert doc["caller_version"] == str(version)
    scripts = {m["role"]: m for m in doc["mounts"] if m["container_path"].startswith("/tools/")}
    entry = "readback_capture_v2" if version == 2 else "readback_capture"
    assert set(scripts) == ({"readback_capture", entry})
    assert doc["command_semantics"]["steps"][1]["arguments"][0] == "{" + entry + "}"
    for mount in scripts.values():
        filename = Path(mount["container_path"]).name
        shipped = path.parent.parent / filename
        assert hashlib.sha256(shipped.read_bytes()).hexdigest() == mount["expected_sha256"]
