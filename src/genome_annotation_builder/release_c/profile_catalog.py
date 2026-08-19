"""Validate and describe the exact caller profiles shipped with Release C."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .caller_profile import load_caller_profile
from .portable_model import PortableEvidenceError

EXPECTED_PROFILE_IDENTITIES: Mapping[str, str] = {
    "cnvpytor-combined-rd-baf-1.3.2.json": (
        "63874b343264e8ae5c858f865e6c53aafd8d5ba0f5668492a8549c40d50e1e33"
    ),
    "delly-short-read-2.5.1.json": (
        "d081743a9c59fb3a848fefbcf9287bd5825f484209d0c711bcaa2c64490d84c8"
    ),
    "expansionhunter-grch38-female-5.0.0.json": (
        "4e2a03346cbbb6c2c694f5df1967c5219ce5549d6d5e17698a6d5ce043e30a26"
    ),
    "expansionhunter-grch38-male-5.0.0.json": (
        "84e93fb1e24cbebaeaee24424642fbbe63090d8c6fdcbbb689cc615849dc0063"
    ),
    "gatk-mitochondria-normal-shifted-4.6.2.0.json": (
        "69910e7de8dd0b7a7d806307455a9b717cbdd2e7faa5d740f47d5a6353a5ff72"
    ),
    "hla-la-short-read-1.0.4.json": (
        "31d97da9f13758bdab1b68eadff1324c80ef46f0ef5ff46cf28478ff111f7815"
    ),
    "manta-germline-1.6.0.json": (
        "aaa137729145e44caac0befe96a272c26aeb94a964c5d20450d7c09c02f8cd33"
    ),
    "readback-bounded-pysam-0.24.0.json": (
        "4a1a8d2e478d1158edebbcc496e3a7dd98b706d2f53967e433cd8c4ec1f79080"
    ),
    "spechla-short-read-bam-1.0.12.json": (
        "85ce0109fd6764337d833dd42b89aecc873c34ef795f15e78cedc4186350e2d1"
    ),
}


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_profile_catalog(
    profile_root: str | Path | None = None,
    *,
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    """Fail closed on profile omission, addition, content drift, or privacy leakage."""

    root = (
        Path(profile_root)
        if profile_root is not None
        else Path(__file__).resolve().with_name("profiles")
    )
    if root.is_symlink() or not root.is_dir():
        raise PortableEvidenceError("shipped caller-profile directory is unavailable")
    paths = sorted(root.glob("*.json"), key=lambda path: path.name)
    actual = {path.name for path in paths}
    expected = set(EXPECTED_PROFILE_IDENTITIES)
    if actual != expected:
        raise PortableEvidenceError("shipped caller-profile set differs from the release boundary")

    profiles: list[dict[str, Any]] = []
    profile_ids: set[str] = set()
    for path in paths:
        if path.is_symlink() or not path.is_file():
            raise PortableEvidenceError("shipped caller profile is not a regular file")
        profile = load_caller_profile(path, forbidden_bindings=forbidden_bindings)
        if profile.identity_sha256 != EXPECTED_PROFILE_IDENTITIES[path.name]:
            raise PortableEvidenceError("shipped caller-profile identity differs")
        summary = profile.portable_summary()
        profile_id = str(summary["profile_id"])
        if profile_id in profile_ids:
            raise PortableEvidenceError("shipped caller-profile ID is duplicated")
        profile_ids.add(profile_id)
        profiles.append(
            {
                "profile_file": path.name,
                **summary,
            }
        )
    payload = {
        "schema_version": 1,
        "status": "PASS",
        "profile_count": len(profiles),
        "profiles": profiles,
        "private_paths_persisted": False,
    }
    payload["catalog_identity_sha256"] = _canonical_sha256(payload)
    return payload


def resolve_shipped_profile(profile_id: str) -> Path:
    """Resolve one validated public profile ID for internal CLI use."""

    catalog = build_profile_catalog()
    matches = [row for row in catalog["profiles"] if row["profile_id"] == profile_id]
    if len(matches) != 1:
        raise PortableEvidenceError("shipped caller-profile ID is unknown or duplicated")
    return Path(__file__).resolve().with_name("profiles") / str(matches[0]["profile_file"])


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate and list the exact read-backed caller profiles in this release."
    )
    parser.parse_args(argv)
    try:
        print(json.dumps(build_profile_catalog(), indent=2, sort_keys=True))
    except (OSError, ValueError, PortableEvidenceError):
        print("error: caller_profile_catalog_invalid", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
