#!/usr/bin/env python3
"""Compare two release builds without masking payload drift."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

CHUNK_SIZE = 1024 * 1024


class ReproducibilityError(RuntimeError):
    """A release build differs at a boundary that must be reproducible."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            digest.update(chunk)
    return digest.hexdigest()


def _single(directory: Path, pattern: str, label: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ReproducibilityError(
            f"expected exactly one {label} in {directory}; found {len(matches)}"
        )
    path = matches[0]
    if path.is_symlink() or not path.is_file():
        raise ReproducibilityError(f"{label} is missing, non-regular, or symlinked: {path}")
    return path


def _member_digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _wheel_payload(path: Path) -> dict[str, str]:
    members: dict[str, str] = {}
    with zipfile.ZipFile(path) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            name = PurePosixPath(info.filename).as_posix()
            if name in members:
                raise ReproducibilityError(f"duplicate wheel member: {name}")
            members[name] = _member_digest(archive.read(info))
    if not members:
        raise ReproducibilityError(f"wheel has no regular payload members: {path}")
    return members


def _sdist_payload(path: Path) -> dict[str, str]:
    members: dict[str, str] = {}
    with tarfile.open(path, "r:gz") as archive:
        for info in archive.getmembers():
            if info.isdir():
                continue
            if not info.isfile():
                raise ReproducibilityError(f"non-regular source-distribution member: {info.name}")
            parts = PurePosixPath(info.name).parts
            if len(parts) < 2 or not parts[0].startswith("genome_annotation_builder-"):
                raise ReproducibilityError(f"unexpected source-distribution root: {info.name}")
            name = PurePosixPath(*parts[1:]).as_posix()
            if name in members:
                raise ReproducibilityError(f"duplicate source-distribution member: {name}")
            extracted = archive.extractfile(info)
            if extracted is None:
                raise ReproducibilityError(f"cannot read source-distribution member: {info.name}")
            members[name] = _member_digest(extracted.read())
    if not members:
        raise ReproducibilityError(f"source distribution has no payload members: {path}")
    return members


def _describe_difference(first: dict[str, str], second: dict[str, str]) -> str:
    first_names = set(first)
    second_names = set(second)
    removed = sorted(first_names - second_names)
    added = sorted(second_names - first_names)
    changed = sorted(name for name in first_names & second_names if first[name] != second[name])
    return f"removed={removed[:5]}, added={added[:5]}, changed={changed[:5]}"


def compare_distributions(first_directory: Path, second_directory: Path) -> dict[str, Any]:
    first_wheel = _single(first_directory, "*.whl", "wheel")
    second_wheel = _single(second_directory, "*.whl", "wheel")
    first_sdist = _single(first_directory, "*.tar.gz", "source distribution")
    second_sdist = _single(second_directory, "*.tar.gz", "source distribution")

    if first_wheel.name != second_wheel.name or first_sdist.name != second_sdist.name:
        raise ReproducibilityError("release archive names differ between builds")

    first_wheel_sha = _sha256(first_wheel)
    second_wheel_sha = _sha256(second_wheel)
    first_wheel_payload = _wheel_payload(first_wheel)
    second_wheel_payload = _wheel_payload(second_wheel)
    if first_wheel_payload != second_wheel_payload:
        raise ReproducibilityError(
            "wheel payload differs: "
            + _describe_difference(first_wheel_payload, second_wheel_payload)
        )
    if first_wheel_sha != second_wheel_sha:
        raise ReproducibilityError("wheel payload matches but archive bytes are not reproducible")

    first_sdist_payload = _sdist_payload(first_sdist)
    second_sdist_payload = _sdist_payload(second_sdist)
    if first_sdist_payload != second_sdist_payload:
        raise ReproducibilityError(
            "source-distribution payload differs: "
            + _describe_difference(first_sdist_payload, second_sdist_payload)
        )

    return {
        "schema_version": 1,
        "status": "PASS",
        "wheel": {
            "archive": first_wheel.name,
            "archive_byte_identical": True,
            "members": len(first_wheel_payload),
            "sha256": first_wheel_sha,
        },
        "source_distribution": {
            "archive": first_sdist.name,
            "archive_sha256_first": _sha256(first_sdist),
            "archive_sha256_second": _sha256(second_sdist),
            "payload_identical": True,
            "payload_members": len(first_sdist_payload),
            "policy": "compare logical member paths and bytes; ignore gzip/tar timestamp metadata",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("first_directory", type=Path)
    parser.add_argument("second_directory", type=Path)
    arguments = parser.parse_args()
    try:
        report = compare_distributions(arguments.first_directory, arguments.second_directory)
    except ReproducibilityError as exc:
        raise SystemExit(f"distribution reproducibility validation failed: {exc}") from exc
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
