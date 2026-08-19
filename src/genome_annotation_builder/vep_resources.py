"""Install and verify the immutable Ensembl VEP 116 GRCh38 cache."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tarfile
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO

from .config import BuildConfig
from .constants import (
    VEP_CACHE_ARCHIVE,
    VEP_CACHE_ARCHIVE_1K_BLOCKS,
    VEP_CACHE_ARCHIVE_BSD_SUM,
    VEP_CACHE_ARCHIVE_BYTES,
    VEP_CACHE_ARCHIVE_SHA256,
    VEP_CACHE_ASSEMBLY,
    VEP_CACHE_CHECKSUMS_SHA256,
    VEP_CACHE_CHECKSUMS_URL,
    VEP_CACHE_URL,
    VEP_CACHE_VERSION,
    VEP_VERSION,
)
from .errors import ResourceError
from .util import atomic_output_directory, atomic_write_json, directory_digest, sha256_file

EXPECTED_CACHE_METADATA = {
    "source_assembly": "GRCh38.p14",
    "source_genebuild": "GENCODE50",
    "source_gencode": "GENCODE 50",
    "source_sift": "6.2.1",
    "source_polyphen": "2.2.3",
    "source_dbSNP": "156",
    "source_ClinVar": "202509",
    "source_gnomADe": "v4.1",
    "source_gnomADg": "v4.1",
}


def _bsd_sum(handle: BinaryIO) -> int:
    """Return the historical 16-bit BSD checksum used by Ensembl CHECKSUMS."""
    checksum = 0
    while chunk := handle.read(8 * 1024 * 1024):
        for byte in chunk:
            checksum = (checksum >> 1) | ((checksum & 1) << 15)
            checksum = (checksum + byte) & 0xFFFF
    return checksum


def _archive_identity(path: Path) -> dict[str, Any]:
    size = path.stat().st_size
    sum_executable = shutil.which("sum")
    if sum_executable is None:
        raise ResourceError("POSIX/GNU sum is required to verify Ensembl's cache checksum")
    completed = subprocess.run(
        [sum_executable, "-r", str(path)], check=False, capture_output=True, text=True
    )
    if completed.returncode != 0:
        raise ResourceError(f"could not calculate Ensembl BSD checksum: {completed.stderr}")
    try:
        bsd_sum, blocks = (int(value) for value in completed.stdout.split()[:2])
    except (ValueError, IndexError) as exc:
        raise ResourceError(f"could not parse BSD checksum output: {completed.stdout!r}") from exc
    return {
        "bytes": size,
        "sha256": sha256_file(path),
        "bsd_sum": bsd_sum,
        "blocks_1k": blocks,
    }


def _verify_archive(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ResourceError(f"VEP cache archive not found: {path}")
    identity = _archive_identity(path)
    errors: list[str] = []
    if identity["bytes"] != VEP_CACHE_ARCHIVE_BYTES:
        errors.append(f"bytes {identity['bytes']} != {VEP_CACHE_ARCHIVE_BYTES}")
    if identity["bsd_sum"] != VEP_CACHE_ARCHIVE_BSD_SUM:
        errors.append(f"BSD sum {identity['bsd_sum']} != {VEP_CACHE_ARCHIVE_BSD_SUM}")
    if identity["blocks_1k"] != VEP_CACHE_ARCHIVE_1K_BLOCKS:
        errors.append(f"1K blocks {identity['blocks_1k']} != {VEP_CACHE_ARCHIVE_1K_BLOCKS}")
    if VEP_CACHE_ARCHIVE_SHA256 and identity["sha256"] != VEP_CACHE_ARCHIVE_SHA256:
        errors.append(f"SHA-256 {identity['sha256']} != {VEP_CACHE_ARCHIVE_SHA256}")
    if errors:
        raise ResourceError("invalid Ensembl VEP cache archive: " + "; ".join(errors))
    return identity


def _resource_root(config: BuildConfig) -> Path:
    return config.vep_cache_manifest.parent


def _archive_path(config: BuildConfig) -> Path:
    return _resource_root(config) / "archive" / VEP_CACHE_ARCHIVE


def _download_archive(config: BuildConfig, destination: Path) -> None:
    if not config.allow_download:
        raise ResourceError(
            "VEP cache archive is absent and resources.allow_download is false; "
            f"place {VEP_CACHE_ARCHIVE} at {destination}"
        )
    curl = shutil.which("curl")
    if curl is None:
        raise ResourceError("curl is required for restart-safe VEP cache download")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f".{destination.name}.partial")
    completed = subprocess.run(
        [
            curl,
            "--fail",
            "--location",
            "--retry",
            "8",
            "--retry-delay",
            "5",
            "--continue-at",
            "-",
            "--output",
            str(partial),
            VEP_CACHE_URL,
        ],
        check=False,
    )
    if completed.returncode != 0:
        raise ResourceError(
            f"VEP cache download failed with exit {completed.returncode}; partial retained at {partial}"
        )
    _verify_archive(partial)
    os.replace(partial, destination)


def _validate_tar_members(archive: Path) -> dict[str, int]:
    files = 0
    directories = 0
    try:
        with tarfile.open(archive, mode="r:gz") as handle:
            for member in handle:
                path = PurePosixPath(member.name)
                if path.is_absolute() or ".." in path.parts:
                    raise ResourceError(f"unsafe path in VEP cache archive: {member.name}")
                if not path.parts or path.parts[0] != "homo_sapiens":
                    raise ResourceError(f"unexpected VEP archive root: {member.name}")
                if member.issym() or member.islnk() or member.isdev() or member.isfifo():
                    raise ResourceError(f"unsafe member type in VEP cache archive: {member.name}")
                if member.isdir():
                    directories += 1
                elif member.isfile():
                    files += 1
                else:
                    raise ResourceError(f"unsupported member type in VEP archive: {member.name}")
    except (tarfile.TarError, OSError) as exc:
        raise ResourceError(f"cannot validate VEP cache archive: {exc}") from exc
    if files == 0:
        raise ResourceError("VEP cache archive contains no files")
    return {"archive_files": files, "archive_directories": directories}


def _extract_archive(archive: Path, destination: Path) -> dict[str, int]:
    members = _validate_tar_members(archive)
    tar = shutil.which("tar")
    if tar is None:
        raise ResourceError("GNU tar is required to extract the VEP cache archive")
    with atomic_output_directory(destination) as staging:
        completed = subprocess.run(
            [
                tar,
                "--extract",
                "--gzip",
                "--file",
                str(archive),
                "--directory",
                str(staging),
                "--no-same-owner",
                "--no-same-permissions",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise ResourceError(f"VEP cache extraction failed: {completed.stderr[-4000:]}")
        expected = staging / "homo_sapiens" / f"{VEP_CACHE_VERSION}_{VEP_CACHE_ASSEMBLY}"
        if not expected.is_dir():
            raise ResourceError(f"extracted VEP cache lacks {expected.relative_to(staging)}")
    return members


def _tree_statistics(path: Path) -> dict[str, int | str]:
    files = [item for item in path.rglob("*") if item.is_file()]
    return {
        "tree_sha256": directory_digest(path),
        "tree_files": len(files),
        "tree_bytes": sum(item.stat().st_size for item in files),
    }


def _cache_metadata(cache_dir: Path) -> dict[str, str]:
    info = cache_dir / "homo_sapiens" / f"{VEP_CACHE_VERSION}_{VEP_CACHE_ASSEMBLY}" / "info.txt"
    try:
        with info.open(encoding="utf-8") as handle:
            rows = [line.rstrip("\n").split("\t", maxsplit=1) for line in handle]
    except OSError as exc:
        raise ResourceError(f"cannot read VEP cache metadata: {exc}") from exc
    metadata = {row[0]: row[1] for row in rows if len(row) == 2 and row[0] and row[1]}
    for field, expected in EXPECTED_CACHE_METADATA.items():
        if metadata.get(field) != expected:
            raise ResourceError(
                f"VEP cache metadata {field} mismatch: {metadata.get(field)!r} != {expected!r}"
            )
    return {field: metadata[field] for field in sorted(EXPECTED_CACHE_METADATA)}


def verify_vep_cache(config: BuildConfig, *, full_digest: bool = True) -> dict[str, Any]:
    try:
        manifest = json.loads(config.vep_cache_manifest.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ResourceError(f"VEP cache manifest not found: {config.vep_cache_manifest}") from exc
    except json.JSONDecodeError as exc:
        raise ResourceError(f"invalid VEP cache manifest: {exc}") from exc
    required = {
        "resource_type",
        "vep_version",
        "cache_version",
        "assembly",
        "archive_sha256",
        "tree_sha256",
        "tree_files",
        "tree_bytes",
        "cache_metadata",
    }
    if not isinstance(manifest, dict) or not required.issubset(manifest):
        raise ResourceError("VEP cache manifest is missing required fields")
    for field in ("archive_sha256", "tree_sha256"):
        if not isinstance(manifest[field], str) or not re.fullmatch(
            r"[0-9a-f]{64}", manifest[field]
        ):
            raise ResourceError(f"VEP cache manifest has invalid {field}")
    expected_values = {
        "resource_type": "ensembl_vep_indexed_cache",
        "vep_version": config.vep_expected_version,
        "cache_version": config.vep_cache_version,
        "assembly": VEP_CACHE_ASSEMBLY,
    }
    for field, expected in expected_values.items():
        if str(manifest[field]) != expected:
            raise ResourceError(
                f"VEP cache manifest {field} mismatch: {manifest[field]!r} != {expected!r}"
            )
    expected_root = (
        config.vep_cache_dir / "homo_sapiens" / (f"{config.vep_cache_version}_{VEP_CACHE_ASSEMBLY}")
    )
    if not expected_root.is_dir():
        raise ResourceError(f"VEP human GRCh38 cache directory not found: {expected_root}")
    if manifest["cache_metadata"] != _cache_metadata(config.vep_cache_dir):
        raise ResourceError("VEP cache semantic metadata differs from the resource manifest")
    if full_digest:
        statistics = _tree_statistics(config.vep_cache_dir)
        for field in ("tree_sha256", "tree_files", "tree_bytes"):
            if statistics[field] != manifest[field]:
                raise ResourceError(
                    f"VEP cache {field} mismatch: {statistics[field]} != {manifest[field]}"
                )
    return manifest


def sync_vep_cache(config: BuildConfig) -> dict[str, Any]:
    if config.vep_cache_version != VEP_CACHE_VERSION:
        raise ResourceError(
            f"pipeline 0.1 supports VEP cache {VEP_CACHE_VERSION}, got {config.vep_cache_version}"
        )
    if config.vep_expected_version != VEP_VERSION:
        raise ResourceError(
            f"pipeline 0.1 supports VEP {VEP_VERSION}, got {config.vep_expected_version}"
        )
    if config.vep_cache_dir.exists() or config.vep_cache_manifest.exists():
        return verify_vep_cache(config)

    archive = _archive_path(config)
    if not archive.is_file():
        _download_archive(config, archive)
    archive_identity = _verify_archive(archive)
    member_counts = _extract_archive(archive, config.vep_cache_dir)
    tree = _tree_statistics(config.vep_cache_dir)
    cache_metadata = _cache_metadata(config.vep_cache_dir)
    manifest: dict[str, Any] = {
        "resource_type": "ensembl_vep_indexed_cache",
        "species": "homo_sapiens",
        "vep_version": VEP_VERSION,
        "cache_version": VEP_CACHE_VERSION,
        "assembly": VEP_CACHE_ASSEMBLY,
        "source_url": VEP_CACHE_URL,
        "checksums_url": VEP_CACHE_CHECKSUMS_URL,
        "checksums_sha256": VEP_CACHE_CHECKSUMS_SHA256,
        "official_bsd_sum": VEP_CACHE_ARCHIVE_BSD_SUM,
        "archive_filename": VEP_CACHE_ARCHIVE,
        "archive_bytes": archive_identity["bytes"],
        "archive_sha256": archive_identity["sha256"],
        "installed_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        **member_counts,
        **tree,
        "cache_metadata": cache_metadata,
    }
    atomic_write_json(config.vep_cache_manifest, manifest)
    return verify_vep_cache(config)
