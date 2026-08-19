"""Freeze, prepare, and verify optional offline VEP predictor resources.

The public package contains only code and immutable release identities.  Large
resource objects, generated indexes, and their manifest live below the private
resource cache selected by the user.  A manifest is written only after every
requested object and derivative has passed its checksum and structural gates.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import sqlite3
import subprocess
import tarfile
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO

import pysam

from .config import BuildConfig, PredictorPluginConfig
from .constants import (
    ALPHAMISSENSE_BYTES,
    ALPHAMISSENSE_FILE,
    ALPHAMISSENSE_MD5,
    ALPHAMISSENSE_SHA256,
    ALPHAMISSENSE_URL,
    ALPHAMISSENSE_VERSION,
    CADD_INDEL_BYTES,
    CADD_INDEL_FILE,
    CADD_INDEL_MD5,
    CADD_INDEL_SHA256,
    CADD_INDEL_TBI_BYTES,
    CADD_INDEL_TBI_MD5,
    CADD_INDEL_TBI_SHA256,
    CADD_ROOT_URL,
    CADD_SNV_BYTES,
    CADD_SNV_FILE,
    CADD_SNV_MD5,
    CADD_SNV_SHA256,
    CADD_SNV_TBI_BYTES,
    CADD_SNV_TBI_MD5,
    CADD_SNV_TBI_SHA256,
    CADD_VERSION,
    LOFTEE_ANCESTOR_BYTES,
    LOFTEE_ANCESTOR_FAI_BYTES,
    LOFTEE_ANCESTOR_FAI_SHA256,
    LOFTEE_ANCESTOR_FILE,
    LOFTEE_ANCESTOR_GZI_BYTES,
    LOFTEE_ANCESTOR_GZI_SHA256,
    LOFTEE_ANCESTOR_SHA256,
    LOFTEE_COMMIT,
    LOFTEE_DATA_ROOT_URL,
    LOFTEE_GERP_BYTES,
    LOFTEE_GERP_FILE,
    LOFTEE_GERP_SHA256,
    LOFTEE_LOF_PM_ORIGINAL_SHA256,
    LOFTEE_LOF_PM_PATCHED_SHA256,
    LOFTEE_PATCH_FILE,
    LOFTEE_PATCH_ID,
    LOFTEE_PATCH_SHA256,
    LOFTEE_PATCHED_SOURCE_DIR,
    LOFTEE_SOURCE_ARCHIVE,
    LOFTEE_SOURCE_BYTES,
    LOFTEE_SOURCE_SHA256,
    LOFTEE_SOURCE_URL,
    LOFTEE_SQL_GZ_BYTES,
    LOFTEE_SQL_GZ_FILE,
    LOFTEE_SQL_GZ_SHA256,
    REVEL_ARCHIVE,
    REVEL_ARCHIVE_BYTES,
    REVEL_ARCHIVE_MD5,
    REVEL_ARCHIVE_SHA256,
    REVEL_MEMBER,
    REVEL_MEMBER_BYTES,
    REVEL_URL,
    REVEL_VERSION,
    SPLICEAI_SCOPE,
    SPLICEAI_SNV_BYTES,
    SPLICEAI_SNV_FILE,
    SPLICEAI_SNV_SHA256,
    SPLICEAI_SNV_TBI_BYTES,
    SPLICEAI_SNV_TBI_SHA256,
    SPLICEAI_SNV_URL,
    SPLICEAI_VERSION,
    VEP_VERSION,
)
from .errors import ResourceError
from .util import (
    CHUNK_SIZE,
    atomic_output_directory,
    atomic_write_json,
    directory_digest,
    fsync_path,
    sha256_file,
)

MANIFEST_SCHEMA_VERSION = 3
PLUGIN_ORDER = ("loftee", "cadd", "revel", "alphamissense", "spliceai")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
MD5_PATTERN = re.compile(r"[0-9a-f]{32}")


@dataclass(frozen=True)
class RemoteArtifact:
    key: str
    plugin: str
    relative_path: str
    url: str
    bytes: int
    sha256: str | None
    md5: str | None = None


def _remote_specs() -> dict[str, RemoteArtifact]:
    specs = (
        RemoteArtifact(
            "loftee_source_archive",
            "loftee",
            f"archives/loftee/{LOFTEE_COMMIT}/{LOFTEE_SOURCE_ARCHIVE}",
            LOFTEE_SOURCE_URL,
            LOFTEE_SOURCE_BYTES,
            LOFTEE_SOURCE_SHA256,
        ),
        RemoteArtifact(
            "loftee_gerp_bigwig",
            "loftee",
            f"archives/loftee/GRCh38/{LOFTEE_GERP_FILE}",
            f"{LOFTEE_DATA_ROOT_URL}/{LOFTEE_GERP_FILE}",
            LOFTEE_GERP_BYTES,
            LOFTEE_GERP_SHA256,
        ),
        RemoteArtifact(
            "loftee_human_ancestor",
            "loftee",
            f"archives/loftee/GRCh38/{LOFTEE_ANCESTOR_FILE}",
            f"{LOFTEE_DATA_ROOT_URL}/{LOFTEE_ANCESTOR_FILE}",
            LOFTEE_ANCESTOR_BYTES,
            LOFTEE_ANCESTOR_SHA256,
        ),
        RemoteArtifact(
            "loftee_human_ancestor_fai",
            "loftee",
            f"archives/loftee/GRCh38/{LOFTEE_ANCESTOR_FILE}.fai",
            f"{LOFTEE_DATA_ROOT_URL}/{LOFTEE_ANCESTOR_FILE}.fai",
            LOFTEE_ANCESTOR_FAI_BYTES,
            LOFTEE_ANCESTOR_FAI_SHA256,
        ),
        RemoteArtifact(
            "loftee_human_ancestor_gzi",
            "loftee",
            f"archives/loftee/GRCh38/{LOFTEE_ANCESTOR_FILE}.gzi",
            f"{LOFTEE_DATA_ROOT_URL}/{LOFTEE_ANCESTOR_FILE}.gzi",
            LOFTEE_ANCESTOR_GZI_BYTES,
            LOFTEE_ANCESTOR_GZI_SHA256,
        ),
        RemoteArtifact(
            "loftee_sql_gz",
            "loftee",
            f"archives/loftee/GRCh38/{LOFTEE_SQL_GZ_FILE}",
            f"{LOFTEE_DATA_ROOT_URL}/{LOFTEE_SQL_GZ_FILE}",
            LOFTEE_SQL_GZ_BYTES,
            LOFTEE_SQL_GZ_SHA256,
        ),
        RemoteArtifact(
            "cadd_snv",
            "cadd",
            f"archives/cadd/{CADD_VERSION}/{CADD_SNV_FILE}",
            f"{CADD_ROOT_URL}/{CADD_SNV_FILE}",
            CADD_SNV_BYTES,
            CADD_SNV_SHA256,
            CADD_SNV_MD5,
        ),
        RemoteArtifact(
            "cadd_snv_tbi",
            "cadd",
            f"archives/cadd/{CADD_VERSION}/{CADD_SNV_FILE}.tbi",
            f"{CADD_ROOT_URL}/{CADD_SNV_FILE}.tbi",
            CADD_SNV_TBI_BYTES,
            CADD_SNV_TBI_SHA256,
            CADD_SNV_TBI_MD5,
        ),
        RemoteArtifact(
            "cadd_indel",
            "cadd",
            f"archives/cadd/{CADD_VERSION}/{CADD_INDEL_FILE}",
            f"{CADD_ROOT_URL}/{CADD_INDEL_FILE}",
            CADD_INDEL_BYTES,
            CADD_INDEL_SHA256,
            CADD_INDEL_MD5,
        ),
        RemoteArtifact(
            "cadd_indel_tbi",
            "cadd",
            f"archives/cadd/{CADD_VERSION}/{CADD_INDEL_FILE}.tbi",
            f"{CADD_ROOT_URL}/{CADD_INDEL_FILE}.tbi",
            CADD_INDEL_TBI_BYTES,
            CADD_INDEL_TBI_SHA256,
            CADD_INDEL_TBI_MD5,
        ),
        RemoteArtifact(
            "revel_archive",
            "revel",
            f"archives/revel/{REVEL_VERSION}/{REVEL_ARCHIVE}",
            REVEL_URL,
            REVEL_ARCHIVE_BYTES,
            REVEL_ARCHIVE_SHA256,
            REVEL_ARCHIVE_MD5,
        ),
        RemoteArtifact(
            "alphamissense_scores",
            "alphamissense",
            f"archives/alphamissense/{ALPHAMISSENSE_VERSION}/{ALPHAMISSENSE_FILE}",
            ALPHAMISSENSE_URL,
            ALPHAMISSENSE_BYTES,
            ALPHAMISSENSE_SHA256,
            ALPHAMISSENSE_MD5,
        ),
        RemoteArtifact(
            "spliceai_snv",
            "spliceai",
            f"archives/spliceai/mane-v1.4/{SPLICEAI_SNV_FILE}",
            SPLICEAI_SNV_URL,
            SPLICEAI_SNV_BYTES,
            SPLICEAI_SNV_SHA256,
        ),
        RemoteArtifact(
            "spliceai_snv_tbi",
            "spliceai",
            f"archives/spliceai/mane-v1.4/{SPLICEAI_SNV_FILE}.tbi",
            f"{SPLICEAI_SNV_URL}.tbi",
            SPLICEAI_SNV_TBI_BYTES,
            SPLICEAI_SNV_TBI_SHA256,
        ),
    )
    return {spec.key: spec for spec in specs}


REMOTE_SPECS = _remote_specs()
EXPECTED_DERIVED_ARTIFACTS: dict[str, frozenset[str]] = {
    "loftee": frozenset({"loftee_source_tree", "loftee_sqlite"}),
    "cadd": frozenset(),
    "revel": frozenset({"revel_scores", "revel_scores_tbi"}),
    "alphamissense": frozenset({"alphamissense_scores_tbi"}),
    "spliceai": frozenset({"spliceai_indel_sentinel", "spliceai_indel_sentinel_tbi"}),
}


def _expected_artifact_keys(plugin: str) -> set[str]:
    return {key for key, spec in REMOTE_SPECS.items() if spec.plugin == plugin} | set(
        EXPECTED_DERIVED_ARTIFACTS[plugin]
    )


def _plugin_config(config: BuildConfig, plugin: str) -> PredictorPluginConfig:
    value = getattr(config.predictors, plugin, None)
    if not isinstance(value, PredictorPluginConfig):
        raise ResourceError(f"unknown predictor plugin: {plugin}")
    return value


def _safe_path(root: Path, relative: str) -> Path:
    relative_path = PurePosixPath(relative)
    if relative_path.is_absolute() or not relative_path.parts or ".." in relative_path.parts:
        raise ResourceError(f"unsafe predictor resource path: {relative!r}")
    candidate = root.joinpath(*relative_path.parts)
    if candidate.is_symlink():
        raise ResourceError(f"predictor resource cannot be a symlink: {relative}")
    try:
        candidate.resolve().relative_to(root.resolve())
    except ValueError as exc:
        raise ResourceError(f"predictor resource escapes root: {relative}") from exc
    return candidate


def _hashes(path: Path, include_md5: bool) -> tuple[str, str | None]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False) if include_md5 else None
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            sha256.update(chunk)
            if md5 is not None:
                md5.update(chunk)
    return sha256.hexdigest(), md5.hexdigest() if md5 is not None else None


def _require_frozen(spec: RemoteArtifact) -> str:
    if spec.sha256 is None or SHA256_PATTERN.fullmatch(spec.sha256) is None:
        raise ResourceError(
            f"predictor release {spec.key} has not yet been frozen to a complete SHA-256"
        )
    if spec.md5 is not None and MD5_PATTERN.fullmatch(spec.md5) is None:
        raise ResourceError(f"invalid pinned MD5 for {spec.key}")
    return spec.sha256


def _verify_remote(path: Path, spec: RemoteArtifact) -> dict[str, Any]:
    expected_sha256 = _require_frozen(spec)
    if not path.is_file() or path.is_symlink():
        raise ResourceError(f"predictor resource is missing or not a regular file: {path}")
    actual_bytes = path.stat().st_size
    if actual_bytes != spec.bytes:
        raise ResourceError(f"{spec.key} byte count mismatch: {actual_bytes} != {spec.bytes}")
    actual_sha256, actual_md5 = _hashes(path, include_md5=spec.md5 is not None)
    if actual_sha256 != expected_sha256:
        raise ResourceError(f"{spec.key} SHA-256 mismatch: {actual_sha256} != {expected_sha256}")
    if spec.md5 is not None and actual_md5 != spec.md5:
        raise ResourceError(f"{spec.key} upstream MD5 mismatch: {actual_md5} != {spec.md5}")
    return {
        "kind": "upstream",
        "plugin": spec.plugin,
        "path": spec.relative_path,
        "url": spec.url,
        "bytes": actual_bytes,
        "sha256": actual_sha256,
        "upstream_md5": actual_md5,
    }


def _download(config: BuildConfig, spec: RemoteArtifact) -> tuple[Path, dict[str, Any]]:
    _require_frozen(spec)
    root = config.predictors.root
    destination = _safe_path(root, spec.relative_path)
    if destination.is_file():
        return destination, _verify_remote(destination, spec)
    if not config.allow_download:
        raise ResourceError(
            f"{spec.key} is absent and resources.allow_download is false; place it at {destination}"
        )
    curl = shutil.which("curl")
    if curl is None:
        raise ResourceError("curl is required for restart-safe predictor downloads")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f".{destination.name}.partial")
    if partial.exists() and not partial.is_file():
        raise ResourceError(f"predictor partial is not a regular file: {partial}")
    if partial.is_file() and partial.stat().st_size > spec.bytes:
        raise ResourceError(
            f"predictor partial exceeds expected size and was retained for review: {partial}"
        )
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
            spec.url,
        ],
        check=False,
    )
    if completed.returncode != 0:
        raise ResourceError(
            f"download failed for {spec.key}; restartable partial retained at {partial}"
        )
    identity = _verify_remote(partial, spec)
    os.replace(partial, destination)
    fsync_path(destination.parent)
    return destination, identity


def _artifact(
    root: Path,
    relative: str,
    *,
    kind: str,
    plugin: str,
    source_sha256: str | None = None,
) -> dict[str, Any]:
    path = _safe_path(root, relative)
    if not path.is_file() or path.is_symlink():
        raise ResourceError(f"prepared predictor artifact missing: {path}")
    value: dict[str, Any] = {
        "kind": kind,
        "plugin": plugin,
        "path": relative,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if source_sha256 is not None:
        value["source_sha256"] = source_sha256
    return value


def _verify_manifest_artifact(
    root: Path, key: str, value: Any, *, full_digest: bool = True
) -> Path:
    if not isinstance(value, dict):
        raise ResourceError(f"predictor manifest artifact {key} must be an object")
    required = {"kind", "plugin", "path", "bytes", "sha256"}
    if not required.issubset(value):
        raise ResourceError(f"predictor manifest artifact {key} lacks required fields")
    if not isinstance(value["path"], str):
        raise ResourceError(f"predictor manifest artifact {key} path must be relative text")
    if not isinstance(value["bytes"], int) or value["bytes"] < 0:
        raise ResourceError(f"predictor manifest artifact {key} has invalid byte count")
    if not isinstance(value["sha256"], str) or not SHA256_PATTERN.fullmatch(value["sha256"]):
        raise ResourceError(f"predictor manifest artifact {key} has invalid SHA-256")
    path = _safe_path(root, value["path"])
    if not path.is_file() or path.stat().st_size != value["bytes"]:
        raise ResourceError(f"predictor artifact {key} is absent or has changed byte count")
    if full_digest and sha256_file(path) != value["sha256"]:
        raise ResourceError(f"predictor artifact {key} SHA-256 does not match manifest")
    return path


def _tree_artifact(root: Path, relative: str, plugin: str) -> dict[str, Any]:
    path = _safe_path(root, relative)
    if not path.is_dir() or path.is_symlink():
        raise ResourceError(f"prepared predictor tree missing: {path}")
    files = [item for item in path.rglob("*") if item.is_file()]
    if any(item.is_symlink() for item in path.rglob("*")):
        raise ResourceError(f"predictor source tree contains a symlink: {path}")
    return {
        "kind": "source_tree",
        "plugin": plugin,
        "path": relative,
        "tree_files": len(files),
        "tree_bytes": sum(item.stat().st_size for item in files),
        "tree_sha256": directory_digest(path),
    }


def _verify_tree_artifact(root: Path, key: str, value: Any, *, full_digest: bool = True) -> Path:
    if not isinstance(value, dict):
        raise ResourceError(f"predictor manifest tree {key} must be an object")
    required = {"kind", "plugin", "path", "tree_files", "tree_bytes", "tree_sha256"}
    if not required.issubset(value) or value.get("kind") != "source_tree":
        raise ResourceError(f"predictor manifest tree {key} lacks required fields")
    path = _safe_path(root, value["path"])
    if not path.is_dir() or path.is_symlink():
        raise ResourceError(f"predictor source tree {key} is missing")
    files = [item for item in path.rglob("*") if item.is_file()]
    if any(item.is_symlink() for item in path.rglob("*")):
        raise ResourceError(f"predictor source tree {key} contains a symlink")
    if len(files) != value["tree_files"]:
        raise ResourceError(f"predictor source tree {key} file count changed")
    if sum(item.stat().st_size for item in files) != value["tree_bytes"]:
        raise ResourceError(f"predictor source tree {key} byte count changed")
    if full_digest and directory_digest(path) != value["tree_sha256"]:
        raise ResourceError(f"predictor source tree {key} digest changed")
    return path


def _copy_stream(source: BinaryIO, destination: BinaryIO) -> int:
    total = 0
    while chunk := source.read(CHUNK_SIZE):
        destination.write(chunk)
        total += len(chunk)
    return total


def _loftee_patch_path() -> Path:
    path = Path(__file__).resolve().parent / "vendor_patches" / LOFTEE_PATCH_FILE
    if not path.is_file() or path.is_symlink():
        raise ResourceError(f"pinned LOFTEE patch is missing or unsafe: {path}")
    if sha256_file(path) != LOFTEE_PATCH_SHA256:
        raise ResourceError("pinned LOFTEE patch SHA-256 does not match the release boundary")
    return path


def _loftee_patch_identity(source_tree: Path) -> dict[str, Any]:
    """Validate the exact patched source and return portable derivation metadata."""

    _loftee_patch_path()
    lof_pm = source_tree / "LoF.pm"
    if not lof_pm.is_file() or lof_pm.is_symlink():
        raise ResourceError("patched LOFTEE source lacks a safe LoF.pm")
    actual = sha256_file(lof_pm)
    if actual != LOFTEE_LOF_PM_PATCHED_SHA256:
        raise ResourceError("patched LOFTEE LoF.pm SHA-256 differs from the pinned transformation")
    try:
        text = lof_pm.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise ResourceError(f"cannot inspect patched LOFTEE LoF.pm: {exc}") from exc
    if "use Bio::DB::HTS::Faidx;" not in text or "Bio::DB::HTS::Faidx->new" not in text:
        raise ResourceError("patched LOFTEE LoF.pm lacks Bio::DB::HTS::Faidx lookup")
    if "`samtools faidx" in text:
        raise ResourceError("patched LOFTEE LoF.pm still invokes external samtools faidx")
    return {
        "source_sha256": LOFTEE_SOURCE_SHA256,
        "source_archive_sha256": LOFTEE_SOURCE_SHA256,
        "patch_id": LOFTEE_PATCH_ID,
        "patch_file": f"vendor_patches/{LOFTEE_PATCH_FILE}",
        "patch_sha256": LOFTEE_PATCH_SHA256,
        "original_lof_pm_sha256": LOFTEE_LOF_PM_ORIGINAL_SHA256,
        "patched_lof_pm_sha256": LOFTEE_LOF_PM_PATCHED_SHA256,
        "patch_command_semantics": [
            "patch",
            "--batch",
            "--forward",
            "--fuzz=0",
            "--no-backup-if-mismatch",
            "--directory",
            "SOURCE_TREE",
            "--strip=1",
            "--input",
            f"PACKAGE/{LOFTEE_PATCH_FILE}",
        ],
    }


def _apply_loftee_patch(source_tree: Path) -> dict[str, Any]:
    """Apply the pinned LOFTEE patch only to the exact expected upstream bytes."""

    lof_pm = source_tree / "LoF.pm"
    if not lof_pm.is_file() or lof_pm.is_symlink():
        raise ResourceError("LOFTEE source lacks a safe LoF.pm before patching")
    actual = sha256_file(lof_pm)
    if actual != LOFTEE_LOF_PM_ORIGINAL_SHA256:
        raise ResourceError(
            "upstream LOFTEE LoF.pm differs from the exact patch preimage: "
            f"expected={LOFTEE_LOF_PM_ORIGINAL_SHA256}, observed={actual}"
        )
    patch_path = _loftee_patch_path()
    executable = shutil.which("patch")
    if executable is None:
        raise ResourceError(
            "GNU patch is required to prepare LOFTEE; use the pinned workflow environment"
        )
    command = [
        executable,
        "--batch",
        "--forward",
        "--fuzz=0",
        "--no-backup-if-mismatch",
        "--directory",
        str(source_tree),
        "--strip=1",
        "--input",
        str(patch_path),
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout)[-4000:]
        raise ResourceError(
            f"pinned LOFTEE patch failed with exit {completed.returncode}: {detail}"
        )
    residue = sorted(
        item.relative_to(source_tree).as_posix()
        for item in source_tree.rglob("*")
        if item.is_file() and item.suffix in {".orig", ".rej"}
    )
    if residue:
        raise ResourceError("LOFTEE patch left backup/reject files: " + ", ".join(residue))
    return _loftee_patch_identity(source_tree)


def _extract_loftee_source(root: Path, archive: Path) -> tuple[str, dict[str, Any]]:
    relative = f"prepared/loftee/{LOFTEE_COMMIT}/{LOFTEE_PATCHED_SOURCE_DIR}"
    destination = _safe_path(root, relative)
    if not destination.exists():
        expected_root = f"loftee-{LOFTEE_COMMIT}"
        try:
            with tarfile.open(archive, mode="r:gz") as handle:
                members = handle.getmembers()
                if not members:
                    raise ResourceError("LOFTEE source archive is empty")
                with atomic_output_directory(destination) as staging:
                    files = 0
                    for member in members:
                        member_path = PurePosixPath(member.name)
                        if (
                            member_path.is_absolute()
                            or ".." in member_path.parts
                            or not member_path.parts
                            or member_path.parts[0] != expected_root
                        ):
                            raise ResourceError(
                                f"unsafe LOFTEE source archive member: {member.name}"
                            )
                        relative_member = PurePosixPath(*member_path.parts[1:])
                        if not relative_member.parts:
                            if not member.isdir():
                                raise ResourceError(
                                    f"invalid LOFTEE archive root member: {member.name}"
                                )
                            continue
                        target = staging.joinpath(*relative_member.parts)
                        if member.isdir():
                            target.mkdir(parents=True, exist_ok=True)
                            continue
                        if not member.isfile():
                            raise ResourceError(
                                f"unsupported LOFTEE archive member type: {member.name}"
                            )
                        source = handle.extractfile(member)
                        if source is None:
                            raise ResourceError(f"cannot read LOFTEE archive member: {member.name}")
                        target.parent.mkdir(parents=True, exist_ok=True)
                        with source, target.open("wb") as output:
                            copied = _copy_stream(source, output)
                        if copied != member.size:
                            raise ResourceError(f"truncated LOFTEE archive member: {member.name}")
                        os.chmod(target, 0o755 if member.mode & 0o111 else 0o644)
                        files += 1
                    if files == 0:
                        raise ResourceError("LOFTEE source archive contains no files")
                    _apply_loftee_patch(staging)
        except (OSError, tarfile.TarError) as exc:
            raise ResourceError(f"invalid LOFTEE source archive: {exc}") from exc
    required = (
        "LoF.pm",
        "LICENSE",
        "maxEntScan/score3.pl",
        "maxEntScan/splicemodels/me2x5",
        "splice_data/donor_model.txt",
        "splice_data/acceptor_model.txt",
    )
    missing = [item for item in required if not (destination / item).is_file()]
    if missing:
        raise ResourceError("LOFTEE source tree lacks required files: " + ", ".join(missing))
    patch_identity = _loftee_patch_identity(destination)
    artifact = _tree_artifact(root, relative, "loftee")
    artifact.update(patch_identity)
    return relative, artifact


def _sqlite_identity(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise ResourceError(f"LOFTEE conservation SQLite file is missing: {path}")
    with path.open("rb") as handle:
        if handle.read(16) != b"SQLite format 3\x00":
            raise ResourceError("LOFTEE conservation file lacks SQLite magic")
    try:
        connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        try:
            quick_check = connection.execute("PRAGMA quick_check").fetchone()
            if quick_check != ("ok",):
                raise ResourceError(f"LOFTEE SQLite quick_check failed: {quick_check!r}")
            tables = {
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
                )
            }
            if "phylocsf_summary" not in tables:
                raise ResourceError("LOFTEE SQLite lacks phylocsf_summary")
            rows = int(connection.execute("SELECT COUNT(*) FROM phylocsf_summary").fetchone()[0])
            if rows < 1:
                raise ResourceError("LOFTEE SQLite phylocsf_summary is empty")
        finally:
            connection.close()
    except sqlite3.Error as exc:
        raise ResourceError(f"cannot validate LOFTEE SQLite: {exc}") from exc
    return {
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "tables": sorted(tables),
        "phylocsf_summary_rows": rows,
    }


def _prepare_loftee_sql(root: Path, compressed: Path) -> tuple[str, dict[str, Any]]:
    relative = "prepared/loftee/GRCh38/loftee.sqlite"
    destination = _safe_path(root, relative)
    if not destination.is_file():
        destination.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{destination.name}.partial.", dir=destination.parent
        )
        temporary = Path(temporary_name)
        try:
            with os.fdopen(descriptor, "wb") as output, gzip.open(compressed, "rb") as source:
                _copy_stream(source, output)
                output.flush()
                os.fsync(output.fileno())
            _sqlite_identity(temporary)
            os.replace(temporary, destination)
            fsync_path(destination.parent)
        except BaseException:
            # Preserve a non-empty failed derivative for forensic review.
            if temporary.exists() and temporary.stat().st_size == 0:
                temporary.unlink()
            raise
    return relative, _sqlite_identity(destination)


def _validate_loftee_ancillary(ancestor: Path, gerp: Path) -> dict[str, Any]:
    with gerp.open("rb") as handle:
        magic = handle.read(4).hex()
    if magic not in {"26fc8f88", "888ffc26"}:
        raise ResourceError(f"LOFTEE GERP file has invalid bigWig magic: {magic}")
    try:
        fasta = pysam.FastaFile(str(ancestor))
        try:
            references = tuple(fasta.references)
            if len(references) < 24:
                raise ResourceError(f"LOFTEE ancestor FASTA has only {len(references)} references")
            first = references[0]
            if not fasta.fetch(first, 0, min(32, fasta.get_reference_length(first))):
                raise ResourceError("LOFTEE ancestor FASTA first reference is empty")
        finally:
            fasta.close()
    except (OSError, ValueError) as exc:
        raise ResourceError(f"cannot open indexed LOFTEE ancestor FASTA: {exc}") from exc
    return {
        "bigwig_magic": magic,
        "ancestor_references": len(references),
        "ancestor_first_reference": first,
    }


REVEL_HEADER = (
    b"chr",
    b"hg19_pos",
    b"grch38_pos",
    b"ref",
    b"alt",
    b"aaref",
    b"aaalt",
    b"REVEL",
    b"Ensembl_transcriptid",
)


def _load_preparation(
    path: Path, source_sha256: str, *, full_digest: bool = True
) -> dict[str, Any] | None:
    report_path = path / "preparation.json"
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except json.JSONDecodeError as exc:
        raise ResourceError(f"invalid predictor preparation report: {report_path}: {exc}") from exc
    if not isinstance(report, dict) or report.get("source_sha256") != source_sha256:
        raise ResourceError(f"predictor preparation source mismatch: {report_path}")
    outputs = report.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ResourceError(f"predictor preparation report has no outputs: {report_path}")
    for relative, identity in outputs.items():
        if not isinstance(relative, str) or not isinstance(identity, dict):
            raise ResourceError(f"invalid predictor preparation output: {report_path}")
        output = path / relative
        if (
            not output.is_file()
            or output.stat().st_size != identity.get("bytes")
            or (full_digest and sha256_file(output) != identity.get("sha256"))
        ):
            raise ResourceError(f"prepared predictor output changed: {output}")
    return report


def _validate_preparation_metadata(
    root: Path,
    artifacts: dict[str, Any],
    preparations: Any,
    plugins: set[str],
    *,
    full_digest: bool,
) -> None:
    expected_plugins = plugins & {"loftee", "revel", "alphamissense", "spliceai"}
    if not isinstance(preparations, dict) or set(preparations) != expected_plugins:
        raise ResourceError("predictor preparation set does not match installed plugins")
    if "loftee" in plugins:
        preparation = preparations["loftee"]
        source_path = _safe_path(root, artifacts["loftee_source_tree"]["path"])
        expected = {
            "source_tree": artifacts["loftee_source_tree"]["path"],
            "source_derivation": _loftee_patch_identity(source_path),
            "sqlite": artifacts["loftee_sqlite"]["path"],
            "sqlite_identity": _sqlite_identity(
                _safe_path(root, artifacts["loftee_sqlite"]["path"])
            ),
        }
        if preparation != expected:
            raise ResourceError("LOFTEE preparation metadata differs from validated outputs")
    if "revel" in plugins:
        preparation = preparations["revel"]
        report = _load_preparation(
            _safe_path(root, f"prepared/revel/{REVEL_VERSION}"),
            REVEL_ARCHIVE_SHA256,
            full_digest=full_digest,
        )
        if report is None or preparation != report:
            raise ResourceError("REVEL preparation report differs from its manifest record")
        source_rows = report.get("source_rows")
        retained_rows = report.get("retained_grch38_rows")
        omitted_rows = report.get("omitted_missing_grch38_position")
        if (
            report.get("status") != "PASS"
            or not isinstance(source_rows, int)
            or not isinstance(retained_rows, int)
            or not isinstance(omitted_rows, int)
            or retained_rows < 1
            or omitted_rows < 0
            or retained_rows + omitted_rows != source_rows
        ):
            raise ResourceError("REVEL preparation row accounting is invalid")
        for key in ("revel_scores", "revel_scores_tbi"):
            artifact = artifacts[key]
            output = report["outputs"].get(Path(artifact["path"]).name)
            if not isinstance(output, dict) or output != {
                "bytes": artifact["bytes"],
                "sha256": artifact["sha256"],
            }:
                raise ResourceError(f"REVEL preparation output identity mismatch for {key}")
    if "alphamissense" in plugins:
        preparation = preparations["alphamissense"]
        expected_fields = {
            "source_sha256": ALPHAMISSENSE_SHA256,
            "index_sha256": artifacts["alphamissense_scores_tbi"]["sha256"],
            "index_semantics": "tabix seq=1 start=2 end=2 one-based",
        }
        if (
            not isinstance(preparation, dict)
            or any(preparation.get(key) != value for key, value in expected_fields.items())
            or not isinstance(preparation.get("pysam_version"), str)
        ):
            raise ResourceError("AlphaMissense preparation metadata is invalid")
    if "spliceai" in plugins:
        preparation = preparations["spliceai"]
        report = _load_preparation(
            _safe_path(root, "prepared/spliceai/mane-v1.4"),
            _require_frozen(REMOTE_SPECS["spliceai_snv"]),
            full_digest=full_digest,
        )
        if report is None or preparation != report:
            raise ResourceError("SpliceAI sentinel report differs from its manifest record")
        if (
            report.get("status") != "PASS"
            or report.get("record_count") != 0
            or report.get("semantic_scope") != "empty operational sentinel; no indel score coverage"
        ):
            raise ResourceError("SpliceAI sentinel preparation semantics are invalid")
        for key in ("spliceai_indel_sentinel", "spliceai_indel_sentinel_tbi"):
            artifact = artifacts[key]
            output = report["outputs"].get(Path(artifact["path"]).name)
            if not isinstance(output, dict) or output != {
                "bytes": artifact["bytes"],
                "sha256": artifact["sha256"],
            }:
                raise ResourceError(f"SpliceAI preparation output identity mismatch for {key}")


def _prepare_revel(root: Path, archive_path: Path) -> tuple[dict[str, str], dict[str, Any]]:
    relative_dir = f"prepared/revel/{REVEL_VERSION}"
    destination = _safe_path(root, relative_dir)
    source_sha256 = REVEL_ARCHIVE_SHA256
    if destination.exists():
        report = _load_preparation(destination, source_sha256)
        if report is None:
            raise ResourceError(
                f"unmanifested REVEL preparation exists and was retained for review: {destination}"
            )
        return {
            "scores": f"{relative_dir}/revel.GRCh38.tsv.gz",
            "index": f"{relative_dir}/revel.GRCh38.tsv.gz.tbi",
        }, report

    sort_executable = shutil.which("sort")
    if sort_executable is None:
        raise ResourceError("GNU/coreutils sort is required to prepare REVEL GRCh38 data")
    expected_zip_root = PurePosixPath(REVEL_MEMBER)
    with atomic_output_directory(destination) as staging:
        unsorted = staging / "revel.GRCh38.unsorted.tsv"
        sorted_path = staging / "revel.GRCh38.sorted.tsv"
        sort_temp = staging / "sort-tmp"
        sort_temp.mkdir()
        source_rows = 0
        retained_rows = 0
        missing_grch38 = 0
        try:
            with zipfile.ZipFile(archive_path) as archive:
                members = archive.infolist()
                if len(members) != 1:
                    raise ResourceError(f"REVEL ZIP contains {len(members)} members, expected one")
                member = members[0]
                member_path = PurePosixPath(member.filename)
                if (
                    member.is_dir()
                    or member_path != expected_zip_root
                    or member.file_size != REVEL_MEMBER_BYTES
                ):
                    raise ResourceError(f"unexpected REVEL ZIP member: {member.filename}")
                with archive.open(member) as source, unsorted.open("wb") as output:
                    header = source.readline().rstrip(b"\r\n").split(b",")
                    if tuple(header) != REVEL_HEADER:
                        raise ResourceError(f"unexpected REVEL header: {header!r}")
                    for line_number, line in enumerate(source, start=2):
                        source_rows += 1
                        fields = line.rstrip(b"\r\n").split(b",")
                        if len(fields) != len(REVEL_HEADER):
                            raise ResourceError(
                                f"REVEL row {line_number} has {len(fields)} columns"
                            )
                        if fields[2] == b".":
                            missing_grch38 += 1
                            continue
                        if fields[0] not in {
                            *(str(value).encode() for value in range(1, 23)),
                            b"X",
                            b"Y",
                        }:
                            raise ResourceError(f"REVEL row {line_number} has invalid chromosome")
                        try:
                            position = int(fields[2])
                        except ValueError as exc:
                            raise ResourceError(
                                f"REVEL row {line_number} has invalid GRCh38 position"
                            ) from exc
                        if position < 1:
                            raise ResourceError(
                                f"REVEL row {line_number} has non-positive GRCh38 position"
                            )
                        output.write(b"\t".join(fields) + b"\n")
                        retained_rows += 1
        except (OSError, zipfile.BadZipFile) as exc:
            raise ResourceError(f"cannot read REVEL ZIP: {exc}") from exc

        completed = subprocess.run(
            [
                sort_executable,
                "--stable",
                "--field-separator=\t",
                "--key=1,1",
                "--key=3,3n",
                f"--temporary-directory={sort_temp}",
                f"--output={sorted_path}",
                str(unsorted),
            ],
            check=False,
            capture_output=True,
            text=True,
            env={**os.environ, "LC_ALL": "C", "TMPDIR": str(sort_temp)},
        )
        if completed.returncode != 0:
            raise ResourceError(f"REVEL external sort failed: {completed.stderr[-4000:]}")

        output_path = staging / "revel.GRCh38.tsv.gz"
        previous_key: tuple[bytes, int] | None = None
        written_rows = 0
        with sorted_path.open("rb") as source, pysam.BGZFile(str(output_path), "wb") as output:
            output.write(b"#" + b"\t".join(REVEL_HEADER) + b"\n")
            for line_number, line in enumerate(source, start=2):
                fields = line.rstrip(b"\n").split(b"\t")
                key = (fields[0], int(fields[2]))
                if previous_key is not None and key < previous_key:
                    raise ResourceError(f"REVEL sorted output regressed at row {line_number}")
                previous_key = key
                output.write(line)
                written_rows += 1
        if written_rows != retained_rows:
            raise ResourceError(f"REVEL row accounting mismatch: {written_rows} != {retained_rows}")
        index_path = staging / "revel.GRCh38.tsv.gz.tbi"
        try:
            pysam.tabix_index(
                str(output_path),
                seq_col=0,
                start_col=2,
                end_col=2,
                meta_char="#",
                zerobased=False,
                force=True,
                index=str(index_path),
            )
        except (OSError, ValueError) as exc:
            raise ResourceError(f"cannot index prepared REVEL data: {exc}") from exc
        for temporary in (unsorted, sorted_path):
            temporary.unlink()
        sort_temp.rmdir()
        output_identities = {
            output_path.name: {
                "bytes": output_path.stat().st_size,
                "sha256": sha256_file(output_path),
            },
            index_path.name: {
                "bytes": index_path.stat().st_size,
                "sha256": sha256_file(index_path),
            },
        }
        sort_version = subprocess.run(
            [sort_executable, "--version"],
            check=False,
            capture_output=True,
            text=True,
        ).stdout.splitlines()[0]
        report = {
            "status": "PASS",
            "source_sha256": source_sha256,
            "source_member": REVEL_MEMBER,
            "source_member_bytes": REVEL_MEMBER_BYTES,
            "source_rows": source_rows,
            "retained_grch38_rows": retained_rows,
            "omitted_missing_grch38_position": missing_grch38,
            "sort_semantics": "LC_ALL=C stable key chromosome lexical, GRCh38 position numeric",
            "sort_version": sort_version,
            "pysam_version": pysam.__version__,
            "outputs": output_identities,
        }
        atomic_write_json(staging / "preparation.json", report)
    return {
        "scores": f"{relative_dir}/revel.GRCh38.tsv.gz",
        "index": f"{relative_dir}/revel.GRCh38.tsv.gz.tbi",
    }, report


def _prepare_tabix_index(
    source: Path,
    *,
    seq_col: int,
    start_col: int,
    end_col: int,
) -> Path:
    index = Path(f"{source}.tbi")
    if index.is_file():
        return index
    temporary = index.with_name(f".{index.name}.partial.{os.getpid()}")
    try:
        pysam.tabix_index(
            str(source),
            seq_col=seq_col,
            start_col=start_col,
            end_col=end_col,
            meta_char="#",
            zerobased=False,
            force=True,
            index=str(temporary),
        )
        if not temporary.is_file():
            raise ResourceError(f"tabix did not create expected index: {temporary}")
        os.replace(temporary, index)
        fsync_path(index.parent)
    except (OSError, ValueError) as exc:
        raise ResourceError(f"cannot index predictor resource {source}: {exc}") from exc
    return index


def _prepare_spliceai_sentinel(root: Path) -> tuple[dict[str, str], dict[str, Any]]:
    relative_dir = "prepared/spliceai/mane-v1.4"
    destination = _safe_path(root, relative_dir)
    source_sha256 = SPLICEAI_SNV_SHA256
    if source_sha256 is None:
        raise ResourceError("SpliceAI SNV source has not yet been frozen to SHA-256")
    if destination.exists():
        report = _load_preparation(destination, source_sha256)
        if report is None:
            raise ResourceError(
                f"unmanifested SpliceAI preparation retained for review: {destination}"
            )
        return {
            "indel_sentinel": f"{relative_dir}/spliceai.indels.NOT_AVAILABLE.vcf.gz",
            "indel_sentinel_index": (f"{relative_dir}/spliceai.indels.NOT_AVAILABLE.vcf.gz.tbi"),
        }, report
    with atomic_output_directory(destination) as staging:
        sentinel = staging / "spliceai.indels.NOT_AVAILABLE.vcf.gz"
        payload = (
            b"##fileformat=VCFv4.2\n"
            b"##reference=GRCh38\n"
            b"##GAB_resource_scope=EMPTY_SENTINEL_NO_SPLICEAI_INDEL_SCORES\n"
            b'##INFO=<ID=SpliceAI,Number=.,Type=String,Description="Empty sentinel; '
            b'no indel scores are present">\n'
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        with pysam.BGZFile(str(sentinel), "wb") as handle:
            handle.write(payload)
        index = staging / "spliceai.indels.NOT_AVAILABLE.vcf.gz.tbi"
        try:
            pysam.tabix_index(str(sentinel), preset="vcf", force=True, index=str(index))
        except (OSError, ValueError) as exc:
            raise ResourceError(f"cannot index empty SpliceAI indel sentinel: {exc}") from exc
        outputs = {
            sentinel.name: {
                "bytes": sentinel.stat().st_size,
                "sha256": sha256_file(sentinel),
            },
            index.name: {"bytes": index.stat().st_size, "sha256": sha256_file(index)},
        }
        report = {
            "status": "PASS",
            "source_sha256": source_sha256,
            "semantic_scope": "empty operational sentinel; no indel score coverage",
            "record_count": 0,
            "pysam_version": pysam.__version__,
            "outputs": outputs,
        }
        atomic_write_json(staging / "preparation.json", report)
    return {
        "indel_sentinel": f"{relative_dir}/spliceai.indels.NOT_AVAILABLE.vcf.gz",
        "indel_sentinel_index": f"{relative_dir}/spliceai.indels.NOT_AVAILABLE.vcf.gz.tbi",
    }, report


def _tabix_identity(
    path: Path,
    index: Path,
    *,
    required_header_tokens: tuple[str, ...],
    allow_empty: bool = False,
) -> dict[str, Any]:
    try:
        tabix = pysam.TabixFile(str(path), index=str(index))
        try:
            header = tuple(tabix.header)
            header_text = "\n".join(header)
            for token in required_header_tokens:
                if token not in header_text:
                    raise ResourceError(f"{path.name} header lacks {token!r}")
            contigs = tuple(tabix.contigs)
            first_record: str | None = None
            for contig in contigs:
                iterator = tabix.fetch(contig)
                first_record = next(iterator, None)
                if first_record is not None:
                    break
            if not allow_empty and first_record is None:
                raise ResourceError(f"tabix resource contains no records: {path}")
            if allow_empty and first_record is not None:
                raise ResourceError(f"empty sentinel unexpectedly contains a record: {path}")
        finally:
            tabix.close()
    except (OSError, ValueError) as exc:
        raise ResourceError(f"cannot open tabix resource {path}: {exc}") from exc
    return {
        "header_lines": len(header),
        "contigs": len(contigs),
        "first_record_present": first_record is not None,
    }


def _validate_plugin_structures(
    root: Path, artifacts: dict[str, Any], plugins: set[str]
) -> dict[str, Any]:
    paths = {
        key: _safe_path(root, value["path"])
        for key, value in artifacts.items()
        if isinstance(value, dict) and isinstance(value.get("path"), str)
    }
    structures: dict[str, Any] = {}
    if "loftee" in plugins:
        source = paths["loftee_source_tree"]
        for required in ("LoF.pm", "extended_splice.pl", "ancestral.pm"):
            if not (source / required).is_file():
                raise ResourceError(f"LOFTEE source tree lacks {required}")
        source_derivation = _loftee_patch_identity(source)
        tree_artifact = artifacts["loftee_source_tree"]
        for field, expected in source_derivation.items():
            if tree_artifact.get(field) != expected:
                raise ResourceError(f"LOFTEE source-tree derivation metadata mismatch for {field}")
        ancillary = _validate_loftee_ancillary(
            paths["loftee_human_ancestor"], paths["loftee_gerp_bigwig"]
        )
        sqlite_identity = _sqlite_identity(paths["loftee_sqlite"])
        structures["loftee"] = {
            "source_derivation": source_derivation,
            "ancestral_lookup_runtime": "Bio::DB::HTS::Faidx in pinned VEP container",
            **ancillary,
            "sqlite": sqlite_identity,
        }
    if "cadd" in plugins:
        structures["cadd"] = {
            "snv": _tabix_identity(
                paths["cadd_snv"],
                paths["cadd_snv_tbi"],
                required_header_tokens=("GRCh38", "#Chrom"),
            ),
            "indel": _tabix_identity(
                paths["cadd_indel"],
                paths["cadd_indel_tbi"],
                required_header_tokens=("GRCh38", "#Chrom"),
            ),
        }
    if "revel" in plugins:
        structures["revel"] = _tabix_identity(
            paths["revel_scores"],
            paths["revel_scores_tbi"],
            required_header_tokens=("#chr", "grch38_pos", "Ensembl_transcriptid"),
        )
    if "alphamissense" in plugins:
        structures["alphamissense"] = _tabix_identity(
            paths["alphamissense_scores"],
            paths["alphamissense_scores_tbi"],
            required_header_tokens=("#CHROM", "am_pathogenicity", "transcript_id"),
        )
    if "spliceai" in plugins:
        structures["spliceai"] = {
            "snv": _tabix_identity(
                paths["spliceai_snv"],
                paths["spliceai_snv_tbi"],
                required_header_tokens=(
                    "SpliceAIv1.3.1",
                    "fileDate=20251015",
                    "reference=GRCh38/hg38",
                ),
            ),
            "indel_sentinel": _tabix_identity(
                paths["spliceai_indel_sentinel"],
                paths["spliceai_indel_sentinel_tbi"],
                required_header_tokens=(
                    "EMPTY_SENTINEL_NO_SPLICEAI_INDEL_SCORES",
                    "no indel scores are present",
                ),
                allow_empty=True,
            ),
        }
    return structures


def _plugin_metadata(config: BuildConfig, plugin: str, artifact_keys: list[str]) -> dict[str, Any]:
    plugin_config = _plugin_config(config, plugin)
    versions = {
        "loftee": LOFTEE_COMMIT,
        "cadd": CADD_VERSION,
        "revel": REVEL_VERSION,
        "alphamissense": ALPHAMISSENSE_VERSION,
        "spliceai": SPLICEAI_VERSION,
    }
    coverage = {
        "loftee": {
            "variant_types": ["predicted loss-of-function"],
            "assembly": "GRCh38",
            "transcript_scope": "Ensembl VEP transcript consequences",
            "ancestral_allele_filter": True,
            "ancestral_lookup_runtime": "Bio::DB::HTS::Faidx in pinned VEP container",
            "source_patch_id": LOFTEE_PATCH_ID,
            "source_patch_sha256": LOFTEE_PATCH_SHA256,
        },
        "cadd": {
            "variant_types": ["SNV", "small indel"],
            "assembly": "GRCh38",
            "pre_scored_only": True,
        },
        "revel": {
            "variant_types": ["missense SNV"],
            "assembly": "GRCh38",
            "transcript_matching": True,
        },
        "alphamissense": {
            "variant_types": ["missense SNV"],
            "assembly": "GRCh38",
            "transcript_matching": True,
        },
        "spliceai": {
            "variant_types": ["SNV"],
            "assembly": "GRCh38",
            "scope": SPLICEAI_SCOPE,
            "indel_scores_available": False,
            "indel_runtime_input": "explicit empty sentinel",
        },
    }
    return {
        "version": versions[plugin],
        "license": plugin_config.license,
        "resource_redistribution_allowed": plugin_config.redistribution_allowed,
        "coverage": coverage[plugin],
        "artifact_keys": sorted(artifact_keys),
    }


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ResourceError(f"predictor resource manifest not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ResourceError(f"invalid predictor resource manifest: {exc}") from exc
    if not isinstance(value, dict):
        raise ResourceError("predictor resource manifest must be an object")
    return value


def _manifest_document(
    config: BuildConfig,
    requested: set[str],
    artifacts: dict[str, Any],
    preparations: dict[str, Any],
    structures: dict[str, Any],
) -> dict[str, Any]:
    plugin_documents = {
        plugin: _plugin_metadata(
            config,
            plugin,
            [key for key, value in artifacts.items() if value.get("plugin") == plugin],
        )
        for plugin in PLUGIN_ORDER
        if plugin in requested
    }
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "resource_type": "vep_predictor_resources",
        "assembly": "GRCh38",
        "vep_version": VEP_VERSION,
        "installed_plugins": [plugin for plugin in PLUGIN_ORDER if plugin in requested],
        "plugins": plugin_documents,
        "artifacts": dict(sorted(artifacts.items())),
        "preparations": preparations,
        "structural_validation": structures,
    }


def sync_predictor_resources(config: BuildConfig) -> dict[str, Any]:
    requested = set(config.predictors.enabled_plugins)
    if not requested:
        return {"status": "not_enabled", "enabled_plugins": []}
    if config.predictors.manifest.is_file():
        result = verify_predictor_resources(config)
        return {**result, "reused": True}

    root = config.predictors.root
    root.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, Any] = {}
    local_paths: dict[str, Path] = {}
    for key, spec in REMOTE_SPECS.items():
        if spec.plugin not in requested:
            continue
        path, identity = _download(config, spec)
        local_paths[key] = path
        artifacts[key] = identity

    preparations: dict[str, Any] = {}
    if "loftee" in requested:
        source_relative, source_artifact = _extract_loftee_source(
            root, local_paths["loftee_source_archive"]
        )
        artifacts["loftee_source_tree"] = source_artifact
        sqlite_relative, sqlite_identity = _prepare_loftee_sql(root, local_paths["loftee_sql_gz"])
        artifacts["loftee_sqlite"] = _artifact(
            root,
            sqlite_relative,
            kind="derived",
            plugin="loftee",
            source_sha256=LOFTEE_SQL_GZ_SHA256,
        )
        preparations["loftee"] = {
            "source_tree": source_relative,
            "source_derivation": _loftee_patch_identity(_safe_path(root, source_relative)),
            "sqlite": sqlite_relative,
            "sqlite_identity": sqlite_identity,
        }
    if "revel" in requested:
        revel_paths, report = _prepare_revel(root, local_paths["revel_archive"])
        artifacts["revel_scores"] = _artifact(
            root,
            revel_paths["scores"],
            kind="derived",
            plugin="revel",
            source_sha256=REVEL_ARCHIVE_SHA256,
        )
        artifacts["revel_scores_tbi"] = _artifact(
            root,
            revel_paths["index"],
            kind="derived_index",
            plugin="revel",
            source_sha256=artifacts["revel_scores"]["sha256"],
        )
        preparations["revel"] = report
    if "alphamissense" in requested:
        alpha = local_paths["alphamissense_scores"]
        alpha_index = _prepare_tabix_index(alpha, seq_col=0, start_col=1, end_col=1)
        alpha_index_relative = alpha_index.relative_to(root).as_posix()
        artifacts["alphamissense_scores_tbi"] = _artifact(
            root,
            alpha_index_relative,
            kind="derived_index",
            plugin="alphamissense",
            source_sha256=ALPHAMISSENSE_SHA256,
        )
        preparations["alphamissense"] = {
            "source_sha256": ALPHAMISSENSE_SHA256,
            "index_sha256": artifacts["alphamissense_scores_tbi"]["sha256"],
            "index_semantics": "tabix seq=1 start=2 end=2 one-based",
            "pysam_version": pysam.__version__,
        }
    if "spliceai" in requested:
        sentinel_paths, report = _prepare_spliceai_sentinel(root)
        artifacts["spliceai_indel_sentinel"] = _artifact(
            root,
            sentinel_paths["indel_sentinel"],
            kind="empty_sentinel",
            plugin="spliceai",
            source_sha256=SPLICEAI_SNV_SHA256,
        )
        artifacts["spliceai_indel_sentinel_tbi"] = _artifact(
            root,
            sentinel_paths["indel_sentinel_index"],
            kind="derived_index",
            plugin="spliceai",
            source_sha256=artifacts["spliceai_indel_sentinel"]["sha256"],
        )
        preparations["spliceai"] = report

    structures = _validate_plugin_structures(root, artifacts, requested)
    manifest = _manifest_document(config, requested, artifacts, preparations, structures)
    atomic_write_json(config.predictors.manifest, manifest)
    return {**verify_predictor_resources(config), "reused": False}


def verify_predictor_resources(config: BuildConfig, *, full_digest: bool = True) -> dict[str, Any]:
    requested = set(config.predictors.enabled_plugins)
    if not requested:
        return {"status": "not_enabled", "enabled_plugins": []}
    manifest = _load_manifest(config.predictors.manifest)
    required_top = {
        "schema_version",
        "resource_type",
        "assembly",
        "vep_version",
        "installed_plugins",
        "plugins",
        "artifacts",
        "preparations",
        "structural_validation",
    }
    if set(manifest) != required_top:
        raise ResourceError("predictor resource manifest field set is invalid")
    if manifest["schema_version"] != MANIFEST_SCHEMA_VERSION:
        raise ResourceError("unsupported predictor resource manifest schema")
    if manifest["resource_type"] != "vep_predictor_resources":
        raise ResourceError("unexpected predictor resource manifest type")
    if manifest["assembly"] != "GRCh38" or manifest["vep_version"] != VEP_VERSION:
        raise ResourceError("predictor resources do not match GRCh38/VEP release boundary")
    installed_raw = manifest["installed_plugins"]
    if (
        not isinstance(installed_raw, list)
        or any(item not in PLUGIN_ORDER for item in installed_raw)
        or len(installed_raw) != len(set(installed_raw))
    ):
        raise ResourceError("predictor manifest installed_plugins is invalid")
    installed = set(installed_raw)
    expected_installed = [plugin for plugin in PLUGIN_ORDER if plugin in requested]
    if installed_raw != expected_installed:
        raise ResourceError(
            "predictor manifest plugin set/order differs from enabled configuration"
        )
    plugin_documents = manifest["plugins"]
    artifacts = manifest["artifacts"]
    if not isinstance(plugin_documents, dict) or not isinstance(artifacts, dict):
        raise ResourceError("predictor manifest plugins/artifacts must be objects")
    if set(plugin_documents) != installed:
        raise ResourceError("predictor manifest plugin documents do not match installed plugins")
    expected_artifacts = {key for plugin in installed for key in _expected_artifact_keys(plugin)}
    if set(artifacts) != expected_artifacts:
        raise ResourceError("predictor manifest contains a missing or unexpected artifact")

    verified_paths: dict[str, Path] = {}
    verified_plugins: dict[str, Any] = {}
    for plugin in sorted(requested):
        document = plugin_documents.get(plugin)
        if not isinstance(document, dict):
            raise ResourceError(f"predictor manifest lacks plugin document: {plugin}")
        plugin_config = _plugin_config(config, plugin)
        if document.get("license") != plugin_config.license:
            raise ResourceError(f"predictor licence mismatch for {plugin}")
        if document.get("resource_redistribution_allowed") is not False:
            raise ResourceError(
                f"predictor manifest may not assert resource redistribution for {plugin}"
            )
        artifact_keys = document.get("artifact_keys")
        if not isinstance(artifact_keys, list) or not artifact_keys:
            raise ResourceError(f"predictor manifest has no artifact keys for {plugin}")
        if artifact_keys != sorted(_expected_artifact_keys(plugin)):
            raise ResourceError(f"predictor artifact set mismatch for {plugin}")
        expected_document = _plugin_metadata(config, plugin, artifact_keys)
        for field in ("version", "license", "resource_redistribution_allowed", "coverage"):
            if document.get(field) != expected_document[field]:
                raise ResourceError(f"predictor manifest {field} mismatch for {plugin}")
        for key in artifact_keys:
            if not isinstance(key, str) or key not in artifacts:
                raise ResourceError(f"predictor manifest references unknown artifact {key!r}")
            value = artifacts[key]
            if not isinstance(value, dict) or value.get("plugin") != plugin:
                raise ResourceError(f"predictor artifact {key} belongs to wrong plugin")
            if value.get("kind") == "source_tree":
                verified_paths[key] = _verify_tree_artifact(
                    config.predictors.root, key, value, full_digest=full_digest
                )
            else:
                verified_paths[key] = _verify_manifest_artifact(
                    config.predictors.root, key, value, full_digest=full_digest
                )
        for key, spec in REMOTE_SPECS.items():
            if spec.plugin != plugin:
                continue
            expected_sha256 = _require_frozen(spec)
            value = artifacts.get(key)
            if (
                not isinstance(value, dict)
                or value.get("kind") != "upstream"
                or value.get("url") != spec.url
                or value.get("path") != spec.relative_path
                or value.get("bytes") != spec.bytes
                or value.get("sha256") != expected_sha256
                or value.get("upstream_md5") != spec.md5
            ):
                raise ResourceError(f"predictor upstream identity mismatch for {key}")
        source_relationships = {
            "loftee_source_tree": LOFTEE_SOURCE_SHA256,
            "loftee_sqlite": LOFTEE_SQL_GZ_SHA256,
            "revel_scores": REVEL_ARCHIVE_SHA256,
            "revel_scores_tbi": artifacts["revel_scores"]["sha256"] if plugin == "revel" else None,
            "alphamissense_scores_tbi": ALPHAMISSENSE_SHA256,
            "spliceai_indel_sentinel": (
                _require_frozen(REMOTE_SPECS["spliceai_snv"]) if plugin == "spliceai" else None
            ),
            "spliceai_indel_sentinel_tbi": artifacts["spliceai_indel_sentinel"]["sha256"]
            if plugin == "spliceai"
            else None,
        }
        for key in _expected_artifact_keys(plugin):
            expected_source = source_relationships.get(key)
            if (
                expected_source is not None
                and artifacts[key].get("source_sha256") != expected_source
            ):
                raise ResourceError(f"predictor derived-source identity mismatch for {key}")
        verified_plugins[plugin] = {
            "version": document.get("version"),
            "license": document["license"],
            "resource_redistribution_allowed": False,
            "coverage": document.get("coverage"),
            "artifact_sha256": {
                key: (
                    artifacts[key].get("tree_sha256")
                    if artifacts[key].get("kind") == "source_tree"
                    else artifacts[key].get("sha256")
                )
                for key in artifact_keys
            },
        }
    _validate_preparation_metadata(
        config.predictors.root,
        artifacts,
        manifest["preparations"],
        requested,
        full_digest=full_digest,
    )
    structures = _validate_plugin_structures(config.predictors.root, artifacts, requested)
    recorded_structures = manifest["structural_validation"]
    if recorded_structures != structures:
        raise ResourceError("predictor structural validation differs from manifest")
    return {
        "status": "verified",
        "assembly": "GRCh38",
        "vep_version": VEP_VERSION,
        "enabled_plugins": [plugin for plugin in PLUGIN_ORDER if plugin in requested],
        "manifest_sha256": sha256_file(config.predictors.manifest),
        "full_digest": full_digest,
        "plugins": verified_plugins,
        "structural_validation": structures,
    }


def resolve_predictor_paths(config: BuildConfig, *, full_digest: bool = True) -> dict[str, Path]:
    """Return validated local paths needed to construct VEP plugin arguments."""

    verify_predictor_resources(config, full_digest=full_digest)
    manifest = _load_manifest(config.predictors.manifest)
    artifacts = manifest["artifacts"]
    key_to_name = {
        "loftee_source_tree": "loftee_path",
        "loftee_gerp_bigwig": "loftee_gerp_bigwig",
        "loftee_human_ancestor": "loftee_human_ancestor",
        "loftee_sqlite": "loftee_conservation_sqlite",
        "cadd_snv": "cadd_snv",
        "cadd_indel": "cadd_indel",
        "revel_scores": "revel_scores",
        "alphamissense_scores": "alphamissense_scores",
        "spliceai_snv": "spliceai_snv",
        "spliceai_indel_sentinel": "spliceai_indel_sentinel",
    }
    requested = set(config.predictors.enabled_plugins)
    return {
        output_name: _safe_path(config.predictors.root, artifacts[key]["path"])
        for key, output_name in key_to_name.items()
        if key in artifacts and artifacts[key].get("plugin") in requested
    }
