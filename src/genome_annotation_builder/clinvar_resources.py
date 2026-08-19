"""Freeze, normalize, and verify an immutable ClinVar GRCh38 snapshot."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pysam

from .config import BuildConfig
from .constants import (
    CLINVAR_ARCHIVE,
    CLINVAR_ARCHIVE_BYTES,
    CLINVAR_ARCHIVE_MD5,
    CLINVAR_ARCHIVE_SHA256,
    CLINVAR_MD5_FILE_SHA256,
    CLINVAR_MD5_URL,
    CLINVAR_SNAPSHOT,
    CLINVAR_TBI_BYTES,
    CLINVAR_TBI_SHA256,
    CLINVAR_TBI_URL,
    CLINVAR_URL,
)
from .errors import ResourceError
from .util import atomic_output_directory, atomic_write_json, sha256_file

CORE_REGIONS = ",".join([*(str(value) for value in range(1, 23)), "X", "Y", "MT"])


def _resource_root(config: BuildConfig) -> Path:
    return config.clinvar_manifest.parent


def _archive_paths(config: BuildConfig) -> tuple[Path, Path, Path]:
    root = _resource_root(config) / "archive"
    vcf = root / CLINVAR_ARCHIVE
    return vcf, Path(f"{vcf}.tbi"), Path(f"{vcf}.md5")


def _md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_archive(config: BuildConfig) -> dict[str, Any]:
    vcf, index, md5_path = _archive_paths(config)
    expected = {
        vcf: (CLINVAR_ARCHIVE_BYTES, CLINVAR_ARCHIVE_SHA256),
        index: (CLINVAR_TBI_BYTES, CLINVAR_TBI_SHA256),
    }
    for path, (size, digest) in expected.items():
        if not path.is_file():
            raise ResourceError(f"ClinVar archive artifact is missing: {path}")
        if path.stat().st_size != size or sha256_file(path) != digest:
            raise ResourceError(f"ClinVar archive artifact failed size/SHA-256: {path}")
    if not md5_path.is_file() or sha256_file(md5_path) != CLINVAR_MD5_FILE_SHA256:
        raise ResourceError("ClinVar publisher MD5 file is missing or changed")
    try:
        published_md5 = md5_path.read_text(encoding="utf-8").split()[0].lower()
    except (OSError, IndexError) as exc:
        raise ResourceError("cannot parse ClinVar publisher MD5 file") from exc
    actual_md5 = _md5_file(vcf)
    if published_md5 != CLINVAR_ARCHIVE_MD5 or actual_md5 != CLINVAR_ARCHIVE_MD5:
        raise ResourceError("ClinVar publisher/local MD5 mismatch")
    handle = pysam.VariantFile(str(vcf))
    try:
        if list(handle.header.samples):
            raise ResourceError("ClinVar resource unexpectedly has sample columns")
        if "fileDate=2026-08-08" not in str(handle.header):
            raise ResourceError("ClinVar VCF fileDate does not match the pinned snapshot")
        if "reference=GRCh38" not in str(handle.header):
            raise ResourceError("ClinVar VCF is not GRCh38")
    finally:
        handle.close()
    return {
        "archive_sha256": CLINVAR_ARCHIVE_SHA256,
        "archive_md5": actual_md5,
        "archive_bytes": vcf.stat().st_size,
        "source_index_sha256": CLINVAR_TBI_SHA256,
        "publisher_md5_file_sha256": CLINVAR_MD5_FILE_SHA256,
    }


def _download_one(config: BuildConfig, destination: Path, url: str) -> None:
    if not config.allow_download:
        raise ResourceError(
            f"ClinVar artifact is absent and resources.allow_download is false: {destination}"
        )
    curl = shutil.which("curl")
    if curl is None:
        raise ResourceError("curl is required for restart-safe ClinVar download")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f".{destination.name}.partial")
    completed = subprocess.run(
        [
            curl,
            "--fail",
            "--location",
            "--retry",
            "8",
            "--continue-at",
            "-",
            "--output",
            str(partial),
            url,
        ],
        check=False,
    )
    if completed.returncode != 0:
        raise ResourceError(
            f"ClinVar download failed with exit {completed.returncode}; retained {partial}"
        )
    os.replace(partial, destination)


def _record_count(config: BuildConfig, path: Path) -> int:
    completed = subprocess.run(
        [config.bcftools, "index", "--nrecords", str(path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise ResourceError(f"cannot count ClinVar records: {completed.stderr}")
    return int(completed.stdout.strip())


def _normalize(config: BuildConfig, archive: Path) -> dict[str, Any]:
    with atomic_output_directory(config.clinvar_vcf.parent) as staging:
        output = staging / config.clinvar_vcf.name
        view_log = staging / "view.stderr.log"
        norm_log = staging / "norm.stderr.log"
        sort_log = staging / "sort.stderr.log"
        with (
            view_log.open("wb") as view_stderr,
            norm_log.open("wb") as norm_stderr,
            sort_log.open("wb") as sort_stderr,
        ):
            view = subprocess.Popen(
                [
                    config.bcftools,
                    "view",
                    "--no-version",
                    "--regions",
                    CORE_REGIONS,
                    "--output-type",
                    "u",
                    str(archive),
                ],
                stdout=subprocess.PIPE,
                stderr=view_stderr,
            )
            assert view.stdout is not None
            norm = subprocess.Popen(
                [
                    config.bcftools,
                    "norm",
                    "--no-version",
                    "--fasta-ref",
                    str(config.reference_fasta),
                    "--check-ref",
                    "e",
                    "--output-type",
                    "u",
                    "-",
                ],
                stdin=view.stdout,
                stdout=subprocess.PIPE,
                stderr=norm_stderr,
            )
            view.stdout.close()
            assert norm.stdout is not None
            sort = subprocess.run(
                [
                    config.bcftools,
                    "sort",
                    "--max-mem",
                    "2G",
                    "--temp-dir",
                    str(staging / "sort-tmp"),
                    "--output-type",
                    "z",
                    "--output-file",
                    str(output),
                    "-",
                ],
                stdin=norm.stdout,
                stderr=sort_stderr,
                check=False,
            )
            norm.stdout.close()
            norm_return = norm.wait()
            view_return = view.wait()
        if (view_return, norm_return, sort.returncode) != (0, 0, 0):
            details = "\n".join(
                path.read_text(encoding="utf-8", errors="replace")[-4000:]
                for path in (view_log, norm_log, sort_log)
            )
            raise ResourceError(
                "ClinVar core filtering/normalization failed: "
                f"view={view_return}, norm={norm_return}, sort={sort.returncode}\n{details}"
            )
        index = subprocess.run(
            [config.bcftools, "index", "--tbi", "--force", str(output)],
            check=False,
            capture_output=True,
            text=True,
        )
        if index.returncode != 0:
            raise ResourceError(f"could not index normalized ClinVar: {index.stderr}")
        record_count = _record_count(config, output)
        output_sha256 = sha256_file(output)
        index_sha256 = sha256_file(Path(f"{output}.tbi"))
    return {
        "normalized_records": record_count,
        "normalized_vcf_sha256": output_sha256,
        "normalized_index_sha256": index_sha256,
    }


def verify_clinvar(config: BuildConfig) -> dict[str, Any]:
    try:
        manifest = json.loads(config.clinvar_manifest.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ResourceError(f"ClinVar manifest not found: {config.clinvar_manifest}") from exc
    except json.JSONDecodeError as exc:
        raise ResourceError(f"invalid ClinVar manifest: {exc}") from exc
    required = {
        "resource_type",
        "snapshot",
        "assembly",
        "archive_sha256",
        "normalized_vcf_sha256",
        "normalized_index_sha256",
        "normalized_records",
        "reference_fasta_sha256",
    }
    if not isinstance(manifest, dict) or not required.issubset(manifest):
        raise ResourceError("ClinVar manifest is missing required fields")
    if manifest["resource_type"] != "ncbi_clinvar_vcf":
        raise ResourceError("unexpected ClinVar resource type")
    if str(manifest["snapshot"]) != config.clinvar_snapshot or manifest["assembly"] != "GRCh38":
        raise ResourceError("ClinVar snapshot/assembly differs from configuration")
    if manifest["archive_sha256"] != CLINVAR_ARCHIVE_SHA256:
        raise ResourceError("ClinVar manifest archive SHA-256 is not pinned")
    if not config.clinvar_vcf.is_file() or not Path(f"{config.clinvar_vcf}.tbi").is_file():
        raise ResourceError("normalized ClinVar VCF or index is missing")
    if sha256_file(config.clinvar_vcf) != manifest["normalized_vcf_sha256"]:
        raise ResourceError("normalized ClinVar VCF SHA-256 mismatch")
    if sha256_file(Path(f"{config.clinvar_vcf}.tbi")) != manifest["normalized_index_sha256"]:
        raise ResourceError("normalized ClinVar index SHA-256 mismatch")
    reference_sha = sha256_file(config.reference_fasta)
    if reference_sha != manifest["reference_fasta_sha256"]:
        raise ResourceError("ClinVar normalization reference SHA-256 mismatch")
    if _record_count(config, config.clinvar_vcf) != manifest["normalized_records"]:
        raise ResourceError("normalized ClinVar record count mismatch")
    handle = pysam.VariantFile(str(config.clinvar_vcf))
    try:
        if list(handle.header.samples) or "fileDate=2026-08-08" not in str(handle.header):
            raise ResourceError("normalized ClinVar header is inconsistent")
    finally:
        handle.close()
    return manifest


def sync_clinvar(config: BuildConfig) -> dict[str, Any]:
    if config.clinvar_snapshot != CLINVAR_SNAPSHOT:
        raise ResourceError(
            f"pipeline 0.1 supports ClinVar {CLINVAR_SNAPSHOT}, got {config.clinvar_snapshot}"
        )
    if config.clinvar_manifest.exists() or config.clinvar_vcf.exists():
        return verify_clinvar(config)
    archive, index, md5_path = _archive_paths(config)
    for destination, url in (
        (archive, CLINVAR_URL),
        (index, CLINVAR_TBI_URL),
        (md5_path, CLINVAR_MD5_URL),
    ):
        if not destination.is_file():
            _download_one(config, destination, url)
    archive_identity = _verify_archive(config)
    normalized = _normalize(config, archive)
    manifest = {
        "resource_type": "ncbi_clinvar_vcf",
        "snapshot": CLINVAR_SNAPSHOT,
        "assembly": "GRCh38",
        "source_url": CLINVAR_URL,
        "source_index_url": CLINVAR_TBI_URL,
        "publisher_checksum_url": CLINVAR_MD5_URL,
        "license_or_terms": "NCBI ClinVar data; see https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/",
        "installed_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "reference_fasta_sha256": sha256_file(config.reference_fasta),
        "normalization": (
            "bcftools view primary contigs; bcftools norm --check-ref e -f GRCh38; "
            "bcftools sort; no-version"
        ),
        **archive_identity,
        **normalized,
    }
    atomic_write_json(config.clinvar_manifest, manifest)
    return verify_clinvar(config)
