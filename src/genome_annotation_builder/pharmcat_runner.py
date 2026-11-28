"""Offline, digest-pinned PharmCAT container execution."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from .config import BuildConfig
from .constants import PHARMCAT_VERSION, STRUCTURALLY_COMPLEX_PGX_GENES
from .errors import InputValidationError, ResourceError
from .util import atomic_output_directory, atomic_write_json, sha256_file

DANGEROUS_PHARMCAT_ARGUMENTS = frozenset(
    {
        "-0",
        "--missing-to-ref",
        "--absent-to-ref",
        "--unspecified-to-ref",
        "-research",
        "--research-mode",
    }
)


def _docker_base(config: BuildConfig) -> list[str]:
    executable = shutil.which(config.pgx.docker_executable)
    if executable is None:
        raise ResourceError(f"PharmCAT container runtime not found: {config.pgx.docker_executable}")
    return [executable]


def _hardened_run_prefix(config: BuildConfig) -> list[str]:
    return [
        *_docker_base(config),
        "run",
        "--rm",
        "--pull",
        "never",
        "--network",
        "none",
        "--read-only",
        "--cap-drop",
        "ALL",
        "--security-opt",
        "no-new-privileges",
        "--pids-limit",
        "512",
        "--tmpfs",
        "/tmp:rw,nosuid,noexec,size=2g",
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "--env",
        "HOME=/tmp",
        "--env",
        "PYTHONDONTWRITEBYTECODE=1",
    ]


def verify_pharmcat_runtime(config: BuildConfig) -> dict[str, Any]:
    """Prove the exact local image and version without enabling networking."""

    inspect = subprocess.run(
        [
            *_docker_base(config),
            "image",
            "inspect",
            config.pgx.container_image,
            "--format",
            "{{json .RepoDigests}}",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if inspect.returncode != 0:
        raise ResourceError(
            "pinned PharmCAT image is not available locally; pull exactly "
            f"{config.pgx.container_image}: {(inspect.stderr or inspect.stdout).strip()}"
        )
    try:
        repo_digests = json.loads(inspect.stdout.strip())
    except json.JSONDecodeError as exc:
        raise ResourceError("could not parse PharmCAT image RepoDigests") from exc
    if config.pgx.container_image not in repo_digests:
        raise ResourceError("local PharmCAT image does not expose the configured immutable digest")

    command = [
        *_hardened_run_prefix(config),
        config.pgx.container_image,
        "pharmcat_pipeline",
        "-V",
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        raise ResourceError(f"PharmCAT version command failed: {output[-2000:]}")
    match = re.search(r"PharmCAT\s+v([0-9]+(?:\.[0-9]+){2})", output)
    if match is None or match.group(1) != PHARMCAT_VERSION:
        raise ResourceError(
            f"PharmCAT version mismatch: expected {PHARMCAT_VERSION}, output={output!r}"
        )
    return {
        "status": "verified",
        "version": match.group(1),
        "container_image": config.pgx.container_image,
        "network": "none",
        "read_only_root": True,
        "runs_as_host_user": True,
    }


def sync_pharmcat_runtime(config: BuildConfig) -> dict[str, Any]:
    """Pull only the configured immutable image, then verify its identity."""

    try:
        return verify_pharmcat_runtime(config)
    except ResourceError:
        if not config.allow_download:
            raise
    completed = subprocess.run(
        [*_docker_base(config), "pull", config.pgx.container_image],
        check=False,
        capture_output=True,
        text=True,
        timeout=config.pgx.timeout_seconds,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise ResourceError(f"pinned PharmCAT image pull failed: {detail[-4000:]}")
    return verify_pharmcat_runtime(config)


def _safe_genes(resource_dir: Path) -> list[str]:
    table = pq.read_table(resource_dir / "pgx_loci.parquet", columns=["genes"])
    genes = {
        str(gene)
        for values in table.column("genes").to_pylist()
        for gene in (values or [])
        if str(gene) not in STRUCTURALLY_COMPLEX_PGX_GENES
    }
    if not genes:
        raise InputValidationError("PharmCAT resource has no non-complex matcher genes")
    return sorted(genes)


def _pipeline_arguments(genes: list[str]) -> list[str]:
    arguments = [
        "pharmcat_pipeline",
        "/input/pgx-input.vcf.gz",
        "--output-dir",
        "/output",
        "--base-filename",
        "gab-pgx",
        "--genes",
        ",".join(genes),
        "--reporter-save-json",
        "--reporter-save-calls-only-tsv",
        "--keep-intermediate-preprocessor-files",
        "--max-concurrent-processes",
        "1",
        "--max-memory",
        "8G",
    ]
    forbidden = set(arguments) & DANGEROUS_PHARMCAT_ARGUMENTS
    if forbidden:
        raise InputValidationError(
            "internal error: unsafe PharmCAT argument requested: " + ", ".join(sorted(forbidden))
        )
    return arguments


def _container_command(
    config: BuildConfig, input_dir: Path, output_dir: Path, genes: list[str]
) -> list[str]:
    return [
        *_hardened_run_prefix(config),
        "--volume",
        f"{input_dir.resolve()}:/input:ro",
        "--volume",
        f"{output_dir.resolve()}:/output:rw",
        "--workdir",
        "/output",
        config.pgx.container_image,
        *_pipeline_arguments(genes),
    ]


def _tail(path: Path, maximum: int = 6000) -> str:
    with path.open("rb") as handle:
        size = path.stat().st_size
        if size > maximum:
            handle.seek(-maximum, os.SEEK_END)
        return handle.read().decode("utf-8", errors="replace")


def _exactly_one(paths: list[Path], label: str) -> Path:
    if len(paths) != 1:
        raise InputValidationError(
            f"PharmCAT produced {len(paths)} {label} files; exactly one was expected"
        )
    return paths[0]


def _canonical_report(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise InputValidationError("PharmCAT report is not valid UTF-8 JSON") from exc
    if not isinstance(raw, dict):
        raise InputValidationError("PharmCAT report root must be an object")
    required = {
        "title",
        "timestamp",
        "pharmcatVersion",
        "dataVersion",
        "genes",
        "drugs",
        "messages",
    }
    if not required.issubset(raw):
        raise InputValidationError("PharmCAT report is missing required top-level fields")
    if str(raw["pharmcatVersion"]) != PHARMCAT_VERSION:
        raise InputValidationError("PharmCAT report version differs from the pinned runtime")
    if not isinstance(raw["genes"], dict) or not isinstance(raw["drugs"], dict):
        raise InputValidationError("PharmCAT report genes/drugs fields have invalid types")
    canonical = dict(raw)
    canonical.pop("timestamp", None)
    return raw, canonical


def run_pharmcat(
    config: BuildConfig,
    *,
    materialization_dir: str | Path,
    resource_dir: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    """Run the pinned container with no network and a read-only root filesystem."""

    if not config.pgx.enabled:
        raise InputValidationError("PharmCAT was invoked while PGx is disabled")
    materialization_dir = Path(materialization_dir)
    resource_dir = Path(resource_dir)
    destination = Path(destination)
    input_vcf = materialization_dir / "pgx-input.vcf.gz"
    input_index = Path(f"{input_vcf}.tbi")
    if not input_vcf.is_file() or not input_index.is_file():
        raise InputValidationError("materialized PharmCAT input VCF or index is missing")
    runtime = verify_pharmcat_runtime(config)
    genes = _safe_genes(resource_dir)

    with atomic_output_directory(destination) as staging:
        stdout_log = staging / "pharmcat.stdout.log"
        stderr_log = staging / "pharmcat.stderr.log"
        command = _container_command(config, materialization_dir, staging, genes)
        with stdout_log.open("wb") as stdout, stderr_log.open("wb") as stderr:
            try:
                completed = subprocess.run(
                    command,
                    check=False,
                    stdout=stdout,
                    stderr=stderr,
                    timeout=config.pgx.timeout_seconds,
                )
            except subprocess.TimeoutExpired as exc:
                raise InputValidationError(
                    f"PharmCAT exceeded the configured {config.pgx.timeout_seconds}s timeout"
                ) from exc
        if completed.returncode != 0:
            detail = _tail(stderr_log) or _tail(stdout_log)
            raise InputValidationError(
                f"PharmCAT failed with exit {completed.returncode}: {detail[-4000:]}"
            )

        raw_report = _exactly_one(sorted(staging.glob("gab-pgx*.report.json")), "JSON report")
        raw_calls = _exactly_one(
            sorted(staging.glob("gab-pgx*.report.tsv")), "calls-only TSV report"
        )
        raw_document, canonical_document = _canonical_report(raw_report)
        stable_raw = staging / "pharmcat.report.raw.json"
        stable_calls = staging / "pharmcat.calls.tsv"
        shutil.copyfile(raw_report, stable_raw)
        shutil.copyfile(raw_calls, stable_calls)
        canonical_path = staging / "pharmcat.report.canonical.json"
        atomic_write_json(canonical_path, canonical_document)

        output_files = sorted(path for path in staging.rglob("*") if path.is_file())
        outputs = {str(path.relative_to(staging)): sha256_file(path) for path in output_files}
        summary = {
            **config.checkpoint_identity,
            "status": "PASS",
            "version": PHARMCAT_VERSION,
            "data_version": str(raw_document["dataVersion"]),
            "container_image": config.pgx.container_image,
            "runtime": runtime,
            "input_vcf_sha256": sha256_file(input_vcf),
            "safe_matcher_genes": genes,
            "excluded_complex_genes": sorted(STRUCTURALLY_COMPLEX_PGX_GENES),
            "dangerous_missing_to_reference_options_used": False,
            "research_mode_used": False,
            "network": "none",
            "read_only_root": True,
            "raw_report_sha256": sha256_file(stable_raw),
            "canonical_report_sha256": sha256_file(canonical_path),
            "calls_tsv_sha256": sha256_file(stable_calls),
            "outputs": outputs,
            "command_semantics": _pipeline_arguments(genes),
            "stdout_tail": _tail(stdout_log),
            "stderr_tail": _tail(stderr_log),
        }
        atomic_write_json(staging / "pharmcat_run_summary.json", summary)
    return summary
