"""Read-only environment, input, reference, and resource preflight."""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pysam

from .acmg_resources import verify_acmg_sf
from .clingen_resources import verify_clingen_resources
from .clinvar_resources import verify_clinvar
from .config import BuildConfig
from .errors import InputValidationError
from .family_evidence import verify_family_resources
from .gwas_resources import verify_gwas_catalog
from .pgs import verify_pgs_resources
from .pgx_resources import verify_pgx_resources
from .pharmcat_runner import verify_pharmcat_runtime
from .predictor_doctor import verify_predictor_doctor
from .predictor_resources import verify_predictor_resources
from .resources import verify_pinned_resources
from .util import atomic_write_json, sha256_file, source_tree_digest
from .vep_resources import verify_vep_cache
from .vep_runner import verify_predictor_plugin_runtime, verify_vep_runtime


def _redacted_file(path: Path, include_path: bool) -> dict[str, Any]:
    stat = path.stat()
    value: dict[str, Any] = {
        "basename": path.name,
        "bytes": stat.st_size,
        "sha256": sha256_file(path),
    }
    if include_path:
        value["path"] = str(path)
    return value


def _sample_label(sample: str, include_identifier: bool) -> str:
    if include_identifier:
        return sample
    return "sha256:" + hashlib.sha256(sample.encode("utf-8")).hexdigest()


def _command_version(executable: str, args: list[str]) -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise InputValidationError(f"required executable not found: {executable}")
    completed = subprocess.run(
        [resolved, *args],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise InputValidationError(f"failed to query {executable} version: {detail}")
    output = (completed.stdout or completed.stderr).strip().splitlines()
    return output[0] if output else resolved


def _index_candidates(vcf: Path) -> list[Path]:
    return [Path(f"{vcf}.tbi"), Path(f"{vcf}.csi")]


def run_preflight(config: BuildConfig, write_report: bool = True) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []

    if not config.input_vcf.is_file():
        errors.append(f"input VCF not found: {config.input_vcf}")
    if not config.reference_fasta.is_file():
        errors.append(f"reference FASTA not found: {config.reference_fasta}")
    fai = Path(f"{config.reference_fasta}.fai")
    if not fai.is_file():
        errors.append(f"reference FASTA index not found: {fai}")
    if config.input_vcf.is_file() and not any(
        path.is_file() for path in _index_candidates(config.input_vcf)
    ):
        errors.append("input VCF has no adjacent .tbi or .csi index")
    if config.target_catalog_enabled and not config.target_catalog_path.is_file():
        errors.append(f"target catalogue not found: {config.target_catalog_path}")
    if errors:
        raise InputValidationError("preflight failed: " + "; ".join(errors))

    assert config.input_vcf.is_file()
    assert config.reference_fasta.is_file()

    try:
        variant_file = pysam.VariantFile(str(config.input_vcf))
    except (OSError, ValueError) as exc:
        raise InputValidationError(f"cannot open input VCF: {exc}") from exc
    try:
        samples = list(variant_file.header.samples)
        if config.sample not in samples:
            raise InputValidationError(
                f"configured sample is absent from VCF; available sample count: {len(samples)}"
            )
        if len(samples) > 1:
            warnings.append("multisample input: only the explicitly configured sample will be read")
        vcf_contigs = {name: record.length for name, record in variant_file.header.contigs.items()}
    finally:
        variant_file.close()

    try:
        fasta = pysam.FastaFile(str(config.reference_fasta))
    except (OSError, ValueError) as exc:
        raise InputValidationError(f"cannot open reference FASTA: {exc}") from exc
    try:
        reference_contigs = dict(zip(fasta.references, fasta.lengths, strict=True))
    finally:
        fasta.close()

    mismatches: list[str] = []
    shared = 0
    for name, declared_length in vcf_contigs.items():
        if name not in reference_contigs:
            mismatches.append(f"{name}: absent from FASTA")
            continue
        shared += 1
        if declared_length is not None and int(declared_length) != reference_contigs[name]:
            mismatches.append(
                f"{name}: VCF length {declared_length}, FASTA length {reference_contigs[name]}"
            )
    if mismatches:
        preview = "; ".join(mismatches[:10])
        raise InputValidationError(f"VCF/FASTA contig mismatch ({len(mismatches)}): {preview}")
    if shared == 0:
        raise InputValidationError("VCF and FASTA have no shared contigs")

    reference_sha = sha256_file(config.reference_fasta)
    if config.reference_expected_sha256 and reference_sha != config.reference_expected_sha256:
        raise InputValidationError(
            "reference FASTA SHA-256 mismatch: "
            f"expected {config.reference_expected_sha256}, got {reference_sha}"
        )

    resources: dict[str, Any] = {"genome_spec": verify_pinned_resources(config)}
    if config.vep_enabled:
        cache_manifest = verify_vep_cache(config)
        vep_version, _version_command = verify_vep_runtime(config)
        resources["vep"] = {
            "status": "verified",
            "vep_version": vep_version,
            "cache_version": cache_manifest["cache_version"],
            "cache_tree_sha256": cache_manifest["tree_sha256"],
            "archive_sha256": cache_manifest["archive_sha256"],
            "container_image": (
                config.vep_container_image if config.vep_mode == "docker" else None
            ),
        }
    else:
        resources["vep"] = {"status": "not_enabled"}
    predictor_resources = verify_predictor_resources(config)
    if config.predictors.enabled:
        predictor_resources["runtime"] = verify_predictor_plugin_runtime(config)
        predictor_resources["semantic_controls"] = verify_predictor_doctor(config)
    resources["predictors"] = predictor_resources
    if config.clinvar_enabled:
        clinvar_manifest = verify_clinvar(config)
        resources["clinvar"] = {
            "status": "verified",
            "snapshot": clinvar_manifest["snapshot"],
            "archive_sha256": clinvar_manifest["archive_sha256"],
            "normalized_vcf_sha256": clinvar_manifest["normalized_vcf_sha256"],
            "normalized_records": clinvar_manifest["normalized_records"],
        }
    else:
        resources["clinvar"] = {"status": "not_enabled"}
    if config.acmg_sf_enabled:
        acmg_manifest = verify_acmg_sf(config)
        resources["acmg_sf"] = {
            "status": "verified",
            "version": acmg_manifest["version"],
            "snapshot": acmg_manifest["snapshot"],
            "canonical_sha256": acmg_manifest["canonical_sha256"],
            "source_workbook_sha256": acmg_manifest["source_workbook_sha256"],
            "policy_pmid": acmg_manifest["policy_pmid"],
            "policy_doi": acmg_manifest["policy_doi"],
            "counts": acmg_manifest["counts"],
        }
    else:
        resources["acmg_sf"] = {"status": "not_enabled"}
    resources["clingen"] = verify_clingen_resources(config)
    resources["pgs"] = verify_pgs_resources(config)
    pgx_resources = verify_pgx_resources(config)
    if config.pgx.enabled:
        pgx_resources["runtime"] = verify_pharmcat_runtime(config)
    resources["pgx"] = pgx_resources
    resources["family_evidence"] = verify_family_resources(config)
    resources["gwas_catalog"] = verify_gwas_catalog(config)
    if config.target_catalog_enabled:
        target_sha256 = sha256_file(config.target_catalog_path)
        if target_sha256 != config.target_catalog_expected_sha256:
            raise InputValidationError(
                "target catalogue SHA-256 mismatch: expected "
                f"{config.target_catalog_expected_sha256}, got {target_sha256}"
            )
        resources["target_catalog"] = {
            "status": "verified",
            **_redacted_file(config.target_catalog_path, config.include_input_paths),
            "assay_scope": config.target_catalog_assay_scope,
        }
    else:
        resources["target_catalog"] = {"status": "not_enabled"}
    bcftools_version = _command_version(config.bcftools, ["--version"])
    disk = shutil.disk_usage(
        config.work_dir.parent if config.work_dir.parent.exists() else config.config_path.parent
    )
    if disk.free < 50 * 1024**3:
        warnings.append("less than 50 GiB free on the work filesystem")

    report: dict[str, Any] = {
        "status": "PASS",
        "checked_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "input": _redacted_file(config.input_vcf, config.include_input_paths),
        "sample": _sample_label(config.sample, config.include_sample_identifier),
        "sample_count_in_vcf": len(samples),
        "reference": _redacted_file(config.reference_fasta, config.include_input_paths),
        "reference_fai_sha256": sha256_file(fai),
        "vcf_contigs": len(vcf_contigs),
        "reference_contigs": len(reference_contigs),
        "shared_contigs": shared,
        "resources": resources,
        "software": {
            "builder_implementation_sha256": source_tree_digest(Path(__file__).resolve().parent),
            "pysam": pysam.__version__,
            "htslib": getattr(pysam, "__samtools_version__", None),
            "bcftools": bcftools_version,
        },
        "filesystem": {
            "free_bytes": disk.free,
            "total_bytes": disk.total,
        },
        "warnings": warnings,
    }
    if write_report:
        config.work_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_json(config.work_dir / "preflight.json", report)
    return report
