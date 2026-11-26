"""Restart-safe orchestration of the Release-A build."""

from __future__ import annotations

import json
import shutil
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .bundle import assemble_bundle
from .config import BuildConfig
from .constants import PIPELINE_VERSION
from .errors import GenomeBuilderError
from .family_evidence import build_family_evidence
from .normalization import normalize_projected_vcf
from .pgs import prepare_pgs_resources, score_pgs
from .pgx_materialize import materialize_pgx_input
from .pgx_project import project_pgx_results
from .pgx_resources import prepare_pgx_resources
from .pharmcat_runner import run_pharmcat
from .preflight import run_preflight
from .target_callability import build_target_callability
from .util import atomic_write_json, directory_digest, sha256_file, sha256_json
from .validator import validate_bundle
from .variant_build import build_variant_datasets
from .vcf_ingest import ingest_reference_confidence_vcf
from .vep_runner import run_vep_annotation


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GenomeBuilderError(f"invalid stage checkpoint {path}: {exc}") from exc
    if not isinstance(value, dict) or value.get("status") != "PASS":
        raise GenomeBuilderError(f"stage checkpoint is not PASS: {path}")
    return value


def _verify_outputs(stage: Path, summary: dict[str, Any]) -> None:
    for relative, expected in summary.get("outputs", {}).items():
        target = stage / relative
        if not target.exists():
            raise GenomeBuilderError(f"checkpoint output is missing: {target}")
        actual = directory_digest(target) if target.is_dir() else sha256_file(target)
        if actual != expected:
            raise GenomeBuilderError(f"checkpoint output digest mismatch: {target}")


def _reuse_stage(stage: Path, summary_name: str) -> dict[str, Any] | None:
    if not stage.exists():
        return None
    summary = _load_json(stage / summary_name)
    _verify_outputs(stage, summary)
    return summary


def _runtime(path: Path, run_key: str, state: str, percent: int, **extra: Any) -> None:
    atomic_write_json(
        path,
        {
            "run_key": run_key,
            "state": state,
            "approximate_percent": percent,
            "updated_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
            **extra,
        },
    )


def build_bundle(config: BuildConfig) -> dict[str, Any]:
    preflight = run_preflight(config)
    run_key = sha256_json(
        {
            "input_sha256": preflight["input"]["sha256"],
            "reference_sha256": preflight["reference"]["sha256"],
            "schema_commit": config.genome_spec_commit,
            "pipeline_version": PIPELINE_VERSION,
            "builder_implementation_sha256": preflight["software"]["builder_implementation_sha256"],
            "sample": config.sample,
            "output_profile": {
                "agent_exploration": config.agent_exploration_profile,
                "reference_confidence": config.reference_confidence_extension,
                "partition_by_chromosome": config.partition_variants_by_chromosome,
                "include_input_paths": config.include_input_paths,
                "include_sample_identifier": config.include_sample_identifier,
            },
            "normalization": config.check_ref,
            "quality": vars(config.quality),
            "vep": {
                "enabled": config.vep_enabled,
                "version": config.vep_expected_version if config.vep_enabled else None,
                "cache_version": config.vep_cache_version if config.vep_enabled else None,
                "cache_tree_sha256": (
                    preflight["resources"]["vep"].get("cache_tree_sha256")
                    if config.vep_enabled
                    else None
                ),
                "runtime": (
                    config.vep_container_image
                    if config.vep_enabled and config.vep_mode == "docker"
                    else config.vep_executable
                    if config.vep_enabled
                    else None
                ),
            },
            "predictors": preflight["resources"]["predictors"],
            "clinvar": {
                "enabled": config.clinvar_enabled,
                "snapshot": config.clinvar_snapshot if config.clinvar_enabled else None,
                "normalized_vcf_sha256": (
                    preflight["resources"]["clinvar"].get("normalized_vcf_sha256")
                    if config.clinvar_enabled
                    else None
                ),
            },
            "clingen": preflight["resources"]["clingen"],
            "acmg_sf": {
                "enabled": config.acmg_sf_enabled,
                "version": config.acmg_sf_version if config.acmg_sf_enabled else None,
                "snapshot": config.acmg_sf_snapshot if config.acmg_sf_enabled else None,
                "canonical_sha256": (
                    preflight["resources"]["acmg_sf"].get("canonical_sha256")
                    if config.acmg_sf_enabled
                    else None
                ),
                "min_clinvar_review_stars": (
                    config.acmg_sf_min_clinvar_review_stars if config.acmg_sf_enabled else None
                ),
            },
            "pgs": preflight["resources"]["pgs"],
            "pgx": preflight["resources"]["pgx"],
            "family_evidence": preflight["resources"]["family_evidence"],
            "gwas_catalog": preflight["resources"]["gwas_catalog"],
            "target_catalog": {
                "enabled": config.target_catalog_enabled,
                "sha256": (
                    preflight["resources"]["target_catalog"].get("sha256")
                    if config.target_catalog_enabled
                    else None
                ),
                "assay_scope": (
                    config.target_catalog_assay_scope if config.target_catalog_enabled else None
                ),
            },
        }
    )[:20]
    run_root = config.work_dir / f"run-{run_key}"
    run_root.mkdir(parents=True, exist_ok=True)
    runtime_path = run_root / "runtime.json"

    if config.bundle_path.exists() and not config.overwrite:
        manifest_path = config.bundle_path / "manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise GenomeBuilderError(f"cannot read existing bundle manifest: {exc}") from exc
        provenance = manifest.get("resource_provenance") if isinstance(manifest, dict) else None
        if not isinstance(provenance, dict) or provenance.get("build_run_key") != run_key:
            raise GenomeBuilderError(
                "existing bundle belongs to a different or unknown build; "
                "use a new output path or explicitly enable overwrite"
            )
        validation = validate_bundle(
            config.bundle_path,
            reference_fasta=config.reference_fasta,
            expected_schema_path=config.resource_cache
            / "genome-spec"
            / config.genome_spec_commit
            / "schema.json",
            run_profile_queries=True,
        )
        return {
            "status": "PASS",
            "reused": True,
            "bundle": str(config.bundle_path),
            "run_key": run_key,
            "validation": validation,
        }

    disk = shutil.disk_usage(run_root)
    if disk.free < 50 * 1024**3:
        raise GenomeBuilderError("less than 50 GiB free; refusing to start a full build")

    _runtime(runtime_path, run_key, "preflight_complete", 5)
    ingest_stage = run_root / "01-ingest"
    ingest = _reuse_stage(ingest_stage, "ingest_summary.json")
    if ingest is None:
        _runtime(runtime_path, run_key, "ingesting_reference_confidence_vcf", 10)
        ingest = ingest_reference_confidence_vcf(config, preflight["input"]["sha256"], ingest_stage)
    _runtime(runtime_path, run_key, "ingest_complete", 45)

    normalization_stage = run_root / "02-normalization"
    normalization = _reuse_stage(normalization_stage, "normalization_summary.json")
    if normalization is None:
        _runtime(runtime_path, run_key, "normalizing_and_atomizing", 50)
        normalization = normalize_projected_vcf(
            config, ingest_stage / "projected.vcf.gz", normalization_stage
        )
    _runtime(runtime_path, run_key, "normalization_complete", 65)

    vep_summary: dict[str, Any] | None = None
    variant_input = normalization_stage / "normalized.vcf.gz"
    if config.vep_enabled:
        vep_stage = run_root / "03-vep"
        vep_summary = _reuse_stage(vep_stage, "vep_summary.json")
        if vep_summary is None:
            _runtime(runtime_path, run_key, "annotating_with_vep", 68)
            vep_summary = run_vep_annotation(config, variant_input, vep_stage)
        if vep_summary["record_count"] != normalization["record_count"]:
            raise GenomeBuilderError(
                "VEP output record count differs from normalized input: "
                f"{vep_summary['record_count']} != {normalization['record_count']}"
            )
        variant_input = vep_stage / "annotated.vcf.gz"
        _runtime(runtime_path, run_key, "vep_annotation_complete", 78)

    variant_stage = run_root / "04-variant-datasets"
    variants = _reuse_stage(variant_stage, "variant_build_summary.json")
    if variants is None:
        _runtime(runtime_path, run_key, "building_typed_parquet", 80)
        variants = build_variant_datasets(
            config,
            variant_input,
            preflight["input"]["sha256"],
            ingest["accepted_evidence_multiset"],
            variant_stage,
        )
    _runtime(runtime_path, run_key, "variant_datasets_complete", 90)

    family_evidence: dict[str, Any] | None = None
    family_evidence_stage: Path | None = None
    if config.family_evidence.enabled:
        family_evidence_stage = run_root / "05-family-evidence"
        family_evidence = _reuse_stage(family_evidence_stage, "family_evidence_summary.json")
        if family_evidence is None:
            _runtime(runtime_path, run_key, "materializing_family_evidence", 91)
            family_evidence = build_family_evidence(
                config,
                variant_dir=variant_stage,
                destination=family_evidence_stage,
            )
        _runtime(runtime_path, run_key, "family_evidence_complete", 92)

    pgs_resources: dict[str, Any] | None = None
    pgs_resource_stage: Path | None = None
    if config.pgs_enabled:
        pgs_resource_stage = run_root / "05-pgs-resources"
        pgs_resources = _reuse_stage(pgs_resource_stage, "pgs_resource_summary.json")
        if pgs_resources is None:
            _runtime(runtime_path, run_key, "preparing_pgs_resources", 91)
            pgs_resources = prepare_pgs_resources(config, pgs_resource_stage)
        _runtime(runtime_path, run_key, "pgs_resources_complete", 92)

    pgx_resources: dict[str, Any] | None = None
    pgx_resource_stage: Path | None = None
    if config.pgx.enabled:
        pgx_resource_stage = run_root / "06-pgx-resources"
        pgx_resources = _reuse_stage(pgx_resource_stage, "pgx_resource_summary.json")
        if pgx_resources is None:
            _runtime(runtime_path, run_key, "preparing_pgx_resources", 92)
            pgx_resources = prepare_pgx_resources(
                config,
                pgx_resource_stage,
                base_target_catalog=(
                    pgs_resource_stage / "combined_targets.tsv"
                    if pgs_resource_stage is not None
                    else None
                ),
                base_target_label=(
                    "configured and PGS target union"
                    if pgs_resource_stage is not None
                    else "configured target catalogue"
                ),
            )
        _runtime(runtime_path, run_key, "pgx_resources_complete", 93)

    target_callability: dict[str, Any] | None = None
    target_callability_stage: Path | None = None
    target_config = config
    if pgx_resources is not None and pgx_resource_stage is not None:
        target_config = replace(
            config,
            target_catalog_enabled=True,
            target_catalog_path=pgx_resource_stage / "combined_targets.tsv",
            target_catalog_expected_sha256=pgx_resources["combined_target_catalog_sha256"],
        )
    elif pgs_resources is not None and pgs_resource_stage is not None:
        target_config = replace(
            config,
            target_catalog_enabled=True,
            target_catalog_path=pgs_resource_stage / "combined_targets.tsv",
            target_catalog_expected_sha256=pgs_resources["combined_target_catalog_sha256"],
        )
    if target_config.target_catalog_enabled:
        target_callability_stage = run_root / "07-target-callability"
        target_callability = _reuse_stage(
            target_callability_stage, "target_callability_summary.json"
        )
        if target_callability is None:
            _runtime(runtime_path, run_key, "materializing_target_callability", 94)
            target_callability = build_target_callability(
                target_config,
                ingest_dir=ingest_stage,
                variant_dir=variant_stage,
                destination=target_callability_stage,
            )
        _runtime(runtime_path, run_key, "target_callability_complete", 95)

    pgs_score: dict[str, Any] | None = None
    pgs_score_stage: Path | None = None
    if config.pgs_enabled:
        if pgs_resource_stage is None or target_callability_stage is None:
            raise GenomeBuilderError("PGS scoring requires prepared resources and target evidence")
        pgs_score_stage = run_root / "08-pgs-scoring"
        pgs_score = _reuse_stage(pgs_score_stage, "pgs_score_summary.json")
        if pgs_score is None:
            _runtime(runtime_path, run_key, "scoring_polygenic_models", 96)
            pgs_score = score_pgs(
                config,
                resource_dir=pgs_resource_stage,
                target_callability_dir=target_callability_stage,
                variant_dir=variant_stage,
                destination=pgs_score_stage,
            )
        _runtime(runtime_path, run_key, "pgs_scoring_complete", 97)

    pgx_materialization: dict[str, Any] | None = None
    pgx_materialization_stage: Path | None = None
    pharmcat_run: dict[str, Any] | None = None
    pharmcat_stage: Path | None = None
    pgx_projection: dict[str, Any] | None = None
    pgx_projection_stage: Path | None = None
    if config.pgx.enabled:
        if pgx_resource_stage is None or target_callability_stage is None:
            raise GenomeBuilderError("PGx requires prepared resources and target evidence")
        pgx_materialization_stage = run_root / "09-pgx-materialization"
        pgx_materialization = _reuse_stage(
            pgx_materialization_stage, "pgx_materialization_summary.json"
        )
        if pgx_materialization is None:
            _runtime(runtime_path, run_key, "materializing_pgx_input", 97)
            pgx_materialization = materialize_pgx_input(
                config,
                resource_dir=pgx_resource_stage,
                target_callability_dir=target_callability_stage,
                variant_dir=variant_stage,
                destination=pgx_materialization_stage,
            )

        pharmcat_stage = run_root / "10-pharmcat"
        pharmcat_run = _reuse_stage(pharmcat_stage, "pharmcat_run_summary.json")
        if pharmcat_run is None:
            _runtime(runtime_path, run_key, "running_pharmcat_offline", 98)
            pharmcat_run = run_pharmcat(
                config,
                materialization_dir=pgx_materialization_stage,
                resource_dir=pgx_resource_stage,
                destination=pharmcat_stage,
            )

        pgx_projection_stage = run_root / "11-pgx-projection"
        pgx_projection = _reuse_stage(pgx_projection_stage, "pgx_projection_summary.json")
        if pgx_projection is None:
            _runtime(runtime_path, run_key, "projecting_pgx_evidence", 98)
            pgx_projection = project_pgx_results(
                config,
                pharmcat_dir=pharmcat_stage,
                resource_dir=pgx_resource_stage,
                destination=pgx_projection_stage,
            )
        _runtime(runtime_path, run_key, "pgx_projection_complete", 99)

    _runtime(runtime_path, run_key, "assembling_and_validating_bundle", 99)
    validation = assemble_bundle(
        config,
        run_key=run_key,
        preflight=preflight,
        ingest_dir=ingest_stage,
        ingest_summary=ingest,
        normalization_summary=normalization,
        vep_summary=vep_summary,
        variant_dir=variant_stage,
        variant_summary=variants,
        target_callability_dir=target_callability_stage,
        target_callability_summary=target_callability,
        pgs_score_dir=pgs_score_stage,
        pgs_score_summary=pgs_score,
        pgx_materialization_dir=pgx_materialization_stage,
        pgx_materialization_summary=pgx_materialization,
        pharmcat_dir=pharmcat_stage,
        pharmcat_summary=pharmcat_run,
        pgx_projection_dir=pgx_projection_stage,
        pgx_projection_summary=pgx_projection,
        family_evidence_dir=family_evidence_stage,
        family_evidence_summary=family_evidence,
    )
    _runtime(
        runtime_path,
        run_key,
        "complete",
        100,
        bundle_sha256=directory_digest(config.bundle_path),
    )
    return {
        "status": "PASS",
        "reused": False,
        "bundle": str(config.bundle_path),
        "run_key": run_key,
        "validation": validation,
    }
