"""Strict YAML configuration loading and path resolution."""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import yaml

from .constants import (
    ACMG_SF_SNAPSHOT,
    ACMG_SF_VERSION,
    ALPHAMISSENSE_LICENSE,
    CADD_LICENSE,
    CLINGEN_GENE_LEVEL_SNAPSHOT,
    CLINGEN_LICENSE,
    CLINVAR_SNAPSHOT,
    GENOME_BUILD,
    GENOME_SPEC_COMMIT,
    GWAS_CATALOG_ARCHIVE,
    GWAS_CATALOG_ARCHIVE_SHA256,
    GWAS_CATALOG_LICENSE,
    GWAS_CATALOG_SNAPSHOT,
    LOFTEE_LICENSE,
    PHARMCAT_DOCKER_IMAGE,
    PHARMCAT_POSITIONS_SHA256,
    PHARMCAT_RELEASE_COMMIT,
    PHARMCAT_VERSION,
    REVEL_LICENSE,
    SPLICEAI_LICENSE,
    VEP_CACHE_VERSION,
    VEP_DOCKER_AMD64,
    VEP_VERSION,
)
from .errors import ConfigurationError


def _only(mapping: dict[str, Any], allowed: set[str], section: str) -> None:
    unknown = sorted(set(mapping) - allowed)
    if unknown:
        raise ConfigurationError(f"unknown keys in {section}: {', '.join(unknown)}")


def _mapping(value: Any, section: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigurationError(f"{section} must be a mapping")
    return value


def _required_text(mapping: dict[str, Any], key: str, section: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{section}.{key} must be a non-empty string")
    return value


def _path(value: str, base: Path) -> Path:
    candidate = Path(value).expanduser()
    return candidate if candidate.is_absolute() else (base / candidate).resolve()


def _sha256(value: Any, field: str, *, required: bool) -> str | None:
    if value is None and not required:
        return None
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", value):
        qualifier = "a" if required else "null or a"
        raise ConfigurationError(f"{field} must be {qualifier} SHA-256 hex digest")
    return value.lower()


def _vcf_tag(value: Any, field: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.]*", value):
        raise ConfigurationError(f"{field} must be null or a simple VCF tag identifier")
    return value


def _string_tuple(value: Any, field: str) -> tuple[str, ...]:
    if not isinstance(value, list) or any(
        not isinstance(item, str) or not item.strip() for item in value
    ):
        raise ConfigurationError(f"{field} must be a list of non-empty strings")
    normalized = tuple(item.strip() for item in value)
    if len(set(normalized)) != len(normalized):
        raise ConfigurationError(f"{field} contains duplicate values")
    return normalized


@dataclass(frozen=True)
class QualityThresholds:
    callable_min_dp: int
    callable_min_gq: int
    medium_min_dp: int
    medium_min_gq: int
    high_min_dp: int
    high_min_gq: int

    def confidence(self, dp: int | None, gq: int | None) -> str:
        if dp is not None and gq is not None:
            if dp >= self.high_min_dp and gq >= self.high_min_gq:
                return "high"
            if dp >= self.medium_min_dp and gq >= self.medium_min_gq:
                return "medium"
        return "low"

    def callable(self, dp: int | None, gq: int | None) -> bool:
        return (
            dp is not None
            and gq is not None
            and dp >= self.callable_min_dp
            and gq >= self.callable_min_gq
        )


@dataclass(frozen=True)
class PgsCalibrationConfig:
    path: Path
    expected_sha256: str


@dataclass(frozen=True)
class PgsModelConfig:
    pgs_id: str
    scoring_file: Path
    expected_sha256: str
    source_version: str
    training_date: date
    trait: str
    model_type: str
    license: str
    redistribution_allowed: bool
    max_missing_weight_fraction: float
    calibration: PgsCalibrationConfig | None


@dataclass(frozen=True)
class PgxConfig:
    enabled: bool
    version: str
    release_commit: str
    positions_vcf: Path
    positions_expected_sha256: str
    container_image: str
    license: str
    redistribution_allowed: bool
    cpic_pairs: Path
    cpic_pairs_expected_sha256: str
    cpic_source_version: str
    cpic_snapshot_date: date
    cpic_license: str
    cpic_redistribution_allowed: bool
    docker_executable: str
    timeout_seconds: int


@dataclass(frozen=True)
class FamilyEvidenceTagConfig:
    """Semantic projection of caller-specific family VCF tags."""

    confidence_format: str | None
    reasons_format: str | None
    caller_concordance_format: str | None
    technical_replicate_concordance_format: str | None
    de_novo_score_format: str | None
    phase_set_format: str | None
    difficult_info_flags: tuple[str, ...]
    de_novo_confidence_values: tuple[str, ...]
    conflict_confidence_values: tuple[str, ...]
    parental_mosaicism_reason_tokens: tuple[str, ...]
    parental_dropout_reason_tokens: tuple[str, ...]
    child_dropout_reason_tokens: tuple[str, ...]


@dataclass(frozen=True)
class FamilyLocoConfig:
    """Optional child-independent leave-one-out evidence boundary."""

    review_tsv: Path
    review_expected_sha256: str
    summary_json: Path
    summary_expected_sha256: str


@dataclass(frozen=True)
class FamilyEvidenceConfig:
    """Optional immutable pedigree-aware evidence source."""

    enabled: bool
    joint_vcf: Path | None
    joint_vcf_expected_sha256: str | None
    pedigree: Path | None
    pedigree_expected_sha256: str | None
    source_version: str
    snapshot_date: date
    tags: FamilyEvidenceTagConfig
    loco: FamilyLocoConfig | None


@dataclass(frozen=True)
class GwasCatalogConfig:
    """Immutable literature-curated GWAS association resource boundary."""

    enabled: bool
    snapshot: str
    archive: Path
    archive_expected_sha256: str
    normalized_parquet: Path
    manifest: Path
    license: str
    redistribution_allowed: bool


@dataclass(frozen=True)
class ClingenConfig:
    """Immutable CC0 ClinGen gene-level evidence boundary."""

    enabled: bool
    snapshot: str
    root: Path
    gene_disease_parquet: Path
    dosage_parquet: Path
    manifest: Path
    license: str
    redistribution_allowed: bool


@dataclass(frozen=True)
class PredictorPluginConfig:
    """Licence and enablement boundary for one offline VEP plugin."""

    enabled: bool
    license: str
    redistribution_allowed: bool


@dataclass(frozen=True)
class PredictorResourcesConfig:
    """Private, checksum-pinned resources used by optional VEP plugins."""

    root: Path
    manifest: Path
    loftee: PredictorPluginConfig
    cadd: PredictorPluginConfig
    revel: PredictorPluginConfig
    alphamissense: PredictorPluginConfig
    spliceai: PredictorPluginConfig

    @property
    def enabled(self) -> bool:
        return bool(self.enabled_plugins)

    @property
    def enabled_plugins(self) -> tuple[str, ...]:
        return tuple(
            name
            for name in ("loftee", "cadd", "revel", "alphamissense", "spliceai")
            if getattr(self, name).enabled
        )


@dataclass(frozen=True)
class BuildConfig:
    config_path: Path
    output_dir: Path
    work_dir: Path
    bundle_name: str
    input_vcf: Path
    sample: str
    genome_build: str
    reference_fasta: Path
    reference_expected_sha256: str | None
    contig_style: str
    resource_cache: Path
    allow_download: bool
    genome_spec_commit: str
    bcftools: str
    check_ref: str
    quality: QualityThresholds
    vep_enabled: bool
    vep_mode: str
    vep_executable: str
    vep_container_image: str
    vep_plugin_dir: str
    vep_cache_dir: Path
    vep_cache_version: str
    vep_cache_manifest: Path
    vep_expected_version: str
    vep_forks: int
    vep_buffer_size: int
    clinvar_enabled: bool
    clinvar_snapshot: str
    clinvar_vcf: Path
    clinvar_manifest: Path
    acmg_sf_enabled: bool
    acmg_sf_version: str
    acmg_sf_snapshot: str
    acmg_sf_dataset: Path
    acmg_sf_manifest: Path
    acmg_sf_min_clinvar_review_stars: int
    pgs_enabled: bool
    pgs_models: tuple[PgsModelConfig, ...]
    pgx: PgxConfig
    family_evidence: FamilyEvidenceConfig
    gwas: GwasCatalogConfig
    clingen: ClingenConfig
    predictors: PredictorResourcesConfig
    target_catalog_enabled: bool
    target_catalog_path: Path
    target_catalog_expected_sha256: str | None
    target_catalog_assay_scope: str
    agent_exploration_profile: bool
    reference_confidence_extension: bool
    partition_variants_by_chromosome: bool
    overwrite: bool
    include_input_paths: bool
    include_sample_identifier: bool

    @property
    def bundle_path(self) -> Path:
        return self.output_dir / self.bundle_name


def load_config(path: str | Path) -> BuildConfig:
    config_path = Path(path).expanduser().resolve()
    try:
        document = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ConfigurationError(f"configuration not found: {config_path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigurationError(f"invalid YAML in {config_path}: {exc}") from exc
    root = _mapping(document, "configuration")
    _only(
        root,
        {
            "project",
            "input",
            "reference",
            "resources",
            "normalization",
            "quality",
            "annotation",
            "callability",
            "output",
            "privacy",
        },
        "configuration",
    )
    base = config_path.parent

    project = _mapping(root.get("project"), "project")
    _only(project, {"output_dir", "work_dir", "bundle_name"}, "project")
    output_dir = _path(_required_text(project, "output_dir", "project"), base)
    work_dir = _path(_required_text(project, "work_dir", "project"), base)
    bundle_name = _required_text(project, "bundle_name", "project")
    if Path(bundle_name).name != bundle_name or not bundle_name.endswith(".genome"):
        raise ConfigurationError("project.bundle_name must be one safe basename ending in .genome")

    input_section = _mapping(root.get("input"), "input")
    _only(input_section, {"vcf", "sample", "genome_build"}, "input")
    input_vcf = _path(_required_text(input_section, "vcf", "input"), base)
    sample = _required_text(input_section, "sample", "input")
    genome_build = input_section.get("genome_build", GENOME_BUILD)
    if genome_build != GENOME_BUILD:
        raise ConfigurationError(f"only {GENOME_BUILD} is supported, got {genome_build!r}")

    reference = _mapping(root.get("reference"), "reference")
    _only(reference, {"fasta", "expected_sha256", "contig_style"}, "reference")
    reference_fasta = _path(_required_text(reference, "fasta", "reference"), base)
    expected_sha = reference.get("expected_sha256")
    if expected_sha is not None and (
        not isinstance(expected_sha, str)
        or len(expected_sha) != 64
        or any(char not in "0123456789abcdef" for char in expected_sha.lower())
    ):
        raise ConfigurationError("reference.expected_sha256 must be null or a SHA-256 hex digest")
    contig_style = reference.get("contig_style", "chr")
    if contig_style not in {"chr", "bare"}:
        raise ConfigurationError("reference.contig_style must be chr or bare")

    resources = _mapping(root.get("resources"), "resources")
    _only(resources, {"cache_dir", "allow_download", "genome_spec_commit"}, "resources")
    resource_cache = _path(_required_text(resources, "cache_dir", "resources"), base)
    allow_download = bool(resources.get("allow_download", False))
    genome_spec_commit = resources.get("genome_spec_commit", GENOME_SPEC_COMMIT)
    if genome_spec_commit != GENOME_SPEC_COMMIT:
        raise ConfigurationError(
            "this release requires genome_spec_commit "
            f"{GENOME_SPEC_COMMIT}, got {genome_spec_commit}"
        )

    normalization = _mapping(root.get("normalization", {}), "normalization")
    _only(normalization, {"bcftools", "check_ref"}, "normalization")
    bcftools = str(normalization.get("bcftools", "bcftools"))
    check_ref = str(normalization.get("check_ref", "error"))
    if check_ref not in {"error", "warn"}:
        raise ConfigurationError("normalization.check_ref must be error or warn")

    quality_section = _mapping(root.get("quality", {}), "quality")
    quality_keys = {
        "callable_min_dp",
        "callable_min_gq",
        "medium_min_dp",
        "medium_min_gq",
        "high_min_dp",
        "high_min_gq",
    }
    _only(quality_section, quality_keys, "quality")
    defaults = {
        "callable_min_dp": 10,
        "callable_min_gq": 20,
        "medium_min_dp": 10,
        "medium_min_gq": 20,
        "high_min_dp": 20,
        "high_min_gq": 30,
    }
    values: dict[str, int] = {}
    for key, default in defaults.items():
        raw = quality_section.get(key, default)
        if not isinstance(raw, int) or isinstance(raw, bool) or raw < 0:
            raise ConfigurationError(f"quality.{key} must be a non-negative integer")
        values[key] = raw
    quality = QualityThresholds(**values)

    annotation = _mapping(root.get("annotation", {}), "annotation")
    _only(
        annotation,
        {
            "vep",
            "predictors",
            "clinvar",
            "acmg_sf",
            "pgs",
            "pgx",
            "family_evidence",
            "gwas",
            "clingen",
        },
        "annotation",
    )
    vep = _mapping(annotation.get("vep", {}), "annotation.vep")
    _only(
        vep,
        {
            "enabled",
            "mode",
            "executable",
            "container_image",
            "plugin_dir",
            "cache_dir",
            "cache_version",
            "cache_manifest",
            "expected_version",
            "forks",
            "buffer_size",
        },
        "annotation.vep",
    )
    vep_enabled = bool(vep.get("enabled", False))
    vep_mode = str(vep.get("mode", "docker"))
    if vep_mode not in {"docker", "executable"}:
        raise ConfigurationError("annotation.vep.mode must be docker or executable")
    vep_executable = str(vep.get("executable", "vep"))
    if not vep_executable:
        raise ConfigurationError("annotation.vep.executable cannot be empty")
    vep_container_image = str(vep.get("container_image", VEP_DOCKER_AMD64))
    vep_plugin_dir = str(vep.get("plugin_dir", "/plugins" if vep_mode == "docker" else "")).strip()
    vep_cache_version = str(vep.get("cache_version", VEP_CACHE_VERSION))
    default_vep_root = resource_cache / "vep" / vep_cache_version
    vep_cache_dir = _path(str(vep.get("cache_dir", default_vep_root / "cache")), base)
    vep_cache_manifest = _path(
        str(vep.get("cache_manifest", vep_cache_dir.parent / "cache_resource_manifest.json")),
        base,
    )
    vep_expected_version = str(vep.get("expected_version", VEP_VERSION))
    vep_forks = vep.get("forks", 8)
    vep_buffer_size = vep.get("buffer_size", 5000)
    for key, value in {"forks": vep_forks, "buffer_size": vep_buffer_size}.items():
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ConfigurationError(f"annotation.vep.{key} must be a positive integer")
    if vep_forks > 128:
        raise ConfigurationError("annotation.vep.forks cannot exceed 128")
    if vep_enabled and vep_cache_version != VEP_CACHE_VERSION:
        raise ConfigurationError(
            f"this release requires VEP cache {VEP_CACHE_VERSION}, got {vep_cache_version!r}"
        )
    if vep_enabled and vep_expected_version != VEP_VERSION:
        raise ConfigurationError(
            f"this release requires VEP {VEP_VERSION}, got {vep_expected_version!r}"
        )
    if (
        vep_enabled
        and vep_mode == "docker"
        and not re.fullmatch(r"[^\s@]+@sha256:[0-9a-f]{64}", vep_container_image)
    ):
        raise ConfigurationError(
            "enabled Docker VEP requires container_image pinned as name@sha256:<64 hex>"
        )

    predictor_section = _mapping(annotation.get("predictors", {}), "annotation.predictors")
    predictor_names = ("loftee", "cadd", "revel", "alphamissense", "spliceai")
    _only(
        predictor_section,
        {"root", "manifest", *predictor_names},
        "annotation.predictors",
    )
    predictor_root = _path(str(predictor_section.get("root", resource_cache / "predictors")), base)
    predictor_manifest = _path(
        str(predictor_section.get("manifest", predictor_root / "predictor_resource_manifest.json")),
        base,
    )
    expected_predictor_licenses = {
        "loftee": LOFTEE_LICENSE,
        "cadd": CADD_LICENSE,
        "revel": REVEL_LICENSE,
        "alphamissense": ALPHAMISSENSE_LICENSE,
        "spliceai": SPLICEAI_LICENSE,
    }
    predictor_plugins: dict[str, PredictorPluginConfig] = {}
    for predictor_name in predictor_names:
        section_name = f"annotation.predictors.{predictor_name}"
        raw_plugin = _mapping(predictor_section.get(predictor_name, {}), section_name)
        _only(
            raw_plugin,
            {"enabled", "license", "redistribution_allowed"},
            section_name,
        )
        enabled_value = raw_plugin.get("enabled", False)
        redistribution_value = raw_plugin.get("redistribution_allowed", False)
        if not isinstance(enabled_value, bool):
            raise ConfigurationError(f"{section_name}.enabled must be boolean")
        if not isinstance(redistribution_value, bool):
            raise ConfigurationError(f"{section_name}.redistribution_allowed must be boolean")
        license_value = str(
            raw_plugin.get("license", expected_predictor_licenses[predictor_name])
        ).strip()
        if enabled_value and not license_value:
            raise ConfigurationError(f"{section_name}.license must be explicit")
        if enabled_value and license_value != expected_predictor_licenses[predictor_name]:
            raise ConfigurationError(
                f"{section_name}.license must match the pinned release boundary"
            )
        if redistribution_value:
            raise ConfigurationError(
                f"{section_name}.redistribution_allowed cannot be asserted by this release; "
                "the public repository distributes code only"
            )
        predictor_plugins[predictor_name] = PredictorPluginConfig(
            enabled=enabled_value,
            license=license_value,
            redistribution_allowed=redistribution_value,
        )
    predictor_config = PredictorResourcesConfig(
        root=predictor_root,
        manifest=predictor_manifest,
        loftee=predictor_plugins["loftee"],
        cadd=predictor_plugins["cadd"],
        revel=predictor_plugins["revel"],
        alphamissense=predictor_plugins["alphamissense"],
        spliceai=predictor_plugins["spliceai"],
    )
    if predictor_config.enabled and not vep_enabled:
        raise ConfigurationError("annotation.predictors requires annotation.vep.enabled=true")
    if predictor_config.enabled and not vep_plugin_dir:
        raise ConfigurationError(
            "annotation.vep.plugin_dir is required when predictors are enabled"
        )
    if predictor_config.enabled and vep_mode != "docker":
        raise ConfigurationError(
            "annotation.predictors requires the pinned Docker VEP runtime in this release"
        )
    if predictor_config.enabled and vep_plugin_dir != "/plugins":
        raise ConfigurationError(
            "annotation.predictors requires the pinned container plugin directory /plugins"
        )
    if predictor_config.enabled and any(
        delimiter in str(predictor_root) for delimiter in (",", ":")
    ):
        raise ConfigurationError(
            "annotation.predictors.root cannot contain comma or colon because VEP plugin "
            "parameters use those delimiters"
        )

    clinvar = _mapping(annotation.get("clinvar", {}), "annotation.clinvar")
    _only(clinvar, {"enabled", "snapshot", "vcf", "manifest"}, "annotation.clinvar")
    clinvar_enabled = bool(clinvar.get("enabled", False))
    clinvar_snapshot = str(clinvar.get("snapshot", CLINVAR_SNAPSHOT))
    default_clinvar_root = resource_cache / "clinvar" / clinvar_snapshot
    clinvar_vcf = _path(
        str(clinvar.get("vcf", default_clinvar_root / "normalized" / "clinvar.GRCh38.vcf.gz")),
        base,
    )
    clinvar_manifest = _path(
        str(clinvar.get("manifest", default_clinvar_root / "resource_manifest.json")), base
    )
    if clinvar_enabled and clinvar_snapshot != CLINVAR_SNAPSHOT:
        raise ConfigurationError(
            f"this release requires ClinVar snapshot {CLINVAR_SNAPSHOT}, got {clinvar_snapshot!r}"
        )

    acmg_sf = _mapping(annotation.get("acmg_sf", {}), "annotation.acmg_sf")
    _only(
        acmg_sf,
        {
            "enabled",
            "version",
            "snapshot",
            "dataset",
            "manifest",
            "min_clinvar_review_stars",
        },
        "annotation.acmg_sf",
    )
    acmg_sf_enabled = bool(acmg_sf.get("enabled", False))
    acmg_sf_version = str(acmg_sf.get("version", ACMG_SF_VERSION))
    acmg_sf_snapshot = str(acmg_sf.get("snapshot", ACMG_SF_SNAPSHOT))
    default_acmg_root = resource_cache / "acmg_sf" / acmg_sf_version
    acmg_sf_dataset = _path(
        str(
            acmg_sf.get(
                "dataset",
                default_acmg_root / f"acmg_sf_v{acmg_sf_version}_{acmg_sf_snapshot}.json",
            )
        ),
        base,
    )
    acmg_sf_manifest = _path(
        str(acmg_sf.get("manifest", default_acmg_root / "resource_manifest.json")), base
    )
    acmg_sf_min_review_stars = acmg_sf.get("min_clinvar_review_stars", 2)
    if (
        not isinstance(acmg_sf_min_review_stars, int)
        or isinstance(acmg_sf_min_review_stars, bool)
        or not 0 <= acmg_sf_min_review_stars <= 4
    ):
        raise ConfigurationError(
            "annotation.acmg_sf.min_clinvar_review_stars must be an integer from 0 to 4"
        )
    if acmg_sf_enabled and acmg_sf_version != ACMG_SF_VERSION:
        raise ConfigurationError(
            f"this release requires ACMG SF {ACMG_SF_VERSION}, got {acmg_sf_version!r}"
        )
    if acmg_sf_enabled and acmg_sf_snapshot != ACMG_SF_SNAPSHOT:
        raise ConfigurationError(
            "this release requires ACMG/ClinGen snapshot "
            f"{ACMG_SF_SNAPSHOT}, got {acmg_sf_snapshot!r}"
        )
    if acmg_sf_enabled and not (vep_enabled and clinvar_enabled):
        raise ConfigurationError(
            "annotation.acmg_sf requires both annotation.vep and annotation.clinvar"
        )

    clingen_section = _mapping(annotation.get("clingen", {}), "annotation.clingen")
    _only(
        clingen_section,
        {
            "enabled",
            "snapshot",
            "root",
            "gene_disease_parquet",
            "dosage_parquet",
            "manifest",
            "license",
            "redistribution_allowed",
        },
        "annotation.clingen",
    )
    clingen_enabled = clingen_section.get("enabled", False)
    if not isinstance(clingen_enabled, bool):
        raise ConfigurationError("annotation.clingen.enabled must be boolean")
    clingen_snapshot = str(clingen_section.get("snapshot", CLINGEN_GENE_LEVEL_SNAPSHOT)).strip()
    default_clingen_root = resource_cache / "clingen" / clingen_snapshot
    clingen_root = _path(str(clingen_section.get("root", default_clingen_root)), base)
    clingen_gene_disease = _path(
        str(
            clingen_section.get(
                "gene_disease_parquet",
                clingen_root / "normalized" / "gene_disease_validity.parquet",
            )
        ),
        base,
    )
    clingen_dosage = _path(
        str(
            clingen_section.get(
                "dosage_parquet",
                clingen_root / "normalized" / "dosage_sensitivity.parquet",
            )
        ),
        base,
    )
    clingen_manifest = _path(
        str(clingen_section.get("manifest", clingen_root / "resource_manifest.json")), base
    )
    clingen_license = str(clingen_section.get("license", CLINGEN_LICENSE)).strip()
    clingen_redistribution = clingen_section.get("redistribution_allowed", True)
    if not isinstance(clingen_redistribution, bool):
        raise ConfigurationError("annotation.clingen.redistribution_allowed must be boolean")
    if clingen_enabled:
        if clingen_snapshot != CLINGEN_GENE_LEVEL_SNAPSHOT:
            raise ConfigurationError(
                "this release requires ClinGen snapshot "
                f"{CLINGEN_GENE_LEVEL_SNAPSHOT}, got {clingen_snapshot!r}"
            )
        if clingen_license != CLINGEN_LICENSE:
            raise ConfigurationError(
                "annotation.clingen.license must match the pinned CC0 release boundary"
            )
        if not clingen_redistribution:
            raise ConfigurationError(
                "the pinned ClinGen snapshot is CC0 and its projected bundle tables are "
                "redistributable; redistribution_allowed must remain true"
            )
        if not vep_enabled:
            raise ConfigurationError(
                "annotation.clingen requires annotation.vep.enabled=true for exact gene mapping"
            )
        expected_clingen_paths = {
            clingen_gene_disease.parent.parent,
            clingen_dosage.parent.parent,
            clingen_manifest.parent,
        }
        if expected_clingen_paths != {clingen_root}:
            raise ConfigurationError(
                "annotation.clingen normalized tables and manifest must be below its root"
            )
    clingen_config = ClingenConfig(
        enabled=clingen_enabled,
        snapshot=clingen_snapshot,
        root=clingen_root,
        gene_disease_parquet=clingen_gene_disease,
        dosage_parquet=clingen_dosage,
        manifest=clingen_manifest,
        license=clingen_license,
        redistribution_allowed=clingen_redistribution,
    )

    gwas_section = _mapping(annotation.get("gwas", {}), "annotation.gwas")
    _only(
        gwas_section,
        {
            "enabled",
            "snapshot",
            "archive",
            "archive_expected_sha256",
            "normalized_parquet",
            "manifest",
            "license",
            "redistribution_allowed",
        },
        "annotation.gwas",
    )
    gwas_enabled = bool(gwas_section.get("enabled", False))
    gwas_snapshot = str(gwas_section.get("snapshot", GWAS_CATALOG_SNAPSHOT)).strip()
    default_gwas_root = resource_cache / "gwas_catalog" / gwas_snapshot
    gwas_archive = _path(
        str(
            gwas_section.get(
                "archive",
                default_gwas_root / "archive" / GWAS_CATALOG_ARCHIVE,
            )
        ),
        base,
    )
    gwas_archive_sha = str(
        gwas_section.get("archive_expected_sha256", GWAS_CATALOG_ARCHIVE_SHA256)
    ).lower()
    gwas_normalized = _path(
        str(
            gwas_section.get(
                "normalized_parquet",
                default_gwas_root / "normalized" / "gwas_catalog.parquet",
            )
        ),
        base,
    )
    gwas_manifest = _path(
        str(gwas_section.get("manifest", default_gwas_root / "resource_manifest.json")),
        base,
    )
    gwas_license = str(gwas_section.get("license", GWAS_CATALOG_LICENSE)).strip()
    gwas_redistribution = gwas_section.get("redistribution_allowed", False)
    if gwas_enabled:
        if gwas_snapshot != GWAS_CATALOG_SNAPSHOT:
            raise ConfigurationError(
                "this release requires GWAS Catalog snapshot "
                f"{GWAS_CATALOG_SNAPSHOT}, got {gwas_snapshot!r}"
            )
        if gwas_archive_sha != GWAS_CATALOG_ARCHIVE_SHA256:
            raise ConfigurationError(
                "annotation.gwas.archive_expected_sha256 does not match the pinned release"
            )
        if not gwas_license:
            raise ConfigurationError("annotation.gwas.license must be explicit")
        if not isinstance(gwas_redistribution, bool):
            raise ConfigurationError("annotation.gwas.redistribution_allowed must be boolean")
    gwas_config = GwasCatalogConfig(
        enabled=gwas_enabled,
        snapshot=gwas_snapshot,
        archive=gwas_archive,
        archive_expected_sha256=gwas_archive_sha,
        normalized_parquet=gwas_normalized,
        manifest=gwas_manifest,
        license=gwas_license,
        redistribution_allowed=(
            bool(gwas_redistribution) if isinstance(gwas_redistribution, bool) else False
        ),
    )

    pgs = _mapping(annotation.get("pgs", {}), "annotation.pgs")
    _only(pgs, {"enabled", "models"}, "annotation.pgs")
    pgs_enabled = bool(pgs.get("enabled", False))
    raw_models = pgs.get("models", [])
    if not isinstance(raw_models, list):
        raise ConfigurationError("annotation.pgs.models must be a list")
    pgs_models: list[PgsModelConfig] = []
    seen_pgs_ids: set[str] = set()
    model_keys = {
        "pgs_id",
        "scoring_file",
        "expected_sha256",
        "source_version",
        "training_date",
        "trait",
        "model_type",
        "license",
        "redistribution_allowed",
        "max_missing_weight_fraction",
        "calibration",
    }
    for index, raw_model in enumerate(raw_models):
        section = f"annotation.pgs.models[{index}]"
        model = _mapping(raw_model, section)
        _only(model, model_keys, section)
        pgs_id = _required_text(model, "pgs_id", section).upper()
        if not re.fullmatch(r"PGS[0-9]{6}", pgs_id):
            raise ConfigurationError(f"{section}.pgs_id must match PGS followed by six digits")
        if pgs_id in seen_pgs_ids:
            raise ConfigurationError(f"duplicate annotation.pgs model: {pgs_id}")
        seen_pgs_ids.add(pgs_id)
        scoring_file = _path(_required_text(model, "scoring_file", section), base)
        scoring_sha = _required_text(model, "expected_sha256", section).lower()
        if not re.fullmatch(r"[0-9a-f]{64}", scoring_sha):
            raise ConfigurationError(f"{section}.expected_sha256 must be a SHA-256 hex digest")
        source_version = _required_text(model, "source_version", section)
        training_date_text = _required_text(model, "training_date", section)
        try:
            training_date = date.fromisoformat(training_date_text)
        except ValueError as exc:
            raise ConfigurationError(f"{section}.training_date must be YYYY-MM-DD") from exc
        trait = _required_text(model, "trait", section)
        model_type = str(model.get("model_type", "linear")).strip()
        if model_type != "linear":
            raise ConfigurationError(
                f"{section}.model_type must be linear; non-additive models are not yet supported"
            )
        license_text = _required_text(model, "license", section)
        redistribution_allowed = model.get("redistribution_allowed")
        if not isinstance(redistribution_allowed, bool):
            raise ConfigurationError(f"{section}.redistribution_allowed must be boolean")
        max_missing = model.get("max_missing_weight_fraction", 0.0)
        if (
            not isinstance(max_missing, (int, float))
            or isinstance(max_missing, bool)
            or not 0.0 <= float(max_missing) <= 1.0
        ):
            raise ConfigurationError(
                f"{section}.max_missing_weight_fraction must be between 0 and 1"
            )
        calibration_raw = model.get("calibration")
        calibration: PgsCalibrationConfig | None = None
        if calibration_raw is not None:
            calibration_section = f"{section}.calibration"
            calibration_mapping = _mapping(calibration_raw, calibration_section)
            _only(
                calibration_mapping,
                {"file", "expected_sha256"},
                calibration_section,
            )
            calibration_path = _path(
                _required_text(calibration_mapping, "file", calibration_section), base
            )
            calibration_sha = _required_text(
                calibration_mapping, "expected_sha256", calibration_section
            ).lower()
            if not re.fullmatch(r"[0-9a-f]{64}", calibration_sha):
                raise ConfigurationError(
                    f"{calibration_section}.expected_sha256 must be a SHA-256 hex digest"
                )
            calibration = PgsCalibrationConfig(calibration_path, calibration_sha)
        pgs_models.append(
            PgsModelConfig(
                pgs_id=pgs_id,
                scoring_file=scoring_file,
                expected_sha256=scoring_sha,
                source_version=source_version,
                training_date=training_date,
                trait=trait,
                model_type=model_type,
                license=license_text,
                redistribution_allowed=redistribution_allowed,
                max_missing_weight_fraction=float(max_missing),
                calibration=calibration,
            )
        )
    if pgs_enabled and not pgs_models:
        raise ConfigurationError("enabled annotation.pgs requires at least one model")

    pgx_section = _mapping(annotation.get("pgx", {}), "annotation.pgx")
    _only(
        pgx_section,
        {"enabled", "pharmcat", "cpic", "docker_executable", "timeout_seconds"},
        "annotation.pgx",
    )
    pgx_enabled = bool(pgx_section.get("enabled", False))
    pharmcat = _mapping(pgx_section.get("pharmcat", {}), "annotation.pgx.pharmcat")
    _only(
        pharmcat,
        {
            "version",
            "release_commit",
            "positions_vcf",
            "expected_sha256",
            "container_image",
            "license",
            "redistribution_allowed",
        },
        "annotation.pgx.pharmcat",
    )
    pharmcat_version = str(pharmcat.get("version", PHARMCAT_VERSION))
    pharmcat_release_commit = str(pharmcat.get("release_commit", PHARMCAT_RELEASE_COMMIT))
    pharmcat_positions_vcf = _path(
        str(
            pharmcat.get(
                "positions_vcf",
                resource_cache
                / "pharmcat"
                / PHARMCAT_VERSION
                / f"pharmcat_positions_{PHARMCAT_VERSION}.vcf.bgz",
            )
        ),
        base,
    )
    pharmcat_positions_sha = str(pharmcat.get("expected_sha256", PHARMCAT_POSITIONS_SHA256)).lower()
    pharmcat_container = str(pharmcat.get("container_image", PHARMCAT_DOCKER_IMAGE))
    pharmcat_license = str(pharmcat.get("license", "")).strip()
    pharmcat_redistribution = pharmcat.get("redistribution_allowed")

    cpic = _mapping(pgx_section.get("cpic", {}), "annotation.pgx.cpic")
    _only(
        cpic,
        {
            "pairs_json",
            "expected_sha256",
            "source_version",
            "snapshot_date",
            "license",
            "redistribution_allowed",
        },
        "annotation.pgx.cpic",
    )
    cpic_pairs = _path(
        str(cpic.get("pairs_json", resource_cache / "cpic" / "cpic_pairs.json")), base
    )
    cpic_sha = str(cpic.get("expected_sha256", "")).lower()
    cpic_source_version = str(cpic.get("source_version", "")).strip()
    cpic_snapshot_text = str(cpic.get("snapshot_date", "1970-01-01"))
    try:
        cpic_snapshot_date = date.fromisoformat(cpic_snapshot_text)
    except ValueError as exc:
        raise ConfigurationError("annotation.pgx.cpic.snapshot_date must be YYYY-MM-DD") from exc
    cpic_license = str(cpic.get("license", "")).strip()
    cpic_redistribution = cpic.get("redistribution_allowed")
    docker_executable = str(pgx_section.get("docker_executable", "docker")).strip()
    timeout_seconds = pgx_section.get("timeout_seconds", 1800)
    if (
        not isinstance(timeout_seconds, int)
        or isinstance(timeout_seconds, bool)
        or not 60 <= timeout_seconds <= 86_400
    ):
        raise ConfigurationError(
            "annotation.pgx.timeout_seconds must be an integer from 60 to 86400"
        )
    if pgx_enabled:
        if pharmcat_version != PHARMCAT_VERSION:
            raise ConfigurationError(
                f"this release requires PharmCAT {PHARMCAT_VERSION}, got {pharmcat_version!r}"
            )
        if pharmcat_release_commit != PHARMCAT_RELEASE_COMMIT:
            raise ConfigurationError(
                "annotation.pgx.pharmcat.release_commit does not match the pinned release"
            )
        if pharmcat_positions_sha != PHARMCAT_POSITIONS_SHA256:
            raise ConfigurationError(
                "annotation.pgx.pharmcat.expected_sha256 does not match the pinned positions asset"
            )
        if pharmcat_container != PHARMCAT_DOCKER_IMAGE:
            raise ConfigurationError(
                "annotation.pgx.pharmcat.container_image must equal the pinned digest"
            )
        if not pharmcat_license:
            raise ConfigurationError("annotation.pgx.pharmcat.license must be explicit")
        if not isinstance(pharmcat_redistribution, bool):
            raise ConfigurationError(
                "annotation.pgx.pharmcat.redistribution_allowed must be boolean"
            )
        if not re.fullmatch(r"[0-9a-f]{64}", cpic_sha):
            raise ConfigurationError(
                "annotation.pgx.cpic.expected_sha256 must be a SHA-256 hex digest"
            )
        if not cpic_source_version:
            raise ConfigurationError("annotation.pgx.cpic.source_version must be explicit")
        if cpic_snapshot_date == date(1970, 1, 1):
            raise ConfigurationError("annotation.pgx.cpic.snapshot_date must be explicit")
        if not cpic_license:
            raise ConfigurationError("annotation.pgx.cpic.license must be explicit")
        if not isinstance(cpic_redistribution, bool):
            raise ConfigurationError("annotation.pgx.cpic.redistribution_allowed must be boolean")
        if not docker_executable:
            raise ConfigurationError("annotation.pgx.docker_executable cannot be empty")
    pgx_config = PgxConfig(
        enabled=pgx_enabled,
        version=pharmcat_version,
        release_commit=pharmcat_release_commit,
        positions_vcf=pharmcat_positions_vcf,
        positions_expected_sha256=pharmcat_positions_sha,
        container_image=pharmcat_container,
        license=pharmcat_license,
        redistribution_allowed=(
            bool(pharmcat_redistribution) if isinstance(pharmcat_redistribution, bool) else False
        ),
        cpic_pairs=cpic_pairs,
        cpic_pairs_expected_sha256=cpic_sha,
        cpic_source_version=cpic_source_version,
        cpic_snapshot_date=cpic_snapshot_date,
        cpic_license=cpic_license,
        cpic_redistribution_allowed=(
            bool(cpic_redistribution) if isinstance(cpic_redistribution, bool) else False
        ),
        docker_executable=docker_executable,
        timeout_seconds=timeout_seconds,
    )

    family_section = _mapping(annotation.get("family_evidence", {}), "annotation.family_evidence")
    _only(
        family_section,
        {
            "enabled",
            "joint_vcf",
            "joint_vcf_expected_sha256",
            "pedigree",
            "pedigree_expected_sha256",
            "source_version",
            "snapshot_date",
            "tags",
            "loco",
        },
        "annotation.family_evidence",
    )
    family_enabled = bool(family_section.get("enabled", False))

    def optional_family_path(key: str) -> Path | None:
        raw = family_section.get(key)
        if raw is None:
            return None
        if not isinstance(raw, str) or not raw.strip():
            raise ConfigurationError(
                f"annotation.family_evidence.{key} must be null or a non-empty path"
            )
        return _path(raw, base)

    family_joint_vcf = optional_family_path("joint_vcf")
    family_pedigree = optional_family_path("pedigree")
    family_joint_sha = _sha256(
        family_section.get("joint_vcf_expected_sha256"),
        "annotation.family_evidence.joint_vcf_expected_sha256",
        required=family_enabled,
    )
    family_pedigree_sha = _sha256(
        family_section.get("pedigree_expected_sha256"),
        "annotation.family_evidence.pedigree_expected_sha256",
        required=family_enabled,
    )
    family_source_version = str(family_section.get("source_version", "")).strip()
    family_snapshot_text = str(family_section.get("snapshot_date", "1970-01-01"))
    try:
        family_snapshot_date = date.fromisoformat(family_snapshot_text)
    except ValueError as exc:
        raise ConfigurationError(
            "annotation.family_evidence.snapshot_date must be YYYY-MM-DD"
        ) from exc

    family_tags_section = _mapping(
        family_section.get("tags", {}), "annotation.family_evidence.tags"
    )
    family_tag_keys = {
        "confidence_format",
        "reasons_format",
        "caller_concordance_format",
        "technical_replicate_concordance_format",
        "de_novo_score_format",
        "phase_set_format",
        "difficult_info_flags",
        "de_novo_confidence_values",
        "conflict_confidence_values",
        "parental_mosaicism_reason_tokens",
        "parental_dropout_reason_tokens",
        "child_dropout_reason_tokens",
    }
    _only(family_tags_section, family_tag_keys, "annotation.family_evidence.tags")
    family_tags = FamilyEvidenceTagConfig(
        confidence_format=_vcf_tag(
            family_tags_section.get("confidence_format"),
            "annotation.family_evidence.tags.confidence_format",
        ),
        reasons_format=_vcf_tag(
            family_tags_section.get("reasons_format"),
            "annotation.family_evidence.tags.reasons_format",
        ),
        caller_concordance_format=_vcf_tag(
            family_tags_section.get("caller_concordance_format"),
            "annotation.family_evidence.tags.caller_concordance_format",
        ),
        technical_replicate_concordance_format=_vcf_tag(
            family_tags_section.get("technical_replicate_concordance_format"),
            "annotation.family_evidence.tags.technical_replicate_concordance_format",
        ),
        de_novo_score_format=_vcf_tag(
            family_tags_section.get("de_novo_score_format"),
            "annotation.family_evidence.tags.de_novo_score_format",
        ),
        phase_set_format=_vcf_tag(
            family_tags_section.get("phase_set_format"),
            "annotation.family_evidence.tags.phase_set_format",
        ),
        difficult_info_flags=_string_tuple(
            family_tags_section.get("difficult_info_flags", []),
            "annotation.family_evidence.tags.difficult_info_flags",
        ),
        de_novo_confidence_values=_string_tuple(
            family_tags_section.get("de_novo_confidence_values", []),
            "annotation.family_evidence.tags.de_novo_confidence_values",
        ),
        conflict_confidence_values=_string_tuple(
            family_tags_section.get("conflict_confidence_values", []),
            "annotation.family_evidence.tags.conflict_confidence_values",
        ),
        parental_mosaicism_reason_tokens=_string_tuple(
            family_tags_section.get("parental_mosaicism_reason_tokens", []),
            "annotation.family_evidence.tags.parental_mosaicism_reason_tokens",
        ),
        parental_dropout_reason_tokens=_string_tuple(
            family_tags_section.get("parental_dropout_reason_tokens", []),
            "annotation.family_evidence.tags.parental_dropout_reason_tokens",
        ),
        child_dropout_reason_tokens=_string_tuple(
            family_tags_section.get("child_dropout_reason_tokens", []),
            "annotation.family_evidence.tags.child_dropout_reason_tokens",
        ),
    )

    loco_section = _mapping(family_section.get("loco", {}), "annotation.family_evidence.loco")
    _only(
        loco_section,
        {
            "review_tsv",
            "review_expected_sha256",
            "summary_json",
            "summary_expected_sha256",
        },
        "annotation.family_evidence.loco",
    )
    family_loco: FamilyLocoConfig | None = None
    if loco_section:
        required_loco_keys = {
            "review_tsv",
            "review_expected_sha256",
            "summary_json",
            "summary_expected_sha256",
        }
        missing_loco = sorted(required_loco_keys - set(loco_section))
        if missing_loco:
            raise ConfigurationError(
                "annotation.family_evidence.loco is all-or-none; missing " + ", ".join(missing_loco)
            )
        review_path = _path(
            _required_text(loco_section, "review_tsv", "annotation.family_evidence.loco"),
            base,
        )
        summary_path = _path(
            _required_text(loco_section, "summary_json", "annotation.family_evidence.loco"),
            base,
        )
        review_sha = _sha256(
            loco_section.get("review_expected_sha256"),
            "annotation.family_evidence.loco.review_expected_sha256",
            required=True,
        )
        summary_sha = _sha256(
            loco_section.get("summary_expected_sha256"),
            "annotation.family_evidence.loco.summary_expected_sha256",
            required=True,
        )
        assert review_sha is not None and summary_sha is not None
        family_loco = FamilyLocoConfig(
            review_tsv=review_path,
            review_expected_sha256=review_sha,
            summary_json=summary_path,
            summary_expected_sha256=summary_sha,
        )
    if family_enabled:
        if family_joint_vcf is None or family_pedigree is None:
            raise ConfigurationError(
                "enabled annotation.family_evidence requires joint_vcf and pedigree"
            )
        if not family_source_version:
            raise ConfigurationError("annotation.family_evidence.source_version must be explicit")
        if family_snapshot_date == date(1970, 1, 1):
            raise ConfigurationError("annotation.family_evidence.snapshot_date must be explicit")
    family_evidence_config = FamilyEvidenceConfig(
        enabled=family_enabled,
        joint_vcf=family_joint_vcf,
        joint_vcf_expected_sha256=family_joint_sha,
        pedigree=family_pedigree,
        pedigree_expected_sha256=family_pedigree_sha,
        source_version=family_source_version,
        snapshot_date=family_snapshot_date,
        tags=family_tags,
        loco=family_loco,
    )

    callability = _mapping(root.get("callability", {}), "callability")
    _only(callability, {"target_catalog"}, "callability")
    target_catalog = _mapping(callability.get("target_catalog", {}), "callability.target_catalog")
    _only(
        target_catalog,
        {"enabled", "path", "expected_sha256", "assay_scope"},
        "callability.target_catalog",
    )
    target_catalog_enabled = bool(target_catalog.get("enabled", False))
    target_catalog_path = _path(
        str(target_catalog.get("path", resource_cache / "targets" / "targets.tsv")), base
    )
    target_catalog_expected_sha256 = target_catalog.get("expected_sha256")
    if target_catalog_expected_sha256 is not None and (
        not isinstance(target_catalog_expected_sha256, str)
        or not re.fullmatch(r"[0-9a-fA-F]{64}", target_catalog_expected_sha256)
    ):
        raise ConfigurationError(
            "callability.target_catalog.expected_sha256 must be null or a SHA-256 hex digest"
        )
    if target_catalog_enabled and target_catalog_expected_sha256 is None:
        raise ConfigurationError("enabled callability.target_catalog requires expected_sha256")
    target_catalog_assay_scope = str(
        target_catalog.get("assay_scope", "small_variant_reference_confidence_vcf")
    ).strip()
    if not target_catalog_assay_scope:
        raise ConfigurationError("callability.target_catalog.assay_scope cannot be empty")

    output = _mapping(root.get("output", {}), "output")
    _only(
        output,
        {
            "agent_exploration_profile",
            "reference_confidence_extension",
            "partition_variants_by_chromosome",
            "overwrite",
        },
        "output",
    )
    privacy = _mapping(root.get("privacy", {}), "privacy")
    _only(privacy, {"include_input_paths", "include_sample_identifier"}, "privacy")

    if output_dir == work_dir:
        raise ConfigurationError("project.output_dir and project.work_dir must be different")

    return BuildConfig(
        config_path=config_path,
        output_dir=output_dir,
        work_dir=work_dir,
        bundle_name=bundle_name,
        input_vcf=input_vcf,
        sample=sample,
        genome_build=genome_build,
        reference_fasta=reference_fasta,
        reference_expected_sha256=expected_sha.lower() if expected_sha else None,
        contig_style=contig_style,
        resource_cache=resource_cache,
        allow_download=allow_download,
        genome_spec_commit=genome_spec_commit,
        bcftools=bcftools,
        check_ref=check_ref,
        quality=quality,
        vep_enabled=vep_enabled,
        vep_mode=vep_mode,
        vep_executable=vep_executable,
        vep_container_image=vep_container_image,
        vep_plugin_dir=vep_plugin_dir,
        vep_cache_dir=vep_cache_dir,
        vep_cache_version=vep_cache_version,
        vep_cache_manifest=vep_cache_manifest,
        vep_expected_version=vep_expected_version,
        vep_forks=vep_forks,
        vep_buffer_size=vep_buffer_size,
        clinvar_enabled=clinvar_enabled,
        clinvar_snapshot=clinvar_snapshot,
        clinvar_vcf=clinvar_vcf,
        clinvar_manifest=clinvar_manifest,
        acmg_sf_enabled=acmg_sf_enabled,
        acmg_sf_version=acmg_sf_version,
        acmg_sf_snapshot=acmg_sf_snapshot,
        acmg_sf_dataset=acmg_sf_dataset,
        acmg_sf_manifest=acmg_sf_manifest,
        acmg_sf_min_clinvar_review_stars=acmg_sf_min_review_stars,
        pgs_enabled=pgs_enabled,
        pgs_models=tuple(pgs_models),
        pgx=pgx_config,
        family_evidence=family_evidence_config,
        gwas=gwas_config,
        clingen=clingen_config,
        predictors=predictor_config,
        target_catalog_enabled=target_catalog_enabled,
        target_catalog_path=target_catalog_path,
        target_catalog_expected_sha256=(
            target_catalog_expected_sha256.lower()
            if target_catalog_expected_sha256 is not None
            else None
        ),
        target_catalog_assay_scope=target_catalog_assay_scope,
        agent_exploration_profile=bool(output.get("agent_exploration_profile", True)),
        reference_confidence_extension=bool(output.get("reference_confidence_extension", True)),
        partition_variants_by_chromosome=bool(output.get("partition_variants_by_chromosome", True)),
        overwrite=bool(output.get("overwrite", False)),
        include_input_paths=bool(privacy.get("include_input_paths", False)),
        include_sample_identifier=bool(privacy.get("include_sample_identifier", False)),
    )
