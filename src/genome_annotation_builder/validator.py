"""Strict core, provenance, absence, and Agent Profile bundle validation."""

from __future__ import annotations

import json
import math
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import duckdb
import jsonschema
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .arrow_schemas import (
    acmg_secondary_findings_evidence_schema,
    callability_schema,
    callable_region_schema,
    clingen_dosage_sensitivity_schema,
    clingen_gene_disease_schema,
    clinical_evidence_schema,
    clinical_finding_schema,
    clinvar_assertion_schema,
    family_evidence_schema,
    gene_index_schema,
    gwas_association_schema,
    gwas_catalog_evidence_schema,
    interpretation_schema,
    non_core_variant_schema,
    observed_partial_variant_schema,
    pgs_locus_audit_schema,
    pgs_model_audit_schema,
    pgx_gene_call_audit_schema,
    pgx_guidance_audit_schema,
    pgx_locus_audit_schema,
    pharmacogenomics_schema,
    prs_schema,
    reference_evidence_schema,
    target_callability_evidence_schema,
    variant_provenance_schema,
    variant_schema,
)
from .clingen_freeze import DOSAGE_SCORE_DESCRIPTIONS, GENE_DISEASE_CLASSIFICATIONS
from .clingen_resources import SEMANTIC_BOUNDARY as CLINGEN_SEMANTIC_BOUNDARY
from .constants import (
    CLINGEN_DOSAGE_ROWS,
    CLINGEN_GENE_DISEASE_ROWS,
    CLINGEN_GENE_LEVEL_CANONICAL_SHA256,
    CLINGEN_GENE_LEVEL_SNAPSHOT,
    CLINGEN_LICENSE,
    GENOME_SPEC_COMMIT,
    PHARMCAT_POSITIONS_ALLELES,
    PHARMCAT_POSITIONS_RECORDS,
    PHARMCAT_VERSION,
    PINNED_UPSTREAM_FILES,
    STRUCTURALLY_COMPLEX_PGX_GENES,
)
from .errors import BundleValidationError
from .pgx_materialize import _resolve_pgx_locus
from .predictor_doctor import (
    CONTROL_DEFINITION_SHA256,
    CONTROL_SET,
    CONTROL_VARIANTS,
)
from .predictor_doctor import SEMANTIC_BOUNDARY as PREDICTOR_DOCTOR_SEMANTIC_BOUNDARY
from .query import query_bundle
from .release_c.extension_stage import validate_release_c_extension_stage
from .release_c.portable_model import DOMAINS as RELEASE_C_DOMAINS
from .release_c.portable_model import PortableEvidenceError
from .release_c.portable_model import schema_for as release_c_schema_for
from .util import MultisetDigest, directory_digest, sha256_file

PROFILE_PATHS = (
    "field_glossary.json",
    "inventory.json",
    "raw_artifacts.json",
    "transformations.json",
    "omissions.json",
    "relationships.json",
    "readmygenome.md",
    "readmygenome/SKILL.md",
    "readmygenome/scripts/query.py",
    "readmygenome/reference/query_catalog.md",
    "readmygenome/reference/synonyms.json",
    "readmygenome/evals/canonical_questions.json",
)

_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_SHA1_RE = re.compile(r"[0-9a-f]{40}")

PREDICTOR_ORDER = ("loftee", "cadd", "revel", "alphamissense", "spliceai")
PREDICTOR_LABELS = {
    "loftee": "LOFTEE",
    "cadd": "CADD",
    "revel": "REVEL",
    "alphamissense": "AlphaMissense",
    "spliceai": "SpliceAI",
}
PREDICTOR_PLUGIN_PREFIXES = {
    "loftee": "LoF,",
    "cadd": "CADD,",
    "revel": "REVEL,",
    "alphamissense": "AlphaMissense,",
    "spliceai": "SpliceAI,",
}
PREDICTOR_PLUGIN_REQUIRED_ARGUMENTS = {
    "loftee": (
        "loftee_path:PREDICTORS/loftee_path",
        "human_ancestor_fa:PREDICTORS/loftee_human_ancestor",
        "gerp_bigwig:PREDICTORS/loftee_gerp_bigwig",
        "conservation_file:PREDICTORS/loftee_conservation_sqlite",
    ),
    "cadd": (
        "snv=PREDICTORS/cadd_snv",
        "indels=PREDICTORS/cadd_indel",
    ),
    "revel": ("file=PREDICTORS/revel_scores",),
    "alphamissense": (
        "file=PREDICTORS/alphamissense_scores",
        "transcript_match=1",
    ),
    "spliceai": (
        "snv=PREDICTORS/spliceai_snv",
        "indel=PREDICTORS/spliceai_indel_sentinel",
        "split_output=1",
    ),
}
PREDICTOR_ARTIFACT_KEYS = {
    "loftee": {
        "loftee_source_archive",
        "loftee_gerp_bigwig",
        "loftee_human_ancestor",
        "loftee_human_ancestor_fai",
        "loftee_human_ancestor_gzi",
        "loftee_sql_gz",
        "loftee_source_tree",
        "loftee_sqlite",
    },
    "cadd": {"cadd_snv", "cadd_snv_tbi", "cadd_indel", "cadd_indel_tbi"},
    "revel": {"revel_archive", "revel_scores", "revel_scores_tbi"},
    "alphamissense": {"alphamissense_scores", "alphamissense_scores_tbi"},
    "spliceai": {
        "spliceai_snv",
        "spliceai_snv_tbi",
        "spliceai_indel_sentinel",
        "spliceai_indel_sentinel_tbi",
    },
}


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and _SHA256_RE.fullmatch(value) is not None


def _validate_manifest_boundary(manifest: dict[str, Any], errors: list[str]) -> None:
    """Validate reproducibility identities not constrained by upstream schema."""

    if manifest.get("genome_build") != "GRCh38":
        errors.append("manifest genome_build is not GRCh38")
    if not _is_sha256(manifest.get("source_bundle_sha256")):
        errors.append("manifest source_bundle_sha256 is invalid")
    pipeline_version = manifest.get("pipeline_version")
    if not isinstance(pipeline_version, str) or not pipeline_version.startswith(
        "genome-annotation-builder-"
    ):
        errors.append("manifest pipeline_version is invalid")
    provenance = manifest.get("resource_provenance")
    if not isinstance(provenance, dict):
        errors.append("manifest resource_provenance is missing")
        return
    if not _is_sha256(provenance.get("builder_implementation_sha256")):
        errors.append("manifest builder implementation digest is invalid")
    if provenance.get("genome_spec_commit") != GENOME_SPEC_COMMIT:
        errors.append("manifest genome-spec commit is not the pinned revision")
    if provenance.get("schema_sha256") != PINNED_UPSTREAM_FILES["schema.json"]:
        errors.append("manifest normative schema digest is not pinned")
    if not _is_sha256(provenance.get("reference_fasta_sha256")):
        errors.append("manifest reference FASTA digest is invalid")
    if not _is_sha256(provenance.get("reference_fai_sha256")):
        errors.append("manifest reference FAI digest is invalid")


def _predictor_integrity_boundary(
    value: Any, plugin: str, errors: list[str]
) -> tuple[str | None, set[str]]:
    prefix = "manifest-sha256:"
    if not isinstance(value, str) or not value.startswith(prefix):
        errors.append(f"{plugin} predictor artifact lacks a manifest SHA-256 boundary")
        return None, set()
    manifest_text, separator, artifact_text = value[len(prefix) :].partition(";artifacts:")
    if not separator or not _is_sha256(manifest_text):
        errors.append(f"{plugin} predictor artifact has an invalid manifest SHA-256")
        return None, set()
    artifact_keys: set[str] = set()
    for token in artifact_text.split(",") if artifact_text else ():
        key, token_separator, digest = token.partition("=sha256:")
        if (
            not token_separator
            or not re.fullmatch(r"[a-z0-9_]+", key)
            or not _is_sha256(digest)
            or key in artifact_keys
        ):
            errors.append(f"{plugin} predictor artifact has a malformed digest entry")
            continue
        artifact_keys.add(key)
    if artifact_keys != PREDICTOR_ARTIFACT_KEYS[plugin]:
        errors.append(
            f"{plugin} predictor checksum coverage mismatch: "
            f"expected={sorted(PREDICTOR_ARTIFACT_KEYS[plugin])}, "
            f"observed={sorted(artifact_keys)}"
        )
    return manifest_text, artifact_keys


def _validate_predictor_metadata(
    documents: dict[str, dict[str, Any] | None], errors: list[str]
) -> dict[str, Any]:
    """Cross-check portable predictor claims without requiring private resources."""

    raw_document = documents.get("raw_artifacts.json")
    transformation_document = documents.get("transformations.json")
    omission_document = documents.get("omissions.json")
    if not all(
        isinstance(document, dict)
        for document in (raw_document, transformation_document, omission_document)
    ):
        return {"status": "metadata_unavailable", "enabled_plugins": []}

    assert isinstance(raw_document, dict)
    assert isinstance(transformation_document, dict)
    assert isinstance(omission_document, dict)
    stages = transformation_document.get("stages")
    if not isinstance(stages, list):
        errors.append("transformations.json stages must be a list")
        stages = []
    vep_stages = [
        stage
        for stage in stages
        if isinstance(stage, dict) and stage.get("stage_id") == "ensembl_vep_annotation"
    ]
    if len(vep_stages) > 1:
        errors.append("transformations.json contains duplicate Ensembl VEP stages")
    vep_stage = vep_stages[0] if vep_stages else None

    enabled_raw = vep_stage.get("enabled_predictors", []) if vep_stage else []
    enabled: list[str]
    if (
        not isinstance(enabled_raw, list)
        or any(not isinstance(value, str) for value in enabled_raw)
        or len(enabled_raw) != len(set(enabled_raw))
        or any(value not in PREDICTOR_ORDER for value in enabled_raw)
    ):
        errors.append("Ensembl VEP stage has an invalid enabled_predictors list")
        enabled = []
    else:
        enabled = list(enabled_raw)
        canonical_order = [value for value in PREDICTOR_ORDER if value in enabled]
        if enabled != canonical_order:
            errors.append("Ensembl VEP enabled_predictors are not in canonical order")

    raw_rows = raw_document.get("artifacts")
    if not isinstance(raw_rows, list):
        errors.append("raw_artifacts.json artifacts must be a list")
        raw_rows = []
    predictor_rows: dict[str, dict[str, Any]] = {}
    manifest_digests: set[str] = set()
    prefix = "offline VEP predictor resource: "
    for row in raw_rows:
        if not isinstance(row, dict):
            continue
        artifact_type = row.get("artifact_type")
        if not isinstance(artifact_type, str) or not artifact_type.startswith(prefix):
            continue
        plugin = artifact_type[len(prefix) :]
        if plugin not in PREDICTOR_ORDER:
            errors.append(f"raw_artifacts.json names an unknown VEP predictor: {plugin!r}")
            continue
        if plugin in predictor_rows:
            errors.append(f"raw_artifacts.json duplicates the {plugin} predictor boundary")
            continue
        predictor_rows[plugin] = row
        if row.get("disposition") != "referenced":
            errors.append(f"{plugin} predictor resource must be referenced, not embedded")
        if row.get("redistribution_allowed") is not False:
            errors.append(f"{plugin} predictor resource must not assert redistribution rights")
        if not isinstance(row.get("license"), str) or not row["license"].strip():
            errors.append(f"{plugin} predictor resource lacks a licence/terms statement")
        if not isinstance(row.get("coverage"), dict):
            errors.append(f"{plugin} predictor resource lacks a typed coverage boundary")
        else:
            coverage = row["coverage"]
            if coverage.get("assembly") != "GRCh38":
                errors.append(f"{plugin} predictor coverage is not GRCh38")
            if plugin == "spliceai" and (
                coverage.get("variant_types") != ["SNV"]
                or coverage.get("indel_scores_available") is not False
                or coverage.get("indel_runtime_input") != "explicit empty sentinel"
            ):
                errors.append("SpliceAI coverage does not preserve its SNV-only boundary")
        manifest_digest, _keys = _predictor_integrity_boundary(
            row.get("integrity_identifier"), plugin, errors
        )
        if manifest_digest is not None:
            manifest_digests.add(manifest_digest)

    raw_plugins = set(predictor_rows)
    if raw_plugins != set(enabled):
        errors.append(
            "enabled VEP predictors and raw resource provenance disagree: "
            f"enabled={sorted(enabled)}, resources={sorted(raw_plugins)}"
        )
    if len(manifest_digests) > 1:
        errors.append("enabled predictors do not share one immutable resource manifest")

    if vep_stage is None and predictor_rows:
        errors.append("predictor resources are declared without an Ensembl VEP transformation")
    if vep_stage is not None:
        stage_coverage = vep_stage.get("predictor_coverage")
        expected_coverage = {
            plugin: predictor_rows[plugin].get("coverage")
            for plugin in enabled
            if plugin in predictor_rows
        }
        if stage_coverage != expected_coverage:
            errors.append("VEP predictor coverage differs from raw resource provenance")
        inputs = vep_stage.get("inputs")
        expected_input = (
            "sha256:" + next(iter(manifest_digests)) if len(manifest_digests) == 1 else None
        )
        if enabled and (
            not isinstance(inputs, list)
            or expected_input is None
            or inputs.count(expected_input) != 1
        ):
            errors.append("VEP stage is not bound to exactly one predictor resource manifest")

        command = vep_stage.get("command_semantics")
        if not isinstance(command, list) or any(
            not isinstance(argument, str) for argument in command
        ):
            if enabled:
                errors.append("VEP predictor command semantics are absent or invalid")
            command = []
        specs: list[str] = []
        for index, argument in enumerate(command):
            if argument != "--plugin":
                continue
            if index + 1 >= len(command):
                errors.append("VEP command contains a dangling --plugin argument")
                continue
            specs.append(command[index + 1])
        matched_specs: dict[str, list[str]] = {plugin: [] for plugin in PREDICTOR_ORDER}
        for spec in specs:
            matches = [
                plugin
                for plugin, plugin_prefix in PREDICTOR_PLUGIN_PREFIXES.items()
                if spec.startswith(plugin_prefix)
            ]
            if len(matches) != 1:
                errors.append(f"VEP command names an unrecognized predictor plugin: {spec}")
                continue
            matched_specs[matches[0]].append(spec)
        for plugin in PREDICTOR_ORDER:
            plugin_specs = matched_specs[plugin]
            expected_count = 1 if plugin in enabled else 0
            if len(plugin_specs) != expected_count:
                errors.append(
                    f"VEP command has {len(plugin_specs)} {plugin} plugin arguments; "
                    f"expected {expected_count}"
                )
                continue
            if plugin_specs and any(
                token not in plugin_specs[0]
                for token in PREDICTOR_PLUGIN_REQUIRED_ARGUMENTS[plugin]
            ):
                errors.append(f"VEP {plugin} plugin arguments do not preserve resource semantics")
        if enabled and "--dir_plugins" not in command:
            errors.append("VEP predictor command lacks an explicit plugin directory")

    omission_rows = omission_document.get("omissions")
    if not isinstance(omission_rows, list):
        errors.append("omissions.json omissions must be a list")
        omission_rows = []
    omission_scopes = [row.get("scope") for row in omission_rows if isinstance(row, dict)]
    for plugin, label in PREDICTOR_LABELS.items():
        scope = f"{label} predictor annotation"
        expected_count = 0 if plugin in enabled else 1
        actual_count = omission_scopes.count(scope)
        if actual_count != expected_count:
            errors.append(
                f"{label} disabled-predictor omission count is {actual_count}; "
                f"expected {expected_count}"
            )
    spliceai_indel_count = omission_scopes.count("SpliceAI indel predictions")
    expected_spliceai_indel_count = 1 if "spliceai" in enabled else 0
    if spliceai_indel_count != expected_spliceai_indel_count:
        errors.append(
            "SpliceAI indel unavailability omission count is "
            f"{spliceai_indel_count}; expected {expected_spliceai_indel_count}"
        )

    return {
        "status": "validated",
        "enabled_plugins": enabled,
        "resource_artifacts": len(predictor_rows),
        "resource_manifest_count": len(manifest_digests),
    }


def _validate_predictor_doctor_metadata(
    documents: dict[str, dict[str, Any] | None],
    manifest: dict[str, Any],
    enabled_plugins: list[str],
    errors: list[str],
) -> dict[str, Any]:
    """Cross-check the external fixed-control report across portable metadata."""

    raw_document = documents.get("raw_artifacts.json")
    transformation_document = documents.get("transformations.json")
    raw_rows = raw_document.get("artifacts", []) if isinstance(raw_document, dict) else []
    stages = (
        transformation_document.get("stages", [])
        if isinstance(transformation_document, dict)
        else []
    )
    raw_matches = [
        row
        for row in raw_rows
        if isinstance(row, dict)
        and row.get("artifact_type") == "offline VEP predictor semantic control report"
    ]
    stage_matches = [
        row
        for row in stages
        if isinstance(row, dict) and row.get("stage_id") == "offline_predictor_semantic_controls"
    ]
    provenance = manifest.get("resource_provenance", {}).get("offline_predictors")
    if not enabled_plugins:
        if raw_matches or stage_matches or provenance is not None:
            errors.append("disabled predictor adapters retain semantic-control metadata")
        return {"status": "not_present", "control_count": 0}
    if len(raw_matches) != 1:
        errors.append("enabled predictors require exactly one semantic-control raw artifact")
    if len(stage_matches) != 1:
        errors.append("enabled predictors require exactly one semantic-control stage")
    if not isinstance(provenance, dict):
        errors.append("enabled predictors lack manifest semantic-control provenance")
        return {"status": "invalid_provenance", "control_count": 0}
    controls = provenance.get("semantic_controls")
    if not isinstance(controls, dict):
        errors.append("predictor provenance lacks semantic_controls")
        controls = {}
    expected_plugins = list(enabled_plugins)
    expected_common = {
        "control_set": CONTROL_SET,
        "control_definition_sha256": CONTROL_DEFINITION_SHA256,
        "record_count": len(CONTROL_VARIANTS),
        "semantic_boundary": PREDICTOR_DOCTOR_SEMANTIC_BOUNDARY,
    }
    if provenance.get("enabled_plugins") != expected_plugins:
        errors.append("predictor semantic-control provenance has a plugin-set mismatch")
    if not _is_sha256(provenance.get("manifest_sha256")):
        errors.append("predictor semantic-control provenance has an invalid manifest digest")
    for key, expected in expected_common.items():
        if controls.get(key) != expected:
            errors.append(f"predictor semantic-control provenance mismatch: {key}")
    for key in (
        "doctor_run_key",
        "report_sha256",
        "record_identity_sha256",
        "reference_fasta_sha256",
        "vep_cache_manifest_sha256",
    ):
        if not _is_sha256(controls.get(key)):
            errors.append(f"predictor semantic-control provenance has invalid {key}")
    if controls.get("reference_fasta_sha256") != manifest.get("resource_provenance", {}).get(
        "reference_fasta_sha256"
    ):
        errors.append("predictor semantic controls used a different reference FASTA")

    result_rows = controls.get("control_results")
    if not isinstance(result_rows, list):
        errors.append("predictor semantic-control results must be a list")
        result_rows = []
    expected_ids = [control["control_id"] for control in CONTROL_VARIANTS]
    if [row.get("control_id") for row in result_rows if isinstance(row, dict)] != expected_ids:
        errors.append("predictor semantic-control result set/order mismatch")
    observed_plugin_controls: set[str] = set()
    for control, row in zip(CONTROL_VARIANTS, result_rows, strict=False):
        if not isinstance(row, dict):
            errors.append("predictor semantic-control result is not an object")
            continue
        checks = row.get("checks")
        if not isinstance(checks, dict):
            errors.append("predictor semantic-control result lacks checks")
            continue
        expected_checks = {
            plugin: value
            for plugin, value in control["checks"].items()
            if plugin in enabled_plugins
        }
        if set(checks) != set(expected_checks):
            errors.append(f"predictor semantic-control check set mismatch: {control['control_id']}")
            continue
        for plugin, expected in expected_checks.items():
            details = checks.get(plugin)
            if (
                not isinstance(details, dict)
                or details.get("status") != "PASS"
                or details.get("expected") != expected
                or details.get("observed") != expected
            ):
                errors.append(
                    f"predictor semantic-control value mismatch: {control['control_id']}/{plugin}"
                )
            observed_plugin_controls.add(plugin)
    if observed_plugin_controls != set(enabled_plugins):
        errors.append("predictor semantic controls do not cover every enabled plugin")

    if raw_matches:
        raw = raw_matches[0]
        expected_raw = {
            "integrity_identifier": "sha256:" + str(controls.get("report_sha256")),
            "control_set": CONTROL_SET,
            "control_definition_sha256": CONTROL_DEFINITION_SHA256,
            "doctor_run_key": controls.get("doctor_run_key"),
            "enabled_plugins": expected_plugins,
            "record_identity_sha256": controls.get("record_identity_sha256"),
            "reference_fasta_sha256": controls.get("reference_fasta_sha256"),
            "vep_cache_manifest_sha256": controls.get("vep_cache_manifest_sha256"),
        }
        for key, expected in expected_raw.items():
            if raw.get(key) != expected:
                errors.append(f"predictor semantic-control raw artifact mismatch: {key}")
    if stage_matches:
        stage = stage_matches[0]
        expected_stage = {
            **expected_common,
            "doctor_run_key": controls.get("doctor_run_key"),
            "report_sha256": controls.get("report_sha256"),
            "record_identity_sha256": controls.get("record_identity_sha256"),
            "reference_fasta_sha256": controls.get("reference_fasta_sha256"),
            "vep_cache_manifest_sha256": controls.get("vep_cache_manifest_sha256"),
            "enabled_plugins": expected_plugins,
            "control_results": result_rows,
            "status": "completed",
        }
        for key, expected in expected_stage.items():
            if stage.get(key) != expected:
                errors.append(f"predictor semantic-control stage mismatch: {key}")
    return {
        "status": "validated",
        "control_count": len(result_rows),
        "enabled_plugins": expected_plugins,
        "covered_plugins": sorted(observed_plugin_controls),
    }


def _validate_predictor_values(
    connection: duckdb.DuckDBPyConnection,
    enabled_plugins: set[str],
    errors: list[str],
) -> dict[str, int]:
    result = connection.execute(
        """
        SELECT
          count(*) FILTER (WHERE predictions.cadd_phred IS NOT NULL) AS cadd_values,
          count(*) FILTER (WHERE predictions.revel IS NOT NULL) AS revel_values,
          count(*) FILTER (WHERE predictions.alphamissense IS NOT NULL) AS alphamissense_values,
          count(*) FILTER (WHERE predictions.spliceai_max IS NOT NULL) AS spliceai_values,
          count(*) FILTER (
            WHERE coalesce(list_count(list_filter(
              consequence.all_annotations, item -> item.lof IS NOT NULL
            )), 0) > 0
          ) AS loftee_values,
          count(*) FILTER (
            WHERE predictions.cadd_phred IS NOT NULL AND
              (NOT isfinite(predictions.cadd_phred) OR predictions.cadd_phred < 0)
          ) AS invalid_cadd_values,
          count(*) FILTER (
            WHERE predictions.revel IS NOT NULL AND
              (NOT isfinite(predictions.revel) OR predictions.revel < 0 OR predictions.revel > 1)
          ) AS invalid_revel_values,
          count(*) FILTER (
            WHERE predictions.alphamissense IS NOT NULL AND
              (NOT isfinite(predictions.alphamissense) OR
               predictions.alphamissense < 0 OR predictions.alphamissense > 1)
          ) AS invalid_alphamissense_values,
          count(*) FILTER (
            WHERE predictions.spliceai_max IS NOT NULL AND
              (NOT isfinite(predictions.spliceai_max) OR
               predictions.spliceai_max < 0 OR predictions.spliceai_max > 1)
          ) AS invalid_spliceai_values,
          count(*) FILTER (
            WHERE coalesce(list_count(list_filter(
              consequence.all_annotations,
              item -> item.lof IS NOT NULL AND item.lof NOT IN ('HC', 'LC')
            )), 0) > 0
          ) AS invalid_loftee_values
        FROM variants
        """
    ).fetchone()
    names = [description[0] for description in connection.description]
    metrics = {key: int(value) for key, value in zip(names, result, strict=True)}
    for plugin in PREDICTOR_ORDER:
        value_count = metrics[f"{plugin}_values"]
        if plugin not in enabled_plugins and value_count:
            errors.append(
                f"variants.parquet has {value_count} {plugin} value(s) without an enabled adapter"
            )
        invalid_count = metrics[f"invalid_{plugin}_values"]
        if invalid_count:
            errors.append(
                f"variants.parquet has {invalid_count} invalid {plugin} predictor value(s)"
            )
    return metrics


def _definition_validator(schema: dict[str, Any], definition: str) -> jsonschema.Validator:
    wrapper = {
        "$schema": schema.get("$schema", "https://json-schema.org/draft/2020-12/schema"),
        "$ref": f"#/$defs/{definition}",
        "$defs": schema["$defs"],
    }
    return jsonschema.Draft202012Validator(wrapper, format_checker=jsonschema.FormatChecker())


def _validate_json_document(
    root: Path,
    relative: str,
    schema: dict[str, Any],
    definition: str,
    errors: list[str],
) -> dict[str, Any] | None:
    path = root / relative
    if not path.is_file():
        errors.append(f"missing JSON document: {relative}")
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"invalid JSON in {relative}: {exc}")
        return None
    for issue in _definition_validator(schema, definition).iter_errors(value):
        location = ".".join(str(part) for part in issue.absolute_path) or "<root>"
        errors.append(f"{relative}:{location}: {issue.message}")
    return value


def _directory_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _parquet_files(path: Path) -> list[Path]:
    return [path] if path.is_file() else sorted(path.rglob("*.parquet"))


def _parquet_rows(path: Path) -> int:
    return sum(pq.ParquetFile(item).metadata.num_rows for item in _parquet_files(path))


def _validate_parquet_physical(
    path: Path,
    expected_schema: pa.Schema,
    label: str,
    errors: list[str],
) -> int:
    files = _parquet_files(path)
    if not files:
        errors.append(f"{label} has no Parquet data files")
        return 0
    rows = 0
    for file_path in files:
        parquet = pq.ParquetFile(file_path)
        rows += parquet.metadata.num_rows
        actual_schema = parquet.schema_arrow
        if not actual_schema.equals(expected_schema, check_metadata=False):
            errors.append(f"{label} Arrow schema mismatch in {file_path.relative_to(path.parent)}")
        if parquet.metadata.format_version != "2.6":
            errors.append(f"{label} requires Parquet 2.6, got {parquet.metadata.format_version}")
        for row_group_index in range(parquet.metadata.num_row_groups):
            row_group = parquet.metadata.row_group(row_group_index)
            for column_index in range(row_group.num_columns):
                column = row_group.column(column_index)
                if column.compression != "ZSTD":
                    errors.append(f"{label} column {column.path_in_schema} is not ZSTD-compressed")
    return rows


def _manifest_coverage(root: Path, manifest: dict[str, Any], errors: list[str]) -> None:
    entries: list[Path] = []
    for relative in manifest.get("files", {}):
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts or relative == "manifest.json":
            errors.append(f"unsafe or self-referential manifest path: {relative}")
            continue
        entries.append(path)
    for index, left in enumerate(entries):
        for right in entries[index + 1 :]:
            if left in right.parents or right in left.parents:
                errors.append(f"overlapping manifest paths: {left} and {right}")
    for file_path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = file_path.relative_to(root)
        if relative == Path("manifest.json"):
            continue
        matches = [entry for entry in entries if relative == entry or entry in relative.parents]
        if len(matches) != 1:
            errors.append(
                f"physical file must be covered by exactly one manifest entry: {relative}"
            )


def _validate_manifest_files(
    root: Path, manifest: dict[str, Any], errors: list[str]
) -> dict[str, int]:
    rows: dict[str, int] = {}
    for relative, metadata in manifest.get("files", {}).items():
        target = root / relative
        if not target.exists():
            errors.append(f"manifest target is missing: {relative}")
            continue
        actual_digest = directory_digest(target) if target.is_dir() else sha256_file(target)
        if actual_digest != metadata.get("sha256"):
            errors.append(f"manifest SHA-256 mismatch: {relative}")
        actual_bytes = _directory_bytes(target) if target.is_dir() else target.stat().st_size
        if metadata.get("bytes") is not None and metadata["bytes"] != actual_bytes:
            errors.append(f"manifest byte count mismatch: {relative}")
        declared_rows = metadata.get("rows")
        if declared_rows is not None:
            try:
                actual_rows = _parquet_rows(target)
            except Exception as exc:
                errors.append(f"cannot verify declared rows for {relative}: {exc}")
            else:
                rows[relative] = actual_rows
                if actual_rows != declared_rows:
                    errors.append(f"manifest row count mismatch: {relative}")
    return rows


def _field_paths(table: str, schema: pa.Schema) -> set[str]:
    output: set[str] = set()

    def visit(field: pa.Field, prefix: str) -> None:
        output.add(f"{table}.{prefix}")
        if pa.types.is_struct(field.type):
            for child in field.type:
                visit(child, f"{prefix}.{child.name}")

    for field in schema:
        visit(field, field.name)
    return output


def _validate_variant_rows(
    connection: duckdb.DuckDBPyConnection,
    manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, int]:
    standard_chroms = ",".join(f"'chr{value}'" for value in [*range(1, 23), "X", "Y", "M"])
    result = connection.execute(
        f"""
        SELECT
          count(*) AS total,
          count(DISTINCT variant_id) AS distinct_ids,
          count(*) FILTER (WHERE
            variant_id IS NULL OR chrom IS NULL OR pos IS NULL OR ref IS NULL OR alt IS NULL OR
            quality IS NULL OR genotype IS NULL OR gene IS NULL OR consequence IS NULL OR
            population_frequency IS NULL OR pathogenicity IS NULL OR predictions IS NULL OR
            pharmacogenomics IS NULL OR trait_associations IS NULL OR clinical_grade IS NULL
          ) AS missing_required,
          count(*) FILTER (WHERE chrom NOT IN ({standard_chroms})) AS bad_chrom,
          count(*) FILTER (WHERE pos < 1 OR NOT regexp_full_match(ref, '^[ACGT]+$') OR
                                  NOT regexp_full_match(alt, '^[ACGT]+$')) AS bad_allele,
          count(*) FILTER (WHERE variant_id != chrom || ':' || CAST(pos AS VARCHAR) || ':' || ref || ':' || alt)
            AS bad_variant_id,
          count(*) FILTER (WHERE variant_type NOT IN ('SNV','INS','DEL','MNV','SV')) AS bad_type,
          count(*) FILTER (WHERE quality.call_confidence NOT IN ('high','medium','low') OR
                                  quality.origin NOT IN ('observed','imputed')) AS bad_quality_enum,
          count(*) FILTER (WHERE quality.dp < 0 OR quality.gq < 0 OR
                                  quality.imputation_r2 < 0 OR quality.imputation_r2 > 1) AS
            bad_quality_value,
          count(*) FILTER (WHERE list_count(genotype.gt) NOT BETWEEN 1 AND 2 OR
                                  list_min(genotype.gt) < 0 OR list_max(genotype.gt) > 1 OR
                                  NOT list_contains(genotype.gt, 1)) AS bad_gt,
          count(*) FILTER (WHERE genotype.zygosity != CASE
              WHEN list_min(genotype.gt)=1 AND list_max(genotype.gt)=1 THEN 'hom_alt'
              WHEN list_contains(genotype.gt,0) AND list_contains(genotype.gt,1) THEN 'het'
              ELSE 'hom_ref' END) AS bad_zygosity,
          count(*) FILTER (WHERE genotype.allele_balance IS NOT NULL AND
                                  (genotype.allele_balance < 0 OR genotype.allele_balance > 1))
            AS bad_allele_balance,
          count(*) FILTER (WHERE population_frequency.is_rare IS NULL OR
                                  pathogenicity.clinvar_has_conflicts IS NULL OR
                                  pathogenicity.is_actionable IS NULL OR
                                  predictions.consensus_damaging IS NULL OR
                                  pharmacogenomics.is_pgx IS NULL OR
                                  trait_associations.is_gwas_hit IS NULL) AS missing_boolean,
          count(*) FILTER (WHERE population_frequency.gnomad_af IS NOT NULL AND
                                  (population_frequency.gnomad_af < 0 OR population_frequency.gnomad_af > 1) OR
                                  population_frequency.gnomad_af_popmax IS NOT NULL AND
                                  (population_frequency.gnomad_af_popmax < 0 OR population_frequency.gnomad_af_popmax > 1))
            AS bad_frequency,
          count(*) FILTER (WHERE pathogenicity.clinvar_review_stars IS NOT NULL AND
                                  (pathogenicity.clinvar_review_stars < 0 OR
                                   pathogenicity.clinvar_review_stars > 4)) AS bad_review_stars,
          count(*) FILTER (WHERE trait_associations.effect_size IS NOT NULL AND
                                  (trait_associations.effect_type IS NULL OR
                                   trait_associations.effect_allele IS NULL)) AS bad_effect_reference
        FROM variants
        """
    ).fetchone()
    names = [description[0] for description in connection.description]
    metrics = dict(zip(names, result, strict=True))
    for name, value in metrics.items():
        if name not in {"total", "distinct_ids"} and value:
            errors.append(f"variants.parquet invariant {name} failed for {value} row(s)")
    if metrics["total"] != metrics["distinct_ids"]:
        errors.append("variants.parquet contains duplicate variant_id values")

    definitions = manifest.get("consensus_definitions", {})
    rare = definitions.get("is_rare", {})
    threshold = rare.get("threshold")
    null_treatment = rare.get("null_treatment")
    if isinstance(threshold, (int, float)) and null_treatment in {"rare", "not_rare"}:
        null_value = "true" if null_treatment == "rare" else "false"
        mismatch = connection.execute(
            "SELECT count(*) FROM variants WHERE population_frequency.is_rare != "
            f"CASE WHEN population_frequency.gnomad_af_popmax IS NULL THEN {null_value} "
            f"ELSE population_frequency.gnomad_af_popmax < {float(threshold)} END"
        ).fetchone()[0]
        if mismatch:
            errors.append(f"is_rare manifest rule mismatch for {mismatch} row(s)")

    damaging = definitions.get("consensus_damaging", {})
    tools = {tool.get("name"): tool.get("threshold") for tool in damaging.get("tools", [])}
    if all(
        isinstance(tools.get(name), (int, float)) for name in ("REVEL", "AlphaMissense", "CADD")
    ):
        mismatch = connection.execute(
            """
            SELECT count(*) FROM variants WHERE predictions.consensus_damaging !=
              (consequence.most_severe = 'missense_variant' AND
               ((CASE WHEN predictions.revel > ? THEN 1 ELSE 0 END) +
                (CASE WHEN predictions.alphamissense > ? THEN 1 ELSE 0 END) +
                (CASE WHEN predictions.cadd_phred > ? THEN 1 ELSE 0 END)) >= ?)
            """,
            [tools["REVEL"], tools["AlphaMissense"], tools["CADD"], damaging.get("n", 2)],
        ).fetchone()[0]
        if mismatch:
            errors.append(f"consensus_damaging manifest rule mismatch for {mismatch} row(s)")
    actionable = definitions.get("is_actionable", {})
    if actionable.get("gene_list", "").startswith("EMPTY_"):
        mismatch = connection.execute(
            "SELECT count(*) FROM variants WHERE pathogenicity.is_actionable"
        ).fetchone()[0]
        if mismatch:
            errors.append("Release-A empty actionable rule has true rows")
    clinical = definitions.get("clinical_grade", {})
    if clinical.get("qualifying_sources") == []:
        mismatch = connection.execute(
            "SELECT count(*) FROM variants WHERE clinical_grade"
        ).fetchone()[0]
        if mismatch:
            errors.append("Release-A empty clinical-grade rule has true rows")
    return {key: int(value) for key, value in metrics.items()}


def _validate_reference_and_callability(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    errors: list[str],
) -> dict[str, int]:
    metrics: dict[str, int] = {}
    if "reference_evidence" in views:
        result = connection.execute(
            """
            SELECT count(*) AS total,
              count(DISTINCT evidence_id) AS distinct_ids,
              count(*) FILTER (WHERE start_pos < 1 OR end_pos < start_pos) AS bad_interval,
              count(*) FILTER (WHERE state NOT IN
                ('reference_block','explicit_reference_site','no_call')) AS bad_state,
              count(*) FILTER (WHERE state IN ('reference_block','explicit_reference_site') AND
                (reference_observed IS DISTINCT FROM true OR gt IS NULL OR
                 list_min(gt) != 0 OR list_max(gt) != 0)) AS bad_reference_semantics,
              count(*) FILTER (WHERE state = 'no_call' AND
                (reference_observed IS NOT NULL OR callable IS DISTINCT FROM false)) AS bad_nocall_semantics,
              count(*) FILTER (WHERE NOT regexp_full_match(source_record_digest,'^[0-9a-f]{64}$') OR
                source_ordinal < 1) AS bad_provenance
            FROM reference_evidence
            """
        ).fetchone()
        names = [description[0] for description in connection.description]
        values = dict(zip(names, result, strict=True))
        metrics.update({f"reference_{key}": int(value) for key, value in values.items()})
        for name, value in values.items():
            if name not in {"total", "distinct_ids"} and value:
                errors.append(f"reference-confidence invariant {name} failed for {value} row(s)")
        if values["total"] != values["distinct_ids"]:
            errors.append("reference-confidence evidence_id values are not unique")
    if "callable_regions" in views:
        result = connection.execute(
            """
            SELECT count(*) AS total, count(DISTINCT region_id) AS distinct_ids,
              count(*) FILTER (WHERE start_pos < 1 OR end_pos < start_pos OR callable IS NULL OR
                                     assay_scope IS NULL OR derivation IS NULL) AS bad_row,
              count(*) FILTER (WHERE minimum_depth IS NOT NULL AND minimum_depth < 0) AS bad_depth
            FROM callable_regions
            """
        ).fetchone()
        names = [description[0] for description in connection.description]
        values = dict(zip(names, result, strict=True))
        metrics.update({f"callable_{key}": int(value) for key, value in values.items()})
        for name, value in values.items():
            if name not in {"total", "distinct_ids"} and value:
                errors.append(f"callable_regions invariant {name} failed for {value} row(s)")
        if values["total"] != values["distinct_ids"]:
            errors.append("callable_regions region_id values are not unique")
    return metrics


def _validate_target_callability(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, int]:
    has_core = "callability" in views
    has_extension = "target_callability_evidence" in views
    if has_core != has_extension:
        errors.append(
            "callability.parquet and target_callability_evidence extension must appear together"
        )
        return {}
    if not has_core:
        return {}

    result = connection.execute(
        """
        SELECT count(*) AS total,
               count(DISTINCT callability_id) AS distinct_ids,
               count(*) FILTER (WHERE pos < 1 OR chrom NOT IN (
                 'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                 'chr20','chr21','chr22','chrX','chrY','chrM')) AS bad_locus,
               count(*) FILTER (WHERE call_confidence IS NOT NULL AND
                                      call_confidence NOT IN ('high','medium','low')) AS bad_confidence,
               count(*) FILTER (WHERE evidence_scope NOT IN
                                      ('site','target_panel','derived_interval')) AS bad_scope,
               count(*) FILTER (WHERE depth < 0 OR genotype_quality < 0) AS bad_quality,
               count(*) FILTER (WHERE variant_id IS NOT NULL AND
                 NOT regexp_full_match(variant_id,
                   '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M):[1-9][0-9]*:[ACGT]+:[ACGT]+$'))
                 AS bad_variant_id
        FROM callability
        """
    ).fetchone()
    names = [description[0] for description in connection.description]
    metrics = {f"callability_{key}": int(value) for key, value in zip(names, result, strict=True)}
    if result[0] != result[1]:
        errors.append("callability.parquet callability_id values are not unique")
    for key, value in zip(names[2:], result[2:], strict=True):
        if value:
            errors.append(f"callability.parquet invariant {key} failed for {value} row(s)")

    extension_result = connection.execute(
        """
        SELECT count(*) AS total,
               count(DISTINCT callability_id) AS distinct_ids,
               count(DISTINCT target_id) AS distinct_targets,
               count(*) FILTER (WHERE end_pos < pos OR catalog_row_number < 2 OR
                 NOT regexp_full_match(catalog_sha256,'^[0-9a-f]{64}$') OR
                 NOT regexp_full_match(catalog_row_digest,'^[0-9a-f]{64}$')) AS bad_catalog,
               count(*) FILTER (WHERE resolution NOT IN (
                 'non_reference_observed','overlapping_variant_observed','no_call',
                 'reference_observed','low_confidence_reference',
                 'callable_without_reference','uncallable_without_reference','unknown'))
                 AS bad_resolution,
               count(*) FILTER (WHERE observed_ploidy IS NOT NULL AND
                                      observed_ploidy NOT IN (1,2)) AS bad_ploidy,
               count(*) FILTER (WHERE callable_min_dp < 0 OR callable_min_gq < 0) AS bad_threshold,
               count(*) FILTER (WHERE
                 (ref IS NULL) != (alt IS NULL) OR
                 (ref IS NULL AND (variant_id IS NOT NULL OR end_pos != pos)) OR
                 (ref IS NOT NULL AND (
                   end_pos != pos + length(ref) - 1 OR
                   variant_id != chrom || ':' || CAST(pos AS VARCHAR) || ':' || ref || ':' || alt
                 ))) AS bad_target_identity,
               count(*) FILTER (WHERE
                 resolution IN ('non_reference_observed','overlapping_variant_observed') AND
                   (reference_observed IS DISTINCT FROM false OR
                    list_count(observed_variant_ids) +
                    list_count(observed_non_core_evidence_ids) = 0) OR
                 resolution = 'no_call' AND
                   (reference_observed IS NOT NULL OR callable IS DISTINCT FROM false OR
                    list_count(overlapping_no_call_evidence_ids) = 0) OR
                 resolution = 'reference_observed' AND
                   (reference_observed IS DISTINCT FROM true OR callable IS DISTINCT FROM true OR
                    list_count(selected_reference_evidence_ids) = 0) OR
                 resolution = 'low_confidence_reference' AND
                   (reference_observed IS DISTINCT FROM true OR callable IS DISTINCT FROM false OR
                    list_count(selected_reference_evidence_ids) = 0) OR
                 resolution = 'callable_without_reference' AND
                   (reference_observed IS NOT NULL OR callable IS DISTINCT FROM true OR
                    list_count(callable_region_ids) = 0) OR
                 resolution = 'uncallable_without_reference' AND
                   (reference_observed IS NOT NULL OR callable IS DISTINCT FROM false OR
                    list_count(callable_region_ids) = 0) OR
                 resolution = 'unknown' AND
                   (reference_observed IS NOT NULL OR callable IS NOT NULL OR
                    list_count(observed_variant_ids) +
                    list_count(observed_non_core_evidence_ids) +
                    list_count(overlapping_no_call_evidence_ids) +
                    list_count(selected_reference_evidence_ids) +
                    list_count(callable_region_ids) != 0)
               ) AS bad_semantics,
               count(*) FILTER (WHERE list_count(observed_variant_ids) !=
                                      list_count(observed_variant_source_digests) OR
                 list_has_any(observed_variant_source_digests,
                   list_filter(observed_variant_source_digests,
                     value -> NOT regexp_full_match(value,'^[0-9a-f]{64}$')))) AS bad_variant_proof
        FROM target_callability_evidence
        """
    ).fetchone()
    extension_names = [description[0] for description in connection.description]
    metrics.update(
        {
            f"target_evidence_{key}": int(value)
            for key, value in zip(extension_names, extension_result, strict=True)
        }
    )
    if extension_result[0] != extension_result[1] or extension_result[0] != extension_result[2]:
        errors.append("target callability IDs and target IDs must both be unique")
    for key, value in zip(extension_names[3:], extension_result[3:], strict=True):
        if value:
            errors.append(f"target_callability_evidence invariant {key} failed for {value} row(s)")

    mismatch = connection.execute(
        """
        SELECT count(*)
        FROM callability c
        FULL OUTER JOIN target_callability_evidence e USING (callability_id)
        WHERE c.callability_id IS NULL OR e.callability_id IS NULL OR
              c.variant_id IS DISTINCT FROM e.variant_id OR
              c.gene_symbol IS DISTINCT FROM e.gene_symbol OR
              c.chrom IS DISTINCT FROM e.chrom OR c.pos IS DISTINCT FROM e.pos OR
              c.reference_observed IS DISTINCT FROM e.reference_observed OR
              c.callable IS DISTINCT FROM e.callable OR
              c.call_confidence IS DISTINCT FROM e.call_confidence OR
              c.depth IS DISTINCT FROM e.depth OR
              c.genotype_quality IS DISTINCT FROM e.genotype_quality OR
              c.assay_scope IS DISTINCT FROM e.assay_scope OR
              c.evidence_scope IS DISTINCT FROM e.evidence_scope
        """
    ).fetchone()[0]
    metrics["target_core_extension_mismatches"] = int(mismatch)
    if mismatch:
        errors.append(f"callability core/extension mismatch for {mismatch} row(s)")

    orphan_variants = connection.execute(
        """
        SELECT count(*) FROM target_callability_evidence e,
          UNNEST(e.observed_variant_ids) AS item(variant_id)
        LEFT JOIN variants v USING (variant_id)
        WHERE v.variant_id IS NULL
        """
    ).fetchone()[0]
    orphan_reference = connection.execute(
        """
        SELECT count(*) FROM (
          SELECT unnest(selected_reference_evidence_ids) AS evidence_id
          FROM target_callability_evidence
          UNION ALL
          SELECT unnest(overlapping_no_call_evidence_ids) AS evidence_id
          FROM target_callability_evidence
        ) e LEFT JOIN reference_evidence r USING (evidence_id)
        WHERE r.evidence_id IS NULL
        """
    ).fetchone()[0]
    orphan_regions = connection.execute(
        """
        SELECT count(*) FROM target_callability_evidence e,
          UNNEST(e.callable_region_ids) AS item(region_id)
        LEFT JOIN callable_regions r USING (region_id)
        WHERE r.region_id IS NULL
        """
    ).fetchone()[0]
    orphan_non_core = 0
    non_core_links = connection.execute(
        "SELECT coalesce(sum(list_count(observed_non_core_evidence_ids)),0) "
        "FROM target_callability_evidence"
    ).fetchone()[0]
    if non_core_links:
        if "non_core_variants" not in views:
            orphan_non_core = int(non_core_links)
        else:
            orphan_non_core = connection.execute(
                """
                SELECT count(*) FROM target_callability_evidence e,
                  UNNEST(e.observed_non_core_evidence_ids) AS item(evidence_id)
                LEFT JOIN non_core_variants n USING (evidence_id)
                WHERE n.evidence_id IS NULL
                """
            ).fetchone()[0]
    metrics.update(
        {
            "target_orphan_variants": int(orphan_variants),
            "target_orphan_reference_evidence": int(orphan_reference),
            "target_orphan_callable_regions": int(orphan_regions),
            "target_orphan_non_core_evidence": int(orphan_non_core),
        }
    )
    if orphan_variants or orphan_reference or orphan_regions or orphan_non_core:
        errors.append("target callability contains orphan source-evidence identifiers")

    expected_catalog_sha = (
        manifest.get("resource_provenance", {}).get("target_catalog", {}).get("sha256")
    )
    catalog_mismatch = connection.execute(
        "SELECT count(*) FROM target_callability_evidence WHERE catalog_sha256 != ?",
        [expected_catalog_sha],
    ).fetchone()[0]
    metrics["target_catalog_manifest_mismatches"] = int(catalog_mismatch)
    if not isinstance(expected_catalog_sha, str) or catalog_mismatch:
        errors.append(
            "target catalogue digest is absent from or disagrees with manifest provenance"
        )
    return metrics


def _validate_provenance(
    connection: duckdb.DuckDBPyConnection,
    root: Path,
    errors: list[str],
) -> dict[str, Any]:
    metrics: dict[str, Any] = {}
    views = {
        row[0]
        for row in connection.execute(
            "SELECT table_name FROM information_schema.tables WHERE table_type='VIEW'"
        ).fetchall()
    }
    if "variant_provenance" not in views:
        errors.append("variant provenance extension is missing")
        return metrics
    total, distinct_variants, bad_digest, bad_source = connection.execute(
        """
        SELECT count(*), count(DISTINCT variant_id),
          count(*) FILTER (WHERE NOT regexp_full_match(source_record_digest,'^[0-9a-f]{64}$')),
          count(*) FILTER (WHERE source_ordinal < 1 OR source_alt_index < 1 OR source_gt IS NULL)
        FROM variant_provenance
        """
    ).fetchone()
    variants = connection.execute("SELECT count(*) FROM variants").fetchone()[0]
    missing_provenance = connection.execute(
        "SELECT count(*) FROM variants v LEFT JOIN variant_provenance p USING (variant_id) "
        "WHERE p.variant_id IS NULL"
    ).fetchone()[0]
    orphan_provenance = connection.execute(
        "SELECT count(*) FROM variant_provenance p LEFT JOIN variants v USING (variant_id) "
        "WHERE v.variant_id IS NULL"
    ).fetchone()[0]
    if (
        total != variants
        or distinct_variants != variants
        or missing_provenance
        or orphan_provenance
    ):
        errors.append("variant/provenance one-to-one relationship failed")
    if bad_digest or bad_source:
        errors.append("variant provenance contains invalid source identifiers")

    sources_sql = "SELECT DISTINCT source_id, source_alt_index FROM variant_provenance"
    partial_rows = 0
    if "observed_partial_variants" in views:
        partial_rows, distinct_partial, bad_partial = connection.execute("""
            SELECT count(*), count(DISTINCT (p.source_id, p.source_alt_index, p.variant_id)),
              count(*) FILTER (WHERE
                r.evidence_id IS NULL OR r.reference_observed IS NOT NULL
                OR r.callable IS DISTINCT FROM false OR p.source_ordinal != r.source_ordinal
                OR p.source_record_digest != r.source_record_digest OR p.source_pos != r.start_pos
                OR p.source_alt_index NOT BETWEEN 1 AND len(r.source_alts)
                OR p.source_alt IS DISTINCT FROM r.source_alts[p.source_alt_index]
                OR p.phased IS DISTINCT FROM r.phased
                OR len(p.gt) != 2 OR list_count(p.gt) != 1 OR NOT list_contains(p.gt,1)
                OR p.gt IS DISTINCT FROM list_transform(r.gt, a -> CASE WHEN a IS NULL THEN NULL ELSE (a=p.source_alt_index)::TINYINT END)
                OR contains(p.source_gt, '|') IS DISTINCT FROM r.phased
                OR replace(p.source_gt,'|','/') IS DISTINCT FROM array_to_string(
                    list_transform(r.gt, a -> coalesce(a::VARCHAR,'.')), '/')
                OR p.variant_id != concat(p.chrom,':',p.pos,':',p.ref,':',p.alt))
            FROM observed_partial_variants p
            LEFT JOIN reference_evidence r ON r.evidence_id=p.source_id AND r.state='no_call'
        """).fetchone()
        if distinct_partial != partial_rows or bad_partial:
            errors.append(
                "observed partial alleles differ from source no-call evidence or repeat a projection"
            )
        sources_sql += " UNION SELECT source_id, source_alt_index FROM observed_partial_variants"
    cursor = connection.execute(sources_sql + " ORDER BY source_id, source_alt_index")
    multiset = MultisetDigest()
    while batch := cursor.fetchmany(100_000):
        for source_id, source_alt_index in batch:
            multiset.add(f"{source_id}:a{source_alt_index}")
    reconciliation_path = root / "extensions" / "reconciliation.json"
    try:
        reconciliation = json.loads(reconciliation_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"invalid reconciliation.json: {exc}")
        reconciliation = {}
    actual_multiset = multiset.as_dict()
    if reconciliation.get("accepted_source_alleles") != actual_multiset:
        errors.append("reconciliation accepted-source multiset mismatch")
    if reconciliation.get("normalized_output_alleles") != actual_multiset:
        errors.append("reconciliation normalized-output multiset mismatch")
    if reconciliation.get("unique_source_alleles") != actual_multiset["count"]:
        errors.append("reconciliation unique-source count mismatch")
    if reconciliation.get("unique_normalized_variants") != variants:
        errors.append("reconciliation normalized-variant count mismatch")
    if reconciliation.get("observed_partial_variant_rows", 0) != partial_rows:
        errors.append("reconciliation observed partial variant count mismatch")
    if reconciliation.get("missing_to_reference_conversions") != 0:
        errors.append("reconciliation reports a missing-to-reference conversion")
    metrics.update(
        {
            "variant_rows": int(variants),
            "provenance_rows": int(total),
            "unique_source_alleles": actual_multiset["count"],
            "source_allele_multiset": actual_multiset,
        }
    )
    return metrics


def _validate_gene_index(
    connection: duckdb.DuckDBPyConnection, views: set[str], errors: list[str]
) -> dict[str, int]:
    if "gene_index" not in views:
        return {"rows": 0}
    row = connection.execute(
        """
        SELECT count(*) AS total,
               count(DISTINCT gene_symbol || '|' || chrom) AS distinct_gene_chrom,
               count(*) FILTER (WHERE gene_symbol IS NULL OR gene_symbol = '' OR
                                      start_pos < 1 OR end_pos < start_pos OR
                                      variant_count < 1 OR actionable_count < 0 OR
                                      actionable_count > variant_count OR pgx_count < 0 OR
                                      pgx_count > variant_count) AS invalid
        FROM gene_index
        """
    ).fetchone()
    total, distinct_gene_chrom, invalid = (int(value) for value in row)
    if total != distinct_gene_chrom:
        errors.append("gene_index contains duplicate gene/chromosome rows")
    if invalid:
        errors.append(f"gene_index contains {invalid} invalid row(s)")

    actionable_relation = (
        "SELECT DISTINCT variant_id, gene_symbol FROM acmg_secondary_findings_evidence "
        "WHERE is_actionable"
        if "acmg_secondary_findings_evidence" in views
        else "SELECT NULL::VARCHAR AS variant_id, NULL::VARCHAR AS gene_symbol WHERE false"
    )
    mismatches = connection.execute(
        f"""
        WITH selected AS (
          SELECT variant_id, chrom, pos, gene.symbol AS gene_symbol,
                 pharmacogenomics.is_pgx AS is_pgx,
                 pharmacogenomics.contributes_to_gene AS pgx_gene
          FROM variants WHERE gene.symbol IS NOT NULL AND gene.symbol != ''
        ), transcript AS (
          SELECT variant_id, chrom, pos, annotation.symbol AS gene_symbol,
                 pharmacogenomics.is_pgx AS is_pgx,
                 pharmacogenomics.contributes_to_gene AS pgx_gene
          FROM variants, UNNEST(consequence.all_annotations) AS nested(annotation)
          WHERE annotation.symbol IS NOT NULL AND annotation.symbol != ''
        ), memberships AS (
          SELECT * FROM selected UNION SELECT * FROM transcript
        ), actionable AS (
          {actionable_relation}
        ), expected AS (
          SELECT m.gene_symbol, m.chrom, min(m.pos) AS start_pos, max(m.pos) AS end_pos,
                 count(*) AS variant_count,
                 count(*) FILTER (WHERE a.variant_id IS NOT NULL) AS actionable_count,
                 count(*) FILTER (WHERE m.is_pgx AND m.pgx_gene = m.gene_symbol) AS pgx_count
          FROM memberships m LEFT JOIN actionable a
            ON m.variant_id = a.variant_id AND m.gene_symbol = a.gene_symbol
          GROUP BY m.gene_symbol, m.chrom
        )
        SELECT count(*) FROM expected e FULL OUTER JOIN gene_index g
          USING (gene_symbol, chrom)
        WHERE e.gene_symbol IS NULL OR g.gene_symbol IS NULL OR
              e.start_pos != g.start_pos OR e.end_pos != g.end_pos OR
              e.variant_count != g.variant_count OR
              e.actionable_count != g.actionable_count OR e.pgx_count != g.pgx_count
        """
    ).fetchone()[0]
    if mismatches:
        errors.append(f"gene_index disagrees with variants for {mismatches} gene row(s)")
    return {
        "rows": total,
        "distinct_gene_chrom": distinct_gene_chrom,
        "invalid_rows": invalid,
        "variant_membership_mismatches": int(mismatches),
    }


def _validate_clinical_evidence(
    connection: duckdb.DuckDBPyConnection, views: set[str], errors: list[str]
) -> dict[str, int]:
    if "clinical_evidence" not in views:
        expected = int(
            connection.execute(
                "SELECT count(*) FROM variants WHERE pathogenicity.clinvar_id IS NOT NULL"
            ).fetchone()[0]
        )
        if "acmg_secondary_findings_evidence" in views:
            expected += int(
                connection.execute(
                    "SELECT count(*) FROM acmg_secondary_findings_evidence"
                ).fetchone()[0]
            )
        if expected:
            errors.append("clinical evidence sources exist but clinical_evidence.parquet is absent")
        return {"rows": 0, "expected_source_rows": expected}

    total, distinct_ids, invalid = (
        int(value)
        for value in connection.execute(
            """
            SELECT count(*), count(DISTINCT evidence_id),
              count(*) FILTER (WHERE evidence_id IS NULL OR evidence_id = '' OR
                (finding_id IS NULL AND variant_id IS NULL) OR source IS NULL OR source = '' OR
                source_version IS NULL OR source_version = '' OR assertion IS NULL OR
                assertion = '' OR retrieved_at IS NULL)
            FROM clinical_evidence
            """
        ).fetchone()
    )
    if total != distinct_ids:
        errors.append("clinical_evidence contains duplicate evidence_id values")
    if invalid:
        errors.append(f"clinical_evidence contains {invalid} invalid row(s)")

    orphan_variants = int(
        connection.execute(
            "SELECT count(*) FROM clinical_evidence e LEFT JOIN variants v USING (variant_id) "
            "WHERE e.variant_id IS NOT NULL AND v.variant_id IS NULL"
        ).fetchone()[0]
    )
    if orphan_variants:
        errors.append(f"clinical_evidence has {orphan_variants} orphan variant link(s)")

    linked_findings = int(
        connection.execute(
            "SELECT count(*) FROM clinical_evidence WHERE finding_id IS NOT NULL"
        ).fetchone()[0]
    )
    orphan_findings = 0
    if linked_findings and "clinical_findings" not in views:
        orphan_findings = linked_findings
    elif "clinical_findings" in views:
        orphan_findings = int(
            connection.execute(
                "SELECT count(*) FROM clinical_evidence e LEFT JOIN clinical_findings f "
                "USING (finding_id) WHERE e.finding_id IS NOT NULL AND f.finding_id IS NULL"
            ).fetchone()[0]
        )
    if orphan_findings:
        errors.append(f"clinical_evidence has {orphan_findings} orphan finding link(s)")

    clinvar_mapping_mismatches = 0
    if "clinvar_assertions" in views:
        clinvar_mapping_mismatches = int(
            connection.execute(
                """
                SELECT count(*)
                FROM clinvar_assertions a FULL OUTER JOIN (
                  SELECT * FROM clinical_evidence WHERE source = 'NCBI ClinVar'
                ) e USING (evidence_id)
                WHERE a.evidence_id IS NULL OR e.evidence_id IS NULL OR
                      e.variant_id IS DISTINCT FROM a.variant_id OR
                      e.source != 'NCBI ClinVar' OR
                      e.source_record_id IS DISTINCT FROM a.clinvar_accession OR
                      e.source_version IS DISTINCT FROM a.snapshot
                """
            ).fetchone()[0]
        )
        if clinvar_mapping_mismatches:
            errors.append(
                "normative clinical evidence disagrees with the lossless ClinVar extension for "
                f"{clinvar_mapping_mismatches} row(s)"
            )
    return {
        "rows": total,
        "distinct_ids": distinct_ids,
        "invalid_rows": invalid,
        "orphan_variant_links": orphan_variants,
        "orphan_finding_links": orphan_findings,
        "clinvar_mapping_mismatches": clinvar_mapping_mismatches,
    }


def _validate_clinical_artifacts(
    connection: duckdb.DuckDBPyConnection, views: set[str], errors: list[str]
) -> dict[str, int]:
    allowed_claims = (
        "'variant_classification','carrier_status','disease_risk','pharmacogenomic',"
        "'trait_association','computational_prediction'"
    )
    finding_rows = 0
    interpretation_rows = 0
    bad_evidence_links = 0
    orphan_links = 0

    if "clinical_findings" in views:
        finding_rows, distinct_findings, invalid_findings = (
            int(value)
            for value in connection.execute(
                f"""
                SELECT count(*), count(DISTINCT finding_id),
                  count(*) FILTER (WHERE finding_id IS NULL OR finding_id = '' OR
                    (variant_id IS NULL AND gene_symbol IS NULL) OR condition IS NULL OR
                    condition = '' OR claim_type NOT IN ({allowed_claims}) OR
                    clinical_grade IS NULL OR evidence_ids IS NULL)
                FROM clinical_findings
                """
            ).fetchone()
        )
        if finding_rows != distinct_findings:
            errors.append("clinical_findings contains duplicate finding_id values")
        if invalid_findings:
            errors.append(f"clinical_findings contains {invalid_findings} invalid row(s)")
        orphan_links += int(
            connection.execute(
                "SELECT count(*) FROM clinical_findings f LEFT JOIN variants v USING (variant_id) "
                "WHERE f.variant_id IS NOT NULL AND v.variant_id IS NULL"
            ).fetchone()[0]
        )
        if "clinical_evidence" in views:
            bad_evidence_links += int(
                connection.execute(
                    """
                    SELECT count(*) FROM clinical_findings f,
                      UNNEST(f.evidence_ids) AS linked(evidence_id)
                    LEFT JOIN clinical_evidence e ON linked.evidence_id = e.evidence_id
                    WHERE e.evidence_id IS NULL
                    """
                ).fetchone()[0]
            )
        elif connection.execute(
            "SELECT count(*) FROM clinical_findings WHERE len(evidence_ids) > 0"
        ).fetchone()[0]:
            bad_evidence_links += finding_rows

    if "interpretations" in views:
        interpretation_rows, distinct_interpretations, invalid_interpretations = (
            int(value)
            for value in connection.execute(
                f"""
                SELECT count(*), count(DISTINCT interpretation_id),
                  count(*) FILTER (WHERE interpretation_id IS NULL OR interpretation_id = '' OR
                    interpretation_status NOT IN
                      ('algorithmic','draft','expert_reviewed','guideline_based') OR
                    claim_type NOT IN ({allowed_claims}) OR confirmation_required IS NULL OR
                    evidence_snapshot IS NULL OR evidence_snapshot = '' OR evidence_ids IS NULL OR
                    (interpretation_status = 'expert_reviewed' AND
                      (reviewed_by IS NULL OR reviewed_at IS NULL)) OR
                    (interpretation_status = 'algorithmic' AND
                      (reviewed_by IS NOT NULL OR reviewed_at IS NOT NULL)))
                FROM interpretations
                """
            ).fetchone()
        )
        if interpretation_rows != distinct_interpretations:
            errors.append("interpretations contains duplicate interpretation_id values")
        if invalid_interpretations:
            errors.append(f"interpretations contains {invalid_interpretations} invalid row(s)")
        orphan_links += int(
            connection.execute(
                "SELECT count(*) FROM interpretations i LEFT JOIN variants v USING (variant_id) "
                "WHERE i.variant_id IS NOT NULL AND v.variant_id IS NULL"
            ).fetchone()[0]
        )
        if "clinical_findings" in views:
            orphan_links += int(
                connection.execute(
                    "SELECT count(*) FROM interpretations i LEFT JOIN clinical_findings f "
                    "USING (finding_id) WHERE i.finding_id IS NOT NULL AND f.finding_id IS NULL"
                ).fetchone()[0]
            )
        elif connection.execute(
            "SELECT count(*) FROM interpretations WHERE finding_id IS NOT NULL"
        ).fetchone()[0]:
            orphan_links += interpretation_rows
        if "clinical_evidence" in views:
            bad_evidence_links += int(
                connection.execute(
                    """
                    SELECT count(*) FROM interpretations i,
                      UNNEST(i.evidence_ids) AS linked(evidence_id)
                    LEFT JOIN clinical_evidence e ON linked.evidence_id = e.evidence_id
                    WHERE e.evidence_id IS NULL
                    """
                ).fetchone()[0]
            )

    if orphan_links:
        errors.append(f"clinical findings/interpretations have {orphan_links} orphan link(s)")
    if bad_evidence_links:
        errors.append(
            f"clinical findings/interpretations have {bad_evidence_links} missing evidence link(s)"
        )

    acmg_materialization_mismatches = 0
    if "acmg_secondary_findings_evidence" in views:
        if "clinical_findings" not in views or "interpretations" not in views:
            acmg_materialization_mismatches = int(
                connection.execute(
                    "SELECT count(*) FROM acmg_secondary_findings_evidence"
                ).fetchone()[0]
            )
        else:
            acmg_materialization_mismatches = int(
                connection.execute(
                    """
                    SELECT count(*) FROM acmg_secondary_findings_evidence a
                    WHERE NOT EXISTS (
                      SELECT 1 FROM clinical_findings f
                      WHERE f.variant_id = a.variant_id AND f.gene_symbol = a.gene_symbol AND
                            f.condition = a.disease AND f.claim_type = 'disease_risk' AND
                            NOT f.clinical_grade
                    ) OR NOT EXISTS (
                      SELECT 1 FROM interpretations i
                      WHERE i.variant_id = a.variant_id AND i.gene_symbol = a.gene_symbol AND
                            i.condition = a.disease AND i.claim_type = 'disease_risk' AND
                            i.interpretation_status = 'algorithmic' AND i.confirmation_required
                    )
                    """
                ).fetchone()[0]
            )
        if acmg_materialization_mismatches:
            errors.append(
                "ACMG decisions lack conservative normative finding/interpretation rows for "
                f"{acmg_materialization_mismatches} decision(s)"
            )
    return {
        "finding_rows": finding_rows,
        "interpretation_rows": interpretation_rows,
        "orphan_links": orphan_links,
        "missing_evidence_links": bad_evidence_links,
        "acmg_materialization_mismatches": acmg_materialization_mismatches,
    }


def _validate_clinvar_assertions(
    connection: duckdb.DuckDBPyConnection, views: set[str], errors: list[str]
) -> dict[str, int]:
    if "clinvar_assertions" not in views:
        matched = connection.execute(
            "SELECT count(*) FROM variants WHERE pathogenicity.clinvar_id IS NOT NULL"
        ).fetchone()[0]
        if matched:
            errors.append(
                "variants contain ClinVar fields but extensions/clinvar_assertions.parquet is absent"
            )
        return {"rows": 0, "matched_variants": int(matched)}
    row = connection.execute(
        """
        SELECT count(*) AS total, count(DISTINCT evidence_id) AS distinct_ids,
          count(*) FILTER (WHERE evidence_id IS NULL OR variant_id IS NULL OR
            source != 'NCBI ClinVar' OR snapshot IS NULL OR clinvar_variation_id IS NULL OR
            clinvar_accession IS NULL OR review_stars NOT BETWEEN 0 AND 4 OR
            raw_significance IS NULL OR review_status IS NULL OR
            conflict_classifications IS NULL OR condition_names IS NULL OR
            condition_identifiers IS NULL OR submitter_accessions IS NULL OR
            gene_info IS NULL OR rsids IS NULL OR origin_codes IS NULL OR hgvs IS NULL OR
            NOT regexp_full_match(source_record_digest, '^[0-9a-f]{64}$')) AS invalid
        FROM clinvar_assertions
        """
    ).fetchone()
    total, distinct_ids, invalid = (int(value) for value in row)
    if total != distinct_ids:
        errors.append("ClinVar assertions contain duplicate evidence_id values")
    if invalid:
        errors.append(f"ClinVar assertions contain {invalid} invalid row(s)")
    orphan = connection.execute(
        "SELECT count(*) FROM clinvar_assertions e LEFT JOIN variants v USING (variant_id) "
        "WHERE v.variant_id IS NULL"
    ).fetchone()[0]
    missing = connection.execute(
        "SELECT count(*) FROM variants v LEFT JOIN clinvar_assertions e USING (variant_id) "
        "WHERE v.pathogenicity.clinvar_id IS NOT NULL AND e.variant_id IS NULL"
    ).fetchone()[0]
    conflict_mismatch = connection.execute(
        """
        WITH evidence AS (
          SELECT variant_id, bool_or(has_conflicts) AS has_conflicts,
                 max(review_stars) AS review_stars
          FROM clinvar_assertions GROUP BY variant_id
        )
        SELECT count(*) FROM variants v JOIN evidence e USING (variant_id)
        WHERE (e.has_conflicts AND NOT v.pathogenicity.clinvar_has_conflicts) OR
              v.pathogenicity.clinvar_review_stars IS DISTINCT FROM e.review_stars
        """
    ).fetchone()[0]
    if orphan:
        errors.append(f"ClinVar assertions have {orphan} orphan row(s)")
    if missing:
        errors.append(f"ClinVar-populated variants lack evidence rows for {missing} variant(s)")
    if conflict_mismatch:
        errors.append(
            f"ClinVar conflict/review aggregation disagrees for {conflict_mismatch} variant(s)"
        )
    matched = connection.execute(
        "SELECT count(DISTINCT variant_id) FROM clinvar_assertions"
    ).fetchone()[0]
    return {
        "rows": total,
        "distinct_ids": distinct_ids,
        "invalid_rows": invalid,
        "matched_variants": int(matched),
        "orphan_rows": int(orphan),
        "missing_evidence_variants": int(missing),
        "aggregation_mismatches": int(conflict_mismatch),
    }


def _validate_acmg_evidence(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, int]:
    definition = manifest.get("consensus_definitions", {}).get("is_actionable", {})
    true_variants = int(
        connection.execute(
            "SELECT count(*) FROM variants WHERE pathogenicity.is_actionable"
        ).fetchone()[0]
    )
    if "acmg_secondary_findings_evidence" not in views:
        if true_variants:
            errors.append("actionable variants exist without ACMG policy evidence")
        return {"rows": 0, "actionable_variants": true_variants}

    expected_gene_digest = str(definition.get("gene_list_sha256", ""))
    row = connection.execute(
        """
        SELECT count(*) AS total, count(DISTINCT decision_id) AS distinct_ids,
          count(*) FILTER (WHERE decision_id IS NULL OR variant_id IS NULL OR
            gene_symbol IS NULL OR policy_key IS NULL OR disease IS NULL OR
            disease_identifiers IS NULL OR inheritance IS NULL OR variants_to_report IS NULL OR
            decision_status IS NULL OR is_actionable IS NULL OR clinical_grade IS NULL OR
            confirmation_required IS DISTINCT FROM true OR analytical_confidence IS NULL OR
            genotype IS NULL OR clinvar_significance IS NULL OR
            clinvar_snapshot IS NULL OR
            clinvar_review_stars NOT BETWEEN 0 AND 4 OR clinvar_has_conflicts IS NULL OR
            clinvar_evidence_ids IS NULL OR clinvar_accessions IS NULL OR
            matched_condition_identifiers IS NULL OR policy_version IS NULL OR
            policy_snapshot IS NULL OR policy_gene_set_sha256 IS NULL OR reasons IS NULL)
            AS invalid,
          count(*) FILTER (WHERE clinical_grade OR
            is_actionable != (decision_status = 'REPORTABLE_CANDIDATE')) AS bad_decision_semantics,
          count(*) FILTER (WHERE is_actionable AND
            (clinvar_has_conflicts OR analytical_confidence != 'high')) AS unsafe_actionable
        FROM acmg_secondary_findings_evidence
        """
    ).fetchone()
    total, distinct_ids, invalid, bad_semantics, unsafe = (int(value) for value in row)
    if total != distinct_ids:
        errors.append("ACMG secondary-findings evidence has duplicate decision IDs")
    if invalid:
        errors.append(f"ACMG secondary-findings evidence contains {invalid} invalid row(s)")
    if bad_semantics:
        errors.append(f"ACMG decision/status semantics failed for {bad_semantics} row(s)")
    if unsafe:
        errors.append(f"unsafe ACMG actionable evidence found for {unsafe} row(s)")
    bad_policy_digest = int(
        connection.execute(
            "SELECT count(*) FROM acmg_secondary_findings_evidence "
            "WHERE policy_gene_set_sha256 != ?",
            [expected_gene_digest],
        ).fetchone()[0]
    )
    if bad_policy_digest:
        errors.append(f"ACMG policy digest mismatch for {bad_policy_digest} row(s)")
    orphan = int(
        connection.execute(
            "SELECT count(*) FROM acmg_secondary_findings_evidence e "
            "LEFT JOIN variants v USING (variant_id) WHERE v.variant_id IS NULL"
        ).fetchone()[0]
    )
    if orphan:
        errors.append(f"ACMG secondary-findings evidence has {orphan} orphan row(s)")
    missing_clinvar_links = 0
    if "clinical_evidence" in views:
        missing_clinvar_links = int(
            connection.execute(
                """
                SELECT count(*) FROM acmg_secondary_findings_evidence a,
                  UNNEST(a.clinvar_evidence_ids) AS linked(evidence_id)
                LEFT JOIN clinical_evidence e ON linked.evidence_id = e.evidence_id
                WHERE e.evidence_id IS NULL OR e.source != 'NCBI ClinVar'
                """
            ).fetchone()[0]
        )
    else:
        missing_clinvar_links = total
    if missing_clinvar_links:
        errors.append(
            f"ACMG evidence has {missing_clinvar_links} unresolved ClinVar evidence link(s)"
        )
    mismatch = int(
        connection.execute(
            """
            WITH expected AS (
              SELECT variant_id, bool_or(is_actionable) AS is_actionable
              FROM acmg_secondary_findings_evidence GROUP BY variant_id
            )
            SELECT count(*) FROM variants v LEFT JOIN expected e USING (variant_id)
            WHERE v.pathogenicity.is_actionable IS DISTINCT FROM coalesce(e.is_actionable, false)
            """
        ).fetchone()[0]
    )
    if mismatch:
        errors.append(f"variant/ACMG actionability flags disagree for {mismatch} variant(s)")
    return {
        "rows": total,
        "distinct_ids": distinct_ids,
        "invalid_rows": invalid,
        "bad_decision_semantics": bad_semantics,
        "unsafe_actionable_rows": unsafe,
        "bad_policy_digest_rows": bad_policy_digest,
        "orphan_rows": orphan,
        "missing_clinvar_evidence_links": missing_clinvar_links,
        "actionability_mismatches": mismatch,
        "actionable_variants": true_variants,
    }


def _validate_pgs(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, Any]:
    has_models = "pgs_model_audit" in views
    has_loci = "pgs_locus_audit" in views
    has_prs = "prs" in views
    if has_models != has_loci:
        errors.append("PGS model and locus audit extensions must appear together")
        return {"model_rows": 0, "locus_rows": 0, "published_rows": 0}
    if not has_models:
        if has_prs:
            errors.append("prs.parquet exists without PGS audit extensions")
        return {"model_rows": 0, "locus_rows": 0, "published_rows": 0}

    model_total, model_distinct, invalid_models = connection.execute(
        """
        SELECT count(*), count(DISTINCT pgs_id),
          count(*) FILTER (WHERE pgs_id IS NULL OR trait IS NULL OR
            NOT regexp_full_match(scoring_file_sha256, '^[0-9a-f]{64}$') OR
            scoring_format_version != '2.0' OR source_version IS NULL OR
            training_date IS NULL OR model_type != 'linear' OR license IS NULL OR
            redistribution_allowed IS NULL OR expected_loci <= 0 OR supported_loci <= 0 OR
            supported_loci > expected_loci OR strand_flipped_loci < 0 OR
            strand_flipped_loci > supported_loci OR alternate_called_loci < 0 OR
            reference_called_loci < 0 OR missing_loci < 0 OR ambiguous_loci < 0 OR
            low_confidence_loci < 0 OR unsupported_loci < 0 OR
            expected_absolute_weight <= 0 OR used_absolute_weight < 0 OR
            missing_absolute_weight < -1e-9 OR weight_completeness NOT BETWEEN 0 AND 1 OR
            max_missing_weight_fraction NOT BETWEEN 0 AND 1 OR score_valid IS NULL OR
            calibrated IS NULL OR publication_status NOT IN
              ('published','valid_uncalibrated','calibration_inapplicable','incomplete') OR
            contributing_variant_count < 0 OR
            contributing_variant_count > alternate_called_loci OR
            contributing_variant_ids_complete IS NULL OR
            (contributing_variant_ids_complete AND
              (contributing_variant_ids IS NULL OR
               list_count(contributing_variant_ids) != contributing_variant_count OR
               list_count(list_distinct(contributing_variant_ids)) !=
                 contributing_variant_count OR
               list_sort(contributing_variant_ids) != contributing_variant_ids)) OR
            (NOT contributing_variant_ids_complete AND contributing_variant_ids IS NOT NULL) OR
            NOT regexp_full_match(target_catalog_sha256, '^[0-9a-f]{64}$'))
        FROM pgs_model_audit
        """
    ).fetchone()
    model_total = int(model_total)
    model_distinct = int(model_distinct)
    invalid_models = int(invalid_models)
    if model_total != model_distinct:
        errors.append("PGS model audit contains duplicate pgs_id values")
    if invalid_models:
        errors.append(f"PGS model audit contains {invalid_models} invalid row(s)")

    locus_total, locus_distinct, invalid_loci = connection.execute(
        """
        SELECT count(*), count(DISTINCT (pgs_id, scoring_row_number)),
          count(*) FILTER (WHERE pgs_id IS NULL OR scoring_row_number <= 0 OR
            NOT regexp_full_match(scoring_row_digest, '^[0-9a-f]{64}$') OR
            effect_allele IS NULL OR source_effect_allele IS NULL OR
            NOT isfinite(effect_weight) OR status NOT IN
              ('alternate_called','reference_called','missing_no_call',
               'missing_interval_only','missing_unknown','ambiguous_other_variant',
               'ambiguous_genotype','low_confidence','unsupported') OR
            observed_variant_ids IS NULL OR reason IS NULL OR
            ((status IN ('alternate_called','reference_called')) !=
             (dosage IS NOT NULL AND contribution IS NOT NULL)) OR
            (dosage IS NOT NULL AND
             (observed_ploidy NOT IN (1,2) OR dosage < 0 OR dosage > observed_ploidy)) OR
            (status = 'alternate_called' AND resolution != 'non_reference_observed') OR
            (status = 'reference_called' AND resolution != 'reference_observed') OR
            (status = 'unsupported' AND
             (target_id IS NOT NULL OR callability_id IS NOT NULL OR
              target_catalog_sha256 IS NOT NULL OR target_catalog_row_digest IS NOT NULL)) OR
            (status != 'unsupported' AND
             (target_id IS NULL OR callability_id IS NULL OR
              strand_flipped IS NULL OR
              NOT regexp_full_match(target_catalog_sha256, '^[0-9a-f]{64}$') OR
              NOT regexp_full_match(target_catalog_row_digest, '^[0-9a-f]{64}$'))))
        FROM pgs_locus_audit
        """
    ).fetchone()
    locus_total = int(locus_total)
    locus_distinct = int(locus_distinct)
    invalid_loci = int(invalid_loci)
    if locus_total != locus_distinct:
        errors.append("PGS locus audit contains duplicate model/row identifiers")
    if invalid_loci:
        errors.append(f"PGS locus audit contains {invalid_loci} invalid row(s)")

    orphan_models = int(
        connection.execute(
            "SELECT count(*) FROM pgs_locus_audit l LEFT JOIN pgs_model_audit m USING (pgs_id) "
            "WHERE m.pgs_id IS NULL"
        ).fetchone()[0]
    )
    if orphan_models:
        errors.append(f"PGS locus audit has {orphan_models} orphan model row(s)")

    target_mismatches = locus_total
    if "target_callability_evidence" in views:
        target_mismatches = int(
            connection.execute(
                """
                SELECT count(*) FROM pgs_locus_audit l
                LEFT JOIN target_callability_evidence e USING (target_id)
                WHERE l.status != 'unsupported' AND
                  (e.target_id IS NULL OR l.callability_id != e.callability_id OR
                   l.target_catalog_sha256 != e.catalog_sha256 OR
                   l.target_catalog_row_digest != e.catalog_row_digest OR
                   l.chrom IS DISTINCT FROM e.chrom OR l.pos IS DISTINCT FROM e.pos OR
                   l.ref IS DISTINCT FROM e.ref OR l.alt IS DISTINCT FROM e.alt OR
                   l.variant_id IS DISTINCT FROM e.variant_id OR
                   l.resolution IS DISTINCT FROM e.resolution)
                """
            ).fetchone()[0]
        )
    if target_mismatches:
        errors.append(f"PGS locus audit has {target_mismatches} target-proof mismatch(es)")

    carried_mismatches = int(
        connection.execute(
            """
            SELECT count(*) FROM pgs_locus_audit l
            LEFT JOIN variants v USING (variant_id)
            WHERE l.status = 'alternate_called' AND v.variant_id IS NULL
            """
        ).fetchone()[0]
    )
    if carried_mismatches:
        errors.append(f"PGS alternate dosage has {carried_mismatches} missing variant link(s)")

    embedded_id_mismatches = int(
        connection.execute(
            """
            WITH contributing AS (
              SELECT DISTINCT pgs_id, variant_id
              FROM pgs_locus_audit
              WHERE status = 'alternate_called' AND contribution IS NOT NULL
                AND variant_id IS NOT NULL
            ), embedded AS (
              SELECT m.pgs_id, item.variant_id
              FROM pgs_model_audit m
              CROSS JOIN unnest(m.contributing_variant_ids) AS item(variant_id)
              WHERE m.contributing_variant_ids_complete
            )
            SELECT count(*)
            FROM embedded e
            LEFT JOIN contributing c
              ON c.pgs_id = e.pgs_id AND c.variant_id = e.variant_id
            WHERE c.variant_id IS NULL
            """
        ).fetchone()[0]
    )
    if embedded_id_mismatches:
        errors.append(
            f"PGS model audit embeds {embedded_id_mismatches} non-contributing variant ID(s)"
        )

    aggregate_rows = connection.execute(
        """
        SELECT pgs_id, count(*) AS expected_loci,
          count(*) FILTER (WHERE status != 'unsupported') AS supported_loci,
          count(*) FILTER (WHERE strand_flipped) AS strand_flipped_loci,
          count(*) FILTER (WHERE status = 'alternate_called') AS alternate_called_loci,
          count(*) FILTER (WHERE status = 'reference_called') AS reference_called_loci,
          count(*) FILTER (WHERE starts_with(status, 'missing_')) AS missing_loci,
          count(*) FILTER (WHERE starts_with(status, 'ambiguous_')) AS ambiguous_loci,
          count(*) FILTER (WHERE status = 'low_confidence') AS low_confidence_loci,
          count(*) FILTER (WHERE status = 'unsupported') AS unsupported_loci,
          fsum(abs(effect_weight)) AS expected_absolute_weight,
          coalesce(fsum(abs(effect_weight)) FILTER (WHERE contribution IS NOT NULL), 0)
            AS used_absolute_weight,
          coalesce(fsum(contribution), 0) AS partial_raw_score,
          count(DISTINCT variant_id) FILTER (
            WHERE status = 'alternate_called' AND contribution IS NOT NULL
              AND variant_id IS NOT NULL
          ) AS contributing_variant_count
        FROM pgs_locus_audit GROUP BY pgs_id ORDER BY pgs_id
        """
    ).fetchall()
    aggregates = {
        str(row[0]): {
            "expected_loci": int(row[1]),
            "supported_loci": int(row[2]),
            "strand_flipped_loci": int(row[3]),
            "alternate_called_loci": int(row[4]),
            "reference_called_loci": int(row[5]),
            "missing_loci": int(row[6]),
            "ambiguous_loci": int(row[7]),
            "low_confidence_loci": int(row[8]),
            "unsupported_loci": int(row[9]),
            "expected_absolute_weight": float(row[10]),
            "used_absolute_weight": float(row[11]),
            "partial_raw_score": float(row[12]),
            "contributing_variant_count": int(row[13]),
        }
        for row in aggregate_rows
    }
    accounting_mismatches = 0
    model_cursor = connection.execute(
        """
        SELECT pgs_id, expected_loci, supported_loci, strand_flipped_loci,
          alternate_called_loci, reference_called_loci, missing_loci, ambiguous_loci,
          low_confidence_loci,
          unsupported_loci, expected_absolute_weight, used_absolute_weight,
          missing_absolute_weight, weight_completeness, max_missing_weight_fraction,
          partial_raw_score, score_valid, calibrated, percentile, reference_population,
          calibration_source, calibration_source_version, calibration_date,
          calibration_sha256, calibration_max_missing_weight_fraction,
          publication_status, contributing_variant_count,
          contributing_variant_ids_complete, target_catalog_sha256
        FROM pgs_model_audit ORDER BY pgs_id
        """
    )
    model_names = [description[0] for description in model_cursor.description]
    for values in model_cursor.fetchall():
        model = dict(zip(model_names, values, strict=True))
        pgs_id = str(model["pgs_id"])
        aggregate = aggregates.get(pgs_id)
        if aggregate is None:
            accounting_mismatches += 1
            continue
        count_names = (
            "expected_loci",
            "supported_loci",
            "strand_flipped_loci",
            "alternate_called_loci",
            "reference_called_loci",
            "missing_loci",
            "ambiguous_loci",
            "low_confidence_loci",
            "unsupported_loci",
        )
        declared_counts = tuple(int(model[name]) for name in count_names)
        observed_counts = tuple(int(aggregate[name]) for name in count_names)
        expected_weight = aggregate["expected_absolute_weight"]
        used_weight = aggregate["used_absolute_weight"]
        raw_score = aggregate["partial_raw_score"]
        missing_weight = max(0.0, expected_weight - used_weight)
        completeness = min(1.0, used_weight / expected_weight)
        valid = (1.0 - completeness) <= float(model["max_missing_weight_fraction"]) + 1e-12
        calibrated = bool(model["calibrated"])
        published = str(model["publication_status"]) == "published"
        calibration_complete = all(
            model[name] is not None
            for name in (
                "percentile",
                "reference_population",
                "calibration_source",
                "calibration_source_version",
                "calibration_date",
                "calibration_sha256",
                "calibration_max_missing_weight_fraction",
            )
        )
        if (
            declared_counts != observed_counts
            or int(model["contributing_variant_count"]) != aggregate["contributing_variant_count"]
            or not math.isclose(
                float(model["expected_absolute_weight"]),
                expected_weight,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
            or not math.isclose(
                float(model["used_absolute_weight"]),
                used_weight,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
            or not math.isclose(
                float(model["missing_absolute_weight"]),
                missing_weight,
                rel_tol=1e-10,
                abs_tol=1e-10,
            )
            or not math.isclose(
                float(model["weight_completeness"]),
                completeness,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
            or not math.isclose(
                float(model["partial_raw_score"]), raw_score, rel_tol=1e-12, abs_tol=1e-12
            )
            or bool(model["score_valid"]) != valid
            or calibrated != published
            or (calibrated and not calibration_complete)
            or (not calibrated and model["percentile"] is not None)
        ):
            accounting_mismatches += 1
    if accounting_mismatches:
        errors.append(f"PGS model accounting mismatch for {accounting_mismatches} model(s)")

    published_models = int(
        connection.execute(
            "SELECT count(*) FROM pgs_model_audit WHERE publication_status='published'"
        ).fetchone()[0]
    )
    prs_mismatches = 0
    if has_prs:
        prs_total, prs_distinct = connection.execute(
            "SELECT count(*), count(DISTINCT training_source) FROM prs"
        ).fetchone()
        if int(prs_total) != int(prs_distinct):
            errors.append("prs.parquet contains duplicate training_source values")
        prs_mismatches = int(
            connection.execute(
                """
                SELECT count(*) FROM prs p
                FULL OUTER JOIN pgs_model_audit m ON p.training_source = m.pgs_id
                  AND m.publication_status = 'published'
                WHERE p.training_source IS NULL OR m.pgs_id IS NULL OR
                  p.trait != m.trait OR abs(p.score_value::DOUBLE - m.partial_raw_score) > 1e-5 OR
                  abs(p.percentile::DOUBLE - m.percentile) > 1e-5 OR
                  p.reference_population != m.reference_population OR
                  p.training_date != m.training_date OR p.model_type != m.model_type OR
                  (m.contributing_variant_ids_complete AND
                    list_sort(coalesce(p.contributing_variants, []::VARCHAR[])) !=
                      list_sort(m.contributing_variant_ids)) OR
                  (NOT m.contributing_variant_ids_complete AND
                    p.contributing_variants IS NOT NULL)
                """
            ).fetchone()[0]
        )
    elif published_models:
        prs_mismatches = published_models
    if prs_mismatches:
        errors.append(f"normative PRS/model audit mismatch for {prs_mismatches} row(s)")

    provenance = manifest.get("resource_provenance", {}).get("polygenic_scores")
    provenance_mismatches = 0
    if not isinstance(provenance, dict):
        provenance_mismatches = model_total
    else:
        expected_target_digest = provenance.get("target_catalog_sha256")
        provenance_mismatches += int(
            connection.execute(
                "SELECT count(*) FROM pgs_model_audit WHERE target_catalog_sha256 != ?",
                [expected_target_digest],
            ).fetchone()[0]
        )
        declared_models = {
            str(row.get("pgs_id")): str(row.get("scoring_file_sha256"))
            for row in provenance.get("models", [])
            if isinstance(row, dict)
        }
        observed_models = {
            str(pgs_id): str(digest)
            for pgs_id, digest in connection.execute(
                "SELECT pgs_id, scoring_file_sha256 FROM pgs_model_audit"
            ).fetchall()
        }
        if declared_models != observed_models:
            provenance_mismatches += 1
    if provenance_mismatches:
        errors.append(f"PGS manifest provenance mismatch ({provenance_mismatches})")

    return {
        "model_rows": model_total,
        "locus_rows": locus_total,
        "published_rows": published_models,
        "invalid_model_rows": invalid_models,
        "invalid_locus_rows": invalid_loci,
        "orphan_model_rows": orphan_models,
        "target_proof_mismatches": target_mismatches,
        "carried_variant_mismatches": carried_mismatches,
        "embedded_contributing_id_mismatches": embedded_id_mismatches,
        "accounting_mismatches": accounting_mismatches,
        "normative_mismatches": prs_mismatches,
        "provenance_mismatches": provenance_mismatches,
    }


def _validate_pgx(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, Any]:
    """Recompute every bundled PGx locus and validate the gated projection."""

    has_loci = "pgx_locus_audit" in views
    has_genes = "pgx_gene_call_audit" in views
    has_guidance = "pgx_guidance_audit" in views
    has_core = "pharmacogenomics" in views
    if not any((has_loci, has_genes, has_guidance, has_core)):
        return {
            "locus_rows": 0,
            "gene_rows": 0,
            "guidance_rows": 0,
            "published_rows": 0,
        }
    if not has_loci or not has_genes:
        errors.append("PGx locus and gene-call audit extensions must appear together")
        return {
            "locus_rows": 0,
            "gene_rows": 0,
            "guidance_rows": 0,
            "published_rows": 0,
        }
    if "target_callability_evidence" not in views or "variants" not in views:
        errors.append("PGx audit requires target callability and carried-variant evidence")
        return {
            "locus_rows": 0,
            "gene_rows": 0,
            "guidance_rows": 0,
            "published_rows": 0,
        }

    locus_total, locus_distinct, allele_total, invalid_loci = (
        int(value)
        for value in connection.execute(
            """
            SELECT count(*), count(DISTINCT locus_id), sum(list_count(alts)),
              count(*) FILTER (WHERE
                locus_id IS NULL OR record_number <= 0 OR
                NOT regexp_full_match(source_record_digest, '^[0-9a-f]{64}$') OR
                NOT regexp_full_match(positions_vcf_sha256, '^[0-9a-f]{64}$') OR
                chrom IS NULL OR pos < 1 OR ref IS NULL OR
                NOT regexp_full_match(ref, '^[ACGT]+$') OR
                alts IS NULL OR list_count(alts) < 1 OR
                list_count(alts) != list_count(target_ids) OR
                list_count(alts) != list_count(target_variant_ids) OR
                list_count(alts) != list_count(normalization_changed) OR
                list_count(alts) != list_count(callability_ids) OR
                list_count(alts) != list_count(evidence_resolutions) OR
                list_count(list_distinct(alts)) != list_count(alts) OR
                list_count(list_distinct(target_ids)) != list_count(target_ids) OR
                list_count(list_distinct(target_variant_ids)) != list_count(target_variant_ids) OR
                list_count(list_filter(target_variant_ids, value ->
                  NOT regexp_full_match(value, '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M):[0-9]+:[ACGT]+:[ACGT]+$'))) > 0 OR
                list_count(list_distinct(callability_ids)) != list_count(callability_ids) OR
                genes IS NULL OR observed_variant_ids IS NULL OR reason IS NULL OR reason = '' OR
                NOT regexp_full_match(materialized_record_digest, '^[0-9a-f]{64}$') OR
                status NOT IN (
                  'variant_call','reference_call','no_call_structurally_complex',
                  'no_call_other_variant','no_call_low_confidence_variant',
                  'no_call_inconsistent_genotype','no_call_low_confidence_reference',
                  'no_call_source_missing','no_call_interval_only',
                  'no_call_inconsistent_reference','no_call_unknown') OR
                (eligible_for_matcher != (status IN ('variant_call','reference_call'))) OR
                ((status IN ('variant_call','reference_call')) != (genotype IS NOT NULL)) OR
                (genotype IS NOT NULL AND
                  (observed_ploidy NOT IN (1,2) OR list_count(genotype) != observed_ploidy OR
                   depth IS NULL OR depth < 0 OR genotype_quality IS NULL OR
                   genotype_quality < 0 OR phased IS NULL)) OR
                (genotype IS NULL AND (phased IS NOT NULL OR phase_set IS NOT NULL)) OR
                (status = 'reference_call' AND
                  (list_min(genotype) != 0 OR list_max(genotype) != 0 OR
                   list_count(list_filter(evidence_resolutions, value -> value != 'reference_observed')) > 0)) OR
                (status = 'variant_call' AND
                  (list_min(genotype) < 0 OR list_max(genotype) < 1 OR
                   list_max(genotype) > list_count(alts) OR
                   NOT list_contains(evidence_resolutions, 'non_reference_observed'))) OR
                (structurally_complex != (status = 'no_call_structurally_complex')))
            FROM pgx_locus_audit
            """
        ).fetchone()
    )
    if locus_total != locus_distinct:
        errors.append("PGx locus audit contains duplicate locus_id values")
    if locus_total != PHARMCAT_POSITIONS_RECORDS or allele_total != PHARMCAT_POSITIONS_ALLELES:
        errors.append("PGx locus audit does not cover the exact pinned PharmCAT positions boundary")
    if invalid_loci:
        errors.append(f"PGx locus audit contains {invalid_loci} invalid row(s)")

    target_mismatches = int(
        connection.execute(
            """
            SELECT count(*)
            FROM pgx_locus_audit l,
              UNNEST(l.target_ids) WITH ORDINALITY AS item(target_id, allele_index)
            LEFT JOIN target_callability_evidence e ON e.target_id = item.target_id
            WHERE e.target_id IS NULL OR
              e.callability_id IS DISTINCT FROM list_extract(l.callability_ids, allele_index) OR
              e.resolution IS DISTINCT FROM list_extract(l.evidence_resolutions, allele_index) OR
              e.variant_id IS DISTINCT FROM list_extract(l.target_variant_ids, allele_index)
            """
        ).fetchone()[0]
    )
    if target_mismatches:
        errors.append(f"PGx locus audit has {target_mismatches} target-proof mismatch(es)")

    orphan_variants = int(
        connection.execute(
            """
            SELECT count(*)
            FROM pgx_locus_audit l,
              UNNEST(l.observed_variant_ids) AS item(variant_id)
            LEFT JOIN variants v ON v.variant_id = item.variant_id
            WHERE v.variant_id IS NULL
            """
        ).fetchone()[0]
    )
    if orphan_variants:
        errors.append(f"PGx locus audit has {orphan_variants} orphan variant identifier(s)")

    audit_cursor = connection.execute("SELECT * FROM pgx_locus_audit ORDER BY record_number")
    audit_names = [description[0] for description in audit_cursor.description]
    audit_rows = {
        str(row["locus_id"]): row
        for values in audit_cursor.fetchall()
        for row in [dict(zip(audit_names, values, strict=True))]
    }
    proof_cursor = connection.execute(
        """
        WITH expanded AS (
          SELECT l.locus_id, l.record_number, l.source_record_digest,
            l.positions_vcf_sha256, l.chrom, l.pos, l.source_id, l.ref, l.genes,
            l.structurally_complex, item.allele_index,
            list_extract(l.alts, item.allele_index) AS alt, item.target_id,
            list_extract(l.target_variant_ids, item.allele_index) AS target_variant_id,
            list_extract(l.normalization_changed, item.allele_index) AS normalization_changed
          FROM pgx_locus_audit l
          CROSS JOIN UNNEST(l.target_ids) WITH ORDINALITY AS item(target_id, allele_index)
        )
        SELECT x.locus_id, x.record_number, x.source_record_digest,
          x.positions_vcf_sha256, x.chrom, x.pos, x.source_id, x.ref, x.genes,
          x.structurally_complex, x.allele_index, x.alt, x.target_id,
          x.target_variant_id, x.normalization_changed,
          e.callability_id, e.resolution, e.reference_observed, e.callable,
          e.call_confidence, e.depth, e.genotype_quality, e.observed_ploidy,
          e.observed_variant_ids, v.genotype.gt AS variant_gt,
          v.genotype.phased AS variant_phased
        FROM expanded x
        LEFT JOIN target_callability_evidence e ON e.target_id = x.target_id
        LEFT JOIN variants v ON v.variant_id = x.target_variant_id
        ORDER BY x.record_number, x.allele_index
        """
    )
    proof_names = [description[0] for description in proof_cursor.description]
    proof_by_locus: dict[str, list[dict[str, Any]]] = {}
    for values in proof_cursor.fetchall():
        row = dict(zip(proof_names, values, strict=True))
        row["position_of_interest"] = False
        proof_by_locus.setdefault(str(row["locus_id"]), []).append(row)

    recomputation_mismatches = 0
    comparable_fields = (
        "status",
        "genotype",
        "phased",
        "phase_set",
        "observed_ploidy",
        "depth",
        "observed_variant_ids",
        "structurally_complex",
        "eligible_for_matcher",
        "reason",
        "materialized_record_digest",
    )
    for locus_id, audit in audit_rows.items():
        proof = proof_by_locus.get(locus_id, [])
        if len(proof) != len(audit["alts"]):
            recomputation_mismatches += 1
            continue
        try:
            recomputed = _resolve_pgx_locus(proof)
        except Exception:
            recomputation_mismatches += 1
            continue
        mismatch = any(audit[field] != recomputed[field] for field in comparable_fields)
        left_gq = audit["genotype_quality"]
        right_gq = recomputed["genotype_quality"]
        if (left_gq is None) != (right_gq is None) or (
            left_gq is not None
            and right_gq is not None
            and not math.isclose(float(left_gq), float(right_gq), rel_tol=1e-6, abs_tol=1e-6)
        ):
            mismatch = True
        recomputation_mismatches += int(mismatch)
    if recomputation_mismatches:
        errors.append(
            f"PGx evidence recomputation disagrees for {recomputation_mismatches} locus/loci"
        )

    gene_total, gene_distinct, invalid_genes, published_rows = (
        int(value)
        for value in connection.execute(
            """
            SELECT count(*), count(DISTINCT gene_symbol),
              count(*) FILTER (WHERE gene_symbol IS NULL OR gene_symbol = '' OR
                source_diplotypes IS NULL OR phenotypes IS NULL OR missing_positions IS NULL OR
                cpic_levels IS NULL OR affected_drugs IS NULL OR guideline_urls IS NULL OR
                caller != 'PharmCAT' OR caller_version != ? OR data_version IS NULL OR
                NOT regexp_full_match(calls_tsv_sha256, '^[0-9a-f]{64}$') OR
                NOT regexp_full_match(report_json_sha256, '^[0-9a-f]{64}$') OR reason IS NULL OR
                call_status NOT IN ('published','complex_requires_validated_outside_call',
                  'ambiguous_diplotype','incomplete_positions','phenotype_unavailable',
                  'no_normative_cpic_level') OR
                list_count(list_distinct(cpic_levels)) != list_count(cpic_levels) OR
                list_count(list_filter(cpic_levels, value -> value NOT IN ('A','B','C','D'))) > 0 OR
                list_sort(cpic_levels) != cpic_levels OR
                list_count(list_distinct(affected_drugs)) != list_count(affected_drugs) OR
                list_sort(affected_drugs) != affected_drugs OR
                (structurally_complex !=
                  (call_status = 'complex_requires_validated_outside_call')) OR
                (call_status = 'published' AND
                  (list_count(source_diplotypes) != 1 OR selected_diplotype IS NULL OR
                   selected_diplotype != list_extract(source_diplotypes, 1) OR
                   list_count(missing_positions) != 0 OR selected_phenotype IS NULL OR
                   list_count(cpic_levels) = 0 OR list_count(affected_drugs) = 0))),
              count(*) FILTER (WHERE call_status = 'published')
            FROM pgx_gene_call_audit
            """,
            [PHARMCAT_VERSION],
        ).fetchone()
    )
    if gene_total != gene_distinct:
        errors.append("PGx gene-call audit contains duplicate gene symbols")
    if invalid_genes:
        errors.append(f"PGx gene-call audit contains {invalid_genes} invalid row(s)")

    normative_mismatches = 0
    if has_core:
        normative_mismatches = int(
            connection.execute(
                """
                SELECT count(*)
                FROM pharmacogenomics p
                FULL OUTER JOIN (
                  SELECT * FROM pgx_gene_call_audit WHERE call_status = 'published'
                ) a ON p.gene_symbol = a.gene_symbol
                WHERE p.gene_symbol IS NULL OR a.gene_symbol IS NULL OR
                  p.diplotype IS DISTINCT FROM a.selected_diplotype OR
                  p.phenotype IS DISTINCT FROM a.selected_phenotype OR
                  p.cpic_level IS DISTINCT FROM list_extract(a.cpic_levels, 1) OR
                  p.caller != 'PharmCAT' OR p.caller_version != ? OR
                  list_count(p.affected_drugs) = 0 OR
                  list_count(list_distinct(p.affected_drugs)) != list_count(p.affected_drugs) OR
                  list_sort(p.affected_drugs) != p.affected_drugs OR
                  NOT list_has_all(a.affected_drugs, p.affected_drugs) OR
                  (p.activity_score IS NULL) != (a.activity_score IS NULL) OR
                  (p.activity_score IS NOT NULL AND
                   abs(p.activity_score::DOUBLE - a.activity_score) > 1e-5)
                """,
                [PHARMCAT_VERSION],
            ).fetchone()[0]
        )
        complex_placeholders = ",".join("?" for _ in STRUCTURALLY_COMPLEX_PGX_GENES)
        unsafe_core = int(
            connection.execute(
                f"""
                SELECT count(*) FROM pharmacogenomics
                WHERE gene_symbol IN ({complex_placeholders}) OR
                  cpic_level NOT IN ('A','B','C','D') OR diplotype IS NULL OR
                  phenotype IS NULL
                """,
                sorted(STRUCTURALLY_COMPLEX_PGX_GENES),
            ).fetchone()[0]
        )
        normative_mismatches += unsafe_core
    elif published_rows:
        normative_mismatches = published_rows
    if normative_mismatches:
        errors.append(
            f"normative pharmacogenomics/gene-audit mismatch for {normative_mismatches} row(s)"
        )

    guidance_rows = 0
    invalid_guidance = 0
    if has_guidance:
        guidance_rows, guidance_distinct, invalid_guidance = (
            int(value)
            for value in connection.execute(
                """
                SELECT count(*), count(DISTINCT guidance_id),
                  count(*) FILTER (WHERE
                    NOT regexp_full_match(guidance_id, '^pgx-guidance:[0-9a-f]{64}$') OR
                    source IS NULL OR gene_symbol IS NULL OR drug IS NULL OR
                    implications IS NULL OR matched IS NULL OR
                    (cpic_level IS NOT NULL AND cpic_level NOT IN ('A','B','C','D')) OR
                    NOT regexp_full_match(report_json_sha256, '^[0-9a-f]{64}$'))
                FROM pgx_guidance_audit
                """
            ).fetchone()
        )
        if guidance_rows != guidance_distinct:
            errors.append("PGx guidance audit contains duplicate guidance identifiers")
        matched_mismatches = int(
            connection.execute(
                """
                SELECT count(*) FROM pgx_guidance_audit g
                LEFT JOIN pgx_gene_call_audit a USING (gene_symbol)
                WHERE g.matched != coalesce(a.call_status = 'published', false)
                """
            ).fetchone()[0]
        )
        invalid_guidance += matched_mismatches
        if invalid_guidance:
            errors.append(f"PGx guidance audit contains {invalid_guidance} invalid row(s)")

    provenance = manifest.get("resource_provenance", {}).get("pharmacogenomics")
    provenance_mismatches = 0
    if not isinstance(provenance, dict):
        provenance_mismatches = 1
    else:
        expected_positions = provenance.get("pharmcat_positions_vcf_sha256")
        expected_report = provenance.get("canonical_report_sha256")
        expected_calls = provenance.get("calls_tsv_sha256")
        cpic_version = provenance.get("cpic_source_version")
        cpic_release = provenance.get("cpic_source_release")
        cpic_api = provenance.get("cpic_api_identity")
        release_valid = (
            isinstance(cpic_version, str)
            and isinstance(cpic_release, dict)
            and cpic_release.get("tag") == f"v{cpic_version}"
            and isinstance(cpic_release.get("commit"), str)
            and _SHA1_RE.fullmatch(cpic_release["commit"]) is not None
            and isinstance(cpic_release.get("github_release_id"), int)
            and isinstance(cpic_release.get("published_at"), str)
            and bool(cpic_release.get("published_at"))
            and cpic_release.get("release_url")
            == f"https://github.com/cpicpgx/cpic-data/releases/tag/v{cpic_version}"
        )
        api_valid = (
            isinstance(cpic_api, dict)
            and _is_sha256(cpic_api.get("openapi_canonical_sha256"))
            and isinstance(cpic_api.get("postgrest_version"), str)
            and bool(cpic_api.get("postgrest_version"))
            and cpic_api.get("swagger_version") == "2.0"
            and cpic_api.get("base_path") == "/v1/"
        )
        provenance_mismatches += int(
            provenance.get("pharmcat_version") != PHARMCAT_VERSION
            or not isinstance(provenance.get("pharmcat_release_commit"), str)
            or _SHA1_RE.fullmatch(provenance.get("pharmcat_release_commit", "")) is None
            or not _is_sha256(expected_positions)
            or provenance.get("pharmcat_positions_records") != PHARMCAT_POSITIONS_RECORDS
            or provenance.get("pharmcat_positions_alleles") != PHARMCAT_POSITIONS_ALLELES
            or not isinstance(provenance.get("pharmcat_normalization_changed_alleles"), int)
            or provenance.get("pharmcat_normalization_changed_alleles", -1) < 0
            or not isinstance(provenance.get("pharmcat_container_image"), str)
            or "@sha256:" not in provenance.get("pharmcat_container_image", "")
            or not isinstance(provenance.get("pharmcat_license"), str)
            or not provenance.get("pharmcat_license")
            or not isinstance(provenance.get("pharmcat_redistribution_allowed"), bool)
            or not isinstance(provenance.get("cpic_license"), str)
            or not provenance.get("cpic_license")
            or not isinstance(provenance.get("cpic_redistribution_allowed"), bool)
            or provenance.get("cpic_source") != "CPIC PostgREST pair_view"
            or not isinstance(provenance.get("cpic_snapshot_date"), str)
            or not _is_sha256(provenance.get("cpic_snapshot_sha256"))
            or not _is_sha256(provenance.get("cpic_rowset_canonical_sha256"))
            or not isinstance(provenance.get("cpic_normative_recommendation_rows"), int)
            or provenance.get("cpic_normative_recommendation_rows", 0) <= 0
            or not release_valid
            or not api_valid
            or not _is_sha256(expected_report)
            or not _is_sha256(expected_calls)
            or provenance.get("missing_to_reference_conversion") is not False
            or provenance.get("network_during_matching") != "none"
            or provenance.get("complex_matcher_genes_excluded")
            != sorted(STRUCTURALLY_COMPLEX_PGX_GENES)
        )
        provenance_mismatches += int(
            connection.execute(
                "SELECT count(*) FROM pgx_locus_audit WHERE positions_vcf_sha256 != ?",
                [expected_positions],
            ).fetchone()[0]
        )
        provenance_mismatches += int(
            connection.execute(
                """
                SELECT count(*) FROM pgx_gene_call_audit
                WHERE report_json_sha256 != ? OR calls_tsv_sha256 != ?
                """,
                [expected_report, expected_calls],
            ).fetchone()[0]
        )
        if has_guidance:
            provenance_mismatches += int(
                connection.execute(
                    "SELECT count(*) FROM pgx_guidance_audit WHERE report_json_sha256 != ?",
                    [expected_report],
                ).fetchone()[0]
            )
    if provenance_mismatches:
        errors.append(f"PGx manifest provenance mismatch ({provenance_mismatches})")

    return {
        "locus_rows": locus_total,
        "allele_rows": allele_total,
        "gene_rows": gene_total,
        "guidance_rows": guidance_rows,
        "published_rows": published_rows,
        "invalid_locus_rows": invalid_loci,
        "target_proof_mismatches": target_mismatches,
        "orphan_variant_ids": orphan_variants,
        "recomputation_mismatches": recomputation_mismatches,
        "invalid_gene_rows": invalid_genes,
        "invalid_guidance_rows": invalid_guidance,
        "normative_mismatches": normative_mismatches,
        "provenance_mismatches": provenance_mismatches,
    }


def _reference_contig(output_chrom: str, references: set[str]) -> str | None:
    token = output_chrom[3:]
    candidates = [output_chrom, token]
    if token == "M":
        candidates.extend(["MT", "chrMT"])
    return next((candidate for candidate in candidates if candidate in references), None)


def _validate_ref_and_partition_order(
    variants_path: Path,
    reference_fasta: Path | None,
    errors: list[str],
) -> int:
    fasta = pysam.FastaFile(str(reference_fasta)) if reference_fasta else None
    checked = 0
    try:
        references = set(fasta.references) if fasta else set()
        for file_path in _parquet_files(variants_path):
            partition = next(
                (
                    part.split("=", maxsplit=1)[1]
                    for part in file_path.relative_to(variants_path).parts
                    if part.startswith("chrom=")
                ),
                None,
            )
            previous_pos = 0
            sequence: str | None = None
            reference_name: str | None = None
            parquet = pq.ParquetFile(file_path)
            for batch in parquet.iter_batches(batch_size=100_000, columns=["chrom", "pos", "ref"]):
                values = batch.to_pydict()
                for chrom, pos, ref in zip(
                    values["chrom"], values["pos"], values["ref"], strict=True
                ):
                    checked += 1
                    if partition is not None and chrom != partition:
                        errors.append(f"partition/value chromosome mismatch in {file_path}")
                        return checked
                    if pos < previous_pos:
                        errors.append(f"variant positions are not sorted in {file_path}")
                        return checked
                    previous_pos = pos
                    if fasta is not None:
                        if reference_name is None:
                            reference_name = _reference_contig(chrom, references)
                            if reference_name is None:
                                errors.append(
                                    f"bundle chromosome absent from validation FASTA: {chrom}"
                                )
                                return checked
                            sequence = fasta.fetch(reference_name).upper()
                        assert sequence is not None
                        observed = sequence[pos - 1 : pos - 1 + len(ref)]
                        if observed != ref:
                            errors.append(
                                f"bundle REF mismatch at {chrom}:{pos}: {ref} != {observed}"
                            )
                            return checked
    finally:
        if fasta is not None:
            fasta.close()
    return checked


def _validate_family_evidence(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, Any]:
    provenance = manifest.get("resource_provenance", {}).get("family_evidence")
    if "family_evidence" not in views:
        if provenance is not None:
            errors.append("family evidence provenance exists without its Parquet extension")
        return {"status": "not_present", "rows": 0}
    if not isinstance(provenance, dict):
        errors.append("family evidence extension lacks resource provenance")
        return {"status": "invalid_provenance", "rows": 0}
    for key in ("joint_vcf_sha256", "pedigree_sha256"):
        if not _is_sha256(provenance.get(key)):
            errors.append(f"family evidence provenance has invalid {key}")
    for key in (
        "relative_identifiers_persisted",
        "relative_genotypes_persisted",
        "genotype_refinements",
    ):
        if provenance.get(key) != 0:
            errors.append(f"family evidence privacy invariant failed: {key}")
    projection = provenance.get("normalization_projection")
    projection_counts: tuple[int, int, int] | None = None
    if not isinstance(projection, dict):
        errors.append("family evidence lacks normalization-projection provenance")
    else:
        count_fields = (
            "source_alleles",
            "output_alleles",
            "extra_output_alleles",
            "family_source_records_projected",
            "incomplete_subject_genotype_records_omitted",
            "incomplete_subject_genotype_records_with_alt_evidence",
        )
        if any(
            not isinstance(projection.get(key), int)
            or isinstance(projection.get(key), bool)
            or projection[key] < 0
            for key in count_fields
        ):
            errors.append("family normalization-projection counts are invalid")
        else:
            source_count = projection["source_alleles"]
            output_count = projection["output_alleles"]
            extra_count = projection["extra_output_alleles"]
            projection_counts = source_count, output_count, extra_count
            if (
                output_count < source_count
                or extra_count != output_count - source_count
                or projection["family_source_records_projected"] != source_count
                or projection["incomplete_subject_genotype_records_with_alt_evidence"]
                > projection["incomplete_subject_genotype_records_omitted"]
            ):
                errors.append("family normalization-projection accounting is inconsistent")
        if not _is_sha256(projection.get("variant_provenance_sha256")):
            errors.append("family normalization projection has an invalid provenance digest")
        if projection.get("incomplete_subject_genotype_policy") != (
            "omit from family variant evidence as no-call; never convert partial GT to "
            "reference or a complete alternate genotype"
        ):
            errors.append("family incomplete-genotype policy is invalid")
        if projection.get("changed_allele_policy") != (
            "project family and LOCO annotations through exact variant provenance while "
            "preserving source and projected genotype/phasing"
        ):
            errors.append("family changed-allele projection policy is invalid")
    if projection_counts is not None and "variant_provenance" in views:
        observed_projection = connection.execute(
            """
            SELECT count(DISTINCT source_id || ':' || CAST(source_alt_index AS VARCHAR)),
                   count(*)
            FROM variant_provenance
            WHERE normalization_changed
            """
        ).fetchone()
        assert observed_projection is not None
        observed_source = int(observed_projection[0])
        observed_output = int(observed_projection[1])
        if projection_counts != (
            observed_source,
            observed_output,
            observed_output - observed_source,
        ):
            errors.append("family normalization projection differs from variant provenance")

    result = connection.execute(
        """
        SELECT
          count(*) AS rows,
          count(DISTINCT evidence_id) AS evidence_ids,
          count(DISTINCT variant_id) AS variant_ids,
          count(DISTINCT joint_vcf_sha256) AS joint_digests,
          count(DISTINCT pedigree_sha256) AS pedigree_digests,
          count(*) FILTER (
            WHERE NOT regexp_full_match(evidence_id, 'family-evidence:[0-9a-f]{64}')
               OR NOT regexp_full_match(subject_evidence_digest, '[0-9a-f]{64}')
               OR NOT regexp_full_match(joint_vcf_sha256, '[0-9a-f]{64}')
               OR NOT regexp_full_match(pedigree_sha256, '[0-9a-f]{64}')
          ) AS invalid_identifiers,
          count(*) FILTER (
            WHERE list_count(subject_gt) NOT IN (1, 2)
               OR NOT list_contains(subject_gt, 1)
               OR list_count(list_filter(subject_gt, allele -> allele IS NULL OR allele NOT IN (0, 1))) > 0
          ) AS invalid_genotypes,
          count(*) FILTER (
            WHERE (subject_phased AND phase_set IS NOT NULL
                   AND phase_confidence != 'phased_with_phase_set')
               OR (subject_phased AND phase_set IS NULL
                   AND phase_confidence != 'phased_without_phase_set')
               OR (NOT subject_phased AND phase_confidence != 'unphased')
          ) AS invalid_phase_semantics,
          count(*) FILTER (
            WHERE complex_region != (list_count(difficult_region_flags) > 0)
          ) AS invalid_complex_semantics,
          count(*) FILTER (
            WHERE loco_status IN (
              'not_individually_persisted',
              'not_assessed_outside_autosomal_loco',
              'not_configured'
            ) AND (
              loco_source_status IS NOT NULL OR
              loco_subject_gt_concordant IS NOT NULL OR
              loco_exact_allele_in_model IS NOT NULL
            ) OR loco_status = 'source_genotype_discordant' AND (
              loco_source_status NOT IN (
                'uninformative_parent_not_callable',
                'de_novo_candidate',
                'mendelian_conflict'
              ) OR loco_subject_gt_concordant IS DISTINCT FROM false OR
              loco_exact_allele_in_model IS NULL
            ) OR loco_status IN (
              'uninformative_parent_not_callable',
              'de_novo_candidate',
              'mendelian_conflict'
            ) AND (
              loco_source_status != loco_status OR
              loco_subject_gt_concordant IS DISTINCT FROM true OR
              loco_exact_allele_in_model IS NULL
            ) OR loco_status NOT IN (
              'not_individually_persisted',
              'not_assessed_outside_autosomal_loco',
              'not_configured',
              'source_genotype_discordant',
              'uninformative_parent_not_callable',
              'de_novo_candidate',
              'mendelian_conflict'
            )
          ) AS invalid_loco_semantics,
          count(*) FILTER (
            WHERE (loco_review_sha256 IS NULL) != (loco_summary_sha256 IS NULL)
               OR (loco_review_sha256 IS NOT NULL AND
                   (NOT regexp_full_match(loco_review_sha256, '[0-9a-f]{64}') OR
                    NOT regexp_full_match(loco_summary_sha256, '[0-9a-f]{64}')))
          ) AS invalid_loco_digests,
          count(*) FILTER (
            WHERE de_novo_candidate AND (
              inheritance_state != 'de_novo_candidate' OR mendelian_status != 'conflict'
            )
          ) AS invalid_de_novo_semantics
        FROM family_evidence
        """
    ).fetchone()
    assert result is not None
    names = [description[0] for description in connection.description]
    metrics = dict(zip(names, result, strict=True))
    if metrics["rows"] != metrics["evidence_ids"] or metrics["rows"] != metrics["variant_ids"]:
        errors.append("family evidence identifiers are not one-to-one")
    expected_digest_count = 1 if metrics["rows"] else 0
    if (
        metrics["joint_digests"] != expected_digest_count
        or metrics["pedigree_digests"] != expected_digest_count
    ):
        errors.append("family evidence does not have one immutable source boundary")
    for key in (
        "invalid_identifiers",
        "invalid_genotypes",
        "invalid_phase_semantics",
        "invalid_complex_semantics",
        "invalid_loco_semantics",
        "invalid_loco_digests",
        "invalid_de_novo_semantics",
    ):
        if metrics[key]:
            errors.append(f"family evidence has {metrics[key]} {key.replace('_', ' ')} row(s)")
    boundary = connection.execute(
        "SELECT min(joint_vcf_sha256), min(pedigree_sha256), min(source_version), "
        "min(CAST(snapshot_date AS VARCHAR)) FROM family_evidence"
    ).fetchone()
    if boundary is not None and metrics["rows"]:
        if boundary[0] != provenance.get("joint_vcf_sha256"):
            errors.append("family joint VCF digest differs from manifest provenance")
        if boundary[1] != provenance.get("pedigree_sha256"):
            errors.append("family PED digest differs from manifest provenance")
        if boundary[2] != provenance.get("source_version"):
            errors.append("family source version differs from manifest provenance")
        if boundary[3] != provenance.get("snapshot_date"):
            errors.append("family snapshot date differs from manifest provenance")

    reconciliation = connection.execute(
        """
        SELECT
          count(*) FILTER (WHERE v.variant_id IS NULL) AS orphan_family_rows,
          count(*) FILTER (WHERE f.variant_id IS NULL) AS variants_without_family_rows,
          count(*) FILTER (
            WHERE f.variant_id IS NOT NULL AND v.variant_id IS NOT NULL AND (
              f.chrom != v.chrom OR f.pos != v.pos OR f.ref != v.ref OR f.alt != v.alt OR
              f.subject_gt != v.genotype.gt OR f.subject_phased != v.genotype.phased
            )
          ) AS subject_ledger_mismatches
        FROM family_evidence f
        FULL OUTER JOIN variants v USING (variant_id)
        """
    ).fetchone()
    assert reconciliation is not None
    metrics.update(
        {
            "orphan_family_rows": int(reconciliation[0]),
            "variants_without_family_rows": int(reconciliation[1]),
            "subject_ledger_mismatches": int(reconciliation[2]),
        }
    )
    for key in (
        "orphan_family_rows",
        "variants_without_family_rows",
        "subject_ledger_mismatches",
    ):
        if metrics[key]:
            errors.append(f"family evidence reconciliation failed: {key}={metrics[key]}")
    expected_rows = provenance.get("family_evidence_rows")
    if not isinstance(expected_rows, int) or expected_rows != metrics["rows"]:
        errors.append("family evidence row count differs from manifest provenance")
    metrics["status"] = "validated"
    return metrics


def _validate_gwas_catalog(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, Any]:
    present = {name for name in ("gwas_associations", "gwas_catalog_evidence") if name in views}
    if not present:
        return {"status": "not_present", "association_rows": 0, "evidence_rows": 0}
    if present != {"gwas_associations", "gwas_catalog_evidence"}:
        errors.append("GWAS compact and lossless evidence tables must be present together")
        return {"status": "invalid_partial_boundary"}

    evidence_rows, evidence_ids, exact_failures, dosage_failures = connection.execute(
        """
        SELECT count(*), count(DISTINCT e.association_id),
               count(*) FILTER (
                 WHERE v.variant_id IS NULL OR e.effect_allele NOT IN (v.ref, v.alt)
               ),
               count(*) FILTER (
                 WHERE v.variant_id IS NOT NULL AND e.effect_allele IN (v.ref, v.alt) AND
                   e.effect_allele_dosage != CASE
                     WHEN e.effect_allele = v.ref THEN
                       list_count(list_filter(v.genotype.gt, allele -> allele = 0))
                     ELSE list_count(list_filter(v.genotype.gt, allele -> allele = 1))
                   END
               )
        FROM gwas_catalog_evidence e
        LEFT JOIN variants v USING (variant_id)
        """
    ).fetchone()
    if evidence_rows != evidence_ids:
        errors.append("GWAS lossless evidence contains duplicate association_id")
    if exact_failures:
        errors.append(f"GWAS evidence has {exact_failures} non-exact variant/allele joins")
    if dosage_failures:
        errors.append(f"GWAS evidence has {dosage_failures} effect-allele dosage mismatches")

    association_rows, association_ids, malformed = connection.execute(
        """
        SELECT count(*), count(DISTINCT association_id),
               count(*) FILTER (
                 WHERE effect_size IS NULL OR effect_type NOT IN
                   ('beta','odds_ratio','hazard_ratio') OR
                   effect_allele IS NULL OR length(effect_allele) = 0 OR
                   len(study_pmids) = 0
               )
        FROM gwas_associations
        """
    ).fetchone()
    if association_rows != association_ids:
        errors.append("gwas_associations contains duplicate association_id")
    if malformed:
        errors.append(f"gwas_associations contains {malformed} malformed normative rows")

    missing_audit, orphan_audit, value_mismatch = connection.execute(
        """
        SELECT
          count(*) FILTER (WHERE e.association_id IS NULL),
          count(*) FILTER (WHERE a.association_id IS NULL AND e.normative_row_emitted),
          count(*) FILTER (
            WHERE a.association_id IS NOT NULL AND e.association_id IS NOT NULL AND
              (a.variant_id != e.variant_id OR a.effect_allele != e.effect_allele OR
               a.effect_size != e.effect_size OR a.effect_type != e.effect_type OR
               a.source_version != e.source_version)
          )
        FROM gwas_associations a
        FULL OUTER JOIN gwas_catalog_evidence e USING (association_id)
        """
    ).fetchone()
    if missing_audit:
        errors.append(f"{missing_audit} GWAS normative rows lack lossless audit evidence")
    if orphan_audit:
        errors.append(f"{orphan_audit} typed GWAS audit rows lack normative rows")
    if value_mismatch:
        errors.append(f"{value_mismatch} GWAS compact/audit rows disagree")

    flag_mismatches = connection.execute(
        """
        WITH matched AS (
          SELECT variant_id, count(*) AS n,
                 coalesce(
                   min(p_value) FILTER (WHERE normative_row_emitted), min(p_value)
                 ) AS selected_p
          FROM gwas_catalog_evidence GROUP BY variant_id
        )
        SELECT count(*)
        FROM variants v LEFT JOIN matched m USING (variant_id)
        WHERE v.trait_associations.is_gwas_hit != (coalesce(m.n, 0) > 0) OR
              (m.n > 0 AND abs(v.trait_associations.p_value - m.selected_p) > 1e-7)
        """
    ).fetchone()[0]
    if flag_mismatches:
        errors.append(f"{flag_mismatches} variant-level GWAS summaries disagree with evidence")

    provenance = manifest.get("resource_provenance", {}).get("gwas_catalog")
    if not isinstance(provenance, dict):
        errors.append("GWAS tables exist without manifest resource provenance")
    else:
        for field in ("archive_sha256", "normalized_parquet_sha256"):
            if not _is_sha256(provenance.get(field)):
                errors.append(f"manifest GWAS provenance has invalid {field}")
        if provenance.get("source_archive_redistributed") is not False:
            errors.append("manifest must state that the raw GWAS source archive is not embedded")

    return {
        "status": "validated",
        "association_rows": int(association_rows),
        "evidence_rows": int(evidence_rows),
        "exact_join_failures": int(exact_failures),
        "dosage_failures": int(dosage_failures),
        "variant_summary_failures": int(flag_mismatches),
    }


def _validate_clingen_gene_context(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    manifest: dict[str, Any],
    documents: dict[str, dict[str, Any] | None],
    root: Path,
    errors: list[str],
) -> dict[str, Any]:
    """Validate ClinGen as immutable gene-level context, never subject evidence."""

    table_names = ("clingen_gene_disease", "clingen_dosage_sensitivity")
    present = {name for name in table_names if name in views}
    provenance = manifest.get("resource_provenance", {}).get("clingen_gene_level")
    annotation_source = manifest.get("annotation_sources", {}).get("clingen_gene_level")

    raw_document = documents.get("raw_artifacts.json")
    transformation_document = documents.get("transformations.json")
    omission_document = documents.get("omissions.json")
    raw_rows = raw_document.get("artifacts", []) if isinstance(raw_document, dict) else []
    stages = (
        transformation_document.get("stages", [])
        if isinstance(transformation_document, dict)
        else []
    )
    omissions = (
        omission_document.get("omissions", []) if isinstance(omission_document, dict) else []
    )
    raw_matches = [
        row
        for row in raw_rows
        if isinstance(row, dict)
        and row.get("artifact_type") == "ClinGen gene-level curation snapshot"
    ]
    stage_matches = [
        row
        for row in stages
        if isinstance(row, dict) and row.get("stage_id") == "clingen_gene_context_projection"
    ]
    disabled_omissions = [
        row
        for row in omissions
        if isinstance(row, dict)
        and row.get("scope") == "ClinGen gene-disease validity and dosage-sensitivity context"
    ]
    table_omission_scopes = {
        "clingen_gene_disease": "ClinGen gene-disease projection table",
        "clingen_dosage_sensitivity": "ClinGen dosage-sensitivity projection table",
    }

    if provenance is None:
        if present:
            errors.append("ClinGen tables exist without manifest resource provenance")
        if raw_matches:
            errors.append("disabled ClinGen adapter retains a raw-artifact claim")
        if stage_matches:
            errors.append("disabled ClinGen adapter retains a transformation stage")
        if len(disabled_omissions) != 1:
            errors.append("disabled ClinGen adapter must have exactly one explicit omission")
        if annotation_source != "not_run":
            errors.append("disabled ClinGen annotation source must be recorded as not_run")
        return {
            "status": "not_present",
            "gene_disease_rows": 0,
            "dosage_rows": 0,
        }

    if not isinstance(provenance, dict):
        errors.append("ClinGen manifest provenance must be an object")
        return {"status": "invalid_provenance"}
    if len(raw_matches) != 1:
        errors.append("enabled ClinGen adapter must have exactly one raw-artifact boundary")
    if len(stage_matches) != 1:
        errors.append("enabled ClinGen adapter must have exactly one transformation stage")
    if disabled_omissions:
        errors.append("enabled ClinGen adapter retains its disabled-state omission")
    if not isinstance(annotation_source, str) or annotation_source == "not_run":
        errors.append("enabled ClinGen annotation source is absent from manifest")

    expected_boundary = {
        "snapshot": CLINGEN_GENE_LEVEL_SNAPSHOT,
        "canonical_sha256": CLINGEN_GENE_LEVEL_CANONICAL_SHA256,
        "license": CLINGEN_LICENSE,
        "redistribution_allowed": True,
        "semantic_boundary": CLINGEN_SEMANTIC_BOUNDARY,
        "variant_pathogenicity_claims": 0,
        "diagnosis_claims": 0,
    }
    for key, expected in expected_boundary.items():
        if provenance.get(key) != expected:
            errors.append(f"ClinGen provenance boundary mismatch: {key}")
    if not _is_sha256(provenance.get("manifest_sha256")):
        errors.append("ClinGen provenance has an invalid resource-manifest digest")

    sources = provenance.get("sources")
    if not isinstance(sources, dict) or set(sources) != {
        "dosage_sensitivity",
        "gene_disease_validity",
    }:
        errors.append("ClinGen provenance has an invalid source set")
        sources = {}
    for name, metadata in sources.items():
        if not isinstance(metadata, dict) or not _is_sha256(metadata.get("sha256")):
            errors.append(f"ClinGen source {name!r} has an invalid digest")

    resource_outputs = provenance.get("resource_outputs")
    expected_resource_counts = {
        "gene_disease_validity": CLINGEN_GENE_DISEASE_ROWS,
        "dosage_sensitivity": CLINGEN_DOSAGE_ROWS,
    }
    if not isinstance(resource_outputs, dict) or set(resource_outputs) != set(
        expected_resource_counts
    ):
        errors.append("ClinGen provenance has an invalid normalized-resource output set")
    else:
        for name, expected_rows in expected_resource_counts.items():
            metadata = resource_outputs.get(name)
            if (
                not isinstance(metadata, dict)
                or metadata.get("rows") != expected_rows
                or not _is_sha256(metadata.get("sha256"))
            ):
                errors.append(f"ClinGen normalized-resource identity mismatch: {name}")

    projection = provenance.get("bundle_projection")
    if not isinstance(projection, dict):
        errors.append("ClinGen provenance lacks bundle projection accounting")
        projection = {}
    projection_counts = {
        key: projection.get(key)
        for key in (
            "input_gene_symbols",
            "matched_gene_symbols",
            "unmatched_gene_symbols",
            "gene_disease_rows",
            "dosage_rows",
        )
    }
    if any(not isinstance(value, int) or value < 0 for value in projection_counts.values()):
        errors.append("ClinGen projection counts are not non-negative integers")
    elif (
        projection_counts["matched_gene_symbols"] + projection_counts["unmatched_gene_symbols"]
        != projection_counts["input_gene_symbols"]
    ):
        errors.append("ClinGen matched/unmatched gene-symbol accounting does not reconcile")

    actual_outputs: dict[str, dict[str, Any]] = {}
    gene_rows = 0
    dosage_rows = 0
    if "clingen_gene_disease" in views:
        path = root / "extensions" / "clingen_gene_disease.parquet"
        gene_rows, distinct_evidence, distinct_assertions, invalid = connection.execute(
            """
            SELECT count(*), count(DISTINCT evidence_id), count(DISTINCT assertion_id),
                   count(*) FILTER (
                     WHERE evidence_id != 'clingen-gene-disease:' || assertion_id OR
                           gene_symbol = '' OR
                           NOT regexp_full_match(hgnc_id, 'HGNC:[1-9][0-9]*') OR
                           NOT regexp_full_match(mondo_id, 'MONDO:[0-9]+') OR
                           classification NOT IN (
                             'Definitive','Strong','Moderate','Limited','Disputed','Refuted',
                             'No Known Disease Relationship'
                           ) OR
                           NOT starts_with(online_report,
                             'https://search.clinicalgenome.org/kb/gene-validity/') OR
                           snapshot != ? OR source_row_number < 1 OR
                           NOT regexp_full_match(source_row_sha256, '[0-9a-f]{64}') OR
                           classification_datetime_source = '' OR gcep = ''
                   )
            FROM clingen_gene_disease
            """,
            [CLINGEN_GENE_LEVEL_SNAPSHOT],
        ).fetchone()
        gene_rows = int(gene_rows)
        if gene_rows != distinct_evidence or gene_rows != distinct_assertions:
            errors.append("ClinGen gene-disease identifiers are not unique")
        if invalid:
            errors.append(f"ClinGen gene-disease table has {invalid} malformed row(s)")
        unknown_classifications = connection.execute(
            "SELECT DISTINCT classification FROM clingen_gene_disease"
        ).fetchall()
        if {row[0] for row in unknown_classifications} - set(GENE_DISEASE_CLASSIFICATIONS):
            errors.append("ClinGen gene-disease table contains an unknown classification")
        actual_outputs[path.name] = {"rows": gene_rows, "sha256": sha256_file(path)}

    if "clingen_dosage_sensitivity" in views:
        path = root / "extensions" / "clingen_dosage_sensitivity.parquet"
        dosage_rows, distinct_evidence, distinct_genes, invalid = connection.execute(
            """
            SELECT count(*), count(DISTINCT evidence_id), count(DISTINCT gene_symbol),
                   count(*) FILTER (
                     WHERE evidence_id != 'clingen-dosage:' || gene_symbol OR
                           gene_symbol = '' OR ncbi_gene_id < 1 OR
                           NOT regexp_full_match(chrom, '([1-9]|1[0-9]|2[0-2]|X|Y|M)') OR
                           start < 1 OR "end" < start OR
                           snapshot != ? OR source_row_number < 1 OR
                           NOT regexp_full_match(source_row_sha256, '[0-9a-f]{64}') OR
                           list_count(list_filter(haploinsufficiency_pmids, value -> value < 1)) > 0 OR
                           list_count(list_filter(triplosensitivity_pmids, value -> value < 1)) > 0
                   )
            FROM clingen_dosage_sensitivity
            """,
            [CLINGEN_GENE_LEVEL_SNAPSHOT],
        ).fetchone()
        dosage_rows = int(dosage_rows)
        if dosage_rows != distinct_evidence or dosage_rows != distinct_genes:
            errors.append("ClinGen dosage identifiers and gene symbols are not unique")
        if invalid:
            errors.append(f"ClinGen dosage table has {invalid} malformed row(s)")
        score_rows = connection.execute(
            """
            SELECT haploinsufficiency_score, haploinsufficiency_source_score,
                   haploinsufficiency_description, triplosensitivity_score,
                   triplosensitivity_source_score, triplosensitivity_description
            FROM clingen_dosage_sensitivity
            """
        ).fetchall()
        invalid_scores = 0
        for row in score_rows:
            for score, source_score, description in (row[:3], row[3:]):
                if score is None:
                    valid = source_score == "Not yet evaluated" and description == source_score
                else:
                    valid = (
                        source_score == str(score)
                        and DOSAGE_SCORE_DESCRIPTIONS.get(score) == description
                    )
                invalid_scores += not valid
        if invalid_scores:
            errors.append(f"ClinGen dosage table has {invalid_scores} invalid score semantic(s)")
        actual_outputs[path.name] = {"rows": dosage_rows, "sha256": sha256_file(path)}

    if present:
        if "gene_index" not in views:
            errors.append("ClinGen projection tables require gene_index.parquet")
        else:
            for table in sorted(present):
                orphan_rows = connection.execute(
                    f"SELECT count(*) FROM {table} c WHERE NOT EXISTS "
                    "(SELECT 1 FROM gene_index g WHERE g.gene_symbol = c.gene_symbol)"
                ).fetchone()[0]
                if orphan_rows:
                    errors.append(
                        f"{table} has {orphan_rows} gene symbol(s) absent from gene_index"
                    )

    matched_symbols = 0
    if present:
        union_sql = " UNION ".join(f"SELECT gene_symbol FROM {name}" for name in sorted(present))
        matched_symbols = int(
            connection.execute(f"SELECT count(*) FROM ({union_sql})").fetchone()[0]
        )
    if projection_counts.get("gene_disease_rows") != gene_rows:
        errors.append("ClinGen gene-disease row count differs from manifest provenance")
    if projection_counts.get("dosage_rows") != dosage_rows:
        errors.append("ClinGen dosage row count differs from manifest provenance")
    if projection_counts.get("matched_gene_symbols") != matched_symbols:
        errors.append("ClinGen matched gene-symbol count differs from projected tables")
    if projection.get("outputs") != actual_outputs:
        errors.append("ClinGen projected output identities differ from manifest provenance")

    for table, scope in table_omission_scopes.items():
        matches = [row for row in omissions if isinstance(row, dict) and row.get("scope") == scope]
        if table in present and matches:
            errors.append(f"present {table} retains an absent-table omission")
        if table not in present and len(matches) != 1:
            errors.append(f"absent {table} must have exactly one explicit omission")

    if raw_matches:
        raw = raw_matches[0]
        source_digests = ",".join(
            f"{name}=sha256:{metadata.get('sha256')}"
            for name, metadata in sorted(sources.items())
            if isinstance(metadata, dict)
        )
        expected_integrity = (
            "canonical-sha256:"
            + str(provenance.get("canonical_sha256"))
            + ";manifest-sha256:"
            + str(provenance.get("manifest_sha256"))
            + ";sources:"
            + source_digests
        )
        if raw.get("integrity_identifier") != expected_integrity:
            errors.append("ClinGen raw-artifact integrity boundary differs from provenance")
        if raw.get("license") != CLINGEN_LICENSE or raw.get("redistribution_allowed") is not True:
            errors.append("ClinGen raw-artifact licence boundary is invalid")
        coverage = raw.get("coverage")
        if not isinstance(coverage, dict) or coverage.get("evidence_level") != "gene":
            errors.append("ClinGen raw artifact does not declare gene-level scope")

    if stage_matches:
        stage = stage_matches[0]
        expected_stage_outputs = sorted("extensions/" + name for name in actual_outputs)
        for key, expected in {
            "status": "completed",
            "matching_policy": "exact case-sensitive HGNC gene symbol",
            "semantic_boundary": CLINGEN_SEMANTIC_BOUNDARY,
            "variant_pathogenicity_claims_created": 0,
            "subject_diagnosis_claims_created": 0,
            "gene_disease_rows": gene_rows,
            "dosage_rows": dosage_rows,
            "matched_gene_symbols": matched_symbols,
            "outputs": expected_stage_outputs,
        }.items():
            if stage.get(key) != expected:
                errors.append(f"ClinGen transformation-stage mismatch: {key}")

    if "clinical_evidence" in views:
        leaked = connection.execute(
            "SELECT count(*) FROM clinical_evidence WHERE lower(source) LIKE '%clingen%'"
        ).fetchone()[0]
        if leaked:
            errors.append(
                f"clinical_evidence contains {leaked} ClinGen gene-context row(s) as variant claims"
            )

    return {
        "status": "validated",
        "gene_disease_rows": gene_rows,
        "dosage_rows": dosage_rows,
        "matched_gene_symbols": matched_symbols,
        "variant_pathogenicity_claims": 0,
        "diagnosis_claims": 0,
    }


def _register_views(connection: duckdb.DuckDBPyConnection, root: Path) -> set[str]:
    paths = {
        "variants": root / "variants.parquet",
        "gene_index": root / "gene_index.parquet",
        "clinical_findings": root / "clinical_findings.parquet",
        "clinical_evidence": root / "clinical_evidence.parquet",
        "interpretations": root / "interpretations.parquet",
        "callability": root / "callability.parquet",
        "prs": root / "prs.parquet",
        "pharmacogenomics": root / "pharmacogenomics.parquet",
        "clinvar_assertions": root / "extensions" / "clinvar_assertions.parquet",
        "acmg_secondary_findings_evidence": (
            root / "extensions" / "acmg_secondary_findings_evidence.parquet"
        ),
        "callable_regions": root / "callable_regions.parquet",
        "variant_provenance": root / "extensions" / "variant_provenance.parquet",
        "reference_evidence": root / "extensions" / "reference_confidence_regions.parquet",
        "non_core_variants": root / "extensions" / "non_core_variants.parquet",
        "observed_partial_variants": root / "extensions" / "observed_partial_variants.parquet",
        "target_callability_evidence": (
            root / "extensions" / "target_callability_evidence.parquet"
        ),
        "pgs_model_audit": root / "extensions" / "pgs_model_audit.parquet",
        "pgs_locus_audit": root / "extensions" / "pgs_locus_audit.parquet",
        "pgx_locus_audit": root / "extensions" / "pgx_locus_audit.parquet",
        "pgx_gene_call_audit": root / "extensions" / "pgx_gene_call_audit.parquet",
        "pgx_guidance_audit": root / "extensions" / "pgx_guidance_audit.parquet",
        "family_evidence": root / "extensions" / "family_evidence.parquet",
        "gwas_associations": root / "gwas_associations.parquet",
        "gwas_catalog_evidence": root / "extensions" / "gwas_catalog_evidence.parquet",
        "clingen_gene_disease": root / "extensions" / "clingen_gene_disease.parquet",
        "clingen_dosage_sensitivity": (root / "extensions" / "clingen_dosage_sensitivity.parquet"),
        **{
            domain: root / "extensions" / "release_c" / f"{domain}.parquet"
            for domain in RELEASE_C_DOMAINS
        },
    }
    views: set[str] = set()
    for name, path in paths.items():
        if not path.exists():
            continue
        parquet_path = path / "**" / "*.parquet" if path.is_dir() else path
        escaped = str(parquet_path).replace("'", "''")
        connection.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{escaped}')")
        views.add(name)
    return views


def _validate_profile(root: Path, manifest: dict[str, Any], errors: list[str]) -> int:
    if not manifest.get("profiles", {}).get("agent_exploration"):
        return 0
    for relative in PROFILE_PATHS:
        if not (root / relative).exists():
            errors.append(f"Agent Exploration Profile missing path: {relative}")
    path = root / "readmygenome" / "evals" / "canonical_questions.json"
    try:
        evaluation = json.loads(path.read_text(encoding="utf-8"))
        cases = evaluation["cases"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        errors.append(f"invalid Agent Profile evaluation suite: {exc}")
        return 0
    if not isinstance(cases, list) or len(cases) < 30:
        errors.append("Agent Profile evaluation suite has fewer than 30 cases")
        return len(cases) if isinstance(cases, list) else 0
    required = {
        "case_id",
        "question",
        "expected_tables_searched",
        "required_findings",
        "claims_that_must_not_be_made",
        "expected_callability_behaviour",
        "expected_citation_fields",
    }
    for case in cases:
        if not isinstance(case, dict) or not required.issubset(case):
            errors.append("Agent Profile evaluation case is missing required semantics")
            break
    return len(cases)


def validate_bundle(
    bundle: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    expected_schema_path: str | Path | None = None,
    run_profile_queries: bool = False,
) -> dict[str, Any]:
    root = Path(bundle).expanduser().resolve()
    errors: list[str] = []
    warnings: list[str] = []
    metrics: dict[str, Any] = {}
    if not root.is_dir():
        raise BundleValidationError(f"bundle directory not found: {root}")
    for item in root.rglob("*"):
        if item.is_symlink() or (item.exists() and not item.is_file() and not item.is_dir()):
            errors.append(f"forbidden bundle entry: {item.relative_to(root)}")

    schema_path = root / "schema.json"
    manifest_path = root / "manifest.json"
    variants_path = root / "variants.parquet"
    if not schema_path.is_file() or not manifest_path.is_file() or not variants_path.is_dir():
        raise BundleValidationError(
            "bundle is missing required schema.json, manifest.json, or variants.parquet/"
        )
    schema_digest = sha256_file(schema_path)
    if schema_digest != PINNED_UPSTREAM_FILES["schema.json"]:
        errors.append("schema.json is not the pinned normative byte sequence")
    if expected_schema_path and schema_path.read_bytes() != Path(expected_schema_path).read_bytes():
        errors.append("schema.json is not byte-identical to the expected normative file")
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise BundleValidationError(f"schema.json is invalid JSON: {exc}") from exc
    manifest = _validate_json_document(root, "manifest.json", schema, "Manifest", errors)
    if manifest is None:
        raise BundleValidationError("manifest.json could not be validated")
    if manifest.get("schema_version") != schema.get("version"):
        errors.append("manifest/schema version mismatch")
    _validate_manifest_boundary(manifest, errors)
    try:
        datetime.fromisoformat(str(manifest.get("generated_at", "")).replace("Z", "+00:00"))
    except ValueError:
        errors.append("manifest.generated_at is not a parseable ISO-8601 datetime")

    metadata_definitions = {
        "field_glossary.json": "FieldGlossary",
        "inventory.json": "Inventory",
        "raw_artifacts.json": "RawArtifacts",
        "transformations.json": "Transformations",
        "omissions.json": "Omissions",
        "relationships.json": "Relationships",
    }
    documents = {
        relative: _validate_json_document(root, relative, schema, definition, errors)
        for relative, definition in metadata_definitions.items()
    }
    predictor_metadata = _validate_predictor_metadata(documents, errors)
    predictor_doctor_metadata = _validate_predictor_doctor_metadata(
        documents,
        manifest,
        list(predictor_metadata.get("enabled_plugins", [])),
        errors,
    )
    _manifest_coverage(root, manifest, errors)
    manifest_rows = _validate_manifest_files(root, manifest, errors)

    variant_rows = _validate_parquet_physical(
        variants_path, variant_schema(), "variants.parquet", errors
    )
    callable_path = root / "callable_regions.parquet"
    callable_rows = 0
    if callable_path.exists():
        callable_rows = _validate_parquet_physical(
            callable_path, callable_region_schema(), "callable_regions.parquet", errors
        )
    target_callability_path = root / "callability.parquet"
    target_callability_rows = 0
    if target_callability_path.exists():
        target_callability_rows = _validate_parquet_physical(
            target_callability_path,
            callability_schema(),
            "callability.parquet",
            errors,
        )
    gene_index_path = root / "gene_index.parquet"
    gene_index_rows = 0
    if gene_index_path.exists():
        gene_index_rows = _validate_parquet_physical(
            gene_index_path, gene_index_schema(), "gene_index.parquet", errors
        )
    clinical_evidence_path = root / "clinical_evidence.parquet"
    clinical_evidence_rows = 0
    if clinical_evidence_path.exists():
        clinical_evidence_rows = _validate_parquet_physical(
            clinical_evidence_path,
            clinical_evidence_schema(),
            "clinical_evidence.parquet",
            errors,
        )
    clinical_findings_path = root / "clinical_findings.parquet"
    clinical_findings_rows = 0
    if clinical_findings_path.exists():
        clinical_findings_rows = _validate_parquet_physical(
            clinical_findings_path,
            clinical_finding_schema(),
            "clinical_findings.parquet",
            errors,
        )
    interpretations_path = root / "interpretations.parquet"
    interpretation_rows = 0
    if interpretations_path.exists():
        interpretation_rows = _validate_parquet_physical(
            interpretations_path,
            interpretation_schema(),
            "interpretations.parquet",
            errors,
        )
    prs_path = root / "prs.parquet"
    prs_rows = 0
    if prs_path.exists():
        prs_rows = _validate_parquet_physical(prs_path, prs_schema(), "prs.parquet", errors)
    pharmacogenomics_path = root / "pharmacogenomics.parquet"
    pharmacogenomics_rows = 0
    if pharmacogenomics_path.exists():
        pharmacogenomics_rows = _validate_parquet_physical(
            pharmacogenomics_path,
            pharmacogenomics_schema(),
            "pharmacogenomics.parquet",
            errors,
        )
    gwas_associations_path = root / "gwas_associations.parquet"
    gwas_association_rows = 0
    if gwas_associations_path.exists():
        gwas_association_rows = _validate_parquet_physical(
            gwas_associations_path,
            gwas_association_schema(),
            "gwas_associations.parquet",
            errors,
        )
    extension_specs = {
        "variant_provenance": (
            root / "extensions" / "variant_provenance.parquet",
            variant_provenance_schema(),
        ),
        "reference_confidence_regions": (
            root / "extensions" / "reference_confidence_regions.parquet",
            reference_evidence_schema(),
        ),
        "observed_partial_variants": (
            root / "extensions" / "observed_partial_variants.parquet",
            observed_partial_variant_schema(),
        ),
        "non_core_variants": (
            root / "extensions" / "non_core_variants.parquet",
            non_core_variant_schema(),
        ),
        "acmg_secondary_findings_evidence": (
            root / "extensions" / "acmg_secondary_findings_evidence.parquet",
            acmg_secondary_findings_evidence_schema(),
        ),
        "clinvar_assertions": (
            root / "extensions" / "clinvar_assertions.parquet",
            clinvar_assertion_schema(),
        ),
        "target_callability_evidence": (
            root / "extensions" / "target_callability_evidence.parquet",
            target_callability_evidence_schema(),
        ),
        "pgs_model_audit": (
            root / "extensions" / "pgs_model_audit.parquet",
            pgs_model_audit_schema(),
        ),
        "pgs_locus_audit": (
            root / "extensions" / "pgs_locus_audit.parquet",
            pgs_locus_audit_schema(),
        ),
        "pgx_locus_audit": (
            root / "extensions" / "pgx_locus_audit.parquet",
            pgx_locus_audit_schema(),
        ),
        "pgx_gene_call_audit": (
            root / "extensions" / "pgx_gene_call_audit.parquet",
            pgx_gene_call_audit_schema(),
        ),
        "pgx_guidance_audit": (
            root / "extensions" / "pgx_guidance_audit.parquet",
            pgx_guidance_audit_schema(),
        ),
        "family_evidence": (
            root / "extensions" / "family_evidence.parquet",
            family_evidence_schema(),
        ),
        "gwas_catalog_evidence": (
            root / "extensions" / "gwas_catalog_evidence.parquet",
            gwas_catalog_evidence_schema(),
        ),
        "clingen_gene_disease": (
            root / "extensions" / "clingen_gene_disease.parquet",
            clingen_gene_disease_schema(),
        ),
        "clingen_dosage_sensitivity": (
            root / "extensions" / "clingen_dosage_sensitivity.parquet",
            clingen_dosage_sensitivity_schema(),
        ),
        **{
            f"release_c/{domain}": (
                root / "extensions" / "release_c" / f"{domain}.parquet",
                release_c_schema_for(domain),
            )
            for domain in RELEASE_C_DOMAINS
        },
    }
    extension_rows: dict[str, int] = {}
    for label, (path, expected) in extension_specs.items():
        if path.exists():
            extension_rows[label] = _validate_parquet_physical(path, expected, label, errors)

    release_c_root = root / "extensions" / "release_c"
    release_c_validation: dict[str, Any] = {"status": "not_present"}
    if release_c_root.exists():
        expected_reference_sha256 = manifest.get("resource_provenance", {}).get(
            "reference_fasta_sha256"
        )
        try:
            release_c_validation = validate_release_c_extension_stage(
                release_c_root,
                expected_reference_sha256=expected_reference_sha256,
            )
        except (PortableEvidenceError, OSError, ValueError) as exc:
            errors.append(f"Release-C extension validation failed: {exc}")

    connection = duckdb.connect(":memory:")
    try:
        views = _register_views(connection, root)
        if "variants" in views:
            metrics["variant_invariants"] = _validate_variant_rows(connection, manifest, errors)
            metrics["predictor_values"] = _validate_predictor_values(
                connection,
                set(predictor_metadata["enabled_plugins"]),
                errors,
            )
        metrics["predictor_metadata"] = predictor_metadata
        metrics["predictor_semantic_controls"] = predictor_doctor_metadata
        metrics.update(_validate_reference_and_callability(connection, views, errors))
        metrics.update(_validate_target_callability(connection, views, manifest, errors))
        metrics["gene_index"] = _validate_gene_index(connection, views, errors)
        metrics["clinical_evidence"] = _validate_clinical_evidence(connection, views, errors)
        metrics["clinical_artifacts"] = _validate_clinical_artifacts(connection, views, errors)
        metrics["clinvar_assertions"] = _validate_clinvar_assertions(connection, views, errors)
        metrics["acmg_secondary_findings"] = _validate_acmg_evidence(
            connection, views, manifest, errors
        )
        metrics["polygenic_scores"] = _validate_pgs(connection, views, manifest, errors)
        metrics["pharmacogenomics"] = _validate_pgx(connection, views, manifest, errors)
        metrics["family_evidence"] = _validate_family_evidence(connection, views, manifest, errors)
        metrics["gwas_catalog"] = _validate_gwas_catalog(connection, views, manifest, errors)
        metrics["clingen_gene_context"] = _validate_clingen_gene_context(
            connection, views, manifest, documents, root, errors
        )
        metrics["provenance"] = _validate_provenance(connection, root, errors)
    finally:
        connection.close()

    reference = Path(reference_fasta).expanduser().resolve() if reference_fasta else None
    if reference is not None and not reference.is_file():
        errors.append(f"reference FASTA not found: {reference}")
        reference = None
    checked_ref_rows = _validate_ref_and_partition_order(variants_path, reference, errors)
    if reference is None:
        warnings.append(
            "GRCh38 REF alleles were not independently checked because no FASTA was supplied"
        )

    glossary = documents.get("field_glossary.json")
    if glossary is not None:
        expected_fields = _field_paths("variants", variant_schema())
        if callable_path.exists():
            expected_fields |= _field_paths("callable_regions", callable_region_schema())
        if target_callability_path.exists():
            expected_fields |= _field_paths("callability", callability_schema())
        if gene_index_path.exists():
            expected_fields |= _field_paths("gene_index", gene_index_schema())
        if clinical_evidence_path.exists():
            expected_fields |= _field_paths("clinical_evidence", clinical_evidence_schema())
        if clinical_findings_path.exists():
            expected_fields |= _field_paths("clinical_findings", clinical_finding_schema())
        if interpretations_path.exists():
            expected_fields |= _field_paths("interpretations", interpretation_schema())
        if prs_path.exists():
            expected_fields |= _field_paths("prs", prs_schema())
        if pharmacogenomics_path.exists():
            expected_fields |= _field_paths("pharmacogenomics", pharmacogenomics_schema())
        if gwas_associations_path.exists():
            expected_fields |= _field_paths("gwas_associations", gwas_association_schema())
        if extension_specs["variant_provenance"][0].exists():
            expected_fields |= _field_paths("variant_provenance", variant_provenance_schema())
        if extension_specs["reference_confidence_regions"][0].exists():
            expected_fields |= _field_paths(
                "reference_confidence_regions", reference_evidence_schema()
            )
        if extension_specs["observed_partial_variants"][0].exists():
            expected_fields |= _field_paths(
                "observed_partial_variants", observed_partial_variant_schema()
            )
        if extension_specs["non_core_variants"][0].exists():
            expected_fields |= _field_paths("non_core_variants", non_core_variant_schema())
        if extension_specs["acmg_secondary_findings_evidence"][0].exists():
            expected_fields |= _field_paths(
                "acmg_secondary_findings_evidence",
                acmg_secondary_findings_evidence_schema(),
            )
        if extension_specs["clinvar_assertions"][0].exists():
            expected_fields |= _field_paths("clinvar_assertions", clinvar_assertion_schema())
        if extension_specs["target_callability_evidence"][0].exists():
            expected_fields |= _field_paths(
                "target_callability_evidence", target_callability_evidence_schema()
            )
        if extension_specs["pgs_model_audit"][0].exists():
            expected_fields |= _field_paths("pgs_model_audit", pgs_model_audit_schema())
        if extension_specs["pgs_locus_audit"][0].exists():
            expected_fields |= _field_paths("pgs_locus_audit", pgs_locus_audit_schema())
        if extension_specs["pgx_locus_audit"][0].exists():
            expected_fields |= _field_paths("pgx_locus_audit", pgx_locus_audit_schema())
        if extension_specs["pgx_gene_call_audit"][0].exists():
            expected_fields |= _field_paths("pgx_gene_call_audit", pgx_gene_call_audit_schema())
        if extension_specs["pgx_guidance_audit"][0].exists():
            expected_fields |= _field_paths("pgx_guidance_audit", pgx_guidance_audit_schema())
        if extension_specs["family_evidence"][0].exists():
            expected_fields |= _field_paths("family_evidence", family_evidence_schema())
        if extension_specs["gwas_catalog_evidence"][0].exists():
            expected_fields |= _field_paths("gwas_catalog_evidence", gwas_catalog_evidence_schema())
        if extension_specs["clingen_gene_disease"][0].exists():
            expected_fields |= _field_paths("clingen_gene_disease", clingen_gene_disease_schema())
        if extension_specs["clingen_dosage_sensitivity"][0].exists():
            expected_fields |= _field_paths(
                "clingen_dosage_sensitivity", clingen_dosage_sensitivity_schema()
            )
        for domain in RELEASE_C_DOMAINS:
            if extension_specs[f"release_c/{domain}"][0].exists():
                expected_fields |= _field_paths(domain, release_c_schema_for(domain))
        actual_fields = set(glossary.get("fields", {}))
        if actual_fields != expected_fields:
            missing = sorted(expected_fields - actual_fields)
            extra = sorted(actual_fields - expected_fields)
            errors.append(
                f"field glossary coverage mismatch; missing={missing[:5]}, extra={extra[:5]}"
            )

    evaluation_cases = _validate_profile(root, manifest, errors)
    if run_profile_queries and manifest.get("profiles", {}).get("agent_exploration"):
        try:
            overview = query_bundle(root, "overview", limit=1)
        except Exception as exc:
            errors.append(f"Agent Profile overview smoke query failed: {exc}")
        else:
            if overview.get("count") != 1:
                errors.append("Agent Profile overview smoke query returned an invalid count")

    metrics.update(
        {
            "variant_rows": variant_rows,
            "callable_region_rows": callable_rows,
            "target_callability_rows": target_callability_rows,
            "gene_index_rows": gene_index_rows,
            "clinical_evidence_rows": clinical_evidence_rows,
            "clinical_findings_rows": clinical_findings_rows,
            "interpretation_rows": interpretation_rows,
            "prs_rows": prs_rows,
            "pharmacogenomics_rows": pharmacogenomics_rows,
            "gwas_association_rows": gwas_association_rows,
            "extension_rows": extension_rows,
            "release_c_extensions": release_c_validation,
            "reference_rows_checked": checked_ref_rows,
            "agent_evaluation_cases": evaluation_cases,
            "manifest_declared_rows": manifest_rows,
        }
    )
    if errors:
        preview = "\n  - ".join(errors[:50])
        suffix = f"\n  ... and {len(errors) - 50} more" if len(errors) > 50 else ""
        raise BundleValidationError(f"bundle validation failed:\n  - {preview}{suffix}")
    return {
        "valid": True,
        "schema_version": manifest["schema_version"],
        "schema_sha256": schema_digest,
        "metrics": metrics,
        "warnings": warnings,
    }
