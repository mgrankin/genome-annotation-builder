"""Build and validate deterministic multi-bundle release reports.

The project report is deliberately generic: private subject identifiers and
filesystem paths are supplied only by the local invocation.  Portable reports
contain the caller-provided bundle IDs, content digests, aggregate counts, and
validated provenance, never host paths.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import os
import re
import shutil
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from statistics import median
from typing import Any

from ..constants import PIPELINE_VERSION
from ..util import sha256_file
from .annotation_qc import (
    QC_IDENTITY_FIELD,
    aggregate_annotation_qc,
    collect_bundle_annotation_qc,
    validate_annotation_qc,
)
from .bundle_extension import validate_enriched_bundle
from .portable_model import (
    CONCORDANCE_STATES,
    FORBIDDEN_PORTABLE_TEXT,
    SOURCE_ROLES,
    PortableEvidenceError,
)

PROJECT_RELEASE_SCHEMA_VERSION = 1
NORMATIVE_SCHEMA_SHA256 = "6d07eaa4bacca4bb2c474eacd8e1a658470475fc4a46438432fb4509d0f11615"
BUNDLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
QC_EXPLANATION_CODE = re.compile(r"^[A-Z][A-Z0-9_]{2,63}$")
OUTLIER_RULES = (
    ("small_variants.total", "relative", 0.12),
    ("small_variants.ti_tv", "relative", 0.05),
    ("small_variants.het_fraction", "relative", 0.08),
    ("vep.gene_fraction", "absolute", 0.05),
    ("gnomad.any_af_fraction", "absolute", 0.03),
    ("difficult_regions.complex_fraction", "absolute", 0.10),
    ("callability.reference_observed_fraction", "absolute", 0.10),
)
REPORT_FILES = (
    "annotation_qc.html",
    "annotation_qc.json",
    "annotation_qc.tsv",
    "conversion_omissions.tsv",
    "family_evidence_summary.html",
    "validation_summary.json",
)
MANIFEST_FILES = (
    "container_digests.tsv",
    "input_manifest.json",
    "input_manifest.tsv",
    "resource_manifest.json",
    "resource_manifest.tsv",
    "software_manifest.json",
)
LOG_FILES = ("project_release.jsonl",)
BUNDLE_VALIDATION_DIRECTORY = "bundle_validation"


def _bundle_validation_relative(bundle_id: str) -> str:
    return f"{BUNDLE_VALIDATION_DIRECTORY}/{bundle_id}.json"


def _expected_report_files(bundle_ids: Sequence[str]) -> set[str]:
    return {
        *REPORT_FILES,
        *(_bundle_validation_relative(bundle_id) for bundle_id in bundle_ids),
    }


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    ).encode("utf-8")


def _canonical_identity(value: object) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _with_identity(value: dict[str, Any], field: str) -> dict[str, Any]:
    if field in value:
        raise PortableEvidenceError("project report identity field already exists")
    return {**value, field: _canonical_identity(value)}


def _canonical_timestamp(value: str) -> str:
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise PortableEvidenceError("generated-at is not ISO-8601") from exc
    if parsed.tzinfo is None:
        raise PortableEvidenceError("generated-at must include a timezone")
    return parsed.astimezone(UTC).isoformat().replace("+00:00", "Z")


def _load_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise PortableEvidenceError(f"{label} is missing or symlinked")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError(f"{label} is unreadable") from exc
    if not isinstance(value, dict):
        raise PortableEvidenceError(f"{label} must be a JSON object")
    return value


def _load_qc_explanations(path: Path | None) -> dict[tuple[str, str], dict[str, str]]:
    if path is None:
        return {}
    raw = Path(path).expanduser()
    if raw.is_symlink() or not raw.is_file():
        raise PortableEvidenceError("QC explanation file is missing or symlinked")
    document = _load_json(raw.resolve(), "QC explanations")
    rows = document.get("explanations")
    if document.get("schema_version") != 1 or not isinstance(rows, list):
        raise PortableEvidenceError("QC explanation document is invalid")
    output: dict[tuple[str, str], dict[str, str]] = {}
    allowed_metrics = {metric for metric, _mode, _threshold in OUTLIER_RULES}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {
            "bundle_id",
            "metric",
            "code",
            "explanation",
        }:
            raise PortableEvidenceError("QC explanation row is invalid")
        bundle_id = row.get("bundle_id")
        metric = row.get("metric")
        code = row.get("code")
        explanation = row.get("explanation")
        if (
            not isinstance(bundle_id, str)
            or BUNDLE_ID.fullmatch(bundle_id) is None
            or metric not in allowed_metrics
            or not isinstance(code, str)
            or QC_EXPLANATION_CODE.fullmatch(code) is None
            or not isinstance(explanation, str)
            or not 1 <= len(explanation) <= 500
            or "\n" in explanation
            or FORBIDDEN_PORTABLE_TEXT.search(explanation)
        ):
            raise PortableEvidenceError("QC explanation value is invalid or non-portable")
        key = bundle_id, metric
        if key in output:
            raise PortableEvidenceError("QC explanation row is duplicated")
        output[key] = {"code": code, "explanation": explanation}
    return output


def _write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def _write_json(path: Path, value: object) -> None:
    _write_bytes(path, _canonical_json_bytes(value))


def _tsv_bytes(rows: Sequence[Sequence[object]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    for row in rows:
        writer.writerow(["" if value is None else value for value in row])
    return stream.getvalue().encode("utf-8")


def _safe_bundle_specs(
    bundles: Mapping[str, str | Path] | Sequence[tuple[str, str | Path]],
) -> list[tuple[str, Path]]:
    items = list(bundles.items()) if isinstance(bundles, Mapping) else list(bundles)
    if not items:
        raise PortableEvidenceError("project release requires at least one bundle")
    seen_ids: set[str] = set()
    seen_paths: set[Path] = set()
    normalized: list[tuple[str, Path]] = []
    for bundle_id, raw_path in items:
        if not isinstance(bundle_id, str) or BUNDLE_ID.fullmatch(bundle_id) is None:
            raise PortableEvidenceError("bundle ID is not portable")
        source = Path(raw_path).expanduser()
        if source.is_symlink() or not source.is_dir():
            raise PortableEvidenceError("bundle root is missing or symlinked")
        resolved = source.resolve()
        if bundle_id in seen_ids or resolved in seen_paths:
            raise PortableEvidenceError("bundle IDs and roots must be unique")
        seen_ids.add(bundle_id)
        seen_paths.add(resolved)
        normalized.append((bundle_id, resolved))
    return sorted(normalized)


def _portable_payload(value: object, label: str) -> None:
    rendered = json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)
    if FORBIDDEN_PORTABLE_TEXT.search(rendered):
        raise PortableEvidenceError(f"{label} contains a host path or secret")


def _selected_metrics(validation: dict[str, Any]) -> dict[str, Any]:
    core = validation.get("core")
    metrics = core.get("metrics") if isinstance(core, dict) else None
    release_c = validation.get("release_c")
    if not isinstance(metrics, dict) or not isinstance(release_c, dict):
        raise PortableEvidenceError("bundle validation report is incomplete")
    names = (
        "variant_rows",
        "callable_region_rows",
        "target_callability_rows",
        "gene_index_rows",
        "clinical_evidence_rows",
        "clinical_findings_rows",
        "interpretation_rows",
        "prs_rows",
        "pharmacogenomics_rows",
        "gwas_association_rows",
        "reference_rows_checked",
    )
    selected: dict[str, Any] = {}
    for name in names:
        value = metrics.get(name)
        if not isinstance(value, int) or value < 0:
            raise PortableEvidenceError(f"bundle metric {name} is invalid")
        selected[name] = value
    for name in ("domains", "rows", "input_gates", "caller_runs", "source_rows_removed"):
        value = release_c.get(name)
        if not isinstance(value, int) or value < 0:
            raise PortableEvidenceError(f"Release-C metric {name} is invalid")
        selected[f"release_c_{name}"] = value
    concordance = release_c.get("concordance")
    source_roles = release_c.get("technical_source_roles")
    if not isinstance(concordance, dict) or not isinstance(source_roles, dict):
        raise PortableEvidenceError("Release-C concordance summary is absent")
    concordance_totals = {state: 0 for state in sorted(CONCORDANCE_STATES)}
    role_totals = {role: 0 for role in sorted(SOURCE_ROLES)}
    for domain, counts in concordance.items():
        if not isinstance(domain, str) or not isinstance(counts, dict):
            raise PortableEvidenceError("Release-C concordance summary is invalid")
        for state, count in counts.items():
            if state not in CONCORDANCE_STATES or not isinstance(count, int) or count < 0:
                raise PortableEvidenceError("Release-C concordance count is invalid")
            concordance_totals[state] += count
    for domain, counts in source_roles.items():
        if not isinstance(domain, str) or not isinstance(counts, dict):
            raise PortableEvidenceError("Release-C technical-source summary is invalid")
        for role, count in counts.items():
            if role not in SOURCE_ROLES or not isinstance(count, int) or count < 0:
                raise PortableEvidenceError("Release-C technical-source count is invalid")
            role_totals[role] += count
    if (
        sum(concordance_totals.values()) != selected["release_c_rows"]
        or sum(role_totals.values()) != selected["release_c_rows"]
    ):
        raise PortableEvidenceError("Release-C concordance cardinality differs")
    selected.update(
        {
            f"release_c_concordance_{state}_rows": count
            for state, count in concordance_totals.items()
        }
    )
    selected.update(
        {f"release_c_source_role_{role}_rows": count for role, count in role_totals.items()}
    )
    if release_c.get("status") != "PASS" or release_c.get("source_rows_removed") != 0:
        raise PortableEvidenceError("Release-C validation did not preserve all source rows")
    return selected


def _bundle_record(
    bundle_id: str,
    root: Path,
    *,
    reference_fasta: Path,
) -> dict[str, Any]:
    validation = validate_enriched_bundle(
        root,
        reference_fasta=reference_fasta,
        run_read_backed_query=True,
    )
    if validation.get("status") != "PASS":
        raise PortableEvidenceError("enriched bundle validation did not return PASS")
    manifest = _load_json(root / "manifest.json", "bundle manifest")
    schema_sha256 = sha256_file(root / "schema.json")
    if schema_sha256 != NORMATIVE_SCHEMA_SHA256:
        raise PortableEvidenceError("bundle does not contain the pinned .genome/1.1 schema")
    raw_artifacts = _load_json(root / "raw_artifacts.json", "raw artifacts")
    transformations = _load_json(root / "transformations.json", "transformations")
    omissions = _load_json(root / "omissions.json", "omissions")
    extension_manifest = _load_json(
        root / "extensions" / "release_c" / "release_c_extensions.json",
        "Release-C extension manifest",
    )
    artifacts = raw_artifacts.get("artifacts")
    stages = transformations.get("stages")
    omission_rows = omissions.get("omissions")
    domains = extension_manifest.get("domains")
    declared_files = manifest.get("files")
    if not isinstance(artifacts, list) or not all(isinstance(row, dict) for row in artifacts):
        raise PortableEvidenceError("raw artifact rows are invalid")
    if not isinstance(stages, list) or not all(isinstance(row, dict) for row in stages):
        raise PortableEvidenceError("transformation rows are invalid")
    if not isinstance(omission_rows, list) or not all(
        isinstance(row, dict) for row in omission_rows
    ):
        raise PortableEvidenceError("omission rows are invalid")
    if not isinstance(domains, dict):
        raise PortableEvidenceError("Release-C domain manifest is invalid")
    if not isinstance(declared_files, dict) or not all(
        isinstance(path, str) and isinstance(row, dict) for path, row in declared_files.items()
    ):
        raise PortableEvidenceError("bundle file inventory is invalid")
    domain_rows: dict[str, int] = {}
    for name, metadata in sorted(domains.items()):
        if not isinstance(name, str) or not isinstance(metadata, dict):
            raise PortableEvidenceError("Release-C domain record is invalid")
        rows = metadata.get("rows")
        if not isinstance(rows, int) or rows < 1:
            raise PortableEvidenceError("Release-C domain row count is invalid")
        domain_rows[name] = rows
    privacy = manifest.get("privacy")
    if not isinstance(privacy, dict) or any(privacy.get(key) is not False for key in privacy):
        raise PortableEvidenceError("bundle privacy declaration is incomplete")
    source_bundle_sha256 = manifest.get("source_bundle_sha256")
    if not isinstance(source_bundle_sha256, str) or SHA256.fullmatch(source_bundle_sha256) is None:
        raise PortableEvidenceError("source bundle identity is invalid")
    annotation_qc = collect_bundle_annotation_qc(root, validation=validation)
    record = {
        "bundle_id": bundle_id,
        "manifest_sha256": sha256_file(root / "manifest.json"),
        "schema_sha256": schema_sha256,
        "source_bundle_sha256": source_bundle_sha256,
        "pipeline_version": manifest.get("pipeline_version"),
        "generated_at": manifest.get("generated_at"),
        "genome_build": manifest.get("genome_build"),
        "metrics": _selected_metrics(validation),
        "annotation_qc": annotation_qc,
        "domain_rows": domain_rows,
        "warnings": validation.get("core", {}).get("warnings", []),
        "query": validation.get("query"),
        "family_evidence": manifest.get("resource_provenance", {}).get("family_evidence", {}),
        "resource_provenance": manifest.get("resource_provenance", {}),
        "declared_files": declared_files,
        "raw_artifacts": artifacts,
        "transformations": stages,
        "omissions": omission_rows,
        "validation": validation,
    }
    _portable_payload(record, "bundle report")
    return record


def _bundle_validation_report(record: dict[str, Any], generated_at: str) -> dict[str, Any]:
    document = {
        "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
        "status": "PASS",
        "generated_at": generated_at,
        "pipeline_version": PIPELINE_VERSION,
        "bundle_id": record["bundle_id"],
        "bundle_manifest_sha256": record["manifest_sha256"],
        "source_bundle_sha256": record["source_bundle_sha256"],
        "normative_schema_sha256": record["schema_sha256"],
        "file_hashes": record["declared_files"],
        "metrics": record["metrics"],
        "annotation_qc": record["annotation_qc"],
        "annotation_sources": record["resource_provenance"],
        "validation": record["validation"],
        "query_regression": record["query"],
        "errors": [],
        "warnings": record["warnings"],
        "privacy": {
            "bundle_path_persisted": False,
            "private_manifest_persisted": False,
            "relative_genotypes_persisted": False,
        },
    }
    return _with_identity(document, "bundle_validation_identity_sha256")


def _outlier_values(document: Mapping[str, Any]) -> dict[str, float | None]:
    total = float(_qc_value(document, "counts", "small_variants.total", "value"))
    het = float(_qc_distribution_value(document, "small_variants.zygosity", "het"))
    return {
        "small_variants.total": total,
        "small_variants.ti_tv": _qc_value(document, "ratios", "small_variants.ti_tv", "value"),
        "small_variants.het_fraction": round(het / total, 12) if total else None,
        "vep.gene_fraction": _qc_value(document, "completeness", "vep.gene", "fraction"),
        "gnomad.any_af_fraction": _qc_value(document, "completeness", "gnomad.any_af", "fraction"),
        "difficult_regions.complex_fraction": _qc_value(
            document, "ratios", "difficult_regions.complex_fraction", "value"
        ),
        "callability.reference_observed_fraction": _qc_value(
            document,
            "ratios",
            "callability.reference_observed_fraction",
            "value",
        ),
    }


def _outlier_analysis(
    records: Sequence[Mapping[str, Any]],
    explanations: Mapping[tuple[str, str], Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    explanation_map = dict(explanations or {})
    rules = [
        {"metric": metric, "deviation": mode, "threshold": threshold}
        for metric, mode, threshold in OUTLIER_RULES
    ]
    if len(records) < 3:
        if explanation_map:
            raise PortableEvidenceError("QC explanations exist for an unevaluable cohort")
        return {
            "status": "PASS",
            "method": "fixed deviation from the cross-bundle median",
            "cohort_evaluable": False,
            "minimum_bundle_count": 3,
            "flag_count": 0,
            "rules": rules,
            "flags": [],
        }
    values_by_bundle: dict[str, dict[str, float | None]] = {}
    for record in records:
        bundle_id = record.get("bundle_id")
        document = record.get("annotation_qc")
        if not isinstance(bundle_id, str) or not isinstance(document, Mapping):
            raise PortableEvidenceError("outlier analysis bundle input is invalid")
        values_by_bundle[bundle_id] = _outlier_values(document)

    flags: list[dict[str, Any]] = []
    for metric, mode, threshold in OUTLIER_RULES:
        available = [
            float(values[metric])
            for values in values_by_bundle.values()
            if values[metric] is not None
        ]
        if len(available) < 3:
            continue
        center = float(median(available))
        for bundle_id, values in sorted(values_by_bundle.items()):
            value = values[metric]
            if value is None:
                continue
            value = float(value)
            absolute = abs(value - center)
            deviation = absolute if mode == "absolute" or center == 0 else absolute / abs(center)
            if deviation <= threshold:
                continue
            key = bundle_id, metric
            explanation = explanation_map.get(key)
            if not isinstance(explanation, Mapping):
                raise PortableEvidenceError(
                    f"annotation QC outlier requires explanation: {bundle_id}/{metric}"
                )
            code = explanation.get("code")
            text = explanation.get("explanation")
            if (
                not isinstance(code, str)
                or QC_EXPLANATION_CODE.fullmatch(code) is None
                or not isinstance(text, str)
                or not 1 <= len(text) <= 500
                or "\n" in text
                or FORBIDDEN_PORTABLE_TEXT.search(text)
            ):
                raise PortableEvidenceError("annotation QC outlier explanation is invalid")
            flags.append(
                {
                    "bundle_id": bundle_id,
                    "metric": metric,
                    "value": round(value, 12),
                    "median": round(center, 12),
                    "deviation": round(deviation, 12),
                    "deviation_kind": mode if center != 0 else "absolute_zero_center",
                    "threshold": threshold,
                    "explanation_code": code,
                    "explanation": text,
                }
            )
    flagged_keys = {(row["bundle_id"], row["metric"]) for row in flags}
    if set(explanation_map) != flagged_keys:
        raise PortableEvidenceError("QC explanations contain stale or non-outlier entries")
    return {
        "status": "PASS",
        "method": "fixed deviation from the cross-bundle median",
        "cohort_evaluable": True,
        "minimum_bundle_count": 3,
        "flag_count": len(flags),
        "rules": rules,
        "flags": flags,
    }


def _annotation_qc(
    records: list[dict[str, Any]],
    generated_at: str,
    *,
    qc_explanations: Mapping[tuple[str, str], Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    metric_names = sorted(records[0]["metrics"])
    aggregate = {
        name: sum(int(record["metrics"][name]) for record in records) for name in metric_names
    }
    document = {
        "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
        "status": "PASS",
        "generated_at": generated_at,
        "pipeline_version": PIPELINE_VERSION,
        "bundle_count": len(records),
        "aggregate_metrics": aggregate,
        "aggregate_qc": aggregate_annotation_qc([record["annotation_qc"] for record in records]),
        "outlier_analysis": _outlier_analysis(records, qc_explanations),
        "bundles": [
            {
                "bundle_id": record["bundle_id"],
                "manifest_sha256": record["manifest_sha256"],
                "metrics": record["metrics"],
                "annotation_qc": record["annotation_qc"],
                "domain_rows": record["domain_rows"],
                "warning_count": len(record["warnings"]),
            }
            for record in records
        ],
    }
    return _with_identity(document, "annotation_qc_identity_sha256")


def _qc_tsv(records: list[dict[str, Any]], outlier_analysis: Mapping[str, Any]) -> bytes:
    aggregate = aggregate_annotation_qc([record["annotation_qc"] for record in records])
    rows: list[list[object]] = [
        [
            "scope",
            "bundle_id",
            "kind",
            "metric",
            "category",
            "value",
            "numerator",
            "denominator",
            "fraction_or_ratio",
            "unit",
            "detail",
        ]
    ]
    documents = [(record["bundle_id"], record["annotation_qc"]) for record in records]
    documents.append(("PROJECT", aggregate))
    for bundle_id, document in documents:
        scope = document["scope"]
        for row in document["counts"]:
            rows.append(
                [
                    scope,
                    bundle_id,
                    "count",
                    row["metric"],
                    "",
                    row["value"],
                    "",
                    "",
                    "",
                    row["unit"],
                    "",
                ]
            )
        for row in document["distributions"]:
            rows.append(
                [
                    scope,
                    bundle_id,
                    "distribution",
                    row["metric"],
                    row["category"],
                    row["count"],
                    "",
                    "",
                    "",
                    "rows",
                    "",
                ]
            )
        for row in document["ratios"]:
            rows.append(
                [
                    scope,
                    bundle_id,
                    "ratio",
                    row["metric"],
                    f"{row['numerator_label']}/{row['denominator_label']}",
                    "",
                    row["numerator"],
                    row["denominator"],
                    row["value"],
                    "ratio",
                    "",
                ]
            )
        for row in document["completeness"]:
            rows.append(
                [
                    scope,
                    bundle_id,
                    "completeness",
                    row["metric"],
                    "populated/applicable",
                    row["populated_rows"],
                    row["populated_rows"],
                    row["applicable_rows"],
                    row["fraction"],
                    "fraction",
                    "",
                ]
            )
    flags = outlier_analysis.get("flags")
    if not isinstance(flags, list):
        raise PortableEvidenceError("annotation QC outlier flags are invalid")
    for row in flags:
        if not isinstance(row, dict):
            raise PortableEvidenceError("annotation QC outlier row is invalid")
        rows.append(
            [
                "project",
                row.get("bundle_id"),
                "outlier",
                row.get("metric"),
                row.get("explanation_code"),
                row.get("value"),
                row.get("deviation"),
                row.get("threshold"),
                "",
                row.get("deviation_kind"),
                row.get("explanation"),
            ]
        )
    return _tsv_bytes(rows)


def _html_table(title: str, headers: Sequence[str], rows: Sequence[Sequence[object]]) -> bytes:
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
        for row in rows
    )
    document = (
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        f"<title>{html.escape(title)}</title>"
        "<style>body{font-family:system-ui,sans-serif;margin:2rem}"
        "table{border-collapse:collapse}th,td{border:1px solid #bbb;padding:.35rem}"
        "th{background:#eee}</style>"
        f"<h1>{html.escape(title)}</h1><table><thead><tr>{head}</tr></thead>"
        f"<tbody>{body}</tbody></table></html>\n"
    )
    return document.encode("utf-8")


def _omission_tsv(records: list[dict[str, Any]]) -> bytes:
    rows: list[list[object]] = [
        [
            "bundle_id",
            "ordinal",
            "scope",
            "pipeline_stage",
            "reason",
            "interpretation_consequence",
            "recovery",
            "row_json",
        ]
    ]
    for record in records:
        for ordinal, row in enumerate(record["omissions"], start=1):
            rows.append(
                [
                    record["bundle_id"],
                    ordinal,
                    row.get("scope"),
                    row.get("pipeline_stage"),
                    row.get("reason"),
                    row.get("interpretation_consequence"),
                    row.get("recovery"),
                    json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False),
                ]
            )
    return _tsv_bytes(rows)


def _input_manifest(records: list[dict[str, Any]], generated_at: str) -> dict[str, Any]:
    document = {
        "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
        "status": "PASS",
        "generated_at": generated_at,
        "bundles": [
            {
                "bundle_id": record["bundle_id"],
                "bundle_manifest_sha256": record["manifest_sha256"],
                "source_bundle_sha256": record["source_bundle_sha256"],
                "artifacts": record["raw_artifacts"],
            }
            for record in records
        ],
    }
    return _with_identity(document, "input_manifest_identity_sha256")


def _input_tsv(records: list[dict[str, Any]]) -> bytes:
    rows: list[list[object]] = [
        ["bundle_id", "ordinal", "artifact_type", "disposition", "integrity_identifier"]
    ]
    for record in records:
        for ordinal, artifact in enumerate(record["raw_artifacts"], start=1):
            rows.append(
                [
                    record["bundle_id"],
                    ordinal,
                    artifact.get("artifact_type"),
                    artifact.get("disposition"),
                    artifact.get("integrity_identifier"),
                ]
            )
    return _tsv_bytes(rows)


def _resource_manifest(records: list[dict[str, Any]], generated_at: str) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for record in records:
        provenance = record["resource_provenance"]
        if not isinstance(provenance, dict):
            raise PortableEvidenceError("resource provenance is invalid")
        for name, payload in sorted(provenance.items()):
            rows.append(
                {
                    "bundle_id": record["bundle_id"],
                    "resource_name": name,
                    "payload_identity_sha256": _canonical_identity(payload),
                    "payload": payload,
                }
            )
    document = {
        "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
        "status": "PASS",
        "generated_at": generated_at,
        "resources": rows,
    }
    return _with_identity(document, "resource_manifest_identity_sha256")


def _resource_tsv(resource_manifest: dict[str, Any]) -> bytes:
    rows: list[list[object]] = [["bundle_id", "resource_name", "payload_identity_sha256"]]
    rows.extend(
        [row["bundle_id"], row["resource_name"], row["payload_identity_sha256"]]
        for row in resource_manifest["resources"]
    )
    return _tsv_bytes(rows)


def _software_manifest(records: list[dict[str, Any]], generated_at: str) -> dict[str, Any]:
    rows: dict[str, dict[str, Any]] = {}
    for record in records:
        for stage in record["transformations"]:
            portable = {
                key: stage.get(key)
                for key in (
                    "stage_id",
                    "software",
                    "software_version",
                    "container_image",
                    "configuration_digest",
                    "status",
                )
                if key in stage
            }
            identity = _canonical_identity(portable)
            rows[identity] = {"identity_sha256": identity, **portable}
    document = {
        "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
        "status": "PASS",
        "generated_at": generated_at,
        "software_boundaries": [rows[key] for key in sorted(rows)],
    }
    return _with_identity(document, "software_manifest_identity_sha256")


def _container_rows(value: object, path: str = "") -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    if isinstance(value, dict):
        for key, child in sorted(value.items()):
            child_path = f"{path}.{key}" if path else str(key)
            rows.extend(_container_rows(child, child_path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            rows.extend(_container_rows(child, f"{path}[{index}]"))
    elif (
        isinstance(value, str)
        and "sha256:" in value
        and (
            "@sha256:" in value
            or any(token in path.lower() for token in ("container", "image", "oci"))
        )
    ):
        digest = value.rsplit("sha256:", 1)[-1]
        if SHA256.fullmatch(digest):
            rows.append((path, value))
    return rows


def _container_tsv(records: list[dict[str, Any]]) -> bytes:
    unique: set[tuple[str, str, str]] = set()
    for record in records:
        for key_path, value in _container_rows(
            {
                "resources": record["resource_provenance"],
                "transformations": record["transformations"],
            }
        ):
            unique.add((record["bundle_id"], key_path, value))
    return _tsv_bytes(
        [["bundle_id", "key_path", "image_or_digest"], *[list(row) for row in sorted(unique)]]
    )


def _validation_summary(
    records: list[dict[str, Any]],
    generated_at: str,
    bundle_reports: Mapping[str, dict[str, Any]],
) -> dict[str, Any]:
    document = {
        "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
        "status": "PASS",
        "generated_at": generated_at,
        "pipeline_version": PIPELINE_VERSION,
        "normative_schema_sha256": NORMATIVE_SCHEMA_SHA256,
        "bundle_count": len(records),
        "privacy": {
            "bundle_paths_persisted": False,
            "private_manifests_persisted": False,
        },
        "bundles": [
            {
                "bundle_id": record["bundle_id"],
                "manifest_sha256": record["manifest_sha256"],
                "schema_sha256": record["schema_sha256"],
                "source_bundle_sha256": record["source_bundle_sha256"],
                "bundle_validation_identity_sha256": bundle_reports[record["bundle_id"]][
                    "bundle_validation_identity_sha256"
                ],
                "annotation_qc_identity_sha256": record["annotation_qc"][QC_IDENTITY_FIELD],
                "metrics": record["metrics"],
                "domain_rows": record["domain_rows"],
                "query": record["query"],
                "warnings": record["warnings"],
                "validation": record["validation"],
            }
            for record in records
        ],
    }
    return _with_identity(document, "validation_summary_identity_sha256")


def _project_log(validation_summary: dict[str, Any]) -> bytes:
    rows = validation_summary.get("bundles")
    generated_at = validation_summary.get("generated_at")
    summary_identity = validation_summary.get("validation_summary_identity_sha256")
    if (
        not isinstance(rows, list)
        or not isinstance(generated_at, str)
        or not isinstance(summary_identity, str)
    ):
        raise PortableEvidenceError("validation summary cannot produce a project log")
    events: list[dict[str, Any]] = [
        {
            "bundle_count": len(rows),
            "event": "project_validation_started",
            "generated_at": generated_at,
            "pipeline_version": PIPELINE_VERSION,
            "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
            "status": "PASS",
        }
    ]
    for row in rows:
        if not isinstance(row, dict):
            raise PortableEvidenceError("validation summary bundle row is invalid")
        events.append(
            {
                "bundle_id": row.get("bundle_id"),
                "event": "bundle_validated",
                "annotation_qc_identity_sha256": row.get("annotation_qc_identity_sha256"),
                "manifest_sha256": row.get("manifest_sha256"),
                "schema_sha256": row.get("schema_sha256"),
                "source_bundle_sha256": row.get("source_bundle_sha256"),
                "status": "PASS",
            }
        )
    events.append(
        {
            "bundle_count": len(rows),
            "event": "project_validation_completed",
            "generated_at": generated_at,
            "status": "PASS",
            "validation_summary_identity_sha256": summary_identity,
        }
    )
    return b"".join(
        (
            json.dumps(
                event,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
        for event in events
    )


def _family_html(records: list[dict[str, Any]]) -> bytes:
    rows: list[list[object]] = []
    for record in records:
        evidence = record["family_evidence"]
        rows.append(
            [
                record["bundle_id"],
                evidence.get("family_evidence_rows") if isinstance(evidence, dict) else None,
                evidence.get("genotype_refinements") if isinstance(evidence, dict) else None,
                evidence.get("relative_identifiers_persisted")
                if isinstance(evidence, dict)
                else None,
                evidence.get("relative_genotypes_persisted")
                if isinstance(evidence, dict)
                else None,
            ]
        )
    return _html_table(
        "Family evidence summary",
        (
            "bundle_id",
            "family_evidence_rows",
            "genotype_refinements",
            "relative_identifiers_persisted",
            "relative_genotypes_persisted",
        ),
        rows,
    )


def _qc_value(document: Mapping[str, Any], section: str, metric: str, key: str) -> Any:
    rows = document.get(section)
    if not isinstance(rows, list):
        raise PortableEvidenceError(f"annotation QC {section} is invalid")
    matches = [row for row in rows if isinstance(row, dict) and row.get("metric") == metric]
    if len(matches) != 1:
        raise PortableEvidenceError(f"annotation QC metric {metric} is absent or duplicated")
    return matches[0].get(key)


def _qc_distribution_value(document: Mapping[str, Any], metric: str, category: str) -> int:
    rows = document.get("distributions")
    if not isinstance(rows, list):
        raise PortableEvidenceError("annotation QC distributions are invalid")
    matches = [
        row
        for row in rows
        if isinstance(row, dict) and row.get("metric") == metric and row.get("category") == category
    ]
    if not matches:
        return 0
    if len(matches) != 1 or not isinstance(matches[0].get("count"), int):
        raise PortableEvidenceError(f"annotation QC distribution {metric}/{category} is invalid")
    return int(matches[0]["count"])


def _qc_distribution_suffix_total(
    document: Mapping[str, Any], metric_prefix: str, category: str
) -> int:
    rows = document.get("distributions")
    if not isinstance(rows, list):
        raise PortableEvidenceError("annotation QC distributions are invalid")
    return sum(
        int(row["count"])
        for row in rows
        if isinstance(row, dict)
        and isinstance(row.get("metric"), str)
        and row["metric"].startswith(metric_prefix)
        and row.get("category") == category
        and isinstance(row.get("count"), int)
    )


def _qc_percent(value: Any) -> str:
    if value is None:
        return "not available"
    if not isinstance(value, (int, float)):
        raise PortableEvidenceError("annotation QC fraction is invalid")
    return f"{100 * float(value):.4f}%"


def _qc_html(records: list[dict[str, Any]], outlier_analysis: Mapping[str, Any]) -> bytes:
    aggregate = aggregate_annotation_qc([record["annotation_qc"] for record in records])
    documents = [(record["bundle_id"], record["annotation_qc"]) for record in records]
    documents.append(("PROJECT", aggregate))
    flags = outlier_analysis.get("flags")
    if not isinstance(flags, list):
        raise PortableEvidenceError("annotation QC outlier flags are invalid")
    outlier_counts: dict[str, int] = defaultdict(int)
    for row in flags:
        if not isinstance(row, dict) or not isinstance(row.get("bundle_id"), str):
            raise PortableEvidenceError("annotation QC outlier row is invalid")
        outlier_counts[row["bundle_id"]] += 1
    outlier_counts["PROJECT"] = len(flags)
    return _html_table(
        "Genome annotation QC",
        (
            "bundle_id",
            "variants",
            "SNV",
            "indel",
            "Ti/Tv",
            "het",
            "hom-alt",
            "haploid",
            "VEP gene",
            "ClinVar matched",
            "gnomAD any AF",
            "difficult region",
            "target reference",
            "target unknown",
            "PGx gene audits",
            "valid PRS models",
            "Release-C rows",
            "Technical-replicate rows",
            "Discordant rows",
            "QC outlier flags",
        ),
        [
            [
                bundle_id,
                _qc_value(document, "counts", "small_variants.total", "value"),
                _qc_distribution_value(document, "small_variants.variant_type", "SNV"),
                _qc_value(document, "counts", "small_variants.indel_total", "value"),
                _qc_value(document, "ratios", "small_variants.ti_tv", "value"),
                _qc_distribution_value(document, "small_variants.zygosity", "het"),
                _qc_distribution_value(document, "small_variants.zygosity", "hom_alt"),
                _qc_value(document, "counts", "small_variants.haploid", "value"),
                _qc_percent(_qc_value(document, "completeness", "vep.gene", "fraction")),
                _qc_distribution_value(document, "clinvar.match", "matched"),
                _qc_percent(_qc_value(document, "completeness", "gnomad.any_af", "fraction")),
                _qc_percent(
                    _qc_value(document, "ratios", "difficult_regions.complex_fraction", "value")
                ),
                _qc_percent(
                    _qc_value(
                        document,
                        "ratios",
                        "callability.reference_observed_fraction",
                        "value",
                    )
                ),
                _qc_percent(_qc_value(document, "ratios", "callability.unknown_fraction", "value")),
                _qc_value(document, "counts", "pgx.gene_rows", "value"),
                _qc_value(document, "counts", "prs.valid_models", "value"),
                _qc_value(document, "counts", "release_c.rows", "value"),
                _qc_distribution_suffix_total(
                    document,
                    "release_c.technical_source_role.",
                    "technical_replicate",
                ),
                _qc_distribution_suffix_total(
                    document,
                    "release_c.concordance.",
                    "comparable_discordant",
                ),
                outlier_counts.get(bundle_id, 0),
            ]
            for bundle_id, document in documents
        ],
    )


def _directory_payloads(root: Path) -> dict[str, bytes]:
    if root.is_symlink() or not root.is_dir():
        raise PortableEvidenceError("project report directory is missing or symlinked")
    payloads: dict[str, bytes] = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise PortableEvidenceError("project report contains a symlink")
        if path.is_file():
            payloads[path.relative_to(root).as_posix()] = path.read_bytes()
        elif not path.is_dir():
            raise PortableEvidenceError("project report contains a special file")
    return payloads


def _publish_directory(source: Path, destination: Path) -> str:
    expected = _directory_payloads(source)
    if destination.exists():
        if destination.is_symlink() or not destination.is_dir():
            raise PortableEvidenceError("project report destination is unsafe")
        if _directory_payloads(destination) != expected:
            raise PortableEvidenceError("published project report directory differs")
        shutil.rmtree(source)
        return "reused"
    os.replace(source, destination)
    return "published"


def _publish_file(source: Path, destination: Path) -> str:
    expected = source.read_bytes()
    if destination.exists():
        if destination.is_symlink() or not destination.is_file():
            raise PortableEvidenceError("project report file destination is unsafe")
        if destination.read_bytes() != expected:
            raise PortableEvidenceError("published project report file differs")
        source.unlink()
        return "reused"
    os.replace(source, destination)
    return "published"


def _checksum_payload(output_root: Path, bundle_ids: Sequence[str]) -> bytes:
    paths: list[Path] = []
    for bundle_id in bundle_ids:
        paths.extend(
            [
                output_root / f"{bundle_id}.genome" / "manifest.json",
                output_root / f"{bundle_id}.genome" / "schema.json",
            ]
        )
    paths.extend((output_root / "reports" / name) for name in REPORT_FILES)
    paths.extend(
        output_root / "reports" / _bundle_validation_relative(bundle_id) for bundle_id in bundle_ids
    )
    paths.extend((output_root / "manifests" / name) for name in MANIFEST_FILES)
    paths.extend((output_root / "logs" / name) for name in LOG_FILES)
    lines: list[str] = []
    for path in sorted(paths):
        if path.is_symlink() or not path.is_file():
            raise PortableEvidenceError("project checksum member is missing or unsafe")
        lines.append(f"{sha256_file(path)}  {path.relative_to(output_root).as_posix()}")
    return ("\n".join(lines) + "\n").encode("utf-8")


def build_project_release(
    bundles: Mapping[str, str | Path] | Sequence[tuple[str, str | Path]],
    output_root: str | Path,
    *,
    reference_fasta: str | Path,
    generated_at: str,
    minimum_free_bytes: int = 0,
    qc_explanations: Mapping[tuple[str, str], Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    specs = _safe_bundle_specs(bundles)
    raw_output = Path(output_root).expanduser()
    if raw_output.is_symlink() or not raw_output.is_dir():
        raise PortableEvidenceError("project output root must already be a regular directory")
    output = raw_output.resolve()
    reference = Path(reference_fasta).expanduser()
    if reference.is_symlink() or not reference.is_file():
        raise PortableEvidenceError("reference FASTA is missing or symlinked")
    reference = reference.resolve()
    if minimum_free_bytes < 0:
        raise PortableEvidenceError("minimum free bytes cannot be negative")
    if shutil.disk_usage(output).free < minimum_free_bytes:
        raise PortableEvidenceError("project report would cross the configured free-space floor")
    for bundle_id, root in specs:
        expected = output / f"{bundle_id}.genome"
        if expected.is_symlink() or root != expected.resolve():
            raise PortableEvidenceError("bundle root does not match its portable output name")
    generated_at = _canonical_timestamp(generated_at)
    records = [
        _bundle_record(bundle_id, root, reference_fasta=reference) for bundle_id, root in specs
    ]
    annotation_qc = _annotation_qc(
        records,
        generated_at,
        qc_explanations=qc_explanations,
    )
    input_manifest = _input_manifest(records, generated_at)
    resource_manifest = _resource_manifest(records, generated_at)
    software_manifest = _software_manifest(records, generated_at)
    bundle_reports = {
        record["bundle_id"]: _bundle_validation_report(record, generated_at) for record in records
    }
    validation_summary = _validation_summary(records, generated_at, bundle_reports)
    for label, document in (
        ("annotation QC", annotation_qc),
        ("input manifest", input_manifest),
        ("resource manifest", resource_manifest),
        ("software manifest", software_manifest),
        ("validation summary", validation_summary),
    ):
        _portable_payload(document, label)
    for bundle_id, document in bundle_reports.items():
        _portable_payload(document, f"bundle validation report {bundle_id}")

    attempt = Path(tempfile.mkdtemp(prefix=".project-release.partial.", dir=output))
    reports = attempt / "reports"
    manifests = attempt / "manifests"
    logs = attempt / "logs"
    reports.mkdir()
    manifests.mkdir()
    logs.mkdir()
    _write_json(reports / "annotation_qc.json", annotation_qc)
    _write_bytes(
        reports / "annotation_qc.tsv",
        _qc_tsv(records, annotation_qc["outlier_analysis"]),
    )
    _write_bytes(
        reports / "annotation_qc.html",
        _qc_html(records, annotation_qc["outlier_analysis"]),
    )
    _write_bytes(reports / "conversion_omissions.tsv", _omission_tsv(records))
    _write_bytes(reports / "family_evidence_summary.html", _family_html(records))
    _write_json(reports / "validation_summary.json", validation_summary)
    for bundle_id, document in bundle_reports.items():
        _write_json(reports / _bundle_validation_relative(bundle_id), document)
    _write_json(manifests / "input_manifest.json", input_manifest)
    _write_bytes(manifests / "input_manifest.tsv", _input_tsv(records))
    _write_json(manifests / "resource_manifest.json", resource_manifest)
    _write_bytes(manifests / "resource_manifest.tsv", _resource_tsv(resource_manifest))
    _write_json(manifests / "software_manifest.json", software_manifest)
    _write_bytes(manifests / "container_digests.tsv", _container_tsv(records))
    _write_bytes(logs / "project_release.jsonl", _project_log(validation_summary))
    for directory in (reports, manifests, logs):
        for payload in _directory_payloads(directory).values():
            if FORBIDDEN_PORTABLE_TEXT.search(payload.decode("utf-8")):
                raise PortableEvidenceError("project metadata contains a host path")
    publication = {
        "reports": _publish_directory(reports, output / "reports"),
        "manifests": _publish_directory(manifests, output / "manifests"),
        "logs": _publish_directory(logs, output / "logs"),
    }
    checksum = attempt / "SHA256SUMS"
    _write_bytes(checksum, _checksum_payload(output, [bundle_id for bundle_id, _ in specs]))
    publication["SHA256SUMS"] = _publish_file(checksum, output / "SHA256SUMS")
    attempt.rmdir()
    validation = validate_project_release(
        specs,
        output,
        reference_fasta=reference,
        revalidate_bundles=False,
    )
    return {**validation, "publication": publication}


def validate_project_release(
    bundles: Mapping[str, str | Path] | Sequence[tuple[str, str | Path]],
    output_root: str | Path,
    *,
    reference_fasta: str | Path,
    revalidate_bundles: bool = True,
) -> dict[str, Any]:
    specs = _safe_bundle_specs(bundles)
    raw_output = Path(output_root).expanduser()
    if raw_output.is_symlink() or not raw_output.is_dir():
        raise PortableEvidenceError("project output root is missing or symlinked")
    output = raw_output.resolve()
    reference = Path(reference_fasta).expanduser()
    if reference.is_symlink() or not reference.is_file():
        raise PortableEvidenceError("reference FASTA is missing or symlinked")
    reference = reference.resolve()
    bundle_ids = [bundle_id for bundle_id, _ in specs]
    expected_reports = _expected_report_files(bundle_ids)
    expected_manifests = set(MANIFEST_FILES)
    expected_logs = set(LOG_FILES)
    if set(_directory_payloads(output / "reports")) != expected_reports:
        raise PortableEvidenceError("project report file set differs")
    if set(_directory_payloads(output / "manifests")) != expected_manifests:
        raise PortableEvidenceError("project manifest file set differs")
    if set(_directory_payloads(output / "logs")) != expected_logs:
        raise PortableEvidenceError("project log file set differs")
    summary = _load_json(output / "reports" / "validation_summary.json", "validation summary")
    supplied_identity = summary.get("validation_summary_identity_sha256")
    unsigned = {
        key: value for key, value in summary.items() if key != "validation_summary_identity_sha256"
    }
    if (
        summary.get("status") != "PASS"
        or summary.get("schema_version") != PROJECT_RELEASE_SCHEMA_VERSION
        or summary.get("normative_schema_sha256") != NORMATIVE_SCHEMA_SHA256
        or supplied_identity != _canonical_identity(unsigned)
    ):
        raise PortableEvidenceError("validation summary identity differs")
    rows = summary.get("bundles")
    if not isinstance(rows, list) or len(rows) != len(specs):
        raise PortableEvidenceError("validation summary bundle cardinality differs")
    by_id = {
        row.get("bundle_id"): row
        for row in rows
        if isinstance(row, dict) and isinstance(row.get("bundle_id"), str)
    }
    if set(by_id) != {bundle_id for bundle_id, _ in specs}:
        raise PortableEvidenceError("validation summary bundle IDs differ")
    if (output / "logs" / "project_release.jsonl").read_bytes() != _project_log(summary):
        raise PortableEvidenceError("project log differs from validation summary")
    identity_documents = (
        (
            output / "reports" / "annotation_qc.json",
            "annotation_qc_identity_sha256",
            "annotation QC",
        ),
        (
            output / "manifests" / "input_manifest.json",
            "input_manifest_identity_sha256",
            "input manifest",
        ),
        (
            output / "manifests" / "resource_manifest.json",
            "resource_manifest_identity_sha256",
            "resource manifest",
        ),
        (
            output / "manifests" / "software_manifest.json",
            "software_manifest_identity_sha256",
            "software manifest",
        ),
    )
    for path, identity_field, label in identity_documents:
        document = _load_json(path, label)
        identity = document.get(identity_field)
        unsigned_document = {key: value for key, value in document.items() if key != identity_field}
        if (
            document.get("schema_version") != PROJECT_RELEASE_SCHEMA_VERSION
            or document.get("status") != "PASS"
            or identity != _canonical_identity(unsigned_document)
        ):
            raise PortableEvidenceError(f"{label} identity differs")

    annotation_document = _load_json(output / "reports" / "annotation_qc.json", "annotation QC")
    annotation_rows = annotation_document.get("bundles")
    if not isinstance(annotation_rows, list) or len(annotation_rows) != len(specs):
        raise PortableEvidenceError("annotation QC bundle cardinality differs")
    annotation_by_id = {
        row.get("bundle_id"): row
        for row in annotation_rows
        if isinstance(row, dict) and isinstance(row.get("bundle_id"), str)
    }
    if set(annotation_by_id) != set(by_id):
        raise PortableEvidenceError("annotation QC bundle IDs differ")
    bundle_qc_documents: list[dict[str, Any]] = []
    for bundle_id, _root in specs:
        annotation_row = annotation_by_id[bundle_id]
        summary_row = by_id[bundle_id]
        bundle_qc = annotation_row.get("annotation_qc")
        if not isinstance(bundle_qc, dict):
            raise PortableEvidenceError("bundle annotation QC document is absent")
        validate_annotation_qc(bundle_qc)
        if (
            annotation_row.get("manifest_sha256") != summary_row.get("manifest_sha256")
            or annotation_row.get("metrics") != summary_row.get("metrics")
            or annotation_row.get("domain_rows") != summary_row.get("domain_rows")
            or bundle_qc.get(QC_IDENTITY_FIELD) != summary_row.get("annotation_qc_identity_sha256")
        ):
            raise PortableEvidenceError("annotation QC bundle binding differs")
        bundle_qc_documents.append(bundle_qc)
    aggregate_qc = annotation_document.get("aggregate_qc")
    if not isinstance(aggregate_qc, dict):
        raise PortableEvidenceError("project aggregate annotation QC is absent")
    validate_annotation_qc(aggregate_qc)
    if aggregate_qc != aggregate_annotation_qc(bundle_qc_documents):
        raise PortableEvidenceError("project aggregate annotation QC differs")
    outlier_analysis = annotation_document.get("outlier_analysis")
    if not isinstance(outlier_analysis, dict) or not isinstance(
        outlier_analysis.get("flags"), list
    ):
        raise PortableEvidenceError("project outlier analysis is absent")
    persisted_explanations: dict[tuple[str, str], dict[str, str]] = {}
    for flag in outlier_analysis["flags"]:
        if not isinstance(flag, dict):
            raise PortableEvidenceError("project outlier flag is invalid")
        key = flag.get("bundle_id"), flag.get("metric")
        if not all(isinstance(value, str) for value in key) or key in persisted_explanations:
            raise PortableEvidenceError("project outlier flag identity is invalid")
        persisted_explanations[key] = {
            "code": flag.get("explanation_code"),
            "explanation": flag.get("explanation"),
        }
    outlier_records = [
        {
            "bundle_id": bundle_id,
            "annotation_qc": annotation_by_id[bundle_id]["annotation_qc"],
        }
        for bundle_id, _ in specs
    ]
    if outlier_analysis != _outlier_analysis(outlier_records, persisted_explanations):
        raise PortableEvidenceError("project outlier analysis differs")
    metric_names = sorted(next(iter(by_id.values()))["metrics"])
    expected_aggregate_metrics = {
        name: sum(int(by_id[bundle_id]["metrics"][name]) for bundle_id, _ in specs)
        for name in metric_names
    }
    if (
        annotation_document.get("bundle_count") != len(specs)
        or annotation_document.get("aggregate_metrics") != expected_aggregate_metrics
    ):
        raise PortableEvidenceError("annotation QC project counters differ")

    for bundle_id, root in specs:
        row = by_id[bundle_id]
        bundle_report = _load_json(
            output / "reports" / _bundle_validation_relative(bundle_id),
            "bundle validation report",
        )
        bundle_report_identity = bundle_report.get("bundle_validation_identity_sha256")
        unsigned_bundle_report = {
            key: value
            for key, value in bundle_report.items()
            if key != "bundle_validation_identity_sha256"
        }
        if (
            bundle_report.get("schema_version") != PROJECT_RELEASE_SCHEMA_VERSION
            or bundle_report.get("status") != "PASS"
            or bundle_report.get("bundle_id") != bundle_id
            or bundle_report.get("bundle_manifest_sha256") != row.get("manifest_sha256")
            or bundle_report.get("source_bundle_sha256") != row.get("source_bundle_sha256")
            or bundle_report.get("normative_schema_sha256") != NORMATIVE_SCHEMA_SHA256
            or bundle_report.get("metrics") != row.get("metrics")
            or bundle_report.get("annotation_qc")
            != annotation_by_id[bundle_id].get("annotation_qc")
            or bundle_report.get("annotation_qc", {}).get(QC_IDENTITY_FIELD)
            != row.get("annotation_qc_identity_sha256")
            or bundle_report.get("validation") != row.get("validation")
            or bundle_report.get("query_regression") != row.get("query")
            or bundle_report.get("errors") != []
            or bundle_report_identity != row.get("bundle_validation_identity_sha256")
            or bundle_report_identity != _canonical_identity(unsigned_bundle_report)
        ):
            raise PortableEvidenceError("bundle validation report identity differs")
        if (
            (output / f"{bundle_id}.genome").is_symlink()
            or root != (output / f"{bundle_id}.genome").resolve()
            or sha256_file(root / "manifest.json") != row.get("manifest_sha256")
            or sha256_file(root / "schema.json") != NORMATIVE_SCHEMA_SHA256
        ):
            raise PortableEvidenceError("project bundle identity differs")
        if revalidate_bundles:
            fresh = _bundle_record(bundle_id, root, reference_fasta=reference)
            if any(
                fresh[key] != expected
                for key, expected in (
                    ("manifest_sha256", row.get("manifest_sha256")),
                    ("schema_sha256", row.get("schema_sha256")),
                    ("source_bundle_sha256", row.get("source_bundle_sha256")),
                    ("metrics", row.get("metrics")),
                    ("annotation_qc", annotation_by_id[bundle_id].get("annotation_qc")),
                    ("domain_rows", row.get("domain_rows")),
                    ("query", row.get("query")),
                    ("warnings", row.get("warnings")),
                    ("validation", row.get("validation")),
                )
            ):
                raise PortableEvidenceError("project bundle revalidation differs from report")
    qc_records = outlier_records
    if (output / "reports" / "annotation_qc.tsv").read_bytes() != _qc_tsv(
        qc_records, outlier_analysis
    ):
        raise PortableEvidenceError("annotation QC TSV differs from JSON")
    if (output / "reports" / "annotation_qc.html").read_bytes() != _qc_html(
        qc_records, outlier_analysis
    ):
        raise PortableEvidenceError("annotation QC HTML differs from JSON")
    expected_checksum = _checksum_payload(output, bundle_ids)
    checksum_path = output / "SHA256SUMS"
    if checksum_path.is_symlink() or not checksum_path.is_file():
        raise PortableEvidenceError("project SHA256SUMS is missing or unsafe")
    if checksum_path.read_bytes() != expected_checksum:
        raise PortableEvidenceError("project SHA256SUMS differs")
    for directory in (output / "reports", output / "manifests", output / "logs"):
        for payload in _directory_payloads(directory).values():
            if FORBIDDEN_PORTABLE_TEXT.search(payload.decode("utf-8")):
                raise PortableEvidenceError("published project metadata contains a host path")
    return {
        "schema_version": PROJECT_RELEASE_SCHEMA_VERSION,
        "status": "PASS",
        "bundle_count": len(specs),
        "validation_summary_identity_sha256": supplied_identity,
        "sha256sums_sha256": sha256_file(checksum_path),
        "bundles_revalidated": revalidate_bundles,
        "private_paths_persisted": False,
    }


def _parse_bundle(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("bundle must use ID=PATH")
    bundle_id, path = value.split("=", 1)
    if not bundle_id or not path:
        raise argparse.ArgumentTypeError("bundle ID and path must be non-empty")
    return bundle_id, Path(path)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build deterministic multi-bundle reports")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("build", "validate"):
        command = subparsers.add_parser(name)
        command.add_argument("--bundle", action="append", type=_parse_bundle, required=True)
        command.add_argument("--output-root", type=Path, required=True)
        command.add_argument("--reference-fasta", type=Path, required=True)
        if name == "build":
            command.add_argument("--generated-at", required=True)
            command.add_argument("--minimum-free-bytes", type=int, default=0)
            command.add_argument("--qc-explanations", type=Path)
    args = parser.parse_args()
    if args.command == "build":
        report = build_project_release(
            args.bundle,
            args.output_root,
            reference_fasta=args.reference_fasta,
            generated_at=args.generated_at,
            minimum_free_bytes=args.minimum_free_bytes,
            qc_explanations=_load_qc_explanations(args.qc_explanations),
        )
    else:
        report = validate_project_release(
            args.bundle,
            args.output_root,
            reference_fasta=args.reference_fasta,
        )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
