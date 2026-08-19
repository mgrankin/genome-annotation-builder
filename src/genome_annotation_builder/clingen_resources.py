"""Materialize and verify the vendored CC0 ClinGen gene-level snapshot."""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .arrow_schemas import (
    clingen_dosage_sensitivity_schema,
    clingen_gene_disease_schema,
)
from .clingen_freeze import load_canonical_snapshot, snapshot_identity
from .config import BuildConfig
from .constants import (
    CLINGEN_ATTRIBUTION,
    CLINGEN_DOSAGE_BYTES,
    CLINGEN_DOSAGE_ROWS,
    CLINGEN_DOSAGE_SHA256,
    CLINGEN_DOSAGE_SNAPSHOT,
    CLINGEN_DOSAGE_URL,
    CLINGEN_GENE_DISEASE_BYTES,
    CLINGEN_GENE_DISEASE_ROWS,
    CLINGEN_GENE_DISEASE_SHA256,
    CLINGEN_GENE_DISEASE_SNAPSHOT,
    CLINGEN_GENE_DISEASE_URL,
    CLINGEN_GENE_LEVEL_CANONICAL_BYTES,
    CLINGEN_GENE_LEVEL_CANONICAL_SHA256,
    CLINGEN_GENE_LEVEL_FILE,
    CLINGEN_GENE_LEVEL_SNAPSHOT,
    CLINGEN_LICENSE,
)
from .errors import ResourceError
from .util import atomic_output_directory, atomic_write_json, sha256_file

RESOURCE_SCHEMA_VERSION = 1
SEMANTIC_BOUNDARY = (
    "Gene-level validity and dosage context only; presence of a variant in a curated gene "
    "is not evidence that the variant is pathogenic or that the subject has the associated "
    "condition."
)


def vendor_clingen_path() -> Path:
    return Path(__file__).resolve().parent / "vendor" / "clingen" / CLINGEN_GENE_LEVEL_FILE


def verify_vendor_clingen_snapshot() -> dict[str, Any]:
    path = vendor_clingen_path()
    if not path.is_file():
        raise ResourceError(f"vendored ClinGen snapshot not found: {path}")
    actual_bytes = path.stat().st_size
    if actual_bytes != CLINGEN_GENE_LEVEL_CANONICAL_BYTES:
        raise ResourceError(
            "vendored ClinGen snapshot byte count mismatch: "
            f"{actual_bytes} != {CLINGEN_GENE_LEVEL_CANONICAL_BYTES}"
        )
    actual_sha256 = sha256_file(path)
    if actual_sha256 != CLINGEN_GENE_LEVEL_CANONICAL_SHA256:
        raise ResourceError("vendored ClinGen snapshot SHA-256 mismatch")
    document = load_canonical_snapshot(path)
    identity = snapshot_identity(document)
    if identity != {
        "snapshot": CLINGEN_GENE_LEVEL_SNAPSHOT,
        "gene_disease_rows": CLINGEN_GENE_DISEASE_ROWS,
        "dosage_rows": CLINGEN_DOSAGE_ROWS,
        "canonical_sha256": CLINGEN_GENE_LEVEL_CANONICAL_SHA256,
    }:
        raise ResourceError(f"vendored ClinGen snapshot identity mismatch: {identity}")
    if document.get("license") != CLINGEN_LICENSE:
        raise ResourceError("vendored ClinGen licence boundary changed")
    if document.get("attribution") != CLINGEN_ATTRIBUTION:
        raise ResourceError("vendored ClinGen attribution changed")
    if document.get("semantic_boundary") != SEMANTIC_BOUNDARY:
        raise ResourceError("vendored ClinGen semantic boundary changed")
    expected_sources = {
        "gene_disease_validity": {
            "url": CLINGEN_GENE_DISEASE_URL,
            "snapshot": CLINGEN_GENE_DISEASE_SNAPSHOT,
            "bytes": CLINGEN_GENE_DISEASE_BYTES,
            "sha256": CLINGEN_GENE_DISEASE_SHA256,
            "row_count": CLINGEN_GENE_DISEASE_ROWS,
        },
        "dosage_sensitivity": {
            "url": CLINGEN_DOSAGE_URL,
            "snapshot": CLINGEN_DOSAGE_SNAPSHOT,
            "bytes": CLINGEN_DOSAGE_BYTES,
            "sha256": CLINGEN_DOSAGE_SHA256,
            "row_count": CLINGEN_DOSAGE_ROWS,
        },
    }
    sources = document.get("sources")
    if not isinstance(sources, dict):
        raise ResourceError("vendored ClinGen source metadata is invalid")
    for name, expected in expected_sources.items():
        actual = sources.get(name)
        if not isinstance(actual, dict) or any(
            actual.get(key) != value for key, value in expected.items()
        ):
            raise ResourceError(f"vendored ClinGen {name} source identity changed")
    return {
        **identity,
        "canonical_bytes": actual_bytes,
        "canonical_path": str(path),
        "sources": sources,
        "license": CLINGEN_LICENSE,
        "attribution": CLINGEN_ATTRIBUTION,
        "semantic_boundary": SEMANTIC_BOUNDARY,
    }


def _gene_disease_rows(document: dict[str, Any]) -> list[dict[str, Any]]:
    output = []
    for row in document["gene_disease_validity"]:
        source_datetime = str(row["classification_date"])
        output.append(
            {
                "evidence_id": "clingen-gene-disease:" + row["assertion_id"],
                "assertion_id": row["assertion_id"],
                "gene_symbol": row["gene_symbol"],
                "hgnc_id": row["hgnc_id"],
                "disease_label": row["disease_label"],
                "mondo_id": row["mondo_id"],
                "mode_of_inheritance": row["mode_of_inheritance"],
                "sop": row["sop"],
                "classification": row["classification"],
                "online_report": row["online_report"],
                "classification_date": date.fromisoformat(source_datetime[:10]),
                "classification_datetime_source": source_datetime,
                "gcep": row["gcep"],
                "snapshot": CLINGEN_GENE_LEVEL_SNAPSHOT,
                "source_row_number": row["source_row_number"],
                "source_row_sha256": row["source_row_sha256"],
            }
        )
    return output


def _dosage_rows(document: dict[str, Any]) -> list[dict[str, Any]]:
    output = []
    for row in document["dosage_sensitivity"]:
        output.append(
            {
                "evidence_id": "clingen-dosage:" + row["gene_symbol"],
                "gene_symbol": row["gene_symbol"],
                "ncbi_gene_id": row["ncbi_gene_id"],
                "cytoband": row["cytoband"],
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "haploinsufficiency_score": row["haploinsufficiency_score"],
                "haploinsufficiency_source_score": row["haploinsufficiency_source_score"],
                "haploinsufficiency_description": row["haploinsufficiency_description"],
                "haploinsufficiency_pmids": row["haploinsufficiency_pmids"],
                "haploinsufficiency_pmid_source_values": row[
                    "haploinsufficiency_pmid_source_values"
                ],
                "haploinsufficiency_disease_id": row["haploinsufficiency_disease_id"],
                "triplosensitivity_score": row["triplosensitivity_score"],
                "triplosensitivity_source_score": row["triplosensitivity_source_score"],
                "triplosensitivity_description": row["triplosensitivity_description"],
                "triplosensitivity_pmids": row["triplosensitivity_pmids"],
                "triplosensitivity_pmid_source_values": row["triplosensitivity_pmid_source_values"],
                "triplosensitivity_disease_id": row["triplosensitivity_disease_id"],
                "last_evaluated": date.fromisoformat(row["last_evaluated"]),
                "snapshot": CLINGEN_GENE_LEVEL_SNAPSHOT,
                "source_row_number": row["source_row_number"],
                "source_row_sha256": row["source_row_sha256"],
            }
        )
    return output


def _write_parquet(rows: list[dict[str, Any]], schema: pa.Schema, path: Path) -> None:
    table = pa.Table.from_pylist(rows, schema=schema)
    pq.write_table(
        table,
        path,
        compression="zstd",
        compression_level=9,
        use_dictionary=True,
        write_statistics=True,
        version="2.6",
    )


def _relative(config: BuildConfig, path: Path, staging: Path) -> Path:
    try:
        relative = path.resolve().relative_to(config.clingen.root.resolve())
    except ValueError as exc:
        raise ResourceError(f"ClinGen output is outside configured root: {path}") from exc
    return staging / relative


def sync_clingen_resources(config: BuildConfig) -> dict[str, Any]:
    if not config.clingen.enabled:
        return {"status": "not_enabled"}
    vendor = verify_vendor_clingen_snapshot()
    if config.clingen.root.exists():
        return verify_clingen_resources(config)
    document = load_canonical_snapshot(vendor_clingen_path())
    with atomic_output_directory(config.clingen.root) as staging:
        gene_path = _relative(config, config.clingen.gene_disease_parquet, staging)
        dosage_path = _relative(config, config.clingen.dosage_parquet, staging)
        manifest_path = _relative(config, config.clingen.manifest, staging)
        gene_path.parent.mkdir(parents=True, exist_ok=True)
        dosage_path.parent.mkdir(parents=True, exist_ok=True)
        _write_parquet(_gene_disease_rows(document), clingen_gene_disease_schema(), gene_path)
        _write_parquet(_dosage_rows(document), clingen_dosage_sensitivity_schema(), dosage_path)
        manifest = {
            "schema_version": RESOURCE_SCHEMA_VERSION,
            "resource_type": "clingen_gene_level_curations",
            "snapshot": CLINGEN_GENE_LEVEL_SNAPSHOT,
            "license": CLINGEN_LICENSE,
            "redistribution_allowed": True,
            "attribution": CLINGEN_ATTRIBUTION,
            "semantic_boundary": SEMANTIC_BOUNDARY,
            "canonical": {
                "file": CLINGEN_GENE_LEVEL_FILE,
                "bytes": vendor["canonical_bytes"],
                "sha256": vendor["canonical_sha256"],
            },
            "sources": vendor["sources"],
            "outputs": {
                "gene_disease_validity": {
                    "path": config.clingen.gene_disease_parquet.resolve()
                    .relative_to(config.clingen.root.resolve())
                    .as_posix(),
                    "rows": CLINGEN_GENE_DISEASE_ROWS,
                    "sha256": sha256_file(gene_path),
                },
                "dosage_sensitivity": {
                    "path": config.clingen.dosage_parquet.resolve()
                    .relative_to(config.clingen.root.resolve())
                    .as_posix(),
                    "rows": CLINGEN_DOSAGE_ROWS,
                    "sha256": sha256_file(dosage_path),
                },
            },
        }
        atomic_write_json(manifest_path, manifest)
    return verify_clingen_resources(config)


def _verify_table(
    path: Path,
    expected_schema: pa.Schema,
    expected_rows: int,
    key: str,
) -> dict[str, Any]:
    if not path.is_file():
        raise ResourceError(f"ClinGen normalized table not found: {path}")
    try:
        table = pq.read_table(path)
    except (OSError, pa.ArrowException) as exc:
        raise ResourceError(f"cannot read ClinGen normalized table {path}: {exc}") from exc
    if not table.schema.equals(expected_schema, check_metadata=False):
        raise ResourceError(f"ClinGen normalized table schema mismatch: {path}")
    if table.num_rows != expected_rows:
        raise ResourceError(
            f"ClinGen normalized table row count mismatch: {table.num_rows} != {expected_rows}"
        )
    values = table[key]
    if values.null_count:
        raise ResourceError(f"ClinGen normalized table has null {key} values")
    unique = pc.count_distinct(values).as_py()
    if unique != expected_rows:
        raise ResourceError(f"ClinGen normalized table has duplicate {key} values")
    return {"rows": table.num_rows, "sha256": sha256_file(path)}


def verify_clingen_resources(config: BuildConfig) -> dict[str, Any]:
    if not config.clingen.enabled:
        return {"status": "not_enabled"}
    vendor = verify_vendor_clingen_snapshot()
    try:
        manifest = json.loads(config.clingen.manifest.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ResourceError(
            f"ClinGen resource manifest not found: {config.clingen.manifest}"
        ) from exc
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ResourceError(f"invalid ClinGen resource manifest: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ResourceError("ClinGen resource manifest must be a JSON object")
    expected_keys = {
        "schema_version",
        "resource_type",
        "snapshot",
        "license",
        "redistribution_allowed",
        "attribution",
        "semantic_boundary",
        "canonical",
        "sources",
        "outputs",
    }
    if set(manifest) != expected_keys:
        raise ResourceError("ClinGen resource manifest has an unexpected field set")
    expected_top = {
        "schema_version": RESOURCE_SCHEMA_VERSION,
        "resource_type": "clingen_gene_level_curations",
        "snapshot": CLINGEN_GENE_LEVEL_SNAPSHOT,
        "license": CLINGEN_LICENSE,
        "redistribution_allowed": True,
        "attribution": CLINGEN_ATTRIBUTION,
        "semantic_boundary": SEMANTIC_BOUNDARY,
    }
    if any(manifest.get(key) != value for key, value in expected_top.items()):
        raise ResourceError("ClinGen resource manifest boundary mismatch")
    expected_canonical = {
        "file": CLINGEN_GENE_LEVEL_FILE,
        "bytes": CLINGEN_GENE_LEVEL_CANONICAL_BYTES,
        "sha256": CLINGEN_GENE_LEVEL_CANONICAL_SHA256,
    }
    if manifest.get("canonical") != expected_canonical:
        raise ResourceError("ClinGen resource manifest canonical identity mismatch")
    if manifest.get("sources") != vendor["sources"]:
        raise ResourceError("ClinGen resource manifest source identities changed")
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict):
        raise ResourceError("ClinGen resource manifest outputs are invalid")
    gene = _verify_table(
        config.clingen.gene_disease_parquet,
        clingen_gene_disease_schema(),
        CLINGEN_GENE_DISEASE_ROWS,
        "evidence_id",
    )
    dosage = _verify_table(
        config.clingen.dosage_parquet,
        clingen_dosage_sensitivity_schema(),
        CLINGEN_DOSAGE_ROWS,
        "evidence_id",
    )
    expected_outputs = {
        "gene_disease_validity": {
            "path": config.clingen.gene_disease_parquet.resolve()
            .relative_to(config.clingen.root.resolve())
            .as_posix(),
            **gene,
        },
        "dosage_sensitivity": {
            "path": config.clingen.dosage_parquet.resolve()
            .relative_to(config.clingen.root.resolve())
            .as_posix(),
            **dosage,
        },
    }
    if outputs != expected_outputs:
        raise ResourceError("ClinGen normalized outputs differ from manifest")
    return {
        "status": "verified",
        "snapshot": CLINGEN_GENE_LEVEL_SNAPSHOT,
        "manifest_sha256": sha256_file(config.clingen.manifest),
        "canonical_sha256": CLINGEN_GENE_LEVEL_CANONICAL_SHA256,
        "license": CLINGEN_LICENSE,
        "redistribution_allowed": True,
        "attribution": CLINGEN_ATTRIBUTION,
        "semantic_boundary": SEMANTIC_BOUNDARY,
        "sources": vendor["sources"],
        "outputs": expected_outputs,
    }
