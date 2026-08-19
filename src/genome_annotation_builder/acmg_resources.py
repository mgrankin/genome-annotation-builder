"""Freeze and verify the current ACMG/ClinGen secondary-findings boundary."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import unicodedata
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

from openpyxl import load_workbook

from .config import BuildConfig
from .constants import (
    ACMG_SF_CANONICAL_SHA256,
    ACMG_SF_CLINGEN_ROW_COUNT,
    ACMG_SF_FREQUENCY_ROW_COUNT,
    ACMG_SF_GENE_COUNT,
    ACMG_SF_GENE_SET_SHA256,
    ACMG_SF_GUIDANCE_ROW_COUNT,
    ACMG_SF_POLICY_DOI,
    ACMG_SF_POLICY_PMID,
    ACMG_SF_POLICY_ROW_COUNT,
    ACMG_SF_SNAPSHOT,
    ACMG_SF_VERSION,
    ACMG_SF_WORKBOOK_BYTES,
    ACMG_SF_WORKBOOK_NAME,
    ACMG_SF_WORKBOOK_SHA256,
    ACMG_SF_WORKBOOK_URL,
)
from .errors import ResourceError
from .util import atomic_write_bytes, atomic_write_json, canonical_json_bytes, sha256_file

POLICY_HEADER = (
    "Gene",
    "Gene MIM",
    "Disease/Phentyope",
    "Disorder MIM",
    "Phenotype Category",
    "Inheritance",
    "SF List Version",
    "Variants to report",
)
GENE_DISEASE_HEADER = (
    "Gene Symbol",
    "HGNC ID",
    "MONDO Disease Name",
    "MONDO ID",
    "MOI",
    "Gene-Disease Validity",
    "Reportable as SF",
)
GUIDANCE_HEADER = ("Gene Symbol", "Reporting Guidance Comment")
VARIANT_GUIDANCE_HEADER = (
    "ClinVar Variant ID",
    "Gene Symbol",
    "Transcript",
    "c.HGVS",
    "p.HGVS",
    "Reportable as SF",
    "Rationale",
)
REQUIRED_SHEETS = {
    "Update Log",
    "Gene-Disease Reporting",
    "ACMG SFClinGen Reporting Guidan",
    "ACMG SF v3.3 List",
    "DAFTs, BA1, BS1 Specifications",
    "Variant-Specific Guidance",
    "Gene Comments",
}


def _vendor_path() -> Path:
    return (
        Path(__file__).parent
        / "vendor"
        / "acmg_sf"
        / f"acmg_sf_v{ACMG_SF_VERSION}_{ACMG_SF_SNAPSHOT}.json"
    )


def _resource_paths(config: BuildConfig) -> tuple[Path, Path, Path]:
    root = config.resource_cache / "acmg_sf" / ACMG_SF_VERSION
    canonical = getattr(
        config,
        "acmg_sf_dataset",
        root / f"acmg_sf_v{ACMG_SF_VERSION}_{ACMG_SF_SNAPSHOT}.json",
    )
    manifest = getattr(config, "acmg_sf_manifest", root / "resource_manifest.json")
    return canonical, manifest, manifest.parent / "archive" / ACMG_SF_WORKBOOK_NAME


def _clean_text(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ResourceError("ACMG workbook contains a non-finite number")
        return str(int(value)) if value.is_integer() else format(value, ".15g")
    text = unicodedata.normalize("NFC", str(value)).replace("\u00a0", " ")
    lines = [re.sub(r"[ \t]+", " ", line).strip() for line in text.splitlines()]
    while lines and not lines[-1]:
        lines.pop()
    return "\n".join(lines).strip() or None


def _identifier(value: Any, prefix: str) -> str:
    text = _clean_text(value)
    if text is None:
        raise ResourceError(f"missing {prefix} identifier in ACMG workbook")
    if text.startswith(prefix + ":"):
        return text
    if not re.fullmatch(r"\d+", text):
        raise ResourceError(f"invalid {prefix} identifier in ACMG workbook: {text!r}")
    return f"{prefix}:{text}"


def _identifiers(value: Any, prefix: str) -> list[str]:
    text = _clean_text(value)
    if text is None or text.lower() in {"n/a", "na"}:
        return []
    identifiers = re.findall(r"\d+", text)
    residual = re.sub(r"[\d,;\s]+", "", text)
    if not identifiers or residual:
        raise ResourceError(f"invalid {prefix} identifier list in ACMG workbook: {text!r}")
    return [f"{prefix}:{identifier}" for identifier in identifiers]


def _version(value: Any) -> str:
    text = _clean_text(value)
    if text is None:
        raise ResourceError("missing ACMG list version")
    try:
        return f"{float(text):.1f}"
    except ValueError as exc:
        raise ResourceError(f"invalid ACMG list version: {text!r}") from exc


def _header(sheet: Any, expected: tuple[str, ...]) -> None:
    row = next(sheet.iter_rows(min_row=1, max_row=1, values_only=True))
    observed = tuple(_clean_text(value) for value in row[: len(expected)])
    if observed != expected:
        raise ResourceError(f"unexpected header in ACMG sheet {sheet.title!r}: {observed!r}")
    if any(_clean_text(value) is not None for value in row[len(expected) :]):
        raise ResourceError(f"unexpected extra header columns in ACMG sheet {sheet.title!r}")


def _policy_rows(sheet: Any) -> list[dict[str, Any]]:
    _header(sheet, POLICY_HEADER)
    output: list[dict[str, Any]] = []
    for values in sheet.iter_rows(min_row=2, values_only=True):
        if not any(value is not None for value in values):
            continue
        if any(_clean_text(value) is not None for value in values[8:]):
            raise ResourceError("unexpected data beyond ACMG policy columns")
        gene, gene_mim, disease, disorder_mim, category, inheritance, version, variants = values[:8]
        output.append(
            {
                "gene_symbol": _clean_text(gene),
                "gene_mim": _identifier(gene_mim, "MIM"),
                "disease": _clean_text(disease),
                "disorder_mim_ids": _identifiers(disorder_mim, "MIM"),
                "disorder_mim_source": _clean_text(disorder_mim),
                "phenotype_category": _clean_text(category),
                "inheritance": _clean_text(inheritance),
                "introduced_version": _version(version),
                "variants_to_report": _clean_text(variants),
            }
        )
    return output


def _gene_disease_rows(sheet: Any) -> list[dict[str, Any]]:
    _header(sheet, GENE_DISEASE_HEADER)
    output: list[dict[str, Any]] = []
    for values in sheet.iter_rows(min_row=2, values_only=True):
        if not any(value is not None for value in values):
            continue
        if any(_clean_text(value) is not None for value in values[7:]):
            raise ResourceError("unexpected data beyond ClinGen gene-disease columns")
        gene, hgnc, disease, mondo, moi, validity, reportable = values[:7]
        output.append(
            {
                "gene_symbol": _clean_text(gene),
                "hgnc_id": _identifier(hgnc, "HGNC"),
                "disease": _clean_text(disease),
                "mondo_id": _clean_text(mondo),
                "mode_of_inheritance": _clean_text(moi),
                "gene_disease_validity": _clean_text(validity),
                "reportable_as_secondary_finding": _clean_text(reportable),
            }
        )
    return output


def _guidance_rows(sheet: Any) -> list[dict[str, Any]]:
    _header(sheet, GUIDANCE_HEADER)
    output: list[dict[str, Any]] = []
    for values in sheet.iter_rows(min_row=2, values_only=True):
        if not any(value is not None for value in values):
            continue
        if any(_clean_text(value) is not None for value in values[2:]):
            raise ResourceError("unexpected data beyond ClinGen reporting-guidance columns")
        output.append({"gene_symbol": _clean_text(values[0]), "comment": _clean_text(values[1])})
    return output


def _frequency_rows(sheet: Any) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    current_gene: str | None = None
    current_transcript: str | None = None
    for row_number, values in enumerate(sheet.iter_rows(min_row=4, values_only=True), start=4):
        values = values[:12]
        if row_number > 107:
            break
        if not any(value is not None for value in values):
            continue
        gene = _clean_text(values[0])
        transcript = _clean_text(values[1])
        if gene is not None:
            current_gene = gene
        if transcript is not None:
            current_transcript = transcript
        if current_gene is None or current_transcript is None:
            raise ResourceError("ACMG frequency row lacks inherited gene/transcript context")
        output.append(
            {
                "gene_symbol": current_gene,
                "mane_transcript": current_transcript,
                "disease": _clean_text(values[2]),
                "current_vcep_ba1": _clean_text(values[3]),
                "current_vcep_bs1": _clean_text(values[4]),
                "prevalence_one_in": _clean_text(values[5]),
                "gene_contribution": _clean_text(values[6]),
                "penetrance": _clean_text(values[7]),
                "calculated_ba1": _clean_text(values[8]),
                "calculated_bs1": _clean_text(values[9]),
                "notes": _clean_text(values[10]),
                "references": _clean_text(values[11]),
            }
        )
    return output


def _simple_rows(sheet: Any, width: int) -> list[list[str | None]]:
    output: list[list[str | None]] = []
    for values in sheet.iter_rows(values_only=True):
        if not any(value is not None for value in values):
            continue
        output.append([_clean_text(value) for value in values[:width]])
    return output


def _validate_snapshot(snapshot: dict[str, Any]) -> dict[str, int]:
    policy = snapshot.get("policy_gene_disease_pairs")
    gene_disease = snapshot.get("clingen_gene_disease_reporting")
    guidance = snapshot.get("clingen_reporting_guidance")
    frequency = snapshot.get("frequency_specifications")
    if not all(isinstance(value, list) for value in (policy, gene_disease, guidance, frequency)):
        raise ResourceError("ACMG canonical snapshot has invalid tables")
    assert isinstance(policy, list)
    assert isinstance(gene_disease, list)
    assert isinstance(guidance, list)
    assert isinstance(frequency, list)
    if len(policy) != ACMG_SF_POLICY_ROW_COUNT:
        raise ResourceError(f"ACMG policy row count is {len(policy)}, expected 100")
    if len(gene_disease) != ACMG_SF_CLINGEN_ROW_COUNT:
        raise ResourceError(f"ClinGen disease row count is {len(gene_disease)}, expected 177")
    if len(guidance) != ACMG_SF_GUIDANCE_ROW_COUNT:
        raise ResourceError(f"ClinGen guidance row count is {len(guidance)}, expected 84")
    if len(frequency) != ACMG_SF_FREQUENCY_ROW_COUNT:
        raise ResourceError(f"ClinGen frequency row count is {len(frequency)}, expected 104")
    policy_genes = sorted({str(row["gene_symbol"]) for row in policy})
    reporting_genes = sorted({str(row["gene_symbol"]) for row in gene_disease})
    guidance_genes = sorted({str(row["gene_symbol"]) for row in guidance})
    if policy_genes != reporting_genes or policy_genes != guidance_genes:
        raise ResourceError("ACMG policy, reporting, and guidance gene sets differ")
    if len(policy_genes) != ACMG_SF_GENE_COUNT:
        raise ResourceError(f"ACMG gene count is {len(policy_genes)}, expected 84")
    gene_digest = hashlib.sha256(("\n".join(policy_genes) + "\n").encode()).hexdigest()
    if gene_digest != ACMG_SF_GENE_SET_SHA256:
        raise ResourceError("ACMG gene-set digest mismatch")
    additions = {
        str(row["gene_symbol"]) for row in policy if row["introduced_version"] == ACMG_SF_VERSION
    }
    if additions != {"ABCD1", "CYP27A1", "PLN"}:
        raise ResourceError(f"unexpected ACMG SF v3.3 additions: {sorted(additions)}")
    allowed = {"Yes", "No", "Pending", "NA", "See Guidance"}
    observed = {str(row["reportable_as_secondary_finding"]) for row in gene_disease}
    if not observed <= allowed:
        raise ResourceError(f"unknown ClinGen reportability values: {sorted(observed - allowed)}")
    return {
        "genes": len(policy_genes),
        "policy_gene_disease_pairs": len(policy),
        "clingen_gene_disease_rows": len(gene_disease),
        "reporting_guidance_rows": len(guidance),
        "frequency_specification_rows": len(frequency),
    }


def extract_acmg_workbook(path: str | Path) -> dict[str, Any]:
    workbook_path = Path(path)
    if not workbook_path.is_file():
        raise ResourceError(f"ACMG workbook is missing: {workbook_path}")
    if (
        workbook_path.stat().st_size != ACMG_SF_WORKBOOK_BYTES
        or sha256_file(workbook_path) != ACMG_SF_WORKBOOK_SHA256
    ):
        raise ResourceError("ACMG workbook size/SHA-256 differs from the pinned snapshot")
    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    try:
        if not set(workbook.sheetnames) >= REQUIRED_SHEETS:
            raise ResourceError("ACMG workbook is missing required sheets")
        variant_rows = _simple_rows(workbook["Variant-Specific Guidance"], 7)
        if not variant_rows or tuple(variant_rows[0]) != VARIANT_GUIDANCE_HEADER:
            raise ResourceError("unexpected variant-specific guidance header")
        if variant_rows[1:] != [["Under construction", None, None, None, None, None, None]]:
            raise ResourceError("variant-specific guidance is no longer the pinned empty state")
        snapshot: dict[str, Any] = {
            "resource_type": "acmg_secondary_findings_policy_and_clingen_curation",
            "version": ACMG_SF_VERSION,
            "snapshot": ACMG_SF_SNAPSHOT,
            "policy": {
                "pmid": ACMG_SF_POLICY_PMID,
                "doi": ACMG_SF_POLICY_DOI,
                "title": "ACMG SF v3.3 list for reporting of secondary findings",
            },
            "source": {
                "url": ACMG_SF_WORKBOOK_URL,
                "workbook_sha256": ACMG_SF_WORKBOOK_SHA256,
                "workbook_bytes": ACMG_SF_WORKBOOK_BYTES,
                "accessed_date": "2026-08-15",
                "license": "CC0-1.0",
                "attribution": "Clinical Genome Resource (ClinGen), accessed 2026-08-15",
            },
            "policy_gene_disease_pairs": _policy_rows(workbook["ACMG SF v3.3 List"]),
            "clingen_gene_disease_reporting": _gene_disease_rows(
                workbook["Gene-Disease Reporting"]
            ),
            "clingen_reporting_guidance": _guidance_rows(
                workbook["ACMG SFClinGen Reporting Guidan"]
            ),
            "frequency_specifications": _frequency_rows(workbook["DAFTs, BA1, BS1 Specifications"]),
            "variant_specific_guidance": {"status": "under_construction", "rows": []},
            "update_log": _simple_rows(workbook["Update Log"], 2),
            "gene_comments": _simple_rows(workbook["Gene Comments"], 3),
        }
    finally:
        workbook.close()
    snapshot["counts"] = _validate_snapshot(snapshot)
    return snapshot


def verify_vendor_acmg_snapshot() -> dict[str, Any]:
    path = _vendor_path()
    try:
        payload = path.read_bytes()
        snapshot = json.loads(payload)
    except FileNotFoundError as exc:
        raise ResourceError(f"vendored ACMG snapshot is missing: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ResourceError(f"vendored ACMG snapshot is invalid JSON: {exc}") from exc
    if payload != canonical_json_bytes(snapshot):
        raise ResourceError("vendored ACMG snapshot is not canonical JSON")
    digest = hashlib.sha256(payload).hexdigest()
    if ACMG_SF_CANONICAL_SHA256 and digest != ACMG_SF_CANONICAL_SHA256:
        raise ResourceError("vendored ACMG canonical SHA-256 mismatch")
    _validate_snapshot(snapshot)
    if snapshot.get("version") != ACMG_SF_VERSION or snapshot.get("snapshot") != ACMG_SF_SNAPSHOT:
        raise ResourceError("vendored ACMG version/snapshot mismatch")
    return snapshot


def _download_workbook(config: BuildConfig, destination: Path) -> None:
    if destination.is_file():
        if (
            destination.stat().st_size != ACMG_SF_WORKBOOK_BYTES
            or sha256_file(destination) != ACMG_SF_WORKBOOK_SHA256
        ):
            raise ResourceError("cached ACMG workbook differs from the pinned source")
        return
    if not config.allow_download:
        return
    curl = shutil.which("curl")
    if curl is None:
        raise ResourceError("curl is required to archive the ACMG source workbook")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f".{destination.name}.partial")
    result = subprocess.run(
        [
            curl,
            "--fail",
            "--location",
            "--retry",
            "8",
            "--proto",
            "=https",
            "--tlsv1.2",
            "--output",
            str(partial),
            ACMG_SF_WORKBOOK_URL,
        ],
        check=False,
    )
    if result.returncode != 0:
        raise ResourceError(
            f"ACMG workbook download failed with exit {result.returncode}; retained {partial}"
        )
    if (
        partial.stat().st_size != ACMG_SF_WORKBOOK_BYTES
        or sha256_file(partial) != ACMG_SF_WORKBOOK_SHA256
    ):
        raise ResourceError("downloaded ACMG workbook is newer/different than the pinned snapshot")
    os.replace(partial, destination)


def verify_acmg_sf(config: BuildConfig) -> dict[str, Any]:
    canonical, manifest_path, _workbook = _resource_paths(config)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        snapshot = json.loads(canonical.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ResourceError("ACMG canonical snapshot or manifest is missing") from exc
    except json.JSONDecodeError as exc:
        raise ResourceError(f"invalid ACMG snapshot/manifest JSON: {exc}") from exc
    digest = sha256_file(canonical)
    if ACMG_SF_CANONICAL_SHA256 and digest != ACMG_SF_CANONICAL_SHA256:
        raise ResourceError("installed ACMG canonical SHA-256 is not pinned")
    if digest != manifest.get("canonical_sha256"):
        raise ResourceError("installed ACMG canonical SHA-256 differs from manifest")
    if manifest.get("version") != ACMG_SF_VERSION:
        raise ResourceError("ACMG manifest version is not pinned")
    if manifest.get("snapshot") != ACMG_SF_SNAPSHOT:
        raise ResourceError("ACMG manifest snapshot is not pinned")
    if manifest.get("source_workbook_sha256") != ACMG_SF_WORKBOOK_SHA256:
        raise ResourceError("ACMG source-workbook identity differs from manifest")
    counts = _validate_snapshot(snapshot)
    if manifest.get("counts") != counts:
        raise ResourceError("ACMG manifest counts differ from canonical data")
    if snapshot.get("version") != ACMG_SF_VERSION or snapshot.get("snapshot") != ACMG_SF_SNAPSHOT:
        raise ResourceError("installed ACMG version/snapshot mismatch")
    return manifest


def load_acmg_sf(config: BuildConfig) -> dict[str, Any]:
    """Return the installed canonical snapshot only after full verification."""

    verify_acmg_sf(config)
    canonical, _manifest_path, _workbook = _resource_paths(config)
    try:
        snapshot = json.loads(canonical.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ResourceError(f"cannot load verified ACMG snapshot: {exc}") from exc
    _validate_snapshot(snapshot)
    return snapshot


def sync_acmg_sf(config: BuildConfig) -> dict[str, Any]:
    canonical, manifest_path, workbook = _resource_paths(config)
    vendor_snapshot = verify_vendor_acmg_snapshot()
    vendor_payload = canonical_json_bytes(vendor_snapshot)
    if manifest_path.exists():
        return verify_acmg_sf(config)
    _download_workbook(config, workbook)
    if workbook.is_file():
        extracted = extract_acmg_workbook(workbook)
        if canonical_json_bytes(extracted) != vendor_payload:
            raise ResourceError("ACMG workbook extraction differs from the vendored snapshot")
    if canonical.exists():
        if canonical.read_bytes() != vendor_payload:
            raise ResourceError("existing ACMG canonical snapshot differs from vendor")
    else:
        atomic_write_bytes(canonical, vendor_payload)
    counts = _validate_snapshot(vendor_snapshot)
    manifest = {
        "resource_type": "acmg_secondary_findings_policy_and_clingen_curation",
        "version": ACMG_SF_VERSION,
        "snapshot": ACMG_SF_SNAPSHOT,
        "canonical_sha256": hashlib.sha256(vendor_payload).hexdigest(),
        "canonical_bytes": len(vendor_payload),
        "source_workbook_url": ACMG_SF_WORKBOOK_URL,
        "source_workbook_sha256": ACMG_SF_WORKBOOK_SHA256,
        "source_workbook_bytes": ACMG_SF_WORKBOOK_BYTES,
        "source_archive_present": workbook.is_file(),
        "policy_pmid": ACMG_SF_POLICY_PMID,
        "policy_doi": ACMG_SF_POLICY_DOI,
        "license": "CC0-1.0",
        "terms_url": "https://www.clinicalgenome.org/docs/terms-of-use/",
        "attribution": "Clinical Genome Resource (ClinGen), accessed 2026-08-15",
        "installed_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "counts": counts,
    }
    atomic_write_json(manifest_path, manifest)
    return verify_acmg_sf(config)
