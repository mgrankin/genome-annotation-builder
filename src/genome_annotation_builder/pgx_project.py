"""Project PharmCAT/CPIC output into normative and lossless PGx tables."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .arrow_schemas import (
    pgx_gene_call_audit_schema,
    pgx_guidance_audit_schema,
    pharmacogenomics_schema,
)
from .config import BuildConfig
from .constants import PHARMCAT_VERSION, STRUCTURALLY_COMPLEX_PGX_GENES
from .errors import InputValidationError
from .pgx_resources import CPIC_NORMATIVE_LEVELS
from .util import atomic_output_directory, atomic_write_json, canonical_json_bytes, sha256_file

CALLS_HEADER = (
    "Gene",
    "Source Diplotype",
    "Phenotype",
    "Activity Score",
    "Haplotype 1",
    "Haplotype 1 Function",
    "Haplotype 1 Activity Value",
    "Haplotype 2",
    "Haplotype 2 Function",
    "Haplotype 2 Activity Value",
    "Outside Call",
    "Match Score",
    "Missing positions",
    "Recommendation Lookup Diplotype",
    "Recommendation Lookup Phenotype",
    "Recommendation Lookup Activity Score",
)
LEVEL_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3}


def _parse_calls_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise InputValidationError("PharmCAT calls-only report is not valid UTF-8 TSV") from exc
    if len(rows) < 3 or rows[0] != [f"PharmCAT {PHARMCAT_VERSION}"]:
        raise InputValidationError("PharmCAT calls-only report version line mismatch")
    if tuple(rows[1]) != CALLS_HEADER:
        raise InputValidationError("PharmCAT calls-only report header boundary mismatch")
    output: list[dict[str, str]] = []
    seen: set[str] = set()
    for row_number, values in enumerate(rows[2:], start=3):
        if not values or all(not value for value in values):
            continue
        if len(values) != len(CALLS_HEADER):
            raise InputValidationError(
                f"PharmCAT calls-only row {row_number} has {len(values)} columns"
            )
        row = {key: value.strip() for key, value in zip(CALLS_HEADER, values, strict=True)}
        gene = row["Gene"]
        if not gene or not re.fullmatch(r"[A-Za-z0-9-]+", gene):
            raise InputValidationError(f"PharmCAT calls-only row {row_number} has invalid gene")
        if gene in seen:
            raise InputValidationError(f"duplicate PharmCAT calls-only gene: {gene}")
        seen.add(gene)
        if row["Outside Call"] not in {"yes", "no"}:
            raise InputValidationError(
                f"PharmCAT calls-only row {row_number} has invalid Outside Call"
            )
        output.append(row)
    if not output:
        raise InputValidationError("PharmCAT calls-only report contains no gene rows")
    return output


def _activity_score(value: str) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _missing_positions(value: str) -> list[str]:
    if value.casefold() == "no":
        return []
    if not value:
        return ["UNSPECIFIED_BY_PHARMCAT"]
    return sorted({token.strip() for token in re.split(r"[;,]", value) if token.strip()})


def _source_diplotypes(value: str) -> list[str]:
    if not value:
        return []
    return [token.strip() for token in value.split(" OR ") if token.strip()]


def _load_cpic_pairs(path: Path) -> list[dict[str, Any]]:
    rows = pq.read_table(path).to_pylist()
    if not rows:
        raise InputValidationError("prepared CPIC pair table is empty")
    return rows


def _strongest_level(levels: set[str]) -> str | None:
    exact = sorted((level for level in levels if level in LEVEL_ORDER), key=LEVEL_ORDER.get)
    return exact[0] if exact else None


def _normative_pairs(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        row
        for row in rows
        if row.get("used_for_recommendation") is True
        and row.get("provisional") is False
        and row.get("cpic_level") in CPIC_NORMATIVE_LEVELS
    ]


def _project_gene_calls(
    calls: list[dict[str, str]],
    pairs: list[dict[str, Any]],
    *,
    calls_sha256: str,
    report_sha256: str,
    data_version: str,
    required_complex_genes: set[str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    pairs_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for pair in _normative_pairs(pairs):
        pairs_by_gene[str(pair["gene_symbol"]).upper()].append(pair)

    core_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    for call in sorted(calls, key=lambda item: item["Gene"]):
        gene = call["Gene"]
        source_diplotypes = _source_diplotypes(call["Source Diplotype"])
        phenotype_values = []
        for value in (call["Phenotype"], call["Recommendation Lookup Phenotype"]):
            if value and value not in phenotype_values:
                phenotype_values.append(value)
        selected_phenotype = phenotype_values[0] if phenotype_values else None
        selected_diplotype = source_diplotypes[0] if len(source_diplotypes) == 1 else None
        missing = _missing_positions(call["Missing positions"])
        outside_call = call["Outside Call"] == "yes"
        complex_gene = gene in STRUCTURALLY_COMPLEX_PGX_GENES
        gene_pairs = sorted(pairs_by_gene.get(gene.upper(), []), key=lambda row: row["pair_id"])
        levels = {str(row["cpic_level"]) for row in gene_pairs}
        selected_level = _strongest_level(levels)
        all_drugs = sorted({str(row["drug_name"]) for row in gene_pairs})
        all_urls = sorted(
            {str(row["guideline_url"]) for row in gene_pairs if row.get("guideline_url")}
        )

        if complex_gene:
            status = "complex_requires_validated_outside_call"
            reason = (
                "structurally complex PGx genes are withheld until a separately validated "
                "CNV/SV/HLA-aware outside-call adapter is configured"
            )
        elif len(source_diplotypes) != 1 or " OR " in call["Source Diplotype"]:
            status = "ambiguous_diplotype"
            reason = "PharmCAT returned zero or multiple possible source diplotypes"
        elif missing:
            status = "incomplete_positions"
            reason = "PharmCAT reports missing required positions"
        elif selected_phenotype is None:
            status = "phenotype_unavailable"
            reason = "PharmCAT did not assign a gene-level phenotype"
        elif selected_level is None:
            status = "no_normative_cpic_level"
            reason = "the frozen CPIC snapshot has no non-provisional Yes row at exact level A-D"
        else:
            status = "published"
            reason = (
                "one complete unambiguous PharmCAT call maps to the strongest exact CPIC "
                "recommendation level; lower-level drug mappings remain in the audit"
            )

        activity_raw = call["Activity Score"] or call["Recommendation Lookup Activity Score"]
        activity = _activity_score(activity_raw)
        audit_rows.append(
            {
                "gene_symbol": gene,
                "source_diplotypes": source_diplotypes,
                "selected_diplotype": selected_diplotype,
                "phenotypes": phenotype_values,
                "selected_phenotype": selected_phenotype,
                "activity_score_raw": activity_raw or None,
                "activity_score": activity,
                "outside_call": outside_call,
                "match_score": call["Match Score"] or None,
                "missing_positions": missing,
                "structurally_complex": complex_gene,
                "call_status": status,
                "cpic_levels": sorted(levels, key=lambda value: LEVEL_ORDER[value]),
                "affected_drugs": all_drugs,
                "guideline_urls": all_urls,
                "caller": "PharmCAT",
                "caller_version": PHARMCAT_VERSION,
                "data_version": data_version,
                "calls_tsv_sha256": calls_sha256,
                "report_json_sha256": report_sha256,
                "reason": reason,
            }
        )
        if status == "published":
            assert selected_diplotype is not None
            assert selected_phenotype is not None
            assert selected_level is not None
            selected_pairs = [row for row in gene_pairs if row.get("cpic_level") == selected_level]
            drugs = sorted({str(row["drug_name"]) for row in selected_pairs})
            urls = sorted(
                {str(row["guideline_url"]) for row in selected_pairs if row.get("guideline_url")}
            )
            core_rows.append(
                {
                    "gene_symbol": gene,
                    "diplotype": selected_diplotype,
                    "phenotype": selected_phenotype,
                    "activity_score": activity,
                    "copy_number": None,
                    "cpic_level": selected_level,
                    "affected_drugs": drugs,
                    "guideline_url": urls[0] if len(urls) == 1 else None,
                    "caller": "PharmCAT",
                    "caller_version": PHARMCAT_VERSION,
                }
            )
    observed_genes = {str(row["gene_symbol"]).upper() for row in audit_rows}
    for gene in sorted((required_complex_genes or set()) - observed_genes):
        gene_pairs = sorted(pairs_by_gene.get(gene.upper(), []), key=lambda row: row["pair_id"])
        levels = {str(row["cpic_level"]) for row in gene_pairs}
        audit_rows.append(
            {
                "gene_symbol": gene,
                "source_diplotypes": [],
                "selected_diplotype": None,
                "phenotypes": [],
                "selected_phenotype": None,
                "activity_score_raw": None,
                "activity_score": None,
                "outside_call": False,
                "match_score": None,
                "missing_positions": [],
                "structurally_complex": True,
                "call_status": "complex_requires_validated_outside_call",
                "cpic_levels": sorted(levels, key=lambda value: LEVEL_ORDER[value]),
                "affected_drugs": sorted({str(row["drug_name"]) for row in gene_pairs}),
                "guideline_urls": sorted(
                    {str(row["guideline_url"]) for row in gene_pairs if row.get("guideline_url")}
                ),
                "caller": "PharmCAT",
                "caller_version": PHARMCAT_VERSION,
                "data_version": data_version,
                "calls_tsv_sha256": calls_sha256,
                "report_json_sha256": report_sha256,
                "reason": (
                    "the gene is intentionally excluded from VCF-only matching and requires a "
                    "separately validated CNV/SV/HLA-aware outside call"
                ),
            }
        )
    audit_rows.sort(key=lambda row: str(row["gene_symbol"]))
    return core_rows, audit_rows


def _annotation_genes(annotation: dict[str, Any]) -> set[str]:
    genes: set[str] = set()
    lookup = annotation.get("lookupKey")
    if isinstance(lookup, list):
        for item in lookup:
            if isinstance(item, dict):
                genes.update(str(key) for key in item if str(key))
    genotypes = annotation.get("genotypes")
    if isinstance(genotypes, list):
        for genotype in genotypes:
            if not isinstance(genotype, dict):
                continue
            diplotypes = genotype.get("diplotypes")
            if not isinstance(diplotypes, list):
                continue
            for diplotype in diplotypes:
                if isinstance(diplotype, dict) and diplotype.get("gene"):
                    genes.add(str(diplotype["gene"]))
    return genes


def _project_guidance(
    report: dict[str, Any],
    pairs: list[dict[str, Any]],
    audit_rows: list[dict[str, Any]],
    report_sha256: str,
) -> list[dict[str, Any]]:
    pairs_by_gene_drug: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    genes_by_drug: dict[str, set[str]] = defaultdict(set)
    for pair in _normative_pairs(pairs):
        key = (str(pair["gene_symbol"]).upper(), str(pair["drug_name"]).casefold())
        pairs_by_gene_drug[key].append(pair)
        genes_by_drug[key[1]].add(key[0])
    published = {
        str(row["gene_symbol"]).upper() for row in audit_rows if row["call_status"] == "published"
    }
    output: list[dict[str, Any]] = []
    drug_groups = report.get("drugs")
    if not isinstance(drug_groups, dict):
        raise InputValidationError("PharmCAT report drugs field must be an object")
    for group_name in sorted(drug_groups):
        drugs = drug_groups[group_name]
        if not isinstance(drugs, dict):
            raise InputValidationError("PharmCAT report drug group must be an object")
        for drug_key in sorted(drugs):
            drug = drugs[drug_key]
            if not isinstance(drug, dict):
                raise InputValidationError("PharmCAT report drug entry must be an object")
            drug_name = str(drug.get("name") or drug_key)
            source = str(drug.get("source") or group_name)
            guidelines = drug.get("guidelines") or []
            if not isinstance(guidelines, list):
                raise InputValidationError("PharmCAT report guidelines must be a list")
            for guideline_index, guideline in enumerate(guidelines):
                if not isinstance(guideline, dict):
                    raise InputValidationError("PharmCAT report guideline must be an object")
                annotations = guideline.get("annotations") or []
                if not isinstance(annotations, list):
                    raise InputValidationError("PharmCAT report annotations must be a list")
                for annotation_index, annotation in enumerate(annotations):
                    if not isinstance(annotation, dict):
                        raise InputValidationError("PharmCAT report annotation must be an object")
                    genes = _annotation_genes(annotation)
                    if not genes:
                        genes = set(genes_by_drug.get(drug_name.casefold(), set()))
                    for gene in sorted(genes):
                        pair_rows = pairs_by_gene_drug.get((gene.upper(), drug_name.casefold()), [])
                        level = _strongest_level({str(row["cpic_level"]) for row in pair_rows})
                        identity = {
                            "source": source,
                            "gene": gene,
                            "drug": drug_name,
                            "guideline": guideline.get("id"),
                            "guideline_index": guideline_index,
                            "annotation_index": annotation_index,
                            "recommendation": annotation.get("drugRecommendation"),
                            "implications": annotation.get("implications") or [],
                        }
                        guidance_id = (
                            "pgx-guidance:"
                            + hashlib.sha256(canonical_json_bytes(identity)).hexdigest()
                        )
                        implications = annotation.get("implications") or []
                        if not isinstance(implications, list):
                            raise InputValidationError(
                                "PharmCAT report annotation implications must be a list"
                            )
                        output.append(
                            {
                                "guidance_id": guidance_id,
                                "source": source,
                                "gene_symbol": gene,
                                "drug": drug_name,
                                "guideline_id": (
                                    str(guideline["id"])
                                    if guideline.get("id") is not None
                                    else None
                                ),
                                "guideline_url": (
                                    str(guideline["url"])
                                    if guideline.get("url") is not None
                                    else None
                                ),
                                "cpic_level": level,
                                "recommendation": (
                                    str(annotation["drugRecommendation"])
                                    if annotation.get("drugRecommendation") is not None
                                    else None
                                ),
                                "implications": [str(value) for value in implications],
                                "classification": (
                                    str(annotation["classification"])
                                    if annotation.get("classification") is not None
                                    else None
                                ),
                                "population": (
                                    str(annotation["population"])
                                    if annotation.get("population") is not None
                                    else None
                                ),
                                "matched": gene.upper() in published,
                                "report_json_sha256": report_sha256,
                            }
                        )
    output.sort(
        key=lambda row: (
            row["source"],
            row["gene_symbol"],
            row["drug"].casefold(),
            row["guidance_id"],
        )
    )
    seen: set[str] = set()
    for row in output:
        if row["guidance_id"] in seen:
            raise InputValidationError("duplicate PharmCAT guidance identity")
        seen.add(row["guidance_id"])
    return output


def _write_table(path: Path, rows: list[dict[str, Any]], schema: pa.Schema) -> None:
    pq.write_table(
        pa.Table.from_pylist(rows, schema=schema),
        path,
        compression="zstd",
        compression_level=9,
        use_dictionary=True,
        write_statistics=True,
        version="2.6",
    )


def project_pgx_results(
    config: BuildConfig,
    *,
    pharmcat_dir: str | Path,
    resource_dir: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    """Create .genome PGx core rows and audit extensions."""

    if not config.pgx.enabled:
        raise InputValidationError("PGx projection was invoked while PGx is disabled")
    pharmcat_dir = Path(pharmcat_dir)
    resource_dir = Path(resource_dir)
    destination = Path(destination)
    calls_path = pharmcat_dir / "pharmcat.calls.tsv"
    report_path = pharmcat_dir / "pharmcat.report.canonical.json"
    if not calls_path.is_file() or not report_path.is_file():
        raise InputValidationError("stable PharmCAT calls/report artifacts are missing")
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise InputValidationError("canonical PharmCAT report is invalid") from exc
    if not isinstance(report, dict) or str(report.get("pharmcatVersion")) != PHARMCAT_VERSION:
        raise InputValidationError("canonical PharmCAT report version mismatch")
    data_version = str(report.get("dataVersion") or "")
    if not data_version:
        raise InputValidationError("canonical PharmCAT report lacks dataVersion")
    calls_sha = sha256_file(calls_path)
    report_sha = sha256_file(report_path)
    calls = _parse_calls_tsv(calls_path)
    pairs = _load_cpic_pairs(resource_dir / "cpic_pairs.parquet")
    locus_rows = pq.read_table(
        resource_dir / "pgx_loci.parquet", columns=["genes", "structurally_complex"]
    ).to_pylist()
    required_complex_genes = {
        str(gene).upper()
        for row in locus_rows
        if row["structurally_complex"]
        for gene in (row["genes"] or [])
        if str(gene).upper() in STRUCTURALLY_COMPLEX_PGX_GENES
    }
    core_rows, audit_rows = _project_gene_calls(
        calls,
        pairs,
        calls_sha256=calls_sha,
        report_sha256=report_sha,
        data_version=data_version,
        required_complex_genes=required_complex_genes,
    )
    guidance_rows = _project_guidance(report, pairs, audit_rows, report_sha)

    with atomic_output_directory(destination) as staging:
        gene_audit_path = staging / "pgx_gene_call_audit.parquet"
        _write_table(gene_audit_path, audit_rows, pgx_gene_call_audit_schema())
        outputs = {"pgx_gene_call_audit.parquet": sha256_file(gene_audit_path)}
        if core_rows:
            core_path = staging / "pharmacogenomics.parquet"
            _write_table(core_path, core_rows, pharmacogenomics_schema())
            outputs["pharmacogenomics.parquet"] = sha256_file(core_path)
        if guidance_rows:
            guidance_path = staging / "pgx_guidance_audit.parquet"
            _write_table(guidance_path, guidance_rows, pgx_guidance_audit_schema())
            outputs["pgx_guidance_audit.parquet"] = sha256_file(guidance_path)
        statuses = Counter(str(row["call_status"]) for row in audit_rows)
        summary = {
            "status": "PASS",
            "caller": "PharmCAT",
            "caller_version": PHARMCAT_VERSION,
            "data_version": data_version,
            "cpic_source_version": config.pgx.cpic_source_version,
            "cpic_snapshot_date": config.pgx.cpic_snapshot_date.isoformat(),
            "calls_tsv_sha256": calls_sha,
            "canonical_report_sha256": report_sha,
            "gene_rows": len(audit_rows),
            "published_gene_rows": len(core_rows),
            "guidance_rows": len(guidance_rows),
            "matched_guidance_rows": sum(bool(row["matched"]) for row in guidance_rows),
            "call_status_counts": dict(sorted(statuses.items())),
            "complex_gene_core_rows": sum(
                row["gene_symbol"] in STRUCTURALLY_COMPLEX_PGX_GENES for row in core_rows
            ),
            "ambiguous_gene_core_rows": 0,
            "outputs": outputs,
        }
        if summary["complex_gene_core_rows"] != 0:
            raise InputValidationError("a structurally complex VCF-only PGx call reached core")
        atomic_write_json(staging / "pgx_projection_summary.json", summary)
    return summary
