"""Bounded read-only queries for local .genome/1.1 bundles.

This module intentionally has no imports from the rest of the package. The
assembler copies it into each Agent Exploration Profile bundle as a portable
offline query script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

import duckdb

MAX_LIMIT = 500
LOCUS_PATTERN = re.compile(
    r"^(chr(?:[1-9]|1\d|2[0-2]|X|Y|M)):(\d+)(?::([ACGT]+):([ACGT]+))?$",
    re.IGNORECASE,
)


def _digest(path: Path) -> str:
    value = hashlib.sha256()
    if path.is_dir():
        for item in sorted(candidate for candidate in path.rglob("*") if candidate.is_file()):
            if item.is_symlink():
                raise ValueError(f"symbolic link is forbidden: {item}")
            value.update(item.relative_to(path).as_posix().encode("utf-8"))
            value.update(b"\0")
            with item.open("rb") as handle:
                while chunk := handle.read(8 * 1024 * 1024):
                    value.update(chunk)
    else:
        with path.open("rb") as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                value.update(chunk)
    return value.hexdigest()


def discover_bundle(start: str | Path | None = None) -> Path:
    candidate = Path(start or Path.cwd()).expanduser().resolve()
    if candidate.is_file():
        candidate = candidate.parent
    search = [candidate, *candidate.parents]
    script = Path(__file__).resolve()
    search.extend([script.parent, *script.parents])
    seen: set[Path] = set()
    for directory in search:
        if directory in seen:
            continue
        seen.add(directory)
        if (directory / "manifest.json").is_file() and (directory / "schema.json").is_file():
            return directory
    raise ValueError("could not discover a .genome bundle root")


def validate_integrity(bundle: str | Path) -> dict[str, Any]:
    root = Path(bundle).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"bundle directory not found: {root}")
    for item in root.rglob("*"):
        if item.is_symlink() or (item.exists() and not item.is_file() and not item.is_dir()):
            raise ValueError(f"forbidden bundle entry: {item.relative_to(root)}")
    manifest_path = root / "manifest.json"
    schema_path = root / "schema.json"
    variants_path = root / "variants.parquet"
    if not manifest_path.is_file() or not schema_path.is_file() or not variants_path.is_dir():
        raise ValueError("bundle is missing manifest.json, schema.json, or variants.parquet/")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.1.0" or schema.get("version") != "1.1.0":
        raise ValueError("unsupported or inconsistent schema version")
    checked = 0
    for relative, metadata in manifest.get("files", {}).items():
        relative_path = Path(relative)
        if relative_path.is_absolute() or ".." in relative_path.parts:
            raise ValueError(f"unsafe path in manifest: {relative}")
        target = root / relative_path
        if not target.exists():
            raise ValueError(f"manifest entry is missing: {relative}")
        actual = _digest(target)
        expected = metadata.get("sha256")
        if actual != expected:
            raise ValueError(f"hash mismatch for {relative}")
        checked += 1
    return {
        "valid": True,
        "schema_version": "1.1.0",
        "files_verified": checked,
        "manifest": manifest,
    }


def _sql_path(path: Path) -> str:
    return str(path).replace("'", "''")


def _dataset_glob(path: Path) -> str:
    return _sql_path(path / "**" / "*.parquet")


def _register_views(connection: duckdb.DuckDBPyConnection, root: Path) -> set[str]:
    views: set[str] = set()
    datasets = {
        "variants": root / "variants.parquet",
        "observed_partial_variants": root / "extensions" / "observed_partial_variants.parquet",
        "gene_index": root / "gene_index.parquet",
        "callability": root / "callability.parquet",
        "callable_regions": root / "callable_regions.parquet",
        "pharmacogenomics": root / "pharmacogenomics.parquet",
        "prs": root / "prs.parquet",
        "clinical_findings": root / "clinical_findings.parquet",
        "clinical_evidence": root / "clinical_evidence.parquet",
        "clinvar_assertions": root / "extensions" / "clinvar_assertions.parquet",
        "acmg_secondary_findings_evidence": (
            root / "extensions" / "acmg_secondary_findings_evidence.parquet"
        ),
        "interpretations": root / "interpretations.parquet",
        "gwas_associations": root / "gwas_associations.parquet",
        "gwas_catalog_evidence": root / "extensions" / "gwas_catalog_evidence.parquet",
        "variant_provenance": root / "extensions" / "variant_provenance.parquet",
        "reference_evidence": root / "extensions" / "reference_confidence_regions.parquet",
        "non_core_variants": root / "extensions" / "non_core_variants.parquet",
        "target_callability_evidence": (
            root / "extensions" / "target_callability_evidence.parquet"
        ),
        "pgs_model_audit": root / "extensions" / "pgs_model_audit.parquet",
        "pgs_locus_audit": root / "extensions" / "pgs_locus_audit.parquet",
        "pgx_locus_audit": root / "extensions" / "pgx_locus_audit.parquet",
        "pgx_gene_call_audit": root / "extensions" / "pgx_gene_call_audit.parquet",
        "pgx_guidance_audit": root / "extensions" / "pgx_guidance_audit.parquet",
        "family_evidence": root / "extensions" / "family_evidence.parquet",
        "clingen_gene_disease": root / "extensions" / "clingen_gene_disease.parquet",
        "clingen_dosage_sensitivity": (root / "extensions" / "clingen_dosage_sensitivity.parquet"),
        "structural_variants": (root / "extensions" / "release_c" / "structural_variants.parquet"),
        "copy_number": root / "extensions" / "release_c" / "copy_number.parquet",
        "repeat_expansions": (root / "extensions" / "release_c" / "repeat_expansions.parquet"),
        "hla_typing": root / "extensions" / "release_c" / "hla_typing.parquet",
        "complex_pharmacogenes": (
            root / "extensions" / "release_c" / "complex_pharmacogenes.parquet"
        ),
        "mitochondrial_heteroplasmy": (
            root / "extensions" / "release_c" / "mitochondrial_heteroplasmy.parquet"
        ),
        "read_backed_evidence": (
            root / "extensions" / "release_c" / "read_backed_evidence.parquet"
        ),
    }
    for name, path in datasets.items():
        if path.is_dir():
            expression = f"read_parquet('{_dataset_glob(path)}')"
        elif path.is_file():
            expression = f"read_parquet('{_sql_path(path)}')"
        else:
            continue
        connection.execute(f"CREATE VIEW {name} AS SELECT * FROM {expression}")
        views.add(name)
    return views


def _json_value(value: Any) -> Any:
    """Make DuckDB values stable and directly JSON serializable."""

    if isinstance(value, datetime):
        if value.tzinfo is not None:
            return value.astimezone(UTC).isoformat().replace("+00:00", "Z")
        return value.isoformat()
    if isinstance(value, date):
        return value.isoformat()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _bounded(
    connection: duckdb.DuckDBPyConnection,
    sql: str,
    parameters: list[Any],
    limit: int,
    table: str,
) -> tuple[list[dict[str, Any]], bool]:
    cursor = connection.execute(f"{sql} LIMIT ?", [*parameters, limit + 1])
    names = [description[0] for description in cursor.description]
    rows = [
        {name: _json_value(value) for name, value in zip(names, values, strict=True)}
        for values in cursor.fetchall()
    ]
    truncated = len(rows) > limit
    rows = rows[:limit]
    for row in rows:
        row["_table"] = table
    return rows, truncated


def _variant_projection() -> str:
    return """
        SELECT variant_id, chrom, pos, ref, alt, variant_type, rsid,
               quality.call_confidence AS call_confidence,
               quality.origin AS origin,
               quality.dp AS dp, quality.gq AS gq,
               genotype.gt AS gt, genotype.phased AS phased,
               genotype.zygosity AS zygosity,
               gene.symbol AS gene_symbol,
               consequence.most_severe AS consequence,
               consequence.hgvsc AS hgvsc, consequence.hgvsp AS hgvsp,
               population_frequency.gnomad_af_popmax AS gnomad_af_popmax,
               population_frequency.is_rare AS is_rare,
               pathogenicity.clinvar_significance AS clinvar_significance,
               pathogenicity.clinvar_has_conflicts AS clinvar_has_conflicts,
               pathogenicity.clinvar_review_stars AS clinvar_review_stars,
               pathogenicity.clinvar_id AS clinvar_id,
               pathogenicity.is_actionable AS is_actionable,
               predictions.consensus_damaging AS consensus_damaging,
               clinical_grade
        FROM variants
    """


def _partial_variants(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    predicate: str,
    parameters: list[Any],
    limit: int,
) -> tuple[list[dict[str, Any]], bool]:
    if "observed_partial_variants" not in views:
        return [], False
    return _bounded(
        connection,
        f"""
        SELECT variant_id, chrom, pos, ref, alt, rsid, gt, phased,
               gene.symbol AS gene_symbol, consequence.most_severe AS consequence,
               consequence.hgvsc AS hgvsc, consequence.hgvsp AS hgvsp,
               source_id, source_alt_index, source_gt
        FROM observed_partial_variants WHERE {predicate}
        ORDER BY chrom, pos, ref, alt, source_ordinal, source_alt_index
        """,
        parameters,
        limit,
        "extensions/observed_partial_variants",
    )


def _parse_locus(term: str) -> tuple[str, int, str | None, str | None]:
    match = LOCUS_PATTERN.fullmatch(term)
    if not match:
        raise ValueError("locus must be chrN:pos or chrN:pos:REF:ALT")
    chrom = match.group(1)
    chrom = "chr" + chrom[3:].upper().replace("MT", "M")
    return chrom, int(match.group(2)), match.group(3), match.group(4)


def _callability(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    term: str,
    limit: int,
) -> dict[str, Any]:
    if term.lower().startswith("rs"):
        tables = ["variants"]
        if "observed_partial_variants" in views:
            tables.append("observed_partial_variants")
        hits = " UNION ".join(
            f"SELECT chrom, pos, ref, alt FROM {table} WHERE lower(rsid) = lower(?)"
            for table in tables
        )
        row = connection.execute(
            f"SELECT DISTINCT * FROM ({hits}) LIMIT 2", [term] * len(tables)
        ).fetchall()
        if len(row) != 1 and "target_callability_evidence" in views:
            row = connection.execute(
                """
                SELECT chrom, pos, ref, alt
                FROM target_callability_evidence
                WHERE lower(target_id) = lower(?) AND ref IS NOT NULL AND alt IS NOT NULL
                LIMIT 2
                """,
                [term],
            ).fetchall()
        if len(row) != 1:
            return {
                "conclusion": "unknown",
                "reason": "rsID did not resolve uniquely inside the bundle",
                "rows": [],
                "count": 0,
                "truncated": False,
            }
        term = f"{row[0][0]}:{row[0][1]}:{row[0][2]}:{row[0][3]}"
    chrom, pos, ref, alt = _parse_locus(term)
    parameters: list[Any] = [chrom, pos]
    predicate = "chrom = ? AND pos = ?"
    if ref and alt:
        predicate += " AND ref = ? AND alt = ?"
        parameters.extend([ref.upper(), alt.upper()])
    variant_rows, variant_truncated = _bounded(
        connection,
        _variant_projection() + f" WHERE {predicate} ORDER BY variant_id",
        parameters,
        limit,
        "variants",
    )
    partial_rows, partial_truncated = _partial_variants(
        connection, views, predicate, parameters, limit
    )
    evidence_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    callability_rows: list[dict[str, Any]] = []
    target_evidence_rows: list[dict[str, Any]] = []
    evidence_truncated = False
    region_truncated = False
    callability_truncated = False
    target_evidence_truncated = False
    target_predicate = "chrom = ? AND pos = ?"
    target_parameters: list[Any] = [chrom, pos]
    if ref and alt:
        target_predicate += " AND variant_id = ?"
        target_parameters.append(f"{chrom}:{pos}:{ref.upper()}:{alt.upper()}")
    if "callability" in views:
        callability_rows, callability_truncated = _bounded(
            connection,
            "SELECT * FROM callability WHERE " + target_predicate + " ORDER BY callability_id",
            target_parameters,
            limit,
            "callability",
        )
    if "target_callability_evidence" in views:
        extension_predicate = "chrom = ? AND pos = ?"
        extension_parameters: list[Any] = [chrom, pos]
        if ref and alt:
            extension_predicate += " AND ref = ? AND alt = ?"
            extension_parameters.extend([ref.upper(), alt.upper()])
        target_evidence_rows, target_evidence_truncated = _bounded(
            connection,
            """
            SELECT callability_id, target_id, target_source, target_source_version,
                   catalog_sha256, catalog_row_digest, catalog_row_number,
                   chrom, pos, end_pos, ref, alt, variant_id, gene_symbol,
                   resolution, reference_observed, callable, call_confidence,
                   depth, genotype_quality, observed_ploidy, assay_scope, evidence_scope,
                   observed_variant_ids, observed_variant_source_digests,
                   observed_non_core_evidence_ids, overlapping_no_call_evidence_ids,
                   selected_reference_evidence_ids, callable_region_ids,
                   callable_min_dp, callable_min_gq, precedence_rule, reason
            FROM target_callability_evidence WHERE
            """
            + extension_predicate
            + " ORDER BY callability_id",
            extension_parameters,
            limit,
            "extensions/target_callability_evidence",
        )
    if "reference_evidence" in views:
        evidence_rows, evidence_truncated = _bounded(
            connection,
            """
            SELECT evidence_id, chrom, start_pos, end_pos, state,
                   reference_observed, callable, call_confidence,
                   ploidy, gt, dp, minimum_dp, gq, reason
            FROM reference_evidence
            WHERE chrom = ? AND start_pos <= ? AND end_pos >= ?
            ORDER BY start_pos, end_pos, evidence_id
            """,
            [chrom, pos, pos],
            limit,
            "extensions/reference_confidence_regions",
        )
    if "callable_regions" in views:
        region_rows, region_truncated = _bounded(
            connection,
            """
            SELECT region_id, chrom, start_pos, end_pos, callable,
                   minimum_depth, assay_scope, derivation
            FROM callable_regions
            WHERE chrom = ? AND start_pos <= ? AND end_pos >= ?
            ORDER BY start_pos, end_pos, region_id
            """,
            [chrom, pos, pos],
            limit,
            "callable_regions",
        )

    resolutions = {str(row["resolution"]) for row in target_evidence_rows}
    materialized_reason = (
        str(target_evidence_rows[0]["reason"]) if len(target_evidence_rows) == 1 else None
    )
    if partial_rows and not variant_rows:
        conclusion = "partial_non_reference_observed"
        reason = "the ALT is present in an incomplete genotype; its exact dosage is unknown"
    elif "non_reference_observed" in resolutions:
        conclusion = "non_reference_observed"
        reason = materialized_reason or "a materialized target has a carried target/site allele"
    elif "overlapping_variant_observed" in resolutions:
        conclusion = "overlapping_variant_observed"
        reason = materialized_reason or "a carried allele overlaps the requested target"
    elif "no_call" in resolutions:
        conclusion = "unknown_no_call"
        reason = materialized_reason or "materialized no-call evidence forbids a negative claim"
    elif "low_confidence_reference" in resolutions:
        conclusion = "low_confidence_reference_observation"
        reason = materialized_reason or "reference evidence did not pass configured thresholds"
    elif "reference_observed" in resolutions:
        conclusion = "reference_observed"
        reason = materialized_reason or "materialized complete reference evidence passed thresholds"
    elif "uncallable_without_reference" in resolutions:
        conclusion = "unknown_uncallable"
        reason = materialized_reason or "the target is interval-uncallable without a genotype"
    elif "callable_without_reference" in resolutions:
        conclusion = "callable_without_explicit_reference"
        reason = (
            materialized_reason or "interval callability does not establish a reference genotype"
        )
    elif target_evidence_rows:
        conclusion = "unknown"
        reason = (
            materialized_reason or "the materialized target has no applicable genotype evidence"
        )
    elif variant_rows:
        conclusion = "non_reference_observed"
        reason = "a matching carried biological ALT is present"
    elif any(row["state"] == "no_call" for row in evidence_rows):
        conclusion = "unknown_no_call"
        reason = "overlapping incomplete genotype evidence forbids a negative claim"
    elif any(
        row["reference_observed"] is True and row["callable"] is True for row in evidence_rows
    ):
        conclusion = "reference_observed"
        reason = "explicit reference evidence covers the locus and passes configured thresholds"
    elif any(row["reference_observed"] is True for row in evidence_rows):
        conclusion = "low_confidence_reference_observation"
        reason = "reference evidence exists but does not pass configured callability thresholds"
    elif any(row["callable"] is True for row in region_rows):
        conclusion = "callable_without_explicit_reference"
        reason = "callability alone is not an explicit reference genotype"
    else:
        conclusion = "unknown"
        reason = "no applicable variant, reference, or callability evidence was found"
    rows = [
        *variant_rows,
        *partial_rows,
        *callability_rows,
        *target_evidence_rows,
        *evidence_rows,
        *region_rows,
    ]
    return {
        "query_locus": {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt},
        "conclusion": conclusion,
        "reason": reason,
        "rows": rows,
        "count": len(rows),
        "truncated": (
            variant_truncated
            or partial_truncated
            or callability_truncated
            or target_evidence_truncated
            or evidence_truncated
            or region_truncated
        ),
    }


def _clinical(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    term: str,
    limit: int,
) -> tuple[list[dict[str, Any]], bool, dict[str, Any]]:
    wildcard = f"%{term}%"
    rows: list[dict[str, Any]] = []
    truncated = False
    searched: list[str] = []

    def append_query(sql: str, parameters: list[Any], table: str) -> None:
        nonlocal truncated
        remaining = max(0, limit - len(rows))
        result, result_truncated = _bounded(connection, sql, parameters, remaining, table)
        rows.extend(result)
        truncated = truncated or result_truncated
        searched.append(table)

    if "clinical_findings" in views:
        append_query(
            """
            SELECT finding_id, variant_id, gene_symbol, condition, claim_type,
                   classification, clinical_grade, evidence_ids
            FROM clinical_findings
            WHERE lower(condition) LIKE lower(?) OR
                  lower(coalesce(classification,'')) LIKE lower(?) OR
                  lower(coalesce(gene_symbol,'')) LIKE lower(?) OR
                  lower(coalesce(variant_id,'')) = lower(?) OR
                  lower(finding_id) = lower(?)
            ORDER BY clinical_grade DESC, finding_id
            """,
            [wildcard, wildcard, wildcard, term, term],
            "clinical_findings",
        )

    if "interpretations" in views:
        append_query(
            """
            SELECT interpretation_id, finding_id, variant_id, gene_symbol,
                   interpretation_status, claim_type, condition, inheritance_model,
                   required_genotype, observed_genotype, genotype_matches_model,
                   analytical_confidence, classification_review_status,
                   actionability_basis, confirmation_required, reviewed_by,
                   reviewed_at, evidence_snapshot, evidence_ids
            FROM interpretations
            WHERE lower(coalesce(condition,'')) LIKE lower(?) OR
                  lower(coalesce(gene_symbol,'')) LIKE lower(?) OR
                  lower(coalesce(actionability_basis,'')) LIKE lower(?) OR
                  lower(coalesce(classification_review_status,'')) LIKE lower(?)
            ORDER BY interpretation_status, interpretation_id
            """,
            [wildcard] * 4,
            "interpretations",
        )

    if "acmg_secondary_findings_evidence" in views:
        append_query(
            """
            SELECT decision_id, variant_id, gene_symbol, disease,
                   disease_identifiers, inheritance, variants_to_report,
                   decision_status, is_actionable, clinical_grade,
                   confirmation_required, analytical_confidence,
                   clinvar_significance, clinvar_review_stars,
                   clinvar_has_conflicts, clinvar_accessions,
                   matched_condition_identifiers, policy_version,
                   policy_snapshot, guidance, reasons
            FROM acmg_secondary_findings_evidence
            WHERE lower(disease) LIKE lower(?) OR
                  lower(gene_symbol) LIKE lower(?) OR
                  lower(decision_status) LIKE lower(?) OR
                  lower(array_to_string(disease_identifiers, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(clinvar_accessions, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(matched_condition_identifiers, ' | ')) LIKE lower(?)
            ORDER BY is_actionable DESC, clinical_grade DESC,
                     clinvar_review_stars DESC, variant_id
            """,
            [wildcard] * 6,
            "extensions/acmg_secondary_findings_evidence",
        )

    if "variants" in views:
        append_query(
            """
            SELECT variant_id, gene.symbol AS gene_symbol,
                   pathogenicity.clinvar_significance AS clinvar_significance,
                   pathogenicity.clinvar_has_conflicts AS clinvar_has_conflicts,
                   pathogenicity.clinvar_conflict_summary AS clinvar_conflict_summary,
                   pathogenicity.clinvar_review_stars AS clinvar_review_stars,
                   pathogenicity.clinvar_id AS clinvar_id,
                   pathogenicity.is_actionable AS is_actionable,
                   clinical_grade
            FROM variants
            WHERE lower(coalesce(pathogenicity.clinvar_significance,'')) LIKE lower(?) OR
                  lower(coalesce(pathogenicity.clinvar_conflict_summary,'')) LIKE lower(?) OR
                  lower(coalesce(pathogenicity.clinvar_id,'')) LIKE lower(?) OR
                  lower(coalesce(gene.symbol,'')) LIKE lower(?)
            ORDER BY clinical_grade DESC, pathogenicity.clinvar_review_stars DESC,
                     variant_id
            """,
            [wildcard, wildcard, wildcard, wildcard],
            "variants",
        )

    if "clinical_evidence" in views:
        append_query(
            """
            SELECT evidence_id, finding_id, variant_id, source, source_record_id,
                   source_version, assertion, review_status, retrieved_at
            FROM clinical_evidence
            WHERE lower(assertion) LIKE lower(?) OR
                  lower(coalesce(review_status,'')) LIKE lower(?) OR
                  lower(source) LIKE lower(?) OR
                  lower(coalesce(source_record_id,'')) = lower(?) OR
                  lower(coalesce(variant_id,'')) = lower(?) OR
                  lower(coalesce(finding_id,'')) = lower(?)
            ORDER BY retrieved_at DESC, source, source_record_id
            """,
            [wildcard, wildcard, wildcard, term, term, term],
            "clinical_evidence",
        )

    if "clinvar_assertions" in views:
        append_query(
            """
            SELECT evidence_id, variant_id, clinvar_accession, raw_significance,
                   aggregate_significance, review_status, review_stars,
                   has_conflicts, conflict_classifications,
                   condition_names, condition_identifiers, gene_info, rsids, hgvs, snapshot
            FROM clinvar_assertions
            WHERE lower(array_to_string(condition_names, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(condition_identifiers, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(raw_significance, ' | ')) LIKE lower(?) OR
                  lower(coalesce(aggregate_significance,'')) LIKE lower(?) OR
                  lower(array_to_string(review_status, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(conflict_classifications, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(gene_info, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(rsids, ' | ')) LIKE lower(?) OR
                  lower(array_to_string(hgvs, ' | ')) LIKE lower(?) OR
                  lower(clinvar_accession) = lower(?) OR
                  lower(clinvar_variation_id) = lower(?)
            ORDER BY review_stars DESC, variant_id, clinvar_accession
            """,
            [wildcard] * 9 + [term, term],
            "extensions/clinvar_assertions",
        )

    if "clingen_gene_disease" in views:
        append_query(
            """
            SELECT evidence_id, assertion_id, gene_symbol, hgnc_id, disease_label,
                   mondo_id, mode_of_inheritance, sop, classification, online_report,
                   classification_date, gcep, snapshot
            FROM clingen_gene_disease
            WHERE lower(gene_symbol) LIKE lower(?) OR
                  lower(disease_label) LIKE lower(?) OR
                  lower(mondo_id) LIKE lower(?) OR
                  lower(mode_of_inheritance) LIKE lower(?) OR
                  lower(classification) LIKE lower(?) OR
                  lower(gcep) LIKE lower(?) OR
                  lower(assertion_id) = lower(?)
            ORDER BY gene_symbol, disease_label, assertion_id
            """,
            [wildcard] * 6 + [term],
            "extensions/clingen_gene_disease",
        )

    if "clingen_dosage_sensitivity" in views:
        append_query(
            """
            SELECT evidence_id, gene_symbol, ncbi_gene_id, cytoband, chrom, start, "end",
                   haploinsufficiency_score, haploinsufficiency_source_score,
                   haploinsufficiency_description, haploinsufficiency_pmids,
                   haploinsufficiency_disease_id, triplosensitivity_score,
                   triplosensitivity_source_score, triplosensitivity_description,
                   triplosensitivity_pmids, triplosensitivity_disease_id,
                   last_evaluated, snapshot
            FROM clingen_dosage_sensitivity
            WHERE lower(gene_symbol) LIKE lower(?) OR
                  lower(cytoband) LIKE lower(?) OR
                  lower(haploinsufficiency_description) LIKE lower(?) OR
                  lower(coalesce(haploinsufficiency_disease_id,'')) LIKE lower(?) OR
                  lower(triplosensitivity_description) LIKE lower(?) OR
                  lower(coalesce(triplosensitivity_disease_id,'')) LIKE lower(?)
            ORDER BY gene_symbol
            """,
            [wildcard] * 6,
            "extensions/clingen_dosage_sensitivity",
        )

    extra: dict[str, Any] = {
        "tables_searched": searched,
        "interpretation": (
            "raw and aggregate external clinical evidence plus any explicitly labelled "
            "policy candidates; no diagnosis was inferred. Review status, reviewer fields, "
            "and confirmation_required describe review state; clinical_grade is only the "
            "bundle policy flag. ClinGen rows are gene-level context only and are never "
            "variant pathogenicity or subject diagnosis claims"
        ),
    }
    if not searched:
        extra["unavailable"] = "no clinical evidence table was built"
    return rows, truncated, extra


RELEASE_C_VIEWS = (
    "structural_variants",
    "copy_number",
    "repeat_expansions",
    "hla_typing",
    "complex_pharmacogenes",
    "mitochondrial_heteroplasmy",
    "read_backed_evidence",
)


def _read_backed(
    connection: duckdb.DuckDBPyConnection,
    views: set[str],
    term: str,
    limit: int,
) -> tuple[list[dict[str, Any]], bool, dict[str, Any]]:
    """Search typed Release-C observations without collapsing their source rows."""

    rows: list[dict[str, Any]] = []
    truncated = False
    searched: list[str] = []
    locus_match = LOCUS_PATTERN.fullmatch(term)

    def append_query(sql: str, parameters: list[Any], table: str) -> None:
        nonlocal truncated
        remaining = max(0, limit - len(rows))
        result, result_truncated = _bounded(
            connection,
            sql,
            parameters,
            remaining,
            "extensions/release_c/" + table,
        )
        rows.extend(result)
        truncated = truncated or result_truncated
        searched.append("extensions/release_c/" + table)

    if locus_match:
        chrom, pos, ref, alt = _parse_locus(term)
        if "structural_variants" in views:
            append_query(
                """
                SELECT * FROM structural_variants
                WHERE (chrom1 = ? AND pos1 = ?) OR (chrom2 = ? AND pos2 = ?)
                ORDER BY chrom1, pos1, chrom2, pos2, evidence_id
                """,
                [chrom, pos, chrom, pos],
                "structural_variants",
            )
        if "copy_number" in views:
            append_query(
                """
                SELECT * FROM copy_number
                WHERE chrom = ? AND start_pos <= ? AND end_pos >= ?
                ORDER BY start_pos, end_pos, evidence_id
                """,
                [chrom, pos, pos],
                "copy_number",
            )
        if "repeat_expansions" in views:
            append_query(
                """
                SELECT * FROM repeat_expansions
                WHERE chrom = ? AND start_pos <= ? AND end_pos >= ?
                ORDER BY start_pos, end_pos, evidence_id
                """,
                [chrom, pos, pos],
                "repeat_expansions",
            )
        for table in ("mitochondrial_heteroplasmy", "read_backed_evidence"):
            if table not in views:
                continue
            predicate = "chrom = ? AND pos = ?"
            parameters: list[Any] = [chrom, pos]
            if ref is not None and alt is not None:
                predicate += " AND ref = ? AND alt = ?"
                parameters.extend([ref.upper(), alt.upper()])
            append_query(
                f"SELECT * FROM {table} WHERE {predicate} ORDER BY evidence_id",
                parameters,
                table,
            )
    else:
        wildcard = f"%{term}%"
        if "structural_variants" in views:
            append_query(
                """
                SELECT * FROM structural_variants
                WHERE lower(event_id) = lower(?) OR lower(evidence_id) = lower(?) OR
                      lower(caller_record_id) = lower(?) OR
                      lower(sv_type) LIKE lower(?) OR lower(caller_name) LIKE lower(?) OR
                      lower(coalesce(concordance_group_id,'')) = lower(?)
                ORDER BY chrom1, pos1, chrom2, pos2, evidence_id
                """,
                [term, term, term, wildcard, wildcard, term],
                "structural_variants",
            )
        if "copy_number" in views:
            append_query(
                """
                SELECT * FROM copy_number
                WHERE lower(event_id) = lower(?) OR lower(evidence_id) = lower(?) OR
                      lower(call_type) LIKE lower(?) OR lower(caller_model) LIKE lower(?) OR
                      lower(caller_name) LIKE lower(?) OR
                      lower(coalesce(concordance_group_id,'')) = lower(?)
                ORDER BY chrom, start_pos, end_pos, evidence_id
                """,
                [term, term, wildcard, wildcard, wildcard, term],
                "copy_number",
            )
        if "repeat_expansions" in views:
            append_query(
                """
                SELECT * FROM repeat_expansions
                WHERE lower(event_id) = lower(?) OR lower(evidence_id) = lower(?) OR
                      lower(locus_id) LIKE lower(?) OR lower(variant_id) LIKE lower(?) OR
                      lower(motif) = lower(?) OR lower(caller_name) LIKE lower(?) OR
                      lower(coalesce(concordance_group_id,'')) = lower(?)
                ORDER BY chrom, start_pos, end_pos, evidence_id
                """,
                [term, term, wildcard, wildcard, term, wildcard, term],
                "repeat_expansions",
            )
        if "hla_typing" in views:
            append_query(
                """
                SELECT * FROM hla_typing
                WHERE lower(event_id) = lower(?) OR lower(evidence_id) = lower(?) OR
                      lower(locus) LIKE lower(?) OR lower(coalesce(allele1,'')) LIKE lower(?) OR
                      lower(coalesce(allele2,'')) LIKE lower(?) OR
                      lower(coalesce(allele1_two_field,'')) LIKE lower(?) OR
                      lower(coalesce(allele2_two_field,'')) LIKE lower(?) OR
                      lower(caller_name) LIKE lower(?) OR
                      lower(coalesce(concordance_group_id,'')) = lower(?)
                ORDER BY locus, caller_name, evidence_id
                """,
                [term, term, wildcard, wildcard, wildcard, wildcard, wildcard, wildcard, term],
                "hla_typing",
            )
        if "complex_pharmacogenes" in views:
            append_query(
                """
                SELECT * FROM complex_pharmacogenes
                WHERE lower(event_id) = lower(?) OR lower(evidence_id) = lower(?) OR
                      lower(gene) LIKE lower(?) OR lower(coalesce(diplotype,'')) LIKE lower(?) OR
                      lower(coalesce(phenotype,'')) LIKE lower(?) OR
                      lower(caller_name) LIKE lower(?) OR
                      lower(coalesce(concordance_group_id,'')) = lower(?)
                ORDER BY gene, caller_name, evidence_id
                """,
                [term, term, wildcard, wildcard, wildcard, wildcard, term],
                "complex_pharmacogenes",
            )
        if "mitochondrial_heteroplasmy" in views:
            append_query(
                """
                SELECT * FROM mitochondrial_heteroplasmy
                WHERE lower(event_id) = lower(?) OR lower(evidence_id) = lower(?) OR
                      lower(caller_name) LIKE lower(?) OR
                      lower(coalesce(concordance_group_id,'')) = lower(?)
                ORDER BY chrom, pos, ref, alt, evidence_id
                """,
                [term, term, wildcard, term],
                "mitochondrial_heteroplasmy",
            )
        if "read_backed_evidence" in views:
            append_query(
                """
                SELECT * FROM read_backed_evidence
                WHERE lower(event_id) = lower(?) OR lower(evidence_id) = lower(?) OR
                      lower(target_id) LIKE lower(?) OR
                      lower(candidate_category) LIKE lower(?) OR
                      lower(caller_name) LIKE lower(?) OR
                      lower(coalesce(concordance_group_id,'')) = lower(?)
                ORDER BY chrom, pos, ref, alt, evidence_id
                """,
                [term, term, wildcard, wildcard, wildcard, term],
                "read_backed_evidence",
            )

    extra: dict[str, Any] = {
        "tables_searched": searched,
        "interpretation": (
            "Each row is a source-specific read-backed observation. Concordance annotations "
            "do not replace, average, or vote away caller or technical-replicate rows; "
            "confirmation_required and limitations must be respected. These rows are not "
            "clinical diagnoses."
        ),
    }
    if not searched:
        extra["unavailable"] = "read-backed Release-C extensions were not built"
    return rows, truncated, extra


def query_bundle(
    bundle: str | Path,
    operation: str,
    *,
    term: str | None = None,
    limit: int = 50,
) -> dict[str, Any]:
    if not 1 <= limit <= MAX_LIMIT:
        raise ValueError(f"limit must be between 1 and {MAX_LIMIT}")
    validation = validate_integrity(bundle)
    root = Path(bundle).expanduser().resolve()
    manifest = validation.pop("manifest")
    connection = duckdb.connect(":memory:", read_only=False)
    try:
        views = _register_views(connection, root)
        rows: list[dict[str, Any]] = []
        truncated = False
        extra: dict[str, Any] = {}
        if operation == "overview":
            rows, truncated = _bounded(
                connection,
                """
                SELECT count(*) AS variant_count,
                       count(*) FILTER (WHERE clinical_grade) AS clinical_grade_count,
                       count(*) FILTER (WHERE pathogenicity.is_actionable) AS actionable_count,
                       count(*) FILTER (WHERE quality.call_confidence = 'high') AS high_confidence_count,
                       count(*) FILTER (WHERE quality.call_confidence = 'medium') AS medium_confidence_count,
                       count(*) FILTER (WHERE quality.call_confidence = 'low') AS low_confidence_count
                FROM variants
                """,
                [],
                1,
                "variants",
            )
            extra["available_tables"] = sorted(views)
            if "observed_partial_variants" in views:
                extra["observed_partial_allele_rows"] = connection.execute(
                    "SELECT count(*) FROM observed_partial_variants"
                ).fetchone()[0]
            extra["read_backed_extension_rows"] = {
                table: connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
                for table in RELEASE_C_VIEWS
                if table in views
            }
        elif operation in {"variant", "gene", "search"}:
            if not term:
                raise ValueError(f"{operation} requires a search term")
            if operation == "variant":
                predicate = "lower(variant_id) = lower(?) OR lower(coalesce(rsid,'')) = lower(?)"
                parameters = [term, term]
            elif operation == "gene":
                predicate = (
                    "lower(coalesce(gene.symbol,'')) = lower(?) OR "
                    "list_contains(list_transform(consequence.all_annotations, "
                    "annotation -> lower(coalesce(annotation.symbol,''))), lower(?))"
                )
                parameters = [term, term]
            else:
                predicate = """
                    lower(variant_id) LIKE lower(?) OR lower(coalesce(rsid,'')) LIKE lower(?) OR
                    lower(coalesce(gene.symbol,'')) LIKE lower(?) OR
                    lower(coalesce(consequence.most_severe,'')) LIKE lower(?)
                """
                wildcard = f"%{term}%"
                parameters = [wildcard] * 4
            core_predicate = predicate
            core_parameters = parameters
            if operation == "search":
                core_predicate += " OR lower(coalesce(pathogenicity.clinvar_id,'')) LIKE lower(?)"
                core_parameters = [*parameters, wildcard]
            rows, truncated = _bounded(
                connection,
                _variant_projection() + f" WHERE {core_predicate} ORDER BY chrom, pos, alt",
                core_parameters,
                limit,
                "variants",
            )
            partial, partial_truncated = _partial_variants(
                connection, views, predicate, parameters, limit - len(rows)
            )
            rows.extend(partial)
            truncated = truncated or partial_truncated
            if operation == "gene":
                context_rows: list[dict[str, Any]] = []
                context_truncated = False
                for table in ("clingen_gene_disease", "clingen_dosage_sensitivity"):
                    if table not in views:
                        continue
                    available = max(0, limit - len(context_rows))
                    result, result_truncated = _bounded(
                        connection,
                        f"SELECT * FROM {table} WHERE lower(gene_symbol) = lower(?) "
                        "ORDER BY evidence_id",
                        [term],
                        available,
                        "extensions/" + table,
                    )
                    context_rows.extend(result)
                    context_truncated = context_truncated or result_truncated
                extra["clingen_gene_context"] = context_rows
                extra["clingen_gene_context_truncated"] = context_truncated
                extra["clingen_semantic_boundary"] = (
                    "Gene-level validity and dosage context only; these rows do not classify "
                    "any subject variant and do not diagnose the associated condition."
                )
        elif operation == "clinical":
            if not term:
                raise ValueError("clinical requires a condition or classification term")
            rows, truncated, clinical_extra = _clinical(connection, views, term, limit)
            extra.update(clinical_extra)
        elif operation == "pgx":
            if not term:
                raise ValueError("pgx requires a drug or gene term")
            audit_rows: list[dict[str, Any]] = []
            audit_truncated = False
            guidance_rows: list[dict[str, Any]] = []
            guidance_truncated = False
            if "pharmacogenomics" in views:
                rows, truncated = _bounded(
                    connection,
                    """
                    SELECT * FROM pharmacogenomics
                    WHERE lower(gene_symbol) = lower(?) OR
                          list_contains(list_transform(affected_drugs, value -> lower(value)), lower(?))
                    ORDER BY gene_symbol
                    """,
                    [term, term],
                    limit,
                    "pharmacogenomics",
                )
            if "pgx_gene_call_audit" in views:
                audit_rows, audit_truncated = _bounded(
                    connection,
                    """
                    SELECT * FROM pgx_gene_call_audit
                    WHERE lower(gene_symbol) = lower(?) OR
                      list_contains(
                        list_transform(affected_drugs, value -> lower(value)), lower(?)
                      )
                    ORDER BY gene_symbol
                    """,
                    [term, term],
                    limit,
                    "extensions/pgx_gene_call_audit",
                )
            if "pgx_guidance_audit" in views:
                guidance_rows, guidance_truncated = _bounded(
                    connection,
                    """
                    SELECT * FROM pgx_guidance_audit
                    WHERE lower(gene_symbol) = lower(?) OR lower(drug) = lower(?)
                    ORDER BY matched DESC, gene_symbol, drug, guidance_id
                    """,
                    [term, term],
                    limit,
                    "extensions/pgx_guidance_audit",
                )
            if rows:
                if audit_rows:
                    extra["gene_call_audit"] = audit_rows
                    extra["gene_call_audit_truncated"] = audit_truncated
                if guidance_rows:
                    extra["guidance_audit"] = guidance_rows
                    extra["guidance_audit_truncated"] = guidance_truncated
            elif audit_rows:
                rows, truncated = audit_rows, audit_truncated
                if guidance_rows:
                    extra["guidance_audit"] = guidance_rows
                    extra["guidance_audit_truncated"] = guidance_truncated
                extra["unavailable"] = (
                    "no normative PGx row passed the complete, unambiguous, non-complex, "
                    "CPIC-level publication gate; audit rows explain the recorded limitation"
                )
            elif guidance_rows:
                rows, truncated = guidance_rows, guidance_truncated
                extra["unavailable"] = (
                    "guidance exists, but no matching normative gene call was published"
                )
            elif "pharmacogenomics" not in views:
                extra["unavailable"] = "pharmacogenomics.parquet was not built"
        elif operation == "prs":
            if not term:
                raise ValueError("prs requires a trait term")
            if "prs" in views:
                rows, truncated = _bounded(
                    connection,
                    "SELECT * FROM prs WHERE lower(trait) LIKE lower(?) ORDER BY trait",
                    [f"%{term}%"],
                    limit,
                    "prs",
                )
            if "pgs_model_audit" in views:
                audit_rows, audit_truncated = _bounded(
                    connection,
                    "SELECT * FROM pgs_model_audit WHERE lower(trait) LIKE lower(?) OR "
                    "lower(pgs_id) = lower(?) ORDER BY trait, pgs_id",
                    [f"%{term}%", term],
                    limit,
                    "extensions/pgs_model_audit",
                )
                if "prs" in views:
                    extra["model_audit"] = audit_rows
                    extra["model_audit_truncated"] = audit_truncated
                else:
                    rows, truncated = audit_rows, audit_truncated
                    extra["unavailable"] = (
                        "no matching calibrated score was published; model audit rows explain "
                        "completeness or calibration limits"
                    )
            elif "prs" not in views:
                extra["unavailable"] = "prs.parquet was not built"
        elif operation == "gwas":
            if not term:
                raise ValueError("gwas requires a trait, rsID, study, PMID, or variant term")
            wildcard = f"%{term}%"
            if "gwas_catalog_evidence" in views:
                rows, truncated = _bounded(
                    connection,
                    """
                    SELECT * FROM gwas_catalog_evidence
                    WHERE lower(variant_id) = lower(?) OR
                          lower(risk_rsid) = lower(?) OR
                          list_contains(
                            list_transform(rsid_aliases, value -> lower(value)), lower(?)
                          ) OR
                          lower(reported_trait) LIKE lower(?) OR
                          lower(array_to_string(mapped_traits, ' | ')) LIKE lower(?) OR
                          lower(study_accession) = lower(?) OR
                          list_contains(study_pmids, ?)
                    ORDER BY p_value, association_id
                    """,
                    [term, term, term, wildcard, wildcard, term, term],
                    limit,
                    "extensions/gwas_catalog_evidence",
                )
            if "gwas_associations" in views:
                typed_sql = """
                    SELECT * FROM gwas_associations a
                    WHERE lower(a.variant_id) = lower(?) OR
                          lower(a.trait) LIKE lower(?) OR
                          list_contains(a.study_pmids, ?)
                """
                typed_parameters: list[object] = [term, wildcard, term]
                if "gwas_catalog_evidence" in views:
                    typed_sql += """
                          OR EXISTS (
                            SELECT 1 FROM gwas_catalog_evidence e
                            WHERE e.association_id = a.association_id AND (
                              lower(e.risk_rsid) = lower(?) OR
                              list_contains(
                                list_transform(e.rsid_aliases, value -> lower(value)), lower(?)
                              ) OR
                              lower(e.reported_trait) LIKE lower(?) OR
                              lower(array_to_string(e.mapped_traits, ' | ')) LIKE lower(?) OR
                              lower(e.study_accession) = lower(?) OR
                              list_contains(e.study_pmids, ?)
                            )
                          )
                    """
                    typed_parameters.extend([term, term, wildcard, wildcard, term, term])
                typed_sql += " ORDER BY a.p_value, a.association_id"
                typed, typed_truncated = _bounded(
                    connection,
                    typed_sql,
                    typed_parameters,
                    limit,
                    "gwas_associations",
                )
                extra["typed_associations"] = typed
                extra["typed_associations_truncated"] = typed_truncated
            if "gwas_catalog_evidence" not in views:
                extra["unavailable"] = "GWAS Catalog evidence was not built"
            extra["interpretation"] = (
                "population association evidence only; no row is a personal diagnosis, "
                "causal claim, or clinical-grade recommendation"
            )
        elif operation == "callability":
            if not term:
                raise ValueError("callability requires chrN:pos[:REF:ALT] or rsID")
            extra.update(_callability(connection, views, term, limit))
            rows = extra.pop("rows")
            truncated = extra.pop("truncated")
        elif operation == "read-backed":
            if not term:
                raise ValueError(
                    "read-backed requires a locus, event, gene, allele, or caller term"
                )
            rows, truncated, read_backed_extra = _read_backed(
                connection,
                views,
                term,
                limit,
            )
            extra.update(read_backed_extra)
        elif operation == "evidence":
            if not term:
                raise ValueError("evidence requires a variant_id")
            evidence_sources: list[str] = []
            if "variant_provenance" in views:
                rows, truncated = _bounded(
                    connection,
                    "SELECT * FROM variant_provenance WHERE lower(variant_id) = lower(?) "
                    "ORDER BY source_ordinal, source_alt_index",
                    [term],
                    limit,
                    "extensions/variant_provenance",
                )
                evidence_sources.append("extensions/variant_provenance")
            else:
                extra["unavailable"] = "variant provenance extension was not built"
            evidence_queries = [
                (
                    "observed_partial_variants",
                    "SELECT * FROM observed_partial_variants WHERE lower(variant_id) = lower(?) "
                    "ORDER BY source_ordinal, source_alt_index",
                    "extensions/observed_partial_variants",
                ),
                (
                    "clinical_evidence",
                    "SELECT * FROM clinical_evidence WHERE lower(variant_id) = lower(?) "
                    "ORDER BY source, source_record_id",
                    "clinical_evidence",
                ),
                (
                    "clinvar_assertions",
                    "SELECT * FROM clinvar_assertions WHERE lower(variant_id) = lower(?) "
                    "ORDER BY clinvar_accession",
                    "extensions/clinvar_assertions",
                ),
                (
                    "acmg_secondary_findings_evidence",
                    "SELECT * FROM acmg_secondary_findings_evidence "
                    "WHERE lower(variant_id) = lower(?) ORDER BY policy_key",
                    "extensions/acmg_secondary_findings_evidence",
                ),
                (
                    "clinical_findings",
                    "SELECT * FROM clinical_findings WHERE lower(variant_id) = lower(?) "
                    "ORDER BY finding_id",
                    "clinical_findings",
                ),
                (
                    "family_evidence",
                    "SELECT * FROM family_evidence WHERE lower(variant_id) = lower(?) "
                    "ORDER BY evidence_id",
                    "extensions/family_evidence",
                ),
                (
                    "gwas_catalog_evidence",
                    "SELECT * FROM gwas_catalog_evidence "
                    "WHERE lower(variant_id) = lower(?) ORDER BY p_value, association_id",
                    "extensions/gwas_catalog_evidence",
                ),
                (
                    "interpretations",
                    "SELECT * FROM interpretations WHERE lower(variant_id) = lower(?) "
                    "ORDER BY interpretation_id",
                    "interpretations",
                ),
            ]
            for view, sql, label in evidence_queries:
                if view not in views or len(rows) >= limit:
                    continue
                additional, additional_truncated = _bounded(
                    connection,
                    sql,
                    [term],
                    limit - len(rows),
                    label,
                )
                rows.extend(additional)
                truncated = truncated or additional_truncated
                evidence_sources.append(label)
            if any(view in views for view in RELEASE_C_VIEWS):
                read_backed_rows, read_backed_truncated, read_backed_extra = _read_backed(
                    connection,
                    views,
                    term,
                    max(0, limit - len(rows)),
                )
                rows.extend(read_backed_rows)
                truncated = truncated or read_backed_truncated
                evidence_sources.extend(read_backed_extra["tables_searched"])
                extra["read_backed_interpretation"] = read_backed_extra["interpretation"]
            extra["tables_searched"] = evidence_sources
        else:
            raise ValueError(f"unknown operation: {operation}")
    finally:
        connection.close()

    return {
        "operation": operation,
        "term": term,
        "count": len(rows),
        "truncated": truncated,
        "limit": limit,
        "rows": rows,
        "snapshot": manifest.get("generated_at"),
        "annotation_sources": manifest.get("annotation_sources", {}),
        "consensus_definitions": manifest.get("consensus_definitions", {}),
        "integrity": validation,
        **extra,
    }


def _table_output(result: dict[str, Any]) -> str:
    rows = result.get("rows") or []
    if not rows:
        return json.dumps(result, indent=2, ensure_ascii=False, default=str)
    columns = list(rows[0])
    widths = {
        column: min(
            40,
            max(len(column), *(len(str(row.get(column, ""))) for row in rows)),
        )
        for column in columns
    }
    header = " | ".join(column.ljust(widths[column]) for column in columns)
    rule = "-+-".join("-" * widths[column] for column in columns)
    body = [
        " | ".join(
            str(row.get(column, ""))[: widths[column]].ljust(widths[column]) for column in columns
        )
        for row in rows
    ]
    return "\n".join([header, rule, *body])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read-only .genome/1.1 query interface")
    parser.add_argument(
        "operation",
        choices=[
            "overview",
            "search",
            "variant",
            "gene",
            "clinical",
            "pgx",
            "prs",
            "gwas",
            "callability",
            "read-backed",
            "evidence",
        ],
    )
    parser.add_argument("term", nargs="?")
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--limit", type=int, default=50)
    parser.add_argument("--format", choices=["json", "table"], default="json")
    arguments = parser.parse_args(argv)
    try:
        bundle = discover_bundle(arguments.bundle)
        result = query_bundle(
            bundle,
            arguments.operation,
            term=arguments.term,
            limit=arguments.limit,
        )
    except (ValueError, OSError, duckdb.Error) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if arguments.format == "table":
        print(_table_output(result))
    else:
        print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
