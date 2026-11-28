"""Materialize target-level callability without inventing reference genotypes."""

from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
import subprocess
import tempfile
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .arrow_schemas import (
    callability_schema,
    target_callability_evidence_schema,
    target_catalog_schema,
)
from .config import BuildConfig
from .errors import InputValidationError
from .parquet_io import PartitionedParquetWriter
from .query import par_query_limit
from .util import (
    atomic_output_directory,
    atomic_write_json,
    canonical_json_bytes,
    sha256_file,
)
from .vcf_ingest import canonical_core_chrom

CATALOG_COLUMNS = (
    "target_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "gene_symbol",
    "target_source",
    "target_source_version",
    "evidence_scope",
)
EVIDENCE_SCOPES = {"site", "target_panel", "derived_interval"}
PRECEDENCE_RULE = (
    "normalized_partial_no_call>carried_variant>overlapping_no_call>complete_reference_evidence>"
    "covering_callable_interval>unknown"
)


def _chrom_rank(chrom: str) -> int:
    token = chrom.removeprefix("chr")
    if token.isdigit() and 1 <= int(token) <= 22:
        return int(token)
    return {"X": 23, "Y": 24, "M": 25}[token]


def _reference_contigs(fasta: pysam.FastaFile) -> dict[str, tuple[str, int]]:
    mapping: dict[str, tuple[str, int]] = {}
    for name, length in zip(fasta.references, fasta.lengths, strict=True):
        canonical = canonical_core_chrom(name)
        if canonical is None:
            continue
        if canonical in mapping:
            raise InputValidationError(
                f"reference contains ambiguous core-contig aliases for {canonical}"
            )
        mapping[canonical] = (name, int(length))
    return mapping


def _record_digest(row: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_json_bytes(row)).hexdigest()


def _write_catalog_table(
    config: BuildConfig,
    catalog_sha256: str,
    destination: Path,
    raw_vcf: Path,
) -> tuple[int, int]:
    try:
        fasta = pysam.FastaFile(str(config.reference_fasta))
    except (OSError, ValueError) as exc:
        raise InputValidationError(
            f"cannot open reference FASTA for target catalogue: {exc}"
        ) from exc
    contigs = _reference_contigs(fasta)
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value="genome-annotation-builder-target-validation")
    for canonical in sorted(contigs, key=_chrom_rank):
        name, length = contigs[canonical]
        header.contigs.add(name, length=length)
    allele_writer = pysam.VariantFile(str(raw_vcf), "w", header=header)
    rows = 0
    allele_rows = 0
    try:
        with (
            config.target_catalog_path.open("r", encoding="utf-8", newline="") as handle,
            PartitionedParquetWriter(
                destination, target_catalog_schema(), partition_field=None
            ) as writer,
        ):
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != CATALOG_COLUMNS:
                raise InputValidationError(
                    "target catalogue header must be exactly: " + "\t".join(CATALOG_COLUMNS)
                )
            for row_number, raw in enumerate(reader, start=2):
                if None in raw or any(value is None for value in raw.values()):
                    raise InputValidationError(
                        f"target catalogue row {row_number} has the wrong number of columns"
                    )
                values = {key: str(raw[key]).strip() for key in CATALOG_COLUMNS}
                for key in (
                    "target_id",
                    "chrom",
                    "pos",
                    "target_source",
                    "target_source_version",
                    "evidence_scope",
                ):
                    if not values[key]:
                        raise InputValidationError(
                            f"target catalogue row {row_number} has empty {key}"
                        )
                if values["evidence_scope"] not in EVIDENCE_SCOPES:
                    raise InputValidationError(
                        f"target catalogue row {row_number} has unsupported evidence_scope"
                    )
                try:
                    pos = int(values["pos"])
                except ValueError as exc:
                    raise InputValidationError(
                        f"target catalogue row {row_number} has a non-integer position"
                    ) from exc
                canonical_chrom = canonical_core_chrom(values["chrom"])
                if canonical_chrom is None or canonical_chrom not in contigs:
                    raise InputValidationError(
                        f"target catalogue row {row_number} is not on a GRCh38 primary contig"
                    )
                reference_name, reference_length = contigs[canonical_chrom]
                if not 1 <= pos <= reference_length:
                    raise InputValidationError(
                        f"target catalogue row {row_number} position is outside the reference"
                    )
                ref = values["ref"].upper() or None
                alt = values["alt"].upper() or None
                if (ref is None) != (alt is None):
                    raise InputValidationError(
                        f"target catalogue row {row_number} must provide both ref and alt or neither"
                    )
                if ref is not None:
                    if not ref.isalpha() or not alt or not alt.isalpha():
                        raise InputValidationError(
                            f"target catalogue row {row_number} has an invalid allele"
                        )
                    if set(ref) - set("ACGT") or set(alt) - set("ACGT"):
                        raise InputValidationError(
                            f"target catalogue row {row_number} alleles must contain only A/C/G/T"
                        )
                    if ref == alt:
                        raise InputValidationError(
                            f"target catalogue row {row_number} REF and ALT are identical"
                        )
                    end_pos = pos + len(ref) - 1
                    if end_pos > reference_length:
                        raise InputValidationError(
                            f"target catalogue row {row_number} REF extends beyond the reference"
                        )
                    observed = fasta.fetch(reference_name, pos - 1, end_pos).upper()
                    if observed != ref:
                        raise InputValidationError(
                            f"target catalogue row {row_number} REF mismatch: {ref} != {observed}"
                        )
                    variant_id = f"{canonical_chrom}:{pos}:{ref}:{alt}"
                else:
                    end_pos = pos
                    variant_id = None

                identity = {
                    "target_id": values["target_id"],
                    "chrom": canonical_chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene_symbol": values["gene_symbol"] or None,
                    "target_source": values["target_source"],
                    "target_source_version": values["target_source_version"],
                    "evidence_scope": values["evidence_scope"],
                }
                callability_id = (
                    "callability:"
                    + hashlib.sha256(
                        canonical_json_bytes(
                            {
                                "catalog_sha256": catalog_sha256,
                                "target": identity,
                            }
                        )
                    ).hexdigest()
                )
                writer.write(
                    {
                        "callability_id": callability_id,
                        **identity,
                        "catalog_sha256": catalog_sha256,
                        "catalog_row_digest": _record_digest(identity),
                        "catalog_row_number": row_number,
                        "chrom_rank": _chrom_rank(canonical_chrom),
                        "end_pos": end_pos,
                        "variant_id": variant_id,
                        "assay_scope": config.target_catalog_assay_scope,
                    }
                )
                rows += 1
                if ref is not None:
                    record = allele_writer.new_record(
                        contig=reference_name,
                        start=pos - 1,
                        id=f"GABT{row_number}",
                        alleles=(ref, alt),
                    )
                    allele_writer.write(record)
                    allele_rows += 1
    except UnicodeDecodeError as exc:
        raise InputValidationError("target catalogue must be UTF-8 TSV") from exc
    finally:
        allele_writer.close()
        fasta.close()
    if rows == 0:
        raise InputValidationError("enabled target catalogue is empty")
    return rows, allele_rows


def _validate_unique_targets(targets: Path) -> None:
    connection = duckdb.connect(":memory:")
    try:
        path = str(targets / "**" / "*.parquet").replace("'", "''")
        duplicate = connection.execute(
            f"SELECT target_id, count(*) AS n FROM read_parquet('{path}') "
            "GROUP BY target_id HAVING n > 1 ORDER BY target_id LIMIT 1"
        ).fetchone()
        duplicate_callability = connection.execute(
            f"SELECT callability_id, count(*) AS n FROM read_parquet('{path}') "
            "GROUP BY callability_id HAVING n > 1 LIMIT 1"
        ).fetchone()
    finally:
        connection.close()
    if duplicate:
        raise InputValidationError(f"target catalogue target_id is not unique: {duplicate[0]}")
    if duplicate_callability:
        raise InputValidationError("target catalogue produces duplicate callability identifiers")


def _validate_normalized_alleles(
    config: BuildConfig,
    raw_vcf: Path,
    normalized_vcf: Path,
    allele_rows: int,
) -> None:
    if allele_rows == 0:
        return
    completed = subprocess.run(
        [
            config.bcftools,
            "norm",
            "--check-ref",
            "e",
            "--fasta-ref",
            str(config.reference_fasta),
            "--output-type",
            "v",
            "--output",
            str(normalized_vcf),
            str(raw_vcf),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise InputValidationError(f"target catalogue allele normalization check failed: {detail}")
    database = sqlite3.connect(raw_vcf.with_suffix(".normalization.sqlite"))
    try:
        database.execute("PRAGMA journal_mode=OFF")
        database.execute("PRAGMA synchronous=OFF")
        database.execute(
            "CREATE TABLE expected "
            "(id TEXT PRIMARY KEY, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT)"
        )
        database.execute(
            "CREATE TABLE observed "
            "(id TEXT PRIMARY KEY, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT)"
        )

        def load(path: Path, table: str) -> int:
            source = pysam.VariantFile(str(path))
            count = 0
            batch: list[tuple[str, str | None, int, str, str]] = []
            try:
                for record in source:
                    batch.append(
                        (
                            str(record.id),
                            canonical_core_chrom(record.contig),
                            record.pos,
                            record.ref.upper(),
                            ",".join(str(value).upper() for value in (record.alts or ())),
                        )
                    )
                    count += 1
                    if len(batch) >= 10_000:
                        database.executemany(f"INSERT INTO {table} VALUES (?,?,?,?,?)", batch)
                        batch.clear()
                if batch:
                    database.executemany(f"INSERT INTO {table} VALUES (?,?,?,?,?)", batch)
            except sqlite3.IntegrityError as exc:
                raise InputValidationError(
                    "target catalogue normalization produced duplicate validation IDs"
                ) from exc
            finally:
                source.close()
            return count

        expected_count = load(raw_vcf, "expected")
        observed_count = load(normalized_vcf, "observed")
        database.commit()
        mismatch = database.execute(
            """
            SELECT e.id, e.chrom, e.pos, e.ref, e.alt,
                   o.chrom, o.pos, o.ref, o.alt
            FROM expected e LEFT JOIN observed o USING (id)
            WHERE o.id IS NULL OR e.chrom != o.chrom OR e.pos != o.pos OR
                  e.ref != o.ref OR e.alt != o.alt
            UNION ALL
            SELECT o.id, NULL, NULL, NULL, NULL, o.chrom, o.pos, o.ref, o.alt
            FROM observed o LEFT JOIN expected e USING (id)
            WHERE e.id IS NULL
            LIMIT 1
            """
        ).fetchone()
        if mismatch:
            row_number = str(mismatch[0]).removeprefix("GABT")
            before_value = tuple(mismatch[1:5])
            after_value = tuple(mismatch[5:9])
            raise InputValidationError(
                "target catalogue contains a non-normalized allele at row "
                f"{row_number}: {before_value} normalizes to {after_value}"
            )
    finally:
        database.close()
    if expected_count != allele_rows or observed_count != allele_rows:
        raise InputValidationError("target catalogue normalization accounting mismatch")


def _register_parquet(connection: duckdb.DuckDBPyConnection, name: str, path: Path) -> None:
    source = path / "**" / "*.parquet" if path.is_dir() else path
    escaped = str(source).replace("'", "''")
    connection.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{escaped}')")


def _build_hit_tables(
    connection: duckdb.DuckDBPyConnection,
    *,
    targets: Path,
    variants: Path,
    provenance: Path,
    partial_variants: Path,
    references: Path,
    regions: Path,
    non_core: Path,
) -> None:
    for name, path in {
        "targets": targets,
        "variants": variants,
        "provenance": provenance,
        "reference_evidence": references,
        "callable_regions": regions,
        "non_core": non_core,
    }.items():
        _register_parquet(connection, name, path)

    connection.execute(
        """
        CREATE TEMP TABLE variant_hits AS
        SELECT t.callability_id,
               count(*)::BIGINT AS variant_hit_count,
               bool_or(t.variant_id IS NOT NULL AND t.variant_id = v.variant_id)
                 AS exact_target_observed,
               bool_and(v.quality.call_confidence IN ('high','medium'))
                 AS variants_callable,
               CASE WHEN bool_or(v.quality.call_confidence = 'low') THEN 'low'
                    WHEN bool_or(v.quality.call_confidence = 'medium') THEN 'medium'
                    ELSE 'high' END AS variant_confidence,
               CASE WHEN count(*) = count(p.source_dp) THEN min(p.source_dp) END AS variant_depth,
               CASE WHEN count(*) = count(p.source_gq) THEN min(p.source_gq) END AS variant_gq,
               CASE WHEN count(DISTINCT list_count(v.genotype.gt)) = 1
                    THEN min(list_count(v.genotype.gt))::TINYINT END AS variant_ploidy,
               list(v.variant_id ORDER BY v.variant_id) AS observed_variant_ids,
               list(p.source_record_digest ORDER BY v.variant_id)
                 AS observed_variant_source_digests
        FROM targets t
        JOIN variants v
          ON t.chrom = v.chrom
         AND v.pos <= t.end_pos
         AND v.pos + length(v.ref) - 1 >= t.pos
        JOIN provenance p ON p.variant_id = v.variant_id
        GROUP BY t.callability_id
        """
    )
    connection.execute(
        """
        CREATE TEMP TABLE non_core_hits AS
        SELECT t.callability_id,
               count(*)::BIGINT AS non_core_hit_count,
               CASE WHEN count(*) = count(n.dp) THEN min(n.dp) END AS non_core_depth,
               CASE WHEN count(*) = count(n.gq) THEN min(n.gq) END AS non_core_gq,
               bool_or(list_count(n.source_gt) < len(n.source_gt)) AS partial_non_core_hit,
               CASE WHEN count(DISTINCT len(n.source_gt)) = 1
                    THEN min(len(n.source_gt))::TINYINT END AS non_core_ploidy,
               list(n.evidence_id ORDER BY n.evidence_id) AS non_core_evidence_ids
        FROM targets t
        JOIN non_core n
          ON t.chrom = n.chrom
         AND n.pos <= t.end_pos
         AND n.pos + length(n.ref) - 1 >= t.pos
        GROUP BY t.callability_id
        """
    )
    # A left-aligned partial allele may fall outside its source no-call interval.
    # Match both coordinates and keep each source record once per target.
    if partial_variants.exists():
        _register_parquet(connection, "partial_variants", partial_variants)
        connection.execute("""
            CREATE TEMP TABLE normalized_no_calls AS
            SELECT DISTINCT t.callability_id, p.source_id AS evidence_id
            FROM targets t JOIN partial_variants p
              ON p.chrom = t.chrom AND p.pos <= t.end_pos AND p.pos + length(p.ref) - 1 >= t.pos
        """)
    else:
        connection.execute("""
            CREATE TEMP TABLE normalized_no_calls(callability_id VARCHAR, evidence_id VARCHAR)
        """)
    connection.execute(
        """
        CREATE TEMP TABLE reference_hits AS
        WITH hits AS (
          SELECT t.callability_id, r.evidence_id
          FROM targets t JOIN reference_evidence r
            ON t.chrom = r.chrom AND r.start_pos <= t.end_pos AND r.end_pos >= t.pos
          UNION
          SELECT callability_id, evidence_id FROM normalized_no_calls
        )
        SELECT t.callability_id,
               count(*) FILTER (WHERE r.state = 'no_call')::BIGINT AS no_call_count,
               list(r.evidence_id ORDER BY r.evidence_id)
                 FILTER (WHERE r.state = 'no_call') AS no_call_evidence_ids,
               CASE WHEN count(*) FILTER (WHERE r.state = 'no_call') =
                              count(r.dp) FILTER (WHERE r.state = 'no_call')
                    THEN min(r.dp) FILTER (WHERE r.state = 'no_call') END AS no_call_depth,
               CASE WHEN count(*) FILTER (WHERE r.state = 'no_call') =
                              count(r.gq) FILTER (WHERE r.state = 'no_call')
                    THEN min(r.gq) FILTER (WHERE r.state = 'no_call') END AS no_call_gq,
               CASE WHEN count(DISTINCT r.ploidy) FILTER (WHERE r.state = 'no_call') = 1
                    THEN min(r.ploidy) FILTER (WHERE r.state = 'no_call') END AS no_call_ploidy,
               count(*) FILTER (
                 WHERE r.state IN ('reference_block','explicit_reference_site')
                   AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
               )::BIGINT AS reference_count,
               list(r.evidence_id ORDER BY r.evidence_id) FILTER (
                 WHERE r.state IN ('reference_block','explicit_reference_site')
                   AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
               ) AS reference_evidence_ids,
               bool_and(r.callable IS TRUE) FILTER (
                 WHERE r.state IN ('reference_block','explicit_reference_site')
                   AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
               ) AS references_callable,
               CASE
                 WHEN bool_or(r.call_confidence = 'low') FILTER (
                   WHERE r.state IN ('reference_block','explicit_reference_site')
                     AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                 ) THEN 'low'
                 WHEN bool_or(r.call_confidence = 'medium') FILTER (
                   WHERE r.state IN ('reference_block','explicit_reference_site')
                     AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                 ) THEN 'medium'
                 WHEN count(*) FILTER (
                   WHERE r.state IN ('reference_block','explicit_reference_site')
                     AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                 ) > 0 THEN 'high'
               END AS reference_confidence,
               CASE WHEN count(*) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            ) = count(coalesce(r.minimum_dp,r.dp)) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            )
                    THEN min(coalesce(r.minimum_dp,r.dp)) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            ) END AS reference_depth,
               CASE WHEN count(*) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            ) = count(r.gq) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            )
                    THEN min(r.gq) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            ) END AS reference_gq,
               CASE WHEN count(DISTINCT r.ploidy) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            ) = 1
                    THEN min(r.ploidy) FILTER (
                              WHERE r.state IN ('reference_block','explicit_reference_site')
                                AND r.start_pos <= t.pos AND r.end_pos >= t.end_pos
                            ) END AS reference_ploidy
        FROM hits
        JOIN targets t USING (callability_id)
        JOIN reference_evidence r USING (evidence_id)
        GROUP BY t.callability_id
        """
    )
    connection.execute(
        """
        CREATE TEMP TABLE region_hits AS
        SELECT t.callability_id,
               count(*)::BIGINT AS region_count,
               bool_and(r.callable IS TRUE) AS regions_callable,
               list(r.region_id ORDER BY r.region_id) AS region_ids
        FROM targets t
        JOIN callable_regions r
          ON t.chrom = r.chrom
         AND r.start_pos <= t.pos
         AND r.end_pos >= t.end_pos
        GROUP BY t.callability_id
        """
    )


def _minimum_complete(left: int | None, right: int | None) -> int | None:
    if left is None or right is None:
        return None
    return min(left, right)


def _resolve(
    row: dict[str, Any], config: BuildConfig, par_scope: Sequence[dict[str, Any]] = ()
) -> tuple[dict[str, Any], dict[str, Any]]:
    variant_count = int(row.get("variant_hit_count") or 0)
    non_core_count = int(row.get("non_core_hit_count") or 0)
    no_call_count = int(row.get("no_call_count") or 0)
    reference_count = int(row.get("reference_count") or 0)
    region_count = int(row.get("region_count") or 0)
    observed_variant_ids = list(row.get("observed_variant_ids") or [])
    source_digests = list(row.get("observed_variant_source_digests") or [])
    non_core_ids = list(row.get("non_core_evidence_ids") or [])
    no_call_ids = list(row.get("no_call_evidence_ids") or [])
    reference_ids = list(row.get("reference_evidence_ids") or [])
    region_ids = list(row.get("region_ids") or [])

    partial_hit = row.get("normalized_no_call") or row.get("partial_non_core_hit")
    if (variant_count or non_core_count) and not partial_hit:
        exact = bool(row.get("exact_target_observed")) or (
            row.get("variant_id") is None and variant_count > 0
        )
        resolution = "non_reference_observed" if exact else "overlapping_variant_observed"
        reference_observed = False
        if non_core_count:
            callable_value = False
            confidence = "low"
        else:
            callable_value = bool(row.get("variants_callable"))
            confidence = str(row["variant_confidence"])
        if variant_count and non_core_count:
            depth = _minimum_complete(row.get("variant_depth"), row.get("non_core_depth"))
            gq = _minimum_complete(row.get("variant_gq"), row.get("non_core_gq"))
            ploidy = (
                row.get("variant_ploidy")
                if row.get("variant_ploidy") == row.get("non_core_ploidy")
                else None
            )
        elif variant_count:
            depth, gq, ploidy = (
                row.get("variant_depth"),
                row.get("variant_gq"),
                row.get("variant_ploidy"),
            )
        else:
            depth, gq, ploidy = (
                row.get("non_core_depth"),
                row.get("non_core_gq"),
                row.get("non_core_ploidy"),
            )
        reason = (
            "the exact target allele or a carried site-level allele is observed"
            if exact
            else "a different or structurally represented carried allele overlaps the target"
        )
    elif no_call_count:
        resolution = "no_call"
        reference_observed = None
        callable_value = False
        confidence = None
        depth, gq, ploidy = (
            row.get("no_call_depth"),
            row.get("no_call_gq"),
            row.get("no_call_ploidy"),
        )
        reason = "overlapping incomplete genotype evidence forbids a reference claim"
    elif reference_count:
        callable_value = bool(row.get("references_callable"))
        resolution = "reference_observed" if callable_value else "low_confidence_reference"
        reference_observed = True
        confidence = str(row["reference_confidence"])
        depth, gq, ploidy = (
            row.get("reference_depth"),
            row.get("reference_gq"),
            row.get("reference_ploidy"),
        )
        reason = (
            "complete explicit reference evidence covers the target and passes thresholds"
            if callable_value
            else "complete reference evidence exists but fails one or more callability thresholds"
        )
    elif region_count:
        callable_value = bool(row.get("regions_callable"))
        resolution = (
            "callable_without_reference" if callable_value else "uncallable_without_reference"
        )
        reference_observed = None
        confidence = None
        depth = gq = ploidy = None
        reason = (
            "a covering interval is callable, but interval evidence is not a reference genotype"
            if callable_value
            else "a covering interval is uncallable and no locus-specific genotype was retained"
        )
    else:
        resolution = "unknown"
        reference_observed = callable_value = confidence = None
        depth = gq = ploidy = None
        reason = "no carried allele, no-call, complete reference, or covering interval was retained"

    limit = par_query_limit(par_scope, row["chrom"], row["pos"], row["end_pos"])
    if limit and (
        limit["scope"] == "all"
        or resolution
        not in {
            "non_reference_observed",
            "overlapping_variant_observed",
        }
    ):
        resolution = "mapping_scope_unassessed"
        reference_observed = confidence = None
        callable_value = False
        reason = limit["reason"]

    core = {
        "callability_id": row["callability_id"],
        "variant_id": row.get("variant_id"),
        "gene_symbol": row.get("gene_symbol"),
        "chrom": row["chrom"],
        "pos": row["pos"],
        "reference_observed": reference_observed,
        "callable": callable_value,
        "call_confidence": confidence,
        "depth": depth,
        "genotype_quality": float(gq) if gq is not None else None,
        "assay_scope": row["assay_scope"],
        "evidence_scope": row["evidence_scope"],
    }
    extension = {
        "callability_id": row["callability_id"],
        "target_id": row["target_id"],
        "target_source": row["target_source"],
        "target_source_version": row["target_source_version"],
        "catalog_sha256": row["catalog_sha256"],
        "catalog_row_digest": row["catalog_row_digest"],
        "catalog_row_number": row["catalog_row_number"],
        "chrom": row["chrom"],
        "pos": row["pos"],
        "end_pos": row["end_pos"],
        "ref": row.get("ref"),
        "alt": row.get("alt"),
        "variant_id": row.get("variant_id"),
        "gene_symbol": row.get("gene_symbol"),
        "resolution": resolution,
        "reference_observed": reference_observed,
        "callable": callable_value,
        "call_confidence": confidence,
        "depth": depth,
        "genotype_quality": float(gq) if gq is not None else None,
        "observed_ploidy": ploidy,
        "assay_scope": row["assay_scope"],
        "evidence_scope": row["evidence_scope"],
        "observed_variant_ids": observed_variant_ids,
        "observed_variant_source_digests": source_digests,
        "observed_non_core_evidence_ids": non_core_ids,
        "overlapping_no_call_evidence_ids": no_call_ids,
        "selected_reference_evidence_ids": reference_ids,
        "callable_region_ids": region_ids,
        "callable_min_dp": config.quality.callable_min_dp,
        "callable_min_gq": config.quality.callable_min_gq,
        "precedence_rule": PRECEDENCE_RULE,
        "reason": reason,
    }
    return core, extension


def _materialize(
    config: BuildConfig,
    connection: duckdb.DuckDBPyConnection,
    core_path: Path,
    extension_path: Path,
    par_scope: list[dict[str, Any]],
) -> Counter[str]:
    cursor = connection.execute(
        """
        SELECT t.*, EXISTS (SELECT 1 FROM normalized_no_calls p
                             WHERE p.callability_id = t.callability_id) AS normalized_no_call,
               v.variant_hit_count, v.exact_target_observed, v.variants_callable,
               v.variant_confidence, v.variant_depth, v.variant_gq, v.variant_ploidy,
               v.observed_variant_ids, v.observed_variant_source_digests,
               n.non_core_hit_count, n.non_core_depth, n.non_core_gq, n.non_core_ploidy,
               n.non_core_evidence_ids, n.partial_non_core_hit,
               r.no_call_count, r.no_call_evidence_ids, r.no_call_depth, r.no_call_gq,
               r.no_call_ploidy, r.reference_count, r.reference_evidence_ids,
               r.references_callable, r.reference_confidence, r.reference_depth,
               r.reference_gq, r.reference_ploidy,
               g.region_count, g.regions_callable, g.region_ids
        FROM targets t
        LEFT JOIN variant_hits v USING (callability_id)
        LEFT JOIN non_core_hits n USING (callability_id)
        LEFT JOIN reference_hits r USING (callability_id)
        LEFT JOIN region_hits g USING (callability_id)
        ORDER BY t.chrom_rank, t.pos, coalesce(t.ref,''), coalesce(t.alt,''), t.target_id
        """
    )
    names = [description[0] for description in cursor.description]
    counts: Counter[str] = Counter()
    with (
        pq.ParquetWriter(
            core_path,
            callability_schema(),
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
            write_statistics=True,
            version="2.6",
        ) as core_writer,
        pq.ParquetWriter(
            extension_path,
            target_callability_evidence_schema(),
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
            write_statistics=True,
            version="2.6",
        ) as extension_writer,
    ):
        while batch := cursor.fetchmany(50_000):
            core_rows: list[dict[str, Any]] = []
            extension_rows: list[dict[str, Any]] = []
            for values in batch:
                resolved, evidence = _resolve(
                    dict(zip(names, values, strict=True)), config, par_scope
                )
                core_rows.append(resolved)
                extension_rows.append(evidence)
                counts[evidence["resolution"]] += 1
            core_writer.write_table(
                pa.Table.from_pylist(core_rows, schema=callability_schema()),
                row_group_size=50_000,
            )
            extension_writer.write_table(
                pa.Table.from_pylist(extension_rows, schema=target_callability_evidence_schema()),
                row_group_size=50_000,
            )
    return counts


def build_target_callability(
    config: BuildConfig,
    *,
    ingest_dir: str | Path,
    variant_dir: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    """Build normative and lossless target-level callability tables."""

    if not config.target_catalog_enabled:
        raise InputValidationError("target callability was invoked while the catalogue is disabled")
    expected = config.target_catalog_expected_sha256
    if expected is None:
        raise InputValidationError("enabled target catalogue has no expected SHA-256")
    if not config.target_catalog_path.is_file():
        raise InputValidationError(f"target catalogue not found: {config.target_catalog_path}")
    actual = sha256_file(config.target_catalog_path)
    if actual != expected:
        raise InputValidationError(
            f"target catalogue SHA-256 mismatch: expected {expected}, got {actual}"
        )

    ingest_dir = Path(ingest_dir)
    variant_dir = Path(variant_dir)
    destination = Path(destination)
    with atomic_output_directory(destination) as staging:
        with tempfile.TemporaryDirectory(prefix=".target-work-", dir=staging) as temporary:
            work = Path(temporary)
            targets = work / "targets.parquet"
            raw_vcf = work / "targets.raw.vcf"
            normalized_vcf = work / "targets.normalized.vcf"
            target_rows, allele_rows = _write_catalog_table(config, actual, targets, raw_vcf)
            _validate_unique_targets(targets)
            _validate_normalized_alleles(config, raw_vcf, normalized_vcf, allele_rows)

            connection = duckdb.connect(str(work / "materialization.duckdb"))
            try:
                connection.execute(
                    f"SET temp_directory='{str(work).replace(chr(39), chr(39) * 2)}'"
                )
                _build_hit_tables(
                    connection,
                    targets=targets,
                    variants=variant_dir / "variants.parquet",
                    provenance=variant_dir / "variant_provenance.parquet",
                    partial_variants=variant_dir / "observed_partial_variants.parquet",
                    references=ingest_dir / "reference_evidence.parquet",
                    regions=ingest_dir / "callable_regions.parquet",
                    non_core=ingest_dir / "non_core_variants.parquet",
                )
                core_path = staging / "callability.parquet"
                extension_path = staging / "target_callability_evidence.parquet"
                ingest = json.loads((ingest_dir / "ingest_summary.json").read_bytes())
                scopes = ingest.get("family_par_query_scope", [])
                counts = _materialize(config, connection, core_path, extension_path, scopes)
            finally:
                connection.close()

        summary = {
            **config.checkpoint_identity,
            "status": "PASS",
            "catalog": {
                "sha256": actual,
                "bytes": config.target_catalog_path.stat().st_size,
                "target_rows": target_rows,
                "allele_rows": allele_rows,
                "normalization_checked": allele_rows,
                "assay_scope": config.target_catalog_assay_scope,
            },
            "resolution_counts": dict(sorted(counts.items())),
            "precedence_rule": PRECEDENCE_RULE,
            "outputs": {
                "callability.parquet": sha256_file(core_path),
                "target_callability_evidence.parquet": sha256_file(extension_path),
            },
        }
        atomic_write_json(staging / "target_callability_summary.json", summary)
    return summary
