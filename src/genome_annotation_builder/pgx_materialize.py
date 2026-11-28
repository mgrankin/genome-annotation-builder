"""Build a conservative, position-complete PharmCAT input VCF.

The input to this module is the already materialized target-evidence table.  A
missing row, an interval-only call, a different overlapping allele, or a
low-confidence call can therefore never become a reference genotype here.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .arrow_schemas import pgx_locus_audit_schema
from .config import BuildConfig
from .constants import PHARMCAT_POSITIONS_ALLELES, PHARMCAT_POSITIONS_RECORDS
from .errors import InputValidationError
from .util import atomic_output_directory, atomic_write_json, canonical_json_bytes, sha256_file
from .vcf_ingest import canonical_core_chrom

PASSING_CONFIDENCE = frozenset({"high", "medium"})


def _canonical_record_digest(values: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_json_bytes(values)).hexdigest()


def _consistent_non_null(values: Iterable[Any]) -> Any | None:
    observed = {value for value in values if value is not None}
    return next(iter(observed)) if len(observed) == 1 else None


def _no_call(
    base: dict[str, Any],
    *,
    status: str,
    reason: str,
    depth: int | None = None,
    genotype_quality: float | None = None,
    observed_ploidy: int | None = None,
) -> dict[str, Any]:
    row = {
        **base,
        "status": status,
        "genotype": None,
        "phased": None,
        "phase_set": None,
        "observed_ploidy": observed_ploidy,
        "depth": depth,
        "genotype_quality": genotype_quality,
        "eligible_for_matcher": False,
        "reason": reason,
    }
    row["materialized_record_digest"] = _canonical_record_digest(
        {
            "locus_id": row["locus_id"],
            "chrom": row["chrom"],
            "pos": row["pos"],
            "ref": row["ref"],
            "alts": row["alts"],
            "target_variant_ids": row["target_variant_ids"],
            "status": status,
            "genotype": None,
            "depth": depth,
            "genotype_quality": genotype_quality,
            "reason": reason,
        }
    )
    return row


def _resolve_pgx_locus(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Resolve all official ALT rows for one PharmCAT locus.

    This pure function is intentionally strict so its truth table can be unit
    tested without invoking PharmCAT or loading a genome-wide VCF.
    """

    if not rows:
        raise InputValidationError("cannot resolve an empty PharmCAT locus")
    ordered = sorted(rows, key=lambda item: int(item["allele_index"]))
    first = ordered[0]
    invariant_fields = (
        "locus_id",
        "record_number",
        "source_record_digest",
        "positions_vcf_sha256",
        "chrom",
        "pos",
        "source_id",
        "ref",
        "genes",
        "position_of_interest",
        "structurally_complex",
    )
    for row in ordered[1:]:
        for field in invariant_fields:
            if row.get(field) != first.get(field):
                raise InputValidationError(
                    f"inconsistent PharmCAT locus resource field {field}: {first['locus_id']}"
                )
    expected_indices = list(range(1, len(ordered) + 1))
    observed_indices = [int(row["allele_index"]) for row in ordered]
    if observed_indices != expected_indices:
        raise InputValidationError(
            f"non-contiguous PharmCAT allele indices for {first['locus_id']}"
        )

    alts = [str(row["alt"]) for row in ordered]
    target_ids = [str(row["target_id"]) for row in ordered]
    target_variant_ids = [str(row["target_variant_id"]) for row in ordered]
    normalization_changed = [bool(row["normalization_changed"]) for row in ordered]
    callability_ids = [str(row["callability_id"]) for row in ordered]
    resolutions = [str(row["resolution"]) for row in ordered]
    observed_variant_ids = sorted(
        {str(value) for row in ordered for value in (row.get("observed_variant_ids") or [])}
    )
    official_variant_ids = set(target_variant_ids)
    base = {
        "locus_id": str(first["locus_id"]),
        "record_number": int(first["record_number"]),
        "source_record_digest": str(first["source_record_digest"]),
        "positions_vcf_sha256": str(first["positions_vcf_sha256"]),
        "chrom": str(first["chrom"]),
        "pos": int(first["pos"]),
        "source_id": first.get("source_id"),
        "ref": str(first["ref"]),
        "alts": alts,
        "genes": list(first.get("genes") or []),
        "target_ids": target_ids,
        "target_variant_ids": target_variant_ids,
        "normalization_changed": normalization_changed,
        "callability_ids": callability_ids,
        "evidence_resolutions": resolutions,
        "observed_variant_ids": observed_variant_ids,
        "structurally_complex": bool(first["structurally_complex"]),
    }
    ploidy_hint = _consistent_non_null(row.get("observed_ploidy") for row in ordered)
    if ploidy_hint not in {1, 2}:
        ploidy_hint = None

    if base["structurally_complex"]:
        return _no_call(
            base,
            status="no_call_structurally_complex",
            reason=(
                "the locus belongs to a gene requiring an outside CNV/SV/HLA-aware call; "
                "SNV/indel VCF matching is deliberately disabled"
            ),
            observed_ploidy=ploidy_hint,
        )

    if "mapping_scope_unassessed" in resolutions:
        return _no_call(
            base,
            status="no_call_unknown",
            reason="the source mapping profile does not assess this target",
            observed_ploidy=ploidy_hint,
        )

    unexpected = sorted(set(observed_variant_ids) - official_variant_ids)
    if unexpected:
        return _no_call(
            base,
            status="no_call_other_variant",
            reason="a carried overlapping allele is not an exact PharmCAT target allele",
            observed_ploidy=ploidy_hint,
        )

    exact_rows = [row for row in ordered if row.get("variant_gt") is not None]
    if exact_rows:
        if any(
            row.get("callable") is not True or row.get("call_confidence") not in PASSING_CONFIDENCE
            for row in exact_rows
        ):
            return _no_call(
                base,
                status="no_call_low_confidence_variant",
                reason="an exact carried PharmCAT allele does not pass the DP/GQ/FILTER gate",
                depth=_consistent_non_null(row.get("depth") for row in exact_rows),
                genotype_quality=_consistent_non_null(
                    row.get("genotype_quality") for row in exact_rows
                ),
                observed_ploidy=ploidy_hint,
            )
        genotypes = [list(row["variant_gt"]) for row in exact_rows]
        lengths = {len(gt) for gt in genotypes}
        if lengths not in ({1}, {2}) or any(
            allele not in {0, 1} for gt in genotypes for allele in gt
        ):
            return _no_call(
                base,
                status="no_call_inconsistent_genotype",
                reason="exact target alleles have unsupported or incomplete projected genotypes",
                observed_ploidy=ploidy_hint,
            )
        ploidy = next(iter(lengths))
        genotype: list[int] = []
        for slot in range(ploidy):
            carried = [
                int(row["allele_index"])
                for row, gt in zip(exact_rows, genotypes, strict=True)
                if gt[slot] == 1
            ]
            if len(carried) > 1:
                return _no_call(
                    base,
                    status="no_call_inconsistent_genotype",
                    reason="multiple target alleles occupy the same genotype slot",
                    observed_ploidy=ploidy,
                )
            genotype.append(carried[0] if carried else 0)
        if not any(genotype):
            raise InputValidationError("exact PharmCAT allele join produced a reference genotype")
        phased_values = {bool(row.get("variant_phased")) for row in exact_rows}
        if len(phased_values) != 1:
            return _no_call(
                base,
                status="no_call_inconsistent_genotype",
                reason="exact target alleles disagree on phasing state",
                observed_ploidy=ploidy,
            )
        depth_values = [row.get("depth") for row in exact_rows]
        gq_values = [row.get("genotype_quality") for row in exact_rows]
        if any(value is None for value in depth_values + gq_values):
            return _no_call(
                base,
                status="no_call_low_confidence_variant",
                reason="an exact carried PharmCAT allele lacks complete DP/GQ evidence",
                observed_ploidy=ploidy,
            )
        depth = min(int(value) for value in depth_values)
        genotype_quality = min(float(value) for value in gq_values)
        result = {
            **base,
            "status": "variant_call",
            "genotype": genotype,
            "phased": next(iter(phased_values)),
            "phase_set": None,
            "observed_ploidy": ploidy,
            "depth": depth,
            "genotype_quality": genotype_quality,
            "eligible_for_matcher": True,
            "reason": "all carried alleles at the locus are exact passing PharmCAT targets",
        }
    else:
        resolution_set = set(resolutions)
        if resolution_set == {"reference_observed"} and all(
            row.get("reference_observed") is True
            and row.get("callable") is True
            and row.get("call_confidence") in PASSING_CONFIDENCE
            for row in ordered
        ):
            ploidy = _consistent_non_null(row.get("observed_ploidy") for row in ordered)
            depths = [row.get("depth") for row in ordered]
            qualities = [row.get("genotype_quality") for row in ordered]
            if ploidy in {1, 2} and not any(value is None for value in [*depths, *qualities]):
                result = {
                    **base,
                    "status": "reference_call",
                    "genotype": [0] * int(ploidy),
                    "phased": False,
                    "phase_set": None,
                    "observed_ploidy": int(ploidy),
                    "depth": min(int(value) for value in depths),
                    "genotype_quality": min(float(value) for value in qualities),
                    "eligible_for_matcher": True,
                    "reason": (
                        "complete explicit reference evidence covers every allele at the locus "
                        "and passes DP/GQ/FILTER thresholds"
                    ),
                }
            else:
                return _no_call(
                    base,
                    status="no_call_inconsistent_reference",
                    reason="reference target rows lack consistent complete ploidy/DP/GQ evidence",
                    observed_ploidy=ploidy_hint,
                )
        else:
            if "low_confidence_reference" in resolution_set:
                status = "no_call_low_confidence_reference"
                reason = "explicit reference evidence exists but fails the DP/GQ/FILTER gate"
            elif "no_call" in resolution_set:
                status = "no_call_source_missing"
                reason = "the source contains an overlapping incomplete genotype"
            elif resolution_set & {
                "callable_without_reference",
                "uncallable_without_reference",
            }:
                status = "no_call_interval_only"
                reason = "interval callability is not a locus-specific reference genotype"
            elif "overlapping_variant_observed" in resolution_set:
                status = "no_call_other_variant"
                reason = "an overlapping variant cannot be reconstructed as an exact target allele"
            else:
                status = "no_call_unknown"
                reason = "no exact variant or complete explicit reference evidence is available"
            return _no_call(
                base,
                status=status,
                reason=reason,
                observed_ploidy=ploidy_hint,
            )

    result["materialized_record_digest"] = _canonical_record_digest(
        {
            "locus_id": result["locus_id"],
            "chrom": result["chrom"],
            "pos": result["pos"],
            "ref": result["ref"],
            "alts": result["alts"],
            "target_variant_ids": result["target_variant_ids"],
            "status": result["status"],
            "genotype": result["genotype"],
            "phased": result["phased"],
            "depth": result["depth"],
            "genotype_quality": result["genotype_quality"],
        }
    )
    return result


def _parquet_source(path: Path) -> str:
    source = path / "**" / "*.parquet" if path.is_dir() else path
    return str(source).replace("'", "''")


def _locus_cursor(
    resource_dir: Path, target_callability_dir: Path, variant_dir: Path
) -> tuple[duckdb.DuckDBPyConnection, duckdb.DuckDBPyConnection]:
    connection = duckdb.connect(":memory:")
    loci = _parquet_source(resource_dir / "pgx_loci.parquet")
    evidence = _parquet_source(target_callability_dir / "target_callability_evidence.parquet")
    variants = _parquet_source(variant_dir / "variants.parquet")
    cursor = connection.execute(
        f"""
        SELECT l.*,
               e.callability_id, e.resolution, e.reference_observed, e.callable,
               e.call_confidence, e.depth, e.genotype_quality, e.observed_ploidy,
               e.observed_variant_ids,
               v.genotype.gt AS variant_gt,
               v.genotype.phased AS variant_phased
        FROM read_parquet('{loci}') l
        JOIN read_parquet('{evidence}') e USING (target_id)
        LEFT JOIN read_parquet('{variants}') v ON v.variant_id = l.target_variant_id
        ORDER BY l.record_number, l.allele_index
        """
    )
    return connection, cursor


def _output_header(config: BuildConfig) -> pysam.VariantHeader:
    fasta = pysam.FastaFile(str(config.reference_fasta))
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value="genome-annotation-builder")
    reference_identity = (
        f"sha256:{config.reference_expected_sha256}"
        if config.reference_expected_sha256
        else "GRCh38"
    )
    header.add_meta("reference", value=reference_identity)
    try:
        for name, length in zip(fasta.references, fasta.lengths, strict=True):
            canonical = canonical_core_chrom(name)
            if canonical is not None and canonical not in header.contigs:
                header.contigs.add(canonical, length=int(length))
    finally:
        fasta.close()
    header.add_meta(
        "FILTER",
        items=[
            ("ID", "GAB_PGX_NO_CALL"),
            ("Description", "No eligible explicit genotype was materialized for PharmCAT"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[("ID", "PX"), ("Number", "."), ("Type", "String"), ("Description", "PharmCAT gene")],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "POI"),
            ("Number", "0"),
            ("Type", "Flag"),
            ("Description", "PharmCAT position of interest"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "GAB_PGX_STATUS"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Evidence materialization status"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "GAB_PGX_EVIDENCE"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "SHA-256 of the materialization decision"),
        ],
    )
    header.formats.add(
        "GT", number=1, type="String", description="Explicit target genotype or no-call"
    )
    header.formats.add("DP", number=1, type="Integer", description="Conservative target depth")
    header.formats.add(
        "GQ", number=1, type="Integer", description="Conservative target genotype quality"
    )
    header.add_sample("SAMPLE")
    return header


def materialize_pgx_input(
    config: BuildConfig,
    *,
    resource_dir: str | Path,
    target_callability_dir: str | Path,
    variant_dir: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    """Create a position-complete VCF and per-locus evidence audit."""

    if not config.pgx.enabled:
        raise InputValidationError("PGx materialization was invoked while PGx is disabled")
    resource_dir = Path(resource_dir)
    target_callability_dir = Path(target_callability_dir)
    variant_dir = Path(variant_dir)
    destination = Path(destination)
    with atomic_output_directory(destination) as staging:
        audit_path = staging / "pgx_locus_audit.parquet"
        vcf_path = staging / "pgx-input.vcf.gz"
        writer = pq.ParquetWriter(
            audit_path,
            pgx_locus_audit_schema(),
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
            write_statistics=True,
            version="2.6",
        )
        vcf = pysam.VariantFile(str(vcf_path), "wz", header=_output_header(config))
        counts: Counter[str] = Counter()
        record_count = 0
        allele_count = 0
        connection: duckdb.DuckDBPyConnection | None = None
        try:
            connection, cursor = _locus_cursor(resource_dir, target_callability_dir, variant_dir)
            names = [description[0] for description in cursor.description]
            pending: list[dict[str, Any]] = []
            current_record: int | None = None

            def emit(rows: list[dict[str, Any]]) -> None:
                nonlocal record_count, allele_count
                audit = _resolve_pgx_locus(rows)
                writer.write_table(pa.Table.from_pylist([audit], schema=pgx_locus_audit_schema()))
                counts[audit["status"]] += 1
                record_count += 1
                allele_count += len(audit["alts"])
                record = vcf.new_record(
                    contig=audit["chrom"],
                    start=int(audit["pos"]) - 1,
                    stop=int(audit["pos"]) - 1 + len(audit["ref"]),
                    id=audit["source_id"],
                    alleles=(audit["ref"], *audit["alts"]),
                )
                if audit["eligible_for_matcher"]:
                    record.filter.add("PASS")
                else:
                    record.filter.add("GAB_PGX_NO_CALL")
                if audit["genes"]:
                    record.info["PX"] = tuple(audit["genes"])
                if any(bool(row.get("position_of_interest")) for row in rows):
                    record.info["POI"] = True
                record.info["GAB_PGX_STATUS"] = audit["status"]
                record.info["GAB_PGX_EVIDENCE"] = audit["materialized_record_digest"]
                sample = record.samples["SAMPLE"]
                if audit["genotype"] is not None:
                    sample["GT"] = tuple(audit["genotype"])
                    sample.phased = bool(audit["phased"])
                else:
                    # VCF has no ploidy-neutral missing GT.  Preserve an observed
                    # one/two-copy hint when available; otherwise use a diploid
                    # missing representation without making a reference claim.
                    missing_ploidy = audit["observed_ploidy"] or 2
                    sample["GT"] = tuple(None for _ in range(int(missing_ploidy)))
                if audit["depth"] is not None:
                    sample["DP"] = int(audit["depth"])
                if audit["genotype_quality"] is not None:
                    sample["GQ"] = int(audit["genotype_quality"])
                vcf.write(record)

            while values := cursor.fetchmany(10_000):
                for value in values:
                    row = dict(zip(names, value, strict=True))
                    record_number = int(row["record_number"])
                    if current_record is None:
                        current_record = record_number
                    if record_number != current_record:
                        emit(pending)
                        pending = []
                        current_record = record_number
                    pending.append(row)
            if pending:
                emit(pending)
        finally:
            vcf.close()
            writer.close()
            if connection is not None:
                connection.close()

        if (record_count, allele_count) != (
            PHARMCAT_POSITIONS_RECORDS,
            PHARMCAT_POSITIONS_ALLELES,
        ):
            raise InputValidationError(
                "PGx materialization accounting mismatch: "
                f"{record_count} records/{allele_count} alleles"
            )
        try:
            pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
        except OSError as exc:
            raise InputValidationError(f"failed to index materialized PGx VCF: {exc}") from exc
        index_path = Path(f"{vcf_path}.tbi")
        summary = {
            **config.checkpoint_identity,
            "status": "PASS",
            "record_count": record_count,
            "allele_count": allele_count,
            "eligible_records": counts["variant_call"] + counts["reference_call"],
            "no_call_records": record_count - counts["variant_call"] - counts["reference_call"],
            "status_counts": dict(sorted(counts.items())),
            "sample_name": "SAMPLE",
            "missing_to_reference_conversion": False,
            "structurally_complex_matcher_enabled": False,
            "outputs": {
                "pgx-input.vcf.gz": sha256_file(vcf_path),
                "pgx-input.vcf.gz.tbi": sha256_file(index_path),
                "pgx_locus_audit.parquet": sha256_file(audit_path),
            },
        }
        atomic_write_json(staging / "pgx_materialization_summary.json", summary)
    return summary
