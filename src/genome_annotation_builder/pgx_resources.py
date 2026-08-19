"""Immutable PharmCAT positions and CPIC pair-snapshot preparation."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import tempfile
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .arrow_schemas import cpic_pair_resource_schema, pgx_locus_resource_schema
from .config import BuildConfig
from .constants import (
    PHARMCAT_POSITIONS_ALLELES,
    PHARMCAT_POSITIONS_BYTES,
    PHARMCAT_POSITIONS_RECORDS,
    PHARMCAT_POSITIONS_URL,
    PHARMCAT_RELEASE_DATE,
    PHARMCAT_VERSION,
    STANDARD_CHROMS,
    STRUCTURALLY_COMPLEX_PGX_GENES,
)
from .errors import ResourceError
from .target_callability import CATALOG_COLUMNS
from .target_union import merge_target_catalogs
from .util import (
    atomic_output_directory,
    atomic_write_bytes,
    atomic_write_json,
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
)
from .vcf_ingest import canonical_core_chrom

PHARMCAT_SOURCE = "PharmCAT"
CPIC_SOURCE = "CPIC PostgREST pair_view"
CPIC_SOURCE_URL = "https://api.cpicpgx.org/v1/pair_view"
CPIC_FIELDS = (
    "pairid",
    "drugid",
    "drugname",
    "genesymbol",
    "guidelinename",
    "guidelineurl",
    "cpiclevel",
    "clinpgxlevel",
    "pgxtesting",
    "pmids",
    "usedforrecommendation",
    "provisional",
)
CPIC_NORMATIVE_LEVELS = frozenset({"A", "B", "C", "D"})
_CHROM_ORDER = {chrom: index for index, chrom in enumerate(STANDARD_CHROMS)}


def _reference_contigs(fasta: pysam.FastaFile) -> dict[str, tuple[str, int]]:
    output: dict[str, tuple[str, int]] = {}
    for name, length in zip(fasta.references, fasta.lengths, strict=True):
        canonical = canonical_core_chrom(name)
        if canonical is None:
            continue
        if canonical in output:
            raise ResourceError(f"reference has ambiguous aliases for {canonical}")
        output[canonical] = (name, int(length))
    return output


def _genes(value: Any) -> list[str]:
    if value is None:
        return []
    values = [value] if isinstance(value, str) else list(value)
    return sorted({str(item).strip() for item in values if str(item).strip()})


def _normalize_position_rows(
    config: BuildConfig,
    rows: list[dict[str, Any]],
    reference_contigs: dict[str, tuple[str, int]],
) -> tuple[list[dict[str, Any]], int]:
    """Bind each source PharmCAT ALT to its exact bcftools-normalized target."""

    with tempfile.TemporaryDirectory(prefix="gab-pgx-normalization-") as temporary_name:
        temporary = Path(temporary_name)
        raw_path = temporary / "source.vcf"
        normalized_path = temporary / "normalized.vcf"
        header = pysam.VariantHeader()
        header.add_meta("fileformat", value="VCFv4.3")
        for chrom in sorted(
            {str(row["chrom"]) for row in rows}, key=lambda value: _CHROM_ORDER[value]
        ):
            reference_name, length = reference_contigs[chrom]
            header.contigs.add(reference_name, length=length)
        writer = pysam.VariantFile(str(raw_path), "w", header=header)
        try:
            for row in rows:
                reference_name = reference_contigs[str(row["chrom"])][0]
                record = writer.new_record(
                    contig=reference_name,
                    start=int(row["pos"]) - 1,
                    id=str(row["target_id"]),
                    alleles=(str(row["ref"]), str(row["alt"])),
                )
                writer.write(record)
        finally:
            writer.close()
        try:
            completed = subprocess.run(
                [
                    config.bcftools,
                    "norm",
                    "--no-version",
                    "--check-ref",
                    "e",
                    "--fasta-ref",
                    str(config.reference_fasta),
                    "--output-type",
                    "v",
                    "--output",
                    str(normalized_path),
                    str(raw_path),
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
            raise ResourceError("failed to execute bcftools for PharmCAT normalization") from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout).strip()
            raise ResourceError(f"PharmCAT allele normalization failed: {detail[-4000:]}")

        observed: dict[str, str] = {}
        normalized = pysam.VariantFile(str(normalized_path))
        try:
            for record in normalized:
                target_id = str(record.id or "")
                chrom = canonical_core_chrom(record.contig)
                alts = tuple(str(alt).upper() for alt in (record.alts or ()))
                if (
                    not target_id
                    or target_id in observed
                    or chrom is None
                    or len(alts) != 1
                    or not record.ref
                ):
                    raise ResourceError("bcftools produced an invalid PharmCAT normalization map")
                observed[target_id] = (
                    f"{chrom}:{int(record.pos)}:{str(record.ref).upper()}:{alts[0]}"
                )
        finally:
            normalized.close()
    if set(observed) != {str(row["target_id"]) for row in rows}:
        raise ResourceError("PharmCAT normalization map has missing or extra target identifiers")
    changed = 0
    for row in rows:
        source_variant = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
        row["target_variant_id"] = observed[str(row["target_id"])]
        row["normalization_changed"] = row["target_variant_id"] != source_variant
        changed += int(row["normalization_changed"])
    return rows, changed


def _parse_positions(
    config: BuildConfig,
    *,
    locus_writer: pq.ParquetWriter | None = None,
    target_writer: csv.writer | None = None,
) -> dict[str, Any]:
    path = config.pgx.positions_vcf
    if not path.is_file():
        raise ResourceError("PharmCAT positions VCF is missing")
    actual_sha256 = sha256_file(path)
    if actual_sha256 != config.pgx.positions_expected_sha256:
        raise ResourceError(
            "PharmCAT positions SHA-256 mismatch: "
            f"expected {config.pgx.positions_expected_sha256}, got {actual_sha256}"
        )
    if path.stat().st_size != PHARMCAT_POSITIONS_BYTES:
        raise ResourceError(
            f"PharmCAT {PHARMCAT_VERSION} positions size mismatch: {path.stat().st_size}"
        )

    fasta = pysam.FastaFile(str(config.reference_fasta))
    reference_contigs = _reference_contigs(fasta)
    source = pysam.VariantFile(str(path))
    rows: list[dict[str, Any]] = []
    record_count = 0
    allele_count = 0
    gene_counts: Counter[str] = Counter()
    complex_loci = 0
    previous_key: tuple[int, int, str, tuple[str, ...]] | None = None
    seen_loci: set[tuple[str, int, str]] = set()
    try:
        if tuple(source.header.samples) != ("PharmCAT",):
            raise ResourceError("PharmCAT positions VCF must contain exactly its PharmCAT sample")
        for record_number, record in enumerate(source, start=1):
            chrom = canonical_core_chrom(record.contig)
            if chrom is None or chrom not in reference_contigs or chrom not in _CHROM_ORDER:
                raise ResourceError(
                    f"PharmCAT positions record {record_number} uses unsupported contig"
                )
            ref = record.ref.upper()
            alts = tuple(str(alt).upper() for alt in (record.alts or ()))
            if not alts or any(not alt or set(alt) - set("ACGT") for alt in alts):
                raise ResourceError(
                    f"PharmCAT positions record {record_number} has non-concrete ALT"
                )
            if not ref or set(ref) - set("ACGT") or any(alt == ref for alt in alts):
                raise ResourceError(
                    f"PharmCAT positions record {record_number} has invalid alleles"
                )
            locus_key = (chrom, int(record.pos), ref)
            if locus_key in seen_loci:
                raise ResourceError(f"duplicate PharmCAT positions locus: {locus_key}")
            seen_loci.add(locus_key)
            sort_key = (_CHROM_ORDER[chrom], int(record.pos), ref, alts)
            if previous_key is not None and sort_key <= previous_key:
                raise ResourceError("PharmCAT positions VCF is not strictly coordinate sorted")
            previous_key = sort_key
            reference_name, reference_length = reference_contigs[chrom]
            if record.pos < 1 or record.pos + len(ref) - 1 > reference_length:
                raise ResourceError(f"PharmCAT positions record {record_number} is outside FASTA")
            observed_ref = fasta.fetch(reference_name, record.pos - 1, record.pos - 1 + len(ref))
            if observed_ref.upper() != ref:
                raise ResourceError(
                    f"PharmCAT positions REF mismatch at {chrom}:{record.pos}: "
                    f"{ref} != {observed_ref.upper()}"
                )
            genes = _genes(record.info.get("PX"))
            position_of_interest = bool(record.info.get("POI", False))
            if not genes and not position_of_interest:
                raise ResourceError(
                    f"PharmCAT positions record {record_number} has neither PX gene nor POI flag"
                )
            for gene in genes:
                if not re.fullmatch(r"[A-Za-z0-9-]+", gene):
                    raise ResourceError(f"invalid PharmCAT gene symbol at record {record_number}")
                gene_counts[gene] += 1
            structurally_complex = bool(set(genes) & STRUCTURALLY_COMPLEX_PGX_GENES)
            complex_loci += int(structurally_complex)
            locus_id = f"pgx:pharmcat:{PHARMCAT_VERSION}:record:{record_number:09d}"
            source_identity = {
                "record_number": record_number,
                "chrom": chrom,
                "pos": int(record.pos),
                "id": None if record.id in {None, "."} else str(record.id),
                "ref": ref,
                "alts": alts,
                "genes": genes,
                "position_of_interest": position_of_interest,
            }
            source_digest = hashlib.sha256(canonical_json_bytes(source_identity)).hexdigest()
            for allele_index, alt in enumerate(alts, start=1):
                target_id = f"{locus_id}:allele:{allele_index:03d}"
                row = {
                    "locus_id": locus_id,
                    "record_number": record_number,
                    "allele_index": allele_index,
                    "source_record_digest": source_digest,
                    "positions_vcf_sha256": actual_sha256,
                    "target_id": target_id,
                    "chrom": chrom,
                    "pos": int(record.pos),
                    "source_id": None if record.id in {None, "."} else str(record.id),
                    "ref": ref,
                    "alt": alt,
                    "genes": genes,
                    "position_of_interest": position_of_interest,
                    "structurally_complex": structurally_complex,
                }
                rows.append(row)
                allele_count += 1
            record_count += 1
    finally:
        source.close()
        fasta.close()
    if (record_count, allele_count) != (
        PHARMCAT_POSITIONS_RECORDS,
        PHARMCAT_POSITIONS_ALLELES,
    ):
        raise ResourceError(
            f"PharmCAT {PHARMCAT_VERSION} positions count mismatch: "
            f"{record_count} records/{allele_count} alleles"
        )
    rows, normalization_changed = _normalize_position_rows(config, rows, reference_contigs)
    if locus_writer is not None:
        locus_writer.write_table(
            pa.Table.from_pylist(rows, schema=pgx_locus_resource_schema()),
            row_group_size=10_000,
        )
    if target_writer is not None:
        for row in rows:
            target_chrom, target_pos, target_ref, target_alt = str(row["target_variant_id"]).split(
                ":", 3
            )
            target_writer.writerow(
                [
                    row["target_id"],
                    target_chrom,
                    target_pos,
                    target_ref,
                    target_alt,
                    ",".join(str(gene) for gene in row["genes"]),
                    PHARMCAT_SOURCE,
                    f"{PHARMCAT_VERSION}@{config.pgx.release_commit}",
                    "target_panel",
                ]
            )
    return {
        "version": PHARMCAT_VERSION,
        "release_commit": config.pgx.release_commit,
        "release_date": PHARMCAT_RELEASE_DATE,
        "positions_vcf_sha256": actual_sha256,
        "positions_vcf_bytes": path.stat().st_size,
        "record_count": record_count,
        "allele_count": allele_count,
        "gene_count": len(gene_counts),
        "genes": sorted(gene_counts),
        "structurally_complex_loci": complex_loci,
        "normalization_changed_alleles": normalization_changed,
        "license": config.pgx.license,
        "redistribution_allowed": config.pgx.redistribution_allowed,
        "container_image": config.pgx.container_image,
    }


def _nullable_text(value: Any) -> str | None:
    if value is None:
        return None
    output = str(value).strip()
    return output or None


def _parse_cpic(config: BuildConfig, *, writer: pq.ParquetWriter | None = None) -> dict[str, Any]:
    path = config.pgx.cpic_pairs
    if not path.is_file():
        raise ResourceError("CPIC pair snapshot is missing")
    actual_sha256 = sha256_file(path)
    if actual_sha256 != config.pgx.cpic_pairs_expected_sha256:
        raise ResourceError(
            "CPIC pair snapshot SHA-256 mismatch: "
            f"expected {config.pgx.cpic_pairs_expected_sha256}, got {actual_sha256}"
        )
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ResourceError("CPIC pair snapshot is not valid UTF-8 JSON") from exc
    if not isinstance(document, dict) or document.get("schema_version") != 2:
        raise ResourceError("unsupported CPIC pair snapshot schema")
    expected_document_fields = {
        "schema_version",
        "source",
        "source_url",
        "source_version",
        "snapshot_date",
        "fields",
        "row_count",
        "rowset_canonical_sha256",
        "stable_retrieval_passes",
        "source_release",
        "api_identity",
        "rows",
    }
    if set(document) != expected_document_fields:
        raise ResourceError("CPIC pair snapshot has an unexpected document boundary")
    if document.get("source") != CPIC_SOURCE:
        raise ResourceError("CPIC pair snapshot source mismatch")
    if document.get("source_url") != CPIC_SOURCE_URL:
        raise ResourceError("CPIC pair snapshot source URL mismatch")
    if document.get("source_version") != config.pgx.cpic_source_version:
        raise ResourceError("CPIC pair snapshot source_version mismatch")
    if document.get("snapshot_date") != config.pgx.cpic_snapshot_date.isoformat():
        raise ResourceError("CPIC pair snapshot date mismatch")
    if tuple(document.get("fields", ())) != CPIC_FIELDS:
        raise ResourceError("CPIC pair snapshot field boundary mismatch")
    raw_rows = document.get("rows")
    if not isinstance(raw_rows, list) or not raw_rows:
        raise ResourceError("CPIC pair snapshot rows must be a non-empty list")
    if document.get("row_count") != len(raw_rows):
        raise ResourceError("CPIC pair snapshot row_count mismatch")
    if document.get("stable_retrieval_passes") != 2:
        raise ResourceError("CPIC pair snapshot was not retrieved in two stable passes")
    rowset_sha256 = document.get("rowset_canonical_sha256")
    if not isinstance(rowset_sha256, str) or re.fullmatch(r"[0-9a-f]{64}", rowset_sha256) is None:
        raise ResourceError("CPIC pair snapshot has an invalid row-set digest")
    actual_rowset_sha256 = sha256_bytes(canonical_json_bytes(raw_rows))
    if actual_rowset_sha256 != rowset_sha256:
        raise ResourceError("CPIC pair snapshot row-set digest mismatch")

    source_release = document.get("source_release")
    expected_release_fields = {
        "tag",
        "commit",
        "github_release_id",
        "published_at",
        "release_url",
    }
    if not isinstance(source_release, dict) or set(source_release) != expected_release_fields:
        raise ResourceError("CPIC pair snapshot has an invalid release identity")
    expected_tag = f"v{config.pgx.cpic_source_version}"
    if source_release.get("tag") != expected_tag:
        raise ResourceError("CPIC pair snapshot release tag mismatch")
    release_commit = source_release.get("commit")
    if not isinstance(release_commit, str) or re.fullmatch(r"[0-9a-f]{40}", release_commit) is None:
        raise ResourceError("CPIC pair snapshot release commit is invalid")
    if not isinstance(source_release.get("github_release_id"), int):
        raise ResourceError("CPIC pair snapshot release ID is invalid")
    if (
        not isinstance(source_release.get("published_at"), str)
        or not source_release["published_at"]
    ):
        raise ResourceError("CPIC pair snapshot release timestamp is invalid")
    expected_release_url = "https://github.com/cpicpgx/cpic-data/releases/tag/" + expected_tag
    if source_release.get("release_url") != expected_release_url:
        raise ResourceError("CPIC pair snapshot release URL mismatch")

    api_identity = document.get("api_identity")
    expected_api_fields = {
        "openapi_canonical_sha256",
        "postgrest_version",
        "swagger_version",
        "base_path",
    }
    if not isinstance(api_identity, dict) or set(api_identity) != expected_api_fields:
        raise ResourceError("CPIC pair snapshot has an invalid API identity")
    openapi_sha256 = api_identity.get("openapi_canonical_sha256")
    if not isinstance(openapi_sha256, str) or re.fullmatch(r"[0-9a-f]{64}", openapi_sha256) is None:
        raise ResourceError("CPIC pair snapshot OpenAPI digest is invalid")
    if api_identity.get("swagger_version") != "2.0" or api_identity.get("base_path") != "/v1/":
        raise ResourceError("CPIC pair snapshot API boundary mismatch")
    if (
        not isinstance(api_identity.get("postgrest_version"), str)
        or not api_identity["postgrest_version"]
    ):
        raise ResourceError("CPIC pair snapshot serving software version is invalid")

    output_rows: list[dict[str, Any]] = []
    seen: set[int] = set()
    previous = -1
    level_counts: Counter[str] = Counter()
    recommendation_rows = 0
    for index, raw in enumerate(raw_rows, start=1):
        if not isinstance(raw, dict) or set(raw) != set(CPIC_FIELDS):
            raise ResourceError(f"CPIC pair row {index} has an unexpected shape")
        try:
            pair_id = int(raw["pairid"])
        except (TypeError, ValueError) as exc:
            raise ResourceError(f"CPIC pair row {index} has invalid pairid") from exc
        if pair_id in seen or pair_id <= previous:
            raise ResourceError("CPIC pair snapshot is not unique and sorted by pairid")
        seen.add(pair_id)
        previous = pair_id
        gene = str(raw["genesymbol"] or "").strip()
        drug = str(raw["drugname"] or "").strip()
        if not gene or not drug:
            raise ResourceError(f"CPIC pair row {index} lacks gene or drug")
        pmids = raw["pmids"] if isinstance(raw["pmids"], list) else []
        used_text = str(raw["usedforrecommendation"] or "").strip()
        if used_text not in {"Yes", "No", "n/a"}:
            raise ResourceError(f"CPIC pair row {index} has invalid recommendation flag")
        if not isinstance(raw["provisional"], bool):
            raise ResourceError(f"CPIC pair row {index} has invalid provisional flag")
        provisional = bool(raw["provisional"])
        if (used_text == "n/a") != provisional:
            raise ResourceError(
                f"CPIC pair row {index} has inconsistent provisional/recommendation status"
            )
        level = _nullable_text(raw["cpiclevel"])
        if level is not None:
            level_counts[level] += 1
        used = {"Yes": True, "No": False, "n/a": None}[used_text]
        recommendation_rows += int(
            used is True and not provisional and level in CPIC_NORMATIVE_LEVELS
        )
        output_rows.append(
            {
                "pair_id": pair_id,
                "drug_id": _nullable_text(raw["drugid"]),
                "drug_name": drug,
                "gene_symbol": gene,
                "guideline_name": _nullable_text(raw["guidelinename"]),
                "guideline_url": _nullable_text(raw["guidelineurl"]),
                "cpic_level": level,
                "clinpgx_level": _nullable_text(raw["clinpgxlevel"]),
                "pgx_testing": _nullable_text(raw["pgxtesting"]),
                "pmids": sorted({str(value) for value in pmids if str(value)}),
                "used_for_recommendation": used,
                "used_for_recommendation_status": used_text,
                "provisional": provisional,
            }
        )
        if writer is not None and len(output_rows) >= 10_000:
            writer.write_table(
                pa.Table.from_pylist(output_rows, schema=cpic_pair_resource_schema()),
                row_group_size=10_000,
            )
            output_rows.clear()
    if writer is not None and output_rows:
        writer.write_table(
            pa.Table.from_pylist(output_rows, schema=cpic_pair_resource_schema()),
            row_group_size=10_000,
        )
    if recommendation_rows == 0:
        raise ResourceError("CPIC pair snapshot has no non-provisional A-D recommendation rows")
    return {
        "source": CPIC_SOURCE,
        "source_version": config.pgx.cpic_source_version,
        "snapshot_date": config.pgx.cpic_snapshot_date.isoformat(),
        "sha256": actual_sha256,
        "rowset_canonical_sha256": rowset_sha256,
        "bytes": path.stat().st_size,
        "row_count": len(raw_rows),
        "normative_recommendation_rows": recommendation_rows,
        "level_counts": dict(sorted(level_counts.items())),
        "source_release": source_release,
        "api_identity": api_identity,
        "license": config.pgx.cpic_license,
        "redistribution_allowed": config.pgx.cpic_redistribution_allowed,
    }


def verify_pgx_resources(config: BuildConfig) -> dict[str, Any]:
    if not config.pgx.enabled:
        return {"status": "not_enabled"}
    return {
        "status": "verified",
        "pharmcat": _parse_positions(config),
        "cpic": _parse_cpic(config),
    }


def sync_pgx_resources(config: BuildConfig) -> dict[str, Any]:
    """Install only the immutable PharmCAT positions asset, then verify all PGx data."""

    if not config.pgx.enabled:
        return {"status": "not_enabled"}
    positions = config.pgx.positions_vcf
    acquisition = "reused"
    if not positions.is_file() or sha256_file(positions) != config.pgx.positions_expected_sha256:
        if not config.allow_download:
            raise ResourceError(
                "PharmCAT positions are absent or invalid and resources.allow_download is false"
            )
        request = urllib.request.Request(
            PHARMCAT_POSITIONS_URL,
            headers={"User-Agent": "genome-annotation-builder-pgx-sync/0.1"},
        )
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                if response.status != 200:
                    raise ResourceError(
                        f"unexpected HTTP status {response.status} for PharmCAT positions"
                    )
                payload = response.read()
        except (urllib.error.URLError, TimeoutError) as exc:
            raise ResourceError(f"failed to download pinned PharmCAT positions: {exc}") from exc
        actual = sha256_bytes(payload)
        if actual != config.pgx.positions_expected_sha256:
            raise ResourceError(
                "downloaded PharmCAT positions digest mismatch: "
                f"expected {config.pgx.positions_expected_sha256}, got {actual}"
            )
        positions.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_bytes(positions, payload)
        acquisition = "downloaded"
    verified = verify_pgx_resources(config)
    verified["positions_acquisition"] = acquisition
    return verified


def prepare_pgx_resources(
    config: BuildConfig,
    destination: str | Path,
    *,
    base_target_catalog: str | Path | None = None,
    base_target_label: str = "upstream target union",
) -> dict[str, Any]:
    if not config.pgx.enabled:
        raise ResourceError("PGx preparation was invoked while annotation.pgx is disabled")
    destination = Path(destination)
    with atomic_output_directory(destination) as staging:
        loci_path = staging / "pgx_loci.parquet"
        cpic_path = staging / "cpic_pairs.parquet"
        pgx_targets = staging / "pgx_targets.tsv"
        with (
            pq.ParquetWriter(
                loci_path,
                pgx_locus_resource_schema(),
                compression="zstd",
                compression_level=9,
                use_dictionary=True,
                write_statistics=True,
                version="2.6",
            ) as locus_writer,
            pq.ParquetWriter(
                cpic_path,
                cpic_pair_resource_schema(),
                compression="zstd",
                compression_level=9,
                use_dictionary=True,
                write_statistics=True,
                version="2.6",
            ) as cpic_writer,
            pgx_targets.open("w", encoding="utf-8", newline="") as target_handle,
        ):
            target_writer = csv.writer(target_handle, delimiter="\t", lineterminator="\n")
            target_writer.writerow(CATALOG_COLUMNS)
            pharmcat_summary = _parse_positions(
                config, locus_writer=locus_writer, target_writer=target_writer
            )
            cpic_summary = _parse_cpic(config, writer=cpic_writer)

        sources: list[tuple[str, Path]] = []
        if base_target_catalog is not None:
            sources.append((base_target_label, Path(base_target_catalog)))
        elif config.target_catalog_enabled:
            sources.append(("configured base", config.target_catalog_path))
        sources.append(("generated PGx", pgx_targets))
        combined_targets = staging / "combined_targets.tsv"
        counts = merge_target_catalogs(sources, combined_targets)
        combined_sha256 = sha256_file(combined_targets)
        summary = {
            "status": "PASS",
            "pharmcat": pharmcat_summary,
            "cpic": cpic_summary,
            "base_target_catalog_sha256": (
                sha256_file(Path(base_target_catalog))
                if base_target_catalog is not None
                else sha256_file(config.target_catalog_path)
                if config.target_catalog_enabled
                else None
            ),
            "target_source_counts": counts,
            "combined_target_rows": sum(counts.values()),
            "combined_target_catalog_sha256": combined_sha256,
            "outputs": {
                "pgx_loci.parquet": sha256_file(loci_path),
                "cpic_pairs.parquet": sha256_file(cpic_path),
                "pgx_targets.tsv": sha256_file(pgx_targets),
                "combined_targets.tsv": combined_sha256,
            },
        }
        atomic_write_json(staging / "pgx_resource_summary.json", summary)
    return summary
