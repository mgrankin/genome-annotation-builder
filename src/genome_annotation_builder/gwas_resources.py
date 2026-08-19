"""Freeze and verify a dated NHGRI-EBI GWAS Catalog association release."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import sqlite3
import subprocess
import tempfile
import zipfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .arrow_schemas import gwas_catalog_resource_schema
from .config import BuildConfig
from .constants import (
    GWAS_CATALOG_ARCHIVE_BYTES,
    GWAS_CATALOG_ARCHIVE_SHA256,
    GWAS_CATALOG_ATTRIBUTION,
    GWAS_CATALOG_MEMBER,
    GWAS_CATALOG_MEMBER_BYTES,
    GWAS_CATALOG_SNAPSHOT,
    GWAS_CATALOG_URL,
)
from .errors import ResourceError
from .util import atomic_write_json, sha256_file

EXPECTED_HEADER = (
    "DATE ADDED TO CATALOG",
    "PUBMEDID",
    "FIRST AUTHOR",
    "DATE",
    "JOURNAL",
    "LINK",
    "STUDY",
    "DISEASE/TRAIT",
    "INITIAL SAMPLE SIZE",
    "REPLICATION SAMPLE SIZE",
    "REGION",
    "CHR_ID",
    "CHR_POS",
    "REPORTED GENE(S)",
    "MAPPED_GENE",
    "UPSTREAM_GENE_ID",
    "DOWNSTREAM_GENE_ID",
    "SNP_GENE_IDS",
    "UPSTREAM_GENE_DISTANCE",
    "DOWNSTREAM_GENE_DISTANCE",
    "STRONGEST SNP-RISK ALLELE",
    "SNPS",
    "MERGED",
    "SNP_ID_CURRENT",
    "CONTEXT",
    "INTERGENIC",
    "RISK ALLELE FREQUENCY",
    "P-VALUE",
    "PVALUE_MLOG",
    "P-VALUE (TEXT)",
    "OR or BETA",
    "95% CI (TEXT)",
    "PLATFORM [SNPS PASSING QC]",
    "CNV",
    "MAPPED_TRAIT",
    "MAPPED_TRAIT_URI",
    "STUDY ACCESSION",
    "GENOTYPING TECHNOLOGY",
)

_CHROM_RANK = {str(index): index for index in range(1, 23)} | {"X": 23, "Y": 24}
_RISK_ALLELE = re.compile(r"^(rs[0-9]+)-([ACGT]+)$", re.IGNORECASE)
_RSID = re.compile(r"rs[0-9]+", re.IGNORECASE)
_CI_NUMBER = r"[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[Ee][-+]?[0-9]+)?"
_CI_RANGE = re.compile(
    rf"[\[(]\s*({_CI_NUMBER})\s*(?:-|\u2013|\u2014|to)\s*({_CI_NUMBER})\s*[\])]",
    re.IGNORECASE,
)


def _archive_identity(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ResourceError(f"GWAS Catalog archive not found: {path}")
    identity = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    if identity["bytes"] != GWAS_CATALOG_ARCHIVE_BYTES:
        raise ResourceError(
            "GWAS Catalog archive byte count mismatch: "
            f"{identity['bytes']} != {GWAS_CATALOG_ARCHIVE_BYTES}"
        )
    if identity["sha256"] != GWAS_CATALOG_ARCHIVE_SHA256:
        raise ResourceError(
            "GWAS Catalog archive SHA-256 mismatch: "
            f"{identity['sha256']} != {GWAS_CATALOG_ARCHIVE_SHA256}"
        )
    try:
        with zipfile.ZipFile(path) as archive:
            members = archive.infolist()
            if len(members) != 1:
                raise ResourceError(
                    f"GWAS Catalog archive has {len(members)} members, expected one"
                )
            member = members[0]
            member_path = PurePosixPath(member.filename)
            if (
                member_path.is_absolute()
                or ".." in member_path.parts
                or member.is_dir()
                or member.filename != GWAS_CATALOG_MEMBER
            ):
                raise ResourceError(f"unexpected GWAS Catalog ZIP member: {member.filename}")
            if member.file_size != GWAS_CATALOG_MEMBER_BYTES:
                raise ResourceError(
                    "GWAS Catalog member byte count mismatch: "
                    f"{member.file_size} != {GWAS_CATALOG_MEMBER_BYTES}"
                )
            bad_member = archive.testzip()
            if bad_member is not None:
                raise ResourceError(f"GWAS Catalog ZIP CRC failed for {bad_member}")
            identity.update(
                {
                    "member": member.filename,
                    "member_bytes": member.file_size,
                    "member_crc32": f"{member.CRC:08x}",
                }
            )
    except (OSError, zipfile.BadZipFile) as exc:
        raise ResourceError(f"invalid GWAS Catalog ZIP archive: {exc}") from exc
    return identity


def _download(config: BuildConfig) -> None:
    destination = config.gwas.archive
    if destination.is_file():
        return
    if not config.allow_download:
        raise ResourceError(
            "GWAS Catalog archive is absent and resources.allow_download is false; "
            f"place the pinned archive at {destination}"
        )
    curl = shutil.which("curl")
    if curl is None:
        raise ResourceError("curl is required for restart-safe GWAS Catalog download")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f".{destination.name}.partial")
    completed = subprocess.run(
        [
            curl,
            "--fail",
            "--location",
            "--retry",
            "8",
            "--retry-delay",
            "5",
            "--continue-at",
            "-",
            "--output",
            str(partial),
            GWAS_CATALOG_URL,
        ],
        check=False,
    )
    if completed.returncode != 0:
        raise ResourceError(
            "GWAS Catalog download failed; restartable partial retained at " + str(partial)
        )
    _archive_identity(partial)
    os.replace(partial, destination)


def _float(value: str) -> float | None:
    text = value.strip()
    if not text or text.upper() in {"NR", "NA", "N/A", "?"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _effect_type(effect_size: float | None, ci_text: str) -> tuple[str | None, str | None]:
    if effect_size is None:
        return None, None
    lowered = ci_text.lower()
    if re.search(r"\bhazard\s+ratio\b|\bhr\b", lowered):
        return "hazard_ratio", "explicit hazard-ratio marker in 95% CI text"
    if re.search(r"\bodds\s+ratio\b", lowered):
        return "odds_ratio", "explicit odds-ratio marker in 95% CI text"
    if re.search(r"\bincrease(?:d)?\b|\bdecrease(?:d)?\b", lowered):
        return "beta", "increase/decrease marker required by GWAS Catalog beta curation"
    return None, None


def _orient_effect(
    effect_size: float | None,
    effect_type: str | None,
    ci_lower: float | None,
    ci_upper: float | None,
    ci_text: str,
) -> tuple[float | None, float | None, float | None]:
    """Apply the GWAS Catalog's textual beta direction without guessing."""

    if effect_size is None or effect_type != "beta":
        return effect_size, ci_lower, ci_upper
    lowered = ci_text.lower()
    if "decrease" in lowered or "decreased" in lowered:
        if ci_upper is not None and ci_upper <= 0:
            # Some catalog rows already carry signed negative confidence bounds.
            oriented_lower, oriented_upper = ci_lower, ci_upper
        else:
            # Other rows carry magnitude-only or unsigned/cross-zero bounds.
            oriented_lower = -ci_upper if ci_upper is not None else None
            oriented_upper = -ci_lower if ci_lower is not None else None
        return -abs(effect_size), oriented_lower, oriented_upper
    if "increase" in lowered or "increased" in lowered:
        return abs(effect_size), ci_lower, ci_upper
    raise ResourceError("typed beta row lacks an explicit increase/decrease direction")


def _ci(text: str) -> tuple[float | None, float | None]:
    match = _CI_RANGE.search(text)
    if match is None:
        return None, None
    lower = _float(match.group(1))
    upper = _float(match.group(2))
    if lower is None or upper is None or lower > upper:
        return None, None
    return lower, upper


def _sample_size(initial: str, replication: str) -> int | None:
    # Free text frequently contains stage, case/control, range, and ancestry
    # counts.  Only a single unambiguous integer is retained; anything more
    # complex remains null rather than silently inventing a cohort size.
    tokens = re.findall(r"(?<![0-9.])([0-9][0-9,]*)(?![0-9.])", f"{initial} {replication}")
    if len(tokens) != 1:
        return None
    value = int(tokens[0].replace(",", ""))
    return value if 0 < value <= 9_223_372_036_854_775_807 else None


def _row_digest(row: dict[str, str]) -> str:
    payload = {key: row.get(key, "") for key in EXPECTED_HEADER}
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _normalized_row(
    row: dict[str, str], source_row_number: int
) -> tuple[dict[str, Any] | None, str | None]:
    chrom_raw = row["CHR_ID"].strip().upper().removeprefix("CHR")
    if chrom_raw not in _CHROM_RANK:
        return None, "non_single_core_chromosome"
    position_raw = row["CHR_POS"].strip()
    if not re.fullmatch(r"[0-9]+", position_raw):
        return None, "non_single_position"
    position = int(position_raw)
    if not 1 <= position <= 2_147_483_647:
        return None, "invalid_position"

    risk = _RISK_ALLELE.fullmatch(row["STRONGEST SNP-RISK ALLELE"].strip())
    if risk is None:
        return None, "unknown_or_complex_risk_allele"
    risk_rsid = risk.group(1).lower()
    effect_allele = risk.group(2).upper()
    aliases = {value.lower() for value in _RSID.findall(row["SNPS"])}
    aliases.add(risk_rsid)
    current = row["SNP_ID_CURRENT"].strip()
    if current.isdigit():
        aliases.add(f"rs{int(current)}")

    reported_trait = row["DISEASE/TRAIT"].strip()
    pmid = row["PUBMEDID"].strip()
    study_accession = row["STUDY ACCESSION"].strip()
    if not reported_trait:
        return None, "missing_reported_trait"
    if not re.fullmatch(r"[0-9]+", pmid):
        return None, "missing_or_invalid_pmid"
    if not re.fullmatch(r"GCST[0-9]+", study_accession):
        return None, "missing_or_invalid_study_accession"
    p_value = _float(row["P-VALUE"])
    if p_value is None or not 0.0 < p_value <= 1.0:
        return None, "missing_or_invalid_p_value"

    raw_effect = row["OR or BETA"].strip() or None
    raw_ci = row["95% CI (TEXT)"].strip() or None
    effect_size = _float(raw_effect or "")
    effect_type, effect_basis = _effect_type(effect_size, raw_ci or "")
    ci_lower, ci_upper = _ci(raw_ci or "")
    effect_size, ci_lower, ci_upper = _orient_effect(
        effect_size,
        effect_type,
        ci_lower,
        ci_upper,
        raw_ci or "",
    )
    mapped_trait = row["MAPPED_TRAIT"].strip()
    mapped_uris = sorted(
        {value.strip() for value in re.split(r"[,;]", row["MAPPED_TRAIT_URI"]) if value.strip()}
    )
    return (
        {
            "chrom": f"chr{chrom_raw}",
            "pos": position,
            "source_row_number": source_row_number,
            "risk_rsid": risk_rsid,
            "rsid_aliases": sorted(aliases),
            "effect_allele": effect_allele,
            "reported_trait": reported_trait,
            "mapped_traits": [mapped_trait] if mapped_trait else [],
            "mapped_trait_uris": mapped_uris,
            "p_value": p_value,
            "effect_size": effect_size,
            "effect_type": effect_type,
            "effect_type_basis": effect_basis,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "sample_size": _sample_size(row["INITIAL SAMPLE SIZE"], row["REPLICATION SAMPLE SIZE"]),
            "pmid": pmid,
            "study_accession": study_accession,
            "raw_effect_text": raw_effect,
            "raw_ci_text": raw_ci,
            "source_record_digest": _row_digest(row),
        },
        None,
    )


_SQL_COLUMNS = (
    "chrom",
    "pos",
    "source_row_number",
    "risk_rsid",
    "rsid_aliases",
    "effect_allele",
    "reported_trait",
    "mapped_traits",
    "mapped_trait_uris",
    "p_value",
    "effect_size",
    "effect_type",
    "effect_type_basis",
    "ci_lower",
    "ci_upper",
    "sample_size",
    "pmid",
    "study_accession",
    "raw_effect_text",
    "raw_ci_text",
    "source_record_digest",
)


def _sql_values(row: dict[str, Any]) -> tuple[Any, ...]:
    values = dict(row)
    for key in ("rsid_aliases", "mapped_traits", "mapped_trait_uris"):
        values[key] = json.dumps(values[key], ensure_ascii=False, separators=(",", ":"))
    return tuple(values[key] for key in _SQL_COLUMNS)


def _arrow_row(values: tuple[Any, ...]) -> dict[str, Any]:
    row = dict(zip(_SQL_COLUMNS, values, strict=True))
    for key in ("rsid_aliases", "mapped_traits", "mapped_trait_uris"):
        row[key] = json.loads(row[key])
    return row


def _normalize_archive(archive_path: Path, destination: Path) -> dict[str, Any]:
    counters: Counter[str] = Counter()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".gwas-normalize-", dir=destination.parent) as temp:
        temp_root = Path(temp)
        database_path = temp_root / "rows.sqlite3"
        output_path = temp_root / "gwas_catalog.parquet"
        connection = sqlite3.connect(database_path)
        try:
            connection.execute("PRAGMA journal_mode=OFF")
            connection.execute("PRAGMA synchronous=OFF")
            connection.execute("PRAGMA temp_store=FILE")
            connection.execute(
                """
                CREATE TABLE resource (
                    chrom TEXT NOT NULL, pos INTEGER NOT NULL,
                    source_row_number INTEGER NOT NULL, risk_rsid TEXT NOT NULL,
                    rsid_aliases TEXT NOT NULL, effect_allele TEXT NOT NULL,
                    reported_trait TEXT NOT NULL, mapped_traits TEXT NOT NULL,
                    mapped_trait_uris TEXT NOT NULL, p_value REAL NOT NULL,
                    effect_size REAL, effect_type TEXT, effect_type_basis TEXT,
                    ci_lower REAL, ci_upper REAL, sample_size INTEGER,
                    pmid TEXT NOT NULL, study_accession TEXT NOT NULL,
                    raw_effect_text TEXT, raw_ci_text TEXT,
                    source_record_digest TEXT NOT NULL
                )
                """
            )
            placeholders = ",".join("?" for _ in _SQL_COLUMNS)
            insert = f"INSERT INTO resource VALUES ({placeholders})"
            pending: list[tuple[Any, ...]] = []
            with (
                zipfile.ZipFile(archive_path) as archive,
                archive.open(GWAS_CATALOG_MEMBER) as binary,
            ):
                text = io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
                reader = csv.DictReader(text, delimiter="\t")
                if tuple(reader.fieldnames or ()) != EXPECTED_HEADER:
                    raise ResourceError("GWAS Catalog TSV header does not match pinned schema")
                for source_row_number, raw in enumerate(reader, start=2):
                    counters["source_rows"] += 1
                    if None in raw:
                        raise ResourceError(
                            f"GWAS Catalog row {source_row_number} has extra columns"
                        )
                    normalized, reason = _normalized_row(raw, source_row_number)
                    if normalized is None:
                        assert reason is not None
                        counters[f"omitted_{reason}"] += 1
                        continue
                    pending.append(_sql_values(normalized))
                    counters["normalized_rows"] += 1
                    counters[
                        "typed_effect_rows"
                        if normalized["effect_type"] is not None
                        else "untyped_effect_rows"
                    ] += 1
                    if len(pending) >= 10_000:
                        connection.executemany(insert, pending)
                        pending.clear()
                if pending:
                    connection.executemany(insert, pending)
            connection.commit()

            order = (
                "CASE chrom "
                + " ".join(f"WHEN 'chr{chrom}' THEN {rank}" for chrom, rank in _CHROM_RANK.items())
                + " END, pos, risk_rsid, effect_allele, source_row_number"
            )
            cursor = connection.execute(
                f"SELECT {','.join(_SQL_COLUMNS)} FROM resource ORDER BY {order}"
            )
            schema = gwas_catalog_resource_schema()
            with pq.ParquetWriter(
                output_path,
                schema,
                compression="zstd",
                use_dictionary=True,
                write_statistics=True,
            ) as writer:
                while values := cursor.fetchmany(50_000):
                    rows = [_arrow_row(tuple(value)) for value in values]
                    writer.write_table(pa.Table.from_pylist(rows, schema=schema))
        finally:
            connection.close()
        os.replace(output_path, destination)
    return dict(sorted(counters.items()))


def sync_gwas_catalog(config: BuildConfig) -> dict[str, Any]:
    if not config.gwas.enabled:
        return {"status": "not_enabled"}
    _download(config)
    archive = _archive_identity(config.gwas.archive)
    counters = _normalize_archive(config.gwas.archive, config.gwas.normalized_parquet)
    manifest = {
        "resource_type": "nhgri_ebi_gwas_catalog_literature_associations",
        "snapshot": GWAS_CATALOG_SNAPSHOT,
        "source_url": GWAS_CATALOG_URL,
        "license": config.gwas.license,
        "attribution": GWAS_CATALOG_ATTRIBUTION,
        "source_archive_redistributed": False,
        "derived_rows_redistribution_allowed": config.gwas.redistribution_allowed,
        "archive_sha256": archive["sha256"],
        "archive_bytes": archive["bytes"],
        "archive_member": archive["member"],
        "archive_member_bytes": archive["member_bytes"],
        "archive_member_crc32": archive["member_crc32"],
        "normalized_parquet_sha256": sha256_file(config.gwas.normalized_parquet),
        "normalized_parquet_bytes": config.gwas.normalized_parquet.stat().st_size,
        "normalized_schema": str(gwas_catalog_resource_schema()),
        "counts": counters,
        "matching_policy": (
            "GRCh38 coordinate + exact rsID alias + exact effect allele equal to REF or ALT"
        ),
        "effect_type_policy": (
            "retain beta/odds_ratio/hazard_ratio only when the legacy bulk row contains "
            "an explicit textual type marker; ambiguous OR-or-BETA rows remain untyped"
        ),
        "generated_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
    }
    config.gwas.manifest.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(config.gwas.manifest, manifest)
    return verify_gwas_catalog(config)


def verify_gwas_catalog(config: BuildConfig) -> dict[str, Any]:
    if not config.gwas.enabled:
        return {"status": "not_enabled"}
    if config.gwas.snapshot != GWAS_CATALOG_SNAPSHOT:
        raise ResourceError("configured GWAS Catalog snapshot is not pinned")
    if config.gwas.archive_expected_sha256 != GWAS_CATALOG_ARCHIVE_SHA256:
        raise ResourceError("configured GWAS Catalog SHA-256 is not pinned")
    archive = _archive_identity(config.gwas.archive)
    try:
        manifest = json.loads(config.gwas.manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ResourceError(f"invalid GWAS Catalog resource manifest: {exc}") from exc
    required = {
        "resource_type",
        "snapshot",
        "source_url",
        "license",
        "attribution",
        "archive_sha256",
        "archive_bytes",
        "archive_member",
        "archive_member_bytes",
        "archive_member_crc32",
        "normalized_parquet_sha256",
        "normalized_parquet_bytes",
        "counts",
        "matching_policy",
        "effect_type_policy",
    }
    if not isinstance(manifest, dict) or not required.issubset(manifest):
        raise ResourceError("GWAS Catalog resource manifest is missing required fields")
    expected = {
        "resource_type": "nhgri_ebi_gwas_catalog_literature_associations",
        "snapshot": GWAS_CATALOG_SNAPSHOT,
        "source_url": GWAS_CATALOG_URL,
        "license": config.gwas.license,
        "attribution": GWAS_CATALOG_ATTRIBUTION,
        "archive_sha256": archive["sha256"],
        "archive_bytes": archive["bytes"],
        "archive_member": archive["member"],
        "archive_member_bytes": archive["member_bytes"],
        "archive_member_crc32": archive["member_crc32"],
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ResourceError(
                f"GWAS Catalog manifest {key} mismatch: {manifest.get(key)!r} != {value!r}"
            )
    if not config.gwas.normalized_parquet.is_file():
        raise ResourceError(
            f"normalized GWAS Catalog Parquet not found: {config.gwas.normalized_parquet}"
        )
    normalized_sha = sha256_file(config.gwas.normalized_parquet)
    if normalized_sha != manifest["normalized_parquet_sha256"]:
        raise ResourceError("normalized GWAS Catalog Parquet SHA-256 mismatch")
    if config.gwas.normalized_parquet.stat().st_size != manifest["normalized_parquet_bytes"]:
        raise ResourceError("normalized GWAS Catalog Parquet byte count mismatch")
    actual_schema = pq.read_schema(config.gwas.normalized_parquet)
    # Parquet canonicalises list child names from ``item`` to ``element``;
    # logical/physical types and nullability must match, field metadata is not
    # used by this resource schema.
    if not actual_schema.equals(gwas_catalog_resource_schema(), check_metadata=False):
        raise ResourceError("normalized GWAS Catalog Parquet physical schema mismatch")
    normalized_rows = int(pq.ParquetFile(config.gwas.normalized_parquet).metadata.num_rows)
    if normalized_rows != int(manifest["counts"].get("normalized_rows", -1)):
        raise ResourceError("normalized GWAS Catalog Parquet row count mismatch")
    return {
        **manifest,
        "status": "verified",
        "manifest_sha256": sha256_file(config.gwas.manifest),
    }


__all__ = [
    "EXPECTED_HEADER",
    "sync_gwas_catalog",
    "verify_gwas_catalog",
]
