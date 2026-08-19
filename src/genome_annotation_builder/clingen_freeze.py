"""Deterministically freeze public ClinGen gene-level curation snapshots.

The frozen document is intentionally gene-level evidence.  It does not claim
that any variant in a curated gene causes the associated condition, and it is
not a substitute for exact-allele variant interpretation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import date, datetime
from pathlib import Path
from typing import Any

from .errors import ResourceError
from .util import canonical_json_bytes, sha256_file

CLINGEN_GENE_DISEASE_HEADER = (
    "GENE SYMBOL",
    "GENE ID (HGNC)",
    "DISEASE LABEL",
    "DISEASE ID (MONDO)",
    "MOI",
    "SOP",
    "CLASSIFICATION",
    "ONLINE REPORT",
    "CLASSIFICATION DATE",
    "GCEP",
)
CLINGEN_DOSAGE_HEADER = (
    "Gene Symbol",
    "Gene ID",
    "cytoBand",
    "Genomic Location",
    "Haploinsufficiency Score",
    "Haploinsufficiency Description",
    "Haploinsufficiency PMID1",
    "Haploinsufficiency PMID2",
    "Haploinsufficiency PMID3",
    "Haploinsufficiency PMID4",
    "Haploinsufficiency PMID5",
    "Haploinsufficiency PMID6",
    "Triplosensitivity Score",
    "Triplosensitivity Description",
    "Triplosensitivity PMID1",
    "Triplosensitivity PMID2",
    "Triplosensitivity PMID3",
    "Triplosensitivity PMID4",
    "Triplosensitivity PMID5",
    "Triplosensitivity PMID6",
    "Date Last Evaluated",
    "Haploinsufficiency Disease ID",
    "Triplosensitivity Disease ID",
)

GENE_DISEASE_CLASSIFICATIONS = frozenset(
    {
        "Definitive",
        "Strong",
        "Moderate",
        "Limited",
        "Disputed",
        "Refuted",
        "No Known Disease Relationship",
    }
)
DOSAGE_SCORE_DESCRIPTIONS = {
    0: "No evidence available",
    1: "Little evidence for dosage pathogenicity",
    2: "Some evidence for dosage pathogenicity",
    3: "Sufficient evidence for dosage pathogenicity",
    30: "Gene associated with autosomal recessive phenotype",
    40: "Dosage sensitivity unlikely",
}

_HGNC = re.compile(r"HGNC:[1-9][0-9]*")
_MONDO = re.compile(r"MONDO:[0-9]+")
_REPORT = re.compile(
    r"https://search\.clinicalgenome\.org/kb/gene-validity/"
    r"((?:CGGV|CGGCIEX):assertion_[A-Za-z0-9_.:-]+)"
)
_LOCATION = re.compile(r"chr([0-9]+|X|Y|M):([0-9]+)-([0-9]+)")
_PMID = re.compile(r"[0-9]+")


def _row_digest(values: list[str]) -> str:
    return hashlib.sha256(canonical_json_bytes(values)).hexdigest()


def _parse_iso_datetime(value: str, field: str) -> str:
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ResourceError(f"invalid {field}: {value!r}") from exc
    if parsed.tzinfo is None:
        # A few historical ClinGen rows omit the suffix.  Preserve that exact
        # source representation rather than inventing a timezone.
        return value
    return value


def _parse_gene_disease(path: Path, snapshot: date) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            source_rows = list(csv.reader(handle))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ResourceError(f"cannot read ClinGen gene-disease CSV: {exc}") from exc
    if len(source_rows) < 7 or any(
        len(row) != len(CLINGEN_GENE_DISEASE_HEADER) for row in source_rows
    ):
        raise ResourceError("ClinGen gene-disease CSV has an invalid row width")
    if source_rows[0][0] != "CLINGEN GENE DISEASE VALIDITY CURATIONS":
        raise ResourceError("unexpected ClinGen gene-disease title")
    expected_created = f"FILE CREATED: {snapshot.isoformat()}"
    if source_rows[1][0] != expected_created:
        raise ResourceError(f"ClinGen gene-disease creation date mismatch: {source_rows[1][0]!r}")
    if tuple(source_rows[4]) != CLINGEN_GENE_DISEASE_HEADER:
        raise ResourceError("ClinGen gene-disease CSV header changed")
    if not all(
        set(row) == {"+" * len(value) for value in row} for row in (source_rows[3], source_rows[5])
    ):
        raise ResourceError("ClinGen gene-disease separator rows changed")

    output: list[dict[str, Any]] = []
    assertion_ids: set[str] = set()
    source_digests: set[str] = set()
    for source_row_number, values in enumerate(source_rows[6:], start=7):
        row = dict(zip(CLINGEN_GENE_DISEASE_HEADER, values, strict=True))
        required = tuple(row[field].strip() for field in CLINGEN_GENE_DISEASE_HEADER)
        if any(not value for value in required):
            raise ResourceError(
                f"ClinGen gene-disease row {source_row_number} has an empty required field"
            )
        if _HGNC.fullmatch(row["GENE ID (HGNC)"]) is None:
            raise ResourceError(f"invalid HGNC ID at ClinGen row {source_row_number}")
        if _MONDO.fullmatch(row["DISEASE ID (MONDO)"]) is None:
            raise ResourceError(f"invalid MONDO ID at ClinGen row {source_row_number}")
        if row["CLASSIFICATION"] not in GENE_DISEASE_CLASSIFICATIONS:
            raise ResourceError(
                f"unknown ClinGen classification at row {source_row_number}: "
                f"{row['CLASSIFICATION']!r}"
            )
        report_match = _REPORT.fullmatch(row["ONLINE REPORT"])
        if report_match is None:
            raise ResourceError(f"invalid ClinGen report URL at row {source_row_number}")
        assertion_id = report_match.group(1)
        if assertion_id in assertion_ids:
            raise ResourceError(f"duplicate ClinGen assertion ID: {assertion_id}")
        assertion_ids.add(assertion_id)
        _parse_iso_datetime(row["CLASSIFICATION DATE"], "ClinGen classification date")
        source_digest = _row_digest(values)
        if source_digest in source_digests:
            raise ResourceError(f"duplicate ClinGen source row digest: {source_digest}")
        source_digests.add(source_digest)
        output.append(
            {
                "assertion_id": assertion_id,
                "gene_symbol": row["GENE SYMBOL"],
                "hgnc_id": row["GENE ID (HGNC)"],
                "disease_label": row["DISEASE LABEL"],
                "mondo_id": row["DISEASE ID (MONDO)"],
                "mode_of_inheritance": row["MOI"],
                "sop": row["SOP"],
                "classification": row["CLASSIFICATION"],
                "online_report": row["ONLINE REPORT"],
                "classification_date": row["CLASSIFICATION DATE"],
                "gcep": row["GCEP"],
                "source_row_number": source_row_number,
                "source_row_sha256": source_digest,
            }
        )
    output.sort(key=lambda row: (row["hgnc_id"], row["mondo_id"], row["assertion_id"]))
    return output, {
        "row_count": len(output),
        "unique_assertions": len(assertion_ids),
        "unique_genes": len({row["hgnc_id"] for row in output}),
        "classifications": sorted({row["classification"] for row in output}),
    }


def _dosage_score(value: str, description: str, source_row_number: int) -> int | None:
    if value == "Not yet evaluated":
        if description != "Not yet evaluated":
            raise ResourceError(f"inconsistent unevaluated dosage row {source_row_number}")
        return None
    try:
        score = int(value)
    except ValueError as exc:
        raise ResourceError(f"invalid dosage score at row {source_row_number}: {value!r}") from exc
    if DOSAGE_SCORE_DESCRIPTIONS.get(score) != description:
        raise ResourceError(
            f"dosage score/description mismatch at row {source_row_number}: "
            f"{value!r}/{description!r}"
        )
    return score


def _pmids(row: dict[str, str], prefix: str) -> tuple[list[int], list[str]]:
    raw = [row[f"{prefix} PMID{index}"] for index in range(1, 7)]
    source_values = [value for value in raw if value]
    parsed: list[int] = []
    seen: set[int] = set()
    for value in source_values:
        for token in _PMID.findall(value):
            pmid = int(token)
            if pmid > 0 and pmid not in seen:
                seen.add(pmid)
                parsed.append(pmid)
    return parsed, source_values


def _parse_dosage(path: Path, snapshot: date) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ResourceError(f"cannot read ClinGen dosage TSV: {exc}") from exc
    comments = [line for line in lines if line.startswith("#")]
    data_lines = [line for line in lines if not line.startswith("#")]
    if len(comments) != 6 or not data_lines:
        raise ResourceError("ClinGen dosage TSV preamble changed")
    try:
        data_date = datetime.strptime(comments[1], "#%d %b,%Y").date()
    except ValueError as exc:
        raise ResourceError(f"invalid ClinGen dosage data date: {comments[1]!r}") from exc
    if data_date != snapshot:
        raise ResourceError(
            f"ClinGen dosage data date mismatch: {data_date.isoformat()} != {snapshot.isoformat()}"
        )
    if "GRCh38 (hg38): GCF_000001405.36" not in comments[2]:
        raise ResourceError("ClinGen dosage assembly declaration changed")
    header = tuple(comments[-1].removeprefix("#").split("\t"))
    if header != CLINGEN_DOSAGE_HEADER:
        raise ResourceError("ClinGen dosage TSV header changed")
    try:
        reader = csv.DictReader(data_lines, fieldnames=header, delimiter="\t", strict=True)
        source_rows = list(reader)
    except csv.Error as exc:
        raise ResourceError(f"invalid ClinGen dosage TSV: {exc}") from exc
    if any(None in row or len(row) != len(header) for row in source_rows):
        raise ResourceError("ClinGen dosage TSV has an invalid row width")

    output: list[dict[str, Any]] = []
    symbols: set[str] = set()
    gene_ids: set[int] = set()
    for source_row_number, row in enumerate(source_rows, start=7):
        symbol = row["Gene Symbol"]
        if not symbol or symbol in symbols:
            raise ResourceError(
                f"missing or duplicate ClinGen dosage symbol at row {source_row_number}"
            )
        symbols.add(symbol)
        try:
            gene_id = int(row["Gene ID"])
        except ValueError as exc:
            raise ResourceError(f"invalid NCBI Gene ID at dosage row {source_row_number}") from exc
        if gene_id <= 0 or gene_id in gene_ids:
            raise ResourceError(f"duplicate/invalid NCBI Gene ID at dosage row {source_row_number}")
        gene_ids.add(gene_id)
        location = _LOCATION.fullmatch(row["Genomic Location"])
        if location is None:
            raise ResourceError(f"invalid GRCh38 location at dosage row {source_row_number}")
        start, end = int(location.group(2)), int(location.group(3))
        if start < 1 or end < start:
            raise ResourceError(f"invalid GRCh38 interval at dosage row {source_row_number}")
        hi_score = _dosage_score(
            row["Haploinsufficiency Score"],
            row["Haploinsufficiency Description"],
            source_row_number,
        )
        ts_score = _dosage_score(
            row["Triplosensitivity Score"],
            row["Triplosensitivity Description"],
            source_row_number,
        )
        hi_pmids, hi_pmid_source_values = _pmids(row, "Haploinsufficiency")
        ts_pmids, ts_pmid_source_values = _pmids(row, "Triplosensitivity")
        last_evaluated = row["Date Last Evaluated"]
        try:
            date.fromisoformat(last_evaluated)
        except ValueError as exc:
            raise ResourceError(
                f"invalid dosage evaluation date at row {source_row_number}: {last_evaluated!r}"
            ) from exc
        for disease_field in (
            "Haploinsufficiency Disease ID",
            "Triplosensitivity Disease ID",
        ):
            disease_id = row[disease_field]
            if disease_id and _MONDO.fullmatch(disease_id) is None:
                raise ResourceError(f"invalid dosage MONDO ID at row {source_row_number}")
        raw_values = [row[field] for field in header]
        output.append(
            {
                "gene_symbol": symbol,
                "ncbi_gene_id": gene_id,
                "cytoband": row["cytoBand"],
                "chrom": location.group(1),
                "start": start,
                "end": end,
                "haploinsufficiency_score": hi_score,
                "haploinsufficiency_source_score": row["Haploinsufficiency Score"],
                "haploinsufficiency_description": row["Haploinsufficiency Description"],
                "haploinsufficiency_pmids": hi_pmids,
                "haploinsufficiency_pmid_source_values": hi_pmid_source_values,
                "haploinsufficiency_disease_id": row["Haploinsufficiency Disease ID"] or None,
                "triplosensitivity_score": ts_score,
                "triplosensitivity_source_score": row["Triplosensitivity Score"],
                "triplosensitivity_description": row["Triplosensitivity Description"],
                "triplosensitivity_pmids": ts_pmids,
                "triplosensitivity_pmid_source_values": ts_pmid_source_values,
                "triplosensitivity_disease_id": row["Triplosensitivity Disease ID"] or None,
                "last_evaluated": last_evaluated,
                "source_row_number": source_row_number,
                "source_row_sha256": _row_digest(raw_values),
            }
        )
    output.sort(key=lambda row: (row["gene_symbol"], row["ncbi_gene_id"]))
    return output, {
        "row_count": len(output),
        "unique_genes": len(symbols),
        "assembly": "GRCh38",
        "assembly_accession": "GCF_000001405.36",
        "coordinate_system": "1-based inclusive source coordinates",
        "score_semantics": {
            str(score): description for score, description in DOSAGE_SCORE_DESCRIPTIONS.items()
        }
        | {"null": "Not yet evaluated"},
    }


def freeze_clingen_snapshot(
    gene_disease_csv: str | Path,
    dosage_tsv: str | Path,
    *,
    gene_disease_snapshot: date,
    dosage_snapshot: date,
    accessed_at: date,
    gene_disease_url: str,
    dosage_url: str,
) -> dict[str, Any]:
    """Return a deterministic, lossless-enough canonical ClinGen snapshot."""

    gene_disease_csv = Path(gene_disease_csv)
    dosage_tsv = Path(dosage_tsv)
    gene_disease, gene_disease_stats = _parse_gene_disease(gene_disease_csv, gene_disease_snapshot)
    dosage, dosage_stats = _parse_dosage(dosage_tsv, dosage_snapshot)
    return {
        "schema_version": 1,
        "resource_type": "clingen_gene_level_curations",
        "snapshot": accessed_at.isoformat(),
        "license": "CC0-1.0",
        "attribution": (
            "Curated content obtained from the Clinical Genome Resource (ClinGen); "
            f"accessed {accessed_at.isoformat()}."
        ),
        "semantic_boundary": (
            "Gene-level validity and dosage context only; presence of a variant in a curated "
            "gene is not evidence that the variant is pathogenic or that the subject has the "
            "associated condition."
        ),
        "sources": {
            "gene_disease_validity": {
                "url": gene_disease_url,
                "snapshot": gene_disease_snapshot.isoformat(),
                "bytes": gene_disease_csv.stat().st_size,
                "sha256": sha256_file(gene_disease_csv),
                **gene_disease_stats,
            },
            "dosage_sensitivity": {
                "url": dosage_url,
                "snapshot": dosage_snapshot.isoformat(),
                "bytes": dosage_tsv.stat().st_size,
                "sha256": sha256_file(dosage_tsv),
                **dosage_stats,
            },
        },
        "gene_disease_validity": gene_disease,
        "dosage_sensitivity": dosage,
    }


def snapshot_identity(document: dict[str, Any]) -> dict[str, Any]:
    """Return stable high-level identity after strict structural checks."""

    if document.get("schema_version") != 1:
        raise ResourceError("unsupported ClinGen canonical snapshot schema")
    if document.get("resource_type") != "clingen_gene_level_curations":
        raise ResourceError("unexpected ClinGen canonical resource type")
    if document.get("license") != "CC0-1.0":
        raise ResourceError("ClinGen canonical snapshot lacks the CC0 boundary")
    gene_rows = document.get("gene_disease_validity")
    dosage_rows = document.get("dosage_sensitivity")
    sources = document.get("sources")
    if not isinstance(gene_rows, list) or not isinstance(dosage_rows, list):
        raise ResourceError("ClinGen canonical snapshot rowsets are invalid")
    if not isinstance(sources, dict):
        raise ResourceError("ClinGen canonical snapshot sources are invalid")
    if sources.get("gene_disease_validity", {}).get("row_count") != len(gene_rows):
        raise ResourceError("ClinGen gene-disease row count mismatch")
    if sources.get("dosage_sensitivity", {}).get("row_count") != len(dosage_rows):
        raise ResourceError("ClinGen dosage row count mismatch")
    assertion_ids = [row.get("assertion_id") for row in gene_rows if isinstance(row, dict)]
    dosage_symbols = [row.get("gene_symbol") for row in dosage_rows if isinstance(row, dict)]
    if len(assertion_ids) != len(gene_rows) or len(set(assertion_ids)) != len(assertion_ids):
        raise ResourceError("ClinGen canonical gene-disease assertions are not unique")
    if len(dosage_symbols) != len(dosage_rows) or len(set(dosage_symbols)) != len(dosage_symbols):
        raise ResourceError("ClinGen canonical dosage genes are not unique")
    return {
        "snapshot": document.get("snapshot"),
        "gene_disease_rows": len(gene_rows),
        "dosage_rows": len(dosage_rows),
        "canonical_sha256": hashlib.sha256(canonical_json_bytes(document)).hexdigest(),
    }


def load_canonical_snapshot(path: str | Path) -> dict[str, Any]:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ResourceError(f"cannot read canonical ClinGen snapshot: {exc}") from exc
    if not isinstance(document, dict):
        raise ResourceError("canonical ClinGen snapshot must be a JSON object")
    snapshot_identity(document)
    return document
