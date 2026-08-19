"""Streaming exact-allele ClinVar matching and lossless clinical evidence."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pysam

from .constants import STANDARD_CHROMS
from .errors import InputValidationError, ResourceError
from .vcf_ingest import canonical_core_chrom

CHROM_RANK = {chrom: index for index, chrom in enumerate(STANDARD_CHROMS)}
CORE_SIGNIFICANCE = {
    "Pathogenic",
    "Likely_pathogenic",
    "Uncertain_significance",
    "Likely_benign",
    "Benign",
    "Conflicting_interpretations",
}


def _raw_values(record: pysam.VariantRecord, key: str) -> list[str]:
    value = record.info.get(key)
    if value is None:
        return []
    values = value if isinstance(value, tuple) else (value,)
    return [str(item) for item in values if item not in {None, "", "."}]


def _split_values(record: pysam.VariantRecord, key: str) -> list[str]:
    return sorted(
        {
            item
            for value in _raw_values(record, key)
            for item in value.split("|")
            if item and item != "."
        }
    )


def _display(value: str) -> str:
    return re.sub(r"_+", " ", value).strip()


def _review_tokens(record: pysam.VariantRecord) -> list[str]:
    return sorted({_display(value).lower() for value in _raw_values(record, "CLNREVSTAT")})


def _review_stars(tokens: list[str]) -> int:
    values = set(tokens)
    if "practice guideline" in values:
        return 4
    if "reviewed by expert panel" in values:
        return 3
    if {"criteria provided", "multiple submitters", "no conflicts"}.issubset(values):
        return 2
    if "criteria provided" in values:
        return 1
    return 0


def _core_significance(raw_values: list[str], conflicts: list[str]) -> tuple[str | None, bool]:
    lowered = {_display(value).lower().replace(" ", "_") for value in raw_values}
    has_conflicts = bool(conflicts) or any("conflict" in value for value in lowered)
    pathogenic = any("pathogenic" in value and "conflict" not in value for value in lowered)
    benign = any("benign" in value and "conflict" not in value for value in lowered)
    if pathogenic and benign:
        has_conflicts = True
    if has_conflicts:
        return "Conflicting_interpretations", True
    if lowered == {"pathogenic"}:
        return "Pathogenic", False
    if any(value in {"likely_pathogenic", "pathogenic/likely_pathogenic"} for value in lowered):
        return "Likely_pathogenic", False
    if "pathogenic" in lowered:
        return "Pathogenic", False
    if "uncertain_significance" in lowered:
        return "Uncertain_significance", False
    if any(value in {"likely_benign", "benign/likely_benign"} for value in lowered):
        return "Likely_benign", False
    if "benign" in lowered:
        return "Benign", False
    return None, False


def _record_key(record: pysam.VariantRecord) -> tuple[int, int, str, str] | None:
    chrom = canonical_core_chrom(record.contig)
    if chrom is None or not record.alts or len(record.alts) != 1:
        return None
    ref = record.ref.upper()
    alt = record.alts[0].upper()
    if not re.fullmatch(r"[ACGT]+", ref) or not re.fullmatch(r"[ACGT]+", alt):
        return None
    return CHROM_RANK[chrom], record.pos, ref, alt


def _variation_accession(identifier: str) -> str:
    return f"VCV{int(identifier):09d}" if identifier.isdigit() else identifier


def _record_digest(record: pysam.VariantRecord) -> str:
    payload = {
        "chrom": record.contig,
        "pos": record.pos,
        "id": record.id,
        "ref": record.ref,
        "alts": list(record.alts or ()),
        "info": {key: record.info[key] for key in sorted(record.info)},
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=list)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _integer_values(record: pysam.VariantRecord, key: str) -> list[int]:
    output: set[int] = set()
    for value in _split_values(record, key):
        for token in value.split(","):
            try:
                output.add(int(token))
            except ValueError:
                continue
    return sorted(output)


def _evidence_row(record: pysam.VariantRecord, variant_id: str, snapshot: str) -> dict[str, Any]:
    identifier = str(record.id or "")
    if not identifier:
        raise InputValidationError(
            f"ClinVar record without Variation ID at {record.contig}:{record.pos}"
        )
    raw_significance = [_display(value) for value in _split_values(record, "CLNSIG")]
    conflict_values = [_display(value) for value in _split_values(record, "CLNSIGCONF")]
    aggregate, has_conflicts = _core_significance(raw_significance, conflict_values)
    review = _review_tokens(record)
    digest = _record_digest(record)
    rsids = [
        value if value.lower().startswith("rs") else f"rs{value}"
        for value in _split_values(record, "RS")
    ]
    allele_id = record.info.get("ALLELEID")
    if isinstance(allele_id, tuple):
        allele_id = allele_id[0] if len(allele_id) == 1 else None
    return {
        "evidence_id": hashlib.sha256(f"{snapshot}\0{digest}".encode()).hexdigest(),
        "variant_id": variant_id,
        "source": "NCBI ClinVar",
        "snapshot": snapshot,
        "clinvar_variation_id": identifier,
        "clinvar_accession": _variation_accession(identifier),
        "allele_id": int(allele_id) if allele_id is not None else None,
        "raw_significance": raw_significance,
        "aggregate_significance": aggregate,
        "review_status": review,
        "review_stars": _review_stars(review),
        "has_conflicts": has_conflicts,
        "conflict_classifications": conflict_values,
        "condition_names": [_display(value) for value in _split_values(record, "CLNDN")],
        "condition_identifiers": _split_values(record, "CLNDISDB"),
        "submitter_accessions": _split_values(record, "CLNSIGSCV"),
        "gene_info": _split_values(record, "GENEINFO"),
        "rsids": rsids,
        "origin_codes": _integer_values(record, "ORIGIN"),
        "hgvs": _split_values(record, "CLNHGVS"),
        "variation_type": next(iter(_raw_values(record, "CLNVC")), None),
        "source_record_digest": digest,
    }


def aggregate_clinvar(
    records: list[pysam.VariantRecord], variant_id: str, snapshot: str
) -> tuple[dict[str, Any], list[dict[str, Any]], str | None]:
    rows = [_evidence_row(record, variant_id, snapshot) for record in records]
    raw_significance = sorted({value for row in rows for value in row["raw_significance"]})
    conflict_values = sorted({value for row in rows for value in row["conflict_classifications"]})
    significance, aggregate_conflict = _core_significance(raw_significance, conflict_values)
    has_conflicts = aggregate_conflict or any(row["has_conflicts"] for row in rows)
    if has_conflicts:
        significance = "Conflicting_interpretations"
    if significance is not None and significance not in CORE_SIGNIFICANCE:
        raise InputValidationError(f"unsupported core ClinVar significance: {significance}")
    accessions = sorted({row["clinvar_accession"] for row in rows})
    submitters = {value for row in rows for value in row["submitter_accessions"]}
    rsids = sorted({value for row in rows for value in row["rsids"]})
    summary_values = sorted(
        {
            *conflict_values,
            *(value for row in rows for value in row["raw_significance"] if row["has_conflicts"]),
        }
    )
    return (
        {
            "clinvar_significance": significance,
            "clinvar_has_conflicts": has_conflicts,
            "clinvar_conflict_summary": "; ".join(summary_values) or None,
            "clinvar_review_stars": max(row["review_stars"] for row in rows),
            "clinvar_submitters_count": len(submitters),
            "clinvar_id": ",".join(accessions),
        },
        rows,
        rsids[0] if rsids else None,
    )


class SortedClinVarLookup:
    """Merge-join exact alleles without loading either genome-wide VCF in memory."""

    def __init__(self, path: str | Path, snapshot: str) -> None:
        self.path = Path(path)
        self.snapshot = snapshot
        try:
            self._handle = pysam.VariantFile(str(self.path))
        except (OSError, ValueError) as exc:
            raise ResourceError(f"cannot open normalized ClinVar VCF: {exc}") from exc
        self._iterator: Iterator[pysam.VariantRecord] = iter(self._handle)
        self._next: tuple[tuple[int, int], str, str, pysam.VariantRecord] | None = None
        self._previous_resource_locus: tuple[int, int] | None = None
        self._previous_query_locus: tuple[int, int] | None = None
        self._seen_query_alleles: set[tuple[str, str]] = set()
        self._current_locus: tuple[int, int] | None = None
        self._current_records: dict[tuple[str, str], list[pysam.VariantRecord]] = {}
        self._advance()

    def _advance(self) -> None:
        for record in self._iterator:
            key = _record_key(record)
            if key is None:
                continue
            locus = key[:2]
            if self._previous_resource_locus is not None and locus < self._previous_resource_locus:
                raise ResourceError("normalized ClinVar VCF is not coordinate sorted")
            self._previous_resource_locus = locus
            self._next = locus, key[2], key[3], record.copy()
            return
        self._next = None

    def _consume_group(
        self,
    ) -> tuple[tuple[int, int], dict[tuple[str, str], list[pysam.VariantRecord]]]:
        assert self._next is not None
        locus = self._next[0]
        records: dict[tuple[str, str], list[pysam.VariantRecord]] = {}
        while self._next is not None and self._next[0] == locus:
            _, ref, alt, record = self._next
            records.setdefault((ref, alt), []).append(record)
            self._advance()
        return locus, records

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> list[pysam.VariantRecord]:
        locus = CHROM_RANK[chrom], pos
        allele = ref, alt
        if self._previous_query_locus is not None and locus < self._previous_query_locus:
            raise InputValidationError(
                "variant stream is not coordinate sorted for ClinVar merge join"
            )
        if locus != self._previous_query_locus:
            self._seen_query_alleles.clear()
        elif allele in self._seen_query_alleles:
            raise InputValidationError("duplicate query allele reached ClinVar merge join")
        self._previous_query_locus = locus
        self._seen_query_alleles.add(allele)

        if self._current_locus is not None and self._current_locus < locus:
            self._current_locus = None
            self._current_records = {}
        while self._next is not None and self._next[0] < locus:
            self._consume_group()
        if self._current_locus != locus and self._next is not None and self._next[0] == locus:
            self._current_locus, self._current_records = self._consume_group()
        if self._current_locus != locus:
            return []
        return list(self._current_records.get(allele, []))

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> SortedClinVarLookup:
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        self.close()
