"""Streaming, exact-allele projection of a frozen GWAS Catalog resource."""

from __future__ import annotations

import hashlib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from .constants import GWAS_CATALOG_SNAPSHOT, STANDARD_CHROMS
from .errors import InputValidationError, ResourceError

CHROM_RANK = {chrom: index for index, chrom in enumerate(STANDARD_CHROMS)}


def _association_id(
    snapshot: str,
    source_digest: str,
    source_row_number: int,
    variant_id: str,
) -> str:
    payload = f"NHGRI-EBI-GWAS\0{snapshot}\0{source_digest}\0{source_row_number}\0{variant_id}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _resource_rows(path: Path) -> Iterator[dict[str, Any]]:
    try:
        parquet = pq.ParquetFile(path)
        for batch in parquet.iter_batches(batch_size=32_768):
            yield from batch.to_pylist()
    except (OSError, ValueError) as exc:
        raise ResourceError(f"cannot stream normalized GWAS Catalog resource: {exc}") from exc


class SortedGwasCatalogLookup:
    """Merge-join coordinate-sorted carried alleles and GWAS resource rows."""

    def __init__(self, path: str | Path, snapshot: str = GWAS_CATALOG_SNAPSHOT) -> None:
        self.path = Path(path)
        self.snapshot = snapshot
        self._iterator = _resource_rows(self.path)
        self._next: tuple[tuple[int, int], dict[str, Any]] | None = None
        self._previous_resource_key: tuple[int, int] | None = None
        self._previous_query_key: tuple[int, int] | None = None
        self._seen_query_alleles: set[tuple[str, str]] = set()
        self._current_key: tuple[int, int] | None = None
        self._current_rows: list[dict[str, Any]] = []
        self._advance()

    def _advance(self) -> None:
        for row in self._iterator:
            chrom = str(row["chrom"])
            if chrom not in CHROM_RANK:
                raise ResourceError(f"non-core chromosome in GWAS resource: {chrom}")
            key = CHROM_RANK[chrom], int(row["pos"])
            if self._previous_resource_key is not None and key < self._previous_resource_key:
                raise ResourceError("normalized GWAS Catalog resource is not coordinate sorted")
            self._previous_resource_key = key
            self._next = key, row
            return
        self._next = None

    def _consume_group(self) -> tuple[tuple[int, int], list[dict[str, Any]]]:
        assert self._next is not None
        key = self._next[0]
        rows: list[dict[str, Any]] = []
        while self._next is not None and self._next[0] == key:
            rows.append(self._next[1])
            self._advance()
        return key, rows

    def lookup(
        self,
        *,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        rsid: str | None,
        gt: tuple[int, ...],
    ) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
        key = CHROM_RANK[chrom], pos
        allele = ref, alt
        if self._previous_query_key is not None and key < self._previous_query_key:
            raise InputValidationError("variant stream is not sorted for GWAS merge join")
        if key != self._previous_query_key:
            self._seen_query_alleles.clear()
        elif allele in self._seen_query_alleles:
            raise InputValidationError("duplicate query allele reached GWAS merge join")
        self._previous_query_key = key
        self._seen_query_alleles.add(allele)

        if self._current_key is not None and self._current_key < key:
            self._current_key = None
            self._current_rows = []
        while self._next is not None and self._next[0] < key:
            self._consume_group()
        if self._current_key != key and self._next is not None and self._next[0] == key:
            self._current_key, self._current_rows = self._consume_group()
        if self._current_key != key or rsid is None:
            return _blank_core(), [], []

        query_rsid = rsid.lower()
        variant_id = f"{chrom}:{pos}:{ref}:{alt}"
        evidence: list[dict[str, Any]] = []
        normative: list[dict[str, Any]] = []
        for resource in self._current_rows:
            aliases = [str(value).lower() for value in resource["rsid_aliases"]]
            effect_allele = str(resource["effect_allele"]).upper()
            if query_rsid not in aliases or effect_allele not in {ref, alt}:
                continue
            effect_index = 0 if effect_allele == ref else 1
            dosage = sum(int(value == effect_index) for value in gt)
            identifier = _association_id(
                self.snapshot,
                str(resource["source_record_digest"]),
                int(resource["source_row_number"]),
                variant_id,
            )
            emitted = resource["effect_size"] is not None and resource["effect_type"] is not None
            omission_reason = None if emitted else "ambiguous_or_missing_effect_type_or_size"
            pmids = [str(resource["pmid"])]
            evidence_row = {
                "association_id": identifier,
                "variant_id": variant_id,
                "risk_rsid": str(resource["risk_rsid"]),
                "rsid_aliases": aliases,
                "effect_allele": effect_allele,
                "effect_allele_dosage": dosage,
                "reported_trait": str(resource["reported_trait"]),
                "mapped_traits": list(resource["mapped_traits"]),
                "mapped_trait_uris": list(resource["mapped_trait_uris"]),
                "p_value": float(resource["p_value"]),
                "effect_size": (
                    float(resource["effect_size"]) if resource["effect_size"] is not None else None
                ),
                "effect_type": resource["effect_type"],
                "effect_type_basis": resource["effect_type_basis"],
                "ci_lower": (
                    float(resource["ci_lower"]) if resource["ci_lower"] is not None else None
                ),
                "ci_upper": (
                    float(resource["ci_upper"]) if resource["ci_upper"] is not None else None
                ),
                "sample_size": (
                    int(resource["sample_size"]) if resource["sample_size"] is not None else None
                ),
                "study_pmids": pmids,
                "study_accession": str(resource["study_accession"]),
                "source_version": self.snapshot,
                "source_row_number": int(resource["source_row_number"]),
                "source_record_digest": str(resource["source_record_digest"]),
                "raw_effect_text": resource["raw_effect_text"],
                "raw_ci_text": resource["raw_ci_text"],
                "normative_row_emitted": emitted,
                "omission_reason": omission_reason,
            }
            evidence.append(evidence_row)
            if emitted:
                normative.append(
                    {
                        "association_id": identifier,
                        "variant_id": variant_id,
                        "trait": str(resource["reported_trait"]),
                        "effect_allele": effect_allele,
                        "effect_size": float(resource["effect_size"]),
                        "effect_type": str(resource["effect_type"]),
                        "p_value": float(resource["p_value"]),
                        "sample_size": evidence_row["sample_size"],
                        "study_pmids": pmids,
                        "source_version": self.snapshot,
                    }
                )
        evidence.sort(key=lambda row: (row["p_value"], row["association_id"]))
        normative.sort(key=lambda row: (row["p_value"], row["association_id"]))
        return _aggregate_core(evidence), evidence, normative

    def close(self) -> None:
        self._iterator.close()

    def __enter__(self) -> SortedGwasCatalogLookup:
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()


def _blank_core() -> dict[str, Any]:
    return {
        "is_gwas_hit": False,
        "traits": None,
        "effect_size": None,
        "effect_type": None,
        "effect_allele": None,
        "ci_lower": None,
        "ci_upper": None,
        "p_value": None,
        "sample_size": None,
        "study_pmids": None,
    }


def _aggregate_core(evidence: list[dict[str, Any]]) -> dict[str, Any]:
    if not evidence:
        return _blank_core()
    traits = sorted(
        {
            str(value)
            for row in evidence
            for value in (row["mapped_traits"] or [row["reported_trait"]])
            if value
        }
    )
    pmids = sorted({value for row in evidence for value in row["study_pmids"]})
    typed = [row for row in evidence if row["normative_row_emitted"]]
    selected = (
        min(typed, key=lambda row: (row["p_value"], row["association_id"])) if typed else None
    )
    p_value = (
        float(selected["p_value"])
        if selected is not None
        else min(float(row["p_value"]) for row in evidence)
    )
    sample_size = selected["sample_size"] if selected is not None else None
    if sample_size is not None and sample_size > 2_147_483_647:
        sample_size = None
    return {
        "is_gwas_hit": True,
        "traits": traits or None,
        "effect_size": selected["effect_size"] if selected is not None else None,
        "effect_type": selected["effect_type"] if selected is not None else None,
        "effect_allele": selected["effect_allele"] if selected is not None else None,
        "ci_lower": selected["ci_lower"] if selected is not None else None,
        "ci_upper": selected["ci_upper"] if selected is not None else None,
        "p_value": p_value,
        "sample_size": sample_size,
        "study_pmids": pmids or None,
    }


__all__ = ["SortedGwasCatalogLookup"]
