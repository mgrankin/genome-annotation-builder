"""Resolution-aware HLA concordance without destructive voting or intersection."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from .portable_model import (
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    schema_for,
    validate_rows,
    write_extension,
)


def _source_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        str(row["caller_name"]),
        str(row["technical_source_role"]),
        str(row["source_alignment_sha256"]),
    )


def _allele_counter(row: dict[str, Any], boundary: str) -> Counter[str] | None:
    if row["call_status"] != "diploid_call":
        return None
    if boundary == "g-group":
        values = (row["allele1_g_group"], row["allele2_g_group"])
    elif boundary == "two-field":
        values = (row["allele1_two_field"], row["allele2_two_field"])
    else:
        raise PortableEvidenceError("unknown HLA concordance boundary")
    if any(value is None for value in values):
        return None
    return Counter(str(value) for value in values)


def _common_boundary(rows: list[dict[str, Any]]) -> str | None:
    if all(_allele_counter(row, "g-group") is not None for row in rows):
        return "g-group"
    if all(_allele_counter(row, "two-field") is not None for row in rows):
        return "two-field"
    return None


def _group_id(locus: str, boundary: str | None) -> str:
    payload = json.dumps(
        {"locus": locus, "boundary": boundary or "not-comparable"},
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return f"hla-concordance:{hashlib.sha256(payload).hexdigest()}"


def _biological_digest(rows: list[dict[str, Any]]) -> str:
    excluded = {
        "evidence_id",
        "confidence_tier",
        "limitations",
        "concordance_group_id",
        "concordance_status",
    }
    payload = [
        {key: value for key, value in sorted(row.items()) if key not in excluded} for row in rows
    ]
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode()).hexdigest()


def annotate_hla_concordance(
    input_rows: list[dict[str, Any]],
    *,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows = validate_rows(
        "hla_typing",
        input_rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    sources = {_source_key(row) for row in rows}
    if len(sources) < 2:
        raise PortableEvidenceError("HLA concordance requires at least two independent sources")
    biological_before = _biological_digest(rows)
    by_locus: dict[str, list[int]] = {}
    seen_source_locus: set[tuple[tuple[str, str, str], str]] = set()
    for index, row in enumerate(rows):
        source_locus = (_source_key(row), str(row["locus"]))
        if source_locus in seen_source_locus:
            raise PortableEvidenceError("one HLA source contains duplicate locus rows")
        seen_source_locus.add(source_locus)
        by_locus.setdefault(str(row["locus"]), []).append(index)
    annotated = [dict(row) for row in rows]
    status_counts: Counter[str] = Counter()
    boundary_counts: Counter[str] = Counter()
    protected_current_rows = 0
    for locus, indices in sorted(by_locus.items()):
        locus_rows = [rows[index] for index in indices]
        boundary = _common_boundary(locus_rows) if len(locus_rows) >= 2 else None
        counters = (
            [_allele_counter(row, boundary) for row in locus_rows] if boundary is not None else []
        )
        if len(locus_rows) == 1:
            status = "source_specific"
        elif boundary is None or any(counter is None for counter in counters):
            status = "not_comparable"
        elif all(counter == counters[0] for counter in counters[1:]):
            status = "comparable_concordant"
        else:
            status = "comparable_discordant"
        group_id = _group_id(locus, boundary)
        current_indices = [
            index for index in indices if not bool(annotated[index]["database_age_limitation"])
        ]
        current_counters = (
            [_allele_counter(annotated[index], boundary) for index in current_indices]
            if boundary is not None
            else []
        )
        current_sources_agree = bool(current_counters) and all(
            counter == current_counters[0] for counter in current_counters[1:]
        )
        for index in indices:
            row = annotated[index]
            row["concordance_group_id"] = group_id
            row["concordance_status"] = status
            limitations = list(row["limitations"])
            if status == "source_specific":
                limitations.append(
                    "no second supplied HLA source covers this locus; absence is not refutation"
                )
            elif status == "not_comparable":
                limitations.append(
                    "HLA sources lack a complete shared G-group or two-field comparison boundary"
                )
            elif status == "comparable_concordant":
                limitations.append(f"independent HLA sources agree at the {boundary} boundary")
            else:
                limitations.append(f"independent HLA sources disagree at the {boundary} boundary")
                if row["database_age_limitation"]:
                    row["confidence_tier"] = "REVIEW"
                    limitations.append(
                        "the historical database result cannot override a newer source"
                    )
                elif current_sources_agree:
                    protected_current_rows += 1
                    limitations.append(
                        "this current-database call is retained; discordant historical audit evidence does not override it"
                    )
                else:
                    row["confidence_tier"] = "REVIEW"
            row["limitations"] = sorted(set(limitations))
            row["evidence_id"] = make_evidence_id("hla_typing", row)
            status_counts[status] += 1
        boundary_counts[boundary or "none"] += len(indices)
    validated = validate_rows(
        "hla_typing",
        annotated,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_after = _biological_digest(validated)
    if biological_before != biological_after:
        raise PortableEvidenceError("HLA concordance changed biological call fields")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "input_rows": len(rows),
        "output_rows": len(validated),
        "sources": len(sources),
        "loci": len(by_locus),
        "row_status_counts": dict(sorted(status_counts.items())),
        "row_boundary_counts": dict(sorted(boundary_counts.items())),
        "current_database_rows_protected_from_historical_override": protected_current_rows,
        "biological_fields_before_sha256": biological_before,
        "biological_fields_after_sha256": biological_after,
        "biological_fields_changed": False,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "private_paths_persisted": False,
    }
    return validated, summary


def read_hla_parquets(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        table = pq.read_table(path)
        if table.schema != schema_for("hla_typing"):
            raise PortableEvidenceError("HLA concordance input has the wrong physical schema")
        rows.extend(table.to_pylist())
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-parquet", type=Path, action="append", required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    args = parser.parse_args()
    rows = read_hla_parquets(args.input_parquet)
    annotated, summary = annotate_hla_concordance(
        rows, expected_reference_sha256=args.reference_sha256
    )
    output = write_extension(
        "hla_typing",
        annotated,
        args.output_parquet,
        expected_reference_sha256=args.reference_sha256,
    )
    print(
        json.dumps(
            {"status": "PASS", "concordance": summary, "output": output},
            sort_keys=True,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
