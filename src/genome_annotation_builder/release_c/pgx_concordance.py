"""Conservative complex-PGx concordance without voting or source-row loss."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
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

STAR_TOKEN = re.compile(
    r"^\*(?P<major>[0-9]+)(?:\.(?P<suballele>[0-9]+))?(?:x(?P<copies>[1-9][0-9]*))?$"
)
MAX_EXPLICIT_COPIES = 50


def _source_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        str(row["caller_name"]),
        str(row["technical_source_role"]),
        str(row["source_alignment_sha256"]),
    )


def _canonical_major_structure(diplotype: str) -> str | None:
    if not isinstance(diplotype, str) or diplotype != diplotype.strip():
        return None
    sides = diplotype.split("/")
    if len(sides) != 2 or any(not side for side in sides):
        return None
    canonical_sides: list[tuple[str, ...]] = []
    for side in sides:
        alleles: list[str] = []
        for token in side.split("+"):
            match = STAR_TOKEN.fullmatch(token)
            if match is None:
                return None
            copies = int(match.group("copies") or "1")
            if copies > MAX_EXPLICIT_COPIES:
                return None
            alleles.extend([f"*{match.group('major')}"] * copies)
        canonical_sides.append(tuple(sorted(alleles)))
    canonical_sides.sort()
    return "/".join("+".join(side) for side in canonical_sides)


def _called_key(row: dict[str, Any]) -> str | None:
    if row["call_status"] != "called" or row["diplotype"] is None:
        return None
    return _canonical_major_structure(str(row["diplotype"]))


def _candidate_keys(row: dict[str, Any]) -> frozenset[str] | None:
    candidates = row["candidate_diplotypes"]
    if not candidates:
        return None
    keys = [_canonical_major_structure(str(candidate)) for candidate in candidates]
    if any(key is None for key in keys):
        return None
    return frozenset(str(key) for key in keys)


def _group_id(gene: str) -> str:
    rendered = json.dumps(
        {"gene": gene, "boundary": "major-star-structural"},
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return f"complex-pgx-concordance:{hashlib.sha256(rendered).hexdigest()}"


def _biological_digest(rows: list[dict[str, Any]]) -> str:
    excluded = {
        "evidence_id",
        "filter_reasons",
        "confidence_tier",
        "limitations",
        "outside_call_eligible",
        "concordance_group_id",
        "concordance_status",
    }
    payload = [
        {key: value for key, value in sorted(row.items()) if key not in excluded} for row in rows
    ]
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode()).hexdigest()


def annotate_pgx_concordance(
    input_rows: list[dict[str, Any]],
    *,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows = validate_rows(
        "complex_pharmacogenes",
        input_rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    sources = {_source_key(row) for row in rows}
    if len(sources) < 2:
        raise PortableEvidenceError(
            "complex-PGX concordance requires at least two independent sources"
        )
    biological_before = _biological_digest(rows)
    by_gene: dict[str, list[int]] = {}
    seen_source_gene: set[tuple[tuple[str, str, str], str]] = set()
    for index, row in enumerate(rows):
        source_gene = (_source_key(row), str(row["gene"]))
        if source_gene in seen_source_gene:
            raise PortableEvidenceError("one complex-PGX source contains duplicate gene rows")
        seen_source_gene.add(source_gene)
        by_gene.setdefault(str(row["gene"]), []).append(index)
    annotated = [dict(row) for row in rows]
    status_counts: Counter[str] = Counter()
    candidate_overlap_groups = 0
    definition_boundary_mismatch_groups = 0
    for gene, indices in sorted(by_gene.items()):
        gene_rows = [rows[index] for index in indices]
        called_keys = [_called_key(row) for row in gene_rows]
        if len(gene_rows) == 1:
            status = "source_specific"
        elif all(row["call_status"] == "called" for row in gene_rows) and all(
            key is not None for key in called_keys
        ):
            status = (
                "comparable_concordant" if len(set(called_keys)) == 1 else "comparable_discordant"
            )
        else:
            status = "not_comparable"
        candidate_sets = [_candidate_keys(row) for row in gene_rows]
        candidate_overlap = (
            len(gene_rows) >= 2
            and all(keys is not None for keys in candidate_sets)
            and bool(set.intersection(*(set(keys) for keys in candidate_sets if keys is not None)))
        )
        if candidate_overlap and status == "not_comparable":
            candidate_overlap_groups += 1
        definitions_differ = len({row["definition_identity_sha256"] for row in gene_rows}) > 1
        if definitions_differ:
            definition_boundary_mismatch_groups += 1
        group_id = _group_id(gene)
        for index in indices:
            row = annotated[index]
            row["concordance_group_id"] = group_id
            row["concordance_status"] = status
            limitations = list(row["limitations"])
            filters = list(row["filter_reasons"])
            if status == "source_specific":
                limitations.append(
                    "no second supplied complex-PGX source covers this gene; absence is not refutation"
                )
            elif status == "not_comparable":
                limitations.append(
                    "complex-PGX sources lack complete unambiguous calls at a shared major-star structural boundary"
                )
                if candidate_overlap:
                    limitations.append(
                        "candidate sets overlap, but ambiguity is retained and no candidate is selected"
                    )
            elif status == "comparable_concordant":
                limitations.append(
                    "independent complex-PGX calls agree at the major-star structural boundary"
                )
            else:
                limitations.append(
                    "independent complex-PGX calls disagree at the major-star structural boundary"
                )
                filters.append("CROSS_CALLER_COMPLEX_PGX_DISCORDANT")
                row["confidence_tier"] = "REVIEW"
                row["outside_call_eligible"] = False
            if definitions_differ:
                limitations.append(
                    "caller definition identities differ; concordance does not imply suballele-definition equivalence"
                )
            row["filter_reasons"] = sorted(set(filters))
            row["limitations"] = sorted(set(limitations))
            row["evidence_id"] = make_evidence_id("complex_pharmacogenes", row)
            status_counts[status] += 1
    validated = validate_rows(
        "complex_pharmacogenes",
        annotated,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_after = _biological_digest(validated)
    if biological_before != biological_after:
        raise PortableEvidenceError("complex-PGX concordance changed biological call fields")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "comparison_boundary": "major-star-structural",
        "input_rows": len(rows),
        "output_rows": len(validated),
        "sources": len(sources),
        "genes": len(by_gene),
        "row_status_counts": dict(sorted(status_counts.items())),
        "candidate_overlap_groups_retained_ambiguous": candidate_overlap_groups,
        "definition_boundary_mismatch_groups": definition_boundary_mismatch_groups,
        "majority_vote_applied": False,
        "source_rows_removed": 0,
        "biological_fields_before_sha256": biological_before,
        "biological_fields_after_sha256": biological_after,
        "biological_fields_changed": False,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "private_paths_persisted": False,
    }
    return validated, summary


def read_pgx_parquets(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        table = pq.read_table(path)
        if table.schema != schema_for("complex_pharmacogenes"):
            raise PortableEvidenceError(
                "complex-PGX concordance input has the wrong physical schema"
            )
        rows.extend(table.to_pylist())
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-parquet", type=Path, action="append", required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    args = parser.parse_args()
    rows = read_pgx_parquets(args.input_parquet)
    annotated, summary = annotate_pgx_concordance(
        rows, expected_reference_sha256=args.reference_sha256
    )
    output = write_extension(
        "complex_pharmacogenes",
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
