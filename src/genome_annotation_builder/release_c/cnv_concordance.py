"""Interval-aware CNV concordance without intersection, voting, or CN rewriting."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, deque
from dataclasses import dataclass
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


@dataclass(frozen=True)
class CnvConcordancePolicy:
    min_reciprocal_overlap: float = 0.5
    max_copy_number_difference: float = 0.75


DEFAULT_CNV_CONCORDANCE_POLICY = CnvConcordancePolicy()


def _configuration_key(row: dict[str, Any]) -> tuple[str, str, str, str, int, str]:
    return (
        str(row["caller_name"]),
        str(row["caller_model"]),
        str(row["technical_source_role"]),
        str(row["source_alignment_sha256"]),
        int(row["bin_size"] or 0),
        str(row["caller_run_identity_sha256"]),
    )


def _caller_key(row: dict[str, Any]) -> tuple[str, str]:
    return str(row["caller_name"]), str(row["caller_version"])


def _reciprocal_overlap(first: dict[str, Any], second: dict[str, Any]) -> float:
    if first["chrom"] != second["chrom"]:
        return 0.0
    start = max(int(first["start_pos"]), int(second["start_pos"]))
    end = min(int(first["end_pos"]), int(second["end_pos"]))
    overlap = max(0, end - start + 1)
    first_size = int(first["end_pos"]) - int(first["start_pos"]) + 1
    second_size = int(second["end_pos"]) - int(second["start_pos"]) + 1
    return min(overlap / first_size, overlap / second_size)


def _biological_digest(rows: list[dict[str, Any]]) -> str:
    excluded = {
        "evidence_id",
        "filter_reasons",
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


def _group_id(rows: list[dict[str, Any]]) -> str:
    payload = sorted(
        (
            str(row["chrom"]),
            int(row["start_pos"]),
            int(row["end_pos"]),
            str(row["call_type"]),
            str(row["caller_name"]),
            str(row["caller_model"]),
            int(row["bin_size"] or 0),
        )
        for row in rows
    )
    rendered = json.dumps(payload, separators=(",", ":")).encode()
    return f"cnv-concordance:{hashlib.sha256(rendered).hexdigest()}"


def annotate_cnv_concordance(
    input_rows: list[dict[str, Any]],
    *,
    policy: CnvConcordancePolicy = DEFAULT_CNV_CONCORDANCE_POLICY,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not 0 < policy.min_reciprocal_overlap <= 1 or policy.max_copy_number_difference < 0:
        raise PortableEvidenceError("invalid CNV concordance policy")
    rows = validate_rows(
        "copy_number",
        input_rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    configurations = {_configuration_key(row) for row in rows}
    if len(configurations) < 2:
        raise PortableEvidenceError("CNV concordance requires at least two configurations")
    biological_before = _biological_digest(rows)
    adjacency: dict[int, set[int]] = {index: set() for index in range(len(rows))}
    for first_index, first in enumerate(rows):
        for second_index in range(first_index + 1, len(rows)):
            second = rows[second_index]
            if _configuration_key(first) == _configuration_key(second):
                continue
            if _reciprocal_overlap(first, second) >= policy.min_reciprocal_overlap:
                adjacency[first_index].add(second_index)
                adjacency[second_index].add(first_index)
    components: list[list[int]] = []
    remaining = set(adjacency)
    while remaining:
        seed = min(remaining)
        queue = deque([seed])
        component: list[int] = []
        remaining.remove(seed)
        while queue:
            current = queue.popleft()
            component.append(current)
            for neighbor in sorted(adjacency[current]):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    queue.append(neighbor)
        components.append(sorted(component))
    annotated = [dict(row) for row in rows]
    status_counts: Counter[str] = Counter()
    cross_caller_concordant_components = 0
    same_caller_cross_configuration_concordant_components = 0
    for indices in components:
        component_rows = [rows[index] for index in indices]
        component_configurations = [_configuration_key(row) for row in component_rows]
        if len(indices) == 1:
            status = "source_specific"
        elif len(set(component_configurations)) != len(component_configurations):
            status = "comparable_discordant"
        else:
            complete = all(
                _reciprocal_overlap(component_rows[first], component_rows[second])
                >= policy.min_reciprocal_overlap
                for first in range(len(component_rows))
                for second in range(first + 1, len(component_rows))
            )
            call_types = {str(row["call_type"]) for row in component_rows}
            copy_numbers = [row["copy_number"] for row in component_rows]
            comparable_copy_numbers = all(value is not None for value in copy_numbers)
            copy_numbers_agree = comparable_copy_numbers and all(
                abs(float(first) - float(second)) <= policy.max_copy_number_difference
                for offset, first in enumerate(copy_numbers)
                for second in copy_numbers[offset + 1 :]
            )
            if not complete or len(call_types) != 1 or not copy_numbers_agree:
                status = "comparable_discordant"
            else:
                status = "comparable_concordant"
                callers = {_caller_key(row) for row in component_rows}
                if len(callers) > 1:
                    cross_caller_concordant_components += 1
                else:
                    same_caller_cross_configuration_concordant_components += 1
        group_id = _group_id(component_rows)
        component_callers = {_caller_key(row) for row in component_rows}
        for index in indices:
            row = annotated[index]
            row["concordance_group_id"] = group_id
            row["concordance_status"] = status
            limitations = list(row["limitations"])
            filters = list(row["filter_reasons"])
            if status == "source_specific":
                limitations.append(
                    "no overlapping call was found in another supplied CNV configuration; absence is not refutation"
                )
            elif status == "comparable_concordant":
                if len(component_callers) > 1:
                    limitations.append(
                        "independent CNV callers agree within declared interval and copy-number tolerances"
                    )
                else:
                    limitations.append(
                        "one CNV caller agrees across configurations; this is not independent caller replication"
                    )
            else:
                limitations.append(
                    "overlapping CNV configurations disagree in type, copy number, interval grouping, or multiplicity"
                )
                filters.append("CROSS_CONFIGURATION_CNV_DISCORDANT")
                row["confidence_tier"] = "REVIEW"
            row["filter_reasons"] = sorted(set(filters))
            row["limitations"] = sorted(set(limitations))
            row["evidence_id"] = make_evidence_id("copy_number", row)
            status_counts[status] += 1
    validated = validate_rows(
        "copy_number",
        annotated,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_after = _biological_digest(validated)
    if biological_before != biological_after:
        raise PortableEvidenceError("CNV concordance changed biological call fields")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "input_rows": len(rows),
        "output_rows": len(validated),
        "callers": len({_caller_key(row) for row in rows}),
        "configurations": len(configurations),
        "components": len(components),
        "candidate_overlap_edges": sum(len(neighbors) for neighbors in adjacency.values()) // 2,
        "row_status_counts": dict(sorted(status_counts.items())),
        "cross_caller_concordant_components": cross_caller_concordant_components,
        "same_caller_cross_configuration_concordant_components": (
            same_caller_cross_configuration_concordant_components
        ),
        "majority_vote_applied": False,
        "source_rows_removed": 0,
        "biological_fields_before_sha256": biological_before,
        "biological_fields_after_sha256": biological_after,
        "biological_fields_changed": False,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "private_paths_persisted": False,
    }
    return validated, summary


def read_cnv_parquets(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        table = pq.read_table(path)
        if table.schema != schema_for("copy_number"):
            raise PortableEvidenceError("CNV concordance input has the wrong physical schema")
        rows.extend(table.to_pylist())
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-parquet", type=Path, action="append", required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--min-reciprocal-overlap", type=float, default=0.5)
    parser.add_argument("--max-copy-number-difference", type=float, default=0.75)
    args = parser.parse_args()
    policy = CnvConcordancePolicy(
        min_reciprocal_overlap=args.min_reciprocal_overlap,
        max_copy_number_difference=args.max_copy_number_difference,
    )
    rows = read_cnv_parquets(args.input_parquet)
    annotated, summary = annotate_cnv_concordance(
        rows,
        policy=policy,
        expected_reference_sha256=args.reference_sha256,
    )
    output = write_extension(
        "copy_number",
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
