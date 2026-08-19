"""Technical-replicate concordance for typed short-read repeat estimates."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter, defaultdict
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
class RepeatEstimate:
    repeat_count: int
    ci_lower: int
    ci_upper: int | None
    count_is_lower_bound: bool

    def exact_identity(self) -> tuple[int, int, int | None, bool]:
        return (
            self.repeat_count,
            self.ci_lower,
            self.ci_upper,
            self.count_is_lower_bound,
        )


DEFINITION_FIELDS = (
    "locus_id",
    "variant_id",
    "variant_subtype",
    "chrom",
    "start_pos",
    "end_pos",
    "motif",
    "ploidy",
    "reference_repeat_count",
)
BIOLOGICAL_FIELDS = (
    *DEFINITION_FIELDS,
    "allele1_repeat_count",
    "allele1_ci_lower",
    "allele1_ci_upper",
    "allele1_count_is_lower_bound",
    "allele2_repeat_count",
    "allele2_ci_lower",
    "allele2_ci_upper",
    "allele2_count_is_lower_bound",
    "source_alignment_sha256",
    "source_record_digest",
)


def _alleles(row: dict[str, Any]) -> tuple[RepeatEstimate, ...]:
    result: list[RepeatEstimate] = []
    for prefix in ("allele1", "allele2"):
        count = row[f"{prefix}_repeat_count"]
        if count is None:
            continue
        result.append(
            RepeatEstimate(
                repeat_count=int(count),
                ci_lower=int(row[f"{prefix}_ci_lower"]),
                ci_upper=(
                    int(row[f"{prefix}_ci_upper"])
                    if row[f"{prefix}_ci_upper"] is not None
                    else None
                ),
                count_is_lower_bound=bool(row[f"{prefix}_count_is_lower_bound"]),
            )
        )
    if len(result) != int(row["ploidy"]):
        raise PortableEvidenceError("repeat allele cardinality differs from ploidy")
    return tuple(result)


def _allele_relation(first: RepeatEstimate, second: RepeatEstimate) -> str:
    if first.exact_identity() == second.exact_identity():
        return "exact"
    if first.ci_upper is None or second.ci_upper is None:
        return "not_comparable"
    overlap = max(first.ci_lower, second.ci_lower) <= min(first.ci_upper, second.ci_upper)
    return "interval_compatible" if overlap else "discordant"


def _genotype_relation(first: dict[str, Any], second: dict[str, Any]) -> str:
    first_alleles = _alleles(first)
    second_alleles = _alleles(second)
    if len(first_alleles) != len(second_alleles):
        return "not_comparable"
    candidates: list[tuple[tuple[int, int, int, int], tuple[str, ...]]] = []
    for permutation in itertools.permutations(second_alleles):
        relations = tuple(
            _allele_relation(left, right)
            for left, right in zip(first_alleles, permutation, strict=True)
        )
        compatibility_rank = (
            3
            if all(relation in {"exact", "interval_compatible"} for relation in relations)
            else 2
            if "discordant" not in relations
            else 1
        )
        score = (
            compatibility_rank,
            relations.count("exact"),
            relations.count("interval_compatible"),
            -relations.count("not_comparable"),
        )
        candidates.append((score, relations))
    _, best = max(candidates, key=lambda item: (item[0], item[1]))
    if all(relation == "exact" for relation in best):
        return "exact"
    if all(relation in {"exact", "interval_compatible"} for relation in best):
        return "interval_compatible"
    if "not_comparable" in best:
        return "not_comparable"
    return "discordant"


def _group_id(rows: list[dict[str, Any]]) -> str:
    payload = [
        {field: row[field] for field in BIOLOGICAL_FIELDS}
        for row in sorted(rows, key=lambda item: item["source_alignment_sha256"])
    ]
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return f"repeat-concordance:{hashlib.sha256(rendered.encode()).hexdigest()}"


def _biological_digest(rows: list[dict[str, Any]]) -> str:
    payload = sorted(
        ({field: row[field] for field in BIOLOGICAL_FIELDS} for row in rows),
        key=lambda item: (item["source_alignment_sha256"], item["source_record_digest"]),
    )
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode()).hexdigest()


def annotate_repeat_replicate_concordance(
    input_rows: list[dict[str, Any]],
    *,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Annotate replicate consistency without voting, deletion, or call replacement."""

    rows = validate_rows(
        "repeat_expansions",
        input_rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    if any(row["caller_name"] != "ExpansionHunter" for row in rows):
        raise PortableEvidenceError("repeat replicate concordance requires ExpansionHunter rows")
    alignments = {str(row["source_alignment_sha256"]) for row in rows}
    if len(alignments) < 2:
        raise PortableEvidenceError(
            "repeat replicate concordance requires two technical alignments"
        )
    biological_before = _biological_digest(rows)
    grouped: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[(str(row["locus_id"]), str(row["variant_id"]))].append(index)

    annotated = [dict(row) for row in rows]
    status_counts: Counter[str] = Counter()
    relation_counts: Counter[str] = Counter()
    for _, indices in sorted(grouped.items()):
        component = [rows[index] for index in indices]
        component_alignments = [str(row["source_alignment_sha256"]) for row in component]
        if len(component_alignments) != len(set(component_alignments)):
            raise PortableEvidenceError(
                "one repeat locus has multiple accepted rows from the same alignment"
            )
        definitions = {tuple(row[field] for field in DEFINITION_FIELDS) for row in component}
        caller_boundaries = {
            (
                str(row["caller_version"]),
                str(row["resource_manifest_sha256"]),
            )
            for row in component
        }
        pair_relations = [
            _genotype_relation(component[first], component[second])
            for first in range(len(component))
            for second in range(first + 1, len(component))
        ]
        if len(component) == 1:
            status = "source_specific"
            group_relation = "single_source"
        elif len(definitions) != 1 or len(caller_boundaries) != 1:
            status = "not_comparable"
            group_relation = "definition_or_resource_mismatch"
        elif "discordant" in pair_relations:
            status = "comparable_discordant"
            group_relation = "discordant"
        elif "not_comparable" in pair_relations:
            status = "not_comparable"
            group_relation = "open_ended_not_comparable"
        elif pair_relations and all(relation == "exact" for relation in pair_relations):
            status = "comparable_concordant"
            group_relation = "exact"
        else:
            status = "comparable_concordant"
            group_relation = "interval_compatible"
        group_id = _group_id(component)
        relation_counts[group_relation] += 1
        for index in indices:
            row = annotated[index]
            row["concordance_group_id"] = group_id
            row["concordance_status"] = status
            limitations = list(row["limitations"])
            filter_reasons = list(row["filter_reasons"])
            if status == "source_specific":
                limitations.append(
                    "no comparable repeat row was supplied by another technical alignment; absence is not refutation"
                )
            elif status == "not_comparable":
                limitations.append(
                    "technical-replicate repeat estimates are not comparable at a finite common boundary"
                )
                filter_reasons.append("REPLICATE:REPEAT_NOT_COMPARABLE")
                row["confidence_tier"] = "REVIEW"
            elif status == "comparable_discordant":
                limitations.append(
                    "technical-replicate repeat confidence intervals are genotype-discordant"
                )
                filter_reasons.append("REPLICATE:REPEAT_GENOTYPE_DISCORDANT")
                row["confidence_tier"] = "REVIEW"
            elif group_relation == "interval_compatible":
                limitations.append(
                    "technical replicates agree only within overlapping repeat confidence intervals"
                )
            row["filter_reasons"] = sorted(set(filter_reasons))
            row["limitations"] = sorted(set(limitations))
            row["evidence_id"] = make_evidence_id("repeat_expansions", row)
            status_counts[status] += 1

    validated = validate_rows(
        "repeat_expansions",
        annotated,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_after = _biological_digest(validated)
    if biological_before != biological_after:
        raise PortableEvidenceError("repeat concordance changed biological call fields")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "input_rows": len(rows),
        "output_rows": len(validated),
        "technical_alignments": len(alignments),
        "locus_variant_groups": len(grouped),
        "pairwise_comparisons": sum(
            len(indices) * (len(indices) - 1) // 2 for indices in grouped.values()
        ),
        "group_relation_counts": dict(sorted(relation_counts.items())),
        "row_status_counts": dict(sorted(status_counts.items())),
        "source_rows_removed": 0,
        "majority_vote_used": False,
        "biological_fields_before_sha256": biological_before,
        "biological_fields_after_sha256": biological_after,
        "biological_fields_changed": False,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def read_repeat_parquets(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        table = pq.read_table(path)
        if table.schema != schema_for("repeat_expansions"):
            raise PortableEvidenceError("repeat concordance input has the wrong physical schema")
        rows.extend(table.to_pylist())
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-parquet", type=Path, action="append", required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    args = parser.parse_args()
    rows = read_repeat_parquets(args.input_parquet)
    annotated, summary = annotate_repeat_replicate_concordance(
        rows, expected_reference_sha256=args.reference_sha256
    )
    output = write_extension(
        "repeat_expansions",
        annotated,
        args.output_parquet,
        expected_reference_sha256=args.reference_sha256,
    )
    print(
        json.dumps(
            {"status": "PASS", "concordance": summary, "output": output},
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
