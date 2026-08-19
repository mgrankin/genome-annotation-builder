"""Replicate concordance for bounded local read-back observations.

Every technical-alignment row is retained.  Concordance is an annotation over
finite REF/ALT fragment counts; it is never a genotype vote and never treats
an absent or unassessed observation as reference evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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
class ReadbackReplicatePolicy:
    minimum_ref_alt_fragments: int = 10
    wilson_z: float = 1.96


DEFAULT_READBACK_REPLICATE_POLICY = ReadbackReplicatePolicy()
POLICY_FIELDS = (
    "minimum_mapping_quality",
    "minimum_base_quality",
    "duplicate_reads_included",
    "overlapping_mate_policy",
    "nearby_indel_window",
)
DEFINITION_FIELDS = (
    "target_manifest_sha256",
    "target_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_type",
    "local_complexity_status",
    "local_haplotype_status",
    "candidate_category",
)
ANNOTATION_FIELDS = frozenset(
    {
        "evidence_id",
        "filter_reasons",
        "confidence_tier",
        "limitations",
        "concordance_group_id",
        "concordance_status",
    }
)


def _wilson_interval(successes: int, total: int, z: float) -> tuple[float, float]:
    proportion = successes / total
    z_squared = z * z
    denominator = 1 + z_squared / total
    centre = (proportion + z_squared / (2 * total)) / denominator
    margin = (
        z
        * math.sqrt(proportion * (1 - proportion) / total + z_squared / (4 * total * total))
        / denominator
    )
    return max(0.0, centre - margin), min(1.0, centre + margin)


def _pair_relation(
    first: dict[str, Any],
    second: dict[str, Any],
    policy: ReadbackReplicatePolicy,
) -> str:
    if first["assessment_status"] != "assessed" or second["assessment_status"] != "assessed":
        return "not_comparable"
    first_total = int(first["ref_reads"]) + int(first["alt_reads"])
    second_total = int(second["ref_reads"]) + int(second["alt_reads"])
    if min(first_total, second_total) < policy.minimum_ref_alt_fragments:
        return "not_comparable"
    first_interval = _wilson_interval(int(first["alt_reads"]), first_total, policy.wilson_z)
    second_interval = _wilson_interval(int(second["alt_reads"]), second_total, policy.wilson_z)
    if first_interval[1] < second_interval[0] or second_interval[1] < first_interval[0]:
        return "discordant"
    return "interval_compatible"


def _biological_digest(rows: list[dict[str, Any]]) -> str:
    payload = sorted(
        (
            {key: value for key, value in sorted(row.items()) if key not in ANNOTATION_FIELDS}
            for row in rows
        ),
        key=lambda row: (
            row["target_manifest_sha256"],
            row["target_id"],
            row["source_alignment_sha256"],
            row["source_record_digest"],
        ),
    )
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode()).hexdigest()


def _group_id(row: dict[str, Any]) -> str:
    payload = {field: row[field] for field in DEFINITION_FIELDS}
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return f"readback-concordance:{hashlib.sha256(rendered.encode()).hexdigest()}"


def _degrade_to_review(row: dict[str, Any]) -> None:
    if row["confidence_tier"] == "SUPPORTED":
        row["confidence_tier"] = "REVIEW"


def annotate_readback_replicate_concordance(
    input_rows: list[dict[str, Any]],
    *,
    policy: ReadbackReplicatePolicy = DEFAULT_READBACK_REPLICATE_POLICY,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Annotate the same bounded targets across independent alignments."""

    if (
        policy.minimum_ref_alt_fragments < 1
        or not math.isfinite(policy.wilson_z)
        or policy.wilson_z <= 0
    ):
        raise PortableEvidenceError("invalid read-back replicate policy")
    rows = validate_rows(
        "read_backed_evidence",
        input_rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    if any(row["caller_name"] != "GAB bounded local pileup" for row in rows):
        raise PortableEvidenceError("read-back concordance requires bounded pileup rows")
    alignments = {str(row["source_alignment_sha256"]) for row in rows}
    if len(alignments) < 2:
        raise PortableEvidenceError("read-back concordance requires two technical alignments")
    target_manifests = {str(row["target_manifest_sha256"]) for row in rows}
    if len(target_manifests) != 1:
        raise PortableEvidenceError("read-back replicates use different target manifests")

    biological_before = _biological_digest(rows)
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[str(row["target_id"])].append(index)

    annotated = [dict(row) for row in rows]
    status_counts: Counter[str] = Counter()
    relation_counts: Counter[str] = Counter()
    for _target_id, indices in sorted(grouped.items()):
        component = [rows[index] for index in indices]
        component_alignments = [str(row["source_alignment_sha256"]) for row in component]
        if len(component_alignments) != len(set(component_alignments)):
            raise PortableEvidenceError(
                "one read-back target has multiple rows from the same alignment"
            )
        definitions = {tuple(row[field] for field in DEFINITION_FIELDS) for row in component}
        if len(definitions) != 1:
            raise PortableEvidenceError("one read-back target ID has conflicting definitions")
        boundaries = {
            (
                str(row["caller_version"]),
                str(row["resource_manifest_sha256"]),
                *(row[field] for field in POLICY_FIELDS),
            )
            for row in component
        }
        pair_relations = [
            _pair_relation(component[first], component[second], policy)
            for first in range(len(component))
            for second in range(first + 1, len(component))
        ]
        if len(component) != len(alignments):
            status = "source_specific"
            relation = "missing_source_observation"
        elif len(boundaries) != 1:
            status = "not_comparable"
            relation = "caller_resource_or_policy_mismatch"
        elif "discordant" in pair_relations:
            status = "comparable_discordant"
            relation = "alt_fraction_interval_discordant"
        elif "not_comparable" in pair_relations:
            status = "not_comparable"
            relation = "unassessed_or_insufficient_fragments"
        else:
            status = "comparable_concordant"
            relation = "alt_fraction_interval_compatible"
        group_id = _group_id(component[0])
        relation_counts[relation] += 1
        for index in indices:
            row = annotated[index]
            row["concordance_group_id"] = group_id
            row["concordance_status"] = status
            limitations = list(row["limitations"])
            filters = list(row["filter_reasons"])
            if status == "source_specific":
                limitations.append(
                    "another supplied alignment lacks this target row; absence is not reference evidence"
                )
            elif status == "not_comparable":
                limitations.append(
                    "technical read-back observations lack a common audited boundary or sufficient finite REF/ALT fragments"
                )
                filters.append("REPLICATE:READBACK_NOT_COMPARABLE")
                _degrade_to_review(row)
            elif status == "comparable_discordant":
                limitations.append(
                    "technical read-back ALT-fraction Wilson intervals are discordant"
                )
                filters.append("REPLICATE:READBACK_ALT_FRACTION_DISCORDANT")
                _degrade_to_review(row)
            else:
                limitations.append(
                    "technical read-back ALT-fraction Wilson intervals overlap; this is not genotype confirmation"
                )
            row["filter_reasons"] = sorted(set(filters))
            row["limitations"] = sorted(set(limitations))
            row["evidence_id"] = make_evidence_id("read_backed_evidence", row)
            status_counts[status] += 1

    validated = validate_rows(
        "read_backed_evidence",
        annotated,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_after = _biological_digest(validated)
    if biological_before != biological_after:
        raise PortableEvidenceError("read-back concordance changed biological observation fields")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "policy": {
            "minimum_ref_alt_fragments": policy.minimum_ref_alt_fragments,
            "wilson_z": policy.wilson_z,
        },
        "input_rows": len(rows),
        "output_rows": len(validated),
        "technical_alignments": len(alignments),
        "targets": len(grouped),
        "pairwise_comparisons": sum(
            len(indices) * (len(indices) - 1) // 2 for indices in grouped.values()
        ),
        "group_relation_counts": dict(sorted(relation_counts.items())),
        "row_status_counts": dict(sorted(status_counts.items())),
        "source_rows_removed": 0,
        "allele_fraction_averaging_used": False,
        "majority_vote_used": False,
        "absence_used_as_reference": False,
        "biological_fields_before_sha256": biological_before,
        "biological_fields_after_sha256": biological_after,
        "biological_fields_changed": False,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def read_readback_parquets(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        table = pq.read_table(path)
        if table.schema != schema_for("read_backed_evidence"):
            raise PortableEvidenceError("read-back concordance input has the wrong schema")
        rows.extend(table.to_pylist())
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-parquet", type=Path, action="append", required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument(
        "--minimum-ref-alt-fragments",
        type=int,
        default=DEFAULT_READBACK_REPLICATE_POLICY.minimum_ref_alt_fragments,
    )
    parser.add_argument(
        "--wilson-z", type=float, default=DEFAULT_READBACK_REPLICATE_POLICY.wilson_z
    )
    args = parser.parse_args()
    rows = read_readback_parquets(args.input_parquet)
    policy = ReadbackReplicatePolicy(
        minimum_ref_alt_fragments=args.minimum_ref_alt_fragments,
        wilson_z=args.wilson_z,
    )
    annotated, summary = annotate_readback_replicate_concordance(
        rows,
        policy=policy,
        expected_reference_sha256=args.reference_sha256,
    )
    output = write_extension(
        "read_backed_evidence",
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
