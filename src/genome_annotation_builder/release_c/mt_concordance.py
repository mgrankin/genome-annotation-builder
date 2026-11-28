"""Technical-replicate concordance for typed mitochondrial heteroplasmy evidence."""

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

from .caller_receipt import canonical_json_sha256
from .portable_model import (
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    schema_for,
    sha256_file,
    validate_rows,
    write_extension,
)


@dataclass(frozen=True)
class MitochondrialReplicatePolicy:
    """Finite sampling boundary used only to assess replicate comparability."""

    minimum_allelic_depth: int = 20
    wilson_z: float = 1.959963984540054


DEFAULT_MT_REPLICATE_POLICY = MitochondrialReplicatePolicy()

MUTABLE_CONCORDANCE_FIELDS = frozenset(
    {
        "evidence_id",
        "filter_reasons",
        "confidence_tier",
        "limitations",
        "concordance_group_id",
        "concordance_status",
    }
)
MT_DEPTH_POLICY = "source_record_ref_plus_all_alt_AD"

BIOLOGICAL_FIELDS = tuple(
    field
    for field in schema_for("mitochondrial_heteroplasmy").names
    if field not in MUTABLE_CONCORDANCE_FIELDS
)


@dataclass(frozen=True)
class AlleleFractionObservation:
    reported_fraction: float
    alt_depth: int
    allelic_depth: int
    lower: float
    upper: float


def _wilson_interval(alt_depth: int, allelic_depth: int, z: float) -> tuple[float, float]:
    proportion = alt_depth / allelic_depth
    z_squared = z * z
    denominator = 1.0 + z_squared / allelic_depth
    centre = (proportion + z_squared / (2.0 * allelic_depth)) / denominator
    radius = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / allelic_depth
            + z_squared / (4.0 * allelic_depth * allelic_depth)
        )
        / denominator
    )
    return max(0.0, centre - radius), min(1.0, centre + radius)


def _source_records(rows: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    records = defaultdict(list)
    for row in rows:
        key = (
            row["source_alignment_sha256"],
            row["raw_artifact_sha256"],
            row["source_record_digest"],
            row["chrom"],
            row["pos"],
            row["ref"],
        )
        records[key].append(row)
    return list(records.values())


def verified_mt_records(rows: list[dict[str, Any]]) -> set[frozenset[str]]:
    """Identify complete sibling sets from an independently verified native partition."""
    return {frozenset(row["evidence_id"] for row in group) for group in _source_records(rows)}


def _allelic_depths(
    rows: list[dict[str, Any]], complete_records: set[frozenset[str]]
) -> dict[str, int | None]:
    """Count REF once across ALT siblings from the same physical Mutect2 record."""
    result = {}
    for siblings in _source_records(rows):
        counts = {(r["depth"], r["ref_depth"]) for r in siblings}
        if len(counts) != 1:
            raise PortableEvidenceError(
                "inconsistent mitochondrial DP or REF depth within one source record"
            )
        if len({r["alt"] for r in siblings}) != len(siblings):
            raise PortableEvidenceError(
                "one mitochondrial allele has multiple accepted rows from the same alignment"
            )
        dp, ref = counts.pop()
        alt_depths = [r["alt_depth"] for r in siblings]
        total = None
        if ref is not None and all(n is not None for n in alt_depths):
            counted = ref + sum(alt_depths)
            # A subset can omit sibling ALTs. A complete adapter partition also
            # resolves the denominator when DP includes unassigned reads.
            verified = frozenset(row["evidence_id"] for row in siblings) in complete_records
            if (verified or counted == dp) and (dp is None or counted <= dp):
                total = counted
        for row in siblings:
            result[row["evidence_id"]] = total
    return result


def _observation(
    row: dict[str, Any], allelic_depth: int | None, policy: MitochondrialReplicatePolicy
) -> AlleleFractionObservation | None:
    if allelic_depth is None or allelic_depth < policy.minimum_allelic_depth:
        return None
    alt = row["alt_depth"]
    lower, upper = _wilson_interval(alt, allelic_depth, policy.wilson_z)
    return AlleleFractionObservation(
        reported_fraction=float(row["allele_fraction"]),
        alt_depth=alt,
        allelic_depth=allelic_depth,
        lower=lower,
        upper=upper,
    )


def _pair_relation(
    first_observation: AlleleFractionObservation | None,
    second_observation: AlleleFractionObservation | None,
) -> str:
    if first_observation is None or second_observation is None:
        return "not_comparable"
    same_reported_fraction = math.isclose(
        first_observation.reported_fraction,
        second_observation.reported_fraction,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    same_observed_fraction = (
        first_observation.alt_depth * second_observation.allelic_depth
        == second_observation.alt_depth * first_observation.allelic_depth
    )
    if same_reported_fraction and same_observed_fraction:
        return "exact_fraction"
    intervals_overlap = max(first_observation.lower, second_observation.lower) <= min(
        first_observation.upper, second_observation.upper
    )
    return "sampling_interval_compatible" if intervals_overlap else "discordant"


def _biological_payload(row: dict[str, Any]) -> dict[str, Any]:
    return {field: row[field] for field in BIOLOGICAL_FIELDS}


def _biological_digest(rows: list[dict[str, Any]]) -> str:
    payload = sorted(
        (_biological_payload(row) for row in rows),
        key=lambda item: (
            item["chrom"],
            item["pos"],
            item["ref"],
            item["alt"],
            item["source_alignment_sha256"],
            item["source_record_digest"],
        ),
    )
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode()).hexdigest()


def _group_id(rows: list[dict[str, Any]]) -> str:
    payload = sorted(
        (_biological_payload(row) for row in rows),
        key=lambda item: (item["source_alignment_sha256"], item["source_record_digest"]),
    )
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return f"mt-concordance:{hashlib.sha256(rendered.encode()).hexdigest()}"


def annotate_mitochondrial_replicate_concordance(
    input_rows: list[dict[str, Any]],
    *,
    policy: MitochondrialReplicatePolicy = DEFAULT_MT_REPLICATE_POLICY,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
    complete_records: set[frozenset[str]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Annotate exact-allele replicate evidence without averaging or row removal."""

    if (
        policy.minimum_allelic_depth < 1
        or not math.isfinite(policy.wilson_z)
        or policy.wilson_z <= 0
    ):
        raise PortableEvidenceError("invalid mitochondrial replicate policy")
    rows = validate_rows(
        "mitochondrial_heteroplasmy",
        input_rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    if any(row["caller_name"] != "GATK Mutect2 mitochondrial mode" for row in rows):
        raise PortableEvidenceError(
            "mitochondrial replicate concordance requires GATK Mutect2 mitochondrial rows"
        )
    references = {str(row["reference_sha256"]) for row in rows}
    if len(references) != 1:
        raise PortableEvidenceError("mitochondrial replicate rows use different references")
    alignments = {str(row["source_alignment_sha256"]) for row in rows}
    if len(alignments) < 2:
        raise PortableEvidenceError(
            "mitochondrial replicate concordance requires two technical alignments"
        )

    biological_before = _biological_digest(rows)
    grouped: dict[tuple[str, int, str, str], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[(str(row["chrom"]), int(row["pos"]), str(row["ref"]), str(row["alt"]))].append(
            index
        )

    depths = _allelic_depths(rows, complete_records or set())
    annotated = [dict(row) for row in rows]
    relation_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    for _, indices in sorted(grouped.items()):
        component = [rows[index] for index in indices]
        component_alignments = [str(row["source_alignment_sha256"]) for row in component]
        if len(component_alignments) != len(set(component_alignments)):
            raise PortableEvidenceError(
                "one mitochondrial allele has multiple accepted rows from the same alignment"
            )
        caller_boundaries = {
            (
                str(row["caller_name"]),
                str(row["caller_version"]),
                str(row["resource_manifest_sha256"]),
            )
            for row in component
        }
        circular_audited = all(bool(row["circular_boundary_audited"]) for row in component)
        observations = [_observation(row, depths[row["evidence_id"]], policy) for row in component]
        pair_relations = [
            _pair_relation(observations[first], observations[second])
            for first in range(len(component))
            for second in range(first + 1, len(component))
        ]
        if len(component) == 1:
            status = "source_specific"
            relation = "single_source"
        elif len(caller_boundaries) != 1 or not circular_audited:
            status = "not_comparable"
            relation = "caller_resource_or_circular_boundary_mismatch"
        elif "discordant" in pair_relations:
            status = "comparable_discordant"
            relation = "allele_fraction_discordant"
        elif "not_comparable" in pair_relations:
            status = "not_comparable"
            relation = (
                "unresolved_allelic_depth"
                if any(depths[row["evidence_id"]] is None for row in component)
                else "insufficient_allelic_depth"
            )
        elif pair_relations and all(value == "exact_fraction" for value in pair_relations):
            status = "comparable_concordant"
            relation = "exact_fraction"
        else:
            status = "comparable_concordant"
            relation = "sampling_interval_compatible"

        group_id = _group_id(component)
        relation_counts[relation] += 1
        for index in indices:
            row = annotated[index]
            row["concordance_group_id"] = group_id
            row["concordance_status"] = status
            limitations = list(row["limitations"])
            filters = list(row["filter_reasons"])
            if status == "source_specific":
                limitations.append(
                    "no exact-allele row was supplied by another technical alignment; variant-only absence is not reference evidence"
                )
            elif status == "not_comparable":
                limitations.append(
                    "mitochondrial technical-replicate evidence lacks a common audited caller/resource boundary or sufficient finite allelic depth"
                )
                filters.append("REPLICATE:MT_NOT_COMPARABLE")
                row["confidence_tier"] = "REVIEW"
            elif status == "comparable_discordant":
                limitations.append(
                    "technical-replicate mitochondrial allele fractions have non-overlapping Wilson intervals"
                )
                filters.append("REPLICATE:MT_ALLELE_FRACTION_DISCORDANT")
                row["confidence_tier"] = "REVIEW"
            elif relation == "sampling_interval_compatible":
                limitations.append(
                    "technical-replicate mitochondrial allele fractions agree only within finite Wilson sampling intervals"
                )
            row["filter_reasons"] = sorted(set(filters))
            row["limitations"] = sorted(set(limitations))
            row["evidence_id"] = make_evidence_id("mitochondrial_heteroplasmy", row)
            status_counts[status] += 1

    validated = validate_rows(
        "mitochondrial_heteroplasmy",
        annotated,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_after = _biological_digest(validated)
    if biological_before != biological_after:
        raise PortableEvidenceError(
            "mitochondrial replicate concordance changed biological call fields"
        )
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "policy": {
            "minimum_allelic_depth": policy.minimum_allelic_depth,
            "denominator": MT_DEPTH_POLICY,
            "wilson_z": policy.wilson_z,
        },
        "input_rows": len(rows),
        "output_rows": len(validated),
        "technical_alignments": len(alignments),
        "exact_allele_groups": len(grouped),
        "pairwise_comparisons": sum(
            len(indices) * (len(indices) - 1) // 2 for indices in grouped.values()
        ),
        "group_relation_counts": dict(sorted(relation_counts.items())),
        "row_status_counts": dict(sorted(status_counts.items())),
        "source_rows_removed": 0,
        "allele_fraction_averaging_used": False,
        "majority_vote_used": False,
        "variant_only_absence_used_as_reference": False,
        "biological_fields_before_sha256": biological_before,
        "biological_fields_after_sha256": biological_after,
        "biological_fields_changed": False,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def complete_mt_partition(
    report: dict[str, Any], item: dict[str, Any], rows: list[dict[str, Any]]
) -> bool:
    """Match a verified output partition to the native parser's complete allele inventory."""
    parsing = report.get("parsing")
    if (
        report.get("schema_version") != 1
        or report.get("profile_id") != "gatk-mitochondria-normal-shifted-4.6.2.0"
        or report.get("status") != "PASS"
        or report.get("source_rows_removed") != 0
        or not isinstance(parsing, list)
        or len(parsing) != 1
        or not isinstance(parsing[0], dict)
        or parsing[0].get("domain") != "mitochondrial_heteroplasmy"
    ):
        return False
    summary = parsing[0].get("summary")
    source = (report.get("source_alignment_sha256"), report.get("raw_artifact_sha256"))
    return (
        isinstance(summary, dict)
        and summary.get("status") == "PASS"
        and summary.get("raw_artifact_sha256") == source[1]
        and summary.get("allele_rows") == len(rows) == item.get("rows")
        and summary.get("canonical_row_multiset") == item.get("canonical_row_multiset")
        and all(
            (row["source_alignment_sha256"], row["raw_artifact_sha256"]) == source for row in rows
        )
    )


def _complete_partition(path: Path, rows: list[dict[str, Any]]) -> bool:
    """Recognize the complete output of one native adapter, including legacy tables."""
    summary = path.parent / "adapter-summary.json"
    if not summary.is_file():
        return False
    try:
        report = json.loads(summary.read_bytes())
    except (ValueError, OSError):
        return False
    if not isinstance(report, dict) or not isinstance(report.get("files"), list):
        return False
    claimed = report.pop("adapter_identity_sha256", None)
    files = [
        item
        for item in report["files"]
        if isinstance(item, dict) and item.get("domain") == "mitochondrial_heteroplasmy"
    ]
    try:
        identity = canonical_json_sha256(report)
    except ValueError:
        return False
    if identity != claimed or len(files) != 1:
        return False
    item = files[0]
    return (
        item.get("path") == path.name
        and item.get("sha256") == sha256_file(path)
        and item.get("canonical_row_multiset") == canonical_row_multiset(rows)
        and complete_mt_partition(report, item, rows)
    )


def read_mitochondrial_parquets(
    paths: list[Path], *, complete_records: set[frozenset[str]] | None = None
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        table = pq.read_table(path)
        if table.schema != schema_for("mitochondrial_heteroplasmy"):
            raise PortableEvidenceError(
                "mitochondrial concordance input has the wrong physical schema"
            )
        partition = table.to_pylist()
        if complete_records is not None and partition and _complete_partition(path, partition):
            complete_records.update(verified_mt_records(partition))
        rows.extend(partition)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-parquet", type=Path, action="append", required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument(
        "--minimum-allelic-depth",
        type=int,
        default=DEFAULT_MT_REPLICATE_POLICY.minimum_allelic_depth,
    )
    parser.add_argument("--wilson-z", type=float, default=DEFAULT_MT_REPLICATE_POLICY.wilson_z)
    args = parser.parse_args()
    complete_records: set[frozenset[str]] = set()
    rows = read_mitochondrial_parquets(args.input_parquet, complete_records=complete_records)
    policy = MitochondrialReplicatePolicy(
        minimum_allelic_depth=args.minimum_allelic_depth,
        wilson_z=args.wilson_z,
    )
    annotated, summary = annotate_mitochondrial_replicate_concordance(
        rows,
        policy=policy,
        complete_records=complete_records,
        expected_reference_sha256=args.reference_sha256,
    )
    output = write_extension(
        "mitochondrial_heteroplasmy",
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
