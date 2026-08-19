"""Type-aware structural-event concordance without destructive intersection."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from .portable_model import (
    CHROM_RANK,
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    schema_for,
    validate_rows,
    write_extension,
)

ORIENTATION_ENDS = {
    "3to5": ("3", "5"),
    "5to3": ("5", "3"),
    "3to3": ("3", "3"),
    "5to5": ("5", "5"),
    "++": ("3", "5"),
    "-+": ("5", "3"),
    "+-": ("3", "3"),
    "--": ("5", "5"),
    "inversion": ("?", "?"),
    # DELLY uses CT=NtoN for insertions, whose breakend orientation is not
    # applicable.  Manta represents the same portable semantics explicitly as
    # ``not_applicable``.
    "NtoN": ("?", "?"),
    "not_applicable": ("?", "?"),
}


@dataclass(frozen=True)
class ConcordancePolicy:
    max_breakpoint_distance: int = 500
    min_reciprocal_overlap: float = 0.5
    min_size_similarity: float = 0.7


DEFAULT_CONCORDANCE_POLICY = ConcordancePolicy()


@dataclass(frozen=True)
class Breakpoint:
    chrom: str
    pos: int
    ci_lower: int
    ci_upper: int
    end_type: str

    @property
    def lower(self) -> int:
        return self.pos + self.ci_lower

    @property
    def upper(self) -> int:
        return self.pos + self.ci_upper


@dataclass
class Event:
    key: tuple[str, str, str, str, str]
    row_indices: list[int]
    sv_type: str
    first: Breakpoint
    second: Breakpoint
    length: int | None
    call_status: str
    genotype: tuple[int, ...] | None
    event_id_collision_partitioned: bool = False


def _canonical_breakpoints(row: dict[str, Any]) -> tuple[Breakpoint, Breakpoint]:
    try:
        end1, end2 = ORIENTATION_ENDS[str(row["orientation"])]
    except KeyError as exc:
        raise PortableEvidenceError("unsupported SV orientation for concordance") from exc
    first = Breakpoint(
        chrom=str(row["chrom1"]),
        pos=int(row["pos1"]),
        ci_lower=int(row["cipos"][0]),
        ci_upper=int(row["cipos"][1]),
        end_type=end1,
    )
    second = Breakpoint(
        chrom=str(row["chrom2"]),
        pos=int(row["pos2"]),
        ci_lower=int(row["ciend"][0]),
        ci_upper=int(row["ciend"][1]),
        end_type=end2,
    )
    if first.lower > first.upper or second.lower > second.upper:
        raise PortableEvidenceError("SV breakpoint confidence interval bounds are reversed")
    first_key = (CHROM_RANK[first.chrom], first.pos, first.end_type)
    second_key = (CHROM_RANK[second.chrom], second.pos, second.end_type)
    return (first, second) if first_key <= second_key else (second, first)


def _event_signature(event: Event) -> dict[str, Any]:
    return {
        "sv_type": event.sv_type,
        "first": event.first.__dict__,
        "second": event.second.__dict__,
        "length": event.length,
    }


def _build_events(rows: list[dict[str, Any]]) -> list[Event]:
    grouped: dict[tuple[str, str, str, str], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        source_key = (
            str(row["caller_name"]),
            str(row["technical_source_role"]),
            str(row["source_alignment_sha256"]),
            str(row["event_id"]),
        )
        grouped[source_key].append(index)
    events: list[Event] = []
    for source_key, indices in sorted(grouped.items()):
        partitions: dict[tuple[str, Breakpoint, Breakpoint], list[int]] = defaultdict(list)
        record_signatures: dict[str, set[tuple[str, Breakpoint, Breakpoint]]] = defaultdict(set)
        for index in indices:
            row = rows[index]
            first, second = _canonical_breakpoints(row)
            signature = (str(row["sv_type"]), first, second)
            partitions[signature].append(index)
            record_signatures[str(row["caller_record_id"])].add(signature)

        collision_partitioned = len(partitions) > 1
        if collision_partitioned:
            signatures = list(partitions)
            point_signatures = {
                (
                    sv_type,
                    first.chrom,
                    first.pos,
                    first.end_type,
                    second.chrom,
                    second.pos,
                    second.end_type,
                )
                for sv_type, first, second in signatures
            }
            confidence_intervals_overlap = all(
                first_a.lower <= first_b.upper
                and first_b.lower <= first_a.upper
                and second_a.lower <= second_b.upper
                and second_b.lower <= second_a.upper
                for offset, (_, first_a, second_a) in enumerate(signatures)
                for _, first_b, second_b in signatures[offset + 1 :]
            )
            if (
                len(point_signatures) != 1
                or not confidence_intervals_overlap
                or any(
                    len(signatures_for_record) != 1
                    for signatures_for_record in record_signatures.values()
                )
            ):
                raise PortableEvidenceError("one caller event has inconsistent breakend records")

        ordered_partitions = sorted(
            partitions.items(),
            key=lambda item: (
                item[0][0],
                CHROM_RANK[item[0][1].chrom],
                item[0][1].pos,
                item[0][1].ci_lower,
                item[0][1].ci_upper,
                item[0][1].end_type,
                CHROM_RANK[item[0][2].chrom],
                item[0][2].pos,
                item[0][2].ci_lower,
                item[0][2].ci_upper,
                item[0][2].end_type,
            ),
        )
        for signature, partition_indices in ordered_partitions:
            sv_type, first, second = signature
            statuses = {str(rows[index]["call_status"]) for index in partition_indices}
            genotypes = {
                tuple(int(value) for value in rows[index]["genotype"])
                if rows[index]["genotype"] is not None
                else None
                for index in partition_indices
            }
            if len(statuses) != 1 or len(genotypes) != 1:
                raise PortableEvidenceError(
                    "one caller event has internally inconsistent genotypes"
                )
            representative = rows[partition_indices[0]]
            raw_length = representative.get("length")
            partition_payload = json.dumps(
                {
                    "sv_type": sv_type,
                    "first": first.__dict__,
                    "second": second.__dict__,
                },
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
            )
            partition_id = hashlib.sha256(partition_payload.encode()).hexdigest()
            events.append(
                Event(
                    key=(*source_key, partition_id),
                    row_indices=partition_indices,
                    sv_type=sv_type,
                    first=first,
                    second=second,
                    length=abs(int(raw_length)) if raw_length is not None else None,
                    call_status=next(iter(statuses)),
                    genotype=next(iter(genotypes)),
                    event_id_collision_partitioned=collision_partitioned,
                )
            )
    return events


def _interval_distance(first: Breakpoint, second: Breakpoint) -> int | None:
    if first.chrom != second.chrom or first.end_type != second.end_type:
        return None
    if first.upper < second.lower:
        return second.lower - first.upper
    if second.upper < first.lower:
        return first.lower - second.upper
    return 0


def _size_similarity(first: int | None, second: int | None) -> float | None:
    if first is None or second is None:
        return None
    if first == second == 0:
        return 1.0
    if first == 0 or second == 0:
        return 0.0
    return min(first, second) / max(first, second)


def _reciprocal_overlap(first: Event, second: Event) -> float:
    first_start, first_end = sorted((first.first.pos, first.second.pos))
    second_start, second_end = sorted((second.first.pos, second.second.pos))
    overlap = max(0, min(first_end, second_end) - max(first_start, second_start) + 1)
    first_span = max(1, first_end - first_start + 1)
    second_span = max(1, second_end - second_start + 1)
    return min(overlap / first_span, overlap / second_span)


def _events_match(first: Event, second: Event, policy: ConcordancePolicy) -> bool:
    if first.sv_type != second.sv_type:
        return False
    first_distance = _interval_distance(first.first, second.first)
    second_distance = _interval_distance(first.second, second.second)
    if first_distance is None or second_distance is None:
        return False
    if (
        first_distance > policy.max_breakpoint_distance
        or second_distance > policy.max_breakpoint_distance
    ):
        return False
    similarity = _size_similarity(first.length, second.length)
    if similarity is not None and similarity < policy.min_size_similarity:
        return False
    if first.sv_type in {"DEL", "DUP", "INV", "CNV"}:
        return _reciprocal_overlap(first, second) >= policy.min_reciprocal_overlap
    return True


def _source(event: Event) -> tuple[str, str, str]:
    return event.key[:3]


def _matching_event_pairs(events: list[Event], policy: ConcordancePolicy) -> list[tuple[int, int]]:
    """Return the exact cross-source matches without an all-events quadratic scan."""
    partitions: dict[tuple[str, str, str, str, str], list[int]] = defaultdict(list)
    for index, event in enumerate(events):
        partitions[
            (
                event.sv_type,
                event.first.chrom,
                event.first.end_type,
                event.second.chrom,
                event.second.end_type,
            )
        ].append(index)

    matches: list[tuple[int, int]] = []
    for partition in sorted(partitions):
        ordered = sorted(
            partitions[partition],
            key=lambda index: (
                events[index].first.lower,
                events[index].first.upper,
                events[index].second.lower,
                events[index].second.upper,
                events[index].key,
            ),
        )
        for offset, first_index in enumerate(ordered):
            first = events[first_index]
            first_window_end = first.first.upper + policy.max_breakpoint_distance
            for second_index in ordered[offset + 1 :]:
                second = events[second_index]
                if second.first.lower > first_window_end:
                    break
                if _source(first) == _source(second):
                    continue
                if _events_match(first, second, policy):
                    matches.append(tuple(sorted((first_index, second_index))))
    return sorted(matches)


def _group_id(events: list[Event]) -> str:
    payload = [_event_signature(event) for event in sorted(events, key=lambda item: item.key)]
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return f"sv-concordance:{hashlib.sha256(rendered.encode()).hexdigest()}"


def _biological_digest(rows: list[dict[str, Any]]) -> str:
    fields = (
        "caller_record_id",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "sv_type",
        "orientation",
        "length",
        "genotype",
        "call_status",
        "copy_number",
        "source_record_digest",
    )
    payload = sorted(
        json.dumps(
            {field: row[field] for field in fields},
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        for row in rows
    )
    rendered = json.dumps(payload, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode()).hexdigest()


def annotate_sv_concordance(
    input_rows: list[dict[str, Any]],
    *,
    policy: ConcordancePolicy = DEFAULT_CONCORDANCE_POLICY,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows = validate_rows(
        "structural_variants",
        input_rows,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_before = _biological_digest(rows)
    events = _build_events(rows)
    sources = {_source(event) for event in events}
    if len(sources) < 2:
        raise PortableEvidenceError("SV concordance requires at least two independent sources")
    adjacency: dict[int, set[int]] = {index: set() for index in range(len(events))}
    for first_index, second_index in _matching_event_pairs(events, policy):
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
    for component_indices in components:
        component = [events[index] for index in component_indices]
        has_edges = any(adjacency[index] for index in component_indices)
        source_counts = Counter(_source(event) for event in component)
        statuses = {event.call_status for event in component}
        genotypes = {event.genotype for event in component}
        if not has_edges:
            status = (
                "source_specific"
                if component[0].call_status == "non_reference"
                else "not_comparable"
            )
        elif any(count > 1 for count in source_counts.values()):
            status = "comparable_discordant"
        elif "no_call" in statuses:
            status = "not_comparable"
        elif len(statuses) == 1 and len(genotypes) == 1:
            status = "comparable_concordant"
        else:
            status = "comparable_discordant"
        group_id = _group_id(component)
        for event in component:
            for row_index in event.row_indices:
                row = annotated[row_index]
                row["concordance_group_id"] = group_id
                row["concordance_status"] = status
                limitations = list(row["limitations"])
                if event.event_id_collision_partitioned:
                    limitations.append(
                        "overlapping caller confidence-interval representations shared an event identifier and were retained as distinct source records"
                    )
                    row["confidence_tier"] = "REVIEW"
                if status == "source_specific":
                    limitations.append(
                        "no comparable event was found in another supplied source; absence is not refutation"
                    )
                elif status == "not_comparable":
                    limitations.append("cross-source structural evidence is not comparable")
                elif status == "comparable_discordant":
                    limitations.append(
                        "comparable structural evidence is ambiguous or genotype-discordant"
                    )
                    row["confidence_tier"] = "REVIEW"
                row["limitations"] = sorted(set(limitations))
                row["evidence_id"] = make_evidence_id("structural_variants", row)
                status_counts[status] += 1
    validated = validate_rows(
        "structural_variants",
        annotated,
        expected_reference_sha256=expected_reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    biological_after = _biological_digest(validated)
    if biological_before != biological_after:
        raise PortableEvidenceError("SV concordance changed biological call fields")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "input_rows": len(rows),
        "output_rows": len(validated),
        "events": len(events),
        "sources": len(sources),
        "event_id_collision_groups_partitioned": len(
            {event.key[:4] for event in events if event.event_id_collision_partitioned}
        ),
        "event_id_collision_rows_partitioned": sum(
            len(event.row_indices) for event in events if event.event_id_collision_partitioned
        ),
        "candidate_edges": sum(len(value) for value in adjacency.values()) // 2,
        "components": len(components),
        "row_status_counts": dict(sorted(status_counts.items())),
        "biological_fields_before_sha256": biological_before,
        "biological_fields_after_sha256": biological_after,
        "biological_fields_changed": False,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "private_paths_persisted": False,
    }
    return validated, summary


def read_sv_parquets(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        table = pq.read_table(path)
        if table.schema != schema_for("structural_variants"):
            raise PortableEvidenceError("SV concordance input has the wrong physical schema")
        rows.extend(table.to_pylist())
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-parquet", type=Path, action="append", required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--max-breakpoint-distance", type=int, default=500)
    parser.add_argument("--min-reciprocal-overlap", type=float, default=0.5)
    parser.add_argument("--min-size-similarity", type=float, default=0.7)
    args = parser.parse_args()
    policy = ConcordancePolicy(
        max_breakpoint_distance=args.max_breakpoint_distance,
        min_reciprocal_overlap=args.min_reciprocal_overlap,
        min_size_similarity=args.min_size_similarity,
    )
    rows = read_sv_parquets(args.input_parquet)
    annotated, summary = annotate_sv_concordance(
        rows, policy=policy, expected_reference_sha256=args.reference_sha256
    )
    output = write_extension(
        "structural_variants",
        annotated,
        args.output_parquet,
        expected_reference_sha256=args.reference_sha256,
    )
    print(
        json.dumps(
            {"status": "PASS", "concordance": summary, "output": output}, indent=2, sort_keys=True
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
