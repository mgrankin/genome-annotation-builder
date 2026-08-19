"""Standalone, path-free raw capture for bounded read-backed evidence.

This module intentionally depends only on the Python standard library and
``pysam`` so that the exact source file can run in a small digest-pinned OCI
image.  It writes no genotype, sample identifier, read name, sequence, or local
path.  The portable adapter validates this raw artifact before creating typed
Parquet evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import statistics
import sys
import tempfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pysam

SCHEMA_VERSION = 1
TARGET_COLUMNS = [
    "target_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "candidate_category",
    "local_complexity_status",
    "local_haplotype_status",
]
STANDARD_CHROMS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY", "chrM"])
CONCRETE_ALLELE = frozenset("ACGT")
CIGAR_SOFT_CLIP = 4
CIGAR_INSERTION = 1
CIGAR_DELETION = 2
CIGAR_REFERENCE_CONSUMING = frozenset({0, 2, 3, 7, 8})


class CaptureError(ValueError):
    """The requested raw read-back capture violates its public contract."""


@dataclass(frozen=True)
class CapturePolicy:
    minimum_mapping_quality: int = 20
    minimum_base_quality: int = 20
    minimum_alt_fragments: int = 3
    nearby_indel_window: int = 10
    include_duplicates: bool = False
    max_depth: int = 1_000_000


DEFAULT_CAPTURE_POLICY = CapturePolicy()


@dataclass(frozen=True)
class Target:
    target_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    candidate_category: str
    local_complexity_status: str
    local_haplotype_status: str


@dataclass(frozen=True)
class Observation:
    allele_class: str
    mapping_quality: int
    base_quality: float | None
    read_position: float | None
    reverse: bool
    orientation: str
    soft_clipped: bool
    nearby_non_target_indel: bool


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode()


def _canonical_chrom(value: str) -> str:
    token = value.strip()
    if token == "MT":
        token = "M"
    if not token.startswith("chr"):
        token = f"chr{token}"
    if token not in STANDARD_CHROMS:
        raise CaptureError("target_chromosome_invalid")
    return token


def _alignment_contig(chrom: str, references: tuple[str, ...]) -> str:
    candidates = [chrom, chrom.removeprefix("chr")]
    if chrom == "chrM":
        candidates.extend(["MT", "M"])
    for candidate in candidates:
        if candidate in references:
            return candidate
    raise CaptureError("target_chromosome_absent")


def _variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(ref) == 1 and alt.startswith(ref):
        return "insertion"
    if len(alt) == 1 and ref.startswith(alt):
        return "deletion"
    return "complex"


def _read_targets(path: Path) -> tuple[list[Target], str]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != TARGET_COLUMNS:
                raise CaptureError("target_schema_invalid")
            raw_rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise CaptureError("target_manifest_unreadable") from exc
    if not raw_rows:
        raise CaptureError("target_manifest_empty")
    targets: list[Target] = []
    seen_ids: set[str] = set()
    seen_alleles: set[tuple[str, int, str, str]] = set()
    for source in raw_rows:
        if None in source:
            raise CaptureError("target_row_width_invalid")
        target_id = source["target_id"].strip()
        if not target_id or "/" in target_id or "\\" in target_id or target_id in seen_ids:
            raise CaptureError("target_id_invalid_or_duplicated")
        seen_ids.add(target_id)
        chrom = _canonical_chrom(source["chrom"])
        try:
            pos = int(source["pos"])
        except ValueError as exc:
            raise CaptureError("target_position_invalid") from exc
        ref = source["ref"].strip().upper()
        alt = source["alt"].strip().upper()
        if (
            pos < 1
            or not ref
            or not alt
            or set(ref) - CONCRETE_ALLELE
            or set(alt) - CONCRETE_ALLELE
            or ref == alt
        ):
            raise CaptureError("target_allele_invalid")
        allele_key = (chrom, pos, ref, alt)
        if allele_key in seen_alleles:
            raise CaptureError("target_allele_duplicated")
        seen_alleles.add(allele_key)
        annotations = [
            source["candidate_category"].strip(),
            source["local_complexity_status"].strip(),
            source["local_haplotype_status"].strip(),
        ]
        if not all(annotations):
            raise CaptureError("target_annotation_incomplete")
        targets.append(
            Target(
                target_id=target_id,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                candidate_category=annotations[0],
                local_complexity_status=annotations[1],
                local_haplotype_status=annotations[2],
            )
        )
    return targets, _sha256_file(path)


def _query_slice(read: pysam.AlignedSegment, start: int, length: int) -> str | None:
    sequence = read.query_sequence
    if sequence is None or start < 0 or start + length > len(sequence):
        return None
    return sequence[start : start + length].upper()


def _base_quality(read: pysam.AlignedSegment, start: int, length: int) -> float | None:
    qualities = read.query_qualities
    if qualities is None or start < 0 or start + length > len(qualities):
        return None
    values = qualities[start : start + length]
    return float(min(values)) if values else None


def _classify(
    pileup_read: pysam.PileupRead, target: Target, variant_type: str
) -> tuple[str, float | None, float | None]:
    read = pileup_read.alignment
    query_position = pileup_read.query_position
    if query_position is None or pileup_read.is_refskip or pileup_read.is_del:
        return "other", None, None
    if variant_type == "SNV":
        observed = _query_slice(read, query_position, 1)
        allele_class = "ref" if observed == target.ref else "alt" if observed == target.alt else "other"
        span = 1
    elif variant_type == "MNV":
        observed = _query_slice(read, query_position, len(target.ref))
        allele_class = (
            "other"
            if pileup_read.indel != 0
            else "ref"
            if observed == target.ref
            else "alt"
            if observed == target.alt
            else "other"
        )
        span = len(target.ref)
    elif variant_type == "insertion":
        length = len(target.alt) - 1
        observed_alt = _query_slice(read, query_position, len(target.alt))
        observed_anchor = _query_slice(read, query_position, 1)
        if pileup_read.indel == length and observed_alt == target.alt:
            allele_class, span = "alt", len(target.alt)
        elif pileup_read.indel == 0 and observed_anchor == target.ref:
            allele_class, span = "ref", 1
        else:
            allele_class, span = "other", 1
    elif variant_type == "deletion":
        length = len(target.ref) - 1
        observed_anchor = _query_slice(read, query_position, 1)
        observed_ref = _query_slice(read, query_position, len(target.ref))
        if pileup_read.indel == -length and observed_anchor == target.alt:
            allele_class, span = "alt", 1
        elif pileup_read.indel == 0 and observed_ref == target.ref:
            allele_class, span = "ref", len(target.ref)
        else:
            allele_class, span = "other", 1
    else:
        return "other", None, None
    return allele_class, _base_quality(read, query_position, span), float(query_position + 1)


def _orientation(read: pysam.AlignedSegment) -> str:
    if (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse):
        return "f1r2"
    if (read.is_read2 and not read.is_reverse) or (read.is_read1 and read.is_reverse):
        return "f2r1"
    return "unclassified_orientation"


def _has_nearby_non_target_indel(
    read: pysam.AlignedSegment,
    target: Target,
    variant_type: str,
    window: int,
) -> bool:
    reference_cursor = int(read.reference_start)
    target_zero = target.pos - 1
    for operation, length in read.cigartuples or ():
        if operation == CIGAR_INSERTION:
            anchor = reference_cursor - 1
            is_target = variant_type == "insertion" and anchor == target_zero and length == len(target.alt) - 1
            if not is_target and abs(anchor - target_zero) <= window:
                return True
        elif operation == CIGAR_DELETION:
            deletion_start = reference_cursor
            deletion_end = reference_cursor + length - 1
            is_target = (
                variant_type == "deletion"
                and deletion_start == target_zero + 1
                and length == len(target.ref) - 1
            )
            distance = min(abs(deletion_start - target_zero), abs(deletion_end - target_zero))
            if not is_target and distance <= window:
                return True
        if operation in CIGAR_REFERENCE_CONSUMING:
            reference_cursor += length
    return False


def _collapse_fragment(observations: list[Observation]) -> tuple[Observation, bool]:
    if len(observations) == 1:
        return observations[0], False
    if len({item.allele_class for item in observations}) != 1:
        representative = max(observations, key=lambda item: (item.mapping_quality, item.base_quality or -1))
        return (
            Observation(
                allele_class="other",
                mapping_quality=representative.mapping_quality,
                base_quality=None,
                read_position=None,
                reverse=representative.reverse,
                orientation="unclassified",
                soft_clipped=any(item.soft_clipped for item in observations),
                nearby_non_target_indel=any(item.nearby_non_target_indel for item in observations),
            ),
            True,
        )
    return max(observations, key=lambda item: (item.base_quality or -1, item.mapping_quality)), False


def _empty_metrics(policy: CapturePolicy) -> dict[str, Any]:
    metrics: dict[str, Any] = {
        "eligible_fragments": 0,
        "ref_reads": 0,
        "alt_reads": 0,
        "other_reads": 0,
        "nearby_indel_window": policy.nearby_indel_window,
        "excluded_low_mapping_quality_reads": 0,
        "excluded_low_base_quality_reads": 0,
        "excluded_duplicate_reads": 0,
        "excluded_nonprimary_reads": 0,
        "excluded_qc_fail_reads": 0,
        "collapsed_overlapping_mates": 0,
        "discordant_overlapping_mates": 0,
        "minimum_mapping_quality": policy.minimum_mapping_quality,
        "minimum_base_quality": policy.minimum_base_quality,
        "duplicate_reads_included": policy.include_duplicates,
        "overlapping_mate_policy": "collapse_by_query_name",
        "nearby_non_target_indel_reads": 0,
    }
    for strand in ("forward", "reverse"):
        for allele in ("ref", "alt"):
            metrics[f"{strand}_{allele}_reads"] = 0
    for orientation in ("f1r2", "f2r1", "unclassified_orientation"):
        for allele in ("ref", "alt"):
            metrics[f"{orientation}_{allele}_reads"] = 0
    for allele in ("ref", "alt"):
        metrics[f"mean_{allele}_mapping_quality"] = None
        metrics[f"mean_{allele}_base_quality"] = None
        metrics[f"mean_{allele}_read_position"] = None
        metrics[f"soft_clipped_{allele}_reads"] = 0
    return metrics


def _mean(values: list[float]) -> float | None:
    return float(statistics.fmean(values)) if values else None


def _observe_target(
    alignment: pysam.AlignmentFile,
    alignment_contig: str,
    target: Target,
    variant_type: str,
    policy: CapturePolicy,
) -> dict[str, Any]:
    metrics = _empty_metrics(policy)
    by_fragment: dict[str, list[Observation]] = defaultdict(list)
    for column in alignment.pileup(
        alignment_contig,
        target.pos - 1,
        target.pos,
        truncate=True,
        stepper="nofilter",
        min_base_quality=0,
        min_mapping_quality=0,
        max_depth=policy.max_depth,
        ignore_overlaps=False,
        ignore_orphans=False,
    ):
        if column.reference_pos != target.pos - 1:
            continue
        for pileup_read in column.pileups:
            read = pileup_read.alignment
            if read.is_secondary or read.is_supplementary:
                metrics["excluded_nonprimary_reads"] += 1
                continue
            if read.is_qcfail:
                metrics["excluded_qc_fail_reads"] += 1
                continue
            if read.is_duplicate and not policy.include_duplicates:
                metrics["excluded_duplicate_reads"] += 1
                continue
            if read.mapping_quality < policy.minimum_mapping_quality:
                metrics["excluded_low_mapping_quality_reads"] += 1
                continue
            allele_class, base_quality, read_position = _classify(pileup_read, target, variant_type)
            if base_quality is not None and base_quality < policy.minimum_base_quality:
                metrics["excluded_low_base_quality_reads"] += 1
                continue
            by_fragment[read.query_name].append(
                Observation(
                    allele_class=allele_class,
                    mapping_quality=int(read.mapping_quality),
                    base_quality=base_quality,
                    read_position=read_position,
                    reverse=bool(read.is_reverse),
                    orientation=_orientation(read),
                    soft_clipped=any(op == CIGAR_SOFT_CLIP for op, _length in read.cigartuples or ()),
                    nearby_non_target_indel=_has_nearby_non_target_indel(
                        read, target, variant_type, policy.nearby_indel_window
                    ),
                )
            )
    collapsed: list[Observation] = []
    for observations in by_fragment.values():
        observation, discordant = _collapse_fragment(observations)
        metrics["collapsed_overlapping_mates"] += len(observations) - 1
        metrics["discordant_overlapping_mates"] += int(discordant)
        collapsed.append(observation)
    metrics["eligible_fragments"] = len(collapsed)
    values: dict[str, list[float]] = {
        f"{allele}_{metric}": []
        for allele in ("ref", "alt")
        for metric in ("mapping", "base", "position")
    }
    for observation in collapsed:
        allele = observation.allele_class
        metrics[f"{allele}_reads"] += 1
        if allele not in {"ref", "alt"}:
            continue
        strand = "reverse" if observation.reverse else "forward"
        metrics[f"{strand}_{allele}_reads"] += 1
        metrics[f"{observation.orientation}_{allele}_reads"] += 1
        values[f"{allele}_mapping"].append(observation.mapping_quality)
        if observation.base_quality is not None:
            values[f"{allele}_base"].append(observation.base_quality)
        if observation.read_position is not None:
            values[f"{allele}_position"].append(observation.read_position)
        metrics[f"soft_clipped_{allele}_reads"] += int(observation.soft_clipped)
        metrics["nearby_non_target_indel_reads"] += int(observation.nearby_non_target_indel)
    for allele in ("ref", "alt"):
        metrics[f"mean_{allele}_mapping_quality"] = _mean(values[f"{allele}_mapping"])
        metrics[f"mean_{allele}_base_quality"] = _mean(values[f"{allele}_base"])
        metrics[f"mean_{allele}_read_position"] = _mean(values[f"{allele}_position"])
    return metrics


def _assessment(
    variant_type: str, metrics: dict[str, Any], policy: CapturePolicy
) -> tuple[str, list[str]]:
    if variant_type == "complex":
        return "unsupported_allele_representation", [
            "COMPLEX_ALLELE_REPRESENTATION_NOT_IMPLEMENTED"
        ]
    if metrics["eligible_fragments"] == 0:
        return "no_eligible_reads", ["NO_ELIGIBLE_READS"]
    reasons: list[str] = []
    if metrics["alt_reads"] == 0:
        reasons.append("NO_ALT_SUPPORT")
    elif metrics["alt_reads"] < policy.minimum_alt_fragments:
        reasons.append("LOW_ALT_SUPPORT")
    if metrics["alt_reads"] > 0 and (
        metrics["forward_alt_reads"] == 0 or metrics["reverse_alt_reads"] == 0
    ):
        reasons.append("ALT_ONE_STRAND_ONLY")
    if metrics["nearby_non_target_indel_reads"]:
        reasons.append("NEARBY_NON_TARGET_INDEL")
    if metrics["discordant_overlapping_mates"]:
        reasons.append("DISCORDANT_OVERLAPPING_MATES")
    return "assessed", sorted(reasons)


def build_capture(
    alignment_path: str | Path,
    reference_fasta: str | Path,
    target_manifest: str | Path,
    *,
    policy: CapturePolicy = DEFAULT_CAPTURE_POLICY,
) -> dict[str, Any]:
    """Compute a deterministic, identifier-free raw evidence payload."""

    if (
        policy.minimum_mapping_quality < 0
        or policy.minimum_base_quality < 0
        or policy.minimum_alt_fragments < 1
        or policy.nearby_indel_window < 0
        or policy.max_depth < 1
    ):
        raise CaptureError("capture_policy_invalid")
    targets, target_sha256 = _read_targets(Path(target_manifest))
    try:
        reference = pysam.FastaFile(str(reference_fasta))
        alignment = pysam.AlignmentFile(
            str(alignment_path), reference_filename=str(reference_fasta)
        )
    except (OSError, ValueError) as exc:
        raise CaptureError("alignment_or_reference_unreadable") from exc
    records: list[dict[str, Any]] = []
    sample_bindings = {
        str(record.get("SM"))
        for record in alignment.header.to_dict().get("RG", [])
        if record.get("SM")
    }
    try:
        if not alignment.has_index():
            raise CaptureError("alignment_index_unusable")
        for target in targets:
            alignment_contig = _alignment_contig(target.chrom, alignment.references)
            reference_contig = _alignment_contig(target.chrom, reference.references)
            observed_ref = reference.fetch(
                reference_contig, target.pos - 1, target.pos - 1 + len(target.ref)
            ).upper()
            if observed_ref != target.ref:
                raise CaptureError("target_reference_mismatch")
            variant_type = _variant_type(target.ref, target.alt)
            metrics = (
                _empty_metrics(policy)
                if variant_type == "complex"
                else _observe_target(alignment, alignment_contig, target, variant_type, policy)
            )
            assessment_status, filter_reasons = _assessment(variant_type, metrics, policy)
            records.append(
                {
                    "target": asdict(target),
                    "variant_type": variant_type,
                    "assessment_status": assessment_status,
                    "metrics": metrics,
                    "filter_reasons": filter_reasons,
                }
            )
    finally:
        alignment.close()
        reference.close()
    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "PASS",
        "assembly": "GRCh38",
        "target_manifest_sha256": target_sha256,
        "policy": asdict(policy),
        "records": records,
        "record_count": len(records),
        "sample_identifier_persisted": False,
        "read_names_persisted": False,
        "read_sequences_persisted": False,
        "private_paths_persisted": False,
    }
    rendered = _canonical_bytes(payload)
    if any(binding.encode() in rendered for binding in sample_bindings):
        raise CaptureError("sample_identifier_persisted")
    return payload


def write_capture(output_path: str | Path, payload: dict[str, Any]) -> dict[str, Any]:
    """Publish canonical raw JSON atomically without replacing an existing file."""

    output = Path(output_path)
    if output.exists():
        raise CaptureError("output_already_exists")
    output.parent.mkdir(parents=True, exist_ok=True)
    data = _canonical_bytes(payload)
    descriptor, raw_name = tempfile.mkstemp(prefix=f".{output.name}.partial.", dir=output.parent)
    temporary = Path(raw_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, output)
        except FileExistsError as exc:
            raise CaptureError("output_already_exists") from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "PASS",
        "records": payload["record_count"],
        "raw_artifact_sha256": hashlib.sha256(data).hexdigest(),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-mapping-quality", type=int, default=20)
    parser.add_argument("--minimum-base-quality", type=int, default=20)
    parser.add_argument("--minimum-alt-fragments", type=int, default=3)
    parser.add_argument("--nearby-indel-window", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=1_000_000)
    parser.add_argument("--include-duplicates", action="store_true")
    args = parser.parse_args()
    try:
        policy = CapturePolicy(
            minimum_mapping_quality=args.minimum_mapping_quality,
            minimum_base_quality=args.minimum_base_quality,
            minimum_alt_fragments=args.minimum_alt_fragments,
            nearby_indel_window=args.nearby_indel_window,
            include_duplicates=args.include_duplicates,
            max_depth=args.max_depth,
        )
        result = write_capture(
            args.output,
            build_capture(
                args.alignment,
                args.reference_fasta,
                args.target_manifest,
                policy=policy,
            ),
        )
        print(json.dumps(result, sort_keys=True))
    except Exception:
        print("read-back capture failed; inspect the authorized private process log", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
