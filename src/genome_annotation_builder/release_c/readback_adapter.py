"""Bounded, exact-allele local read-back adapter for prioritized WGS targets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from .portable_model import (
    CHROM_RANK,
    CONCRETE_ALLELE,
    SHA256,
    SOURCE_ROLES,
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    sha256_file,
    validate_rows,
    write_extension,
)
from .readback_capture import TARGET_COLUMNS
from .readback_capture_v2 import EMPTY_TARGET_SHA256

CIGAR_SOFT_CLIP = 4
CIGAR_INSERTION = 1
CIGAR_DELETION = 2
CIGAR_REFERENCE_CONSUMING = frozenset({0, 2, 3, 7, 8})


@dataclass(frozen=True)
class ReadbackPolicy:
    minimum_mapping_quality: int = 20
    minimum_base_quality: int = 20
    minimum_alt_fragments: int = 3
    nearby_indel_window: int = 10
    include_duplicates: bool = False
    max_depth: int = 1_000_000


DEFAULT_READBACK_POLICY = ReadbackPolicy()


@dataclass(frozen=True)
class ReadbackProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    caller_name: str = "GAB bounded local pileup"
    caller_version: str = "1"


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


def _canonical_chrom(value: str) -> str:
    token = value.strip()
    if token == "MT":
        token = "M"
    if not token.startswith("chr"):
        token = f"chr{token}"
    if token not in CHROM_RANK:
        raise PortableEvidenceError("read-back target uses a non-canonical chromosome")
    return token


def _alignment_contig(chrom: str, references: tuple[str, ...]) -> str:
    candidates = [chrom, chrom.removeprefix("chr")]
    if chrom == "chrM":
        candidates.extend(["MT", "M"])
    for candidate in candidates:
        if candidate in references:
            return candidate
    raise PortableEvidenceError("read-back target chromosome is absent from the alignment")


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


def _read_targets(path: str | Path) -> tuple[list[Target], str]:
    target_sha256 = sha256_file(Path(path))
    try:
        with Path(path).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != TARGET_COLUMNS:
                raise PortableEvidenceError("read-back target manifest column schema differs")
            raw_rows = list(reader)
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("cannot read the read-back target manifest") from exc
    except csv.Error as exc:
        raise PortableEvidenceError("malformed read-back target manifest") from exc
    if not raw_rows and target_sha256 != EMPTY_TARGET_SHA256:
        raise PortableEvidenceError("empty read-back target manifest identity differs")
    targets: list[Target] = []
    seen_ids: set[str] = set()
    seen_alleles: set[tuple[str, int, str, str]] = set()
    for source in raw_rows:
        if None in source:
            raise PortableEvidenceError("read-back target row width differs")
        target_id = source["target_id"].strip()
        if not target_id or "/" in target_id or "\\" in target_id:
            raise PortableEvidenceError("read-back target ID is empty or path-like")
        if target_id in seen_ids:
            raise PortableEvidenceError("duplicate read-back target ID")
        seen_ids.add(target_id)
        chrom = _canonical_chrom(source["chrom"])
        try:
            pos = int(source["pos"])
        except ValueError as exc:
            raise PortableEvidenceError("invalid read-back target position") from exc
        ref = source["ref"].strip().upper()
        alt = source["alt"].strip().upper()
        if pos < 1 or CONCRETE_ALLELE.fullmatch(ref) is None:
            raise PortableEvidenceError("invalid read-back target REF")
        if CONCRETE_ALLELE.fullmatch(alt) is None or ref == alt:
            raise PortableEvidenceError("invalid read-back target ALT")
        allele_key = (chrom, pos, ref, alt)
        if allele_key in seen_alleles:
            raise PortableEvidenceError("duplicate read-back target allele")
        seen_alleles.add(allele_key)
        candidate_category = source["candidate_category"].strip()
        local_complexity = source["local_complexity_status"].strip()
        local_haplotype = source["local_haplotype_status"].strip()
        if not candidate_category or not local_complexity or not local_haplotype:
            raise PortableEvidenceError("read-back target annotations are incomplete")
        targets.append(
            Target(
                target_id=target_id,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                candidate_category=candidate_category,
                local_complexity_status=local_complexity,
                local_haplotype_status=local_haplotype,
            )
        )
    return targets, target_sha256


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
        allele_class = (
            "ref" if observed == target.ref else "alt" if observed == target.alt else "other"
        )
        span = 1
    elif variant_type == "MNV":
        observed = _query_slice(read, query_position, len(target.ref))
        if pileup_read.indel != 0:
            allele_class = "other"
        else:
            allele_class = (
                "ref" if observed == target.ref else "alt" if observed == target.alt else "other"
            )
        span = len(target.ref)
    elif variant_type == "insertion":
        insertion_length = len(target.alt) - 1
        observed_alt = _query_slice(read, query_position, len(target.alt))
        observed_anchor = _query_slice(read, query_position, 1)
        if pileup_read.indel == insertion_length and observed_alt == target.alt:
            allele_class = "alt"
            span = len(target.alt)
        elif pileup_read.indel == 0 and observed_anchor == target.ref:
            allele_class = "ref"
            span = 1
        else:
            allele_class = "other"
            span = 1
    elif variant_type == "deletion":
        deletion_length = len(target.ref) - 1
        observed_anchor = _query_slice(read, query_position, 1)
        observed_ref = _query_slice(read, query_position, len(target.ref))
        if pileup_read.indel == -deletion_length and observed_anchor == target.alt:
            allele_class = "alt"
            span = 1
        elif pileup_read.indel == 0 and observed_ref == target.ref:
            allele_class = "ref"
            span = len(target.ref)
        else:
            allele_class = "other"
            span = 1
    else:
        return "other", None, None
    quality = _base_quality(read, query_position, span)
    return allele_class, quality, float(query_position + 1)


def _orientation(read: pysam.AlignedSegment) -> str:
    if (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse):
        return "f1r2"
    if (read.is_read2 and not read.is_reverse) or (read.is_read1 and read.is_reverse):
        return "f2r1"
    return "unclassified"


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
            is_target = (
                variant_type == "insertion"
                and anchor == target_zero
                and length == len(target.alt) - 1
            )
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
    allele_classes = {observation.allele_class for observation in observations}
    if len(allele_classes) != 1:
        representative = max(
            observations, key=lambda item: (item.mapping_quality, item.base_quality or -1)
        )
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
    return (
        max(observations, key=lambda item: (item.base_quality or -1, item.mapping_quality)),
        False,
    )


def _mean(values: list[float]) -> float | None:
    return float(statistics.fmean(values)) if values else None


def _empty_metrics(policy: ReadbackPolicy) -> dict[str, Any]:
    return {
        "eligible_fragments": 0,
        "ref_reads": 0,
        "alt_reads": 0,
        "other_reads": 0,
        "forward_ref_reads": 0,
        "reverse_ref_reads": 0,
        "forward_alt_reads": 0,
        "reverse_alt_reads": 0,
        "f1r2_ref_reads": 0,
        "f2r1_ref_reads": 0,
        "unclassified_orientation_ref_reads": 0,
        "f1r2_alt_reads": 0,
        "f2r1_alt_reads": 0,
        "unclassified_orientation_alt_reads": 0,
        "mean_ref_mapping_quality": None,
        "mean_alt_mapping_quality": None,
        "mean_ref_base_quality": None,
        "mean_alt_base_quality": None,
        "mean_ref_read_position": None,
        "mean_alt_read_position": None,
        "soft_clipped_ref_reads": 0,
        "soft_clipped_alt_reads": 0,
        "nearby_non_target_indel_reads": 0,
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
    }


def _observe_target(
    alignment: pysam.AlignmentFile,
    alignment_contig: str,
    target: Target,
    variant_type: str,
    policy: ReadbackPolicy,
) -> dict[str, Any]:
    metrics = _empty_metrics(policy)
    observations_by_fragment: dict[str, list[Observation]] = defaultdict(list)
    columns = alignment.pileup(
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
    )
    for column in columns:
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
            observations_by_fragment[read.query_name].append(
                Observation(
                    allele_class=allele_class,
                    mapping_quality=int(read.mapping_quality),
                    base_quality=base_quality,
                    read_position=read_position,
                    reverse=bool(read.is_reverse),
                    orientation=_orientation(read),
                    soft_clipped=any(
                        operation == CIGAR_SOFT_CLIP for operation, _ in read.cigartuples or ()
                    ),
                    nearby_non_target_indel=_has_nearby_non_target_indel(
                        read,
                        target,
                        variant_type,
                        policy.nearby_indel_window,
                    ),
                )
            )
    collapsed: list[Observation] = []
    for fragment_observations in observations_by_fragment.values():
        observation, discordant = _collapse_fragment(fragment_observations)
        metrics["collapsed_overlapping_mates"] += len(fragment_observations) - 1
        metrics["discordant_overlapping_mates"] += int(discordant)
        collapsed.append(observation)
    metrics["eligible_fragments"] = len(collapsed)
    metric_values: dict[str, list[float]] = {
        "ref_mapping": [],
        "alt_mapping": [],
        "ref_base": [],
        "alt_base": [],
        "ref_position": [],
        "alt_position": [],
    }
    for observation in collapsed:
        allele_class = observation.allele_class
        metrics[f"{allele_class}_reads"] += 1
        if allele_class not in {"ref", "alt"}:
            continue
        strand = "reverse" if observation.reverse else "forward"
        metrics[f"{strand}_{allele_class}_reads"] += 1
        metrics[f"{observation.orientation}_{allele_class}_reads"] += 1
        metric_values[f"{allele_class}_mapping"].append(observation.mapping_quality)
        if observation.base_quality is not None:
            metric_values[f"{allele_class}_base"].append(observation.base_quality)
        if observation.read_position is not None:
            metric_values[f"{allele_class}_position"].append(observation.read_position)
        metrics[f"soft_clipped_{allele_class}_reads"] += int(observation.soft_clipped)
        metrics["nearby_non_target_indel_reads"] += int(observation.nearby_non_target_indel)
    for allele_class in ("ref", "alt"):
        metrics[f"mean_{allele_class}_mapping_quality"] = _mean(
            metric_values[f"{allele_class}_mapping"]
        )
        metrics[f"mean_{allele_class}_base_quality"] = _mean(metric_values[f"{allele_class}_base"])
        metrics[f"mean_{allele_class}_read_position"] = _mean(
            metric_values[f"{allele_class}_position"]
        )
    return metrics


def parse_readback_targets(
    alignment_path: str | Path,
    reference_fasta: str | Path,
    target_manifest: str | Path,
    provenance: ReadbackProvenance,
    *,
    policy: ReadbackPolicy = DEFAULT_READBACK_POLICY,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if policy.minimum_mapping_quality < 0 or policy.minimum_base_quality < 0:
        raise PortableEvidenceError("read-back quality thresholds must be non-negative")
    if policy.minimum_alt_fragments < 1 or policy.nearby_indel_window < 0:
        raise PortableEvidenceError("read-back evidence thresholds are invalid")
    if sha256_file(Path(alignment_path)) != provenance.source_alignment_sha256:
        raise PortableEvidenceError("read-back alignment SHA-256 differs from provenance")
    if sha256_file(Path(reference_fasta)) != provenance.reference_sha256:
        raise PortableEvidenceError("read-back reference SHA-256 differs from provenance")
    targets, target_manifest_sha256 = _read_targets(target_manifest)
    if not targets and not Path(str(reference_fasta) + ".fai").is_file():
        raise PortableEvidenceError("empty selection requires an indexed read-back reference")
    try:
        reference = pysam.FastaFile(str(reference_fasta))
        alignment = pysam.AlignmentFile(
            str(alignment_path), reference_filename=str(reference_fasta)
        )
    except (OSError, ValueError) as exc:
        raise PortableEvidenceError("cannot open indexed read-back alignment/reference") from exc
    rows_without_artifact: list[dict[str, Any]] = []
    sample_bindings = {
        str(record.get("SM"))
        for record in alignment.header.to_dict().get("RG", [])
        if record.get("SM")
    }
    try:
        if not alignment.has_index():
            raise PortableEvidenceError("read-back alignment lacks a usable index")
        if not targets and not reference.references:
            raise PortableEvidenceError("empty selection requires a usable read-back reference")
        for target in targets:
            alignment_contig = _alignment_contig(target.chrom, alignment.references)
            reference_contig = _alignment_contig(target.chrom, reference.references)
            observed_ref = reference.fetch(
                reference_contig, target.pos - 1, target.pos - 1 + len(target.ref)
            ).upper()
            if observed_ref != target.ref:
                raise PortableEvidenceError("read-back target REF differs from the exact reference")
            variant_type = _variant_type(target.ref, target.alt)
            if variant_type == "complex":
                metrics = _empty_metrics(policy)
                assessment_status = "unsupported_allele_representation"
                confidence = "UNRESOLVED"
                filter_reasons = ["COMPLEX_ALLELE_REPRESENTATION_NOT_IMPLEMENTED"]
            else:
                metrics = _observe_target(alignment, alignment_contig, target, variant_type, policy)
                if metrics["eligible_fragments"] == 0:
                    assessment_status = "no_eligible_reads"
                    confidence = "UNRESOLVED"
                    filter_reasons = ["NO_ELIGIBLE_READS"]
                else:
                    assessment_status = "assessed"
                    filter_reasons = []
                    if metrics["alt_reads"] == 0:
                        filter_reasons.append("NO_ALT_SUPPORT")
                    elif metrics["alt_reads"] < policy.minimum_alt_fragments:
                        filter_reasons.append("LOW_ALT_SUPPORT")
                    if metrics["alt_reads"] > 0 and (
                        metrics["forward_alt_reads"] == 0 or metrics["reverse_alt_reads"] == 0
                    ):
                        filter_reasons.append("ALT_ONE_STRAND_ONLY")
                    if metrics["nearby_non_target_indel_reads"]:
                        filter_reasons.append("NEARBY_NON_TARGET_INDEL")
                    if metrics["discordant_overlapping_mates"]:
                        filter_reasons.append("DISCORDANT_OVERLAPPING_MATES")
                    confidence = (
                        "SUPPORTED"
                        if metrics["alt_reads"] >= policy.minimum_alt_fragments
                        and not filter_reasons
                        else "REVIEW"
                    )
            semantic_source = {
                "target": target.__dict__,
                "variant_type": variant_type,
                "assessment_status": assessment_status,
                "metrics": metrics,
                "filter_reasons": sorted(filter_reasons),
            }
            row: dict[str, Any] = {
                "target_id": target.target_id,
                "target_manifest_sha256": target_manifest_sha256,
                "chrom": target.chrom,
                "pos": target.pos,
                "ref": target.ref,
                "alt": target.alt,
                "variant_type": variant_type,
                "assessment_status": assessment_status,
                **metrics,
                "local_complexity_status": target.local_complexity_status,
                "local_haplotype_status": target.local_haplotype_status,
                "candidate_category": target.candidate_category,
                "evidence_id": "pending",
                "event_id": f"readback:{target.chrom}:{target.pos}:{target.ref}:{target.alt}",
                "assembly": "GRCh38",
                "reference_sha256": provenance.reference_sha256,
                "technical_source_role": provenance.technical_source_role,
                "caller_name": provenance.caller_name,
                "caller_version": provenance.caller_version,
                "caller_run_identity_sha256": provenance.caller_run_identity_sha256,
                "resource_manifest_sha256": provenance.resource_manifest_sha256,
                "source_alignment_sha256": provenance.source_alignment_sha256,
                "raw_artifact_id": provenance.raw_artifact_id,
                "raw_artifact_sha256": "0" * 64,
                "source_record_digest": hashlib.sha256(
                    json.dumps(
                        semantic_source,
                        sort_keys=True,
                        separators=(",", ":"),
                    ).encode()
                ).hexdigest(),
                "filter_reasons": sorted(filter_reasons),
                "confidence_tier": confidence,
                "limitations": [
                    "bounded short-read pileup evidence is not a clinical confirmation",
                    "absence of ALT reads does not establish a reference genotype",
                    "local haplotype assembly and orthogonal confirmation may still be required",
                ],
                "confirmation_required": True,
                "concordance_group_id": None,
                "concordance_status": "not_assessed",
            }
            rows_without_artifact.append(row)
    finally:
        alignment.close()
        reference.close()
    raw_identity_payload = {
        "target_manifest_sha256": target_manifest_sha256,
        "policy": policy.__dict__,
        "record_digests": sorted(row["source_record_digest"] for row in rows_without_artifact),
    }
    raw_artifact_sha256 = hashlib.sha256(
        json.dumps(raw_identity_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    for row in rows_without_artifact:
        row["raw_artifact_sha256"] = raw_artifact_sha256
        row["evidence_id"] = make_evidence_id("read_backed_evidence", row)
    private_bindings = (*forbidden_bindings, *sorted(sample_bindings))
    validated = _validate_readback_rows(rows_without_artifact, provenance, private_bindings)
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "targets": len(targets),
        "typed_rows": len(validated),
        "record_reconciliation_delta": len(targets) - len(validated),
        "assessment_status_counts": {
            status: sum(row["assessment_status"] == status for row in validated)
            for status in (
                "assessed",
                "no_eligible_reads",
                "unsupported_allele_representation",
            )
        },
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "unresolved_rows": sum(row["confidence_tier"] == "UNRESOLVED" for row in validated),
        "target_manifest_sha256": target_manifest_sha256,
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": raw_artifact_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    rendered = json.dumps(summary, sort_keys=True)
    if any(binding in rendered for binding in sample_bindings):
        raise PortableEvidenceError("read-back summary persisted an alignment sample identifier")
    return validated, summary


def _validate_readback_rows(
    rows: list[dict[str, Any]], provenance: ReadbackProvenance, forbidden: tuple[str, ...]
) -> list[dict[str, Any]]:
    if rows:
        return validate_rows(
            "read_backed_evidence",
            rows,
            expected_reference_sha256=provenance.reference_sha256,
            forbidden_bindings=forbidden,
        )
    for digest in (
        provenance.reference_sha256,
        provenance.caller_run_identity_sha256,
        provenance.resource_manifest_sha256,
        provenance.source_alignment_sha256,
    ):
        if not isinstance(digest, str) or not SHA256.fullmatch(digest):
            raise PortableEvidenceError("empty read-back selection lacks valid provenance")
    raw_id = provenance.raw_artifact_id
    if (
        provenance.technical_source_role not in SOURCE_ROLES
        or not isinstance(raw_id, str)
        or not raw_id.strip()
        or "/" in raw_id
        or "\\" in raw_id
        or any(value and value in raw_id for value in forbidden)
    ):
        raise PortableEvidenceError(
            "empty read-back selection has an invalid source or artifact ID"
        )
    return []


def _materialize(
    rows: list[dict[str, Any]],
    parsing: dict[str, Any],
    output_parquet: str | Path,
    provenance: ReadbackProvenance,
    forbidden: tuple[str, ...],
) -> dict[str, Any]:
    output = None
    if rows:
        output = write_extension(
            "read_backed_evidence",
            rows,
            output_parquet,
            expected_reference_sha256=provenance.reference_sha256,
            forbidden_bindings=forbidden,
        )
    elif Path(output_parquet).exists() or Path(output_parquet).is_symlink():
        raise PortableEvidenceError("existing extension conflicts with empty read-back selection")
    return {"status": "PASS", "parsing": parsing, "output": output}


def materialize_readback_targets(
    alignment_path: str | Path,
    reference_fasta: str | Path,
    target_manifest: str | Path,
    output_parquet: str | Path,
    provenance: ReadbackProvenance,
    *,
    policy: ReadbackPolicy = DEFAULT_READBACK_POLICY,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_readback_targets(
        alignment_path,
        reference_fasta,
        target_manifest,
        provenance,
        policy=policy,
        forbidden_bindings=forbidden_bindings,
    )
    return _materialize(rows, parsing, output_parquet, provenance, forbidden_bindings)


def _capture_policy(value: Any) -> ReadbackPolicy:
    expected = set(ReadbackPolicy.__dataclass_fields__)
    if not isinstance(value, dict) or set(value) != expected:
        raise PortableEvidenceError("read-back capture policy schema differs")
    try:
        policy = ReadbackPolicy(**value)
    except TypeError as exc:
        raise PortableEvidenceError("read-back capture policy is invalid") from exc
    if (
        isinstance(policy.include_duplicates, bool) is False
        or any(
            isinstance(item, bool) or not isinstance(item, int)
            for item in (
                policy.minimum_mapping_quality,
                policy.minimum_base_quality,
                policy.minimum_alt_fragments,
                policy.nearby_indel_window,
                policy.max_depth,
            )
        )
        or policy.minimum_mapping_quality < 0
        or policy.minimum_base_quality < 0
        or policy.minimum_alt_fragments < 1
        or policy.nearby_indel_window < 0
        or policy.max_depth < 1
    ):
        raise PortableEvidenceError("read-back capture policy is invalid")
    return policy


def _captured_assessment(
    variant_type: str,
    metrics: dict[str, Any],
    policy: ReadbackPolicy,
) -> tuple[str, list[str], str]:
    if variant_type == "complex":
        return (
            "unsupported_allele_representation",
            ["COMPLEX_ALLELE_REPRESENTATION_NOT_IMPLEMENTED"],
            "UNRESOLVED",
        )
    if metrics["eligible_fragments"] == 0:
        return "no_eligible_reads", ["NO_ELIGIBLE_READS"], "UNRESOLVED"
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
    reasons.sort()
    confidence = (
        "SUPPORTED"
        if metrics["alt_reads"] >= policy.minimum_alt_fragments and not reasons
        else "REVIEW"
    )
    return "assessed", reasons, confidence


def parse_readback_capture(
    capture_json: str | Path,
    provenance: ReadbackProvenance,
    *,
    allow_empty: bool = True,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Convert one receipt-bound OCI capture into portable typed rows."""

    capture_path = Path(capture_json)
    if capture_path.is_symlink() or not capture_path.is_file():
        raise PortableEvidenceError("read-back capture is not a regular file")
    try:
        payload = json.loads(capture_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError("read-back capture is unreadable") from exc
    expected_payload_fields = {
        "schema_version",
        "status",
        "assembly",
        "target_manifest_sha256",
        "policy",
        "records",
        "record_count",
        "sample_identifier_persisted",
        "read_names_persisted",
        "read_sequences_persisted",
        "private_paths_persisted",
    }
    if (
        not isinstance(payload, dict)
        or set(payload) != expected_payload_fields
        or payload.get("schema_version") != 1
        or payload.get("status") != "PASS"
        or payload.get("assembly") != "GRCh38"
        or payload.get("sample_identifier_persisted") is not False
        or payload.get("read_names_persisted") is not False
        or payload.get("read_sequences_persisted") is not False
        or payload.get("private_paths_persisted") is not False
    ):
        raise PortableEvidenceError("read-back capture boundary is invalid")
    target_manifest_sha256 = payload.get("target_manifest_sha256")
    if (
        not isinstance(target_manifest_sha256, str)
        or SHA256.fullmatch(target_manifest_sha256) is None
    ):
        raise PortableEvidenceError("read-back capture target digest is invalid")
    policy = _capture_policy(payload.get("policy"))
    records = payload.get("records")
    record_count = payload.get("record_count")
    if (
        not isinstance(records, list)
        or isinstance(record_count, bool)
        or not isinstance(record_count, int)
        or record_count != len(records)
    ):
        raise PortableEvidenceError("read-back capture record cardinality differs")
    if (not records) != (target_manifest_sha256 == EMPTY_TARGET_SHA256) or (
        not records and not allow_empty
    ):
        raise PortableEvidenceError(
            "read-back capture target digest or profile conflicts with its record count"
        )
    raw_artifact_sha256 = sha256_file(capture_path)
    expected_metric_fields = set(_empty_metrics(policy))
    expected_record_fields = {
        "target",
        "variant_type",
        "assessment_status",
        "metrics",
        "filter_reasons",
    }
    expected_target_fields = set(Target.__dataclass_fields__)
    seen_target_ids: set[str] = set()
    seen_alleles: set[tuple[str, int, str, str]] = set()
    rows: list[dict[str, Any]] = []
    for captured in records:
        if not isinstance(captured, dict) or set(captured) != expected_record_fields:
            raise PortableEvidenceError("read-back capture record schema differs")
        raw_target = captured.get("target")
        if not isinstance(raw_target, dict) or set(raw_target) != expected_target_fields:
            raise PortableEvidenceError("read-back capture target schema differs")
        try:
            target = Target(**raw_target)
        except TypeError as exc:
            raise PortableEvidenceError("read-back capture target is invalid") from exc
        if (
            not target.target_id
            or "/" in target.target_id
            or "\\" in target.target_id
            or target.target_id in seen_target_ids
            or isinstance(target.pos, bool)
            or not isinstance(target.pos, int)
            or target.pos < 1
        ):
            raise PortableEvidenceError("read-back capture target identity is invalid")
        seen_target_ids.add(target.target_id)
        chrom = _canonical_chrom(target.chrom)
        ref = target.ref.upper()
        alt = target.alt.upper()
        if (
            ref != target.ref
            or alt != target.alt
            or CONCRETE_ALLELE.fullmatch(ref) is None
            or CONCRETE_ALLELE.fullmatch(alt) is None
            or ref == alt
            or not target.candidate_category
            or not target.local_complexity_status
            or not target.local_haplotype_status
        ):
            raise PortableEvidenceError("read-back capture target allele is invalid")
        allele_key = (chrom, target.pos, ref, alt)
        if allele_key in seen_alleles:
            raise PortableEvidenceError("read-back capture target allele is duplicated")
        seen_alleles.add(allele_key)
        variant_type = _variant_type(ref, alt)
        if captured.get("variant_type") != variant_type:
            raise PortableEvidenceError("read-back capture variant classification differs")
        metrics = captured.get("metrics")
        if not isinstance(metrics, dict) or set(metrics) != expected_metric_fields:
            raise PortableEvidenceError("read-back capture metric schema differs")
        assessment, filter_reasons, confidence = _captured_assessment(variant_type, metrics, policy)
        if (
            captured.get("assessment_status") != assessment
            or captured.get("filter_reasons") != filter_reasons
        ):
            raise PortableEvidenceError("read-back capture assessment differs from its metrics")
        semantic_source = {
            "target": target.__dict__,
            "variant_type": variant_type,
            "assessment_status": assessment,
            "metrics": metrics,
            "filter_reasons": filter_reasons,
        }
        row: dict[str, Any] = {
            "target_id": target.target_id,
            "target_manifest_sha256": target_manifest_sha256,
            "chrom": chrom,
            "pos": target.pos,
            "ref": ref,
            "alt": alt,
            "variant_type": variant_type,
            "assessment_status": assessment,
            **metrics,
            "local_complexity_status": target.local_complexity_status,
            "local_haplotype_status": target.local_haplotype_status,
            "candidate_category": target.candidate_category,
            "evidence_id": "pending",
            "event_id": f"readback:{chrom}:{target.pos}:{ref}:{alt}",
            "assembly": "GRCh38",
            "reference_sha256": provenance.reference_sha256,
            "technical_source_role": provenance.technical_source_role,
            "caller_name": provenance.caller_name,
            "caller_version": provenance.caller_version,
            "caller_run_identity_sha256": provenance.caller_run_identity_sha256,
            "resource_manifest_sha256": provenance.resource_manifest_sha256,
            "source_alignment_sha256": provenance.source_alignment_sha256,
            "raw_artifact_id": provenance.raw_artifact_id,
            "raw_artifact_sha256": raw_artifact_sha256,
            "source_record_digest": hashlib.sha256(
                json.dumps(
                    semantic_source,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
            "filter_reasons": filter_reasons,
            "confidence_tier": confidence,
            "limitations": [
                "bounded short-read pileup evidence is not a clinical confirmation",
                "absence of ALT reads does not establish a reference genotype",
                "local haplotype assembly and orthogonal confirmation may still be required",
            ],
            "confirmation_required": True,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
        row["evidence_id"] = make_evidence_id("read_backed_evidence", row)
        rows.append(row)
    validated = _validate_readback_rows(rows, provenance, forbidden_bindings)
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "capture_records": record_count,
        "typed_rows": len(validated),
        "record_reconciliation_delta": record_count - len(validated),
        "assessment_status_counts": {
            status: sum(row["assessment_status"] == status for row in validated)
            for status in (
                "assessed",
                "no_eligible_reads",
                "unsupported_allele_representation",
            )
        },
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "unresolved_rows": sum(row["confidence_tier"] == "UNRESOLVED" for row in validated),
        "target_manifest_sha256": target_manifest_sha256,
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": raw_artifact_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "read_names_persisted": False,
        "read_sequences_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def materialize_readback_capture(
    capture_json: str | Path,
    output_parquet: str | Path,
    provenance: ReadbackProvenance,
    *,
    allow_empty: bool = True,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_readback_capture(
        capture_json,
        provenance,
        allow_empty=allow_empty,
        forbidden_bindings=forbidden_bindings,
    )
    return _materialize(rows, parsing, output_parquet, provenance, forbidden_bindings)


def main() -> int:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--alignment", type=Path)
    source.add_argument("--capture-json", type=Path)
    parser.add_argument("--reference-fasta", type=Path)
    parser.add_argument("--target-manifest", type=Path)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--technical-source-role", required=True)
    parser.add_argument("--caller-run-identity-sha256", required=True)
    parser.add_argument("--resource-manifest-sha256", required=True)
    parser.add_argument("--source-alignment-sha256", required=True)
    parser.add_argument("--raw-artifact-id", required=True)
    parser.add_argument("--caller-version", default="1")
    parser.add_argument("--minimum-mapping-quality", type=int, default=20)
    parser.add_argument("--minimum-base-quality", type=int, default=20)
    parser.add_argument("--minimum-alt-fragments", type=int, default=3)
    parser.add_argument("--nearby-indel-window", type=int, default=10)
    parser.add_argument("--include-duplicates", action="store_true")
    args = parser.parse_args()
    provenance = ReadbackProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        caller_version=args.caller_version,
    )
    if args.capture_json is not None:
        if args.reference_fasta is not None or args.target_manifest is not None:
            parser.error("--capture-json cannot be combined with direct read inputs")
        result = materialize_readback_capture(
            args.capture_json,
            args.output_parquet,
            provenance,
        )
    else:
        if args.reference_fasta is None or args.target_manifest is None:
            parser.error("direct read-back requires --reference-fasta and --target-manifest")
        policy = ReadbackPolicy(
            minimum_mapping_quality=args.minimum_mapping_quality,
            minimum_base_quality=args.minimum_base_quality,
            minimum_alt_fragments=args.minimum_alt_fragments,
            nearby_indel_window=args.nearby_indel_window,
            include_duplicates=args.include_duplicates,
        )
        result = materialize_readback_targets(
            args.alignment,
            args.reference_fasta,
            args.target_manifest,
            args.output_parquet,
            provenance,
            policy=policy,
        )
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
