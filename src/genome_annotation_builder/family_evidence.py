"""Privacy-minimized pedigree evidence without genotype refinement.

The adapter intentionally keeps only the bundle subject's accepted genotype and
small derived family statuses.  Relative identifiers and relative genotypes are
used transiently and are never written to the portable bundle.
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq
import pysam

from .arrow_schemas import family_evidence_schema
from .config import BuildConfig, FamilyEvidenceConfig
from .constants import STANDARD_CHROMS
from .errors import InputValidationError, ReconciliationError
from .parquet_io import PartitionedParquetWriter
from .pedigree import PedigreeMember, parse_pedigree
from .util import (
    atomic_output_directory,
    atomic_write_json,
    directory_digest,
    sha256_file,
    sha256_json,
)
from .vcf_ingest import canonical_core_chrom

CONCRETE_ALLELE = re.compile(r"[ACGT]+")
AUTOSOMES = frozenset(f"chr{value}" for value in range(1, 23))
TRANSMISSION_STATE = {
    (1, 0): "inherited_paternal",
    (0, 1): "inherited_maternal",
    (1, 1): "inherited_biparental",
    (0, 0): "not_assessed",  # No carried ALT to assign.
}
CHROM_RANK = {chrom: rank for rank, chrom in enumerate(STANDARD_CHROMS)}
SourceAlleleKey = tuple[str, int, str, str]
LOCO_REVIEW_COLUMNS = frozenset(
    {
        "chrom",
        "pos",
        "ref",
        "alt",
        "status",
        "child_gt",
        "exact_allele_in_loco_model",
    }
)
LOCO_CHILD_UNINFORMATIVE = "uninformative_child_genotype"
LOCO_EXCEPTION_STATUSES = frozenset(
    {
        LOCO_CHILD_UNINFORMATIVE,
        "uninformative_parent_not_callable",
        "de_novo_candidate",
        "mendelian_conflict",
    }
)


@dataclass(frozen=True)
class CoreVariant:
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: tuple[int, ...]
    phased: bool


@dataclass(frozen=True)
class NormalizationProjection:
    source_key: SourceAlleleKey
    source_identity: tuple[str, int]
    source_gt: tuple[int, ...]
    source_phased: bool
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    projected_gt: tuple[int, ...]
    projected_phased: bool


@dataclass(frozen=True)
class FamilyVariant:
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    subject_gt: tuple[int, ...]
    subject_phased: bool
    father_gt: tuple[int | None, ...] | None
    mother_gt: tuple[int | None, ...] | None
    phase_set: str | None
    family_confidence: str | None
    family_reasons: tuple[str, ...]
    caller_support: int | None
    replicate_support: int | None
    de_novo_score: float | None
    difficult_region_flags: tuple[str, ...]


@dataclass(frozen=True)
class LocoReviewRow:
    key: tuple[int, int, str, str]
    status: str
    child_gt: tuple[int, ...] | None
    exact_allele_in_model: bool

    def concordance(self, accepted: tuple[int, ...]) -> bool | None:
        if self.child_gt is None:
            return None
        return sorted(self.child_gt) == sorted(accepted)


def _index_path(vcf: Path) -> Path | None:
    for candidate in (Path(f"{vcf}.tbi"), Path(f"{vcf}.csi")):
        if candidate.is_file():
            return candidate
    return None


def _contains_pedigree_identifier(text: str, identifiers: tuple[str, ...]) -> bool:
    folded = text.casefold()
    for identifier in identifiers:
        if len(identifier) >= 3 and identifier in folded:
            return True
        if len(identifier) < 3 and re.search(
            rf"(?<![a-z0-9]){re.escape(identifier)}(?![a-z0-9])", folded
        ):
            return True
    return False


def _verify_loco_identity(document: dict[str, Any], sample: str, review_sha256: str) -> None:
    if document.get("heldout_child") != sample:
        raise InputValidationError("LOCO summary held-out subject does not match input sample")
    if "status" in document and document["status"] != "PASS":
        raise InputValidationError("LOCO measurement did not pass")
    unit = document.get("validated_unit")
    if (unit is not None and unit != sample) or "excluded_person_units" in document:
        # A technical unit can have its own VCF sample name. The producer
        # explicitly binds it to the held-out person and excludes it from the model.
        excluded = document.get("excluded_person_units")
        model = document.get("model_samples")
        if (
            not isinstance(unit, str)
            or not unit
            or not isinstance(excluded, list)
            or unit not in excluded
            or not isinstance(model, list)
            or not model
            or sample in model
            or any(name in model for name in excluded)
            or document.get("model_header_excludes_all_person_units") is not True
        ):
            raise InputValidationError(
                "LOCO summary validated unit or model exclusions disagree with the held-out subject"
            )
    crosscheck = document.get("parental_haplotype_crosscheck")
    if not isinstance(crosscheck, dict) or crosscheck.get("model_is_child_independent") is not True:
        raise InputValidationError("LOCO summary does not prove a child-independent model")
    if "artifacts_sha256" in document:
        artifacts = document["artifacts_sha256"]
        if not isinstance(artifacts, dict) or artifacts.get("review.tsv") != review_sha256:
            raise InputValidationError("LOCO summary does not identify the configured review TSV")


def _load_loco_summary(config: BuildConfig) -> dict[str, Any] | None:
    loco = config.family_evidence.loco
    if loco is None:
        return None
    try:
        document = json.loads(loco.summary_json.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputValidationError(f"invalid LOCO summary JSON: {exc}") from exc
    if not isinstance(document, dict) or not isinstance(document.get("counters"), dict):
        raise InputValidationError("LOCO summary must contain a counters mapping")
    _verify_loco_identity(document, config.sample, loco.review_expected_sha256)
    counters = document["counters"]
    # The native producer serializes dict(Counter): zero bins are normally absent.
    # It can also observe source records without any carried concrete autosomal ALT.
    if "heldout_carried_alleles" not in counters and not (
        type(counters.get("source_records")) is int and counters["source_records"] > 0
    ):
        raise InputValidationError("LOCO summary lacks a carried-allele or source-record count")
    normalized: dict[str, int] = {}
    for key in {"heldout_carried_alleles", "mendelian_consistent", *LOCO_EXCEPTION_STATUSES}:
        value = counters.get(key, 0)
        if type(value) is not int or value < 0:
            raise InputValidationError(f"LOCO counter {key} must be a non-negative integer")
        normalized[key] = value
    review_rows = document.get("review_rows")
    informative = document.get("informative_alleles")
    if type(review_rows) is not int or review_rows < 0:
        raise InputValidationError("LOCO review_rows must be a non-negative integer")
    if type(informative) is not int or informative < 0:
        raise InputValidationError("LOCO informative_alleles must be a non-negative integer")
    exception_total = sum(normalized[status] for status in LOCO_EXCEPTION_STATUSES)
    if review_rows != exception_total:
        raise InputValidationError("LOCO review row accounting does not match exception counters")
    if normalized["heldout_carried_alleles"] != (normalized["mendelian_consistent"] + review_rows):
        raise InputValidationError("LOCO carried-allele accounting is inconsistent")
    expected_informative = (
        normalized["mendelian_consistent"]
        + normalized["de_novo_candidate"]
        + normalized["mendelian_conflict"]
    )
    if informative != expected_informative:
        raise InputValidationError("LOCO informative-allele accounting is inconsistent")
    fraction = document.get("mendelian_consistency_fraction")
    expected_fraction = normalized["mendelian_consistent"] / informative if informative else None
    if expected_fraction is None:
        if fraction is not None:
            raise InputValidationError(
                "LOCO consistency fraction must be null for zero denominator"
            )
    elif type(fraction) not in (int, float) or not math.isclose(
        float(fraction), expected_fraction, rel_tol=1e-12, abs_tol=1e-12
    ):
        raise InputValidationError("LOCO consistency fraction is inconsistent with counters")
    return {
        "review_rows": review_rows,
        "informative_alleles": informative,
        "mendelian_consistency_fraction": expected_fraction,
        "counters": normalized,
    }


def verify_family_resources(config: BuildConfig) -> dict[str, Any]:
    family = config.family_evidence
    if not family.enabled:
        return {"status": "not_enabled"}
    assert family.joint_vcf is not None
    assert family.joint_vcf_expected_sha256 is not None
    assert family.pedigree is not None
    assert family.pedigree_expected_sha256 is not None
    if not family.joint_vcf.is_file() or not family.pedigree.is_file():
        raise InputValidationError("family joint VCF or PED is missing")
    index = _index_path(family.joint_vcf)
    if index is None:
        raise InputValidationError("family joint VCF has no adjacent TBI/CSI index")
    joint_sha = sha256_file(family.joint_vcf)
    pedigree_sha = sha256_file(family.pedigree)
    if joint_sha != family.joint_vcf_expected_sha256:
        raise InputValidationError("family joint VCF SHA-256 mismatch")
    if pedigree_sha != family.pedigree_expected_sha256:
        raise InputValidationError("family PED SHA-256 mismatch")
    members = parse_pedigree(family.pedigree)
    if config.sample not in members:
        raise InputValidationError("configured sample is absent from the family PED")
    subject = members[config.sample]
    if family.loco is not None and len(subject.parents) != 2:
        raise InputValidationError("LOCO evidence requires two PED parents for the subject")

    try:
        joint = pysam.VariantFile(str(family.joint_vcf))
    except (OSError, ValueError) as exc:
        raise InputValidationError(f"cannot open family joint VCF: {exc}") from exc
    try:
        samples = set(joint.header.samples)
        required_samples = {config.sample, *subject.parents}
        if not required_samples.issubset(samples):
            raise InputValidationError("family joint VCF lacks the subject or a required parent")
        tags = family.tags
        for tag in (
            tags.confidence_format,
            tags.reasons_format,
            tags.caller_concordance_format,
            tags.technical_replicate_concordance_format,
            tags.de_novo_score_format,
            tags.phase_set_format,
        ):
            if tag is not None and tag not in joint.header.formats:
                raise InputValidationError(f"configured family FORMAT tag is absent: {tag}")
        for tag in tags.difficult_info_flags:
            if tag not in joint.header.info:
                raise InputValidationError(f"configured family INFO flag is absent: {tag}")
    finally:
        joint.close()

    pedigree_identifiers = tuple(
        sorted(
            {
                value.casefold()
                for member in members.values()
                for value in (member.individual, member.family_id)
            }
        )
    )
    if _contains_pedigree_identifier(family.source_version, pedigree_identifiers):
        raise InputValidationError("family source_version must not contain a pedigree identifier")

    loco_summary = None
    loco_resource = None
    if family.loco is not None:
        loco = family.loco
        if not loco.review_tsv.is_file() or not loco.summary_json.is_file():
            raise InputValidationError("configured LOCO review or summary file is missing")
        review_sha = sha256_file(loco.review_tsv)
        summary_sha = sha256_file(loco.summary_json)
        if review_sha != loco.review_expected_sha256:
            raise InputValidationError("LOCO review SHA-256 mismatch")
        if summary_sha != loco.summary_expected_sha256:
            raise InputValidationError("LOCO summary SHA-256 mismatch")
        loco_summary = _load_loco_summary(config)
        loco_resource = {
            "review_sha256": review_sha,
            "review_bytes": loco.review_tsv.stat().st_size,
            "summary_sha256": summary_sha,
            "summary_bytes": loco.summary_json.stat().st_size,
            "summary": loco_summary,
        }
    return {
        "status": "verified",
        "source_version": family.source_version,
        "snapshot_date": family.snapshot_date.isoformat(),
        "joint_vcf": {
            "sha256": joint_sha,
            "bytes": family.joint_vcf.stat().st_size,
            "index_sha256": sha256_file(index),
        },
        "pedigree": {
            "sha256": pedigree_sha,
            "bytes": family.pedigree.stat().st_size,
            "member_count": len(members),
            "subject_role": "child" if subject.parents else "founder",
            "parent_count": len(subject.parents),
        },
        "loco": loco_resource,
    }


def _parse_review_gt(value: str, status: str) -> tuple[int, ...] | None:
    if status == LOCO_CHILD_UNINFORMATIVE:
        # child_gt is the projection. The producer writes "." if it cannot
        # score the native child_source_gt against this specific ALT.
        if value != ".":
            raise InputValidationError(
                "uninformative LOCO child genotype must have a missing projection"
            )
        return None
    if not isinstance(value, str) or not re.fullmatch(r"[01](?:[|/][01])?", value):
        raise InputValidationError("LOCO review contains an invalid projected subject genotype")
    genotype = tuple(int(field) for field in re.split(r"[|/]", value))
    if 1 not in genotype:
        raise InputValidationError("LOCO review genotype does not carry the reviewed ALT")
    return genotype


def _parse_provenance_gt(value: str) -> tuple[tuple[int, ...], bool]:
    separators = {character for character in value if character in {"/", "|"}}
    if len(separators) > 1:
        raise ReconciliationError("variant provenance contains mixed genotype separators")
    fields = re.split(r"[|/]", value)
    if (
        len(fields) not in {1, 2}
        or any(not field.isdigit() for field in fields)
        or any(int(field) not in {0, 1} for field in fields)
    ):
        raise ReconciliationError("variant provenance contains an invalid source genotype")
    genotype = tuple(int(field) for field in fields)
    if 1 not in genotype:
        raise ReconciliationError("normalization projection does not carry its source ALT")
    return genotype, separators == {"|"}


def _projection_variant_identity(variant_id: str) -> tuple[str, int, str, str]:
    fields = variant_id.split(":")
    if len(fields) != 4:
        raise ReconciliationError("variant provenance contains an invalid normalized identity")
    chrom, pos_text, ref, alt = fields
    try:
        pos = int(pos_text)
    except ValueError as exc:
        raise ReconciliationError(
            "variant provenance contains a non-integer normalized position"
        ) from exc
    ref = ref.upper()
    alt = alt.upper()
    if (
        chrom not in CHROM_RANK
        or pos < 1
        or not CONCRETE_ALLELE.fullmatch(ref)
        or not CONCRETE_ALLELE.fullmatch(alt)
    ):
        raise ReconciliationError("variant provenance contains an invalid normalized allele")
    return chrom, pos, ref, alt


def _load_normalization_projections(
    variant_dir: Path,
) -> dict[SourceAlleleKey, tuple[NormalizationProjection, ...]]:
    root = variant_dir / "variant_provenance.parquet"
    files = sorted(root.rglob("*.parquet")) if root.is_dir() else []
    if not files:
        raise ReconciliationError("variant provenance dataset is missing")
    grouped: dict[SourceAlleleKey, list[NormalizationProjection]] = {}
    variant_ids: set[str] = set()
    for path in files:
        parquet = pq.ParquetFile(path)
        for batch in parquet.iter_batches(
            batch_size=100_000,
            columns=[
                "variant_id",
                "source_id",
                "source_alt_index",
                "source_chrom",
                "source_pos",
                "source_ref",
                "source_alt",
                "source_gt",
                "projected_gt",
                "phased",
                "normalization_changed",
            ],
        ):
            for raw in batch.to_pylist():
                if not raw["normalization_changed"]:
                    continue
                source_chrom = canonical_core_chrom(str(raw["source_chrom"]))
                source_pos = int(raw["source_pos"])
                source_ref = str(raw["source_ref"]).upper()
                source_alt = str(raw["source_alt"]).upper()
                if (
                    source_chrom not in CHROM_RANK
                    or source_pos < 1
                    or not CONCRETE_ALLELE.fullmatch(source_ref)
                    or not CONCRETE_ALLELE.fullmatch(source_alt)
                ):
                    raise ReconciliationError(
                        "variant provenance contains an invalid source allele"
                    )
                variant_id = str(raw["variant_id"])
                chrom, pos, ref, alt = _projection_variant_identity(variant_id)
                if chrom != source_chrom:
                    raise ReconciliationError(
                        "normalization projection moved an allele across chromosomes"
                    )
                source_gt, source_phased = _parse_provenance_gt(str(raw["source_gt"]))
                projected_raw = raw["projected_gt"]
                if (
                    not isinstance(projected_raw, list)
                    or len(projected_raw) not in {1, 2}
                    or any(allele not in {0, 1} for allele in projected_raw)
                    or 1 not in projected_raw
                ):
                    raise ReconciliationError(
                        "variant provenance contains an invalid projected genotype"
                    )
                projected_gt = tuple(int(allele) for allele in projected_raw)
                projected_phased = bool(raw["phased"])
                if len(projected_gt) != len(source_gt):
                    raise ReconciliationError(
                        "normalization projection changed source genotype ploidy"
                    )
                if len(source_gt) == 2 and projected_phased != source_phased:
                    raise ReconciliationError(
                        "normalization projection changed source phasing semantics"
                    )
                if len(source_gt) == 1:
                    # A separatorless haploid GT cannot encode phase. htslib
                    # nevertheless exposes a representational phase bit for
                    # it, so retain the normalized/core representation without
                    # pretending it can be compared to '/' or '|' in source GT.
                    source_phased = projected_phased
                source_identity = (str(raw["source_id"]), int(raw["source_alt_index"]))
                source_key = (source_chrom, source_pos, source_ref, source_alt)
                if variant_id in variant_ids:
                    raise ReconciliationError(
                        "variant provenance contains a duplicate normalized identity"
                    )
                variant_ids.add(variant_id)
                grouped.setdefault(source_key, []).append(
                    NormalizationProjection(
                        source_key=source_key,
                        source_identity=source_identity,
                        source_gt=source_gt,
                        source_phased=source_phased,
                        variant_id=variant_id,
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        projected_gt=projected_gt,
                        projected_phased=projected_phased,
                    )
                )
    result: dict[SourceAlleleKey, tuple[NormalizationProjection, ...]] = {}
    for source_key, projections in grouped.items():
        identities = {projection.source_identity for projection in projections}
        source_genotypes = {
            (projection.source_gt, projection.source_phased) for projection in projections
        }
        if len(identities) != 1 or len(source_genotypes) != 1:
            raise ReconciliationError(
                "normalization projections disagree on their source allele identity"
            )
        result[source_key] = tuple(
            sorted(
                projections,
                key=lambda projection: (
                    CHROM_RANK[projection.chrom],
                    projection.pos,
                    projection.ref,
                    projection.alt,
                ),
            )
        )
    return result


def _review_key(chrom: str, pos: str, ref: str, alt: str) -> tuple[int, int, str, str]:
    canonical = canonical_core_chrom(chrom)
    if canonical not in AUTOSOMES:
        raise InputValidationError("LOCO review must contain only canonical autosomes")
    try:
        position = int(pos)
    except ValueError as exc:
        raise InputValidationError("LOCO review contains a non-integer position") from exc
    ref = ref.upper()
    alt = alt.upper()
    if position < 1 or not CONCRETE_ALLELE.fullmatch(ref) or not CONCRETE_ALLELE.fullmatch(alt):
        raise InputValidationError("LOCO review contains an invalid allele")
    return CHROM_RANK[canonical], position, ref, alt


def _boolean_text(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise InputValidationError("LOCO exact-allele flag must be True or False")


class LocoReviewStream:
    """Sorted join over the exception-only LOCO ledger and normalized alleles."""

    def __init__(
        self,
        path: Path,
        expected_summary: dict[str, Any],
        normalization_projections: dict[SourceAlleleKey, tuple[NormalizationProjection, ...]],
    ) -> None:
        self._expected_summary = expected_summary
        source_previous: tuple[int, int, str, str] | None = None
        source_keys: set[tuple[int, int, str, str]] = set()
        projected_rows: list[LocoReviewRow] = []
        self._current: LocoReviewRow | None = None
        self.status_counts: Counter[str] = Counter()
        self.rows = 0
        self.normalization_projected_source_rows = 0
        self.normalization_projection_source_gt_discordant = 0
        self.source_order_inversions_normalized = 0
        self.matched = 0
        self.unmatched = 0
        self.matched_genotype_unassessable = 0
        self.matched_genotype_concordant = 0
        self.matched_genotype_discordant = 0
        self.discordant_status_counts: Counter[str] = Counter()
        with path.open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or not LOCO_REVIEW_COLUMNS.issubset(reader.fieldnames):
                raise InputValidationError("LOCO review header is missing required columns")
            for raw in reader:
                key = _review_key(raw["chrom"], raw["pos"], raw["ref"], raw["alt"])
                if key in source_keys:
                    raise InputValidationError("LOCO review contains a duplicate source allele")
                source_keys.add(key)
                if source_previous is not None and key < source_previous:
                    self.source_order_inversions_normalized += 1
                source_previous = key
                status = raw["status"]
                if status not in LOCO_EXCEPTION_STATUSES:
                    raise InputValidationError(f"unsupported LOCO review status: {status}")
                child_gt = _parse_review_gt(raw["child_gt"], status)
                exact = _boolean_text(raw["exact_allele_in_loco_model"])
                source_key = (STANDARD_CHROMS[key[0]], key[1], key[2], key[3])
                projections = normalization_projections.get(source_key)
                if projections is None:
                    projected_rows.append(
                        LocoReviewRow(
                            key=key,
                            status=status,
                            child_gt=child_gt,
                            exact_allele_in_model=exact,
                        )
                    )
                else:
                    expected_gt = projections[0].source_gt
                    if child_gt is not None and sorted(child_gt) != sorted(expected_gt):
                        self.normalization_projection_source_gt_discordant += 1
                    self.normalization_projected_source_rows += 1
                    projected_rows.extend(
                        LocoReviewRow(
                            key=(
                                CHROM_RANK[projection.chrom],
                                projection.pos,
                                projection.ref,
                                projection.alt,
                            ),
                            status=status,
                            child_gt=child_gt,
                            exact_allele_in_model=exact,
                        )
                        for projection in projections
                    )
                self.rows += 1
                self.status_counts[status] += 1
        projected_rows.sort(key=lambda row: row.key)
        if any(current.key <= previous.key for previous, current in pairwise(projected_rows)):
            raise ReconciliationError("normalized LOCO review contains a duplicate allele")
        self._rows = iter(projected_rows)
        self.projected_rows = len(projected_rows)
        self._advance()

    def _advance(self) -> None:
        try:
            self._current = next(self._rows)
        except StopIteration:
            self._current = None

    def match(
        self, key: tuple[int, int, str, str], subject_gt: tuple[int, ...]
    ) -> LocoReviewRow | None:
        while self._current is not None and self._current.key < key:
            self.unmatched += 1
            self._advance()
        if self._current is None or self._current.key > key:
            return None
        row = self._current
        concordance = row.concordance(subject_gt)
        if concordance is None:
            self.matched_genotype_unassessable += 1
        elif concordance:
            self.matched_genotype_concordant += 1
        else:
            self.matched_genotype_discordant += 1
            self.discordant_status_counts[row.status] += 1
        self.matched += 1
        self._advance()
        return row

    def close(self) -> dict[str, Any]:
        while self._current is not None:
            self.unmatched += 1
            self._advance()
        expected = self._expected_summary
        if self.rows != expected["review_rows"]:
            raise ReconciliationError("LOCO review row count disagrees with its summary")
        for status in LOCO_EXCEPTION_STATUSES:
            if self.status_counts[status] != expected["counters"][status]:
                raise ReconciliationError(
                    f"LOCO review status count disagrees with summary: {status}"
                )
        return {
            "review_rows": self.rows,
            "projected_review_rows": self.projected_rows,
            "normalization_projected_source_rows": self.normalization_projected_source_rows,
            "normalization_projection_source_gt_discordant": (
                self.normalization_projection_source_gt_discordant
            ),
            "source_order_inversions_normalized": self.source_order_inversions_normalized,
            "matched_subject_variants": self.matched,
            "matched_genotype_concordant": self.matched_genotype_concordant,
            "matched_genotype_unassessable": self.matched_genotype_unassessable,
            "matched_genotype_discordant": self.matched_genotype_discordant,
            "genotype_discordant_source_status_counts": dict(
                sorted(self.discordant_status_counts.items())
            ),
            "projected_review_rows_outside_subject_core": self.unmatched,
            "status_counts": dict(sorted(self.status_counts.items())),
        }


def _iter_core_groups(
    variant_dir: Path, chrom: str
) -> Iterator[tuple[int, dict[str, CoreVariant]]]:
    root = variant_dir / "variants.parquet" / f"chrom={chrom}"
    files = sorted(root.glob("*.parquet")) if root.is_dir() else []
    previous_pos = 0
    group_pos: int | None = None
    group: dict[str, CoreVariant] = {}
    for path in files:
        parquet = pq.ParquetFile(path)
        for batch in parquet.iter_batches(
            batch_size=100_000,
            columns=[
                "variant_id",
                "chrom",
                "pos",
                "ref",
                "alt",
                "genotype.gt",
                "genotype.phased",
            ],
        ):
            for raw in batch.to_pylist():
                if raw["chrom"] != chrom:
                    raise ReconciliationError("core variant partition/chromosome mismatch")
                pos = int(raw["pos"])
                if pos < previous_pos:
                    raise ReconciliationError("core variants are not position sorted")
                previous_pos = pos
                if group_pos is not None and pos != group_pos:
                    yield group_pos, group
                    group = {}
                group_pos = pos
                ref = str(raw["ref"]).upper()
                alt = str(raw["alt"]).upper()
                variant_id = str(raw["variant_id"])
                expected_id = f"{chrom}:{pos}:{ref}:{alt}"
                if variant_id != expected_id or variant_id in group:
                    raise ReconciliationError("core variant identity is invalid or duplicated")
                genotype = raw["genotype"]
                group[variant_id] = CoreVariant(
                    variant_id=variant_id,
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    gt=tuple(int(value) for value in genotype["gt"]),
                    phased=bool(genotype["phased"]),
                )
    if group_pos is not None:
        yield group_pos, group


def _sample_gt(sample: pysam.libcbcf.VariantRecordSample) -> tuple[int | None, ...] | None:
    raw = sample.get("GT")
    if raw is None:
        return None
    return tuple(None if allele is None else int(allele) for allele in raw)


def _format_value(sample: pysam.libcbcf.VariantRecordSample, tag: str | None) -> Any:
    return sample.get(tag) if tag is not None else None


def _reasons(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    raw_values = value if isinstance(value, tuple) else (value,)
    return tuple(str(item) for item in raw_values if item not in {None, "", "."})


def _optional_int(value: Any, label: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, tuple):
        if len(value) != 1:
            raise InputValidationError(f"family {label} FORMAT value is not scalar")
        value = value[0]
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise InputValidationError(f"family {label} FORMAT value is not an integer") from exc


def _optional_float(value: Any, label: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, tuple):
        if len(value) != 1:
            raise InputValidationError(f"family {label} FORMAT value is not scalar")
        value = value[0]
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise InputValidationError(f"family {label} FORMAT value is not numeric") from exc
    if not math.isfinite(result):
        raise InputValidationError(f"family {label} FORMAT value is not finite")
    return result


def _source_contig(header: pysam.VariantHeader, chrom: str) -> str | None:
    token = chrom[3:]
    candidates = [chrom, token]
    if chrom == "chrM":
        candidates.extend(["MT", "chrMT"])
    return next((value for value in candidates if value in header.contigs), None)


def _carried_source_key(
    record: pysam.VariantRecord, chrom: str, counters: Counter[str]
) -> SourceAlleleKey | None:
    if not record.alts or len(record.alts) != 1:
        raise InputValidationError("carried family records must be source-biallelic")
    ref = record.ref.upper()
    alt = record.alts[0].upper()
    if not CONCRETE_ALLELE.fullmatch(ref) or not CONCRETE_ALLELE.fullmatch(alt):
        counters["carried_symbolic_or_non_acgt_records_omitted"] += 1
        return None
    observed_chrom = canonical_core_chrom(record.contig)
    if observed_chrom != chrom:
        raise ReconciliationError("family joint contig does not match canonical partition")
    return chrom, record.pos, ref, alt


def _family_variant_from_record(
    record: pysam.VariantRecord,
    subject: PedigreeMember,
    config: FamilyEvidenceConfig,
    forbidden_identifiers: tuple[str, ...],
    *,
    projection: NormalizationProjection | None,
) -> FamilyVariant:
    tags = config.tags
    sample = record.samples[subject.individual]
    raw_subject_gt = _sample_gt(sample)
    if (
        raw_subject_gt is None
        or len(raw_subject_gt) not in {1, 2}
        or any(allele not in {0, 1} for allele in raw_subject_gt)
        or not any(allele == 1 for allele in raw_subject_gt)
    ):
        raise ReconciliationError(
            "family source allele does not have a complete carried biallelic genotype"
        )
    source_gt = tuple(int(allele) for allele in raw_subject_gt)
    source_phased = bool(sample.phased)
    if projection is None:
        chrom = canonical_core_chrom(record.contig)
        if chrom not in CHROM_RANK or not record.alts or len(record.alts) != 1:
            raise ReconciliationError("family source allele identity is invalid")
        pos = record.pos
        ref = record.ref.upper()
        alt = record.alts[0].upper()
        subject_gt = source_gt
        subject_phased = source_phased
    else:
        same_source_gt = (
            source_gt == projection.source_gt
            if source_phased or projection.source_phased
            else sorted(source_gt) == sorted(projection.source_gt)
        )
        source_phase_matches = len(source_gt) == 1 or source_phased == projection.source_phased
        if not same_source_gt or not source_phase_matches:
            raise ReconciliationError(
                "family source genotype/phasing differs from variant provenance"
            )
        chrom = projection.chrom
        pos = projection.pos
        ref = projection.ref
        alt = projection.alt
        subject_gt = projection.projected_gt
        subject_phased = projection.projected_phased
    family_confidence_value = _format_value(sample, tags.confidence_format)
    confidence = (
        str(family_confidence_value) if family_confidence_value not in {None, "", "."} else None
    )
    reasons = _reasons(_format_value(sample, tags.reasons_format))
    for value in (*reasons, confidence or ""):
        if _contains_pedigree_identifier(value.casefold(), forbidden_identifiers):
            raise InputValidationError(
                "family evidence tag contains a pedigree identifier; refusing privacy leak"
            )
    phase_set_value = _format_value(sample, tags.phase_set_format)
    phase_set = str(phase_set_value) if phase_set_value not in {None, "", "."} else None
    caller_support = _optional_int(
        _format_value(sample, tags.caller_concordance_format), "caller concordance"
    )
    replicate_support = _optional_int(
        _format_value(sample, tags.technical_replicate_concordance_format),
        "technical replicate concordance",
    )
    de_novo_score = _optional_float(
        _format_value(sample, tags.de_novo_score_format), "de-novo score"
    )
    difficult = tuple(flag for flag in tags.difficult_info_flags if bool(record.info.get(flag)))
    father_gt = _sample_gt(record.samples[subject.father]) if subject.father is not None else None
    mother_gt = _sample_gt(record.samples[subject.mother]) if subject.mother is not None else None
    return FamilyVariant(
        variant_id=f"{chrom}:{pos}:{ref}:{alt}",
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        subject_gt=subject_gt,
        subject_phased=subject_phased,
        father_gt=father_gt,
        mother_gt=mother_gt,
        phase_set=phase_set,
        family_confidence=confidence,
        family_reasons=reasons,
        caller_support=caller_support,
        replicate_support=replicate_support,
        de_novo_score=de_novo_score,
        difficult_region_flags=difficult,
    )


def _projected_family_exception_groups(
    joint: pysam.VariantFile,
    chrom: str,
    subject: PedigreeMember,
    config: FamilyEvidenceConfig,
    forbidden_identifiers: tuple[str, ...],
    normalization_projections: dict[SourceAlleleKey, tuple[NormalizationProjection, ...]],
) -> list[tuple[int, dict[str, FamilyVariant]]]:
    by_position: dict[int, dict[str, FamilyVariant]] = {}
    for source_key, projections in normalization_projections.items():
        if source_key[0] != chrom:
            continue
        contig = _source_contig(joint.header, source_key[0])
        if contig is None:
            raise ReconciliationError("family joint lacks a normalization source contig")
        matches = [
            record
            for record in joint.fetch(contig, source_key[1] - 1, source_key[1])
            if record.pos == source_key[1]
            and record.ref.upper() == source_key[2]
            and record.alts
            and len(record.alts) == 1
            and record.alts[0].upper() == source_key[3]
        ]
        if len(matches) != 1:
            raise ReconciliationError(
                "family joint does not contain exactly one normalization source allele"
            )
        record = matches[0]
        for projection in projections:
            variant = _family_variant_from_record(
                record,
                subject,
                config,
                forbidden_identifiers,
                projection=projection,
            )
            group = by_position.setdefault(variant.pos, {})
            if variant.variant_id in group:
                raise ReconciliationError(
                    "normalization projection produced a duplicate family allele"
                )
            group[variant.variant_id] = variant
    return sorted(by_position.items())


def _iter_family_groups(
    joint: pysam.VariantFile,
    chrom: str,
    subject: PedigreeMember,
    config: FamilyEvidenceConfig,
    forbidden_identifiers: tuple[str, ...],
    counters: Counter[str],
    normalization_projections: dict[SourceAlleleKey, tuple[NormalizationProjection, ...]],
) -> Iterator[tuple[int, dict[str, FamilyVariant]]]:
    contig = _source_contig(joint.header, chrom)
    if contig is None:
        return
    exception_groups = _projected_family_exception_groups(
        joint,
        chrom,
        subject,
        config,
        forbidden_identifiers,
        normalization_projections,
    )

    def ordinary_groups() -> Iterator[tuple[int, dict[str, FamilyVariant]]]:
        previous_pos = 0
        group_pos: int | None = None
        group: dict[str, FamilyVariant] = {}
        for record in joint.fetch(contig):
            counters["joint_records_scanned"] += 1
            if record.pos < previous_pos:
                raise InputValidationError("family joint VCF is not position sorted")
            previous_pos = record.pos
            sample = record.samples[subject.individual]
            raw_subject_gt = _sample_gt(sample)
            if raw_subject_gt is None:
                counters["subject_absent_genotype_records_omitted"] += 1
                continue
            if any(allele is None for allele in raw_subject_gt):
                counters["subject_incomplete_genotype_records_omitted"] += 1
                if any(allele is not None and allele > 0 for allele in raw_subject_gt):
                    counters["subject_incomplete_genotype_records_with_alt_evidence"] += 1
                continue
            if not any(allele > 0 for allele in raw_subject_gt):
                continue
            source_key = _carried_source_key(record, chrom, counters)
            if source_key is None:
                continue
            if source_key in normalization_projections:
                counters["normalization_projected_family_source_records"] += 1
                continue
            variant = _family_variant_from_record(
                record,
                subject,
                config,
                forbidden_identifiers,
                projection=None,
            )
            if group_pos is not None and variant.pos != group_pos:
                yield group_pos, group
                group = {}
            group_pos = variant.pos
            if variant.variant_id in group:
                raise InputValidationError(
                    "family joint VCF contains a duplicate source-biallelic allele"
                )
            group[variant.variant_id] = variant
        if group_pos is not None:
            yield group_pos, group

    ordinary = iter(ordinary_groups())
    exceptions = iter(exception_groups)
    ordinary_group = _advance_group(ordinary)
    exception_group = _advance_group(exceptions)
    while ordinary_group is not None or exception_group is not None:
        if ordinary_group is None:
            assert exception_group is not None
            yield exception_group
            exception_group = _advance_group(exceptions)
        elif exception_group is None or ordinary_group[0] < exception_group[0]:
            yield ordinary_group
            ordinary_group = _advance_group(ordinary)
        elif exception_group[0] < ordinary_group[0]:
            yield exception_group
            exception_group = _advance_group(exceptions)
        else:
            if set(ordinary_group[1]) & set(exception_group[1]):
                raise ReconciliationError(
                    "ordinary and projected family alleles contain a duplicate identity"
                )
            yield ordinary_group[0], ordinary_group[1] | exception_group[1]
            ordinary_group = _advance_group(ordinary)
            exception_group = _advance_group(exceptions)


def _advance_group(
    iterator: Iterator[tuple[int, dict[str, Any]]],
) -> tuple[int, dict[str, Any]] | None:
    try:
        return next(iterator)
    except StopIteration:
        return None


def _complete_diploid(gt: tuple[int | None, ...] | None) -> tuple[int, int] | None:
    if gt is None or len(gt) != 2 or any(allele not in {0, 1} for allele in gt):
        return None
    return int(gt[0]), int(gt[1])  # type: ignore[arg-type]


def _pedigree_status(variant: FamilyVariant, subject: PedigreeMember) -> tuple[str, str, str]:
    if not subject.parents:
        return "not_applicable_founder", "not_applicable", "founder_observation"
    if variant.chrom not in AUTOSOMES or len(subject.parents) != 2:
        return "not_assessed", "unknown", "not_assessed"
    if len(variant.subject_gt) != 2:
        return "not_assessed", "unknown", "not_assessed"
    father = _complete_diploid(variant.father_gt)
    mother = _complete_diploid(variant.mother_gt)
    child = _complete_diploid(variant.subject_gt)
    if father is None or mother is None or child is None:
        return "parent_or_subject_no_call", "unknown", "unresolved_parent_no_call"
    father_alt = 1 in father
    mother_alt = 1 in mother
    if father_alt and mother_alt:
        origin = "both_parents_observed"
    elif father_alt:
        origin = "paternal_observed"
    elif mother_alt:
        origin = "maternal_observed"
    else:
        origin = "neither_parent_observed"
    # Each valid pair is ordered by PED role. The ordering of the child's GT
    # does not label a paternal or maternal chromosome, even when GT is phased.
    child_alt_count = sum(child)
    transmissions = {
        (paternal, maternal)
        for paternal in father
        for maternal in mother
        if paternal + maternal == child_alt_count
    }
    if not transmissions:
        return "conflict", origin, "mendelian_conflict"
    if len(transmissions) > 1:
        return "consistent", origin, "inherited_parent_ambiguous"
    return "consistent", origin, TRANSMISSION_STATE[transmissions.pop()]


def _caller_status(value: int | None) -> str:
    if value is None:
        return "not_assessed"
    if value == 1:
        return "exact_genotype_concordant"
    if value == 0:
        return "discordant"
    raise InputValidationError("caller concordance FORMAT value must be 0, 1, or missing")


def _replicate_status(value: int | None) -> str:
    if value is None:
        return "not_assessed"
    if value == 2:
        return "both_sources_exact_genotype_concordant"
    if value == 1:
        return "single_source_support"
    if value == 0:
        return "discordant"
    raise InputValidationError(
        "technical replicate concordance FORMAT value must be 0, 1, 2, or missing"
    )


def _phase_confidence(variant: FamilyVariant) -> str:
    if variant.subject_phased and variant.phase_set is not None:
        return "phased_with_phase_set"
    if variant.subject_phased:
        return "phased_without_phase_set"
    return "unphased"


def _read_evidence_status(confidence: str | None) -> str:
    if confidence == "HIGH_CONFIDENCE":
        return "family_policy_high"
    if confidence == "SUPPORTED":
        return "family_policy_supported"
    if confidence is None:
        return "not_available"
    return "review_required"


def _token_present(reasons: tuple[str, ...], configured: tuple[str, ...]) -> bool:
    return bool(set(reasons) & set(configured))


def _make_row(
    variant: FamilyVariant,
    subject: PedigreeMember,
    config: FamilyEvidenceConfig,
    joint_sha256: str,
    pedigree_sha256: str,
    loco_row: LocoReviewRow | None,
) -> dict[str, Any]:
    mendelian_status, parental_origin, inheritance_state = _pedigree_status(variant, subject)
    loco_source_status = loco_row.status if loco_row is not None else None
    loco_subject_gt_concordant = loco_row.concordance(variant.subject_gt) if loco_row else None
    if loco_row is not None:
        loco_status = (
            "source_genotype_discordant" if loco_subject_gt_concordant is False else loco_row.status
        )
    elif config.loco is None:
        loco_status = "not_configured"
    elif variant.chrom in AUTOSOMES:
        loco_status = "not_individually_persisted"
    else:
        loco_status = "not_assessed_outside_autosomal_loco"
    tags = config.tags
    de_novo_candidate = (
        variant.family_confidence in tags.de_novo_confidence_values
        or loco_status == "de_novo_candidate"
    )
    source_conflict = variant.family_confidence in tags.conflict_confidence_values
    if source_conflict and mendelian_status == "consistent":
        raise ReconciliationError(
            "family confidence marks a conflict but direct autosomal genotypes are consistent"
        )
    if source_conflict or de_novo_candidate:
        mendelian_status = "conflict"
        if subject.parents:
            inheritance_state = "de_novo_candidate" if de_novo_candidate else "mendelian_conflict"
    subject_payload = {
        "variant_id": variant.variant_id,
        "subject_gt": list(variant.subject_gt),
        "subject_phased": variant.subject_phased,
        "phase_set": variant.phase_set,
        "family_confidence": variant.family_confidence,
        "family_reasons": list(variant.family_reasons),
        "inheritance_state": inheritance_state,
        "parental_origin": parental_origin,
        "mendelian_status": mendelian_status,
        "loco_status": loco_status,
        "loco_source_status": loco_source_status,
        "loco_subject_gt_concordant": loco_subject_gt_concordant,
        "loco_exact_allele_in_model": (
            loco_row.exact_allele_in_model if loco_row is not None else None
        ),
        "caller_support": variant.caller_support,
        "replicate_support": variant.replicate_support,
        "difficult_region_flags": list(variant.difficult_region_flags),
    }
    subject_digest = sha256_json(subject_payload)
    evidence_id = "family-evidence:" + sha256_json(
        {
            "variant_id": variant.variant_id,
            "joint_vcf_sha256": joint_sha256,
            "subject_evidence_digest": subject_digest,
            "source_version": config.source_version,
        }
    )
    return {
        "evidence_id": evidence_id,
        "variant_id": variant.variant_id,
        "chrom": variant.chrom,
        "pos": variant.pos,
        "ref": variant.ref,
        "alt": variant.alt,
        "source_version": config.source_version,
        "snapshot_date": config.snapshot_date,
        "joint_vcf_sha256": joint_sha256,
        "pedigree_sha256": pedigree_sha256,
        "loco_review_sha256": (
            config.loco.review_expected_sha256 if config.loco is not None else None
        ),
        "loco_summary_sha256": (
            config.loco.summary_expected_sha256 if config.loco is not None else None
        ),
        "subject_gt": list(variant.subject_gt),
        "subject_phased": variant.subject_phased,
        "phase_set": variant.phase_set,
        "phase_confidence": _phase_confidence(variant),
        "family_confidence": variant.family_confidence,
        "family_reasons": list(variant.family_reasons),
        "inheritance_state": inheritance_state,
        "parental_origin": parental_origin,
        "mendelian_status": mendelian_status,
        "loco_status": loco_status,
        "loco_source_status": loco_source_status,
        "loco_subject_gt_concordant": loco_subject_gt_concordant,
        "loco_exact_allele_in_model": (
            loco_row.exact_allele_in_model if loco_row is not None else None
        ),
        "de_novo_score": variant.de_novo_score,
        "de_novo_candidate": de_novo_candidate,
        "parental_mosaicism_candidate": _token_present(
            variant.family_reasons, tags.parental_mosaicism_reason_tokens
        ),
        "parental_dropout_suspected": _token_present(
            variant.family_reasons, tags.parental_dropout_reason_tokens
        ),
        "child_dropout_suspected": _token_present(
            variant.family_reasons, tags.child_dropout_reason_tokens
        ),
        "caller_concordance": _caller_status(variant.caller_support),
        "technical_replicate_concordance": _replicate_status(variant.replicate_support),
        "complex_region": bool(variant.difficult_region_flags),
        "difficult_region_flags": list(variant.difficult_region_flags),
        "read_evidence_status": _read_evidence_status(variant.family_confidence),
        "subject_evidence_digest": subject_digest,
    }


def build_family_evidence(
    config: BuildConfig,
    *,
    run_key: str,
    variant_dir: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    """Reconcile family evidence to every accepted core variant."""

    family = config.family_evidence
    if not family.enabled:
        raise InputValidationError("family evidence materialization was called while disabled")
    assert family.joint_vcf is not None
    assert family.joint_vcf_expected_sha256 is not None
    assert family.pedigree is not None
    assert family.pedigree_expected_sha256 is not None
    members = parse_pedigree(family.pedigree)
    subject = members[config.sample]
    forbidden_identifiers = tuple(
        sorted(
            {member.individual.casefold() for member in members.values()}
            | {member.family_id.casefold() for member in members.values()}
        )
    )
    variant_dir = Path(variant_dir)
    normalization_projections = _load_normalization_projections(variant_dir)
    normalization_projection_summary = {
        "source_alleles": len(normalization_projections),
        "output_alleles": sum(len(values) for values in normalization_projections.values()),
        "extra_output_alleles": sum(
            len(values) - 1 for values in normalization_projections.values()
        ),
        "variant_provenance_sha256": directory_digest(variant_dir / "variant_provenance.parquet"),
        "incomplete_subject_genotype_policy": (
            "omit from family variant evidence as no-call; never convert partial GT to "
            "reference or a complete alternate genotype"
        ),
        "changed_allele_policy": (
            "project family and LOCO annotations through exact variant provenance while "
            "preserving source and projected genotype/phasing"
        ),
    }
    loco_summary = _load_loco_summary(config)
    loco_stream = (
        LocoReviewStream(
            family.loco.review_tsv,
            loco_summary,
            normalization_projections,
        )
        if family.loco is not None and loco_summary is not None
        else None
    )
    counters: Counter[str] = Counter()
    counters["normalization_projection_source_alleles"] = len(normalization_projections)
    counters["normalization_projection_output_alleles"] = sum(
        len(values) for values in normalization_projections.values()
    )
    destination = Path(destination)
    with atomic_output_directory(destination) as staging:
        output = staging / "family_evidence.parquet"
        joint = pysam.VariantFile(str(family.joint_vcf))
        try:
            with PartitionedParquetWriter(output, family_evidence_schema()) as writer:
                for chrom in STANDARD_CHROMS:
                    core_iterator = iter(_iter_core_groups(variant_dir, chrom))
                    family_iterator = iter(
                        _iter_family_groups(
                            joint,
                            chrom,
                            subject,
                            family,
                            forbidden_identifiers,
                            counters,
                            normalization_projections,
                        )
                    )
                    core_group = _advance_group(core_iterator)
                    family_group = _advance_group(family_iterator)
                    while core_group is not None or family_group is not None:
                        if core_group is None:
                            raise ReconciliationError(
                                "family joint contains a carried allele absent from the core ledger"
                            )
                        if family_group is None:
                            raise ReconciliationError(
                                "core ledger contains an allele absent from the family joint"
                            )
                        core_pos, core_rows = core_group
                        family_pos, family_rows = family_group
                        if core_pos != family_pos:
                            relation = (
                                "core allele absent from family joint"
                                if core_pos < family_pos
                                else "family carried allele absent from core ledger"
                            )
                            raise ReconciliationError(f"{relation} on {chrom}")
                        if set(core_rows) != set(family_rows):
                            raise ReconciliationError(
                                f"core/family exact-allele mismatch at {chrom}:{core_pos}"
                            )
                        for variant_id in sorted(core_rows):
                            core = core_rows[variant_id]
                            source = family_rows[variant_id]
                            same_gt = (
                                core.gt == source.subject_gt
                                if core.phased or source.subject_phased
                                else sorted(core.gt) == sorted(source.subject_gt)
                            )
                            if not same_gt or core.phased != source.subject_phased:
                                raise ReconciliationError(
                                    "family evidence subject genotype/phasing differs from core"
                                )
                            loco_row = None
                            if loco_stream is not None and chrom in AUTOSOMES:
                                loco_row = loco_stream.match(
                                    (CHROM_RANK[chrom], core.pos, core.ref, core.alt), core.gt
                                )
                            row = _make_row(
                                source,
                                subject,
                                family,
                                family.joint_vcf_expected_sha256,
                                family.pedigree_expected_sha256,
                                loco_row,
                            )
                            writer.write(row)
                            counters["family_evidence_rows"] += 1
                            counters[f"inheritance_{row['inheritance_state']}"] += 1
                            counters[f"mendelian_{row['mendelian_status']}"] += 1
                            counters[f"loco_{row['loco_status']}"] += 1
                            counters[f"caller_{row['caller_concordance']}"] += 1
                            counters["replicate_" + row["technical_replicate_concordance"]] += 1
                            counters["complex_region_rows"] += int(row["complex_region"])
                        core_group = _advance_group(core_iterator)
                        family_group = _advance_group(family_iterator)
        finally:
            joint.close()
        loco_accounting = loco_stream.close() if loco_stream is not None else None
        normalization_projection_summary.update(
            {
                "family_source_records_projected": counters.get(
                    "normalization_projected_family_source_records", 0
                ),
                "incomplete_subject_genotype_records_omitted": counters.get(
                    "subject_incomplete_genotype_records_omitted", 0
                ),
                "incomplete_subject_genotype_records_with_alt_evidence": counters.get(
                    "subject_incomplete_genotype_records_with_alt_evidence", 0
                ),
            }
        )
        summary = {
            "status": "PASS",
            "build_run_key": run_key,
            "source_version": family.source_version,
            "snapshot_date": family.snapshot_date.isoformat(),
            "joint_vcf_sha256": family.joint_vcf_expected_sha256,
            "pedigree_sha256": family.pedigree_expected_sha256,
            "subject_role": "child" if subject.parents else "founder",
            "relative_identifiers_persisted": 0,
            "relative_genotypes_persisted": 0,
            "genotype_refinements": 0,
            "normalization_projection": normalization_projection_summary,
            "counters": dict(sorted(counters.items())),
            "loco": (
                {
                    "review_sha256": family.loco.review_expected_sha256,
                    "summary_sha256": family.loco.summary_expected_sha256,
                    "summary_metrics": loco_summary,
                    "materialization": loco_accounting,
                    "absence_semantics": (
                        "an absent exception row is not promoted to per-variant consistency; "
                        "it remains not_individually_persisted"
                    ),
                }
                if family.loco is not None
                else None
            ),
            "outputs": {"family_evidence.parquet": directory_digest(output)},
        }
        atomic_write_json(staging / "family_evidence_summary.json", summary)
    return summary
