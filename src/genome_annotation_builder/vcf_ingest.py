"""Lossless classification and allele projection for reference-confidence VCFs."""

from __future__ import annotations

import base64
import hashlib
import heapq
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

import pysam

from .arrow_schemas import (
    callable_region_schema,
    non_core_variant_schema,
    reference_evidence_schema,
)
from .config import BuildConfig
from .constants import SYMBOLIC_REFERENCE_ALLELES
from .errors import InputValidationError
from .parquet_io import PartitionedParquetWriter
from .util import (
    MultisetDigest,
    atomic_output_directory,
    atomic_write_json,
    canonical_json_bytes,
    directory_digest,
    sha256_file,
)

CONCRETE_ALLELE = re.compile(r"^[ACGT]+$")


class _CallableRegionMerger:
    def __init__(self, writer: PartitionedParquetWriter) -> None:
        self.writer = writer
        self._pending_chrom: str | None = None
        self._pending_start: int | None = None
        self._pending: list[tuple[int, bool, int | None]] = []
        self._chrom: str | None = None
        self._cursor: int | None = None
        self._active: dict[int, tuple[int, bool, int | None]] = {}
        self._end_events: list[tuple[int, int]] = []
        self._next_interval_id = 0
        self._output_active: dict[str, Any] | None = None

    def add(
        self,
        *,
        chrom: str,
        start_pos: int,
        end_pos: int,
        callable_value: bool,
        minimum_depth: int | None,
    ) -> None:
        core_chrom = canonical_core_chrom(chrom)
        if core_chrom is None:
            return
        if start_pos < 1 or end_pos < start_pos:
            raise InputValidationError(
                f"invalid callability interval {core_chrom}:{start_pos}-{end_pos}"
            )

        if self._pending_chrom is None:
            self._set_pending(core_chrom, start_pos, end_pos, callable_value, minimum_depth)
            return

        assert self._pending_start is not None
        if core_chrom == self._pending_chrom and start_pos == self._pending_start:
            self._pending.append((end_pos, callable_value, minimum_depth))
            return
        if core_chrom == self._pending_chrom and start_pos < self._pending_start:
            raise InputValidationError("callability intervals are not coordinate sorted")

        previous_chrom = self._pending_chrom
        self._commit_pending()
        if core_chrom != previous_chrom:
            self._finish_chrom()
        self._set_pending(core_chrom, start_pos, end_pos, callable_value, minimum_depth)

    def _set_pending(
        self,
        chrom: str,
        start_pos: int,
        end_pos: int,
        callable_value: bool,
        minimum_depth: int | None,
    ) -> None:
        self._pending_chrom = chrom
        self._pending_start = start_pos
        self._pending = [(end_pos, callable_value, minimum_depth)]

    def _expire(self) -> None:
        assert self._cursor is not None
        while self._end_events and self._end_events[0][0] <= self._cursor:
            _boundary, interval_id = heapq.heappop(self._end_events)
            self._active.pop(interval_id, None)

    def _emit_segment(self, start_pos: int, end_pos: int) -> None:
        if not self._active or self._chrom is None or end_pos < start_pos:
            return
        callable_value = all(item[1] for item in self._active.values())
        thresholds = [item[2] for item in self._active.values() if item[2] is not None]
        minimum_depth = max(thresholds) if thresholds else None
        candidate = {
            "chrom": self._chrom,
            "start_pos": start_pos,
            "end_pos": end_pos,
            "callable": callable_value,
            "minimum_depth": minimum_depth,
            "assay_scope": "small_variant_reference_confidence_vcf",
            "derivation": "conservative_sweep_of_record_FILTER_FT_DP_GQ_thresholds",
        }
        active = self._output_active
        if (
            active is not None
            and active["chrom"] == candidate["chrom"]
            and active["callable"] == candidate["callable"]
            and active["minimum_depth"] == candidate["minimum_depth"]
            and candidate["start_pos"] == active["end_pos"] + 1
        ):
            active["end_pos"] = candidate["end_pos"]
            return
        self._flush_output()
        self._output_active = candidate

    def _advance_to(self, target: int) -> None:
        if self._cursor is None:
            self._cursor = target
        while self._cursor < target:
            self._expire()
            if not self._active:
                self._cursor = target
                break
            next_boundary = min(target, self._end_events[0][0])
            self._emit_segment(self._cursor, next_boundary - 1)
            self._cursor = next_boundary
        self._expire()

    def _commit_pending(self) -> None:
        if self._pending_chrom is None or self._pending_start is None:
            return
        if self._chrom is None:
            self._chrom = self._pending_chrom
            self._cursor = self._pending_start
        elif self._chrom != self._pending_chrom:
            raise InputValidationError("internal callability chromosome transition error")

        self._advance_to(self._pending_start)
        for end_pos, callable_value, minimum_depth in self._pending:
            interval_id = self._next_interval_id
            self._next_interval_id += 1
            self._active[interval_id] = (end_pos, callable_value, minimum_depth)
            heapq.heappush(self._end_events, (end_pos + 1, interval_id))
        self._pending_chrom = None
        self._pending_start = None
        self._pending = []

    def _flush_output(self) -> None:
        if self._output_active is None:
            return
        identity = dict(self._output_active)
        row = dict(self._output_active)
        row["region_id"] = "region:" + hashlib.sha256(canonical_json_bytes(identity)).hexdigest()
        self.writer.write(row)
        self._output_active = None

    def _finish_chrom(self) -> None:
        if self._chrom is None:
            return
        if self._active:
            final_boundary = max(item[0] for item in self._active.values()) + 1
            self._advance_to(final_boundary)
        self._flush_output()
        self._chrom = None
        self._cursor = None
        self._active.clear()
        self._end_events.clear()

    def flush(self) -> None:
        self._commit_pending()
        self._finish_chrom()


def canonical_core_chrom(contig: str) -> str | None:
    token = contig[3:] if contig.lower().startswith("chr") else contig
    token = token.upper()
    if token == "MT":
        token = "M"
    if token in {str(value) for value in range(1, 23)} | {"X", "Y", "M"}:
        return f"chr{token}"
    return None


def _is_concrete(allele: str | None) -> bool:
    return bool(allele and CONCRETE_ALLELE.fullmatch(allele.upper()))


def _integer(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        if len(value) != 1:
            return None
        value = value[0]
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _integer_list(value: Any) -> list[int] | None:
    if value is None:
        return None
    raw = value if isinstance(value, (tuple, list)) else (value,)
    converted = [_integer(item) for item in raw]
    return [item for item in converted if item is not None] if converted else None


def _json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        # VCF permits non-finite floats. Keep their spelling in valid JSON.
        return str(value)
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _canonical_json_text(value: Any) -> str | None:
    if not value:
        return None
    return json.dumps(
        _json_value(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _encode_text(value: str | None) -> str | None:
    if value is None:
        return None
    return base64.urlsafe_b64encode(value.encode("utf-8")).decode("ascii").rstrip("=")


def _filters(value: Any) -> list[str]:
    if value is None:
        return []
    if hasattr(value, "keys"):
        raw = list(value.keys())
    elif isinstance(value, str):
        raw = re.split(r"[;,]", value)
    elif isinstance(value, (tuple, list, set)):
        raw = list(value)
    else:
        raw = [value]
    return [str(item) for item in raw if item not in {None, "", "."}]


def _filters_pass(record_filters: list[str], genotype_filters: list[str]) -> bool:
    accepted = {"PASS"}
    return all(item in accepted for item in record_filters + genotype_filters)


def _source_gt_text(gt: tuple[int | None, ...], phased: bool) -> str:
    separator = "|" if phased else "/"
    return separator.join("." if allele is None else str(allele) for allele in gt)


def source_id_from_digest(input_sha256: str, ordinal: int, record_digest: str) -> str:
    return hashlib.sha256(
        canonical_json_bytes(
            {
                "input_sha256": input_sha256,
                "source_ordinal": ordinal,
                "source_record_digest": record_digest,
            }
        )
    ).hexdigest()


def _source_identifiers(input_sha256: str, ordinal: int, record_text: str) -> tuple[str, str]:
    record_digest = hashlib.sha256(record_text.encode("utf-8")).hexdigest()
    source_id = source_id_from_digest(input_sha256, ordinal, record_digest)
    return source_id, record_digest


def _build_projection_header(input_header: pysam.VariantHeader) -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value="genome-annotation-builder")
    for name, contig in input_header.contigs.items():
        length = contig.length
        header.contigs.add(name, length=length if length is not None else None)
    header.add_meta(
        "FILTER",
        items=[
            ("ID", "GAB_SOURCE_FILTERED"),
            ("Description", "Source record or genotype was not PASS"),
        ],
    )
    info_fields = {
        "GAB_SOURCE_ID": ("1", "String", "Stable source record identifier"),
        "GAB_SOURCE_DIGEST": ("1", "String", "SHA-256 of the source VCF row"),
        "GAB_SOURCE_ORDINAL": ("1", "Integer", "One-based source record ordinal"),
        "GAB_ALT_INDEX": ("1", "Integer", "One-based ALT index in the source record"),
        "GAB_ORIG_CONTIG": ("1", "String", "Source contig before core canonicalization"),
        "GAB_ORIG_POS": ("1", "Integer", "Source one-based position before normalization"),
        "GAB_ORIG_REF": ("1", "String", "Source REF before normalization"),
        "GAB_ORIG_ALT": ("1", "String", "Selected source ALT before normalization"),
        "GAB_ORIG_GT": ("1", "String", "Complete source genotype"),
        "GAB_PROJECTED_GT": (
            "1",
            "String",
            "Selected-ALT biallelic genotype before normalization and atomization",
        ),
        "GAB_RSID": ("1", "String", "Source rs identifier when unambiguous"),
        "GAB_SOURCE_FILTERS": ("1", "String", "Pipe-separated source FILTER/FT values"),
        "GAB_REF_AD": ("1", "Integer", "Source reference allele depth"),
        "GAB_ALT_AD": ("1", "Integer", "Source selected ALT allele depth"),
        "GAB_TOTAL_AD": ("1", "Integer", "Sum of all available source allele depths"),
        "GAB_ORIG_AD": (".", "Integer", "Complete source allele depth vector"),
        "GAB_ORIG_PL": (".", "Integer", "Complete source phred likelihood vector"),
        "GAB_SOURCE_INFO_B64": ("1", "String", "Base64url canonical source INFO JSON"),
        "GAB_OTHER_FORMAT_B64": (
            "1",
            "String",
            "Base64url canonical non-core FORMAT JSON",
        ),
    }
    for identifier, (number, field_type, description) in info_fields.items():
        header.add_meta(
            "INFO",
            items=[
                ("ID", identifier),
                ("Number", number),
                ("Type", field_type),
                ("Description", description),
            ],
        )
    header.formats.add("GT", number=1, type="String", description="Projected biallelic genotype")
    header.formats.add("DP", number=1, type="Integer", description="Source read depth")
    header.formats.add("GQ", number=1, type="Integer", description="Source genotype quality")
    header.formats.add(
        "AD", number="R", type="Integer", description="Source REF and selected ALT depths"
    )
    header.add_sample("SAMPLE")
    return header


def _reference_row(
    *,
    evidence_id: str,
    chrom: str,
    start_pos: int,
    end_pos: int,
    state: str,
    reference_observed: bool | None,
    callable_value: bool | None,
    confidence: str | None,
    gt: tuple[int | None, ...],
    phased: bool,
    dp: int | None,
    minimum_dp: int | None,
    gq: int | None,
    ad: list[int] | None,
    pl: list[int] | None,
    record_filters: list[str],
    genotype_filters: list[str],
    source_info_json: str | None,
    other_format_json: str | None,
    source_record_digest: str,
    source_ordinal: int,
    source_alts: tuple[str, ...],
    reason: str | None,
) -> dict[str, Any]:
    return {
        "evidence_id": evidence_id,
        "chrom": canonical_core_chrom(chrom) or chrom,
        "start_pos": start_pos,
        "end_pos": end_pos,
        "state": state,
        "reference_observed": reference_observed,
        "callable": callable_value,
        "call_confidence": confidence,
        "ploidy": len(gt) if gt else None,
        "gt": [allele for allele in gt] if gt else None,
        "phased": phased if gt else None,
        "dp": dp,
        "minimum_dp": minimum_dp,
        "gq": gq,
        "ad": ad,
        "pl": pl,
        "record_filters": record_filters,
        "genotype_filters": genotype_filters,
        "source_info_json": source_info_json,
        "other_format_json": other_format_json,
        "source_record_digest": source_record_digest,
        "source_ordinal": source_ordinal,
        "source_alts": list(source_alts),
        "reason": reason,
    }


def ingest_reference_confidence_vcf(
    config: BuildConfig,
    input_sha256: str,
    destination: str | Path,
) -> dict[str, Any]:
    """Fully scan input and produce a projected VCF plus lossless evidence tables."""
    destination = Path(destination)
    with atomic_output_directory(destination) as staging:
        projected_path = staging / "projected.vcf.gz"
        reference_path = staging / "reference_evidence.parquet"
        non_core_path = staging / "non_core_variants.parquet"
        callable_regions_path = staging / "callable_regions.parquet"

        try:
            source = pysam.VariantFile(str(config.input_vcf))
        except (OSError, ValueError) as exc:
            raise InputValidationError(f"cannot open input VCF for full scan: {exc}") from exc
        if config.sample not in source.header.samples:
            source.close()
            raise InputValidationError("configured sample disappeared from input header")
        projection_header = _build_projection_header(source.header)
        projected = pysam.VariantFile(str(projected_path), "wz", header=projection_header)

        reference = pysam.FastaFile(str(config.reference_fasta))
        contig_order = {name: index for index, name in enumerate(source.header.contigs)}
        previous_order = -1
        previous_pos = -1
        current_contig: str | None = None
        current_reference = ""
        counters: Counter[str] = Counter()
        accepted_evidence_digest = MultisetDigest()

        try:
            with (
                PartitionedParquetWriter(
                    reference_path, reference_evidence_schema(), partition_field=None
                ) as evidence_writer,
                PartitionedParquetWriter(
                    non_core_path, non_core_variant_schema(), partition_field=None
                ) as non_core_writer,
                PartitionedParquetWriter(
                    callable_regions_path, callable_region_schema(), partition_field=None
                ) as callable_writer,
            ):
                callable_merger = _CallableRegionMerger(callable_writer)
                for ordinal, record in enumerate(source, start=1):
                    counters["records_total"] += 1
                    order = contig_order.get(record.contig)
                    if order is None:
                        raise InputValidationError(
                            f"record contig absent from header: {record.contig}"
                        )
                    if order < previous_order or (
                        order == previous_order and record.pos < previous_pos
                    ):
                        raise InputValidationError(
                            f"input is not coordinate sorted at record {ordinal}: "
                            f"{record.contig}:{record.pos}"
                        )
                    previous_order, previous_pos = order, record.pos

                    if record.contig != current_contig:
                        if record.contig not in reference.references:
                            raise InputValidationError(
                                f"record contig absent from FASTA: {record.contig}"
                            )
                        current_contig = record.contig
                        current_reference = reference.fetch(record.contig).upper()
                    ref = record.ref.upper()
                    end_pos = int(record.stop)
                    if record.pos < 1 or end_pos < record.pos or end_pos > len(current_reference):
                        raise InputValidationError(
                            f"invalid interval at source record {ordinal}: "
                            f"{record.contig}:{record.pos}-{end_pos}"
                        )
                    observed_ref = current_reference[record.pos - 1 : record.pos - 1 + len(ref)]
                    if observed_ref != ref:
                        raise InputValidationError(
                            f"REF mismatch at source record {ordinal} "
                            f"{record.contig}:{record.pos}: "
                            f"VCF={ref}, FASTA={observed_ref}"
                        )

                    record_text = str(record).rstrip("\n")
                    source_id, source_digest = _source_identifiers(
                        input_sha256, ordinal, record_text
                    )
                    sample = record.samples[config.sample]
                    raw_gt = sample.get("GT")
                    gt = tuple(raw_gt) if raw_gt is not None else ()
                    if len(gt) not in {0, 1, 2}:
                        raise InputValidationError(
                            f"unsupported ploidy {len(gt)} at source record {ordinal}"
                        )
                    phased = bool(sample.phased) if gt else False
                    complete = bool(gt) and all(
                        isinstance(allele, int) and allele >= 0 for allele in gt
                    )
                    alts = tuple(str(alt).upper() for alt in (record.alts or ()))
                    if any(allele is not None and not 0 <= allele <= len(alts) for allele in gt):
                        raise InputValidationError(
                            f"genotype allele index exceeds ALT count at source record {ordinal}"
                        )

                    dp = _integer(sample.get("DP"))
                    gq = _integer(sample.get("GQ"))
                    minimum_dp = _integer(sample.get("MIN_DP"))
                    ad = _integer_list(sample.get("AD"))
                    pl = _integer_list(sample.get("PL"))
                    source_info_json = _canonical_json_text(dict(record.info.items()))
                    other_format_json = _canonical_json_text(
                        {
                            key: sample.get(key)
                            for key in sample
                            if key not in {"GT", "DP", "MIN_DP", "GQ", "AD", "PL", "FT"}
                            and sample.get(key) is not None
                        }
                    )
                    record_filters = _filters(record.filter)
                    genotype_filters = _filters(sample.get("FT"))
                    passed = _filters_pass(record_filters, genotype_filters)
                    confidence = config.quality.confidence(dp, gq) if passed else "low"
                    effective_dp = minimum_dp if minimum_dp is not None else dp
                    callable_value = (
                        passed and complete and config.quality.callable(effective_dp, gq)
                    )
                    callable_merger.add(
                        chrom=record.contig,
                        start_pos=record.pos,
                        end_pos=end_pos,
                        callable_value=callable_value,
                        minimum_depth=config.quality.callable_min_dp,
                    )

                    called_indices = sorted(
                        {allele for allele in gt if allele is not None and allele > 0}
                    )
                    if not complete:
                        counters["records_no_call"] += 1
                        evidence_writer.write(
                            _reference_row(
                                evidence_id=source_id,
                                chrom=record.contig,
                                start_pos=record.pos,
                                end_pos=end_pos,
                                state="no_call",
                                reference_observed=None,
                                callable_value=False,
                                confidence=None,
                                gt=gt,
                                phased=phased,
                                dp=dp,
                                minimum_dp=minimum_dp,
                                gq=gq,
                                ad=ad,
                                pl=pl,
                                record_filters=record_filters,
                                genotype_filters=genotype_filters,
                                source_info_json=source_info_json,
                                other_format_json=other_format_json,
                                source_record_digest=source_digest,
                                source_ordinal=ordinal,
                                source_alts=alts,
                                reason="incomplete_or_absent_GT",
                            )
                        )
                        # Retain the locus-level no-call. Only an explicit, concrete GT
                        # allele has a partial variant annotation; native FORMAT evidence
                        # alone cannot select a genotype or a different ALT.
                        called_indices = [
                            i
                            for i in called_indices
                            if canonical_core_chrom(record.contig)
                            and _is_concrete(ref)
                            and _is_concrete(alts[i - 1])
                        ]
                        if not called_indices:
                            continue
                    if not called_indices:
                        symbolic_reference = bool(alts) and all(
                            alt in SYMBOLIC_REFERENCE_ALLELES for alt in alts
                        )
                        state = (
                            "reference_block" if symbolic_reference else "explicit_reference_site"
                        )
                        counters[f"records_{state}"] += 1
                        evidence_writer.write(
                            _reference_row(
                                evidence_id=source_id,
                                chrom=record.contig,
                                start_pos=record.pos,
                                end_pos=end_pos,
                                state=state,
                                reference_observed=True,
                                callable_value=callable_value,
                                confidence=confidence,
                                gt=gt,
                                phased=phased,
                                dp=dp,
                                minimum_dp=minimum_dp,
                                gq=gq,
                                ad=ad,
                                pl=pl,
                                record_filters=record_filters,
                                genotype_filters=genotype_filters,
                                source_info_json=source_info_json,
                                other_format_json=other_format_json,
                                source_record_digest=source_digest,
                                source_ordinal=ordinal,
                                source_alts=alts,
                                reason=None,
                            )
                        )
                        continue

                    ad_tuple = tuple(ad or ())
                    ad_total = sum(ad_tuple) if ad_tuple else None
                    source_gt = _source_gt_text(gt, phased)
                    source_filter_text = "|".join(record_filters + genotype_filters) or "PASS"
                    core_chrom = canonical_core_chrom(record.contig)
                    for alt_index in called_indices:
                        alt = alts[alt_index - 1]
                        evidence_id = f"{source_id}:a{alt_index}"
                        concrete = _is_concrete(ref) and _is_concrete(alt)
                        category: str | None = None
                        if not concrete:
                            category = "symbolic_or_non_acgt_allele"
                        elif core_chrom is None:
                            category = "non_core_contig"

                        if category is not None:
                            counters[f"called_alleles_{category}"] += 1
                            non_core_writer.write(
                                {
                                    "evidence_id": evidence_id,
                                    "chrom": core_chrom or record.contig,
                                    "pos": record.pos,
                                    "ref": ref,
                                    "alt": alt,
                                    "category": category,
                                    "source_alt_index": alt_index,
                                    "source_gt": [int(allele) for allele in gt],
                                    "phased": phased,
                                    "dp": dp,
                                    "gq": gq,
                                    "ad": ad,
                                    "pl": pl,
                                    "source_info_json": source_info_json,
                                    "other_format_json": other_format_json,
                                    "source_record_digest": source_digest,
                                    "source_ordinal": ordinal,
                                }
                            )
                            continue

                        projected_gt = tuple(
                            None if allele is None else int(allele == alt_index) for allele in gt
                        )
                        if 1 not in projected_gt:
                            raise InputValidationError(
                                "internal allele projection lost selected ALT"
                            )
                        output_record = projected.new_record(
                            contig=record.contig,
                            start=record.pos - 1,
                            stop=record.pos - 1 + len(ref),
                            id=f"gab{source_id[:20]}a{alt_index}",
                            alleles=(ref, alt),
                            qual=record.qual,
                        )
                        if passed:
                            output_record.filter.add("PASS")
                        else:
                            output_record.filter.add("GAB_SOURCE_FILTERED")
                        output_record.info["GAB_SOURCE_ID"] = source_id
                        output_record.info["GAB_SOURCE_DIGEST"] = source_digest
                        output_record.info["GAB_SOURCE_ORDINAL"] = ordinal
                        output_record.info["GAB_ALT_INDEX"] = alt_index
                        output_record.info["GAB_ORIG_CONTIG"] = record.contig
                        output_record.info["GAB_ORIG_POS"] = record.pos
                        output_record.info["GAB_ORIG_REF"] = ref
                        output_record.info["GAB_ORIG_ALT"] = alt
                        output_record.info["GAB_ORIG_GT"] = source_gt
                        output_record.info["GAB_PROJECTED_GT"] = _source_gt_text(
                            projected_gt, phased
                        )
                        output_record.info["GAB_SOURCE_FILTERS"] = source_filter_text
                        if record.id and re.fullmatch(r"rs\d+", record.id):
                            output_record.info["GAB_RSID"] = record.id
                        ref_ad = ad_tuple[0] if len(ad_tuple) > 0 else None
                        alt_ad = ad_tuple[alt_index] if len(ad_tuple) > alt_index else None
                        if ref_ad is not None:
                            output_record.info["GAB_REF_AD"] = ref_ad
                        if alt_ad is not None:
                            output_record.info["GAB_ALT_AD"] = alt_ad
                        if ad_total is not None:
                            output_record.info["GAB_TOTAL_AD"] = ad_total
                        if ad:
                            output_record.info["GAB_ORIG_AD"] = ad
                        if pl:
                            output_record.info["GAB_ORIG_PL"] = pl
                        encoded_info = _encode_text(source_info_json)
                        encoded_format = _encode_text(other_format_json)
                        if encoded_info:
                            output_record.info["GAB_SOURCE_INFO_B64"] = encoded_info
                        if encoded_format:
                            output_record.info["GAB_OTHER_FORMAT_B64"] = encoded_format
                        projected_sample = output_record.samples["SAMPLE"]
                        projected_sample["GT"] = projected_gt
                        projected_sample.phased = phased
                        if dp is not None:
                            projected_sample["DP"] = dp
                        if gq is not None:
                            projected_sample["GQ"] = gq
                        if ref_ad is not None and alt_ad is not None:
                            projected_sample["AD"] = (ref_ad, alt_ad)
                        projected.write(output_record)
                        accepted_evidence_digest.add(evidence_id)
                        counters["called_concrete_alleles_projected"] += 1
                        if not complete:
                            counters["observed_partial_source_alleles_projected"] += 1

                    uncalled = len(alts) - len(called_indices)
                    if uncalled > 0:
                        counters["uncalled_alt_definitions"] += uncalled
                callable_merger.flush()
        except (OSError, ValueError) as exc:
            if isinstance(exc, InputValidationError):
                raise
            raise InputValidationError(f"VCF full scan failed: {exc}") from exc
        finally:
            projected.close()
            source.close()
            reference.close()

        try:
            pysam.tabix_index(str(projected_path), preset="vcf", force=True)
        except OSError as exc:
            raise InputValidationError(f"failed to index projected VCF: {exc}") from exc

        summary = {
            "status": "PASS",
            "input_sha256": input_sha256,
            "input_basename": config.input_vcf.name,
            "accepted_evidence_multiset": accepted_evidence_digest.as_dict(),
            "counters": dict(sorted(counters.items())),
            "outputs": {
                "projected.vcf.gz": sha256_file(projected_path),
                "projected.vcf.gz.tbi": sha256_file(Path(f"{projected_path}.tbi")),
                "reference_evidence.parquet": directory_digest(reference_path),
                "non_core_variants.parquet": directory_digest(non_core_path),
                "callable_regions.parquet": directory_digest(callable_regions_path),
            },
        }
        atomic_write_json(staging / "ingest_summary.json", summary)
    return summary
