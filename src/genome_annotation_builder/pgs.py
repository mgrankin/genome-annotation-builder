"""Checksum-pinned PGS Catalog ingestion and evidence-preserving scoring."""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import date
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, TextIO

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .arrow_schemas import (
    pgs_locus_audit_schema,
    pgs_locus_resource_schema,
    pgs_model_audit_schema,
    prs_schema,
)
from .config import BuildConfig, PgsModelConfig
from .errors import InputValidationError, ResourceError
from .target_callability import CATALOG_COLUMNS
from .target_union import merge_target_catalogs
from .util import (
    atomic_output_directory,
    atomic_write_json,
    canonical_json_bytes,
    sha256_file,
)
from .vcf_ingest import canonical_core_chrom

PGS_CATALOG_SOURCE = "PGS Catalog"
SUPPORTED_FORMAT_VERSION = "2.0"
REQUIRED_SCORING_COLUMNS = {
    "effect_allele",
    "other_allele",
    "effect_weight",
    "hm_chr",
    "hm_pos",
}
CALIBRATION_COLUMNS = ("score", "percentile")
CALIBRATION_METADATA = {
    "pgs_id",
    "scoring_file_sha256",
    "reference_population",
    "calibration_source",
    "calibration_source_version",
    "calibration_date",
    "max_missing_weight_fraction",
}
DNA_COMPLEMENT = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A"})
PALINDROMIC_SNV_PAIRS = {frozenset({"A", "T"}), frozenset({"C", "G"})}
MAX_EMBEDDED_CONTRIBUTING_VARIANT_IDS = 10_000


@contextmanager
def _text_reader(path: Path) -> Iterator[TextIO]:
    raw = path.open("rb")
    try:
        magic = raw.read(2)
        raw.seek(0)
        binary: io.BufferedIOBase | gzip.GzipFile
        binary = gzip.GzipFile(fileobj=raw) if magic == b"\x1f\x8b" else raw
        text = io.TextIOWrapper(binary, encoding="utf-8", newline="")
        try:
            yield text
        finally:
            text.detach()
            if binary is not raw:
                binary.close()
    except UnicodeDecodeError as exc:
        raise ResourceError(f"resource is not UTF-8 text: {path.name}") from exc
    finally:
        raw.close()


def _reference_contigs(fasta: pysam.FastaFile) -> dict[str, tuple[str, int]]:
    output: dict[str, tuple[str, int]] = {}
    for name, length in zip(fasta.references, fasta.lengths, strict=True):
        canonical = canonical_core_chrom(name)
        if canonical is None:
            continue
        if canonical in output:
            raise InputValidationError(
                f"reference contains ambiguous core-contig aliases for {canonical}"
            )
        output[canonical] = (name, int(length))
    return output


def _metadata_and_reader(handle: TextIO, label: str) -> tuple[dict[str, str], csv.DictReader]:
    metadata: dict[str, str] = {}
    header_line: str | None = None
    while True:
        line = handle.readline()
        if line == "":
            break
        stripped = line.rstrip("\r\n")
        if stripped.startswith("#"):
            payload = stripped.lstrip("#")
            if "=" in payload:
                key, value = payload.split("=", 1)
                key = key.strip()
                if key in metadata:
                    raise ResourceError(f"duplicate metadata key {key!r} in {label}")
                metadata[key] = value.strip()
            continue
        if not stripped:
            continue
        header_line = stripped
        break
    if header_line is None:
        raise ResourceError(f"missing tabular header in {label}")
    names = header_line.split("\t")
    if len(names) != len(set(names)):
        raise ResourceError(f"duplicate tabular column in {label}")
    return metadata, csv.DictReader(handle, fieldnames=names, delimiter="\t")


def _truthy(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes", "y"}


def _false_qc(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"false", "0", "no", "n"}


def _cell(row: dict[str | None, str | list[str] | None], key: str) -> str:
    value = row.get(key)
    return str(value or "").strip()


def _nullable(value: str) -> str | None:
    return value or None


def _parse_weight(value: str, pgs_id: str, row_number: int) -> tuple[float, Decimal]:
    try:
        decimal_value = Decimal(value)
    except InvalidOperation as exc:
        raise ResourceError(
            f"{pgs_id} row {row_number} has invalid effect_weight {value!r}"
        ) from exc
    if not decimal_value.is_finite():
        raise ResourceError(f"{pgs_id} row {row_number} has non-finite effect_weight")
    floating = float(decimal_value)
    if not math.isfinite(floating):
        raise ResourceError(f"{pgs_id} row {row_number} effect_weight exceeds float64")
    return floating, decimal_value


def _parse_model(
    model: PgsModelConfig,
    reference_fasta: Path,
    *,
    locus_writer: pq.ParquetWriter | None = None,
    target_writer: csv.writer | None = None,
) -> dict[str, Any]:
    actual_sha256 = sha256_file(model.scoring_file)
    if actual_sha256 != model.expected_sha256:
        raise ResourceError(
            f"{model.pgs_id} scoring file SHA-256 mismatch: "
            f"expected {model.expected_sha256}, got {actual_sha256}"
        )
    fasta = pysam.FastaFile(str(reference_fasta))
    contigs = _reference_contigs(fasta)
    row_buffer: list[dict[str, Any]] = []
    unsupported = Counter[str]()
    total_weight = Decimal(0)
    supported_weight = Decimal(0)
    rows = 0
    supported_rows = 0
    strand_flipped_rows = 0
    try:
        with _text_reader(model.scoring_file) as handle:
            metadata, reader = _metadata_and_reader(handle, model.scoring_file.name)
            fields = tuple(reader.fieldnames or ())
            missing = sorted(REQUIRED_SCORING_COLUMNS - set(fields))
            if missing:
                raise ResourceError(
                    f"{model.pgs_id} scoring file lacks required columns: {', '.join(missing)}"
                )
            if metadata.get("format_version") != SUPPORTED_FORMAT_VERSION:
                raise ResourceError(
                    f"{model.pgs_id} requires scoring format {SUPPORTED_FORMAT_VERSION}, "
                    f"got {metadata.get('format_version')!r}"
                )
            if metadata.get("pgs_id", "").upper() != model.pgs_id:
                raise ResourceError(f"configured {model.pgs_id} does not match scoring-file pgs_id")
            if metadata.get("HmPOS_build") != "GRCh38":
                raise ResourceError(
                    f"{model.pgs_id} must be an official harmonized GRCh38 scoring file"
                )
            declared_rows_text = metadata.get("variants_number")
            try:
                declared_rows = int(str(declared_rows_text))
            except (TypeError, ValueError) as exc:
                raise ResourceError(f"{model.pgs_id} has invalid variants_number metadata") from exc
            for row_number, raw in enumerate(reader, start=1):
                if None in raw or any(value is None for value in raw.values()):
                    raise ResourceError(
                        f"{model.pgs_id} scoring row {row_number} has the wrong column count"
                    )
                raw_row = {str(key): str(value) for key, value in raw.items()}
                source_effect_allele = _cell(raw, "effect_allele").upper()
                source_other_allele = _cell(raw, "other_allele").upper()
                effect_allele = source_effect_allele
                other_allele = source_other_allele
                effect_weight, decimal_weight = _parse_weight(
                    _cell(raw, "effect_weight"), model.pgs_id, row_number
                )
                total_weight += abs(decimal_weight)
                reason: str | None = None
                chrom: str | None = None
                pos: int | None = None
                ref: str | None = None
                alt: str | None = None
                variant_id: str | None = None
                effect_is_reference: bool | None = None
                strand_flipped: bool | None = None
                target_id: str | None = None

                if _truthy(_cell(raw, "is_haplotype")):
                    reason = "haplotype_model_row"
                elif _truthy(_cell(raw, "is_diplotype")):
                    reason = "diplotype_model_row"
                elif _truthy(_cell(raw, "is_dominant")):
                    reason = "dominant_model_row"
                elif _truthy(_cell(raw, "is_recessive")):
                    reason = "recessive_model_row"
                elif any(
                    _cell(raw, key)
                    for key in ("dosage_0_weight", "dosage_1_weight", "dosage_2_weight")
                ):
                    reason = "non_additive_dosage_weights"
                elif _truthy(_cell(raw, "hm_inferOtherAllele")):
                    reason = "inferred_other_allele"
                elif _false_qc(_cell(raw, "hm_match_chr")) or _false_qc(_cell(raw, "hm_match_pos")):
                    reason = "harmonization_qc_failed"
                elif (
                    len(effect_allele) != 1
                    or len(other_allele) != 1
                    or effect_allele not in "ACGT"
                    or other_allele not in "ACGT"
                    or effect_allele == other_allele
                ):
                    reason = "unsupported_non_biallelic_snv_alleles"
                elif frozenset({effect_allele, other_allele}) in PALINDROMIC_SNV_PAIRS:
                    reason = "strand_ambiguous_palindromic_snv"
                else:
                    chrom = canonical_core_chrom(_cell(raw, "hm_chr"))
                    try:
                        pos = int(_cell(raw, "hm_pos"))
                    except ValueError:
                        reason = "invalid_harmonized_coordinate"
                    if reason is None and (chrom is None or chrom not in contigs):
                        reason = "non_primary_or_missing_harmonized_contig"
                    if reason is None:
                        assert chrom is not None and pos is not None
                        reference_name, reference_length = contigs[chrom]
                        if not 1 <= pos <= reference_length:
                            reason = "harmonized_coordinate_outside_reference"
                        else:
                            observed_ref = fasta.fetch(reference_name, pos - 1, pos).upper()
                            if effect_allele == observed_ref:
                                ref, alt, effect_is_reference = (
                                    effect_allele,
                                    other_allele,
                                    True,
                                )
                            elif other_allele == observed_ref:
                                ref, alt, effect_is_reference = (
                                    other_allele,
                                    effect_allele,
                                    False,
                                )
                            else:
                                complemented_effect = effect_allele.translate(DNA_COMPLEMENT)
                                complemented_other = other_allele.translate(DNA_COMPLEMENT)
                                if complemented_effect == observed_ref:
                                    effect_allele = complemented_effect
                                    other_allele = complemented_other
                                    ref, alt, effect_is_reference = (
                                        effect_allele,
                                        other_allele,
                                        True,
                                    )
                                    strand_flipped = True
                                elif complemented_other == observed_ref:
                                    effect_allele = complemented_effect
                                    other_allele = complemented_other
                                    ref, alt, effect_is_reference = (
                                        other_allele,
                                        effect_allele,
                                        False,
                                    )
                                    strand_flipped = True
                                else:
                                    reason = "alleles_do_not_match_grch38_reference"
                            if reason is None and strand_flipped is None:
                                strand_flipped = False
                    if reason is None:
                        assert chrom is not None and pos is not None and ref and alt
                        variant_id = f"{chrom}:{pos}:{ref}:{alt}"
                        target_id = f"pgs:{model.pgs_id}:row:{row_number:09d}"
                        supported_rows += 1
                        supported_weight += abs(decimal_weight)
                        if strand_flipped:
                            strand_flipped_rows += 1

                supported = reason is None
                if reason is not None:
                    unsupported[reason] += 1
                row_identity = {
                    "pgs_id": model.pgs_id,
                    "scoring_row_number": row_number,
                    "columns": raw_row,
                }
                locus = {
                    "pgs_id": model.pgs_id,
                    "scoring_row_number": row_number,
                    "scoring_row_digest": hashlib.sha256(
                        canonical_json_bytes(row_identity)
                    ).hexdigest(),
                    "scoring_file_sha256": actual_sha256,
                    "supported": supported,
                    "unsupported_reason": reason,
                    "target_id": target_id,
                    "chrom": chrom if supported else None,
                    "pos": pos if supported else None,
                    "ref": ref,
                    "alt": alt,
                    "variant_id": variant_id,
                    "effect_allele": effect_allele,
                    "other_allele": _nullable(other_allele),
                    "source_effect_allele": source_effect_allele,
                    "source_other_allele": _nullable(source_other_allele),
                    "strand_flipped": strand_flipped,
                    "effect_weight": effect_weight,
                    "effect_is_reference": effect_is_reference,
                    "source_rsid": _nullable(_cell(raw, "rsID")),
                    "source_locus_name": _nullable(_cell(raw, "locus_name")),
                    "harmonization_source": _nullable(_cell(raw, "hm_source")),
                    "harmonization_rsid": _nullable(_cell(raw, "hm_rsID")),
                }
                if locus_writer is not None:
                    row_buffer.append(locus)
                    if len(row_buffer) >= 20_000:
                        locus_writer.write_table(
                            pa.Table.from_pylist(row_buffer, schema=pgs_locus_resource_schema()),
                            row_group_size=20_000,
                        )
                        row_buffer.clear()
                if supported and target_writer is not None:
                    target_writer.writerow(
                        [
                            target_id,
                            chrom,
                            pos,
                            ref,
                            alt,
                            "",
                            PGS_CATALOG_SOURCE,
                            f"{model.pgs_id}@{model.source_version}",
                            "target_panel",
                        ]
                    )
                rows += 1
    finally:
        fasta.close()
    if locus_writer is not None and row_buffer:
        locus_writer.write_table(
            pa.Table.from_pylist(row_buffer, schema=pgs_locus_resource_schema()),
            row_group_size=20_000,
        )
    if rows != declared_rows:
        raise ResourceError(
            f"{model.pgs_id} variants_number={declared_rows}, but parsed {rows} rows"
        )
    if rows > 2_147_483_647:
        raise ResourceError(f"{model.pgs_id} has too many scoring rows for the typed schema")
    if total_weight == 0:
        raise ResourceError(f"{model.pgs_id} has zero total absolute weight")
    if supported_rows == 0:
        raise ResourceError(f"{model.pgs_id} has no supported biallelic GRCh38 SNV scoring loci")
    return {
        "pgs_id": model.pgs_id,
        "scoring_file_sha256": actual_sha256,
        "scoring_file_bytes": model.scoring_file.stat().st_size,
        "format_version": metadata["format_version"],
        "catalog_pgs_name": metadata.get("pgs_name"),
        "catalog_trait_reported": metadata.get("trait_reported"),
        "catalog_trait_mapped": metadata.get("trait_mapped"),
        "catalog_trait_efo": metadata.get("trait_efo"),
        "harmonized_build": metadata["HmPOS_build"],
        "harmonized_date": metadata.get("HmPOS_date"),
        "expected_loci": rows,
        "supported_loci": supported_rows,
        "strand_flipped_loci": strand_flipped_rows,
        "unsupported_loci": rows - supported_rows,
        "unsupported_reason_counts": dict(sorted(unsupported.items())),
        "expected_absolute_weight": float(total_weight),
        "supported_absolute_weight": float(supported_weight),
        "source_version": model.source_version,
        "training_date": model.training_date.isoformat(),
        "trait": model.trait,
        "model_type": model.model_type,
        "license": model.license,
        "redistribution_allowed": model.redistribution_allowed,
        "max_missing_weight_fraction": model.max_missing_weight_fraction,
    }


def _parse_calibration(model: PgsModelConfig, *, include_points: bool) -> dict[str, Any] | None:
    configured = model.calibration
    if configured is None:
        return None
    if not configured.path.is_file():
        raise ResourceError(f"{model.pgs_id} calibration file not found")
    actual = sha256_file(configured.path)
    if actual != configured.expected_sha256:
        raise ResourceError(
            f"{model.pgs_id} calibration SHA-256 mismatch: "
            f"expected {configured.expected_sha256}, got {actual}"
        )
    with _text_reader(configured.path) as handle:
        metadata, reader = _metadata_and_reader(handle, configured.path.name)
        if tuple(reader.fieldnames or ()) != CALIBRATION_COLUMNS:
            raise ResourceError(
                f"{model.pgs_id} calibration header must be exactly: "
                + "\t".join(CALIBRATION_COLUMNS)
            )
        missing = sorted(CALIBRATION_METADATA - set(metadata))
        if missing:
            raise ResourceError(
                f"{model.pgs_id} calibration metadata missing: {', '.join(missing)}"
            )
        if metadata["pgs_id"].upper() != model.pgs_id:
            raise ResourceError(f"{model.pgs_id} calibration pgs_id mismatch")
        if metadata["scoring_file_sha256"].lower() != model.expected_sha256:
            raise ResourceError(f"{model.pgs_id} calibration scoring-file digest mismatch")
        try:
            calibration_date = date.fromisoformat(metadata["calibration_date"])
        except ValueError as exc:
            raise ResourceError(f"{model.pgs_id} calibration_date must be YYYY-MM-DD") from exc
        try:
            maximum_missing = float(metadata["max_missing_weight_fraction"])
        except ValueError as exc:
            raise ResourceError(
                f"{model.pgs_id} calibration max_missing_weight_fraction is invalid"
            ) from exc
        if not math.isfinite(maximum_missing) or not 0 <= maximum_missing <= 1:
            raise ResourceError(
                f"{model.pgs_id} calibration max_missing_weight_fraction must be 0..1"
            )
        scores: list[float] = []
        percentiles: list[float] = []
        for row_number, row in enumerate(reader, start=1):
            if None in row or any(value is None for value in row.values()):
                raise ResourceError(
                    f"{model.pgs_id} calibration row {row_number} has wrong column count"
                )
            try:
                score = float(str(row["score"]).strip())
                percentile = float(str(row["percentile"]).strip())
            except ValueError as exc:
                raise ResourceError(
                    f"{model.pgs_id} calibration row {row_number} is not numeric"
                ) from exc
            if not math.isfinite(score) or not math.isfinite(percentile):
                raise ResourceError(f"{model.pgs_id} calibration row {row_number} is non-finite")
            if scores and score <= scores[-1]:
                raise ResourceError(f"{model.pgs_id} calibration scores are not increasing")
            if not 0 <= percentile <= 100 or (percentiles and percentile < percentiles[-1]):
                raise ResourceError(
                    f"{model.pgs_id} calibration percentiles must be monotone in 0..100"
                )
            scores.append(score)
            percentiles.append(percentile)
    if len(scores) < 2:
        raise ResourceError(f"{model.pgs_id} calibration requires at least two points")
    result: dict[str, Any] = {
        "sha256": actual,
        "bytes": configured.path.stat().st_size,
        "reference_population": metadata["reference_population"],
        "source": metadata["calibration_source"],
        "source_version": metadata["calibration_source_version"],
        "date": calibration_date.isoformat(),
        "max_missing_weight_fraction": maximum_missing,
        "points": len(scores),
    }
    if include_points:
        result["scores"] = scores
        result["percentiles"] = percentiles
    return result


def verify_pgs_resources(config: BuildConfig) -> dict[str, Any]:
    if not config.pgs_enabled:
        return {"status": "not_enabled", "models": []}
    if not config.reference_fasta.is_file():
        raise ResourceError("PGS verification requires the configured GRCh38 FASTA")
    models: list[dict[str, Any]] = []
    for model in config.pgs_models:
        if not model.scoring_file.is_file():
            raise ResourceError(f"{model.pgs_id} scoring file not found")
        summary = _parse_model(model, config.reference_fasta)
        summary["calibration"] = _parse_calibration(model, include_points=False)
        models.append(summary)
    return {
        "status": "verified",
        "source": PGS_CATALOG_SOURCE,
        "models": models,
        "model_count": len(models),
        "scoring_file_format": SUPPORTED_FORMAT_VERSION,
    }


def _merge_targets(
    config: BuildConfig, pgs_targets: Path, combined_targets: Path
) -> tuple[int, int]:
    sources: list[tuple[str, Path]] = []
    if config.target_catalog_enabled:
        sources.append(("configured base", config.target_catalog_path))
    sources.append(("generated PGS", pgs_targets))
    counts = merge_target_catalogs(sources, combined_targets)
    return counts.get("configured base", 0), counts["generated PGS"]


def prepare_pgs_resources(config: BuildConfig, destination: str | Path) -> dict[str, Any]:
    if not config.pgs_enabled:
        raise ResourceError("PGS preparation was invoked while annotation.pgs is disabled")
    destination = Path(destination)
    with atomic_output_directory(destination) as staging:
        loci_path = staging / "pgs_loci.parquet"
        pgs_targets = staging / "pgs_targets.tsv"
        model_summaries: list[dict[str, Any]] = []
        with (
            pq.ParquetWriter(
                loci_path,
                pgs_locus_resource_schema(),
                compression="zstd",
                compression_level=9,
                use_dictionary=True,
                write_statistics=True,
                version="2.6",
            ) as locus_writer,
            pgs_targets.open("w", encoding="utf-8", newline="") as target_handle,
        ):
            target_writer = csv.writer(target_handle, delimiter="\t", lineterminator="\n")
            target_writer.writerow(CATALOG_COLUMNS)
            for model in config.pgs_models:
                summary = _parse_model(
                    model,
                    config.reference_fasta,
                    locus_writer=locus_writer,
                    target_writer=target_writer,
                )
                summary["calibration"] = _parse_calibration(model, include_points=False)
                model_summaries.append(summary)

        combined_targets = staging / "combined_targets.tsv"
        base_rows, pgs_rows = _merge_targets(config, pgs_targets, combined_targets)
        combined_sha256 = sha256_file(combined_targets)
        summary = {
            "status": "PASS",
            "source": PGS_CATALOG_SOURCE,
            "models": model_summaries,
            "model_count": len(model_summaries),
            "base_target_rows": base_rows,
            "pgs_target_rows": pgs_rows,
            "combined_target_rows": base_rows + pgs_rows,
            "combined_target_catalog_sha256": combined_sha256,
            "outputs": {
                "pgs_loci.parquet": sha256_file(loci_path),
                "pgs_targets.tsv": sha256_file(pgs_targets),
                "combined_targets.tsv": combined_sha256,
            },
        }
        atomic_write_json(staging / "pgs_resource_summary.json", summary)
    return summary


def _interpolate_percentile(calibration: dict[str, Any], score: float) -> float:
    scores = calibration["scores"]
    percentiles = calibration["percentiles"]
    if score <= scores[0]:
        return float(percentiles[0])
    if score >= scores[-1]:
        return float(percentiles[-1])
    right = bisect.bisect_right(scores, score)
    left = right - 1
    fraction = (score - scores[left]) / (scores[right] - scores[left])
    return float(percentiles[left] + fraction * (percentiles[right] - percentiles[left]))


def _resolve_locus(row: dict[str, Any]) -> dict[str, Any]:
    base = {
        "pgs_id": row["pgs_id"],
        "scoring_row_number": row["scoring_row_number"],
        "scoring_row_digest": row["scoring_row_digest"],
        "target_id": row.get("target_id"),
        "callability_id": row.get("callability_id"),
        "target_catalog_sha256": row.get("target_catalog_sha256"),
        "target_catalog_row_digest": row.get("target_catalog_row_digest"),
        "chrom": row.get("chrom"),
        "pos": row.get("pos"),
        "ref": row.get("ref"),
        "alt": row.get("alt"),
        "variant_id": row.get("variant_id"),
        "effect_allele": row["effect_allele"],
        "other_allele": row.get("other_allele"),
        "source_effect_allele": row["source_effect_allele"],
        "source_other_allele": row.get("source_other_allele"),
        "strand_flipped": row.get("strand_flipped"),
        "effect_weight": row["effect_weight"],
        "effect_is_reference": row.get("effect_is_reference"),
        "resolution": row.get("resolution"),
        "dosage": None,
        "observed_ploidy": row.get("observed_ploidy"),
        "contribution": None,
        "observed_variant_ids": list(row.get("observed_variant_ids") or []),
        "depth": row.get("depth"),
        "genotype_quality": row.get("genotype_quality"),
    }
    if not row["supported"]:
        return {
            **base,
            "status": "unsupported",
            "reason": str(row.get("unsupported_reason") or "unsupported scoring row"),
        }
    if row.get("callability_id") is None:
        raise InputValidationError(
            f"PGS target evidence is missing for {row['pgs_id']} row {row['scoring_row_number']}"
        )
    resolution = row.get("resolution")
    if resolution == "non_reference_observed":
        observed_ids = list(row.get("observed_variant_ids") or [])
        if observed_ids != [row["variant_id"]]:
            return {
                **base,
                "status": "ambiguous_other_variant",
                "reason": "the target allele co-occurs with or is obscured by another overlapping allele",
            }
        if row.get("callable") is not True or row.get("call_confidence") not in {
            "high",
            "medium",
        }:
            return {
                **base,
                "status": "low_confidence",
                "reason": "the exact carried allele does not pass configured genotype confidence",
            }
        gt = list(row.get("variant_gt") or [])
        if len(gt) not in {1, 2} or any(value not in {0, 1} for value in gt):
            return {
                **base,
                "status": "ambiguous_genotype",
                "reason": "the exact carried allele has unsupported or incomplete ploidy/genotype",
            }
        ploidy = len(gt)
        alt_dosage = sum(value == 1 for value in gt)
        dosage = ploidy - alt_dosage if row["effect_is_reference"] else alt_dosage
        return {
            **base,
            "status": "alternate_called",
            "dosage": dosage,
            "observed_ploidy": ploidy,
            "contribution": float(dosage) * float(row["effect_weight"]),
            "reason": "exact carried target allele with complete passing genotype evidence",
        }
    if resolution == "reference_observed":
        ploidy = row.get("observed_ploidy")
        if (
            row.get("reference_observed") is not True
            or row.get("callable") is not True
            or ploidy not in {1, 2}
        ):
            return {
                **base,
                "status": "low_confidence",
                "reason": "reference state lacks complete passing ploidy/callability evidence",
            }
        dosage = int(ploidy) if row["effect_is_reference"] else 0
        return {
            **base,
            "status": "reference_called",
            "dosage": dosage,
            "contribution": float(dosage) * float(row["effect_weight"]),
            "reason": "explicit complete reference evidence covers the exact scoring locus",
        }
    if resolution == "overlapping_variant_observed":
        return {
            **base,
            "status": "ambiguous_other_variant",
            "reason": "a different carried allele overlaps the scoring locus",
        }
    if resolution == "low_confidence_reference":
        return {
            **base,
            "status": "low_confidence",
            "reason": "reference evidence exists but fails configured confidence thresholds",
        }
    if resolution == "no_call":
        return {
            **base,
            "status": "missing_no_call",
            "reason": "an explicit no-call supersedes reference or interval evidence",
        }
    if resolution in {"callable_without_reference", "uncallable_without_reference"}:
        return {
            **base,
            "status": "missing_interval_only",
            "reason": "interval callability is not an explicit genotype and is not scored",
        }
    return {
        **base,
        "status": "missing_unknown",
        "reason": "no applicable explicit variant or reference genotype evidence was retained",
    }


def _register_parquet(connection: duckdb.DuckDBPyConnection, name: str, path: Path) -> None:
    source = path / "**" / "*.parquet" if path.is_dir() else path
    escaped = str(source).replace("'", "''")
    connection.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{escaped}')")


def score_pgs(
    config: BuildConfig,
    *,
    resource_dir: str | Path,
    target_callability_dir: str | Path,
    variant_dir: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    if not config.pgs_enabled:
        raise InputValidationError("PGS scoring was invoked while annotation.pgs is disabled")
    resource_dir = Path(resource_dir)
    target_callability_dir = Path(target_callability_dir)
    variant_dir = Path(variant_dir)
    destination = Path(destination)
    resource_summary = json.loads(
        (resource_dir / "pgs_resource_summary.json").read_text(encoding="utf-8")
    )
    target_catalog_sha256 = resource_summary["combined_target_catalog_sha256"]
    model_resources = {row["pgs_id"]: row for row in resource_summary["models"]}
    model_configs = {model.pgs_id: model for model in config.pgs_models}
    calibrations = {
        model.pgs_id: _parse_calibration(model, include_points=True) for model in config.pgs_models
    }

    with atomic_output_directory(destination) as staging:
        locus_audit_path = staging / "pgs_locus_audit.parquet"
        model_audit_path = staging / "pgs_model_audit.parquet"
        model_aggregates: dict[str, dict[str, Any]] = {}
        connection = duckdb.connect(str(staging / ".pgs-score.duckdb"))
        try:
            connection.execute(f"SET temp_directory='{str(staging).replace(chr(39), chr(39) * 2)}'")
            _register_parquet(connection, "pgs_loci", resource_dir / "pgs_loci.parquet")
            _register_parquet(
                connection,
                "target_evidence",
                target_callability_dir / "target_callability_evidence.parquet",
            )
            _register_parquet(connection, "variants", variant_dir / "variants.parquet")
            duplicate = connection.execute(
                "SELECT target_id FROM target_evidence GROUP BY target_id "
                "HAVING count(*) != 1 LIMIT 1"
            ).fetchone()
            if duplicate:
                raise InputValidationError(
                    f"target evidence is not unique for target_id {duplicate[0]}"
                )
            cursor = connection.execute(
                """
                SELECT l.pgs_id, l.scoring_row_number, l.scoring_row_digest,
                       l.supported, l.unsupported_reason, l.target_id,
                       l.chrom, l.pos, l.ref, l.alt, l.variant_id,
                       l.effect_allele, l.other_allele, l.effect_weight,
                       l.source_effect_allele, l.source_other_allele,
                       l.strand_flipped, l.effect_is_reference,
                       e.callability_id,
                       e.catalog_sha256 AS target_catalog_sha256,
                       e.catalog_row_digest AS target_catalog_row_digest,
                       e.resolution, e.reference_observed, e.callable,
                       e.call_confidence, e.depth, e.genotype_quality,
                       e.observed_ploidy, e.observed_variant_ids,
                       v.genotype.gt AS variant_gt
                FROM pgs_loci l
                LEFT JOIN target_evidence e ON l.target_id = e.target_id
                LEFT JOIN variants v ON l.variant_id = v.variant_id
                ORDER BY l.pgs_id, l.scoring_row_number
                """
            )
            names = [description[0] for description in cursor.description]
            with pq.ParquetWriter(
                locus_audit_path,
                pgs_locus_audit_schema(),
                compression="zstd",
                compression_level=9,
                use_dictionary=True,
                write_statistics=True,
                version="2.6",
            ) as locus_writer:
                while batch := cursor.fetchmany(20_000):
                    audit_batch: list[dict[str, Any]] = []
                    for values in batch:
                        source = dict(zip(names, values, strict=True))
                        resolved = _resolve_locus(source)
                        if (
                            resolved["target_catalog_sha256"] is not None
                            and resolved["target_catalog_sha256"] != target_catalog_sha256
                        ):
                            raise InputValidationError(
                                "PGS locus points to a different target-catalogue digest"
                            )
                        audit_batch.append(resolved)
                    locus_writer.write_table(
                        pa.Table.from_pylist(audit_batch, schema=pgs_locus_audit_schema()),
                        row_group_size=20_000,
                    )

            # Aggregate from the on-disk audit instead of retaining every scoring
            # locus in Python. Peak memory is therefore independent of locus count.
            _register_parquet(connection, "pgs_locus_audit", locus_audit_path)
            aggregate_cursor = connection.execute(
                """
                SELECT pgs_id,
                  count(*) AS expected_loci,
                  count(*) FILTER (WHERE status != 'unsupported') AS supported_loci,
                  count(*) FILTER (WHERE strand_flipped) AS strand_flipped_loci,
                  count(*) FILTER (WHERE status = 'alternate_called')
                    AS alternate_called_loci,
                  count(*) FILTER (WHERE status = 'reference_called')
                    AS reference_called_loci,
                  count(*) FILTER (WHERE starts_with(status, 'missing_')) AS missing_loci,
                  count(*) FILTER (WHERE starts_with(status, 'ambiguous_'))
                    AS ambiguous_loci,
                  count(*) FILTER (WHERE status = 'low_confidence')
                    AS low_confidence_loci,
                  count(*) FILTER (WHERE status = 'unsupported') AS unsupported_loci,
                  coalesce(fsum(abs(effect_weight)), 0.0) AS expected_absolute_weight,
                  coalesce(
                    fsum(abs(effect_weight)) FILTER (WHERE contribution IS NOT NULL), 0.0
                  ) AS used_absolute_weight,
                  coalesce(fsum(contribution), 0.0) AS partial_raw_score,
                  count(DISTINCT variant_id) FILTER (
                    WHERE status = 'alternate_called' AND contribution IS NOT NULL
                      AND variant_id IS NOT NULL
                  ) AS contributing_variant_count
                FROM pgs_locus_audit
                GROUP BY pgs_id
                ORDER BY pgs_id
                """
            )
            aggregate_names = [description[0] for description in aggregate_cursor.description]
            for values in aggregate_cursor.fetchall():
                aggregate = dict(zip(aggregate_names, values, strict=True))
                pgs_id = str(aggregate["pgs_id"])
                contributing_count = int(aggregate["contributing_variant_count"])
                if contributing_count <= MAX_EMBEDDED_CONTRIBUTING_VARIANT_IDS:
                    contributing = connection.execute(
                        """
                        SELECT list_sort(list(DISTINCT variant_id))
                        FROM pgs_locus_audit
                        WHERE pgs_id = ? AND status = 'alternate_called'
                          AND contribution IS NOT NULL AND variant_id IS NOT NULL
                        """,
                        [pgs_id],
                    ).fetchone()[0]
                    aggregate["contributing_variant_ids"] = list(contributing or [])
                    aggregate["contributing_variant_ids_complete"] = True
                else:
                    aggregate["contributing_variant_ids"] = None
                    aggregate["contributing_variant_ids_complete"] = False
                model_aggregates[pgs_id] = aggregate
        finally:
            connection.close()
            database_path = staging / ".pgs-score.duckdb"
            if database_path.exists():
                database_path.unlink()

        model_audits: list[dict[str, Any]] = []
        normative_rows: list[dict[str, Any]] = []
        status_counts: Counter[str] = Counter()
        expected_models = set(model_configs)
        if set(model_aggregates) != expected_models:
            missing = sorted(expected_models - set(model_aggregates))
            unexpected = sorted(set(model_aggregates) - expected_models)
            raise InputValidationError(
                "PGS locus/model aggregation mismatch: "
                f"missing={missing!r}, unexpected={unexpected!r}"
            )
        for pgs_id in sorted(model_aggregates):
            aggregate = model_aggregates[pgs_id]
            model = model_configs[pgs_id]
            resource = model_resources[pgs_id]
            expected_loci = int(aggregate["expected_loci"])
            supported_loci = int(aggregate["supported_loci"])
            strand_flipped_loci = int(aggregate["strand_flipped_loci"])
            expected_absolute_weight = float(aggregate["expected_absolute_weight"])
            used_absolute_weight = float(aggregate["used_absolute_weight"])
            if (
                expected_loci != int(resource["expected_loci"])
                or supported_loci != int(resource["supported_loci"])
                or strand_flipped_loci != int(resource["strand_flipped_loci"])
                or not math.isclose(
                    expected_absolute_weight,
                    float(resource["expected_absolute_weight"]),
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
            ):
                raise InputValidationError(
                    f"{pgs_id} scored locus aggregate disagrees with its pinned resource summary"
                )
            missing_absolute_weight = max(0.0, expected_absolute_weight - used_absolute_weight)
            weight_completeness = min(1.0, used_absolute_weight / expected_absolute_weight)
            missing_weight_fraction = 1.0 - weight_completeness
            partial_raw_score = float(aggregate["partial_raw_score"])
            score_valid = missing_weight_fraction <= model.max_missing_weight_fraction + 1e-12
            calibration = calibrations[pgs_id]
            calibrated = False
            percentile: float | None = None
            if (
                score_valid
                and calibration is not None
                and missing_weight_fraction
                <= float(calibration["max_missing_weight_fraction"]) + 1e-12
            ):
                calibrated = True
                percentile = _interpolate_percentile(calibration, partial_raw_score)
            if not score_valid:
                publication_status = "incomplete"
                reason = (
                    "missing/ambiguous/low-confidence weight exceeds the model's declared threshold"
                )
            elif calibration is None:
                publication_status = "valid_uncalibrated"
                reason = "raw score is complete enough, but no population calibration was supplied"
            elif not calibrated:
                publication_status = "calibration_inapplicable"
                reason = "raw score passes the model threshold but exceeds the calibration's missingness scope"
            else:
                publication_status = "published"
                reason = "raw score and percentile pass locus, weight, and calibration gates"
            status_counts[publication_status] += 1
            contributing = aggregate["contributing_variant_ids"]
            audit = {
                "pgs_id": pgs_id,
                "trait": model.trait,
                "scoring_file_sha256": model.expected_sha256,
                "scoring_format_version": resource["format_version"],
                "source_version": model.source_version,
                "training_date": model.training_date,
                "model_type": model.model_type,
                "license": model.license,
                "redistribution_allowed": model.redistribution_allowed,
                "expected_loci": expected_loci,
                "supported_loci": supported_loci,
                "strand_flipped_loci": strand_flipped_loci,
                "alternate_called_loci": int(aggregate["alternate_called_loci"]),
                "reference_called_loci": int(aggregate["reference_called_loci"]),
                "missing_loci": int(aggregate["missing_loci"]),
                "ambiguous_loci": int(aggregate["ambiguous_loci"]),
                "low_confidence_loci": int(aggregate["low_confidence_loci"]),
                "unsupported_loci": int(aggregate["unsupported_loci"]),
                "expected_absolute_weight": expected_absolute_weight,
                "used_absolute_weight": used_absolute_weight,
                "missing_absolute_weight": missing_absolute_weight,
                "weight_completeness": weight_completeness,
                "max_missing_weight_fraction": model.max_missing_weight_fraction,
                "partial_raw_score": partial_raw_score,
                "score_valid": score_valid,
                "calibrated": calibrated,
                "percentile": percentile,
                "reference_population": (
                    calibration["reference_population"] if calibrated else None
                ),
                "calibration_source": calibration["source"] if calibration else None,
                "calibration_source_version": (
                    calibration["source_version"] if calibration else None
                ),
                "calibration_date": (
                    date.fromisoformat(calibration["date"]) if calibration else None
                ),
                "calibration_sha256": calibration["sha256"] if calibration else None,
                "calibration_max_missing_weight_fraction": (
                    calibration["max_missing_weight_fraction"] if calibration else None
                ),
                "publication_status": publication_status,
                "reason": reason,
                "contributing_variant_count": int(aggregate["contributing_variant_count"]),
                "contributing_variant_ids_complete": bool(
                    aggregate["contributing_variant_ids_complete"]
                ),
                "contributing_variant_ids": contributing,
                "target_catalog_sha256": target_catalog_sha256,
            }
            model_audits.append(audit)
            if publication_status == "published":
                assert percentile is not None and calibration is not None
                normative_rows.append(
                    {
                        "trait": model.trait,
                        "score_value": partial_raw_score,
                        "percentile": percentile,
                        "reference_population": calibration["reference_population"],
                        "training_source": pgs_id,
                        "training_date": model.training_date,
                        "contributing_variants": (
                            contributing if aggregate["contributing_variant_ids_complete"] else None
                        ),
                        "model_type": model.model_type,
                    }
                )

        pq.write_table(
            pa.Table.from_pylist(model_audits, schema=pgs_model_audit_schema()),
            model_audit_path,
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
            write_statistics=True,
            version="2.6",
        )
        outputs = {
            "pgs_locus_audit.parquet": sha256_file(locus_audit_path),
            "pgs_model_audit.parquet": sha256_file(model_audit_path),
        }
        if normative_rows:
            prs_path = staging / "prs.parquet"
            pq.write_table(
                pa.Table.from_pylist(normative_rows, schema=prs_schema()),
                prs_path,
                compression="zstd",
                compression_level=9,
                use_dictionary=True,
                write_statistics=True,
                version="2.6",
            )
            outputs["prs.parquet"] = sha256_file(prs_path)
        summary = {
            "status": "PASS",
            "model_count": len(model_audits),
            "published_model_count": len(normative_rows),
            "publication_status_counts": dict(sorted(status_counts.items())),
            "target_catalog_sha256": target_catalog_sha256,
            "models": [
                {
                    "pgs_id": row["pgs_id"],
                    "publication_status": row["publication_status"],
                    "score_valid": row["score_valid"],
                    "calibrated": row["calibrated"],
                    "weight_completeness": row["weight_completeness"],
                    "expected_loci": row["expected_loci"],
                    "strand_flipped_loci": row["strand_flipped_loci"],
                    "alternate_called_loci": row["alternate_called_loci"],
                    "reference_called_loci": row["reference_called_loci"],
                    "missing_loci": row["missing_loci"],
                    "ambiguous_loci": row["ambiguous_loci"],
                    "low_confidence_loci": row["low_confidence_loci"],
                    "unsupported_loci": row["unsupported_loci"],
                }
                for row in model_audits
            ],
            "outputs": outputs,
        }
        atomic_write_json(staging / "pgs_score_summary.json", summary)
    return summary
