"""Allele-aware GATK mitochondrial FilterMutectCalls portable adapter."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from .portable_model import (
    CONCRETE_ALLELE,
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    sha256_file,
    validate_rows,
    write_extension,
)


@dataclass(frozen=True)
class MitochondrialProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    caller_name: str = "GATK Mutect2 mitochondrial mode"
    caller_version: str = "4.6.2.0"
    circular_boundary_audited: bool = True


def _tuple(value: Any) -> tuple[Any, ...]:
    if value is None:
        return ()
    return value if isinstance(value, tuple) else (value,)


def _info_value(record: pysam.VariantRecord, key: str) -> Any:
    """Read an optional INFO value without asking htslib for an unknown header key."""
    if key not in record.header.info:
        return None
    return record.info.get(key)


def _integer(value: Any, label: str) -> int:
    if value is None:
        raise PortableEvidenceError(f"missing required mitochondrial {label}")
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid mitochondrial {label}") from exc


def _float(value: Any, label: str) -> float:
    if value is None:
        raise PortableEvidenceError(f"missing required mitochondrial {label}")
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid mitochondrial {label}") from exc


def _allele_value(
    record: pysam.VariantRecord, key: str, alt_index: int, *, required: bool = False
) -> Any:
    value = _info_value(record, key)
    if value is None:
        if required:
            raise PortableEvidenceError(f"missing required allele-aware INFO/{key}")
        return None
    values = _tuple(value)
    number = record.header.info[key].number
    offset = alt_index if number == "R" else alt_index - 1 if number == "A" else 0
    if offset >= len(values):
        raise PortableEvidenceError(f"INFO/{key} cardinality does not match ALT alleles")
    return values[offset]


def _format_allele_value(
    sample: pysam.libcbcf.VariantRecordSample,
    key: str,
    alt_index: int,
    *,
    number: str,
    required: bool = False,
) -> Any:
    value = sample.get(key)
    if value is None:
        if required:
            raise PortableEvidenceError(f"missing required allele-aware FORMAT/{key}")
        return None
    values = _tuple(value)
    offset = alt_index if number == "R" else alt_index - 1
    if offset >= len(values):
        raise PortableEvidenceError(f"FORMAT/{key} cardinality does not match ALT alleles")
    return values[offset]


def _strand_counts(record: pysam.VariantRecord, alt_index: int) -> tuple[int | None, int | None]:
    raw = _info_value(record, "AS_SB_TABLE")
    if raw is not None:
        # htslib may split this nominal String field at commas.  Rejoining the
        # tuple recovers GATK's REF_F,REF_R|ALT_F,ALT_R representation.
        rendered = ",".join(str(item) for item in _tuple(raw))
        allele_groups = rendered.split("|")
        if len(allele_groups) != len(record.alleles):
            raise PortableEvidenceError("INFO/AS_SB_TABLE cardinality is invalid")
        tokens = allele_groups[alt_index].split(",")
        if len(tokens) != 2:
            raise PortableEvidenceError("INFO/AS_SB_TABLE allele entry is invalid")
        return _integer(tokens[0], "forward alt count"), _integer(tokens[1], "reverse alt count")
    if len(record.alts or ()) == 1:
        sample = next(iter(record.samples.values()))
        values = _tuple(sample.get("SB"))
        if len(values) == 4:
            return _integer(values[2], "forward alt count"), _integer(
                values[3], "reverse alt count"
            )
    return None, None


def _allele_filter_status(record: pysam.VariantRecord, alt_index: int) -> list[str]:
    value = _info_value(record, "AS_FilterStatus")
    if value is None:
        return ["SITE"]
    values = _tuple(value)
    alt_count = len(record.alts or ())
    # FilterMutectCalls serializes multiple filters within one allele with
    # commas and separates allele groups with ``|``.  Because the header still
    # declares Number=A, htslib first splits those internal commas, so rebuild
    # the serialized value before applying GATK's allele delimiter.  Retain
    # support for ordinary VCF Number=A comma-separated values as well.
    rendered = ",".join(str(item) for item in values)
    if "|" in rendered:
        allele_groups = rendered.split("|")
    elif alt_count == 1:
        allele_groups = [rendered]
    else:
        allele_groups = [str(item) for item in values]
    if len(allele_groups) != alt_count:
        raise PortableEvidenceError("INFO/AS_FilterStatus cardinality does not match ALT alleles")
    statuses = [token for token in re_split_filter(allele_groups[alt_index - 1]) if token]
    return statuses or ["SITE"]


def re_split_filter(value: str) -> list[str]:
    # GATK versions have used both comma and semicolon inside an allele status.
    return [token.strip() for group in value.split(";") for token in group.split(",")]


def _record_filters(record: pysam.VariantRecord, allele_status: list[str]) -> list[str]:
    filters = [str(value) for value in record.filter if str(value) not in {"PASS", "."}]
    filters.extend(
        f"ALLELE:{value}" for value in allele_status if value.upper() not in {"PASS", "SITE"}
    )
    return sorted(set(filters))


def _record_digest(record: pysam.VariantRecord) -> str:
    return hashlib.sha256(str(record).rstrip("\n").encode()).hexdigest()


def parse_filtered_mutect_vcf(
    input_vcf: str | Path,
    provenance: MitochondrialProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    path = Path(input_vcf)
    raw_artifact_sha256 = sha256_file(path)
    source = pysam.VariantFile(str(path))
    rows: list[dict[str, Any]] = []
    source_records = 0
    try:
        if len(source.header.samples) != 1:
            raise PortableEvidenceError("mitochondrial VCF must contain exactly one sample")
        for record in source:
            source_records += 1
            if record.contig not in {"MT", "chrM"}:
                raise PortableEvidenceError("mitochondrial VCF contains a non-MT record")
            ref = record.ref.upper()
            if CONCRETE_ALLELE.fullmatch(ref) is None:
                raise PortableEvidenceError("mitochondrial REF is not concrete uppercase DNA")
            alts = tuple(str(alt).upper() for alt in record.alts or ())
            if not alts or any(CONCRETE_ALLELE.fullmatch(alt) is None for alt in alts):
                raise PortableEvidenceError("mitochondrial ALT is not concrete uppercase DNA")
            sample = next(iter(record.samples.values()))
            depth = _integer(sample.get("DP"), "FORMAT/DP")
            source_digest = _record_digest(record)
            for alt_index, alt in enumerate(alts, start=1):
                allele_fraction = _float(
                    _format_allele_value(sample, "AF", alt_index, number="A", required=True),
                    "FORMAT/AF",
                )
                ref_depth = _integer(
                    _format_allele_value(sample, "AD", 0, number="R", required=True),
                    "FORMAT/AD reference depth",
                )
                alt_depth = _integer(
                    _format_allele_value(sample, "AD", alt_index, number="R", required=True),
                    "FORMAT/AD alternate depth",
                )
                forward_alt, reverse_alt = _strand_counts(record, alt_index)
                allele_status = _allele_filter_status(record, alt_index)
                filter_reasons = _record_filters(record, allele_status)
                ocm = _info_value(record, "OCM")
                original_mismatch = int(ocm) if ocm is not None else 0
                possible_numt = (
                    any(reason.lower().endswith("possible_numt") for reason in filter_reasons)
                    or original_mismatch > 0
                )
                limitations = [
                    "short-read mitochondrial evidence is not a diagnosis",
                    "low allele fractions and NUMT context require orthogonal confirmation",
                ]
                if not provenance.circular_boundary_audited:
                    limitations.append("circular control-region boundary was not audited")
                if possible_numt:
                    limitations.append("possible NUMT or paralogous alignment evidence")
                confidence = (
                    "SUPPORTED"
                    if not filter_reasons
                    and provenance.circular_boundary_audited
                    and depth > 0
                    and allele_fraction > 0
                    else "REVIEW"
                )
                rpa = _tuple(_info_value(record, "RPA"))
                short_tandem_repeat = bool(_info_value(record, "STR") or False)
                repeat_unit = _info_value(record, "RU")
                event_id = f"mt:chrM:{record.pos}:{ref}:{alt}"
                row: dict[str, Any] = {
                    "chrom": "chrM",
                    "pos": int(record.pos),
                    "ref": ref,
                    "alt": alt,
                    "allele_fraction": allele_fraction,
                    "depth": depth,
                    "ref_depth": ref_depth,
                    "alt_depth": alt_depth,
                    "forward_alt": forward_alt,
                    "reverse_alt": reverse_alt,
                    "f1r2_alt": _optional_int(
                        _format_allele_value(sample, "F1R2", alt_index, number="R")
                    ),
                    "f2r1_alt": _optional_int(
                        _format_allele_value(sample, "F2R1", alt_index, number="R")
                    ),
                    "median_base_quality": _optional_float(_allele_value(record, "MBQ", alt_index)),
                    "median_fragment_length": _optional_float(
                        _allele_value(record, "MFRL", alt_index)
                    ),
                    "median_mapping_quality": _optional_float(
                        _allele_value(record, "MMQ", alt_index)
                    ),
                    "median_read_position": _optional_float(
                        _allele_value(record, "MPOS", alt_index)
                    ),
                    "tumor_lod": _optional_float(_allele_value(record, "TLOD", alt_index)),
                    "population_af_neg_log10": _optional_float(
                        _allele_value(record, "POPAF", alt_index)
                    ),
                    "unique_alt_read_count": _optional_int(
                        _allele_value(record, "AS_UNIQ_ALT_READ_COUNT", alt_index)
                    ),
                    "original_alignment_mismatch_reads": original_mismatch,
                    "allele_filter_status": allele_status,
                    "short_tandem_repeat": short_tandem_repeat,
                    "repeat_unit": str(repeat_unit) if repeat_unit is not None else None,
                    "repeat_count_ref": _optional_int(rpa[0]) if len(rpa) > 0 else None,
                    "repeat_count_alt": (
                        _optional_int(rpa[alt_index]) if len(rpa) > alt_index else None
                    ),
                    "numt_or_paralogous_context": possible_numt,
                    "circular_boundary_audited": provenance.circular_boundary_audited,
                    "evidence_id": "pending",
                    "event_id": event_id,
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
                    "source_record_digest": source_digest,
                    "filter_reasons": filter_reasons,
                    "confidence_tier": confidence,
                    "limitations": limitations,
                    "confirmation_required": True,
                    "concordance_group_id": None,
                    "concordance_status": "not_assessed",
                }
                row["evidence_id"] = make_evidence_id("mitochondrial_heteroplasmy", row)
                rows.append(row)
    finally:
        source.close()
    validated = validate_rows(
        "mitochondrial_heteroplasmy",
        rows,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "source_records": source_records,
        "allele_rows": len(validated),
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "numt_context_rows": sum(row["numt_or_paralogous_context"] for row in validated),
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": raw_artifact_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def _optional_int(value: Any) -> int | None:
    return None if value is None else _integer(value, "optional integer metric")


def _optional_float(value: Any) -> float | None:
    return None if value is None else _float(value, "optional floating metric")


def materialize_filtered_mutect_vcf(
    input_vcf: str | Path,
    output_parquet: str | Path,
    provenance: MitochondrialProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_filtered_mutect_vcf(
        input_vcf, provenance, forbidden_bindings=forbidden_bindings
    )
    output = write_extension(
        "mitochondrial_heteroplasmy",
        rows,
        output_parquet,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    return {"status": "PASS", "parsing": parsing, "output": output}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--technical-source-role", required=True)
    parser.add_argument("--caller-run-identity-sha256", required=True)
    parser.add_argument("--resource-manifest-sha256", required=True)
    parser.add_argument("--source-alignment-sha256", required=True)
    parser.add_argument("--raw-artifact-id", required=True)
    parser.add_argument("--caller-version", default="4.6.2.0")
    parser.add_argument("--circular-boundary-audited", action="store_true")
    args = parser.parse_args()
    provenance = MitochondrialProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        caller_version=args.caller_version,
        circular_boundary_audited=args.circular_boundary_audited,
    )
    result = materialize_filtered_mutect_vcf(args.input_vcf, args.output_parquet, provenance)
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
