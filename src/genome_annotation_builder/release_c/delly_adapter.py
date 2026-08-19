"""Single-sample DELLY 2.5.1 BCF/VCF to portable structural evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from .portable_model import (
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    sha256_file,
    validate_rows,
    write_extension,
)

BND_PARTNER = re.compile(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]")
SV_TYPES = frozenset({"DEL", "DUP", "INV", "INS", "BND", "CNV"})


@dataclass(frozen=True)
class DellyProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    caller_name: str = "DELLY"
    caller_version: str = "2.5.1"


def _canonical_chrom(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise PortableEvidenceError("invalid DELLY chromosome")
    core = value[3:] if value.startswith("chr") else value
    if core in {"M", "MT"}:
        return "chrM"
    if core in {str(index) for index in range(1, 23)} | {"X", "Y"}:
        return f"chr{core}"
    raise PortableEvidenceError("DELLY record uses a non-canonical chromosome")


def _info(record: pysam.VariantRecord, key: str) -> Any:
    if key not in record.header.info:
        return None
    return record.info.get(key)


def _optional_int(value: Any, label: str) -> int | None:
    if value is None:
        return None
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid DELLY {label}") from exc
    return result


def _optional_nonnegative_int(value: Any, label: str) -> int | None:
    result = _optional_int(value, label)
    if result is not None and result < 0:
        raise PortableEvidenceError(f"negative DELLY {label}")
    return result


def _optional_float(value: Any, label: str) -> float | None:
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid DELLY {label}") from exc
    if result != result or result in {float("inf"), float("-inf")}:
        raise PortableEvidenceError(f"non-finite DELLY {label}")
    return result


def _copy_number(value: Any) -> tuple[float | None, bool]:
    """Map DELLY's pinned ``RDCN=-1`` missing sentinel to portable null."""

    result = _optional_float(value, "FORMAT/RDCN")
    if result is None:
        return None, False
    if result == -1:
        return None, True
    if result < 0:
        raise PortableEvidenceError("invalid negative DELLY FORMAT/RDCN")
    return result, False


def _confidence_interval(value: Any, label: str, *, imprecise: bool) -> list[int]:
    if value is None:
        if imprecise:
            raise PortableEvidenceError(f"imprecise DELLY record lacks INFO/{label}")
        return [0, 0]
    values = tuple(value) if isinstance(value, tuple) else (value,)
    if len(values) != 2:
        raise PortableEvidenceError(f"DELLY INFO/{label} requires two bounds")
    result = [_optional_int(item, f"INFO/{label}") for item in values]
    if any(item is None for item in result):
        raise PortableEvidenceError(f"DELLY INFO/{label} contains a missing bound")
    assert result[0] is not None and result[1] is not None
    if result[0] > result[1]:
        raise PortableEvidenceError(f"DELLY INFO/{label} bounds are reversed")
    return [result[0], result[1]]


def _sample_filters(value: Any) -> list[str]:
    if value is None:
        return ["GENOTYPE:MISSING_FT"]
    values = value if isinstance(value, tuple) else (value,)
    result: list[str] = []
    for item in values:
        for token in re.split(r"[,;]", str(item)):
            token = token.strip()
            if token and token not in {"PASS", "."}:
                result.append(f"GENOTYPE:{token}")
    return sorted(set(result))


def _genotype(value: Any) -> tuple[list[int] | None, str]:
    if value is None:
        return None, "no_call"
    values = tuple(value) if isinstance(value, tuple) else (value,)
    if not values or any(item is None for item in values):
        return None, "no_call"
    try:
        genotype = [int(item) for item in values]
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError("invalid DELLY genotype") from exc
    if any(item not in {0, 1} for item in genotype):
        raise PortableEvidenceError("DELLY genotype is not biallelic")
    return genotype, "non_reference" if any(item == 1 for item in genotype) else "reference"


def _likelihoods(value: Any) -> list[int]:
    if value is None:
        return []
    values = tuple(value) if isinstance(value, tuple) else (value,)
    result: list[int] = []
    for item in values:
        if item is None:
            continue
        parsed = _optional_nonnegative_int(item, "FORMAT/PL")
        assert parsed is not None
        result.append(parsed)
    return result


def _bnd_partner(record: pysam.VariantRecord) -> tuple[str, int]:
    chr2 = _info(record, "CHR2")
    pos2 = _optional_int(_info(record, "POS2"), "INFO/POS2")
    if chr2 is None or pos2 is None or pos2 < 1:
        raise PortableEvidenceError("DELLY BND lacks CHR2/POS2")
    chrom2 = _canonical_chrom(str(chr2))
    alts = tuple(record.alts or ())
    if len(alts) != 1:
        raise PortableEvidenceError("DELLY BND must have exactly one ALT")
    partner = BND_PARTNER.search(str(alts[0]))
    if partner is None:
        raise PortableEvidenceError("DELLY BND ALT lacks partner coordinates")
    if (_canonical_chrom(partner.group(1)), int(partner.group(2))) != (chrom2, pos2):
        raise PortableEvidenceError("DELLY BND ALT and CHR2/POS2 disagree")
    return chrom2, pos2


def parse_delly_calls(
    input_vcf: str | Path,
    provenance: DellyProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    path = Path(input_vcf)
    raw_artifact_sha256 = sha256_file(path)
    try:
        source = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise PortableEvidenceError("invalid DELLY BCF/VCF artifact") from exc
    rows: list[dict[str, Any]] = []
    source_records = 0
    rdcn_missing_sentinel_rows = 0
    seen_record_ids: set[str] = set()
    try:
        if len(source.header.samples) != 1:
            raise PortableEvidenceError("DELLY artifact must contain exactly one sample")
        for record in source:
            source_records += 1
            record_id = record.id
            if not record_id or record_id == ".":
                raise PortableEvidenceError("DELLY record lacks a caller ID")
            if record_id in seen_record_ids:
                raise PortableEvidenceError("duplicate DELLY caller record ID")
            seen_record_ids.add(record_id)
            sv_type = _info(record, "SVTYPE")
            if sv_type not in SV_TYPES:
                raise PortableEvidenceError("unsupported DELLY SVTYPE")
            chrom1 = _canonical_chrom(record.contig)
            pos1 = int(record.pos)
            if pos1 < 1:
                raise PortableEvidenceError("invalid DELLY primary breakpoint")
            if sv_type == "BND":
                chrom2, pos2 = _bnd_partner(record)
            else:
                chrom2 = chrom1
                pos2 = int(record.stop)
                if pos2 < pos1:
                    raise PortableEvidenceError("invalid DELLY interval endpoint")
            precise = bool(_info(record, "PRECISE"))
            imprecise = bool(_info(record, "IMPRECISE"))
            if precise == imprecise:
                raise PortableEvidenceError("DELLY precision flags are ambiguous")
            sample = next(iter(record.samples.values()))
            genotype, call_status = _genotype(sample.get("GT"))
            copy_number, rdcn_missing_sentinel = _copy_number(sample.get("RDCN"))
            if rdcn_missing_sentinel:
                rdcn_missing_sentinel_rows += 1
            site_filters = [
                f"SITE:{item}" for item in record.filter if str(item) not in {"PASS", "."}
            ]
            filter_reasons = sorted(set(site_filters + _sample_filters(sample.get("FT"))))
            sex_chromosome_model_boundary = chrom1 in {"chrX", "chrY"} or chrom2 in {
                "chrX",
                "chrY",
            }
            if sex_chromosome_model_boundary:
                filter_reasons = sorted({*filter_reasons, "MODEL:DIPLOID_ONLY_SEX_CHROMOSOME"})
            gq = _optional_float(sample.get("GQ"), "FORMAT/GQ")
            paired_alt = _optional_nonnegative_int(sample.get("DV"), "FORMAT/DV")
            split_alt = _optional_nonnegative_int(sample.get("RV"), "FORMAT/RV")
            paired_total = _optional_nonnegative_int(_info(record, "PE"), "INFO/PE")
            split_total = _optional_nonnegative_int(_info(record, "SR"), "INFO/SR")
            evidence_support = (paired_alt or 0) + (split_alt or 0)
            confidence = (
                "SUPPORTED"
                if not filter_reasons
                and call_status in {"non_reference", "reference"}
                and (gq is None or gq >= 20)
                and (call_status == "reference" or evidence_support > 0)
                else "REVIEW"
            )
            limitations = [
                "short-read structural-variant evidence is not a diagnosis",
                "breakpoints, copy number, and complex loci may require orthogonal confirmation",
                "INFO/PE and INFO/SR are caller-site discovery metrics, not per-sample support denominators",
            ]
            if imprecise:
                limitations.append("caller reports imprecise breakpoint intervals")
            if sex_chromosome_model_boundary:
                limitations.append(
                    "DELLY 2.5.1 genotyping assumes diploidy; sex-chromosome GT requires a separate ploidy-aware review"
                )
            if rdcn_missing_sentinel:
                limitations.append(
                    "DELLY RDCN=-1 missing/invalid sentinel is represented as null copy number"
                )
            if call_status != "non_reference":
                limitations.append(
                    "this row is a caller candidate observation, not a positive non-reference call"
                )
            mate_id = _info(record, "MATEID")
            if sv_type == "BND":
                mate_link_status = (
                    "paired_mate_id" if mate_id is not None else "single_record_partner_coordinates"
                )
            else:
                mate_link_status = "not_applicable"
                mate_id = None
            orientation = _info(record, "CT")
            if not isinstance(orientation, str) or not orientation:
                raise PortableEvidenceError("DELLY record lacks INFO/CT orientation")
            event_id = f"sv:{sv_type}:{chrom1}:{pos1}:{chrom2}:{pos2}:{orientation}"
            row: dict[str, Any] = {
                "caller_record_id": record_id,
                "chrom1": chrom1,
                "pos1": pos1,
                "chrom2": chrom2,
                "pos2": pos2,
                "sv_type": sv_type,
                "orientation": orientation,
                "length": _optional_int(_info(record, "SVLEN"), "INFO/SVLEN"),
                "genotype": genotype,
                "call_status": call_status,
                "genotype_quality": gq,
                "genotype_likelihoods": _likelihoods(sample.get("PL")),
                "copy_number": copy_number,
                "quality": _optional_float(record.qual, "QUAL"),
                "imprecise": imprecise,
                "cipos": _confidence_interval(_info(record, "CIPOS"), "CIPOS", imprecise=imprecise),
                "ciend": _confidence_interval(_info(record, "CIEND"), "CIEND", imprecise=imprecise),
                "paired_end_support": paired_alt,
                "caller_site_paired_end_support": paired_total,
                "paired_end_reference_reads": _optional_nonnegative_int(
                    sample.get("DR"), "FORMAT/DR"
                ),
                "split_read_support": split_alt,
                "caller_site_split_read_support": split_total,
                "split_read_reference_reads": _optional_nonnegative_int(
                    sample.get("RR"), "FORMAT/RR"
                ),
                "read_depth_support": None,
                "regional_depth": _optional_nonnegative_int(sample.get("RC"), "FORMAT/RC"),
                "mapping_quality": _optional_float(_info(record, "MAPQ"), "INFO/MAPQ"),
                "split_read_mapping_quality": _optional_float(
                    _info(record, "SRMAPQ"), "INFO/SRMAPQ"
                ),
                "insertion_length": _optional_nonnegative_int(
                    _info(record, "INSLEN"), "INFO/INSLEN"
                ),
                "homology_length": _optional_nonnegative_int(
                    _info(record, "HOMLEN"), "INFO/HOMLEN"
                ),
                "bnd_mate_id": str(mate_id) if mate_id is not None else None,
                "mate_link_status": mate_link_status,
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
                "source_record_digest": hashlib.sha256(
                    str(record).rstrip("\n").encode()
                ).hexdigest(),
                "filter_reasons": filter_reasons,
                "confidence_tier": confidence,
                "limitations": limitations,
                "confirmation_required": True,
                "concordance_group_id": None,
                "concordance_status": "not_assessed",
            }
            row["evidence_id"] = make_evidence_id("structural_variants", row)
            rows.append(row)
    finally:
        source.close()
    validated = validate_rows(
        "structural_variants",
        rows,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    status_counts = {
        status: sum(row["call_status"] == status for row in validated)
        for status in ("non_reference", "reference", "no_call")
    }
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "source_records": source_records,
        "typed_rows": len(validated),
        "record_reconciliation_delta": source_records - len(validated),
        "rdcn_missing_sentinel_rows": rdcn_missing_sentinel_rows,
        "call_status_counts": status_counts,
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "bnd_rows": sum(row["sv_type"] == "BND" for row in validated),
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": raw_artifact_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def materialize_delly_calls(
    input_vcf: str | Path,
    output_parquet: str | Path,
    provenance: DellyProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_delly_calls(input_vcf, provenance, forbidden_bindings=forbidden_bindings)
    output = write_extension(
        "structural_variants",
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
    parser.add_argument("--caller-version", default="2.5.1")
    args = parser.parse_args()
    provenance = DellyProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        caller_version=args.caller_version,
    )
    result = materialize_delly_calls(args.input_vcf, args.output_parquet, provenance)
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
