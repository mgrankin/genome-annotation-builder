"""Manta 1.6.0 germline SV adapter with reciprocal BND validation."""

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
BND_ENDS = {
    "++": ("3", "5"),
    "+-": ("3", "3"),
    "-+": ("5", "3"),
    "--": ("5", "5"),
}


@dataclass(frozen=True)
class MantaProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    caller_name: str = "Manta"
    caller_version: str = "1.6.0"


def _canonical_chrom(value: Any) -> str:
    core = _contig_key(value)
    if core == "MT":
        return "chrM"
    if core in {str(index) for index in range(1, 23)} | {"X", "Y"}:
        return f"chr{core}"
    raise PortableEvidenceError("Manta record uses a non-canonical chromosome")


def _contig_key(value: Any) -> str:
    """Return a comparison key without restricting the contig to the portable core."""

    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
        or any(character.isspace() or character in "[]:" for character in value)
    ):
        raise PortableEvidenceError("invalid Manta chromosome")
    core = value[3:] if value.startswith("chr") else value
    if not core:
        raise PortableEvidenceError("invalid Manta chromosome")
    return "MT" if core in {"M", "MT"} else core


def _portable_chrom(value: Any) -> str | None:
    core = _contig_key(value)
    if core == "MT":
        return "chrM"
    if core in {str(index) for index in range(1, 23)} | {"X", "Y"}:
        return f"chr{core}"
    return None


def _info(record: pysam.VariantRecord, key: str) -> Any:
    if key not in record.header.info:
        return None
    return record.info.get(key)


def _first(value: Any) -> Any:
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def _optional_int(value: Any, label: str, *, nonnegative: bool = False) -> int | None:
    value = _first(value)
    if value is None:
        return None
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid Manta {label}") from exc
    if nonnegative and result < 0:
        raise PortableEvidenceError(f"negative Manta {label}")
    return result


def _optional_float(value: Any, label: str) -> float | None:
    value = _first(value)
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid Manta {label}") from exc
    if result != result or result in {float("inf"), float("-inf")}:
        raise PortableEvidenceError(f"non-finite Manta {label}")
    return result


def _ci(record: pysam.VariantRecord, key: str) -> list[int]:
    value = _info(record, key)
    if value is None:
        return [0, 0]
    values = tuple(value) if isinstance(value, tuple) else (value,)
    if len(values) != 2:
        raise PortableEvidenceError(f"Manta INFO/{key} requires two bounds")
    result = [_optional_int(item, f"INFO/{key}") for item in values]
    if any(item is None for item in result):
        raise PortableEvidenceError(f"Manta INFO/{key} has a missing bound")
    assert result[0] is not None and result[1] is not None
    if result[0] > result[1]:
        raise PortableEvidenceError(f"Manta INFO/{key} bounds are reversed")
    return [result[0], result[1]]


def _mate_id(record: pysam.VariantRecord) -> str | None:
    value = _info(record, "MATEID")
    if value is None:
        return None
    values = tuple(value) if isinstance(value, tuple) else (value,)
    if len(values) != 1 or not values[0]:
        raise PortableEvidenceError("Manta BND requires exactly one MATEID")
    return str(values[0])


def _raw_bnd_partner(record: pysam.VariantRecord) -> tuple[str, int, str]:
    alts = tuple(record.alts or ())
    if len(alts) != 1:
        raise PortableEvidenceError("Manta BND requires exactly one ALT")
    alt = str(alts[0])
    match = BND_PARTNER.search(alt)
    if match is None:
        raise PortableEvidenceError("Manta BND ALT lacks partner coordinates")
    if alt.startswith("["):
        orientation = "--"
    elif alt.startswith("]"):
        orientation = "-+"
    elif "[" in alt:
        orientation = "++"
    elif "]" in alt:
        orientation = "+-"
    else:
        raise PortableEvidenceError("Manta BND ALT has no bracket orientation")
    partner_pos = int(match.group(2))
    if partner_pos < 1:
        raise PortableEvidenceError("Manta BND ALT has an invalid partner position")
    return _contig_key(match.group(1)), partner_pos, orientation


def _bnd_partner(record: pysam.VariantRecord) -> tuple[str, int, str]:
    partner_chrom, partner_pos, orientation = _raw_bnd_partner(record)
    return _canonical_chrom(partner_chrom), partner_pos, orientation


def _record_map(source: pysam.VariantFile) -> dict[str, pysam.VariantRecord]:
    records: dict[str, pysam.VariantRecord] = {}
    for record in source:
        if not record.id or record.id == ".":
            raise PortableEvidenceError("Manta record lacks a caller ID")
        if record.id in records:
            raise PortableEvidenceError("duplicate Manta caller record ID")
        records[record.id] = record
    return records


def _validate_bnd_pairs(
    records: dict[str, pysam.VariantRecord],
) -> dict[str, tuple[pysam.VariantRecord, str]]:
    mates: dict[str, tuple[pysam.VariantRecord, str]] = {}
    for record_id, record in records.items():
        if _info(record, "SVTYPE") != "BND":
            continue
        mate_id = _mate_id(record)
        if mate_id is None or mate_id not in records:
            raise PortableEvidenceError("Manta BND mate is missing")
        if mate_id == record_id:
            raise PortableEvidenceError("Manta BND cannot be its own mate")
        mate = records[mate_id]
        if _info(mate, "SVTYPE") != "BND" or _mate_id(mate) != record_id:
            raise PortableEvidenceError("Manta BND mate link is not reciprocal")
        partner_chrom, partner_pos, orientation = _raw_bnd_partner(record)
        mate_chrom = _contig_key(mate.contig)
        mate_ci = _ci(mate, "CIPOS")
        if partner_chrom != mate_chrom or not (
            int(mate.pos) + mate_ci[0] <= partner_pos <= int(mate.pos) + mate_ci[1]
        ):
            raise PortableEvidenceError("Manta BND ALT falls outside its mate interval")
        mate_orientation = _raw_bnd_partner(mate)[2]
        local_end, remote_end = BND_ENDS[orientation]
        mate_local_end, mate_remote_end = BND_ENDS[mate_orientation]
        if (local_end, remote_end) != (mate_remote_end, mate_local_end):
            raise PortableEvidenceError("Manta BND mate orientation is not reciprocal")
        mates[record_id] = (mate, orientation)
    return mates


def validate_manta_bnd_pairing(input_vcf: str | Path) -> dict[str, Any]:
    path = Path(input_vcf)
    artifact_sha256 = sha256_file(path)
    try:
        source = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise PortableEvidenceError("invalid Manta VCF artifact") from exc
    try:
        records = _record_map(source)
        mates = _validate_bnd_pairs(records)
    finally:
        source.close()
    pair_ids = {
        tuple(sorted((record_id, _mate_id(records[record_id]) or ""))) for record_id in mates
    }
    return {
        "status": "PASS",
        "source_records": len(records),
        "bnd_records": len(mates),
        "reciprocal_bnd_pairs": len(pair_ids),
        "artifact_sha256": artifact_sha256,
        "sample_identifiers_persisted": False,
        "private_paths_persisted": False,
    }


def _sample_filters(value: Any) -> list[str]:
    if value is None:
        return ["GENOTYPE:MISSING_FT"]
    values = value if isinstance(value, tuple) else (value,)
    filters: list[str] = []
    for value in values:
        for token in re.split(r"[,;]", str(value)):
            token = token.strip()
            if token and token not in {"PASS", "."}:
                filters.append(f"GENOTYPE:{token}")
    return sorted(set(filters))


def _genotype(value: Any) -> tuple[list[int] | None, str]:
    if value is None:
        return None, "no_call"
    values = tuple(value) if isinstance(value, tuple) else (value,)
    if not values or any(item is None for item in values):
        return None, "no_call"
    genotype = [int(item) for item in values]
    if any(item not in {0, 1} for item in genotype):
        raise PortableEvidenceError("Manta genotype is not biallelic")
    return genotype, "non_reference" if any(item == 1 for item in genotype) else "reference"


def _likelihoods(value: Any) -> list[int]:
    if value is None:
        return []
    values = tuple(value) if isinstance(value, tuple) else (value,)
    result: list[int] = []
    for value in values:
        if value is None:
            continue
        parsed = _optional_int(value, "FORMAT/PL", nonnegative=True)
        assert parsed is not None
        result.append(parsed)
    return result


def _allelic_support(value: Any, label: str) -> tuple[int | None, int | None]:
    if value is None:
        return None, None
    values = tuple(value) if isinstance(value, tuple) else (value,)
    if len(values) < 2:
        raise PortableEvidenceError(f"Manta {label} lacks REF/ALT support")
    ref = _optional_int(values[0], label, nonnegative=True)
    alt = _optional_int(values[1], label, nonnegative=True)
    return ref, alt


def _non_bnd_orientation(sv_type: str) -> str:
    return {
        "DEL": "3to5",
        "DUP": "5to3",
        "INV": "inversion",
        "INS": "not_applicable",
        "CNV": "not_applicable",
    }[sv_type]


def _bnd_event_id(
    record: pysam.VariantRecord,
    mate: pysam.VariantRecord,
    orientation: str,
) -> str:
    mate_orientation = _bnd_partner(mate)[2]
    endpoints = sorted(
        (
            (_canonical_chrom(record.contig), int(record.pos), *_ci(record, "CIPOS")),
            (_canonical_chrom(mate.contig), int(mate.pos), *_ci(mate, "CIPOS")),
        )
    )
    endpoint_text = ":".join(str(value) for endpoint in endpoints for value in endpoint)
    orientation_text = ":".join(sorted((orientation, mate_orientation)))
    return f"sv:BND:{endpoint_text}:{orientation_text}"


def parse_manta_germline_calls(
    input_vcf: str | Path,
    provenance: MantaProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    path = Path(input_vcf)
    artifact_sha256 = sha256_file(path)
    try:
        source = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise PortableEvidenceError("invalid Manta germline VCF artifact") from exc
    rows: list[dict[str, Any]] = []
    try:
        if len(source.header.samples) != 1:
            raise PortableEvidenceError("Manta germline VCF must contain exactly one sample")
        records = _record_map(source)
        bnd_mates = _validate_bnd_pairs(records)
        omitted_bnd_ids: set[str] = set()
        omitted_non_bnd_ids: set[str] = set()
        core_to_noncanonical_bnd_pairs = 0
        noncanonical_to_noncanonical_bnd_pairs = 0
        source_bnd_pair_ids = {
            tuple(sorted((record_id, _mate_id(record) or "")))
            for record_id, record in records.items()
            if _info(record, "SVTYPE") == "BND"
        }
        for left_id, right_id in source_bnd_pair_ids:
            left_core = _portable_chrom(records[left_id].contig) is not None
            right_core = _portable_chrom(records[right_id].contig) is not None
            if left_core and right_core:
                continue
            omitted_bnd_ids.update((left_id, right_id))
            if left_core or right_core:
                core_to_noncanonical_bnd_pairs += 1
            else:
                noncanonical_to_noncanonical_bnd_pairs += 1
        for record_id, record in records.items():
            sv_type = _info(record, "SVTYPE")
            if sv_type not in SV_TYPES:
                raise PortableEvidenceError("unsupported Manta SVTYPE")
            if sv_type != "BND" and _portable_chrom(record.contig) is None:
                omitted_non_bnd_ids.add(record_id)
        omitted_ids = omitted_bnd_ids | omitted_non_bnd_ids
        for record_id, record in records.items():
            sv_type = _info(record, "SVTYPE")
            if record_id in omitted_ids:
                continue
            chrom1 = _canonical_chrom(record.contig)
            pos1 = int(record.pos)
            if sv_type == "BND":
                mate, orientation = bnd_mates[record_id]
                chrom2 = _canonical_chrom(mate.contig)
                pos2 = int(mate.pos)
                ciend = _ci(mate, "CIPOS")
                event_id = _bnd_event_id(record, mate, orientation)
                mate_id = mate.id
                mate_link_status = "paired_mate_id"
            else:
                chrom2 = chrom1
                pos2 = int(record.stop)
                if pos2 < pos1:
                    raise PortableEvidenceError("invalid Manta interval endpoint")
                orientation = _non_bnd_orientation(str(sv_type))
                ciend = _ci(record, "CIEND")
                event_id = f"sv:{sv_type}:{chrom1}:{pos1}:{chrom2}:{pos2}:{orientation}"
                mate_id = None
                mate_link_status = "not_applicable"
            sample = next(iter(record.samples.values()))
            genotype, call_status = _genotype(sample.get("GT"))
            paired_ref, paired_alt = _allelic_support(sample.get("PR"), "FORMAT/PR")
            split_ref, split_alt = _allelic_support(sample.get("SR"), "FORMAT/SR")
            site_filters = [
                f"SITE:{item}" for item in record.filter if str(item) not in {"PASS", "."}
            ]
            filter_reasons = sorted(set(site_filters + _sample_filters(sample.get("FT"))))
            gq = _optional_float(sample.get("GQ"), "FORMAT/GQ")
            evidence_support = (paired_alt or 0) + (split_alt or 0)
            confidence = (
                "SUPPORTED"
                if not filter_reasons
                and call_status in {"non_reference", "reference"}
                and (gq is None or gq >= 15)
                and (call_status == "reference" or evidence_support > 0)
                else "REVIEW"
            )
            limitations = [
                "short-read assembly and breakpoint evidence is not a diagnosis",
                "complex or repetitive SVs may require orthogonal confirmation",
            ]
            if bool(_info(record, "IMPRECISE")):
                limitations.append("caller reports an imprecise breakpoint")
            if call_status != "non_reference":
                limitations.append(
                    "this row is a caller candidate observation, not a positive non-reference call"
                )
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
                "copy_number": None,
                "quality": _optional_float(record.qual, "QUAL"),
                "imprecise": bool(_info(record, "IMPRECISE")),
                "cipos": _ci(record, "CIPOS"),
                "ciend": ciend,
                "paired_end_support": paired_alt,
                "caller_site_paired_end_support": _optional_int(
                    _info(record, "PAIR_COUNT"), "INFO/PAIR_COUNT", nonnegative=True
                ),
                "paired_end_reference_reads": paired_ref,
                "split_read_support": split_alt,
                "caller_site_split_read_support": None,
                "split_read_reference_reads": split_ref,
                "read_depth_support": None,
                "regional_depth": _optional_int(
                    _info(record, "BND_DEPTH"), "INFO/BND_DEPTH", nonnegative=True
                ),
                "mapping_quality": None,
                "split_read_mapping_quality": None,
                "insertion_length": _optional_int(
                    _info(record, "SVINSLEN"), "INFO/SVINSLEN", nonnegative=True
                ),
                "homology_length": _optional_int(
                    _info(record, "HOMLEN"), "INFO/HOMLEN", nonnegative=True
                ),
                "bnd_mate_id": mate_id,
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
                "raw_artifact_sha256": artifact_sha256,
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
    pair_ids = {
        tuple(sorted((row["caller_record_id"], row["bnd_mate_id"] or "")))
        for row in validated
        if row["sv_type"] == "BND"
    }
    omitted_records = len(omitted_ids)
    reconciliation_delta = len(records) - len(validated)
    if (
        len(validated) + omitted_records != len(records)
        or reconciliation_delta != omitted_records
        or len(omitted_bnd_ids)
        != 2 * (core_to_noncanonical_bnd_pairs + noncanonical_to_noncanonical_bnd_pairs)
    ):
        raise PortableEvidenceError("Manta portable-scope reconciliation failed")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "portable_contig_scope": "grch38_core_chr1_22_x_y_m",
        "source_records": len(records),
        "typed_rows": len(validated),
        "record_reconciliation_delta": reconciliation_delta,
        "omitted_noncanonical_records": omitted_records,
        "omitted_noncanonical_non_bnd_records": len(omitted_non_bnd_ids),
        "omitted_noncanonical_bnd_records": len(omitted_bnd_ids),
        "omitted_noncanonical_bnd_pairs": len(omitted_bnd_ids) // 2,
        "core_to_noncanonical_bnd_pairs": core_to_noncanonical_bnd_pairs,
        "noncanonical_to_noncanonical_bnd_pairs": noncanonical_to_noncanonical_bnd_pairs,
        "source_bnd_records": len(bnd_mates),
        "source_reciprocal_bnd_pairs": len(source_bnd_pair_ids),
        "bnd_records": sum(row["sv_type"] == "BND" for row in validated),
        "reciprocal_bnd_pairs": len(pair_ids),
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_sha256": artifact_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def materialize_manta_germline_calls(
    input_vcf: str | Path,
    output_parquet: str | Path,
    provenance: MantaProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_manta_germline_calls(
        input_vcf, provenance, forbidden_bindings=forbidden_bindings
    )
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
    parser.add_argument("--caller-version", default="1.6.0")
    args = parser.parse_args()
    provenance = MantaProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        caller_version=args.caller_version,
    )
    result = materialize_manta_germline_calls(args.input_vcf, args.output_parquet, provenance)
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
