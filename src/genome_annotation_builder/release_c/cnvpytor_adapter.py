"""Reconcile CNVpytor mean-shift TSV calls against audited HDF5 rows."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .portable_model import (
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    sha256_file,
    validate_rows,
    write_extension,
)

REGION = re.compile(r"^([^:]+):([0-9]+)-([0-9]+)$")


@dataclass(frozen=True)
class CnvpytorPolicy:
    min_size: int = 50_000
    max_e_value_1: float = 1e-4
    max_q0_fraction: float = 0.5
    max_reference_gap_fraction: float = 0.5
    min_distance_to_large_gap: int = 100_000


DEFAULT_CNVPYTOR_POLICY = CnvpytorPolicy()


@dataclass(frozen=True)
class CnvpytorProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    bin_size: int
    flags_suffix: str
    rd_mask_applied: bool
    x_ploidy: int | None = None
    y_ploidy: int | None = None
    caller_name: str = "CNVpytor"
    caller_version: str = "1.3.2"


@dataclass(frozen=True)
class RdCall:
    call_type: str
    chrom: str
    start: int
    end: int
    size: int
    normalized_read_depth: float
    e_value_1: float
    e_value_2: float
    e_value_3: float
    e_value_4: float
    q0_fraction: float
    reference_gap_fraction: float
    distance_to_large_gap: int
    source_line: str

    def key(self) -> tuple[str, str, int, int, int]:
        return self.call_type, self.chrom, self.start, self.end, self.size


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PortableEvidenceError(f"duplicate CNVpytor audit JSON key: {key}")
        result[key] = value
    return result


def _load_audit(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError("invalid CNVpytor HDF5 audit manifest") from exc
    if not isinstance(payload, dict):
        raise PortableEvidenceError("CNVpytor HDF5 audit root must be an object")
    return payload


def _canonical_chrom(value: str) -> str:
    core = value[3:] if value.startswith("chr") else value
    if core in {str(index) for index in range(1, 23)} | {"X", "Y"}:
        return f"chr{core}"
    raise PortableEvidenceError("CNVpytor call uses a non-canonical nuclear chromosome")


def _finite_float(value: str | float, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid CNVpytor {label}") from exc
    if not math.isfinite(result):
        raise PortableEvidenceError(f"non-finite CNVpytor {label}")
    return result


def _nonnegative_float(value: str | float, label: str) -> float:
    result = _finite_float(value, label)
    if result < 0:
        raise PortableEvidenceError(f"negative CNVpytor {label}")
    return result


def _parse_tsv(path: Path) -> list[RdCall]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("unreadable CNVpytor calls TSV") from exc
    calls: list[RdCall] = []
    seen_lines: set[str] = set()
    for ordinal, line in enumerate(text.splitlines(), start=1):
        if not line:
            raise PortableEvidenceError("blank line in CNVpytor calls TSV")
        if line in seen_lines:
            raise PortableEvidenceError("duplicate line in CNVpytor calls TSV")
        seen_lines.add(line)
        fields = line.split("\t")
        if len(fields) != 11:
            raise PortableEvidenceError(
                f"CNVpytor mean-shift TSV line {ordinal} does not have 11 columns"
            )
        call_type = fields[0]
        if call_type not in {"deletion", "duplication"}:
            raise PortableEvidenceError("unsupported CNVpytor mean-shift call type")
        match = REGION.fullmatch(fields[1])
        if match is None:
            raise PortableEvidenceError("invalid CNVpytor call interval")
        chrom = _canonical_chrom(match.group(1))
        start, end = int(match.group(2)), int(match.group(3))
        size = int(fields[2])
        distance = int(fields[10])
        if start < 1 or end < start or size != end - start + 1:
            raise PortableEvidenceError("CNVpytor call size and interval disagree")
        if distance < 0:
            raise PortableEvidenceError("CNVpytor call lacks exact-reference gap distance")
        call = RdCall(
            call_type=call_type,
            chrom=chrom,
            start=start,
            end=end,
            size=size,
            normalized_read_depth=_nonnegative_float(fields[3], "normalized depth"),
            e_value_1=_nonnegative_float(fields[4], "e-value 1"),
            e_value_2=_nonnegative_float(fields[5], "e-value 2"),
            e_value_3=_nonnegative_float(fields[6], "e-value 3"),
            e_value_4=_nonnegative_float(fields[7], "e-value 4"),
            q0_fraction=_nonnegative_float(fields[8], "Q0 fraction"),
            reference_gap_fraction=_nonnegative_float(fields[9], "reference-gap fraction"),
            distance_to_large_gap=distance,
            source_line=line,
        )
        if call.q0_fraction > 1 or call.reference_gap_fraction > 1:
            raise PortableEvidenceError("CNVpytor fractions must be within [0,1]")
        calls.append(call)
    return calls


def _audit_calls(
    payload: dict[str, Any], provenance: CnvpytorProvenance
) -> dict[tuple[str, str, int, int, int], list[tuple[RdCall, list[float]]]]:
    if (
        payload.get("status") != "PASS"
        or payload.get("schema_version") != 1
        or payload.get("cnvpytor_version") != provenance.caller_version
        or payload.get("bin_size") != provenance.bin_size
        or payload.get("flags_suffix") != provenance.flags_suffix
        or payload.get("volatile_attribute_values_persisted") is not False
        or payload.get("private_paths_persisted") is not False
    ):
        raise PortableEvidenceError("CNVpytor HDF5 audit boundary is invalid")
    datasets = payload.get("call_datasets")
    if not isinstance(datasets, dict) or not datasets:
        raise PortableEvidenceError("CNVpytor HDF5 audit lacks call datasets")
    calls: dict[tuple[str, str, int, int, int], list[tuple[RdCall, list[float]]]] = {}
    for dataset in datasets.values():
        if not isinstance(dataset, dict):
            raise PortableEvidenceError("invalid CNVpytor HDF5 call dataset audit")
        contig = _canonical_chrom(str(dataset.get("contig")))
        if (
            dataset.get("bin_size") != provenance.bin_size
            or dataset.get("flags_suffix") != provenance.flags_suffix
        ):
            raise PortableEvidenceError("CNVpytor HDF5 call dataset run key is wrong")
        rows = dataset.get("rows")
        if not isinstance(rows, list):
            raise PortableEvidenceError("CNVpytor HDF5 audit call rows are missing")
        for raw_row in rows:
            if not isinstance(raw_row, list) or len(raw_row) != 13 or int(raw_row[0]) != 12:
                raise PortableEvidenceError("invalid CNVpytor HDF5 mean-shift row")
            row = [_finite_float(value, "HDF5 call value") for value in raw_row]
            type_code = int(row[1])
            if type_code not in {-1, 1}:
                raise PortableEvidenceError("invalid CNVpytor HDF5 event type")
            call = RdCall(
                call_type="deletion" if type_code == -1 else "duplication",
                chrom=contig,
                start=int(row[2]),
                end=int(row[3]),
                size=int(row[4]),
                normalized_read_depth=row[5],
                e_value_1=row[6],
                e_value_2=row[7],
                e_value_3=row[8],
                e_value_4=row[9],
                q0_fraction=row[10],
                reference_gap_fraction=row[11],
                distance_to_large_gap=int(row[12]),
                source_line="",
            )
            calls.setdefault(call.key(), []).append((call, row))
    return calls


def _metrics_match(tsv: RdCall, hdf5: RdCall) -> bool:
    return (
        all(
            math.isclose(first, second, rel_tol=5e-6, abs_tol=5e-5)
            for first, second in (
                (tsv.normalized_read_depth, hdf5.normalized_read_depth),
                (tsv.e_value_1, hdf5.e_value_1),
                (tsv.e_value_2, hdf5.e_value_2),
                (tsv.e_value_3, hdf5.e_value_3),
                (tsv.e_value_4, hdf5.e_value_4),
                (tsv.q0_fraction, hdf5.q0_fraction),
                (tsv.reference_gap_fraction, hdf5.reference_gap_fraction),
            )
        )
        and tsv.distance_to_large_gap == hdf5.distance_to_large_gap
    )


def _ploidy(chrom: str, provenance: CnvpytorProvenance) -> int:
    if chrom == "chrX":
        value = provenance.x_ploidy
    elif chrom == "chrY":
        value = provenance.y_ploidy
    else:
        value = 2
    if value not in {1, 2}:
        raise PortableEvidenceError("CNVpytor sex-chromosome ploidy is not declared")
    return value


def parse_cnvpytor_mean_shift_calls(
    calls_tsv: str | Path,
    pytor_file: str | Path,
    hdf5_audit_json: str | Path,
    provenance: CnvpytorProvenance,
    *,
    policy: CnvpytorPolicy = DEFAULT_CNVPYTOR_POLICY,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if provenance.bin_size <= 0:
        raise PortableEvidenceError("CNVpytor bin size must be positive")
    tsv_path, pytor_path, audit_path = map(Path, (calls_tsv, pytor_file, hdf5_audit_json))
    tsv_sha256 = sha256_file(tsv_path)
    pytor_sha256 = sha256_file(pytor_path)
    audit_sha256 = sha256_file(audit_path)
    audit = _load_audit(audit_path)
    if audit.get("pytor_sha256") != pytor_sha256:
        raise PortableEvidenceError("CNVpytor HDF5 audit does not bind the supplied .pytor file")
    tsv_calls = _parse_tsv(tsv_path)
    hdf5_calls = _audit_calls(audit, provenance)
    if len(tsv_calls) != sum(len(values) for values in hdf5_calls.values()):
        raise PortableEvidenceError("CNVpytor TSV/HDF5 call cardinality mismatch")
    artifact_set_sha256 = hashlib.sha256(
        json.dumps(
            {
                "calls_tsv_sha256": tsv_sha256,
                "hdf5_audit_sha256": audit_sha256,
                "pytor_sha256": pytor_sha256,
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    rows: list[dict[str, Any]] = []
    for call in tsv_calls:
        candidates = hdf5_calls.get(call.key(), [])
        match_index = next(
            (
                index
                for index, (candidate, _) in enumerate(candidates)
                if _metrics_match(call, candidate)
            ),
            None,
        )
        if match_index is None:
            raise PortableEvidenceError("CNVpytor TSV call does not match its HDF5 signal row")
        _, hdf5_row = candidates.pop(match_index)
        if call.size % provenance.bin_size:
            raise PortableEvidenceError("CNVpytor call size is not an integral bin count")
        ploidy = _ploidy(call.chrom, provenance)
        copy_number = call.normalized_read_depth * ploidy
        filter_reasons: list[str] = []
        if call.size < policy.min_size:
            filter_reasons.append(f"SIZE_LT_{policy.min_size}")
        if call.e_value_1 > policy.max_e_value_1:
            filter_reasons.append(f"E_VALUE_1_GT_{policy.max_e_value_1:g}")
        if call.q0_fraction > policy.max_q0_fraction:
            filter_reasons.append(f"Q0_GT_{policy.max_q0_fraction:g}")
        if call.reference_gap_fraction > policy.max_reference_gap_fraction:
            filter_reasons.append(
                f"REFERENCE_GAP_FRACTION_GT_{policy.max_reference_gap_fraction:g}"
            )
        if call.distance_to_large_gap < policy.min_distance_to_large_gap:
            filter_reasons.append(f"DISTANCE_TO_GAP_LT_{policy.min_distance_to_large_gap}")
        confidence = "SUPPORTED" if not filter_reasons else "REVIEW"
        limitations = [
            "read-depth copy-number evidence is not a diagnosis",
            "mean-shift calls do not provide exact sequence breakpoints",
            "copy number is derived from normalized depth and declared chromosome ploidy",
        ]
        event_id = (
            f"cnv:{call.chrom}:{call.start}:{call.end}:{call.call_type}:bin{provenance.bin_size}"
        )
        source_record_digest = hashlib.sha256(
            json.dumps(
                {
                    "hdf5_row": [float(value).hex() for value in hdf5_row],
                    "tsv_line": call.source_line,
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
        ).hexdigest()
        row: dict[str, Any] = {
            "call_type": call.call_type,
            "caller_model": "mean_shift_read_depth",
            "chrom": call.chrom,
            "start_pos": call.start,
            "end_pos": call.end,
            "copy_number": copy_number,
            "copy_number_basis": "normalized_read_depth_times_declared_ploidy",
            "read_depth_copy_number": copy_number,
            "major_copy_number": None,
            "minor_copy_number": None,
            "normalized_read_depth": call.normalized_read_depth,
            "log2_ratio": (
                math.log2(call.normalized_read_depth) if call.normalized_read_depth > 0 else None
            ),
            "baf": None,
            "baf_deviation_from_half": None,
            "combined_p_value": None,
            "caller_reported_rd_p_value": None,
            "baf_p_value": None,
            "caller_reported_rd_p_value_duplicates_combined": None,
            "likelihood_deletion": None,
            "likelihood_copy_neutral_loh": None,
            "likelihood_duplication": None,
            "e_value_1": call.e_value_1,
            "e_value_2": call.e_value_2,
            "e_value_3": call.e_value_3,
            "e_value_4": call.e_value_4,
            "q0_fraction": call.q0_fraction,
            "reference_gap_fraction": call.reference_gap_fraction,
            "segment_gap_fraction": None,
            "strict_mask_fraction": None,
            "distance_to_large_gap": call.distance_to_large_gap,
            "rd_mask_applied": provenance.rd_mask_applied,
            "baf_mask_applied": False,
            "bin_size": provenance.bin_size,
            "bins_observed": call.size // provenance.bin_size,
            "heterozygous_sites": None,
            "homozygous_sites": None,
            "copy_number_models": [],
            "source_baf_vcf_sha256": None,
            "sex_ploidy": ploidy,
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
            "raw_artifact_sha256": artifact_set_sha256,
            "source_record_digest": source_record_digest,
            "filter_reasons": sorted(filter_reasons),
            "confidence_tier": confidence,
            "limitations": limitations,
            "confirmation_required": True,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
        row["evidence_id"] = make_evidence_id("copy_number", row)
        rows.append(row)
    if any(hdf5_calls.values()):
        raise PortableEvidenceError("CNVpytor HDF5 contains unmatched mean-shift calls")
    validated = (
        validate_rows(
            "copy_number",
            rows,
            expected_reference_sha256=provenance.reference_sha256,
            forbidden_bindings=forbidden_bindings,
        )
        if rows
        else []
    )
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "tsv_calls": len(tsv_calls),
        "hdf5_calls": len(tsv_calls),
        "typed_rows": len(validated),
        "record_reconciliation_delta": len(tsv_calls) - len(validated),
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "calls_tsv_sha256": tsv_sha256,
        "pytor_sha256": pytor_sha256,
        "hdf5_audit_sha256": audit_sha256,
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_set_sha256": artifact_set_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "extension_required": bool(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def materialize_cnvpytor_mean_shift_calls(
    calls_tsv: str | Path,
    pytor_file: str | Path,
    hdf5_audit_json: str | Path,
    output_parquet: str | Path,
    provenance: CnvpytorProvenance,
    *,
    policy: CnvpytorPolicy = DEFAULT_CNVPYTOR_POLICY,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_cnvpytor_mean_shift_calls(
        calls_tsv,
        pytor_file,
        hdf5_audit_json,
        provenance,
        policy=policy,
        forbidden_bindings=forbidden_bindings,
    )
    if not rows:
        return {"status": "PASS", "parsing": parsing, "output": None}
    output = write_extension(
        "copy_number",
        rows,
        output_parquet,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    return {"status": "PASS", "parsing": parsing, "output": output}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calls-tsv", type=Path, required=True)
    parser.add_argument("--pytor-file", type=Path, required=True)
    parser.add_argument("--hdf5-audit-json", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--technical-source-role", required=True)
    parser.add_argument("--caller-run-identity-sha256", required=True)
    parser.add_argument("--resource-manifest-sha256", required=True)
    parser.add_argument("--source-alignment-sha256", required=True)
    parser.add_argument("--raw-artifact-id", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--flags-suffix", required=True)
    parser.add_argument("--rd-mask-applied", action="store_true")
    parser.add_argument("--x-ploidy", type=int)
    parser.add_argument("--y-ploidy", type=int)
    args = parser.parse_args()
    provenance = CnvpytorProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        bin_size=args.bin_size,
        flags_suffix=args.flags_suffix,
        rd_mask_applied=args.rd_mask_applied,
        x_ploidy=args.x_ploidy,
        y_ploidy=args.y_ploidy,
    )
    result = materialize_cnvpytor_mean_shift_calls(
        args.calls_tsv,
        args.pytor_file,
        args.hdf5_audit_json,
        args.output_parquet,
        provenance,
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
