"""Reconcile CNVpytor 1.3.2 combined RD+BAF TSV and audited HDF5 calls."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Any

from .cnvpytor_combined_hdf5_probe import COMBINED_KEYS, validate_combined_row
from .portable_model import (
    SHA256,
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    sha256_file,
    validate_rows,
    write_extension,
)

REGION = re.compile(r"^([^:]+):([0-9]+)-([0-9]+)$")
GENOTYPE = re.compile(r"^CN([0-9]+)/CN([0-9]+)$")
RD_SUFFIX = re.compile(r"^(?:_AT)?(?:_GC)?(?:_mask)?$")
SNP_SUFFIX = re.compile(r"^(?:_mask)?(?:_id)?$")


@dataclass(frozen=True)
class CnvpytorCombinedPolicy:
    min_size: int = 50_000
    max_combined_p_value: float = 1e-4
    max_q0_fraction: float = 0.5
    max_reference_gap_fraction: float = 0.5
    max_segment_gap_fraction: float = 0.5
    min_strict_mask_fraction: float = 0.5
    min_heterozygous_sites: int = 5


DEFAULT_COMBINED_POLICY = CnvpytorCombinedPolicy()


@dataclass(frozen=True)
class CnvpytorCombinedProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    raw_artifact_sha256: str
    source_baf_vcf_sha256: str
    bin_size: int
    rd_flags_suffix: str
    snp_flags_suffix: str
    rd_mask_applied: bool
    baf_mask_applied: bool
    caller_name: str = "CNVpytor"
    caller_version: str = "1.3.2"


@dataclass(frozen=True)
class Model:
    total: int
    first: int
    second: int
    likelihood: float
    cell_fraction: float


@dataclass(frozen=True)
class CombinedCall:
    call_type: str
    chrom: str
    start: int
    end: int
    size: int
    normalized_depth: float
    combined_p_value: float
    likelihood_deletion: float
    likelihood_loh: float
    likelihood_duplication: float
    q0_fraction: float
    reference_gap_fraction: float
    segment_gap_fraction: float
    strict_mask_fraction: float
    bin_size: int
    bins: int
    baf_deviation: float
    reported_rd_p_value: float
    baf_p_value: float
    heterozygous_sites: int
    homozygous_sites: int
    models: tuple[Model, ...]
    source_line: str

    def key(self) -> tuple[str, str, int, int, int, int]:
        return self.call_type, self.chrom, self.start, self.end, self.size, self.bins


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PortableEvidenceError(f"duplicate CNVpytor combined audit key: {key}")
        result[key] = value
    return result


def _load_audit(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError("invalid CNVpytor combined HDF5 audit") from exc
    if not isinstance(payload, dict):
        raise PortableEvidenceError("CNVpytor combined HDF5 audit root is not an object")
    return payload


def _canonical_chrom(value: str) -> str:
    core = value[3:] if value.startswith("chr") else value
    if core in {str(index) for index in range(1, 23)}:
        return f"chr{core}"
    raise PortableEvidenceError("CNVpytor combined call must be autosomal")


def _finite(value: str | float, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid CNVpytor combined {label}") from exc
    if not math.isfinite(number):
        raise PortableEvidenceError(f"non-finite CNVpytor combined {label}")
    return number


def _nonnegative_float(value: str | float, label: str) -> float:
    number = _finite(value, label)
    if number < 0:
        raise PortableEvidenceError(f"negative CNVpytor combined {label}")
    return number


def _nonnegative_int(value: str | float, label: str) -> int:
    number = _finite(value, label)
    if not number.is_integer() or number < 0:
        raise PortableEvidenceError(f"invalid CNVpytor combined {label}")
    return int(number)


def _parse_model(
    total_value: str | float,
    genotype_value: str,
    likelihood_value: str | float,
    cell_fraction_value: str | float,
) -> Model:
    total = _nonnegative_int(total_value, "model total copy number")
    match = GENOTYPE.fullmatch(str(genotype_value))
    if match is None:
        raise PortableEvidenceError("invalid CNVpytor combined model genotype")
    first, second = int(match.group(1)), int(match.group(2))
    if first + second != total:
        raise PortableEvidenceError("CNVpytor combined model genotype does not sum to total")
    likelihood = _nonnegative_float(likelihood_value, "model likelihood")
    cell_fraction = _nonnegative_float(cell_fraction_value, "model cell fraction")
    if cell_fraction > 1:
        raise PortableEvidenceError("CNVpytor combined model cell fraction exceeds one")
    return Model(total, first, second, likelihood, cell_fraction)


def _parse_tsv(path: Path, provenance: CnvpytorCombinedProvenance) -> list[CombinedCall]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("unreadable CNVpytor combined calls TSV") from exc
    calls: list[CombinedCall] = []
    seen: set[str] = set()
    for ordinal, line in enumerate(text.splitlines(), start=1):
        if not line:
            raise PortableEvidenceError("blank line in CNVpytor combined calls TSV")
        if line in seen:
            raise PortableEvidenceError("duplicate CNVpytor combined TSV line")
        seen.add(line)
        fields = line.split("\t")
        if len(fields) != 27:
            raise PortableEvidenceError(
                f"CNVpytor combined TSV line {ordinal} does not have 27 columns"
            )
        if fields[0] not in {"deletion", "duplication", "cnnloh"}:
            raise PortableEvidenceError("unsupported CNVpytor combined call type")
        match = REGION.fullmatch(fields[1])
        if match is None:
            raise PortableEvidenceError("invalid CNVpytor combined interval")
        chrom = _canonical_chrom(match.group(1))
        start, end = int(match.group(2)), int(match.group(3))
        size = _nonnegative_int(fields[2], "size")
        if start < 1 or end < start or size != end - start + 1:
            raise PortableEvidenceError("CNVpytor combined interval and size disagree")
        bin_size = _nonnegative_int(fields[12], "bin size")
        bins = _nonnegative_int(fields[13], "bin count")
        if bin_size != provenance.bin_size or bins <= 0:
            raise PortableEvidenceError("CNVpytor combined bin configuration differs")
        if fields[4] != fields[15]:
            raise PortableEvidenceError(
                "CNVpytor 1.3.2 reported RD p-value is not the expected combined duplicate"
            )
        model_one = _parse_model(fields[19], fields[20], fields[21], fields[22])
        model_two = _parse_model(fields[23], fields[24], fields[25], fields[26])
        expected_type = (
            "deletion"
            if model_one.total < 2
            else "duplication"
            if model_one.total > 2
            else "cnnloh"
        )
        if fields[0] != expected_type:
            raise PortableEvidenceError("CNVpytor combined type differs from its top model")
        call = CombinedCall(
            call_type=fields[0],
            chrom=chrom,
            start=start,
            end=end,
            size=size,
            normalized_depth=_nonnegative_float(fields[3], "normalized depth"),
            combined_p_value=_nonnegative_float(fields[4], "combined p-value"),
            likelihood_deletion=_nonnegative_float(fields[5], "deletion likelihood"),
            likelihood_loh=_nonnegative_float(fields[6], "LOH likelihood"),
            likelihood_duplication=_nonnegative_float(fields[7], "duplication likelihood"),
            q0_fraction=_nonnegative_float(fields[8], "Q0 fraction"),
            reference_gap_fraction=_nonnegative_float(fields[9], "reference-gap fraction"),
            segment_gap_fraction=_nonnegative_float(fields[10], "segment-gap fraction"),
            strict_mask_fraction=_nonnegative_float(fields[11], "strict-mask fraction"),
            bin_size=bin_size,
            bins=bins,
            baf_deviation=_nonnegative_float(fields[14], "BAF deviation"),
            reported_rd_p_value=_nonnegative_float(fields[15], "reported RD p-value"),
            baf_p_value=_nonnegative_float(fields[16], "BAF p-value"),
            heterozygous_sites=_nonnegative_int(fields[17], "heterozygous count"),
            homozygous_sites=_nonnegative_int(fields[18], "homozygous count"),
            models=(model_one, model_two),
            source_line=line,
        )
        for value in (
            call.q0_fraction,
            call.reference_gap_fraction,
            call.segment_gap_fraction,
            call.strict_mask_fraction,
        ):
            if value > 1:
                raise PortableEvidenceError("CNVpytor combined fraction exceeds one")
        if call.baf_deviation > 0.5:
            raise PortableEvidenceError("CNVpytor combined BAF deviation exceeds one half")
        calls.append(call)
    return calls


def _models_from_hdf5(row: list[float]) -> tuple[Model, ...]:
    models: list[Model] = []
    for index in range(10):
        offset = 1 + COMBINED_KEYS + index * 5
        models.append(
            Model(
                int(row[offset]),
                int(row[offset + 1]),
                int(row[offset + 2]),
                row[offset + 3],
                row[offset + 4],
            )
        )
    return tuple(models)


def _call_from_hdf5(contig: str, raw_row: list[float]) -> tuple[CombinedCall, list[float]]:
    try:
        row = validate_combined_row(raw_row)
    except ValueError as exc:
        raise PortableEvidenceError("invalid audited CNVpytor combined HDF5 row") from exc
    event_type = int(row[1])
    call_type = "deletion" if event_type == -1 else "duplication" if event_type == 1 else "cnnloh"
    return (
        CombinedCall(
            call_type=call_type,
            chrom=_canonical_chrom(contig),
            start=int(row[2]),
            end=int(row[3]),
            size=int(row[4]),
            normalized_depth=row[5],
            combined_p_value=row[6],
            likelihood_deletion=row[7],
            likelihood_loh=row[8],
            likelihood_duplication=row[9],
            q0_fraction=row[10],
            reference_gap_fraction=row[18],
            segment_gap_fraction=row[19],
            strict_mask_fraction=row[20],
            bin_size=0,
            bins=int(row[11]),
            baf_deviation=row[12],
            reported_rd_p_value=row[13],
            baf_p_value=row[14],
            heterozygous_sites=int(row[16]),
            homozygous_sites=int(row[17]),
            models=_models_from_hdf5(row),
            source_line="",
        ),
        row,
    )


def _audit_calls(
    payload: dict[str, Any], provenance: CnvpytorCombinedProvenance
) -> dict[tuple[str, str, int, int, int, int], list[tuple[CombinedCall, list[float]]]]:
    if (
        payload.get("schema_version") != 1
        or payload.get("status") != "PASS"
        or payload.get("caller_model") != "combined_rd_baf"
        or payload.get("cnvpytor_version") != provenance.caller_version
        or payload.get("bin_size") != provenance.bin_size
        or payload.get("rd_flags_suffix") != provenance.rd_flags_suffix
        or payload.get("snp_flags_suffix") != provenance.snp_flags_suffix
        or payload.get("volatile_attribute_values_persisted") is not False
        or payload.get("private_paths_persisted") is not False
    ):
        raise PortableEvidenceError("CNVpytor combined HDF5 audit boundary is invalid")
    datasets = payload.get("call_datasets")
    if not isinstance(datasets, dict) or not datasets:
        raise PortableEvidenceError("CNVpytor combined HDF5 audit lacks call datasets")
    calls: dict[tuple[str, str, int, int, int, int], list[tuple[CombinedCall, list[float]]]] = {}
    for dataset in datasets.values():
        if not isinstance(dataset, dict):
            raise PortableEvidenceError("invalid CNVpytor combined call dataset audit")
        if (
            dataset.get("bin_size") != provenance.bin_size
            or dataset.get("rd_flags_suffix") != provenance.rd_flags_suffix
            or dataset.get("snp_flags_suffix") != provenance.snp_flags_suffix
            or not isinstance(dataset.get("rows"), list)
        ):
            raise PortableEvidenceError("CNVpytor combined call dataset run key differs")
        contig = str(dataset.get("contig"))
        for raw_row in dataset["rows"]:
            if not isinstance(raw_row, list):
                raise PortableEvidenceError("invalid CNVpytor combined HDF5 call row")
            call, row = _call_from_hdf5(contig, raw_row)
            call = CombinedCall(**{**call.__dict__, "bin_size": provenance.bin_size})
            calls.setdefault(call.key(), []).append((call, row))
    return calls


def _close(first: float, second: float) -> bool:
    # CNVpytor renders several HDF5 fractions with four decimal places.  Treat
    # exactly half a rendered unit as inclusive without letting binary-float
    # subtraction move that decimal boundary by a few ulps.
    return math.isclose(first, second, rel_tol=5e-6, abs_tol=0.0) or abs(
        Decimal(str(first)) - Decimal(str(second))
    ) <= Decimal("0.00005")


def _calls_match(tsv: CombinedCall, hdf5: CombinedCall) -> bool:
    numeric = (
        "normalized_depth",
        "combined_p_value",
        "likelihood_deletion",
        "likelihood_loh",
        "likelihood_duplication",
        "q0_fraction",
        "reference_gap_fraction",
        "segment_gap_fraction",
        "strict_mask_fraction",
        "baf_deviation",
        "reported_rd_p_value",
        "baf_p_value",
    )
    if any(not _close(getattr(tsv, field), getattr(hdf5, field)) for field in numeric):
        return False
    if (
        tsv.heterozygous_sites != hdf5.heterozygous_sites
        or tsv.homozygous_sites != hdf5.homozygous_sites
    ):
        return False
    for first, second in zip(tsv.models, hdf5.models[:2], strict=True):
        if (
            first.total != second.total
            or first.first != second.first
            or first.second != second.second
            or not _close(first.likelihood, second.likelihood)
            or not _close(first.cell_fraction, second.cell_fraction)
        ):
            return False
    return True


def parse_cnvpytor_combined_calls(
    calls_tsv: str | Path,
    pytor_file: str | Path,
    hdf5_audit_json: str | Path,
    source_baf_vcf: str | Path,
    provenance: CnvpytorCombinedProvenance,
    *,
    policy: CnvpytorCombinedPolicy = DEFAULT_COMBINED_POLICY,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if provenance.bin_size <= 0:
        raise PortableEvidenceError("CNVpytor combined bin size must be positive")
    if SHA256.fullmatch(provenance.raw_artifact_sha256) is None:
        raise PortableEvidenceError("CNVpytor caller artifact SHA-256 is invalid")
    if (
        RD_SUFFIX.fullmatch(provenance.rd_flags_suffix) is None
        or SNP_SUFFIX.fullmatch(provenance.snp_flags_suffix) is None
    ):
        raise PortableEvidenceError("CNVpytor combined signal suffix is invalid")
    if provenance.rd_mask_applied != provenance.rd_flags_suffix.endswith("_mask"):
        raise PortableEvidenceError("CNVpytor combined RD mask provenance differs")
    if not provenance.baf_mask_applied or not provenance.snp_flags_suffix.startswith("_mask"):
        raise PortableEvidenceError("CNVpytor combined BAF strict mask is mandatory")
    tsv_path, pytor_path, audit_path, baf_path = map(
        Path, (calls_tsv, pytor_file, hdf5_audit_json, source_baf_vcf)
    )
    digests = {
        "calls_tsv_sha256": sha256_file(tsv_path),
        "pytor_sha256": sha256_file(pytor_path),
        "hdf5_audit_sha256": sha256_file(audit_path),
        "source_baf_vcf_sha256": sha256_file(baf_path),
    }
    if digests["source_baf_vcf_sha256"] != provenance.source_baf_vcf_sha256:
        raise PortableEvidenceError("CNVpytor combined BAF VCF differs from provenance")
    audit = _load_audit(audit_path)
    if audit.get("pytor_sha256") != digests["pytor_sha256"]:
        raise PortableEvidenceError("CNVpytor combined audit does not bind the .pytor file")
    tsv_calls = _parse_tsv(tsv_path, provenance)
    hdf5_calls = _audit_calls(audit, provenance)
    if len(tsv_calls) != sum(len(values) for values in hdf5_calls.values()):
        raise PortableEvidenceError("CNVpytor combined TSV/HDF5 cardinality differs")
    artifact_set_sha256 = hashlib.sha256(
        json.dumps(digests, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    rows: list[dict[str, Any]] = []
    for call in tsv_calls:
        candidates = hdf5_calls.get(call.key(), [])
        match_index = next(
            (
                index
                for index, (candidate, _) in enumerate(candidates)
                if _calls_match(call, candidate)
            ),
            None,
        )
        if match_index is None:
            raise PortableEvidenceError("CNVpytor combined TSV call does not match its HDF5 row")
        hdf5_call, hdf5_row = candidates.pop(match_index)
        if not _close(call.combined_p_value, call.reported_rd_p_value) or not _close(
            hdf5_call.combined_p_value, hdf5_call.reported_rd_p_value
        ):
            raise PortableEvidenceError(
                "CNVpytor combined reported RD p-value is not a combined duplicate"
            )
        model_rows = []
        for rank, model in enumerate(hdf5_call.models, start=1):
            padding = rank > 2 and (
                model.total,
                model.first,
                model.second,
                model.likelihood,
                model.cell_fraction,
            ) == (0, 0, 0, 0.0, 0.0)
            model_rows.append(
                {
                    "rank": rank,
                    "total_copy_number": model.total,
                    "haplotype1_copy_number": model.first,
                    "haplotype2_copy_number": model.second,
                    "likelihood": model.likelihood,
                    "cell_fraction": model.cell_fraction,
                    "reported_in_tsv": rank <= 2,
                    "padding": padding,
                }
            )
        top = hdf5_call.models[0]
        second = hdf5_call.models[1]
        filter_reasons: list[str] = []
        if call.size < policy.min_size:
            filter_reasons.append(f"SIZE_LT_{policy.min_size}")
        if call.combined_p_value > policy.max_combined_p_value:
            filter_reasons.append(f"COMBINED_P_GT_{policy.max_combined_p_value:g}")
        if call.q0_fraction > policy.max_q0_fraction:
            filter_reasons.append(f"Q0_GT_{policy.max_q0_fraction:g}")
        if call.reference_gap_fraction > policy.max_reference_gap_fraction:
            filter_reasons.append(
                f"REFERENCE_GAP_FRACTION_GT_{policy.max_reference_gap_fraction:g}"
            )
        if call.segment_gap_fraction > policy.max_segment_gap_fraction:
            filter_reasons.append(f"SEGMENT_GAP_FRACTION_GT_{policy.max_segment_gap_fraction:g}")
        if call.strict_mask_fraction < policy.min_strict_mask_fraction:
            filter_reasons.append(f"STRICT_MASK_FRACTION_LT_{policy.min_strict_mask_fraction:g}")
        if call.heterozygous_sites < policy.min_heterozygous_sites:
            filter_reasons.append(f"HET_SITES_LT_{policy.min_heterozygous_sites}")
        if not provenance.rd_mask_applied:
            filter_reasons.append("RD_STRICT_MASK_NOT_APPLIED")
        if top.likelihood <= second.likelihood or math.isclose(
            top.likelihood, second.likelihood, rel_tol=1e-9, abs_tol=0.0
        ):
            filter_reasons.append("TOP_TWO_MODELS_NOT_SEPARATED")
        confidence = "SUPPORTED" if not filter_reasons else "REVIEW"
        source_record_digest = hashlib.sha256(
            json.dumps(
                {
                    "hdf5_row": [value.hex() for value in hdf5_row],
                    "tsv_line": call.source_line,
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
        ).hexdigest()
        event_id = (
            f"cnv:{call.chrom}:{call.start}:{call.end}:{call.call_type}:"
            f"combined-bin{provenance.bin_size}"
        )
        row: dict[str, Any] = {
            "call_type": call.call_type,
            "caller_model": "combined_rd_baf",
            "chrom": call.chrom,
            "start_pos": call.start,
            "end_pos": call.end,
            "copy_number": float(top.total),
            "copy_number_basis": "combined_rd_baf_top_model_total_copy_number",
            "read_depth_copy_number": call.normalized_depth * 2,
            "major_copy_number": float(max(top.first, top.second)),
            "minor_copy_number": float(min(top.first, top.second)),
            "normalized_read_depth": call.normalized_depth,
            "log2_ratio": math.log2(call.normalized_depth) if call.normalized_depth > 0 else None,
            "baf": None,
            "baf_deviation_from_half": call.baf_deviation,
            "combined_p_value": call.combined_p_value,
            "caller_reported_rd_p_value": call.reported_rd_p_value,
            "baf_p_value": call.baf_p_value,
            "caller_reported_rd_p_value_duplicates_combined": True,
            "likelihood_deletion": call.likelihood_deletion,
            "likelihood_copy_neutral_loh": call.likelihood_loh,
            "likelihood_duplication": call.likelihood_duplication,
            "e_value_1": None,
            "e_value_2": None,
            "e_value_3": None,
            "e_value_4": None,
            "q0_fraction": call.q0_fraction,
            "reference_gap_fraction": call.reference_gap_fraction,
            "segment_gap_fraction": call.segment_gap_fraction,
            "strict_mask_fraction": call.strict_mask_fraction,
            "distance_to_large_gap": None,
            "rd_mask_applied": provenance.rd_mask_applied,
            "baf_mask_applied": provenance.baf_mask_applied,
            "bin_size": provenance.bin_size,
            "bins_observed": call.bins,
            "heterozygous_sites": call.heterozygous_sites,
            "homozygous_sites": call.homozygous_sites,
            "copy_number_models": model_rows,
            "source_baf_vcf_sha256": provenance.source_baf_vcf_sha256,
            "sex_ploidy": 2,
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
            "raw_artifact_sha256": provenance.raw_artifact_sha256,
            "source_record_digest": source_record_digest,
            "filter_reasons": sorted(filter_reasons),
            "confidence_tier": confidence,
            "limitations": [
                "CNVpytor combined RD+BAF is a prototype caller and not a diagnosis",
                "combined segments do not provide exact sequence breakpoints",
                "unphased BAF does not establish long-range haplotype phase",
                "CNVpytor 1.3.2 reports the combined p-value again in its RD p-value field",
                "copy-number model likelihoods are retained as caller scores, not normalized probabilities",
            ],
            "confirmation_required": True,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
        row["evidence_id"] = make_evidence_id("copy_number", row)
        rows.append(row)
    if any(hdf5_calls.values()):
        raise PortableEvidenceError("CNVpytor combined HDF5 contains unmatched calls")
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
        "caller_model": "combined_rd_baf",
        "tsv_calls": len(tsv_calls),
        "hdf5_calls": len(tsv_calls),
        "typed_rows": len(validated),
        "record_reconciliation_delta": len(tsv_calls) - len(validated),
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "reported_rd_p_value_duplicates_combined": True,
        **digests,
        "raw_artifact_id": provenance.raw_artifact_id,
        "caller_raw_artifact_sha256": provenance.raw_artifact_sha256,
        "adapter_input_set_sha256": artifact_set_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "extension_required": bool(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def materialize_cnvpytor_combined_calls(
    calls_tsv: str | Path,
    pytor_file: str | Path,
    hdf5_audit_json: str | Path,
    source_baf_vcf: str | Path,
    output_parquet: str | Path,
    provenance: CnvpytorCombinedProvenance,
    *,
    policy: CnvpytorCombinedPolicy = DEFAULT_COMBINED_POLICY,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_cnvpytor_combined_calls(
        calls_tsv,
        pytor_file,
        hdf5_audit_json,
        source_baf_vcf,
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
    parser.add_argument("--source-baf-vcf", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--technical-source-role", required=True)
    parser.add_argument("--caller-run-identity-sha256", required=True)
    parser.add_argument("--resource-manifest-sha256", required=True)
    parser.add_argument("--source-alignment-sha256", required=True)
    parser.add_argument("--raw-artifact-id", required=True)
    parser.add_argument("--raw-artifact-sha256", required=True)
    parser.add_argument("--source-baf-vcf-sha256", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--rd-flags-suffix", required=True)
    parser.add_argument("--snp-flags-suffix", required=True)
    parser.add_argument("--rd-mask-applied", action="store_true")
    parser.add_argument("--baf-mask-applied", action="store_true")
    args = parser.parse_args()
    provenance = CnvpytorCombinedProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        raw_artifact_sha256=args.raw_artifact_sha256,
        source_baf_vcf_sha256=args.source_baf_vcf_sha256,
        bin_size=args.bin_size,
        rd_flags_suffix=args.rd_flags_suffix,
        snp_flags_suffix=args.snp_flags_suffix,
        rd_mask_applied=args.rd_mask_applied,
        baf_mask_applied=args.baf_mask_applied,
    )
    result = materialize_cnvpytor_combined_calls(
        args.calls_tsv,
        args.pytor_file,
        args.hdf5_audit_json,
        args.source_baf_vcf,
        args.output_parquet,
        provenance,
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
