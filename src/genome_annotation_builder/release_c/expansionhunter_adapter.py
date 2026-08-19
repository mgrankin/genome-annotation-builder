"""Reconcile ExpansionHunter JSON/VCF output into portable repeat evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from .portable_model import (
    CONCRETE_ALLELE,
    IUPAC_DNA_MOTIF,
    PortableEvidenceError,
    canonical_row_multiset,
    make_evidence_id,
    sha256_file,
    validate_rows,
    write_extension,
)

REGION = re.compile(r"^(?P<chrom>[^:]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$")
COUNT = re.compile(r"^(?P<count>[0-9]+)(?P<lower_bound>\+)?$")
FINITE_CI = re.compile(r"^(?P<lower>[0-9]+)-(?P<upper>[0-9]+)$")
OPEN_CI = re.compile(r"^(?P<lower>[0-9]+)\+$")
HISTOGRAM_PAIR = re.compile(r"\(([0-9]+),\s*([0-9]+)\)")


@dataclass(frozen=True)
class ExpansionHunterProvenance:
    reference_sha256: str
    technical_source_role: str
    caller_run_identity_sha256: str
    resource_manifest_sha256: str
    source_alignment_sha256: str
    raw_artifact_id: str
    caller_name: str = "ExpansionHunter"
    caller_version: str = "5.0.0"


@dataclass(frozen=True)
class RepeatAllele:
    repeat_count: int
    ci_lower: int
    ci_upper: int | None
    count_is_lower_bound: bool
    support_type: str | None = None
    spanning_reads: int | None = None
    flanking_reads: int | None = None
    in_repeat_reads: int | None = None

    def identity(self) -> tuple[int, int, int | None, bool]:
        return (
            self.repeat_count,
            self.ci_lower,
            self.ci_upper,
            self.count_is_lower_bound,
        )


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PortableEvidenceError(f"duplicate ExpansionHunter JSON key: {key}")
        result[key] = value
    return result


def _load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(
                handle,
                object_pairs_hook=_reject_duplicate_json_keys,
                parse_constant=lambda token: (_ for _ in ()).throw(
                    PortableEvidenceError(f"non-finite ExpansionHunter JSON number: {token}")
                ),
            )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError("invalid ExpansionHunter JSON artifact") from exc
    if not isinstance(value, dict):
        raise PortableEvidenceError("ExpansionHunter JSON root must be an object")
    return value


def _sha256_json(value: Any) -> str:
    rendered = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode()).hexdigest()


def _canonical_chrom(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise PortableEvidenceError("invalid ExpansionHunter chromosome")
    core = value[3:] if value.startswith("chr") else value
    if core in {"M", "MT"}:
        return "chrM"
    if core in {str(index) for index in range(1, 23)} | {"X", "Y"}:
        return f"chr{core}"
    raise PortableEvidenceError("ExpansionHunter repeat uses a non-canonical chromosome")


def _split(value: Any, label: str, *, expected: int) -> list[str]:
    if not isinstance(value, str):
        raise PortableEvidenceError(f"missing ExpansionHunter {label}")
    tokens = value.split("/")
    if len(tokens) != expected or any(not token for token in tokens):
        raise PortableEvidenceError(f"ExpansionHunter {label} ploidy is inconsistent")
    return tokens


def _parse_count(value: str, label: str) -> tuple[int, bool]:
    match = COUNT.fullmatch(value)
    if match is None:
        raise PortableEvidenceError(f"invalid ExpansionHunter {label}: {value}")
    return int(match.group("count")), match.group("lower_bound") is not None


def _parse_ci(value: str, label: str) -> tuple[int, int | None]:
    finite = FINITE_CI.fullmatch(value)
    if finite is not None:
        lower = int(finite.group("lower"))
        upper = int(finite.group("upper"))
        if upper < lower:
            raise PortableEvidenceError(f"invalid ExpansionHunter {label}: {value}")
        return lower, upper
    opened = OPEN_CI.fullmatch(value)
    if opened is not None:
        return int(opened.group("lower")), None
    raise PortableEvidenceError(f"invalid ExpansionHunter {label}: {value}")


def _optional_read_count(value: str, label: str) -> int | None:
    if value == ".":
        return None
    try:
        parsed = int(value)
    except ValueError as exc:
        raise PortableEvidenceError(f"invalid ExpansionHunter {label}: {value}") from exc
    if parsed < 0:
        raise PortableEvidenceError(f"negative ExpansionHunter {label}")
    return parsed


def _histogram_total(value: Any, label: str) -> int:
    if not isinstance(value, str):
        raise PortableEvidenceError(f"missing ExpansionHunter {label}")
    pairs = HISTOGRAM_PAIR.findall(value)
    residual = HISTOGRAM_PAIR.sub("", value).replace(",", "").strip()
    if value == "()":
        return 0
    if not pairs or residual:
        raise PortableEvidenceError(f"invalid ExpansionHunter {label}")
    return sum(int(reads) for _, reads in pairs)


def _info(record: pysam.VariantRecord, key: str) -> Any:
    if key not in record.header.info:
        return None
    return record.info.get(key)


def _read_vcf(path: Path) -> tuple[dict[str, dict[str, Any]], int]:
    try:
        source = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise PortableEvidenceError("invalid ExpansionHunter VCF artifact") from exc
    records: dict[str, dict[str, Any]] = {}
    try:
        if len(source.header.samples) != 1:
            raise PortableEvidenceError("ExpansionHunter VCF must contain exactly one sample")
        for record in source:
            varid = _info(record, "VARID")
            if not isinstance(varid, str) or not varid:
                raise PortableEvidenceError("ExpansionHunter VCF record lacks INFO/VARID")
            if varid in records:
                raise PortableEvidenceError("duplicate ExpansionHunter VCF VARID")
            sample = next(iter(record.samples.values()))
            records[varid] = {
                "chrom": _canonical_chrom(record.contig),
                "pos": int(record.pos),
                "end": int(record.stop),
                "ref_repeat_count": _info(record, "REF"),
                "repeat_unit": _info(record, "RU"),
                "repeat_id": _info(record, "REPID"),
                "filters": sorted(
                    str(item) for item in record.filter if str(item) not in {"PASS", "."}
                ),
                "repcn": sample.get("REPCN"),
                "repci": sample.get("REPCI"),
                "support_type": sample.get("SO"),
                "spanning": sample.get("ADSP"),
                "flanking": sample.get("ADFL"),
                "in_repeat": sample.get("ADIR"),
                "coverage": sample.get("LC"),
                "rendered": str(record).rstrip("\n"),
            }
    except (KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, PortableEvidenceError):
            raise
        raise PortableEvidenceError("malformed ExpansionHunter VCF record") from exc
    finally:
        source.close()
    return records, len(records)


def _parse_region(value: Any) -> tuple[str, int, int]:
    if not isinstance(value, str):
        raise PortableEvidenceError("ExpansionHunter repeat lacks ReferenceRegion")
    match = REGION.fullmatch(value)
    if match is None:
        raise PortableEvidenceError("invalid ExpansionHunter ReferenceRegion")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start < 1 or end < start:
        raise PortableEvidenceError("invalid ExpansionHunter ReferenceRegion interval")
    return _canonical_chrom(match.group("chrom")), start, end


def _alleles_from_call(
    counts: Any,
    intervals: Any,
    *,
    ploidy: int,
    support_types: Any = None,
    spanning: Any = None,
    flanking: Any = None,
    in_repeat: Any = None,
) -> list[RepeatAllele]:
    count_tokens = _split(counts, "repeat genotype", expected=ploidy)
    interval_tokens = _split(intervals, "repeat confidence interval", expected=ploidy)
    support_tokens = (
        [None] * ploidy
        if support_types is None
        else _split(support_types, "support type", expected=ploidy)
    )
    support_fields: list[list[str | None]] = []
    for value, label in (
        (spanning, "spanning-read count"),
        (flanking, "flanking-read count"),
        (in_repeat, "in-repeat-read count"),
    ):
        support_fields.append(
            [None] * ploidy if value is None else _split(value, label, expected=ploidy)
        )
    alleles: list[RepeatAllele] = []
    for index, (count_token, interval_token) in enumerate(
        zip(count_tokens, interval_tokens, strict=True)
    ):
        repeat_count, lower_bound = _parse_count(count_token, "repeat count")
        ci_lower, ci_upper = _parse_ci(interval_token, "repeat confidence interval")
        if repeat_count < ci_lower or (ci_upper is not None and repeat_count > ci_upper):
            raise PortableEvidenceError("ExpansionHunter count lies outside its interval")
        if ci_upper is None:
            lower_bound = True
        alleles.append(
            RepeatAllele(
                repeat_count=repeat_count,
                ci_lower=ci_lower,
                ci_upper=ci_upper,
                count_is_lower_bound=lower_bound,
                support_type=support_tokens[index],
                spanning_reads=(
                    None
                    if support_fields[0][index] is None
                    else _optional_read_count(str(support_fields[0][index]), "spanning-read count")
                ),
                flanking_reads=(
                    None
                    if support_fields[1][index] is None
                    else _optional_read_count(str(support_fields[1][index]), "flanking-read count")
                ),
                in_repeat_reads=(
                    None
                    if support_fields[2][index] is None
                    else _optional_read_count(str(support_fields[2][index]), "in-repeat-read count")
                ),
            )
        )
    return sorted(
        alleles,
        key=lambda allele: (
            allele.repeat_count,
            allele.ci_lower,
            allele.ci_upper if allele.ci_upper is not None else 2**31,
            allele.count_is_lower_bound,
            allele.support_type or "",
        ),
    )


def _row_allele(row: dict[str, Any], prefix: str, allele: RepeatAllele | None) -> None:
    row[f"{prefix}_repeat_count"] = None if allele is None else allele.repeat_count
    row[f"{prefix}_ci_lower"] = None if allele is None else allele.ci_lower
    row[f"{prefix}_ci_upper"] = None if allele is None else allele.ci_upper
    row[f"{prefix}_count_is_lower_bound"] = False if allele is None else allele.count_is_lower_bound
    row[f"{prefix}_support_type"] = None if allele is None else allele.support_type
    row[f"{prefix}_spanning_reads"] = None if allele is None else allele.spanning_reads
    row[f"{prefix}_flanking_reads"] = None if allele is None else allele.flanking_reads
    row[f"{prefix}_in_repeat_reads"] = None if allele is None else allele.in_repeat_reads


def parse_expansionhunter_calls(
    input_json: str | Path,
    input_vcf: str | Path,
    provenance: ExpansionHunterProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    json_path = Path(input_json)
    vcf_path = Path(input_vcf)
    json_sha256 = sha256_file(json_path)
    vcf_sha256 = sha256_file(vcf_path)
    artifact_pair_sha256 = _sha256_json(
        {"expansionhunter_json_sha256": json_sha256, "expansionhunter_vcf_sha256": vcf_sha256}
    )
    payload = _load_json(json_path)
    loci = payload.get("LocusResults")
    if not isinstance(loci, dict):
        raise PortableEvidenceError("ExpansionHunter JSON lacks LocusResults")
    vcf_records, vcf_record_count = _read_vcf(vcf_path)
    json_variants: dict[str, tuple[str, dict[str, Any], dict[str, Any]]] = {}
    for locus_key, locus in loci.items():
        if not isinstance(locus, dict):
            raise PortableEvidenceError("invalid ExpansionHunter locus result")
        locus_id = locus.get("LocusId", locus_key)
        if not isinstance(locus_id, str) or not locus_id:
            raise PortableEvidenceError("invalid ExpansionHunter locus ID")
        variants = locus.get("Variants")
        if not isinstance(variants, dict):
            raise PortableEvidenceError("ExpansionHunter locus lacks Variants")
        for variant_key, variant in variants.items():
            if not isinstance(variant, dict):
                raise PortableEvidenceError("invalid ExpansionHunter variant result")
            variant_id = variant.get("VariantId", variant_key)
            if not isinstance(variant_id, str) or not variant_id:
                raise PortableEvidenceError("invalid ExpansionHunter variant ID")
            if variant_id in json_variants:
                raise PortableEvidenceError("duplicate ExpansionHunter JSON variant ID")
            json_variants[variant_id] = (locus_id, locus, variant)
    if set(json_variants) != set(vcf_records):
        missing_in_vcf = sorted(set(json_variants) - set(vcf_records))[:5]
        missing_in_json = sorted(set(vcf_records) - set(json_variants))[:5]
        raise PortableEvidenceError(
            "ExpansionHunter JSON/VCF variant identity mismatch: "
            f"missing_in_vcf={missing_in_vcf}, missing_in_json={missing_in_json}"
        )

    rows: list[dict[str, Any]] = []
    omissions: list[tuple[str, str]] = []
    for variant_id in sorted(json_variants):
        locus_id, locus, variant = json_variants[variant_id]
        variant_type = variant.get("VariantType")
        if variant_type != "Repeat":
            omissions.append((variant_id, f"OUTSIDE_REPEAT_EXTENSION:{variant_type}"))
            continue
        if variant.get("Genotype") is None:
            omissions.append((variant_id, "NO_REPEAT_CALL"))
            continue
        try:
            ploidy = int(locus.get("AlleleCount"))
        except (TypeError, ValueError) as exc:
            raise PortableEvidenceError("invalid ExpansionHunter AlleleCount") from exc
        if ploidy not in {1, 2}:
            raise PortableEvidenceError("ExpansionHunter repeat ploidy must be one or two")
        vcf_record = vcf_records[variant_id]
        chrom, start, end = _parse_region(variant.get("ReferenceRegion"))
        if (chrom, start, end) != (
            vcf_record["chrom"],
            vcf_record["pos"],
            vcf_record["end"],
        ):
            raise PortableEvidenceError("ExpansionHunter JSON/VCF interval mismatch")
        motif = variant.get("RepeatUnit")
        if not isinstance(motif, str) or motif.upper() != str(vcf_record["repeat_unit"]).upper():
            raise PortableEvidenceError("ExpansionHunter JSON/VCF repeat motif mismatch")
        motif = motif.upper()
        if IUPAC_DNA_MOTIF.fullmatch(motif) is None:
            raise PortableEvidenceError(
                "ExpansionHunter repeat motif contains invalid IUPAC DNA symbols"
            )
        motif_is_iupac_ambiguous = CONCRETE_ALLELE.fullmatch(motif) is None
        json_alleles = _alleles_from_call(
            variant.get("Genotype"),
            variant.get("GenotypeConfidenceInterval"),
            ploidy=ploidy,
        )
        vcf_alleles = _alleles_from_call(
            vcf_record["repcn"],
            vcf_record["repci"],
            ploidy=ploidy,
            support_types=vcf_record["support_type"],
            spanning=vcf_record["spanning"],
            flanking=vcf_record["flanking"],
            in_repeat=vcf_record["in_repeat"],
        )
        if [allele.identity() for allele in json_alleles] != [
            allele.identity() for allele in vcf_alleles
        ]:
            raise PortableEvidenceError("ExpansionHunter JSON/VCF repeat genotype mismatch")
        coverage = locus.get("Coverage")
        try:
            coverage = float(coverage)
            vcf_coverage = float(vcf_record["coverage"])
        except (TypeError, ValueError) as exc:
            raise PortableEvidenceError("missing ExpansionHunter locus coverage") from exc
        if (
            not math.isfinite(coverage)
            or coverage < 0
            or not math.isclose(coverage, vcf_coverage, rel_tol=1e-6, abs_tol=1e-5)
        ):
            raise PortableEvidenceError("ExpansionHunter JSON/VCF coverage mismatch")
        filter_reasons = [f"VCF:{value}" for value in vcf_record["filters"]]
        if motif_is_iupac_ambiguous:
            filter_reasons.append("CATALOG:IUPAC_AMBIGUOUS_MOTIF")
        short_read_resolvable = all(
            allele.ci_upper is not None
            and allele.repeat_count == allele.ci_lower == allele.ci_upper
            and not allele.count_is_lower_bound
            and allele.support_type == "SPANNING"
            for allele in vcf_alleles
        )
        limitations = [
            "short-read repeat estimates are evidence, not a diagnosis",
            "repeat thresholds and pathogenicity require locus-specific review",
            "off-target evidence status is unavailable unless separately catalog-audited",
        ]
        if not short_read_resolvable:
            limitations.append("at least one repeat allele is not exactly short-read-resolved")
        if motif_is_iupac_ambiguous:
            limitations.append(
                "repeat motif uses IUPAC ambiguity codes; sequence composition is not uniquely specified"
            )
        confidence = (
            "SUPPORTED"
            if not filter_reasons and coverage >= 10 and short_read_resolvable
            else "REVIEW"
        )
        source_record_digest = _sha256_json(
            {"json_variant": variant, "vcf_record": vcf_record["rendered"]}
        )
        event_id = f"repeat:{chrom}:{start}:{end}:{variant_id}"
        row: dict[str, Any] = {
            "locus_id": locus_id,
            "variant_id": variant_id,
            "variant_subtype": str(variant.get("VariantSubtype", "Repeat")),
            "chrom": chrom,
            "start_pos": start,
            "end_pos": end,
            "motif": motif,
            "ploidy": ploidy,
            "reference_repeat_count": (
                int(vcf_record["ref_repeat_count"])
                if vcf_record["ref_repeat_count"] is not None
                else None
            ),
            "locus_coverage": coverage,
            "spanning_reads": _histogram_total(
                variant.get("CountsOfSpanningReads"), "CountsOfSpanningReads"
            ),
            "flanking_reads": _histogram_total(
                variant.get("CountsOfFlankingReads"), "CountsOfFlankingReads"
            ),
            "in_repeat_reads": _histogram_total(
                variant.get("CountsOfInrepeatReads"), "CountsOfInrepeatReads"
            ),
            "off_target_evidence": None,
            "short_read_resolvable": short_read_resolvable,
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
            "raw_artifact_sha256": artifact_pair_sha256,
            "source_record_digest": source_record_digest,
            "filter_reasons": filter_reasons,
            "confidence_tier": confidence,
            "limitations": limitations,
            "confirmation_required": True,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
        _row_allele(row, "allele1", vcf_alleles[0])
        _row_allele(row, "allele2", vcf_alleles[1] if ploidy == 2 else None)
        row["evidence_id"] = make_evidence_id("repeat_expansions", row)
        rows.append(row)

    validated = (
        validate_rows(
            "repeat_expansions",
            rows,
            expected_reference_sha256=provenance.reference_sha256,
            forbidden_bindings=forbidden_bindings,
        )
        if rows
        else []
    )
    omission_lines = [f"{variant_id}\t{reason}" for variant_id, reason in sorted(omissions)]
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "json_loci": len(loci),
        "json_variants": len(json_variants),
        "vcf_records": vcf_record_count,
        "repeat_rows": len(validated),
        "supported_rows": sum(row["confidence_tier"] == "SUPPORTED" for row in validated),
        "review_rows": sum(row["confidence_tier"] == "REVIEW" for row in validated),
        "motif_class_counts": {
            "concrete": sum(
                "CATALOG:IUPAC_AMBIGUOUS_MOTIF" not in row["filter_reasons"] for row in validated
            ),
            "iupac_ambiguous": sum(
                "CATALOG:IUPAC_AMBIGUOUS_MOTIF" in row["filter_reasons"] for row in validated
            ),
        },
        "omitted_records": len(omissions),
        "omission_reason_counts": {
            reason: sum(item_reason == reason for _, item_reason in omissions)
            for reason in sorted({reason for _, reason in omissions})
        },
        "omission_identity_sha256": hashlib.sha256("\n".join(omission_lines).encode()).hexdigest(),
        "expansionhunter_json_sha256": json_sha256,
        "expansionhunter_vcf_sha256": vcf_sha256,
        "raw_artifact_id": provenance.raw_artifact_id,
        "raw_artifact_pair_sha256": artifact_pair_sha256,
        "canonical_row_multiset": canonical_row_multiset(validated),
        "extension_required": bool(validated),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    return validated, summary


def materialize_expansionhunter_calls(
    input_json: str | Path,
    input_vcf: str | Path,
    output_parquet: str | Path,
    provenance: ExpansionHunterProvenance,
    *,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    rows, parsing = parse_expansionhunter_calls(
        input_json,
        input_vcf,
        provenance,
        forbidden_bindings=forbidden_bindings,
    )
    if not rows:
        return {"status": "PASS", "parsing": parsing, "output": None}
    output = write_extension(
        "repeat_expansions",
        rows,
        output_parquet,
        expected_reference_sha256=provenance.reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    return {"status": "PASS", "parsing": parsing, "output": output}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--input-vcf", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--technical-source-role", required=True)
    parser.add_argument("--caller-run-identity-sha256", required=True)
    parser.add_argument("--resource-manifest-sha256", required=True)
    parser.add_argument("--source-alignment-sha256", required=True)
    parser.add_argument("--raw-artifact-id", required=True)
    parser.add_argument("--caller-version", default="5.0.0")
    args = parser.parse_args()
    provenance = ExpansionHunterProvenance(
        reference_sha256=args.reference_sha256,
        technical_source_role=args.technical_source_role,
        caller_run_identity_sha256=args.caller_run_identity_sha256,
        resource_manifest_sha256=args.resource_manifest_sha256,
        source_alignment_sha256=args.source_alignment_sha256,
        raw_artifact_id=args.raw_artifact_id,
        caller_version=args.caller_version,
    )
    result = materialize_expansionhunter_calls(
        args.input_json, args.input_vcf, args.output_parquet, provenance
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
