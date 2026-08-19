"""Reference-backed allele normalization with explicit sort and validation."""

from __future__ import annotations

import os
import re
import subprocess
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pysam

from .config import BuildConfig
from .constants import STANDARD_CHROMS
from .errors import InputValidationError
from .util import atomic_output_directory, atomic_write_json, sha256_file
from .vcf_ingest import canonical_core_chrom

CONCRETE_ALLELE = re.compile(r"^[ACGT]+$")
CHROM_RANK = {chrom: rank for rank, chrom in enumerate(STANDARD_CHROMS)}


def _run(command: list[str], environment: dict[str, str]) -> dict[str, Any]:
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    result = {
        "executable": Path(command[0]).name,
        "arguments": [
            Path(value).name if value.startswith("/") else value for value in command[1:]
        ],
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise InputValidationError(
            f"{Path(command[0]).name} failed with exit {completed.returncode}: {detail[-4000:]}"
        )
    return result


def _complete_gt(value: Any, field: str) -> tuple[tuple[int, ...], bool]:
    text = str(value)
    if not text or text == "." or ("/" in text and "|" in text):
        raise InputValidationError(f"invalid {field} genotype provenance")
    separator = "|" if "|" in text else "/" if "/" in text else None
    tokens = text.split(separator) if separator is not None else [text]
    if len(tokens) not in {1, 2} or any(not token.isdigit() for token in tokens):
        raise InputValidationError(f"invalid {field} genotype provenance")
    genotype = tuple(int(token) for token in tokens)
    return genotype, separator == "|"


def _projected_genotype(record: pysam.VariantRecord) -> tuple[tuple[int, ...], bool]:
    source_gt, source_phased = _complete_gt(record.info.get("GAB_ORIG_GT"), "source")
    projected_gt, projected_phased = _complete_gt(record.info.get("GAB_PROJECTED_GT"), "projected")
    try:
        source_alt_index = int(record.info["GAB_ALT_INDEX"])
    except (KeyError, TypeError, ValueError) as exc:
        raise InputValidationError("invalid selected ALT provenance during normalization") from exc
    expected = tuple(1 if allele == source_alt_index else 0 for allele in source_gt)
    if source_alt_index < 1 or 1 not in expected:
        raise InputValidationError("selected ALT is absent from source genotype provenance")
    if projected_gt != expected or projected_phased != source_phased:
        raise InputValidationError("projected GT disagrees with source GT/ALT provenance")
    return projected_gt, projected_phased


def _concrete_atom(record: pysam.VariantRecord) -> tuple[str, bool]:
    alts = tuple(str(alt).upper() for alt in record.alts or ())
    if len(alts) == 1 and CONCRETE_ALLELE.fullmatch(alts[0]):
        return alts[0], False
    if len(alts) == 2 and CONCRETE_ALLELE.fullmatch(alts[0]) and alts[1] == "*":
        return alts[0], True
    raise InputValidationError(
        "bcftools atomization produced an unsupported ALT representation at "
        f"{record.contig}:{record.pos}: {','.join(alts) or '.'}"
    )


def _canonicalized_header(
    source: pysam.VariantHeader,
) -> tuple[pysam.VariantHeader, list[str], list[str]]:
    """Return a header whose core contigs use the public canonical order.

    Coordinate sorting in BCF/VCF follows header RID order, not numeric-looking
    contig names.  Reference FASTAs legitimately use other orders (for example
    ``1,10,...,19,2``), while the streaming clinical resources use
    ``1,2,...,22,X,Y,MT``.  Rebuild the header before the final external sort so
    every downstream merge join has one explicit, reference-independent order.
    Non-core contigs retain their original relative order after the core set.
    """

    by_canonical: dict[str, str] = {}
    non_core: list[str] = []
    for name in source.contigs:
        canonical = canonical_core_chrom(name)
        if canonical is None:
            non_core.append(name)
            continue
        previous = by_canonical.get(canonical)
        if previous is not None:
            raise InputValidationError(
                "multiple contig aliases map to the same canonical chromosome: "
                f"{previous}, {name} -> {canonical}"
            )
        by_canonical[canonical] = name

    canonical_core = [chrom for chrom in STANDARD_CHROMS if chrom in by_canonical]
    ordered_names = [by_canonical[chrom] for chrom in canonical_core] + non_core
    contig_records = {
        str(record.get("ID")): record
        for record in source.records
        if record.key == "contig" and record.get("ID") is not None
    }

    header = pysam.VariantHeader()
    for record in source.records:
        if record.key in {"fileformat", "contig"}:
            continue
        header.add_record(record)
    for name in ordered_names:
        record = contig_records.get(name)
        if record is not None:
            header.add_record(record)
        else:
            length = source.contigs[name].length
            header.contigs.add(name, length=length if length is not None else None)
    for sample in source.samples:
        header.add_sample(sample)
    return header, canonical_core, ordered_names


def _validate_canonical_order(path: Path) -> dict[str, Any]:
    source = pysam.VariantFile(str(path))
    try:
        canonical_header_order: list[str] = []
        seen_header_core: dict[str, str] = {}
        non_core_header_contigs = 0
        for name in source.header.contigs:
            canonical = canonical_core_chrom(name)
            if canonical is None:
                non_core_header_contigs += 1
                continue
            previous = seen_header_core.get(canonical)
            if previous is not None:
                raise InputValidationError(
                    "normalized header contains duplicate canonical aliases: "
                    f"{previous}, {name} -> {canonical}"
                )
            seen_header_core[canonical] = name
            canonical_header_order.append(canonical)
        expected_header_order = [chrom for chrom in STANDARD_CHROMS if chrom in seen_header_core]
        if canonical_header_order != expected_header_order:
            raise InputValidationError("normalized header does not use canonical core contig order")

        previous_locus: tuple[int, int] | None = None
        current_locus: tuple[int, int] | None = None
        seen_alleles: set[tuple[str, str]] = set()
        records = 0
        first_locus: str | None = None
        last_locus: str | None = None
        for record in source:
            canonical = canonical_core_chrom(record.contig)
            if canonical is None:
                raise InputValidationError(
                    "non-core record reached normalized concrete variant stream"
                )
            if not record.alts or len(record.alts) != 1:
                raise InputValidationError(
                    "non-biallelic record reached canonical order validation"
                )
            ref = record.ref.upper()
            alt = record.alts[0].upper()
            if not CONCRETE_ALLELE.fullmatch(ref) or not CONCRETE_ALLELE.fullmatch(alt):
                raise InputValidationError("non-concrete allele reached canonical order validation")
            locus = CHROM_RANK[canonical], record.pos
            if previous_locus is not None and locus < previous_locus:
                raise InputValidationError(
                    "normalized record stream does not use canonical coordinate order"
                )
            if locus != current_locus:
                current_locus = locus
                seen_alleles.clear()
            allele = ref, alt
            if allele in seen_alleles:
                raise InputValidationError(
                    "duplicate normalized query allele would violate streaming merge joins"
                )
            seen_alleles.add(allele)
            previous_locus = locus
            rendered = f"{canonical}:{record.pos}"
            if first_locus is None:
                first_locus = rendered
            last_locus = rendered
            records += 1
    finally:
        source.close()
    return {
        "status": "PASS",
        "canonical_core_header_order": canonical_header_order,
        "header_contigs": len(seen_header_core) + non_core_header_contigs,
        "non_core_header_contigs": non_core_header_contigs,
        "records": records,
        "first_locus": first_locus,
        "last_locus": last_locus,
        "duplicate_query_alleles": 0,
    }


def _biallelicize_atomized_bcf(source_path: Path, destination: Path) -> dict[str, int]:
    """Remove only bcftools' spanning-star companion and restore projected GT.

    ``bcftools norm --atomize`` represents each component of a complex carried
    allele as ``CONCRETE,*`` and may encode the carried source allele as index
    2.  The portable format requires one concrete ALT per row.  The selected
    source genotype is therefore restored from the immutable projection tags;
    every relationship is checked before a row is emitted.
    """

    counters: Counter[str] = Counter()
    source = pysam.VariantFile(str(source_path))
    try:
        if list(source.header.samples) != ["SAMPLE"]:
            raise InputValidationError("atomized projection must contain exactly SAMPLE")
        canonical_header, _canonical_core, _ordered_names = _canonicalized_header(source.header)
        output = pysam.VariantFile(str(destination), "wb", header=canonical_header)
        try:
            for record in source:
                counters["atomized_records"] += 1
                concrete_alt, spanning_star = _concrete_atom(record)
                projected_gt, projected_phased = _projected_genotype(record)
                sample = record.samples["SAMPLE"]
                raw_gt_value = sample.get("GT")
                raw_gt = tuple(raw_gt_value) if raw_gt_value is not None else ()
                if len(raw_gt) != len(projected_gt) or any(
                    allele is None or int(allele) < 0 or int(allele) >= len(record.alleles)
                    for allele in raw_gt
                ):
                    raise InputValidationError(
                        f"invalid atomized GT at {record.contig}:{record.pos}"
                    )
                semantic_gt = tuple(1 if int(allele) > 0 else 0 for allele in raw_gt)
                phase_disagrees = len(projected_gt) == 2 and bool(sample.phased) != projected_phased
                if semantic_gt != projected_gt or phase_disagrees:
                    raise InputValidationError(
                        "atomized GT/phase disagrees with projected provenance at "
                        f"{record.contig}:{record.pos}"
                    )

                projected = output.new_record(
                    contig=record.contig,
                    start=record.start,
                    stop=record.stop,
                    id=record.id,
                    alleles=(record.ref.upper(), concrete_alt),
                    qual=record.qual,
                )
                for item in record.filter:
                    projected.filter.add(item)
                for key, value in record.info.items():
                    projected.info[key] = value

                projected_sample = projected.samples["SAMPLE"]
                projected_sample["GT"] = projected_gt
                # A haploid VCF GT has no separator, so htslib reports it as
                # phased regardless of an attempted phase-bit assignment.
                # Preserve that representation; phase equality is meaningful
                # and enforced only for diploid genotypes.
                projected_sample.phased = (
                    projected_phased if len(projected_gt) == 2 else bool(sample.phased)
                )
                for key in ("DP", "GQ"):
                    value = sample.get(key)
                    if value is not None:
                        projected_sample[key] = value
                ref_ad = record.info.get("GAB_REF_AD")
                alt_ad = record.info.get("GAB_ALT_AD")
                if ref_ad is not None and alt_ad is not None:
                    projected_sample["AD"] = (int(ref_ad), int(alt_ad))
                output.write(projected)
                counters["biallelic_records"] += 1
                if spanning_star:
                    counters["spanning_star_companion_records_removed"] += 1
        finally:
            output.close()
    finally:
        source.close()
    if counters["atomized_records"] != counters["biallelic_records"]:
        raise InputValidationError("atomized/biallelic record reconciliation failed")
    return dict(sorted(counters.items()))


def normalize_projected_vcf(
    config: BuildConfig,
    projected_vcf: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    projected_vcf = Path(projected_vcf)
    destination = Path(destination)
    started = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    with atomic_output_directory(destination) as staging:
        atomized = staging / "atomized.unsorted.bcf"
        intermediate = staging / "normalized.unsorted.bcf"
        normalized = staging / "normalized.vcf.gz"
        environment = dict(os.environ)
        environment["TMPDIR"] = str(staging)
        check_ref = "e" if config.check_ref == "error" else "w"
        norm_command = [
            config.bcftools,
            "norm",
            "--no-version",
            "--fasta-ref",
            str(config.reference_fasta),
            "--atomize",
            "--atom-overlaps",
            "*",
            "--check-ref",
            check_ref,
            "--output-type",
            "b",
            "--output",
            str(atomized),
            str(projected_vcf),
        ]
        sort_command = [
            config.bcftools,
            "sort",
            "--temp-dir",
            str(staging / "bcftools-sort"),
            "--output-type",
            "z",
            "--output",
            str(normalized),
            str(intermediate),
        ]
        norm_result = _run(norm_command, environment)
        projection = _biallelicize_atomized_bcf(atomized, intermediate)
        sort_result = _run(sort_command, environment)
        index_result = _run(
            [config.bcftools, "index", "--tbi", "--force", str(normalized)], environment
        )
        canonical_order = _validate_canonical_order(normalized)
        count_result = _run([config.bcftools, "index", "--nrecords", str(normalized)], environment)
        try:
            record_count = int(count_result["stdout"].strip())
        except ValueError as exc:
            raise InputValidationError("bcftools did not return a normalized record count") from exc
        report = {
            "status": "PASS",
            "started_at": started,
            "completed_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
            "normalization": (
                "bcftools norm -f GRCh38 --atomize --atom-overlaps '*' followed by "
                "provenance-checked spanning-star removal/GT restoration, canonical "
                "core-contig header reconstruction and bcftools sort"
            ),
            "record_count": record_count,
            "projection": projection,
            "canonical_order": canonical_order,
            "input_sha256": sha256_file(projected_vcf),
            "outputs": {
                "normalized.vcf.gz": sha256_file(normalized),
                "normalized.vcf.gz.tbi": sha256_file(Path(f"{normalized}.tbi")),
                "atomized.unsorted.bcf": sha256_file(atomized),
                "normalized.unsorted.bcf": sha256_file(intermediate),
            },
            "commands": [norm_result, sort_result, index_result],
        }
        atomic_write_json(staging / "normalization_summary.json", report)
    return report
