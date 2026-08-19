"""Pinned, no-shell mitochondrial realignment and Mutect2 workflow.

This module re-expresses the scientifically relevant normal/8,000-base-shifted
branches of the GATK mitochondrial WDL with one immutable OCI runtime.  It is
intended to be invoked by :mod:`oci_runner`; host paths, sample identifiers, and
shell interpolation are deliberately absent from the public caller profile.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

MT_LENGTH = 16_569
SHIFT_BASES = 8_000
MT_MD5 = "c68f52674c9fb33aef52dcf399755519"
SHIFTED_MT_MD5 = "8ca1f1339b6185c6f63cf3752285cc37"
NORMAL_INTERVAL = "chrM:576-16024"
SHIFTED_INTERVAL = "chrM:8025-9144"
SOURCE_CONTIG = re.compile(r"^(?:MT|chrM)$")
EXPECTED_POLICY_KEYS = {
    "schema_version",
    "source_contig",
    "threads",
    "max_reads_per_alignment_start",
    "max_alt_allele_count",
    "min_allele_fraction",
    "contamination_estimate",
}
RESOURCE_FILES = {
    "normal_fasta": "mt.fa",
    "normal_fai": "mt.fa.fai",
    "normal_dict": "mt.dict",
    "normal_amb": "mt.fa.amb",
    "normal_ann": "mt.fa.ann",
    "normal_bwt": "mt.fa.bwt",
    "normal_pac": "mt.fa.pac",
    "normal_sa": "mt.fa.sa",
    "shifted_fasta": "mt.shifted.fa",
    "shifted_fai": "mt.shifted.fa.fai",
    "shifted_dict": "mt.shifted.dict",
    "shifted_amb": "mt.shifted.fa.amb",
    "shifted_ann": "mt.shifted.fa.ann",
    "shifted_bwt": "mt.shifted.fa.bwt",
    "shifted_pac": "mt.shifted.fa.pac",
    "shifted_sa": "mt.shifted.fa.sa",
    "shift_back_chain": "ShiftBack.chain",
    "blacklist": "blacklist_sites.hg38.chrM.bed",
}


class MitochondrialWorkflowError(RuntimeError):
    """A path-free workflow failure suitable for the outer OCI runner log."""


def _fail(code: str) -> None:
    raise MitochondrialWorkflowError(code)


def _regular_file(path: Path, code: str) -> None:
    if path.is_symlink() or not path.is_file():
        _fail(code)


def _integer(value: Any, code: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        _fail(code)
    if value < minimum or value > maximum:
        _fail(code)
    return value


def _fraction(value: Any, code: str, *, nullable: bool = False) -> float | None:
    if nullable and value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(code)
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0 or parsed > 1:
        _fail(code)
    return parsed


def load_policy(path: str | Path) -> dict[str, Any]:
    """Load one identifier-free per-alignment MT filtering policy."""

    policy_path = Path(path)
    _regular_file(policy_path, "policy_missing_or_symlinked")
    try:
        value = json.loads(policy_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise MitochondrialWorkflowError("policy_unreadable") from exc
    if not isinstance(value, dict) or set(value) != EXPECTED_POLICY_KEYS:
        _fail("policy_schema_invalid")
    source_contig = value.get("source_contig")
    if not isinstance(source_contig, str) or SOURCE_CONTIG.fullmatch(source_contig) is None:
        _fail("policy_source_contig_invalid")
    if value.get("schema_version") != 1:
        _fail("policy_schema_invalid")
    return {
        "schema_version": 1,
        "source_contig": source_contig,
        "threads": _integer(value.get("threads"), "policy_threads_invalid", minimum=1, maximum=16),
        "max_reads_per_alignment_start": _integer(
            value.get("max_reads_per_alignment_start"),
            "policy_read_cap_invalid",
            minimum=1,
            maximum=10_000,
        ),
        "max_alt_allele_count": _integer(
            value.get("max_alt_allele_count"),
            "policy_alt_cap_invalid",
            minimum=1,
            maximum=20,
        ),
        "min_allele_fraction": _fraction(value.get("min_allele_fraction"), "policy_min_af_invalid"),
        "contamination_estimate": _fraction(
            value.get("contamination_estimate"),
            "policy_contamination_invalid",
            nullable=True,
        ),
    }


def _read_fasta(path: Path) -> tuple[str, str]:
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeError) as exc:
        raise MitochondrialWorkflowError("mt_reference_unreadable") from exc
    if not lines or not lines[0].startswith(">"):
        _fail("mt_reference_invalid")
    name = lines[0][1:].split()[0]
    sequence = "".join(lines[1:]).upper()
    if name != "chrM" or len(sequence) != MT_LENGTH or set(sequence) - set("ACGTN"):
        _fail("mt_reference_invalid")
    return name, sequence


def validate_resource_root(root: str | Path) -> dict[str, Path]:
    """Verify required exact-reference files and circular rotation semantics."""

    resource_root = Path(root)
    if resource_root.is_symlink() or not resource_root.is_dir():
        _fail("mt_resource_root_invalid")
    paths = {role: resource_root / name for role, name in RESOURCE_FILES.items()}
    for path in paths.values():
        _regular_file(path, "mt_resource_member_missing_or_symlinked")
    _, normal = _read_fasta(paths["normal_fasta"])
    _, shifted = _read_fasta(paths["shifted_fasta"])
    if hashlib.md5(normal.encode(), usedforsecurity=False).hexdigest() != MT_MD5:
        _fail("mt_reference_md5_mismatch")
    if hashlib.md5(shifted.encode(), usedforsecurity=False).hexdigest() != SHIFTED_MT_MD5:
        _fail("shifted_mt_reference_md5_mismatch")
    if shifted != normal[SHIFT_BASES:] + normal[:SHIFT_BASES]:
        _fail("shifted_mt_rotation_invalid")
    return paths


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _vcf_record_stream(path: Path) -> tuple[str, int]:
    """Hash the exact ordered record stream while excluding volatile VCF headers."""

    try:
        records = [
            line
            for line in path.read_bytes().splitlines(keepends=True)
            if not line.startswith(b"#")
        ]
    except OSError as exc:
        raise MitochondrialWorkflowError("final_mt_vcf_unreadable") from exc
    if any(not line.endswith((b"\n", b"\r")) for line in records):
        _fail("final_mt_vcf_record_truncated")
    digest = hashlib.sha256()
    for line in records:
        digest.update(line)
    return digest.hexdigest(), len(records)


def _run(
    name: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    logs: Path,
    stdout_path: Path | None = None,
) -> None:
    """Execute one exact argv vector, retaining full private logs without a shell."""

    log_path = logs / f"{name}.log"
    with log_path.open("wb") as stderr_handle:
        if stdout_path is None:
            stdout_handle: Any = stderr_handle
            completed = subprocess.run(
                list(arguments),
                cwd=cwd,
                check=False,
                stdin=subprocess.DEVNULL,
                stdout=stdout_handle,
                stderr=stderr_handle,
            )
        else:
            with stdout_path.open("xb") as stdout_handle:
                completed = subprocess.run(
                    list(arguments),
                    cwd=cwd,
                    check=False,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                )
    if completed.returncode != 0:
        _fail(f"step_failed:{name}:{completed.returncode}")


def _gatk(java_memory: str, tool: str, *arguments: str) -> list[str]:
    return [
        "gatk",
        "--java-options",
        f"-Xmx{java_memory} -Djava.io.tmpdir=/work/tmp",
        tool,
        *arguments,
    ]


def _align(
    label: str,
    *,
    reference: Path,
    interleaved_fastq: Path,
    unmapped_bam: Path,
    calls: Path,
    logs: Path,
    threads: int,
) -> tuple[Path, Path]:
    aligned_sam = calls / f"{label}.aligned.sam"
    merged_bam = calls / f"{label}.merged.bam"
    marked_bam = calls / f"{label}.marked.bam"
    sorted_bam = calls / f"{label}.bam"
    metrics = calls / f"{label}.duplicate-metrics.txt"
    _run(
        f"{label}-bwa-mem",
        [
            "bwa",
            "mem",
            "-K",
            "100000000",
            "-p",
            "-v",
            "3",
            "-t",
            str(threads),
            "-Y",
            str(reference),
            str(interleaved_fastq),
        ],
        cwd=calls,
        logs=logs,
        stdout_path=aligned_sam,
    )
    _run(
        f"{label}-merge-alignment",
        _gatk(
            "12g",
            "MergeBamAlignment",
            "--VALIDATION_STRINGENCY",
            "SILENT",
            "--EXPECTED_ORIENTATIONS",
            "FR",
            "--ATTRIBUTES_TO_RETAIN",
            "X0",
            "--ATTRIBUTES_TO_REMOVE",
            "NM",
            "--ATTRIBUTES_TO_REMOVE",
            "MD",
            "--ALIGNED_BAM",
            str(aligned_sam),
            "--UNMAPPED_BAM",
            str(unmapped_bam),
            "--OUTPUT",
            str(merged_bam),
            "--REFERENCE_SEQUENCE",
            str(reference),
            "--PAIRED_RUN",
            "true",
            "--SORT_ORDER",
            "unsorted",
            "--ALIGNED_READS_ONLY",
            "false",
            "--CLIP_ADAPTERS",
            "false",
            "--ADD_MATE_CIGAR",
            "true",
            "--MAX_INSERTIONS_OR_DELETIONS",
            "-1",
            "--PRIMARY_ALIGNMENT_STRATEGY",
            "MostDistant",
            "--UNMAPPED_READ_STRATEGY",
            "COPY_TO_TAG",
            "--ALIGNER_PROPER_PAIR_FLAGS",
            "true",
            "--UNMAP_CONTAMINANT_READS",
            "true",
            "--ADD_PG_TAG_TO_READS",
            "false",
        ),
        cwd=calls,
        logs=logs,
    )
    _run(
        f"{label}-mark-duplicates",
        _gatk(
            "12g",
            "MarkDuplicates",
            "--INPUT",
            str(merged_bam),
            "--OUTPUT",
            str(marked_bam),
            "--METRICS_FILE",
            str(metrics),
            "--VALIDATION_STRINGENCY",
            "SILENT",
            "--OPTICAL_DUPLICATE_PIXEL_DISTANCE",
            "2500",
            "--ASSUME_SORT_ORDER",
            "queryname",
            "--CLEAR_DT",
            "false",
            "--ADD_PG_TAG_TO_READS",
            "false",
        ),
        cwd=calls,
        logs=logs,
    )
    _run(
        f"{label}-sort",
        _gatk(
            "12g",
            "SortSam",
            "--INPUT",
            str(marked_bam),
            "--OUTPUT",
            str(sorted_bam),
            "--SORT_ORDER",
            "coordinate",
            "--CREATE_INDEX",
            "true",
        ),
        cwd=calls,
        logs=logs,
    )
    index = sorted_bam.with_suffix(".bai")
    _regular_file(sorted_bam, "mt_alignment_missing")
    _regular_file(index, "mt_alignment_index_missing")
    return sorted_bam, index


def run_workflow(
    *,
    alignment: str | Path,
    alignment_index: str | Path,
    full_reference: str | Path,
    full_reference_fai: str | Path,
    full_reference_dict: str | Path,
    mt_resource_root: str | Path,
    policy_path: str | Path,
    workspace: str | Path,
) -> dict[str, Any]:
    """Run the exact MT pipeline and return a path-free scientific summary."""

    source = Path(alignment)
    source_index = Path(alignment_index)
    reference = Path(full_reference)
    reference_fai = Path(full_reference_fai)
    reference_dict = Path(full_reference_dict)
    for path in (source, source_index, reference, reference_fai, reference_dict):
        _regular_file(path, "primary_input_missing_or_symlinked")
    policy = load_policy(policy_path)
    resources = validate_resource_root(mt_resource_root)
    root = Path(workspace)
    if root.is_symlink() or not root.is_dir() or any(root.iterdir()):
        _fail("workspace_not_empty_or_invalid")
    home = root / "home"
    temporary = root / "tmp"
    calls = root / "mt"
    logs = calls / "logs"
    for directory in (home, temporary, calls, logs):
        directory.mkdir(parents=True, exist_ok=False)

    subset = calls / "source-contig.bam"
    reverted = calls / "source-contig.unmapped.bam"
    interleaved = calls / "source-contig.interleaved.fastq"
    _run(
        "subset-source-contig",
        _gatk(
            "8g",
            "PrintReads",
            "--reference",
            str(reference),
            "--intervals",
            policy["source_contig"],
            "--read-filter",
            "MateOnSameContigOrNoMappedMateReadFilter",
            "--read-filter",
            "MateUnmappedAndUnmappedReadFilter",
            "--input",
            str(source),
            "--read-index",
            str(source_index),
            "--output",
            str(subset),
        ),
        cwd=calls,
        logs=logs,
    )
    _run(
        "revert-source-contig",
        _gatk(
            "8g",
            "RevertSam",
            "--INPUT",
            str(subset),
            "--OUTPUT_BY_READGROUP",
            "false",
            "--OUTPUT",
            str(reverted),
            "--VALIDATION_STRINGENCY",
            "LENIENT",
            "--ATTRIBUTE_TO_CLEAR",
            "FT",
            "--ATTRIBUTE_TO_CLEAR",
            "CO",
            "--SORT_ORDER",
            "queryname",
            "--RESTORE_ORIGINAL_QUALITIES",
            "false",
        ),
        cwd=calls,
        logs=logs,
    )
    _run(
        "sam-to-interleaved-fastq",
        _gatk(
            "8g",
            "SamToFastq",
            "--INPUT",
            str(reverted),
            "--FASTQ",
            str(interleaved),
            "--INTERLEAVE",
            "true",
            "--NON_PF",
            "true",
        ),
        cwd=calls,
        logs=logs,
    )
    normal_bam, _ = _align(
        "normal",
        reference=resources["normal_fasta"],
        interleaved_fastq=interleaved,
        unmapped_bam=reverted,
        calls=calls,
        logs=logs,
        threads=policy["threads"],
    )
    shifted_bam, _ = _align(
        "shifted",
        reference=resources["shifted_fasta"],
        interleaved_fastq=interleaved,
        unmapped_bam=reverted,
        calls=calls,
        logs=logs,
        threads=policy["threads"],
    )

    common = [
        "--read-filter",
        "MateOnSameContigOrNoMappedMateReadFilter",
        "--read-filter",
        "MateUnmappedAndUnmappedReadFilter",
        "--annotation",
        "StrandBiasBySample",
        "--mitochondria-mode",
        "--max-reads-per-alignment-start",
        str(policy["max_reads_per_alignment_start"]),
        "--max-mnp-distance",
        "0",
    ]
    normal_vcf = calls / "normal.raw.vcf.gz"
    shifted_vcf = calls / "shifted.raw.vcf.gz"
    for label, input_bam, mt_reference, interval, output_vcf in (
        ("normal", normal_bam, resources["normal_fasta"], NORMAL_INTERVAL, normal_vcf),
        (
            "shifted",
            shifted_bam,
            resources["shifted_fasta"],
            SHIFTED_INTERVAL,
            shifted_vcf,
        ),
    ):
        _run(
            f"{label}-mutect2",
            _gatk(
                "16g",
                "Mutect2",
                "--reference",
                str(mt_reference),
                "--input",
                str(input_bam),
                "--output",
                str(output_vcf),
                "--intervals",
                interval,
                *common,
            ),
            cwd=calls,
            logs=logs,
        )

    shifted_back = calls / "shifted-back.vcf"
    rejected = calls / "shifted-rejected.vcf"
    merged = calls / "merged.raw.vcf"
    merged_stats = calls / "merged.raw.vcf.stats"
    _run(
        "liftover-shifted-calls",
        _gatk(
            "8g",
            "LiftoverVcf",
            "--INPUT",
            str(shifted_vcf),
            "--OUTPUT",
            str(shifted_back),
            "--CHAIN",
            str(resources["shift_back_chain"]),
            "--REJECT",
            str(rejected),
            "--REFERENCE_SEQUENCE",
            str(resources["normal_fasta"]),
            "--VALIDATION_STRINGENCY",
            "STRICT",
            "--CREATE_INDEX",
            "true",
        ),
        cwd=calls,
        logs=logs,
    )
    _run(
        "merge-normal-shifted-vcfs",
        _gatk(
            "8g",
            "MergeVcfs",
            "--INPUT",
            str(shifted_back),
            "--INPUT",
            str(normal_vcf),
            "--OUTPUT",
            str(merged),
            "--CREATE_INDEX",
            "true",
        ),
        cwd=calls,
        logs=logs,
    )
    _run(
        "merge-mutect-stats",
        _gatk(
            "8g",
            "MergeMutectStats",
            "--stats",
            f"{shifted_vcf}.stats",
            "--stats",
            f"{normal_vcf}.stats",
            "--output",
            str(merged_stats),
        ),
        cwd=calls,
        logs=logs,
    )
    filtered_intermediate = calls / "filtered.intermediate.vcf"
    filter_arguments = [
        "--reference",
        str(resources["normal_fasta"]),
        "--variant",
        str(merged),
        "--stats",
        str(merged_stats),
        "--output",
        str(filtered_intermediate),
        "--max-alt-allele-count",
        str(policy["max_alt_allele_count"]),
        "--mitochondria-mode",
        "--min-allele-fraction",
        format(policy["min_allele_fraction"], ".17g"),
    ]
    if policy["contamination_estimate"] is not None:
        filter_arguments.extend(
            [
                "--contamination-estimate",
                format(policy["contamination_estimate"], ".17g"),
            ]
        )
    _run(
        "filter-mutect-calls",
        _gatk("8g", "FilterMutectCalls", *filter_arguments),
        cwd=calls,
        logs=logs,
    )
    final_vcf = calls / "filtered.vcf"
    _run(
        "mask-blacklisted-sites",
        _gatk(
            "8g",
            "VariantFiltration",
            "--variant",
            str(filtered_intermediate),
            "--output",
            str(final_vcf),
            "--apply-allele-specific-filters",
            "--mask",
            str(resources["blacklist"]),
            "--mask-name",
            "blacklisted_site",
        ),
        cwd=calls,
        logs=logs,
    )
    _run(
        "validate-final-reference",
        [
            "bcftools",
            "norm",
            "--check-ref",
            "e",
            "--fasta-ref",
            str(resources["normal_fasta"]),
            "--output-type",
            "v",
            "--output",
            "/dev/null",
            str(final_vcf),
        ],
        cwd=calls,
        logs=logs,
    )
    try:
        rejected_records = sum(
            1
            for line in rejected.read_text(encoding="ascii").splitlines()
            if line and not line.startswith("#")
        )
    except (OSError, UnicodeError) as exc:
        raise MitochondrialWorkflowError("liftover_reject_file_unreadable") from exc
    if rejected_records:
        _fail("shifted_liftover_rejected_records")
    _regular_file(final_vcf, "final_mt_vcf_missing")
    record_stream_sha256, record_count = _vcf_record_stream(final_vcf)
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "caller": "GATK Mutect2 mitochondrial mode",
        "caller_version": "4.6.2.0",
        "source_contig": policy["source_contig"],
        "mt_length": MT_LENGTH,
        "shift_bases": SHIFT_BASES,
        "normal_interval": NORMAL_INTERVAL,
        "shifted_interval": SHIFTED_INTERVAL,
        "normal_reference_md5": MT_MD5,
        "shifted_reference_md5": SHIFTED_MT_MD5,
        "contamination_estimate_supplied": policy["contamination_estimate"] is not None,
        "liftover_rejected_records": 0,
        "final_vcf_sha256": _sha256_file(final_vcf),
        "record_stream_sha256": record_stream_sha256,
        "records": record_count,
        "private_paths_persisted": False,
        "subject_identifier_persisted": False,
    }
    report = calls / "workflow-report.json"
    report.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--alignment-index", type=Path, required=True)
    parser.add_argument("--full-reference", type=Path, required=True)
    parser.add_argument("--full-reference-fai", type=Path, required=True)
    parser.add_argument("--full-reference-dict", type=Path, required=True)
    parser.add_argument("--mt-resource-root", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    arguments = parser.parse_args(argv)
    try:
        summary = run_workflow(
            alignment=arguments.alignment,
            alignment_index=arguments.alignment_index,
            full_reference=arguments.full_reference,
            full_reference_fai=arguments.full_reference_fai,
            full_reference_dict=arguments.full_reference_dict,
            mt_resource_root=arguments.mt_resource_root,
            policy_path=arguments.policy,
            workspace=arguments.workspace,
        )
    except MitochondrialWorkflowError as exc:
        print(f"mitochondrial workflow failed: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
