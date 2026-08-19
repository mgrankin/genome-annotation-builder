#!/usr/bin/env python3
"""Create a non-biological all-reference VCF for an end-to-end PGx smoke test.

This utility is deliberately test-only: every official PharmCAT position is
written as an explicit, passing diploid reference observation.  It must never
be used as a substitute for a real sample or as evidence about any person.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import pysam


def _reference_contig(source: str, references: set[str]) -> str | None:
    token = source[3:] if source.startswith("chr") else source
    candidates = [source, token, f"chr{token}"]
    if token in {"M", "MT"}:
        candidates.extend(["M", "MT", "chrM", "chrMT"])
    return next((candidate for candidate in candidates if candidate in references), None)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positions", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    positions = args.positions.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    output = args.output.expanduser().resolve()
    index = Path(f"{output}.tbi")
    if output.suffixes[-2:] != [".vcf", ".gz"]:
        raise ValueError("--output must end in .vcf.gz")
    if not positions.is_file() or not reference.is_file() or not Path(f"{reference}.fai").is_file():
        raise FileNotFoundError("positions, reference FASTA, or FASTA index is missing")
    if not args.force and (output.exists() or index.exists()):
        raise FileExistsError("output or index already exists; pass --force to replace test output")
    output.parent.mkdir(parents=True, exist_ok=True)

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.partial.", suffix=".vcf.gz", dir=output.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    temporary.unlink()
    temporary_index = Path(f"{temporary}.tbi")
    fasta = pysam.FastaFile(str(reference))
    source = pysam.VariantFile(str(positions))
    try:
        header = pysam.VariantHeader()
        header.add_meta("fileformat", value="VCFv4.3")
        header.add_meta("source", value="genome-annotation-builder synthetic PGx smoke fixture")
        header.add_meta(
            "GAB_SYNTHETIC",
            value="ALL_REFERENCE_NON_BIOLOGICAL_TEST_DATA_DO_NOT_INTERPRET",
        )
        reference_names = set(fasta.references)
        contig_map: dict[str, str] = {}
        for contig in source.header.contigs:
            reference_contig = _reference_contig(contig, reference_names)
            if reference_contig is None:
                raise ValueError(f"positions contig is absent from reference: {contig}")
            contig_map[contig] = reference_contig
            if reference_contig not in header.contigs:
                header.contigs.add(
                    reference_contig, length=fasta.get_reference_length(reference_contig)
                )
        header.add_meta(
            "ALT",
            items=[
                ("ID", "NON_REF"),
                ("Description", "Represents any possible alternative allele at this location"),
            ],
        )
        header.formats.add("GT", 1, "String", "Synthetic explicit reference genotype")
        header.formats.add("DP", 1, "Integer", "Synthetic smoke-test depth")
        header.formats.add("GQ", 1, "Integer", "Synthetic smoke-test genotype quality")
        header.formats.add("MIN_DP", 1, "Integer", "Synthetic smoke-test minimum depth")
        header.add_sample("SYNTHETIC_PGX_REFERENCE")

        writer = pysam.VariantFile(str(temporary), "wz", header=header)
        try:
            count = 0
            for source_record in source:
                output_contig = contig_map[source_record.contig]
                ref = str(source_record.ref).upper()
                observed = fasta.fetch(
                    output_contig,
                    source_record.pos - 1,
                    source_record.pos - 1 + len(ref),
                ).upper()
                if observed != ref:
                    raise ValueError(
                        f"REF mismatch at {source_record.contig}:{source_record.pos}: {ref} != {observed}"
                    )
                record = writer.new_record(
                    contig=output_contig,
                    start=source_record.pos - 1,
                    stop=source_record.pos - 1 + len(ref),
                    alleles=(ref, "<NON_REF>"),
                )
                record.filter.add("PASS")
                sample = record.samples["SYNTHETIC_PGX_REFERENCE"]
                sample["GT"] = (0, 0)
                sample["DP"] = 30
                sample["GQ"] = 60
                sample["MIN_DP"] = 30
                writer.write(record)
                count += 1
        finally:
            writer.close()
        if count == 0:
            raise ValueError("positions resource contains no records")
        pysam.tabix_index(str(temporary), preset="vcf", force=True)
        os.replace(temporary, output)
        os.replace(temporary_index, index)
    finally:
        source.close()
        fasta.close()
        if temporary.exists():
            temporary.unlink()
        if temporary_index.exists():
            temporary_index.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
