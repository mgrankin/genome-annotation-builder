"""Audit CNVpytor 1.3.2 combined RD+BAF scientific HDF5 signals in-container."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

COMBINED_KEYS = 20
COMBINED_MODELS = 10
COMBINED_ROW_WIDTH = 1 + COMBINED_KEYS + 5 * COMBINED_MODELS


class CombinedProbeError(ValueError):
    """The CNVpytor combined-call HDF5 boundary is invalid."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _decode(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _integer(value: float, label: str) -> int:
    if not math.isfinite(value) or not value.is_integer():
        raise CombinedProbeError(f"CNVpytor combined {label} is not an integer")
    return int(value)


def validate_combined_row(raw_row: Iterable[float]) -> list[float]:
    row = [float(value) for value in raw_row]
    if len(row) != COMBINED_ROW_WIDTH or _integer(row[0], "key count") != COMBINED_KEYS:
        raise CombinedProbeError("CNVpytor combined-call row has the wrong width")
    if any(not math.isfinite(value) for value in row):
        raise CombinedProbeError("CNVpytor combined-call row contains a non-finite value")
    event_type = _integer(row[1], "event type")
    if event_type not in {-1, 0, 1}:
        raise CombinedProbeError("CNVpytor combined-call event type is invalid")
    start = _integer(row[2], "start")
    end = _integer(row[3], "end")
    size = _integer(row[4], "size")
    if start < 1 or end < start or size != end - start + 1:
        raise CombinedProbeError("CNVpytor combined-call interval is inconsistent")
    if any(value < 0 for value in row[5:10]):
        raise CombinedProbeError("CNVpytor combined-call depth/p-value/likelihood is negative")
    if any(not 0 <= value <= 1 for value in (row[10], row[18], row[19], row[20])):
        raise CombinedProbeError("CNVpytor combined-call mask/gap fraction is invalid")
    bins = _integer(row[11], "bin count")
    baf = row[12]
    if bins <= 0 or not 0 <= baf <= 0.5:
        raise CombinedProbeError("CNVpytor combined-call BAF/bin evidence is invalid")
    if any(value < 0 for value in (row[13], row[14])):
        raise CombinedProbeError("CNVpytor combined-call reported p-value is negative")
    if _integer(row[15], "segment") < 0:
        raise CombinedProbeError("CNVpytor combined-call segment index is negative")
    if _integer(row[16], "heterozygous count") < 0 or _integer(row[17], "homozygous count") < 0:
        raise CombinedProbeError("CNVpytor combined-call SNP count is negative")
    for index in range(COMBINED_MODELS):
        offset = 1 + COMBINED_KEYS + index * 5
        total = _integer(row[offset], "model total copy number")
        first = _integer(row[offset + 1], "model first haplotype")
        second = _integer(row[offset + 2], "model second haplotype")
        likelihood = row[offset + 3]
        cell_fraction = row[offset + 4]
        if min(total, first, second, likelihood, cell_fraction) < 0:
            raise CombinedProbeError("CNVpytor combined-call model contains a negative value")
        if first + second != total or cell_fraction > 1:
            raise CombinedProbeError("CNVpytor combined-call model genotype is inconsistent")
        if index == 0:
            expected_type = -1 if total < 2 else 1 if total > 2 else 0
            if event_type != expected_type:
                raise CombinedProbeError(
                    "CNVpytor combined-call type differs from its top copy-number model"
                )
    return row


def _canonical_numeric_row(row: list[float]) -> str:
    return "\t".join(value.hex() for value in row)


def required_combined_signal_names(
    contig: str,
    *,
    bin_size: int,
    rd_flags_suffix: str,
    snp_flags_suffix: str,
) -> tuple[set[str], str]:
    """Return the exact CNVpytor 1.3.2 signals consumed by ``call_2d``."""
    names = {
        f"{contig}_rd_p",
        f"{contig}_rd_u",
        f"his_rd_p_{contig}_{bin_size}{rd_flags_suffix}",
        # CNVpytor 1.3.2 call_2d consumes the GC/mask-corrected parity
        # histogram for its RD signal, while qrd_p and qrd_u are read from
        # the unflagged binned parity/unique signals. It does not create a
        # GC-corrected unique-RD histogram.
        f"his_rd_p_{contig}_{bin_size}",
        f"his_rd_u_{contig}_{bin_size}",
        f"snp_bafc_01_{contig}_{bin_size}{snp_flags_suffix}",
        f"snp_bafc_10_{contig}_{bin_size}{snp_flags_suffix}",
        f"snp_bafc_11_{contig}_{bin_size}{snp_flags_suffix}",
        f"snp_baf_{contig}_{bin_size}{snp_flags_suffix}",
    }
    call_name = f"calls_2d_{contig}_{bin_size}{snp_flags_suffix}{rd_flags_suffix}"
    names.add(call_name)
    return names, call_name


def audit_combined_hdf5(
    input_path: Path,
    *,
    bin_size: int,
    rd_flags_suffix: str,
    snp_flags_suffix: str,
    expected_contigs: list[str],
) -> dict[str, Any]:
    try:
        import h5py
    except ImportError as exc:
        raise CombinedProbeError("h5py is required in the pinned CNVpytor runtime") from exc
    if bin_size <= 0 or not expected_contigs or len(set(expected_contigs)) != len(expected_contigs):
        raise CombinedProbeError("invalid CNVpytor combined probe configuration")
    artifact_sha256 = sha256_file(input_path)
    with h5py.File(input_path, "r") as source:
        version = _decode(source.attrs.get("Version", ""))
        if version != "1.3.2":
            raise CombinedProbeError("CNVpytor HDF5 Version is not exactly 1.3.2")
        if "rd_chromosomes" not in source or "snp_chromosomes" not in source:
            raise CombinedProbeError("CNVpytor HDF5 lacks RD/SNP chromosome inventories")
        rd_chromosomes = {_decode(value) for value in source["rd_chromosomes"][...]}
        snp_chromosomes = {_decode(value) for value in source["snp_chromosomes"][...]}
        expected = set(expected_contigs)
        if not expected <= rd_chromosomes or not expected <= snp_chromosomes:
            raise CombinedProbeError("CNVpytor HDF5 does not cover every declared RD+BAF contig")
        required_names: set[str] = set()
        call_datasets: dict[str, dict[str, Any]] = {}
        for contig in expected_contigs:
            names, call_name = required_combined_signal_names(
                contig,
                bin_size=bin_size,
                rd_flags_suffix=rd_flags_suffix,
                snp_flags_suffix=snp_flags_suffix,
            )
            required_names.update(names)
            missing = sorted(names - set(source.keys()))
            if missing:
                raise CombinedProbeError(
                    f"CNVpytor HDF5 lacks required combined signals: {missing[:5]}"
                )
            dataset = source[call_name]
            rows = [validate_combined_row(raw_row) for raw_row in dataset[...]]
            canonical = sorted(_canonical_numeric_row(row) for row in rows)
            call_datasets[call_name] = {
                "contig": contig,
                "bin_size": bin_size,
                "rd_flags_suffix": rd_flags_suffix,
                "snp_flags_suffix": snp_flags_suffix,
                "shape": list(dataset.shape),
                "dtype": str(dataset.dtype),
                "rows": rows,
                "canonical_rows_sha256": hashlib.sha256("\n".join(canonical).encode()).hexdigest(),
            }
        inventory = {
            name: {"shape": list(source[name].shape), "dtype": str(source[name].dtype)}
            for name in sorted(required_names)
        }
        return {
            "schema_version": 1,
            "status": "PASS",
            "caller_model": "combined_rd_baf",
            "cnvpytor_version": version,
            "pytor_sha256": artifact_sha256,
            "bin_size": bin_size,
            "rd_flags_suffix": rd_flags_suffix,
            "snp_flags_suffix": snp_flags_suffix,
            "expected_contigs": sorted(expected_contigs),
            "call_datasets": call_datasets,
            "required_signal_inventory": inventory,
            "root_attribute_names": sorted(str(name) for name in source.attrs),
            "bam_attribute_present": "BAM" in source.attrs,
            "volatile_attribute_values_persisted": False,
            "sample_identifier_persisted": False,
            "private_paths_persisted": False,
        }


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    if path.exists():
        raise CombinedProbeError("refusing to replace an existing CNVpytor HDF5 audit")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{path.name}.partial.", suffix=".json", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            os.fchmod(handle.fileno(), 0o644)
        try:
            os.link(name, path)
        except FileExistsError as exc:
            raise CombinedProbeError("refusing to replace an existing CNVpytor HDF5 audit") from exc
        Path(name).unlink()
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        Path(name).unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--rd-flags-suffix", required=True)
    parser.add_argument("--snp-flags-suffix", required=True)
    parser.add_argument("--expected-contig", action="append", required=True)
    args = parser.parse_args()
    try:
        payload = audit_combined_hdf5(
            args.input,
            bin_size=args.bin_size,
            rd_flags_suffix=args.rd_flags_suffix,
            snp_flags_suffix=args.snp_flags_suffix,
            expected_contigs=args.expected_contig,
        )
        _write_json_atomic(args.output, payload)
        print(json.dumps(payload, sort_keys=True, indent=2))
        return 0
    except Exception:
        print(
            "CNVpytor HDF5 audit failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
