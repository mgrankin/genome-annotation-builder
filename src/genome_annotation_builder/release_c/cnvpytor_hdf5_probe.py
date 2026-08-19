"""Run inside the pinned CNVpytor image to audit scientific HDF5 signals."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any

CALL_DATASET = re.compile(r"^calls_([A-Za-z0-9.]+)_([0-9]+)(.*)$")
DEFAULT_REQUIRED_TEMPLATES = (
    "{contig}_rd_p",
    "{contig}_rd_u",
    "his_rd_p_{contig}_{bin_size}_partition{flags_suffix}_merge",
    "calls_{contig}_{bin_size}{flags_suffix}",
)


class ProbeError(ValueError):
    """The CNVpytor HDF5 artifact violates the declared run boundary."""


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


def _canonical_numeric_row(row: list[float]) -> str:
    return "\t".join(float(value).hex() for value in row)


def audit_hdf5(
    input_path: Path,
    *,
    bin_size: int,
    flags_suffix: str,
    expected_contigs: list[str],
    required_templates: list[str],
) -> dict[str, Any]:
    try:
        import h5py
    except ImportError as exc:
        raise ProbeError("h5py is required inside the pinned CNVpytor runtime") from exc
    if bin_size <= 0 or not expected_contigs:
        raise ProbeError("positive bin size and expected contigs are required")
    if len(set(expected_contigs)) != len(expected_contigs):
        raise ProbeError("duplicate expected CNVpytor contig")
    artifact_sha256 = sha256_file(input_path)
    with h5py.File(input_path, "r") as source:
        version = _decode(source.attrs.get("Version", ""))
        if version != "1.3.2":
            raise ProbeError("CNVpytor HDF5 Version is not exactly 1.3.2")
        if "rd_chromosomes" not in source:
            raise ProbeError("CNVpytor HDF5 lacks rd_chromosomes")
        rd_chromosomes = [_decode(value) for value in source["rd_chromosomes"][...]]
        if set(rd_chromosomes) != set(expected_contigs):
            raise ProbeError("CNVpytor HDF5 RD chromosome boundary is unexpected")
        required_names = {
            template.format(contig=contig, bin_size=bin_size, flags_suffix=flags_suffix)
            for contig in expected_contigs
            for template in required_templates
        }
        missing = sorted(required_names - set(source.keys()))
        if missing:
            raise ProbeError(f"CNVpytor HDF5 lacks required scientific signals: {missing[:5]}")
        call_datasets: dict[str, dict[str, Any]] = {}
        for contig in expected_contigs:
            name = f"calls_{contig}_{bin_size}{flags_suffix}"
            dataset = source[name]
            rows: list[list[float]] = []
            if dataset.size:
                for raw_row in dataset[...]:
                    row = [float(value) for value in raw_row]
                    if len(row) != 13 or int(row[0]) != 12:
                        raise ProbeError("CNVpytor mean-shift call dataset has wrong width")
                    if any(not math.isfinite(value) for value in row):
                        raise ProbeError("CNVpytor call dataset contains a non-finite value")
                    if int(row[1]) not in {-1, 1}:
                        raise ProbeError("CNVpytor mean-shift call has an invalid event type")
                    rows.append(row)
            canonical_rows = sorted(_canonical_numeric_row(row) for row in rows)
            call_datasets[name] = {
                "contig": contig,
                "bin_size": bin_size,
                "flags_suffix": flags_suffix,
                "shape": list(dataset.shape),
                "dtype": str(dataset.dtype),
                "rows": rows,
                "canonical_rows_sha256": hashlib.sha256(
                    "\n".join(canonical_rows).encode()
                ).hexdigest(),
            }
        signal_inventory = {
            name: {"shape": list(source[name].shape), "dtype": str(source[name].dtype)}
            for name in sorted(required_names)
        }
        return {
            "schema_version": 1,
            "status": "PASS",
            "cnvpytor_version": version,
            "pytor_sha256": artifact_sha256,
            "bin_size": bin_size,
            "flags_suffix": flags_suffix,
            "rd_chromosomes": sorted(rd_chromosomes),
            "call_datasets": call_datasets,
            "required_signal_inventory": signal_inventory,
            "root_attribute_names": sorted(str(name) for name in source.attrs),
            "bam_attribute_present": "BAM" in source.attrs,
            "volatile_attribute_values_persisted": False,
            "sample_identifier_persisted": False,
            "private_paths_persisted": False,
        }


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
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
        os.replace(name, path)
    except BaseException:
        Path(name).unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--flags-suffix", required=True)
    parser.add_argument("--expected-contig", action="append", required=True)
    parser.add_argument("--required-dataset-template", action="append")
    args = parser.parse_args()
    payload = audit_hdf5(
        args.input,
        bin_size=args.bin_size,
        flags_suffix=args.flags_suffix,
        expected_contigs=args.expected_contig,
        required_templates=args.required_dataset_template or list(DEFAULT_REQUIRED_TEMPLATES),
    )
    _write_json_atomic(args.output, payload)
    print(json.dumps(payload, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
