"""Run CNVpytor combined RD+BAF calling and capture only validated call rows."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REGION = re.compile(r"^[^:\s]+:([0-9]+)-([0-9]+)$")
GENOTYPE = re.compile(r"^CN([0-9]+)/CN([0-9]+)$")
ERROR_TEXT = re.compile(r"\berror\b", re.IGNORECASE)
EVENT_TYPES = frozenset({"deletion", "duplication", "cnnloh"})


class CnvpytorCaptureError(ValueError):
    """The CNVpytor process or its textual call boundary is invalid."""


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _finite(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise CnvpytorCaptureError(f"invalid CNVpytor combined {label}") from exc
    if not math.isfinite(parsed):
        raise CnvpytorCaptureError(f"non-finite CNVpytor combined {label}")
    return parsed


def _integer(value: str, label: str, *, minimum: int = 0) -> int:
    parsed = _finite(value, label)
    if not parsed.is_integer() or parsed < minimum:
        raise CnvpytorCaptureError(f"invalid CNVpytor combined {label}")
    return int(parsed)


def validate_call_line(line: str, *, expected_bin_size: int) -> str:
    """Validate the exact 27-column CNVpytor 1.3.2 combined-call text row."""

    if not line or "\r" in line or "\x00" in line:
        raise CnvpytorCaptureError("invalid CNVpytor combined call line")
    fields = line.split("\t")
    if len(fields) != 27:
        raise CnvpytorCaptureError("CNVpytor combined call does not have 27 columns")
    if fields[0] not in EVENT_TYPES:
        raise CnvpytorCaptureError("CNVpytor combined call type is invalid")
    region = REGION.fullmatch(fields[1])
    if region is None:
        raise CnvpytorCaptureError("CNVpytor combined interval is invalid")
    start, end = int(region.group(1)), int(region.group(2))
    size = _integer(fields[2], "size", minimum=1)
    if start < 1 or end < start or size != end - start + 1:
        raise CnvpytorCaptureError("CNVpytor combined interval is inconsistent")
    for index, label in (
        (3, "normalized depth"),
        (4, "combined p-value"),
        (5, "deletion likelihood"),
        (6, "LOH likelihood"),
        (7, "duplication likelihood"),
        (8, "Q0 fraction"),
        (9, "reference-gap fraction"),
        (10, "segment-gap fraction"),
        (11, "strict-mask fraction"),
        (14, "BAF deviation"),
        (15, "reported RD p-value"),
        (16, "BAF p-value"),
        (21, "first-model likelihood"),
        (22, "first-model cell fraction"),
        (25, "second-model likelihood"),
        (26, "second-model cell fraction"),
    ):
        if _finite(fields[index], label) < 0:
            raise CnvpytorCaptureError(f"negative CNVpytor combined {label}")
    for index in (8, 9, 10, 11, 22, 26):
        if float(fields[index]) > 1:
            raise CnvpytorCaptureError("CNVpytor combined fraction exceeds one")
    if float(fields[14]) > 0.5:
        raise CnvpytorCaptureError("CNVpytor combined BAF deviation exceeds one half")
    bin_size = _integer(fields[12], "bin size", minimum=1)
    _integer(fields[13], "bin count", minimum=1)
    _integer(fields[17], "heterozygous count")
    _integer(fields[18], "homozygous count")
    if bin_size != expected_bin_size:
        raise CnvpytorCaptureError("CNVpytor combined bin size differs from the request")
    models: list[int] = []
    for total_index, genotype_index, label in ((19, 20, "first"), (23, 24, "second")):
        total = _integer(fields[total_index], f"{label}-model copy number")
        genotype = GENOTYPE.fullmatch(fields[genotype_index])
        if genotype is None:
            raise CnvpytorCaptureError("CNVpytor combined model genotype is invalid")
        if int(genotype.group(1)) + int(genotype.group(2)) != total:
            raise CnvpytorCaptureError("CNVpytor combined model genotype is inconsistent")
        models.append(total)
    expected_type = "deletion" if models[0] < 2 else "duplication" if models[0] > 2 else "cnnloh"
    if fields[0] != expected_type:
        raise CnvpytorCaptureError("CNVpytor combined event type differs from its top model")
    return line


def _regular_input(path: Path, label: str) -> None:
    if not path.is_absolute() or path.is_symlink() or not path.is_file():
        raise CnvpytorCaptureError(label)
    try:
        if path.resolve(strict=True) != path:
            raise CnvpytorCaptureError(label)
    except OSError as exc:
        raise CnvpytorCaptureError(label) from exc


def _output_parent(path: Path) -> None:
    if not path.is_absolute() or path.exists():
        raise CnvpytorCaptureError("CNVpytor combined output already exists or is not absolute")
    parent = path.parent
    if parent.is_symlink() or not parent.is_dir():
        raise CnvpytorCaptureError("CNVpytor combined output parent is invalid")
    try:
        if parent.resolve(strict=True) != parent:
            raise CnvpytorCaptureError("CNVpytor combined output parent is symlinked")
    except OSError as exc:
        raise CnvpytorCaptureError("CNVpytor combined output parent is unreadable") from exc


def _atomic_no_replace(path: Path, payload: bytes) -> None:
    descriptor, raw_name = tempfile.mkstemp(prefix=f".{path.name}.partial.", dir=path.parent)
    temporary = Path(raw_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as exc:
            raise CnvpytorCaptureError("refusing to replace CNVpytor combined output") from exc
        temporary.unlink()
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def capture_combined_calls(
    pytor_path: str | Path,
    config_path: str | Path,
    output_path: str | Path,
    *,
    bin_size: int,
    max_cores: int,
    expected_contigs: list[str],
) -> dict[str, object]:
    """Execute one exact combined call and atomically retain validated TSV rows."""

    pytor = Path(pytor_path)
    config = Path(config_path)
    output = Path(output_path)
    _regular_input(pytor, "CNVpytor HDF5 input is invalid")
    _regular_input(config, "CNVpytor reference configuration is invalid")
    _output_parent(output)
    if bin_size <= 0 or max_cores <= 0:
        raise CnvpytorCaptureError("CNVpytor combined numeric policy is invalid")
    if (
        not expected_contigs
        or len(set(expected_contigs)) != len(expected_contigs)
        or any(not re.fullmatch(r"[0-9]+", contig) for contig in expected_contigs)
    ):
        raise CnvpytorCaptureError("CNVpytor expected contig policy is invalid")
    command = [
        "cnvpytor",
        "-conf",
        str(config),
        "-root",
        str(pytor),
        "-j",
        str(max_cores),
        "-chrom",
        *expected_contigs,
        "-rd_use_mask",
        "-call",
        "combined",
        str(bin_size),
    ]
    process = subprocess.run(command, check=False, capture_output=True)
    if (
        process.returncode != 0
        or ERROR_TEXT.search(process.stdout.decode("utf-8", errors="replace"))
        or ERROR_TEXT.search(process.stderr.decode("utf-8", errors="replace"))
    ):
        raise CnvpytorCaptureError("CNVpytor combined caller failed")
    try:
        stdout = process.stdout.decode("utf-8", errors="strict")
        stderr = process.stderr.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise CnvpytorCaptureError("CNVpytor combined caller emitted non-UTF-8 text") from exc
    lines = stdout.splitlines()
    validated = [validate_call_line(line, expected_bin_size=bin_size) for line in lines]
    payload = (("\n".join(validated) + "\n") if validated else "").encode("ascii")
    _atomic_no_replace(output, payload)
    if stderr:
        sys.stderr.write(stderr)
    return {
        "schema_version": 1,
        "status": "PASS",
        "caller": "CNVpytor",
        "caller_version": "1.3.2",
        "caller_model": "combined_rd_baf",
        "bin_size": bin_size,
        "expected_contigs": len(expected_contigs),
        "call_rows": len(validated),
        "stdout_sha256": _sha256_bytes(process.stdout),
        "stderr_sha256": _sha256_bytes(process.stderr),
        "output_sha256": _sha256_file(output),
        "output_bytes": output.stat().st_size,
        "private_paths_persisted": False,
        "sample_identifier_persisted": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pytor", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--max-cores", type=int, required=True)
    parser.add_argument("--expected-contig", action="append", required=True)
    args = parser.parse_args()
    try:
        report = capture_combined_calls(
            args.pytor,
            args.config,
            args.output,
            bin_size=args.bin_size,
            max_cores=args.max_cores,
            expected_contigs=args.expected_contig,
        )
        print(json.dumps(report, sort_keys=True, indent=2))
        return 0
    except Exception:
        print(
            "CNVpytor combined-call capture failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
