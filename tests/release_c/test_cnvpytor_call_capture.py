from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import cnvpytor_call_capture
from genome_annotation_builder.release_c.cnvpytor_call_capture import (
    CnvpytorCaptureError,
    capture_combined_calls,
    validate_call_line,
)


def valid_line(*, bin_size: int = 10000) -> str:
    fields = [
        "deletion",
        "1:100001-160000",
        "60000",
        "0.5",
        "1e-8",
        "0.8",
        "0.1",
        "0.1",
        "0.1",
        "0.0",
        "0.0",
        "0.9",
        str(bin_size),
        "6",
        "0.2",
        "1e-8",
        "1e-6",
        "20",
        "5",
        "1",
        "CN1/CN0",
        "0.9",
        "1.0",
        "2",
        "CN1/CN1",
        "0.1",
        "1.0",
    ]
    assert len(fields) == 27
    return "\t".join(fields)


def fixture_paths(tmp_path: Path) -> tuple[Path, Path, Path]:
    pytor = tmp_path / "calls.pytor"
    config = tmp_path / "reference_genomes_conf.py"
    output = tmp_path / "combined.tsv"
    pytor.write_bytes(b"synthetic HDF5 boundary")
    config.write_text("import_reference_genomes = {}\n", encoding="utf-8")
    return pytor, config, output


def completed(
    stdout: bytes, *, stderr: bytes = b"", returncode: int = 0
) -> subprocess.CompletedProcess[bytes]:
    return subprocess.CompletedProcess(["cnvpytor"], returncode, stdout, stderr)


def test_exact_row_boundary_accepts_expected_combined_call() -> None:
    line = valid_line()
    assert validate_call_line(line, expected_bin_size=10000) == line


@pytest.mark.parametrize(
    "mutation",
    ["wrong_width", "wrong_bin", "wrong_type", "wrong_model", "non_finite"],
)
def test_row_boundary_rejects_malformed_or_inconsistent_call(mutation: str) -> None:
    fields = valid_line().split("\t")
    if mutation == "wrong_width":
        fields.pop()
    elif mutation == "wrong_bin":
        fields[12] = "100000"
    elif mutation == "wrong_type":
        fields[0] = "duplication"
    elif mutation == "wrong_model":
        fields[20] = "CN1/CN1"
    else:
        fields[4] = "nan"
    with pytest.raises(CnvpytorCaptureError):
        validate_call_line("\t".join(fields), expected_bin_size=10000)


def test_capture_writes_validated_rows_atomically_and_reports_only_hashes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytor, config, output = fixture_paths(tmp_path)
    observed: list[list[str]] = []

    def run(command: list[str], **kwargs) -> subprocess.CompletedProcess[bytes]:
        observed.append(command)
        return completed((valid_line() + "\n").encode(), stderr=b"diagnostic\n")

    monkeypatch.setattr(cnvpytor_call_capture.subprocess, "run", run)
    report = capture_combined_calls(
        pytor,
        config,
        output,
        bin_size=10000,
        max_cores=8,
        expected_contigs=["1", "2"],
    )
    assert output.read_text(encoding="ascii") == valid_line() + "\n"
    assert report["status"] == "PASS"
    assert report["call_rows"] == 1
    assert report["private_paths_persisted"] is False
    retained = json.dumps(report, sort_keys=True)
    assert str(tmp_path) not in retained
    command = observed[0]
    assert command[:7] == [
        "cnvpytor",
        "-conf",
        str(config),
        "-root",
        str(pytor),
        "-j",
        "8",
    ]
    assert command[-3:] == ["-call", "combined", "10000"]

    with pytest.raises(CnvpytorCaptureError, match="already exists"):
        capture_combined_calls(
            pytor,
            config,
            output,
            bin_size=10000,
            max_cores=8,
            expected_contigs=["1", "2"],
        )


def test_empty_callset_is_explicit_valid_zero_row_artifact(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytor, config, output = fixture_paths(tmp_path)
    monkeypatch.setattr(
        cnvpytor_call_capture.subprocess,
        "run",
        lambda *args, **kwargs: completed(b""),
    )
    report = capture_combined_calls(
        pytor,
        config,
        output,
        bin_size=100000,
        max_cores=4,
        expected_contigs=["1"],
    )
    assert output.read_bytes() == b""
    assert report["call_rows"] == 0


@pytest.mark.parametrize(
    ("process", "message"),
    [
        (completed(b"", returncode=1), "caller failed"),
        (completed(b"ERROR: missing resource\n"), "caller failed"),
        (completed(b"", stderr=b"Error opening file\n"), "caller failed"),
        (completed(b"not a call row\n"), "27 columns"),
    ],
)
def test_false_zero_nonzero_and_non_call_output_are_rejected_without_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    process: subprocess.CompletedProcess[bytes],
    message: str,
) -> None:
    pytor, config, output = fixture_paths(tmp_path)
    monkeypatch.setattr(
        cnvpytor_call_capture.subprocess,
        "run",
        lambda *args, **kwargs: process,
    )
    with pytest.raises(CnvpytorCaptureError, match=message):
        capture_combined_calls(
            pytor,
            config,
            output,
            bin_size=10000,
            max_cores=8,
            expected_contigs=["1"],
        )
    assert not output.exists()


def test_cli_failure_is_sanitized_and_does_not_echo_private_path(tmp_path: Path) -> None:
    private_token = "privateid"
    private_root = tmp_path / private_token
    private_root.mkdir()
    script = Path(cnvpytor_call_capture.__file__)
    failed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--pytor",
            str(private_root / "missing.pytor"),
            "--config",
            str(private_root / "missing.py"),
            "--output",
            str(private_root / "calls.tsv"),
            "--bin-size",
            "10000",
            "--max-cores",
            "8",
            "--expected-contig",
            "1",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert failed.returncode == 2
    assert failed.stdout == ""
    assert "Traceback" not in failed.stderr
    assert private_token not in failed.stderr
    assert str(private_root) not in failed.stderr
