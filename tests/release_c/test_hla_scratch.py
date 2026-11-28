"""The matrix may retire only generated HLA*LA mate scratch after typed readback."""

import json
import os
import shutil
import subprocess
from contextlib import ExitStack
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.cli import main
from genome_annotation_builder.release_c import (
    caller_lifecycle,
    hla_scratch,
    matrix_runner,
    oci_runner,
    profile_adapter,
)

from .test_caller_profile import file_binding, make_hla_fixture, write_json
from .test_hla_la_adapter import control_table
from .test_matrix_runner import _inputs, _install_fake_execution, _small_matrix


@pytest.mark.parametrize("valid_hla_options", [True, False])
def test_matrix_cli_hla_validates_before_cleanup_and_replays(
    tmp_path, monkeypatch, valid_hla_options
):
    manifest, resources, freeze, options = _inputs(tmp_path, two_datasets=False)
    config = _small_matrix(tmp_path)
    hla = "hla-la-short-read-1.0.4"
    adapter_options = json.loads(options.read_text())
    adapter_options["profiles"][hla]["database_age_limitation"] = valid_hla_options
    write_json(options, adapter_options)
    cfg = json.loads(config.read_text())
    cfg["fixed_profiles"] = [hla]
    write_json(config, cfg)
    fixture = make_hla_fixture(tmp_path / "source")
    root = tmp_path / "matrix"
    attempt = root / "job-0001/attempt-0001"
    subject = attempt / "caller-workspace/SUBJECT"
    table = control_table(tmp_path)
    original = fixture["alignment"].read_bytes()
    executions = []

    # Use the established stand-in for the unrelated karyotype profile. HLA uses
    # its real profile lifecycle, caller receipt, parser, Parquet writer and validator.
    _install_fake_execution(monkeypatch)
    other_lifecycle = matrix_runner.run_caller_lifecycle
    other_materialize = matrix_runner.materialize_profile_adapter
    other_validate = matrix_runner.validate_profile_adapter_output

    def lifecycle(*args, profile_id, **kwargs):
        if profile_id != hla:
            return other_lifecycle(*args, profile_id=profile_id, **kwargs)
        bindings = json.loads(fixture["bindings"].read_text())
        bindings["mounts"]["workspace"]["host_path"] = str(subject.parent)
        bindings["mounts"]["resource_manifest"] = file_binding(freeze)
        write_json(args[2], bindings)
        work = args[4]
        work.mkdir(exist_ok=True)
        if not (work / "input-gate.json").exists():
            shutil.copyfile(fixture["gate"], work / "input-gate.json")
        return caller_lifecycle.run_caller_lifecycle(
            *args, profile_path=fixture["profile"], **kwargs
        )

    def materialize(*args, profile_id, **kwargs):
        if profile_id == hla:
            assert all((subject / n).is_file() for n in hla_scratch.MATES)
            return profile_adapter.materialize_profile_adapter(*args, profile_id=hla, **kwargs)
        return other_materialize(*args, profile_id=profile_id, **kwargs)

    def validate(*args, profile_id, **kwargs):
        if profile_id == hla:
            return profile_adapter.validate_profile_adapter_output(*args, profile_id=hla, **kwargs)
        return other_validate(*args, profile_id=profile_id, **kwargs)

    image = json.loads(fixture["profile"].read_text())["image"]

    def capture(argv):
        result = "Docker synthetic" if argv[1:] == ["--version"] else json.dumps([image])
        return subprocess.CompletedProcess(argv, 0, result + "\n", "")

    def execute(argv, stdout, stderr):
        executions.append(argv)
        if "HLA-LA.pl" in argv:
            (subject / "hla").mkdir(parents=True)
            shutil.copyfile(table, subject / "hla/R1_bestguess_G.txt")
            shutil.copyfile(fixture["alignment"], subject / "extraction.bam")
            for name in hla_scratch.MATES:
                (subject / name).write_bytes(b"@pair\nAC\n+\nII\n")
            (subject / "R_U.fastq").touch()
        return 0

    monkeypatch.setattr(matrix_runner, "run_caller_lifecycle", lifecycle)
    monkeypatch.setattr(matrix_runner, "materialize_profile_adapter", materialize)
    monkeypatch.setattr(matrix_runner, "validate_profile_adapter_output", validate)
    monkeypatch.setattr(profile_adapter, "resolve_shipped_profile", lambda _: fixture["profile"])
    monkeypatch.setattr(oci_runner, "_capture", capture)
    monkeypatch.setattr(oci_runner, "_execute", execute)
    argv = [
        "release-c",
        "matrix-runner",
        "--manifest",
        str(manifest),
        "--resource-bindings",
        str(resources),
        "--resource-manifest",
        str(freeze),
        "--work-root",
        str(root),
        "--adapter-options",
        str(options),
        "--matrix-config",
        str(config),
        "--user",
        "1000:1000",
        "--minimum-free-bytes",
        "1",
        "--indexed-probe-bases",
        "1",
    ]
    if not valid_hla_options:
        assert main(argv) == 2
        state = json.loads((root / matrix_runner.RUNTIME_NAME).read_text())["jobs"][0]
        assert state["last_completed_stage"] == "caller_lifecycle"
        assert "database-age limitation" in state["error_message"]
        assert all((subject / n).is_file() for n in hla_scratch.MATES)
        return
    assert main(argv) == 0
    assert not any((subject / n).exists() for n in hla_scratch.MATES)
    assert (
        fixture["alignment"].read_bytes() == (subject / "extraction.bam").read_bytes() == original
    )
    assert (subject / "R_U.fastq").is_file()
    assert (subject / "hla/R1_bestguess_G.txt").read_bytes() == table.read_bytes()
    assert pq.read_table(attempt / "typed/hla-la.parquet").num_rows == 13
    # A later file in a retained workspace belongs to its user.
    (subject / hla_scratch.MATES[0]).write_bytes(b"keep on validation")
    assert main([*argv, "--validate-only"]) == 0
    assert (subject / hla_scratch.MATES[0]).read_bytes() == b"keep on validation"
    assert len(executions) == 2


@pytest.mark.parametrize(
    "change", ["symlink", "hardlink", "replaced", "changed-bam", "moved", "declared"]
)
def test_keeps_uncertain_or_changed_scratch(tmp_path: Path, change: str):
    workspace = tmp_path / "workspace"
    original = tmp_path / "retained"
    original.write_bytes(b"retained")
    identity = hla_scratch.new_workspace(workspace)
    assert identity is not None and hla_scratch.new_workspace(workspace) is None
    with ExitStack() as handles:
        subject = workspace / "SUBJECT"
        subject.mkdir()
        source = subject / "extraction.bam"
        source.write_bytes(b"bam")
        first, second = (subject / name for name in hla_scratch.MATES)
        first.write_bytes(b"first")
        if change == "symlink":
            second.symlink_to(original)
        elif change == "hardlink":
            os.link(original, second)
        else:
            second.write_bytes(b"second")
        components = {"bestguess_g_group": str(subject / "hla/R1_bestguess_G.txt")}
        if change == "declared":
            components["retained_reads"] = str(first)
        captured = hla_scratch.snapshot_mates(
            workspace, identity, {"components": components}, handles
        )
        if change == "replaced":
            second.unlink()
            second.write_bytes(b"replacement")
        elif change == "changed-bam":
            source.unlink()
        elif change == "moved":
            subject.rename(tmp_path / "moved")
            subject.symlink_to(tmp_path / "moved", target_is_directory=True)
        hla_scratch.remove_validated_mates(captured)
    assert first.is_file() and second.is_file() and original.read_bytes() == b"retained"
