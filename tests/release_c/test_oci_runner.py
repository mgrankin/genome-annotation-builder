from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import oci_runner
from genome_annotation_builder.release_c.caller_receipt import (
    canonical_json_sha256,
    load_prepared_caller_run_boundary,
    write_caller_run_receipt,
    write_prepared_caller_run_boundary,
)
from genome_annotation_builder.release_c.oci_runner import (
    OciExecutionFailure,
    load_oci_execution_report,
    run_oci_caller,
    validate_oci_execution_report_payload,
)
from genome_annotation_builder.release_c.resource_tree import create_resource_tree_manifest

TOOL_IDENTITY = "4" * 64
GATE_RUN = "5" * 64
PRIVATE_BINDING = "privateid"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def file_mount(role: str, host: Path, container: str) -> dict[str, object]:
    return {
        "role": role,
        "kind": "file",
        "host_path": str(host),
        "container_path": container,
        "mode": "ro",
        "bytes": host.stat().st_size,
        "sha256": sha256_file(host),
        "tree_manifest_path": None,
        "tree_manifest_bytes": None,
        "tree_manifest_sha256": None,
    }


def workspace_mount(workspace: Path) -> dict[str, object]:
    return {
        "role": "workspace",
        "kind": "directory",
        "host_path": str(workspace),
        "container_path": "/work",
        "mode": "rw",
        "bytes": None,
        "sha256": None,
        "tree_manifest_path": None,
        "tree_manifest_bytes": None,
        "tree_manifest_sha256": None,
    }


def tree_mount(tmp_path: Path) -> dict[str, object]:
    root = tmp_path / "catalog"
    (root / "nested").mkdir(parents=True)
    (root / "catalog.json").write_text('{"loci": []}\n', encoding="utf-8")
    (root / "nested" / "index.bin").write_bytes(b"synthetic-index\n")
    manifest_path = tmp_path / "catalog-tree-manifest.json"
    create_resource_tree_manifest(root, manifest_path)
    return {
        "role": "resource_tree",
        "kind": "directory",
        "host_path": str(root),
        "container_path": "/resources/catalog",
        "mode": "ro",
        "bytes": None,
        "sha256": None,
        "tree_manifest_path": str(manifest_path),
        "tree_manifest_bytes": manifest_path.stat().st_size,
        "tree_manifest_sha256": sha256_file(manifest_path),
    }


def prepared_fixture(tmp_path: Path) -> dict[str, object]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    alignment = tmp_path / "alignment.bam"
    reference = tmp_path / "reference.fa"
    resources = tmp_path / "resource-manifest.json"
    alignment.write_bytes(b"synthetic-alignment\n")
    reference.write_text(">1\nACGTACGT\n", encoding="ascii")
    write_json(resources, {"schema_version": 1, "status": "PASS"})
    gate = tmp_path / "input-gate.json"
    write_json(
        gate,
        {
            "schema_version": 1,
            "status": "PASS",
            "caller": "manta",
            "selected_alignment": {"sha256": sha256_file(alignment)},
            "reference": {"sha256": sha256_file(reference)},
            "run_identity_sha256": GATE_RUN,
            "logical_subject_count": 1,
            "subject_identifier_persisted": False,
            "private_paths_persisted": False,
            "private_manifest_digest_persisted": False,
        },
    )
    command = {
        "steps": [
            {
                "name": "configure",
                "entrypoint": "configManta.py",
                "arguments": [
                    "--bam",
                    "{selected_alignment}",
                    "--referenceFasta",
                    "{reference_fasta}",
                    "--runDir",
                    "{working_directory}",
                ],
            },
            {
                "name": "call",
                "entrypoint": "runWorkflow.py",
                "arguments": ["--mode", "local", "--jobs", "2"],
            },
        ],
        "environment": {"LC_ALL": "C", "TMPDIR": "{working_directory}/tmp"},
        "parameters": {"threads": 2},
    }
    boundary = tmp_path / "prepared-boundary.json"
    write_prepared_caller_run_boundary(
        boundary,
        gate,
        resources,
        caller_name="Manta",
        caller_version="1.6.0",
        gate_caller="manta",
        expected_resource_manifest_sha256=sha256_file(resources),
        tool={
            "execution_kind": "oci_image",
            "logical_name": "manta",
            "identity_sha256": TOOL_IDENTITY,
            "network_disabled": True,
        },
        command_semantics=command,
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    mounts = [
        file_mount("selected_alignment", alignment, "/inputs/alignment.bam"),
        file_mount("reference_fasta", reference, "/inputs/reference.fa"),
        file_mount("resource_manifest", resources, "/provenance/resources.json"),
        tree_mount(tmp_path),
        workspace_mount(workspace),
    ]
    specification = {
        "schema_version": 1,
        "runtime": "docker",
        "image": f"example/manta@sha256:{TOOL_IDENTITY}",
        "user": "1000:1000",
        "mounts": mounts,
        "placeholder_bindings": {
            "{selected_alignment}": "/inputs/alignment.bam",
            "{reference_fasta}": "/inputs/reference.fa",
            "{working_directory}": "/work",
        },
        "tmpfs": [{"container_path": "/tmp", "bytes": 64 * 1024**2, "executable": False}],
        "limits": {
            "pids": 256,
            "memory_bytes": 2 * 1024**3,
            "cpus": 2.0,
            "shm_bytes": 128 * 1024**2,
        },
        "minimum_free_bytes": 1,
    }
    spec_path = tmp_path / "execution-spec.json"
    write_json(spec_path, specification)
    return {
        "boundary": boundary,
        "gate": gate,
        "alignment": alignment,
        "reference": reference,
        "resources": resources,
        "workspace": workspace,
        "specification": specification,
        "spec_path": spec_path,
    }


def install_fake_runtime(
    monkeypatch: pytest.MonkeyPatch,
    commands: list[list[str]],
    *,
    exit_codes: list[int] | None = None,
    mutate_after_last: Path | None = None,
) -> None:
    def capture(command: list[str]) -> subprocess.CompletedProcess[str]:
        if command[1:] == ["--version"]:
            return subprocess.CompletedProcess(command, 0, "Docker version synthetic\n", "")
        return subprocess.CompletedProcess(
            command,
            0,
            json.dumps([f"example/manta@sha256:{TOOL_IDENTITY}"]) + "\n",
            "",
        )

    codes = iter(exit_codes or [0, 0])

    def execute(command: list[str], stdout_handle, stderr_handle) -> int:
        commands.append(command)
        stdout_handle.write(b"synthetic stdout\n")
        stderr_handle.write(b"")
        code = next(codes)
        if mutate_after_last is not None and len(commands) == 2:
            mutate_after_last.write_bytes(mutate_after_last.read_bytes() + b"changed")
        return code

    monkeypatch.setattr(oci_runner, "_capture", capture)
    monkeypatch.setattr(oci_runner, "_execute", execute)


def test_oci_runner_enforces_security_and_emits_portable_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = prepared_fixture(tmp_path)
    commands: list[list[str]] = []
    install_fake_runtime(monkeypatch, commands)
    report_path = tmp_path / "execution-report.json"
    log_dir = tmp_path / "private-logs"
    report = run_oci_caller(
        fixture["boundary"],
        fixture["spec_path"],
        report_path,
        log_dir,
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    assert report["status"] == "PASS"
    assert report["error_code"] is None
    assert report["runtime"]["image_digest_verified"] is True
    assert report["input_verification"]["post_execution_stat_unchanged"] is True
    assert len(report["steps"]) == 2
    assert len(commands) == 2
    for command in commands:
        assert command[:4] == ["docker", "run", "--rm", "--pull"]
        for required_pair in (
            ("--pull", "never"),
            ("--network", "none"),
            ("--hostname", "gab-caller"),
            ("--add-host", "gab-caller:127.0.0.1"),
            ("--cap-drop", "ALL"),
            ("--security-opt", "no-new-privileges"),
            ("--user", "1000:1000"),
        ):
            index = command.index(required_pair[0])
            assert command[index : index + 2] == list(required_pair)
        assert "--read-only" in command
        assert subprocess.DEVNULL not in command
        assert "/bin/sh" not in command
        assert "-c" not in command
        assert command.count("--mount") == 5
        assert any(
            value.endswith(",readonly") for value in command if value.startswith("type=bind")
        )
    loaded, report_sha256 = load_oci_execution_report(
        report_path, forbidden_bindings=(PRIVATE_BINDING,)
    )
    assert loaded == report
    assert report_sha256 == sha256_file(report_path)
    assert len(list(log_dir.glob("*.log"))) == 4
    payload = report_path.read_text(encoding="utf-8")
    forbidden_values = {
        PRIVATE_BINDING,
        *(str(value) for value in fixture.values() if isinstance(value, Path)),
        str(report_path),
        str(log_dir),
    }
    assert not any(value in payload for value in forbidden_values)
    checkpoint, _ = load_prepared_caller_run_boundary(fixture["boundary"])
    raw_output = fixture["workspace"] / "diploidSV.vcf.gz"
    raw_output.write_bytes(b"synthetic caller output\n")
    receipt = write_caller_run_receipt(
        tmp_path / "caller-receipt.json",
        checkpoint["run_boundary"],
        fixture["gate"],
        fixture["resources"],
        execution_report=report_path,
        raw_artifact_id="manta-diploid-sv-vcf",
        identity_scheme="file_sha256",
        component_paths={"primary_vcf": raw_output},
    )
    assert (
        receipt["execution"]["report"]["execution_identity_sha256"]
        == report["execution_identity_sha256"]
    )


def test_negative_signal_exit_is_retained_without_invalid_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = prepared_fixture(tmp_path)
    commands: list[list[str]] = []
    install_fake_runtime(monkeypatch, commands, exit_codes=[0, -9])
    report_path = tmp_path / "failed-report.json"
    with pytest.raises(OciExecutionFailure, match="caller_step_failed") as raised:
        run_oci_caller(
            fixture["boundary"],
            fixture["spec_path"],
            report_path,
            tmp_path / "failed-logs",
        )
    assert raised.value.report_written is True
    report, _ = load_oci_execution_report(report_path)
    assert report["status"] == "FAIL_RETAINED"
    assert report["error_code"] == "caller_step_failed"
    assert [step["exit_code"] for step in report["steps"]] == [0, -9]
    assert report["input_verification"]["post_execution_stat_unchanged"] is False


def test_input_change_during_execution_is_detected_and_retained(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = prepared_fixture(tmp_path)
    commands: list[list[str]] = []
    install_fake_runtime(
        monkeypatch,
        commands,
        mutate_after_last=fixture["alignment"],
    )
    report_path = tmp_path / "changed-input-report.json"
    with pytest.raises(OciExecutionFailure, match="read_only_input_changed"):
        run_oci_caller(
            fixture["boundary"],
            fixture["spec_path"],
            report_path,
            tmp_path / "changed-input-logs",
        )
    report, _ = load_oci_execution_report(report_path)
    assert report["status"] == "FAIL_RETAINED"
    assert report["error_code"] == "read_only_input_changed_during_execution"
    assert len(report["steps"]) == 2


def test_mount_digest_directory_inventory_and_privacy_fail_before_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = prepared_fixture(tmp_path)
    commands: list[list[str]] = []
    install_fake_runtime(monkeypatch, commands)
    tampered_spec = json.loads(fixture["spec_path"].read_text(encoding="utf-8"))
    tampered_spec["mounts"][0]["sha256"] = "0" * 64
    write_json(fixture["spec_path"], tampered_spec)
    with pytest.raises(OciExecutionFailure, match="digest_or_size_changed"):
        run_oci_caller(
            fixture["boundary"],
            fixture["spec_path"],
            tmp_path / "tampered-report.json",
            tmp_path / "tampered-logs",
        )
    assert commands == []
    assert not (tmp_path / "tampered-report.json").exists()
    assert not (tmp_path / "tampered-logs").exists()

    second = prepared_fixture(tmp_path / "tree")
    tree_file = tmp_path / "tree" / "catalog" / "catalog.json"
    tree_file.write_text('{"loci": ["changed"]}\n', encoding="utf-8")
    with pytest.raises(OciExecutionFailure, match="read_only_directory_manifest_validation_failed"):
        run_oci_caller(
            second["boundary"],
            second["spec_path"],
            tmp_path / "tree" / "tampered-report.json",
            tmp_path / "tree" / "tampered-logs",
        )

    third = prepared_fixture(tmp_path / "privacy")
    private_spec = json.loads(third["spec_path"].read_text(encoding="utf-8"))
    private_spec["mounts"][3]["role"] = f"artifact-{PRIVATE_BINDING}"
    write_json(third["spec_path"], private_spec)
    with pytest.raises(OciExecutionFailure, match="private_binding"):
        run_oci_caller(
            third["boundary"],
            third["spec_path"],
            tmp_path / "privacy" / "private-report.json",
            tmp_path / "privacy" / "private-logs",
            forbidden_bindings=(PRIVATE_BINDING,),
        )


def test_independent_execution_report_validator_detects_semantic_tamper(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = prepared_fixture(tmp_path)
    install_fake_runtime(monkeypatch, [])
    report_path = tmp_path / "report.json"
    report = run_oci_caller(
        fixture["boundary"],
        fixture["spec_path"],
        report_path,
        tmp_path / "logs",
    )
    tampered = json.loads(report_path.read_text(encoding="utf-8"))
    tampered["security"]["network_none"] = False
    without_identity = {
        key: value for key, value in tampered.items() if key != "execution_identity_sha256"
    }
    tampered["execution_identity_sha256"] = canonical_json_sha256(without_identity)
    with pytest.raises(OciExecutionFailure, match="security_invalid"):
        validate_oci_execution_report_payload(tampered)

    digest_tamper = dict(report)
    digest_tamper["steps"] = [dict(step) for step in report["steps"]]
    digest_tamper["steps"][0]["stdout_sha256"] = "0" * 64
    with pytest.raises(OciExecutionFailure, match="identity_changed"):
        validate_oci_execution_report_payload(digest_tamper)
