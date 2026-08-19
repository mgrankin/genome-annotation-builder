from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pytest

from genome_annotation_builder.release_c.caller_receipt import (
    canonical_json_sha256,
    load_caller_run_receipt,
    load_prepared_caller_run_boundary,
    prepare_caller_run_boundary,
    validate_caller_run_receipt_payload,
    validate_prepared_caller_run_boundary_payload,
    write_caller_run_receipt,
    write_prepared_caller_run_boundary,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError

REFERENCE = "1" * 64
ALIGNMENT = "2" * 64
GATE_RUN = "3" * 64
TOOL_IDENTITY = "4" * 64
PRIVATE_SUBJECT = "SYNTHETIC_PRIVATE_SUBJECT"
MODULE = "genome_annotation_builder.release_c.caller_receipt"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def gate_report(tmp_path: Path, *, caller: str = "manta") -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / f"gate-{caller}.json"
    write_json(
        path,
        {
            "schema_version": 1,
            "status": "PASS",
            "caller": caller,
            "selected_alignment": {"sha256": ALIGNMENT},
            "reference": {"sha256": REFERENCE},
            "run_identity_sha256": GATE_RUN,
            "logical_subject_count": 1,
            "subject_identifier_persisted": False,
            "private_paths_persisted": False,
            "private_manifest_digest_persisted": False,
        },
    )
    return path


def resource_manifest(tmp_path: Path) -> Path:
    path = tmp_path / "resources.json"
    write_json(
        path,
        {
            "schema_version": 1,
            "status": "PASS",
            "resource": {"logical_id": "synthetic-reference-track", "sha256": "5" * 64},
        },
    )
    return path


def tool() -> dict[str, object]:
    return {
        "execution_kind": "oci_image",
        "logical_name": "manta",
        "identity_sha256": TOOL_IDENTITY,
        "network_disabled": True,
    }


def command() -> dict[str, object]:
    return {
        "steps": [
            {
                "name": "discover",
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
                "arguments": ["--mode", "local", "--jobs", "1"],
            },
        ],
        "environment": {"LC_ALL": "C"},
        "parameters": {"threads": 1, "output_role": "primary_vcf"},
    }


def boundary(tmp_path: Path) -> tuple[dict[str, object], Path, Path]:
    gate = gate_report(tmp_path)
    resources = resource_manifest(tmp_path)
    prepared = prepare_caller_run_boundary(
        gate,
        resources,
        caller_name="Manta",
        caller_version="1.6.0",
        gate_caller="manta",
        expected_resource_manifest_sha256=sha256_file(resources),
        tool=tool(),
        command_semantics=command(),
        forbidden_bindings=(PRIVATE_SUBJECT,),
    )
    return prepared, gate, resources


def run_cli(arguments: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", MODULE, *arguments],
        check=False,
        capture_output=True,
        text=True,
    )


def write_execution_report(
    tmp_path: Path,
    prepared_boundary: dict[str, object],
    *,
    filename: str = "execution-report.json",
) -> Path:
    empty_sha256 = hashlib.sha256(b"").hexdigest()
    steps = [
        {
            "ordinal": ordinal,
            "name": step["name"],
            "semantic_step_sha256": canonical_json_sha256(step),
            "exit_code": 0,
            "stdout_sha256": empty_sha256,
            "stdout_bytes": 0,
            "stderr_sha256": empty_sha256,
            "stderr_bytes": 0,
        }
        for ordinal, step in enumerate(prepared_boundary["command_semantics"]["steps"], start=1)
    ]
    report = {
        "schema_version": 1,
        "status": "PASS",
        "error_code": None,
        "caller_run_identity_sha256": prepared_boundary["caller_run_identity_sha256"],
        "source_alignment_sha256": prepared_boundary["source_alignment_sha256"],
        "reference_sha256": prepared_boundary["reference_sha256"],
        "resource_manifest_sha256": prepared_boundary["resource_manifest_sha256"],
        "tool": prepared_boundary["tool"],
        "runtime": {
            "name": "docker",
            "version_output_sha256": "a" * 64,
            "image_digest_sha256": prepared_boundary["tool"]["identity_sha256"],
            "image_digest_verified": True,
        },
        "security": {
            "pull_never": True,
            "network_none": True,
            "fixed_loopback_hostname": True,
            "read_only_root": True,
            "cap_drop_all": True,
            "no_new_privileges": True,
            "non_root_numeric_user": True,
            "stdin_disabled": True,
            "tmpfs_count": 1,
            "tmpfs_total_bytes": 64 * 1024**2,
            "tmpfs_all_noexec": True,
            "pids_limit": 256,
            "memory_limit_bytes": 2 * 1024**3,
            "cpu_limit": 2.0,
            "shm_limit_bytes": 128 * 1024**2,
        },
        "input_verification": {
            "mounts": [
                {
                    "role": "reference_fasta",
                    "kind": "file",
                    "mode": "ro",
                    "sha256": prepared_boundary["reference_sha256"],
                    "bytes": 1,
                    "tree_manifest_sha256": None,
                    "tree_identity_sha256": None,
                    "tree_files": None,
                    "tree_directories": None,
                },
                {
                    "role": "resource_manifest",
                    "kind": "file",
                    "mode": "ro",
                    "sha256": prepared_boundary["resource_manifest_sha256"],
                    "bytes": prepared_boundary["resource_manifest_bytes"],
                    "tree_manifest_sha256": None,
                    "tree_identity_sha256": None,
                    "tree_files": None,
                    "tree_directories": None,
                },
                {
                    "role": "selected_alignment",
                    "kind": "file",
                    "mode": "ro",
                    "sha256": prepared_boundary["source_alignment_sha256"],
                    "bytes": 1,
                    "tree_manifest_sha256": None,
                    "tree_identity_sha256": None,
                    "tree_files": None,
                    "tree_directories": None,
                },
                {
                    "role": "workspace",
                    "kind": "directory",
                    "mode": "rw",
                    "sha256": None,
                    "bytes": None,
                    "tree_manifest_sha256": None,
                    "tree_identity_sha256": None,
                    "tree_files": None,
                    "tree_directories": None,
                },
            ],
            "minimum_free_bytes": 1,
            "free_space_gate_passed": True,
            "post_execution_stat_unchanged": True,
        },
        "declared_steps": len(steps),
        "steps": steps,
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    report["execution_identity_sha256"] = canonical_json_sha256(report)
    path = tmp_path / filename
    write_json(path, report)
    return path


def test_file_receipt_is_deterministic_portable_and_exact(tmp_path: Path) -> None:
    prepared, gate, resources = boundary(tmp_path)
    repeated = prepare_caller_run_boundary(
        gate,
        resources,
        caller_name="Manta",
        caller_version="1.6.0",
        gate_caller="manta",
        expected_resource_manifest_sha256=sha256_file(resources),
        tool=tool(),
        command_semantics=command(),
        forbidden_bindings=(PRIVATE_SUBJECT,),
    )
    assert repeated == prepared
    raw = tmp_path / "diploidSV.vcf.gz"
    raw.write_bytes(b"synthetic raw VCF control\n")
    execution = write_execution_report(tmp_path, prepared)
    output = tmp_path / "caller-receipt.json"
    receipt = write_caller_run_receipt(
        output,
        prepared,
        gate,
        resources,
        execution_report=execution,
        raw_artifact_id="manta-diploid-sv-vcf",
        identity_scheme="file_sha256",
        component_paths={"primary_vcf": raw},
        forbidden_bindings=(PRIVATE_SUBJECT,),
    )
    assert receipt["status"] == "PASS"
    assert (
        receipt["run_boundary"]["caller_run_identity_sha256"]
        == prepared["caller_run_identity_sha256"]
    )
    assert receipt["raw_artifact"]["raw_artifact_sha256"] == sha256_file(raw)
    loaded, file_digest = load_caller_run_receipt(output, forbidden_bindings=(PRIVATE_SUBJECT,))
    assert loaded == receipt
    assert file_digest == sha256_file(output)
    payload = output.read_text(encoding="utf-8")
    assert PRIVATE_SUBJECT not in payload
    assert str(gate) not in payload
    assert str(resources) not in payload
    assert str(raw) not in payload
    assert str(output) not in payload


def test_receipt_rejects_raw_artifact_through_symbolic_parent(tmp_path: Path) -> None:
    prepared, gate, resources = boundary(tmp_path)
    real_directory = tmp_path / "real-output"
    real_directory.mkdir()
    raw = real_directory / "diploidSV.vcf.gz"
    raw.write_bytes(b"synthetic raw VCF control\n")
    alias = tmp_path / "aliased-output"
    alias.symlink_to(real_directory, target_is_directory=True)
    execution = write_execution_report(tmp_path, prepared)

    with pytest.raises(PortableEvidenceError, match="symbolic path component"):
        write_caller_run_receipt(
            tmp_path / "caller-receipt.json",
            prepared,
            gate,
            resources,
            execution_report=execution,
            raw_artifact_id="manta-diploid-sv-vcf",
            identity_scheme="file_sha256",
            component_paths={"primary_vcf": alias / raw.name},
        )


def test_prepared_boundary_is_restart_safe_atomic_and_tamper_evident(tmp_path: Path) -> None:
    gate = gate_report(tmp_path)
    resources = resource_manifest(tmp_path)
    output = tmp_path / "prepared-boundary.json"
    checkpoint = write_prepared_caller_run_boundary(
        output,
        gate,
        resources,
        caller_name="Manta",
        caller_version="1.6.0",
        gate_caller="manta",
        expected_resource_manifest_sha256=sha256_file(resources),
        tool=tool(),
        command_semantics=command(),
        forbidden_bindings=(PRIVATE_SUBJECT,),
    )
    assert checkpoint["status"] == "PREPARED"
    loaded, checkpoint_file_sha256 = load_prepared_caller_run_boundary(
        output, forbidden_bindings=(PRIVATE_SUBJECT,)
    )
    assert loaded == checkpoint
    assert checkpoint_file_sha256 == sha256_file(output)
    payload = output.read_text(encoding="utf-8")
    assert PRIVATE_SUBJECT not in payload
    assert str(gate) not in payload
    assert str(resources) not in payload
    original = output.read_bytes()
    with pytest.raises(PortableEvidenceError, match="refusing to replace"):
        write_prepared_caller_run_boundary(
            output,
            gate,
            resources,
            caller_name="Manta",
            caller_version="1.6.0",
            gate_caller="manta",
            expected_resource_manifest_sha256=sha256_file(resources),
            tool=tool(),
            command_semantics=command(),
        )
    assert output.read_bytes() == original

    tampered = json.loads(output.read_text(encoding="utf-8"))
    tampered["run_boundary"]["command_semantics"]["parameters"]["threads"] = 2
    boundary_without_identity = {
        key: value
        for key, value in tampered["run_boundary"].items()
        if key != "caller_run_identity_sha256"
    }
    tampered["run_boundary"]["caller_run_identity_sha256"] = canonical_json_sha256(
        boundary_without_identity
    )
    with pytest.raises(PortableEvidenceError, match="checkpoint identity changed"):
        validate_prepared_caller_run_boundary_payload(tampered)


def test_prepare_finalize_validate_cli_is_path_free_and_sanitized(tmp_path: Path) -> None:
    gate = gate_report(tmp_path)
    resources = resource_manifest(tmp_path)
    tool_json = tmp_path / "tool.json"
    command_json = tmp_path / "command.json"
    write_json(tool_json, tool())
    write_json(command_json, command())
    forbidden = tmp_path / "private-bindings.txt"
    forbidden.write_text(f"# private runtime values\n{PRIVATE_SUBJECT}\n", encoding="utf-8")
    boundary_output = tmp_path / "prepared.json"
    prepared = run_cli(
        [
            "prepare",
            "--input-gate",
            str(gate),
            "--resource-manifest",
            str(resources),
            "--expected-resource-manifest-sha256",
            sha256_file(resources),
            "--caller-name",
            "Manta",
            "--caller-version",
            "1.6.0",
            "--gate-caller",
            "manta",
            "--tool-json",
            str(tool_json),
            "--command-json",
            str(command_json),
            "--forbidden-bindings-file",
            str(forbidden),
            "--output",
            str(boundary_output),
        ]
    )
    assert prepared.returncode == 0, prepared.stderr
    assert json.loads(prepared.stdout)["status"] == "PREPARED"
    assert boundary_output.is_file()
    boundary_validated = run_cli(
        [
            "validate-boundary",
            "--boundary",
            str(boundary_output),
            "--forbidden-bindings-file",
            str(forbidden),
        ]
    )
    assert boundary_validated.returncode == 0, boundary_validated.stderr
    assert json.loads(boundary_validated.stdout) == json.loads(prepared.stdout)
    checkpoint = json.loads(boundary_output.read_text(encoding="utf-8"))
    execution = write_execution_report(tmp_path, checkpoint["run_boundary"])

    raw = tmp_path / "raw.vcf"
    raw.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    artifact_spec = tmp_path / "artifact.json"
    write_json(
        artifact_spec,
        {
            "schema_version": 1,
            "raw_artifact_id": "primary-vcf",
            "identity_scheme": "file_sha256",
            "components": {"primary_vcf": str(raw)},
            "derived_payload": None,
        },
    )
    receipt_output = tmp_path / "receipt.json"
    finalized = run_cli(
        [
            "finalize",
            "--boundary",
            str(boundary_output),
            "--input-gate",
            str(gate),
            "--resource-manifest",
            str(resources),
            "--execution-report",
            str(execution),
            "--artifact-spec",
            str(artifact_spec),
            "--forbidden-bindings-file",
            str(forbidden),
            "--output",
            str(receipt_output),
        ]
    )
    assert finalized.returncode == 0, finalized.stderr
    assert json.loads(finalized.stdout)["status"] == "PASS"
    validated = run_cli(
        [
            "validate",
            "--receipt",
            str(receipt_output),
            "--forbidden-bindings-file",
            str(forbidden),
        ]
    )
    assert validated.returncode == 0, validated.stderr
    assert json.loads(validated.stdout) == json.loads(finalized.stdout)

    portable_payload = (
        boundary_output.read_text(encoding="utf-8")
        + receipt_output.read_text(encoding="utf-8")
        + prepared.stdout
        + boundary_validated.stdout
        + finalized.stdout
        + validated.stdout
    )
    forbidden_values = {
        PRIVATE_SUBJECT,
        str(gate),
        str(resources),
        str(tool_json),
        str(command_json),
        str(forbidden),
        str(boundary_output),
        str(execution),
        str(raw),
        str(artifact_spec),
        str(receipt_output),
    }
    assert not any(value in portable_payload for value in forbidden_values)

    repeated = run_cli(
        [
            "prepare",
            "--input-gate",
            str(gate),
            "--resource-manifest",
            str(resources),
            "--expected-resource-manifest-sha256",
            sha256_file(resources),
            "--caller-name",
            "Manta",
            "--caller-version",
            "1.6.0",
            "--gate-caller",
            "manta",
            "--tool-json",
            str(tool_json),
            "--command-json",
            str(command_json),
            "--forbidden-bindings-file",
            str(forbidden),
            "--output",
            str(boundary_output),
        ]
    )
    assert repeated.returncode == 2
    assert repeated.stdout == ""
    assert repeated.stderr == (
        "caller receipt operation failed; inspect the authorized private process log\n"
    )
    assert not any(value in repeated.stderr for value in forbidden_values)


def test_named_artifact_set_matches_adapter_hash_convention(tmp_path: Path) -> None:
    prepared, gate, resources = boundary(tmp_path)
    execution = write_execution_report(tmp_path, prepared)
    json_output = tmp_path / "expansion.json"
    vcf_output = tmp_path / "expansion.vcf"
    json_output.write_text('{"LocusResults": {}}\n', encoding="utf-8")
    vcf_output.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    receipt = write_caller_run_receipt(
        tmp_path / "artifact-set-receipt.json",
        prepared,
        gate,
        resources,
        execution_report=execution,
        raw_artifact_id="expansionhunter-json-vcf",
        identity_scheme="named_sha256_set",
        component_paths={
            "expansionhunter_vcf_sha256": vcf_output,
            "expansionhunter_json_sha256": json_output,
        },
    )
    expected = canonical_json_sha256(
        {
            "expansionhunter_json_sha256": sha256_file(json_output),
            "expansionhunter_vcf_sha256": sha256_file(vcf_output),
        }
    )
    assert receipt["raw_artifact"]["raw_artifact_sha256"] == expected
    assert [component["role"] for component in receipt["raw_artifact"]["components"]] == [
        "expansionhunter_json_sha256",
        "expansionhunter_vcf_sha256",
    ]


def test_derived_evidence_receipt_binds_full_canonical_payload(tmp_path: Path) -> None:
    prepared, gate, resources = boundary(tmp_path)
    execution = write_execution_report(tmp_path, prepared)
    derived = {
        "target_manifest_sha256": "6" * 64,
        "policy": {"minimum_mapping_quality": 20, "minimum_base_quality": 20},
        "record_digests": ["7" * 64, "8" * 64],
    }
    receipt = write_caller_run_receipt(
        tmp_path / "derived-receipt.json",
        prepared,
        gate,
        resources,
        execution_report=execution,
        raw_artifact_id="bounded-readback-evidence",
        identity_scheme="derived_canonical_json",
        derived_payload=derived,
    )
    assert receipt["raw_artifact"]["components"] == []
    assert receipt["raw_artifact"]["derived_payload"] == derived
    assert receipt["raw_artifact"]["raw_artifact_sha256"] == canonical_json_sha256(derived)


def test_boundary_rejects_changed_resource_and_nonisolated_tool(tmp_path: Path) -> None:
    gate = gate_report(tmp_path)
    resources = resource_manifest(tmp_path)
    with pytest.raises(PortableEvidenceError, match="resource manifest digest changed"):
        prepare_caller_run_boundary(
            gate,
            resources,
            caller_name="Manta",
            caller_version="1.6.0",
            gate_caller="manta",
            expected_resource_manifest_sha256="0" * 64,
            tool=tool(),
            command_semantics=command(),
        )
    unsafe_tool = {**tool(), "network_disabled": False}
    with pytest.raises(PortableEvidenceError, match="network isolation"):
        prepare_caller_run_boundary(
            gate,
            resources,
            caller_name="Manta",
            caller_version="1.6.0",
            gate_caller="manta",
            expected_resource_manifest_sha256=sha256_file(resources),
            tool=unsafe_tool,
            command_semantics=command(),
        )
    native_tool = {**tool(), "execution_kind": "binary"}
    with pytest.raises(PortableEvidenceError, match="execution kind is unsupported"):
        prepare_caller_run_boundary(
            gate,
            resources,
            caller_name="Manta",
            caller_version="1.6.0",
            gate_caller="manta",
            expected_resource_manifest_sha256=sha256_file(resources),
            tool=native_tool,
            command_semantics=command(),
        )


def test_receipt_rehashes_gate_and_resource_after_execution(tmp_path: Path) -> None:
    prepared, gate, resources = boundary(tmp_path)
    execution = write_execution_report(tmp_path, prepared)
    raw = tmp_path / "raw.vcf"
    raw.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    resources.write_text("changed after caller preparation\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="resource manifest changed"):
        write_caller_run_receipt(
            tmp_path / "changed-resource-receipt.json",
            prepared,
            gate,
            resources,
            execution_report=execution,
            raw_artifact_id="primary-vcf",
            identity_scheme="file_sha256",
            component_paths={"primary_vcf": raw},
        )

    prepared, gate, resources = boundary(tmp_path / "second")
    execution = write_execution_report(tmp_path / "second", prepared)
    raw = tmp_path / "second-raw.vcf"
    raw.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    report = json.loads(gate.read_text(encoding="utf-8"))
    report["technical_note"] = "changed"
    write_json(gate, report)
    with pytest.raises(PortableEvidenceError, match="input-gate report changed"):
        write_caller_run_receipt(
            tmp_path / "changed-gate-receipt.json",
            prepared,
            gate,
            resources,
            execution_report=execution,
            raw_artifact_id="primary-vcf",
            identity_scheme="file_sha256",
            component_paths={"primary_vcf": raw},
        )


def test_receipt_rejects_failed_or_different_oci_execution(tmp_path: Path) -> None:
    prepared, gate, resources = boundary(tmp_path)
    raw = tmp_path / "raw.vcf"
    raw.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    failed_execution = write_execution_report(tmp_path, prepared, filename="failed-execution.json")
    failed = json.loads(failed_execution.read_text(encoding="utf-8"))
    failed["status"] = "FAIL_RETAINED"
    failed["error_code"] = "caller_step_failed"
    failed["input_verification"]["post_execution_stat_unchanged"] = False
    failed["execution_identity_sha256"] = canonical_json_sha256(
        {key: value for key, value in failed.items() if key != "execution_identity_sha256"}
    )
    write_json(failed_execution, failed)
    with pytest.raises(PortableEvidenceError, match="not canonical PASS"):
        write_caller_run_receipt(
            tmp_path / "failed-execution-receipt.json",
            prepared,
            gate,
            resources,
            execution_report=failed_execution,
            raw_artifact_id="primary-vcf",
            identity_scheme="file_sha256",
            component_paths={"primary_vcf": raw},
        )

    different_execution = write_execution_report(
        tmp_path, prepared, filename="different-execution.json"
    )
    different = json.loads(different_execution.read_text(encoding="utf-8"))
    different["caller_run_identity_sha256"] = "0" * 64
    different["execution_identity_sha256"] = canonical_json_sha256(
        {key: value for key, value in different.items() if key != "execution_identity_sha256"}
    )
    write_json(different_execution, different)
    with pytest.raises(PortableEvidenceError, match="differs from its run boundary"):
        write_caller_run_receipt(
            tmp_path / "different-execution-receipt.json",
            prepared,
            gate,
            resources,
            execution_report=different_execution,
            raw_artifact_id="primary-vcf",
            identity_scheme="file_sha256",
            component_paths={"primary_vcf": raw},
        )


def test_receipt_rejects_paths_private_tokens_and_overwrite(tmp_path: Path) -> None:
    gate = gate_report(tmp_path)
    resources = resource_manifest(tmp_path)
    leaking_command = deepcopy(command())
    leaking_command["steps"][0]["arguments"].append("/tmp/private-input.bam")
    with pytest.raises(PortableEvidenceError, match="local path"):
        prepare_caller_run_boundary(
            gate,
            resources,
            caller_name="Manta",
            caller_version="1.6.0",
            gate_caller="manta",
            expected_resource_manifest_sha256=sha256_file(resources),
            tool=tool(),
            command_semantics=leaking_command,
        )
    unsafe_step_name = deepcopy(command())
    unsafe_step_name["steps"][0]["name"] = "nested/step"
    with pytest.raises(PortableEvidenceError, match="step name is invalid"):
        prepare_caller_run_boundary(
            gate,
            resources,
            caller_name="Manta",
            caller_version="1.6.0",
            gate_caller="manta",
            expected_resource_manifest_sha256=sha256_file(resources),
            tool=tool(),
            command_semantics=unsafe_step_name,
        )

    prepared, gate, resources = boundary(tmp_path / "portable")
    execution = write_execution_report(tmp_path / "portable", prepared)
    raw = tmp_path / "portable-raw.vcf"
    raw.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    output = tmp_path / "portable-receipt.json"
    write_caller_run_receipt(
        output,
        prepared,
        gate,
        resources,
        execution_report=execution,
        raw_artifact_id="primary-vcf",
        identity_scheme="file_sha256",
        component_paths={"primary_vcf": raw},
    )
    original = output.read_bytes()
    with pytest.raises(PortableEvidenceError, match="refusing to replace"):
        write_caller_run_receipt(
            output,
            prepared,
            gate,
            resources,
            execution_report=execution,
            raw_artifact_id="primary-vcf",
            identity_scheme="file_sha256",
            component_paths={"primary_vcf": raw},
        )
    assert output.read_bytes() == original

    receipt = json.loads(output.read_text(encoding="utf-8"))
    receipt["raw_artifact"]["raw_artifact_id"] = "artifact-A"
    without_identity = {
        key: value for key, value in receipt.items() if key != "receipt_identity_sha256"
    }
    receipt["receipt_identity_sha256"] = canonical_json_sha256(without_identity)
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        validate_caller_run_receipt_payload(receipt, forbidden_bindings=("A",))


def test_independent_validator_detects_receipt_tampering(tmp_path: Path) -> None:
    prepared, gate, resources = boundary(tmp_path)
    execution = write_execution_report(tmp_path, prepared)
    raw = tmp_path / "raw.vcf"
    raw.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    output = tmp_path / "receipt.json"
    write_caller_run_receipt(
        output,
        prepared,
        gate,
        resources,
        execution_report=execution,
        raw_artifact_id="primary-vcf",
        identity_scheme="file_sha256",
        component_paths={"primary_vcf": raw},
    )
    receipt = json.loads(output.read_text(encoding="utf-8"))
    receipt["raw_artifact"]["components"][0]["bytes"] += 1
    write_json(output, receipt)
    with pytest.raises(PortableEvidenceError, match="receipt identity changed"):
        load_caller_run_receipt(output)
