from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path

import pytest

from genome_annotation_builder.release_c.caller_profile import load_caller_profile
from genome_annotation_builder.release_c.caller_receipt import (
    canonical_json_sha256,
    prepare_caller_run_boundary,
    write_caller_run_receipt,
)
from genome_annotation_builder.release_c.extension_stage import (
    MANIFEST_FILENAME,
    PortableEvidenceError,
    stage_release_c_extensions,
    subject_binding_from_private_manifest,
    validate_release_c_extension_stage,
)
from genome_annotation_builder.release_c.portable_model import (
    canonical_row_multiset,
    make_evidence_id,
    write_extension,
)
from genome_annotation_builder.release_c.profile_catalog import resolve_shipped_profile

from .test_caller_receipt import write_execution_report
from .test_portable_model import rows as example_rows

REFERENCE = "1" * 64
PRIVATE_SUBJECT = "SYNTHETIC_LOGICAL_SUBJECT"
PRIVATE_LAB = "SYNTHETIC_LAB"
ALIGNMENT_A = "a" * 64
ALIGNMENT_B = "b" * 64


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def refresh_receipt_identity(receipt: dict[str, object], *, boundary: bool = False) -> None:
    if boundary:
        run_boundary = receipt["run_boundary"]
        assert isinstance(run_boundary, dict)
        run_boundary["caller_run_identity_sha256"] = canonical_json_sha256(
            {
                key: value
                for key, value in run_boundary.items()
                if key != "caller_run_identity_sha256"
            }
        )
        execution = receipt["execution"]
        assert isinstance(execution, dict)
        report = execution["report"]
        assert isinstance(report, dict)
        report.update(
            {
                "caller_run_identity_sha256": run_boundary["caller_run_identity_sha256"],
                "source_alignment_sha256": run_boundary["source_alignment_sha256"],
                "reference_sha256": run_boundary["reference_sha256"],
                "resource_manifest_sha256": run_boundary["resource_manifest_sha256"],
                "tool": run_boundary["tool"],
            }
        )
        for step, semantic in zip(
            report["steps"], run_boundary["command_semantics"]["steps"], strict=True
        ):
            step["name"] = semantic["name"]
            step["semantic_step_sha256"] = canonical_json_sha256(semantic)
        report["execution_identity_sha256"] = canonical_json_sha256(
            {key: value for key, value in report.items() if key != "execution_identity_sha256"}
        )
        execution["execution_report_file_sha256"] = hashlib.sha256(
            (json.dumps(report, indent=2, sort_keys=True) + "\n").encode()
        ).hexdigest()
    receipt["receipt_identity_sha256"] = canonical_json_sha256(
        {key: value for key, value in receipt.items() if key != "receipt_identity_sha256"}
    )


def private_manifest(tmp_path: Path, *, mixed_subjects: bool = False) -> Path:
    manifest = {
        "schema_version": 1,
        "status": "PASS",
        "reference": {
            "assembly": "GRCh38 synthetic control",
            "fasta_path": str(tmp_path / "private-reference.fa"),
            "fasta_sha256": REFERENCE,
        },
        "datasets": [
            {
                "dataset_id": "synthetic-source-a",
                "logical_subject": PRIVATE_SUBJECT,
                "laboratory": PRIVATE_LAB,
                "original_markdup_bam": {
                    "path": str(tmp_path / "private-a.bam"),
                    "index_path": str(tmp_path / "private-a.bam.bai"),
                    "sha256": ALIGNMENT_A,
                },
                "baf_vcf": {"sample": "SYNTHETIC_TECHNICAL_A"},
            },
            {
                "dataset_id": "synthetic-source-b",
                "logical_subject": "ANOTHER_SYNTHETIC_SUBJECT"
                if mixed_subjects
                else PRIVATE_SUBJECT,
                "laboratory": PRIVATE_LAB,
                "original_markdup_bam": {
                    "path": str(tmp_path / "private-b.bam"),
                    "index_path": str(tmp_path / "private-b.bam.bai"),
                    "sha256": ALIGNMENT_B,
                },
                "baf_vcf": {"sample": "SYNTHETIC_TECHNICAL_B"},
            },
        ],
    }
    path = tmp_path / "private-input-manifest.json"
    write_json(path, manifest)
    return path


def gate_report(tmp_path: Path, caller: str, alignment: str) -> Path:
    report = {
        "schema_version": 1,
        "status": "PASS",
        "caller": caller,
        "selected_alignment": {"sha256": alignment},
        "reference": {"sha256": REFERENCE},
        "run_identity_sha256": hashlib.sha256(f"{caller}:{alignment}".encode()).hexdigest(),
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    path = tmp_path / f"gate-{caller}-{alignment[0]}.json"
    write_json(path, report)
    return path


def receipt_command(name: str) -> dict[str, object]:
    return {
        "steps": [
            {
                "name": "call",
                "entrypoint": name,
                "arguments": ["--input", "{selected_alignment}", "--output", "{raw_artifact}"],
            }
        ],
        "environment": {"LC_ALL": "C"},
        "parameters": {"synthetic_control": True},
    }


def caller_receipts(
    tmp_path: Path, input_gates: list[Path]
) -> tuple[list[Path], dict[str, dict[str, str]]]:
    specifications = [
        ("Manta", "1.6.0", "manta", "8", "synthetic-manta-vcf"),
        (
            "GATK Mutect2 mitochondrial mode",
            "4.6.2.0",
            "gatk_mitochondria",
            "9",
            "synthetic-mitochondrial-vcf",
        ),
    ]
    gate_by_caller = {
        str(json.loads(path.read_text(encoding="utf-8"))["caller"]): path for path in input_gates
    }
    receipt_paths: list[Path] = []
    provenance: dict[str, dict[str, str]] = {}
    for caller_name, caller_version, gate_caller, identity_digit, artifact_id in specifications:
        gate = gate_by_caller[gate_caller]
        resources = tmp_path / f"resources-{gate_caller}.json"
        write_json(
            resources,
            {
                "schema_version": 1,
                "status": "PASS",
                "resource_identity_sha256": identity_digit * 64,
            },
        )
        raw = tmp_path / f"raw-{gate_caller}.vcf"
        raw.write_text(f"synthetic {gate_caller} raw artifact\n", encoding="utf-8")
        boundary = prepare_caller_run_boundary(
            gate,
            resources,
            caller_name=caller_name,
            caller_version=caller_version,
            gate_caller=gate_caller,
            expected_resource_manifest_sha256=sha256_file(resources),
            tool={
                "execution_kind": "oci_image",
                "logical_name": gate_caller,
                "identity_sha256": identity_digit * 64,
                "network_disabled": True,
            },
            command_semantics=receipt_command(gate_caller),
            forbidden_bindings=(PRIVATE_SUBJECT, PRIVATE_LAB),
        )
        execution = write_execution_report(
            tmp_path, boundary, filename=f"execution-{gate_caller}.json"
        )
        receipt_path = tmp_path / f"receipt-{gate_caller}.json"
        receipt = write_caller_run_receipt(
            receipt_path,
            boundary,
            gate,
            resources,
            execution_report=execution,
            raw_artifact_id=artifact_id,
            identity_scheme="file_sha256",
            component_paths={"primary_output": raw},
            forbidden_bindings=(PRIVATE_SUBJECT, PRIVATE_LAB),
        )
        receipt_paths.append(receipt_path)
        provenance[caller_name] = {
            "caller_run_identity_sha256": str(boundary["caller_run_identity_sha256"]),
            "resource_manifest_sha256": str(boundary["resource_manifest_sha256"]),
            "raw_artifact_id": artifact_id,
            "raw_artifact_sha256": str(receipt["raw_artifact"]["raw_artifact_sha256"]),
        }
    return receipt_paths, provenance


def extension_inputs(
    tmp_path: Path, receipt_provenance: dict[str, dict[str, str]]
) -> dict[str, Path]:
    structural = deepcopy(example_rows()["structural_variants"])
    structural.update(
        {
            "caller_name": "Manta",
            "caller_version": "1.6.0",
            "source_alignment_sha256": ALIGNMENT_A,
            **receipt_provenance["Manta"],
            "source_record_digest": "4" * 64,
        }
    )
    structural["evidence_id"] = make_evidence_id("structural_variants", structural)
    mitochondrial = deepcopy(example_rows()["mitochondrial_heteroplasmy"])
    mitochondrial.update(
        {
            "caller_name": "GATK Mutect2 mitochondrial mode",
            "caller_version": "4.6.2.0",
            "source_alignment_sha256": ALIGNMENT_B,
            **receipt_provenance["GATK Mutect2 mitochondrial mode"],
            "source_record_digest": "7" * 64,
        }
    )
    mitochondrial["evidence_id"] = make_evidence_id("mitochondrial_heteroplasmy", mitochondrial)
    paths = {
        "structural_variants": tmp_path / "input-structural.parquet",
        "mitochondrial_heteroplasmy": tmp_path / "input-mitochondrial.parquet",
    }
    write_extension(
        "structural_variants",
        [structural],
        paths["structural_variants"],
        expected_reference_sha256=REFERENCE,
    )
    write_extension(
        "mitochondrial_heteroplasmy",
        [mitochondrial],
        paths["mitochondrial_heteroplasmy"],
        expected_reference_sha256=REFERENCE,
    )
    return paths


def write_adapter_summary(
    path: Path,
    receipt_path: Path,
    *,
    profile_id: str,
    technical_source_role: str,
    domain_files: dict[str, Path],
    manta_omitted_non_bnd_records: int = 0,
) -> Path:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    boundary = receipt["run_boundary"]
    artifact = receipt["raw_artifact"]
    profile = load_caller_profile(resolve_shipped_profile(profile_id))
    files: list[dict[str, object]] = []
    parsing: list[dict[str, object]] = []
    import pyarrow.parquet as pq

    for domain, source in sorted(domain_files.items()):
        rows = pq.read_table(source).to_pylist()
        multiset = canonical_row_multiset(rows)
        files.append(
            {
                "domain": domain,
                "path": source.name,
                "sha256": sha256_file(source),
                "bytes": source.stat().st_size,
                "rows": len(rows),
                "canonical_row_multiset": multiset,
            }
        )
        parse_summary: dict[str, object] = {
            "schema_version": 1,
            "status": "PASS",
            "source_records": len(rows),
            "typed_rows": len(rows),
            "record_reconciliation_delta": 0,
            "raw_artifact_id": artifact["raw_artifact_id"],
            "raw_artifact_sha256": artifact["raw_artifact_sha256"],
            "canonical_row_multiset": multiset,
            "sample_identifier_persisted": False,
            "private_paths_persisted": False,
        }
        if profile_id == "manta-germline-1.6.0":
            omitted = manta_omitted_non_bnd_records
            parse_summary.update(
                {
                    "source_records": len(rows) + omitted,
                    "record_reconciliation_delta": omitted,
                    "portable_contig_scope": "grch38_core_chr1_22_x_y_m",
                    "omitted_noncanonical_records": omitted,
                    "omitted_noncanonical_non_bnd_records": omitted,
                    "omitted_noncanonical_bnd_records": 0,
                    "omitted_noncanonical_bnd_pairs": 0,
                    "core_to_noncanonical_bnd_pairs": 0,
                    "noncanonical_to_noncanonical_bnd_pairs": 0,
                    "source_bnd_records": 0,
                    "source_reciprocal_bnd_pairs": 0,
                    "bnd_records": 0,
                    "reciprocal_bnd_pairs": 0,
                }
            )
        parsing.append({"domain": domain, "summary": parse_summary})
    summary: dict[str, object] = {
        "schema_version": 1,
        "status": "PASS",
        "profile_id": profile_id,
        "profile_identity_sha256": profile.identity_sha256,
        "adapter_options_identity_sha256": "c" * 64,
        "technical_source_role": technical_source_role,
        "reference_sha256": boundary["reference_sha256"],
        "source_alignment_sha256": boundary["source_alignment_sha256"],
        "resource_manifest_sha256": boundary["resource_manifest_sha256"],
        "caller_run_identity_sha256": boundary["caller_run_identity_sha256"],
        "caller_receipt_identity_sha256": receipt["receipt_identity_sha256"],
        "caller_receipt_file_sha256": sha256_file(receipt_path),
        "lifecycle_summary_sha256": "d" * 64,
        "raw_artifact_id": artifact["raw_artifact_id"],
        "raw_artifact_sha256": artifact["raw_artifact_sha256"],
        "raw_component_roles": sorted(row["role"] for row in artifact["components"]),
        "files": files,
        "parsing": parsing,
        "extension_required": bool(files),
        "source_rows_removed": 0,
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "raw_artifacts_embedded": False,
    }
    summary["adapter_identity_sha256"] = canonical_json_sha256(summary)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, summary)
    return path


def adapter_summaries(
    tmp_path: Path,
    receipts: list[Path],
    inputs: dict[str, Path],
) -> list[Path]:
    by_caller = {
        json.loads(path.read_text(encoding="utf-8"))["run_boundary"]["caller_name"]: path
        for path in receipts
    }
    return [
        write_adapter_summary(
            tmp_path / "adapter-manta" / "adapter-summary.json",
            by_caller["Manta"],
            profile_id="manta-germline-1.6.0",
            technical_source_role="primary",
            domain_files={"structural_variants": inputs["structural_variants"]},
            manta_omitted_non_bnd_records=1,
        ),
        write_adapter_summary(
            tmp_path / "adapter-mt" / "adapter-summary.json",
            by_caller["GATK Mutect2 mitochondrial mode"],
            profile_id="gatk-mitochondria-normal-shifted-4.6.2.0",
            technical_source_role="primary",
            domain_files={"mitochondrial_heteroplasmy": inputs["mitochondrial_heteroplasmy"]},
        ),
    ]


def binding(tmp_path: Path):
    return subject_binding_from_private_manifest(
        private_manifest(tmp_path), ["synthetic-source-a", "synthetic-source-b"]
    )


def gates(tmp_path: Path) -> list[Path]:
    return [
        gate_report(tmp_path, "manta", ALIGNMENT_A),
        gate_report(tmp_path, "gatk_mitochondria", ALIGNMENT_B),
    ]


def complete_stage_inputs(tmp_path: Path):
    subject_binding = binding(tmp_path)
    input_gates = gates(tmp_path)
    receipts, provenance = caller_receipts(tmp_path, input_gates)
    inputs = extension_inputs(tmp_path, provenance)
    adapters = adapter_summaries(tmp_path, receipts, inputs)
    return subject_binding, inputs, input_gates, receipts, adapters


def test_atomic_extension_stage_is_subject_bound_and_portable(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    output = tmp_path / "release-c-extensions"
    result = stage_release_c_extensions(
        output, inputs, input_gates, receipts, adapters, subject_binding
    )
    assert result["status"] == "PASS"
    assert result["validation"]["domains"] == 2
    assert result["validation"]["rows"] == 2
    assert result["validation"]["input_gates"] == 2
    assert result["validation"]["caller_runs"] == 2
    assert result["validation"]["adapter_runs"] == 2
    assert result["validation"]["adapter_reconciliation"]["portable_core_scope"] == {
        "policy": "grch38_core_chr1_22_x_y_m",
        "adapter_runs": 1,
        "source_records": 2,
        "typed_rows": 1,
        "omitted_records": 1,
        "omitted_non_bnd_records": 1,
        "omitted_bnd_records": 0,
        "omitted_bnd_pairs": 0,
        "core_to_noncanonical_bnd_pairs": 0,
        "noncanonical_to_noncanonical_bnd_pairs": 0,
    }
    assert result["validation"]["concordance"] == {
        "mitochondrial_heteroplasmy": {"not_assessed": 1},
        "structural_variants": {"not_assessed": 1},
    }
    assert result["validation"]["technical_source_roles"] == {
        "mitochondrial_heteroplasmy": {"primary": 1},
        "structural_variants": {"primary": 1},
    }
    assert result["manifest"]["logical_subject_count"] == 1
    assert result["manifest"]["source_dataset_count"] == 2
    assert result["manifest"]["subject_identifier_persisted"] is False
    assert result["manifest"]["private_manifest_digest_persisted"] is False
    assert result["manifest"]["source_rows_removed"] == 0
    assert result["manifest"]["empty_optional_files"] == 0
    for domain, source in inputs.items():
        copied = output / f"{domain}.parquet"
        assert copied.read_bytes() == source.read_bytes()
    payload = (output / MANIFEST_FILENAME).read_text(encoding="utf-8")
    forbidden = {
        PRIVATE_SUBJECT,
        PRIVATE_LAB,
        "synthetic-source-a",
        "synthetic-source-b",
        "SYNTHETIC_TECHNICAL_A",
        "SYNTHETIC_TECHNICAL_B",
        str(output),
        *(str(path) for path in inputs.values()),
        *(str(path) for path in input_gates),
        *(str(path) for path in receipts),
        *(str(path) for path in adapters),
    }
    forbidden.update(
        hashlib.sha256(value.encode()).hexdigest()
        for value in (
            PRIVATE_SUBJECT,
            PRIVATE_LAB,
            "synthetic-source-a",
            "synthetic-source-b",
        )
    )
    assert not any(value in payload for value in forbidden)


def test_private_binding_rejects_mixed_subjects(tmp_path: Path) -> None:
    with pytest.raises(PortableEvidenceError, match="one logical subject"):
        subject_binding_from_private_manifest(
            private_manifest(tmp_path, mixed_subjects=True),
            ["synthetic-source-a", "synthetic-source-b"],
        )


def test_private_binding_rejects_duplicate_technical_alignment(tmp_path: Path) -> None:
    manifest_path = private_manifest(tmp_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["datasets"][1]["original_markdup_bam"]["sha256"] = ALIGNMENT_A
    write_json(manifest_path, manifest)
    with pytest.raises(PortableEvidenceError, match="duplicates a technical alignment"):
        subject_binding_from_private_manifest(
            manifest_path, ["synthetic-source-a", "synthetic-source-b"]
        )


def test_missing_or_unused_gate_fails_before_publication(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    output = tmp_path / "missing-gate-output"
    with pytest.raises(PortableEvidenceError, match="lacks its exact"):
        stage_release_c_extensions(
            output, inputs, input_gates[:1], receipts, adapters, subject_binding
        )
    assert not output.exists()

    extra_gate = gate_report(tmp_path, "cnvpytor", ALIGNMENT_A)
    with pytest.raises(PortableEvidenceError, match=r"input gates.*strict-adapter"):
        stage_release_c_extensions(
            output,
            inputs,
            [*input_gates, extra_gate],
            receipts,
            adapters,
            subject_binding,
        )
    assert not output.exists()


def test_gate_privacy_claim_must_be_explicitly_pass(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    report = json.loads(input_gates[0].read_text(encoding="utf-8"))
    report["subject_identifier_persisted"] = True
    write_json(input_gates[0], report)
    output = tmp_path / "privacy-invalid-gate-output"
    with pytest.raises(PortableEvidenceError, match="privacy/subject boundary"):
        stage_release_c_extensions(output, inputs, input_gates, receipts, adapters, subject_binding)
    assert not output.exists()


def test_zero_row_adapter_run_remains_explicitly_receipted(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    mitochondrial_adapter = next(path for path in adapters if path.parent.name == "adapter-mt")
    payload = json.loads(mitochondrial_adapter.read_text(encoding="utf-8"))
    empty_multiset = canonical_row_multiset([])
    payload["files"] = []
    payload["extension_required"] = False
    parsing = payload["parsing"][0]["summary"]
    parsing["source_records"] = 0
    parsing["typed_rows"] = 0
    parsing["record_reconciliation_delta"] = 0
    parsing["canonical_row_multiset"] = empty_multiset
    payload["adapter_identity_sha256"] = canonical_json_sha256(
        {key: value for key, value in payload.items() if key != "adapter_identity_sha256"}
    )
    write_json(mitochondrial_adapter, payload)
    nonempty_inputs = {"structural_variants": inputs["structural_variants"]}
    output = tmp_path / "zero-row-adapter-stage"
    result = stage_release_c_extensions(
        output,
        nonempty_inputs,
        input_gates,
        receipts,
        adapters,
        subject_binding,
    )
    assert result["validation"]["rows"] == 1
    assert result["validation"]["caller_runs"] == 2
    assert result["validation"]["adapter_runs"] == 2
    assert result["validation"]["adapter_reconciliation"]["typed_rows"] == 1


def test_multiple_files_for_one_domain_require_exact_parsing_bindings(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    import pyarrow.parquet as pq

    first_row = pq.read_table(inputs["mitochondrial_heteroplasmy"]).to_pylist()[0]
    second_row = deepcopy(first_row)
    second_row.update(
        {
            "event_id": "mt:chrM:101:C:T",
            "pos": 101,
            "ref": "C",
            "alt": "T",
            "source_record_digest": "8" * 64,
        }
    )
    second_row["evidence_id"] = make_evidence_id("mitochondrial_heteroplasmy", second_row)
    first_file = tmp_path / "mt-first.parquet"
    second_file = tmp_path / "mt-second.parquet"
    combined_file = tmp_path / "mt-combined.parquet"
    write_extension(
        "mitochondrial_heteroplasmy", [first_row], first_file, expected_reference_sha256=REFERENCE
    )
    write_extension(
        "mitochondrial_heteroplasmy",
        [second_row],
        second_file,
        expected_reference_sha256=REFERENCE,
    )
    write_extension(
        "mitochondrial_heteroplasmy",
        [first_row, second_row],
        combined_file,
        expected_reference_sha256=REFERENCE,
    )

    mitochondrial_adapter = next(path for path in adapters if path.parent.name == "adapter-mt")
    payload = json.loads(mitochondrial_adapter.read_text(encoding="utf-8"))
    template = payload["parsing"][0]["summary"]
    files: list[dict[str, object]] = []
    parsing: list[dict[str, object]] = []
    for source, rows in ((first_file, [first_row]), (second_file, [second_row])):
        multiset = canonical_row_multiset(rows)
        files.append(
            {
                "domain": "mitochondrial_heteroplasmy",
                "path": source.name,
                "sha256": sha256_file(source),
                "bytes": source.stat().st_size,
                "rows": 1,
                "canonical_row_multiset": multiset,
            }
        )
        summary = deepcopy(template)
        summary.update(
            {
                "source_records": 1,
                "typed_rows": 1,
                "record_reconciliation_delta": 0,
                "canonical_row_multiset": multiset,
            }
        )
        parsing.append({"domain": "mitochondrial_heteroplasmy", "summary": summary})
    payload["files"] = sorted(files, key=lambda row: (str(row["domain"]), str(row["path"])))
    payload["parsing"] = sorted(
        parsing,
        key=lambda row: (str(row["domain"]), canonical_json_sha256(row["summary"])),
    )
    payload["adapter_identity_sha256"] = canonical_json_sha256(
        {key: value for key, value in payload.items() if key != "adapter_identity_sha256"}
    )
    write_json(mitochondrial_adapter, payload)
    expanded_inputs = dict(inputs)
    expanded_inputs["mitochondrial_heteroplasmy"] = combined_file
    output = tmp_path / "multi-file-domain-stage"
    result = stage_release_c_extensions(
        output,
        expanded_inputs,
        input_gates,
        receipts,
        adapters,
        subject_binding,
    )
    assert result["validation"]["rows"] == 3
    assert result["validation"]["adapter_reconciliation"]["typed_rows"] == 3
    assert result["validation"]["adapter_reconciliation"]["parsing_summaries"] == 3

    invalid_payload = deepcopy(payload)
    invalid_payload["parsing"][1]["summary"]["canonical_row_multiset"] = deepcopy(
        invalid_payload["parsing"][0]["summary"]["canonical_row_multiset"]
    )
    invalid_payload["adapter_identity_sha256"] = canonical_json_sha256(
        {key: value for key, value in invalid_payload.items() if key != "adapter_identity_sha256"}
    )
    invalid_adapter = tmp_path / "invalid-multi-file-adapter.json"
    write_json(invalid_adapter, invalid_payload)
    invalid_adapters = [
        invalid_adapter if path == mitochondrial_adapter else path for path in adapters
    ]
    with pytest.raises(PortableEvidenceError, match="exact typed-file metadata"):
        stage_release_c_extensions(
            tmp_path / "invalid-multi-file-domain-stage",
            expanded_inputs,
            input_gates,
            receipts,
            invalid_adapters,
            subject_binding,
        )


def test_missing_or_unused_caller_receipt_fails_before_publication(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    output = tmp_path / "missing-receipt-output"
    with pytest.raises(PortableEvidenceError, match="lacks its exact caller-run"):
        stage_release_c_extensions(
            output, inputs, input_gates, receipts[:1], adapters, subject_binding
        )
    assert not output.exists()

    extra_payload = json.loads(receipts[0].read_text(encoding="utf-8"))
    extra_payload["raw_artifact"]["raw_artifact_id"] = "unused-raw-artifact"
    refresh_receipt_identity(extra_payload)
    extra_receipt = tmp_path / "receipt-unused.json"
    write_json(extra_receipt, extra_payload)
    with pytest.raises(PortableEvidenceError, match="caller receipts and strict-adapter"):
        stage_release_c_extensions(
            output,
            inputs,
            input_gates,
            [*receipts, extra_receipt],
            adapters,
            subject_binding,
        )
    assert not output.exists()


def test_receipt_must_match_exact_gate_and_private_boundary(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    wrong_gate_payload = json.loads(receipts[0].read_text(encoding="utf-8"))
    wrong_gate_payload["run_boundary"]["input_gate"]["report_sha256"] = "0" * 64
    refresh_receipt_identity(wrong_gate_payload, boundary=True)
    wrong_gate = tmp_path / "receipt-wrong-gate.json"
    write_json(wrong_gate, wrong_gate_payload)
    output = tmp_path / "wrong-receipt-boundary-output"
    with pytest.raises(PortableEvidenceError, match="input-gate boundary differ"):
        stage_release_c_extensions(
            output,
            inputs,
            input_gates,
            [wrong_gate, receipts[1]],
            adapters,
            subject_binding,
        )
    assert not output.exists()

    leaking_payload = json.loads(receipts[0].read_text(encoding="utf-8"))
    leaking_payload["run_boundary"]["command_semantics"]["parameters"]["note"] = (
        f"artifact-{PRIVATE_SUBJECT}"
    )
    refresh_receipt_identity(leaking_payload, boundary=True)
    leaking_receipt = tmp_path / "receipt-private-leak.json"
    write_json(leaking_receipt, leaking_payload)
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        stage_release_c_extensions(
            output,
            inputs,
            input_gates,
            [leaking_receipt, receipts[1]],
            adapters,
            subject_binding,
        )
    assert not output.exists()


def test_wrong_physical_schema_and_no_overwrite_fail_closed(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    wrong = dict(inputs)
    wrong["copy_number"] = inputs["structural_variants"]
    output = tmp_path / "wrong-schema-output"
    with pytest.raises(PortableEvidenceError, match="wrong physical schema"):
        stage_release_c_extensions(output, wrong, input_gates, receipts, adapters, subject_binding)
    assert not output.exists()

    valid_output = tmp_path / "valid-output"
    stage_release_c_extensions(
        valid_output, inputs, input_gates, receipts, adapters, subject_binding
    )
    original_manifest = (valid_output / MANIFEST_FILENAME).read_bytes()
    with pytest.raises(PortableEvidenceError, match="refusing to replace"):
        stage_release_c_extensions(
            valid_output, inputs, input_gates, receipts, adapters, subject_binding
        )
    assert (valid_output / MANIFEST_FILENAME).read_bytes() == original_manifest


def test_independent_validator_detects_table_and_manifest_tampering(tmp_path: Path) -> None:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    output = tmp_path / "tamper-output"
    stage_release_c_extensions(output, inputs, input_gates, receipts, adapters, subject_binding)
    table = output / "structural_variants.parquet"
    table.write_bytes(table.read_bytes() + b"tamper")
    with pytest.raises(PortableEvidenceError, match="digest changed"):
        validate_release_c_extension_stage(output, expected_reference_sha256=REFERENCE)

    second_output = tmp_path / "manifest-tamper-output"
    stage_release_c_extensions(
        second_output, inputs, input_gates, receipts, adapters, subject_binding
    )
    manifest_path = second_output / MANIFEST_FILENAME
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source_rows_removed"] = 1
    write_json(manifest_path, manifest)
    with pytest.raises(PortableEvidenceError, match="privacy/losslessness"):
        validate_release_c_extension_stage(second_output, expected_reference_sha256=REFERENCE)

    third_output = tmp_path / "receipt-tamper-output"
    stage_release_c_extensions(
        third_output, inputs, input_gates, receipts, adapters, subject_binding
    )
    third_manifest_path = third_output / MANIFEST_FILENAME
    third_manifest = json.loads(third_manifest_path.read_text(encoding="utf-8"))
    third_manifest["caller_runs"][0]["receipt"]["raw_artifact"]["components"][0]["bytes"] += 1
    third_manifest["extension_set_identity_sha256"] = canonical_json_sha256(
        {
            key: value
            for key, value in third_manifest.items()
            if key != "extension_set_identity_sha256"
        }
    )
    write_json(third_manifest_path, third_manifest)
    with pytest.raises(PortableEvidenceError, match="caller receipt identity changed"):
        validate_release_c_extension_stage(third_output, expected_reference_sha256=REFERENCE)

    fourth_output = tmp_path / "adapter-arithmetic-tamper-output"
    stage_release_c_extensions(
        fourth_output, inputs, input_gates, receipts, adapters, subject_binding
    )
    fourth_manifest_path = fourth_output / MANIFEST_FILENAME
    fourth_manifest = json.loads(fourth_manifest_path.read_text(encoding="utf-8"))
    manta = next(
        row["adapter_summary"]
        for row in fourth_manifest["adapter_runs"]
        if row["adapter_summary"]["profile_id"] == "manta-germline-1.6.0"
    )
    manta["parsing"][0]["summary"]["omitted_noncanonical_records"] += 1
    manta["adapter_identity_sha256"] = canonical_json_sha256(
        {key: value for key, value in manta.items() if key != "adapter_identity_sha256"}
    )
    fourth_manifest["extension_set_identity_sha256"] = canonical_json_sha256(
        {
            key: value
            for key, value in fourth_manifest.items()
            if key != "extension_set_identity_sha256"
        }
    )
    write_json(fourth_manifest_path, fourth_manifest)
    with pytest.raises(PortableEvidenceError, match="Manta raw/typed/scope-omission"):
        validate_release_c_extension_stage(fourth_output, expected_reference_sha256=REFERENCE)
