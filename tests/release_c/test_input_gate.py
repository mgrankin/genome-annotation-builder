from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from genome_annotation_builder.release_c.input_gate import (
    InputGateFailure,
    InputGatePolicy,
    _assert_portable,
    run_input_gate,
)

SYNTHETIC_SUBJECT = "SYNTHETIC_SUBJECT"
SYNTHETIC_DATASET = "synthetic-dataset-001"
SYNTHETIC_LAB = "SYNTHETIC_LAB"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def run(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(f"synthetic setup failed with exit {completed.returncode}")


def file_spec(path: Path, index: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "index_path": str(index),
        "index_bytes": index.stat().st_size,
        "index_sha256": sha256_file(index),
    }


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def make_inputs(
    tmp_path: Path, *, baf_sample: str = SYNTHETIC_SUBJECT
) -> tuple[Path, dict[str, object]]:
    reference = tmp_path / "reference.fa"
    reference.write_text(">1\n" + "A" * 200 + "\n", encoding="ascii")
    run(["samtools", "faidx", str(reference)])
    fai = tmp_path / "reference.fa.fai"
    dictionary = tmp_path / "reference.dict"
    run(["samtools", "dict", "-o", str(dictionary), str(reference)])

    sam = tmp_path / "reads.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:1\tLN:200\n"
        f"@RG\tID:rg1\tSM:{SYNTHETIC_SUBJECT}\tPL:ILLUMINA\n"
        "pair1\t99\t1\t10\t60\t10M\t=\t30\t30\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg1\n"
        "pair1\t147\t1\t30\t60\t10M\t=\t10\t-30\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg1\n",
        encoding="ascii",
    )
    alignment = tmp_path / "alignment.bam"
    alignment_index = tmp_path / "alignment.bam.bai"
    run(["samtools", "view", "-b", "-o", str(alignment), str(sam)])
    run(["samtools", "index", "-o", str(alignment_index), str(alignment)])

    raw_vcf = tmp_path / "baf.vcf"
    raw_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=200>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        f"{baf_sample}\n"
        "1\t20\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        encoding="ascii",
    )
    baf = tmp_path / "baf.vcf.gz"
    baf_index = tmp_path / "baf.vcf.gz.tbi"
    run(["bcftools", "view", "--no-version", "-Oz", "-o", str(baf), str(raw_vcf)])
    run(["bcftools", "index", "--tbi", "-o", str(baf_index), str(baf)])

    manifest = {
        "schema_version": 1,
        "status": "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH",
        "caller_alignment_policy": {
            "cnvpytor": "original_markdup_bam",
            "manta": "original_markdup_bam",
        },
        "reference": {
            "assembly": "GRCh38 synthetic control",
            "fasta_path": str(reference),
            "fasta_bytes": reference.stat().st_size,
            "fasta_sha256": sha256_file(reference),
            "fai_path": str(fai),
            "fai_bytes": fai.stat().st_size,
            "fai_sha256": sha256_file(fai),
            "dictionary_path": str(dictionary),
            "dictionary_bytes": dictionary.stat().st_size,
            "dictionary_sha256": sha256_file(dictionary),
        },
        "datasets": [
            {
                "dataset_id": SYNTHETIC_DATASET,
                "logical_subject": SYNTHETIC_SUBJECT,
                "laboratory": SYNTHETIC_LAB,
                "original_markdup_bam": file_spec(alignment, alignment_index),
                "baf_vcf": {
                    **file_spec(baf, baf_index),
                    "sample": baf_sample,
                },
            }
        ],
    }
    manifest_path = tmp_path / "private-manifest.json"
    write_json(manifest_path, manifest)
    return manifest_path, manifest


def forbidden_values(
    manifest_path: Path, output_path: Path, manifest: dict[str, object]
) -> set[str]:
    dataset = manifest["datasets"][0]
    assert isinstance(dataset, dict)
    reference = manifest["reference"]
    assert isinstance(reference, dict)
    values = {
        SYNTHETIC_SUBJECT,
        SYNTHETIC_DATASET,
        SYNTHETIC_LAB,
        str(dataset["baf_vcf"]["sample"]),
        str(manifest_path),
        str(output_path),
        str(reference["fasta_path"]),
        str(dataset["original_markdup_bam"]["path"]),
        str(dataset["baf_vcf"]["path"]),
    }
    values.update(
        hashlib.sha256(value.encode()).hexdigest()
        for value in (SYNTHETIC_SUBJECT, SYNTHETIC_DATASET, SYNTHETIC_LAB)
    )
    return values


def test_generic_input_gate_rehashes_and_emits_portable_report(tmp_path: Path) -> None:
    manifest_path, manifest = make_inputs(tmp_path)
    output = tmp_path / "portable-input-gate.json"
    report = run_input_gate(
        manifest_path,
        SYNTHETIC_DATASET,
        "cnvpytor",
        output,
        policy=InputGatePolicy(minimum_free_bytes=1, indexed_probe_bases=100),
    )
    assert report["status"] == "PASS"
    assert report["input_manifest_state"].startswith("PARTIAL_PASS")
    assert report["selected_alignment"]["sort_order"] == "coordinate"
    assert report["selected_alignment"]["indexed_probe_records"] == 2
    assert report["reference"]["indexed_probe_succeeded"] is True
    assert report["baf_vcf"]["sample_identity_match"] is True
    assert report["baf_vcf"]["records"] == 1
    assert report["logical_subject_count"] == 1
    assert report["subject_identifier_persisted"] is False
    assert report["private_manifest_digest_persisted"] is False
    assert "manifest_sha256" not in report
    payload = output.read_text(encoding="utf-8")
    assert not any(value in payload for value in forbidden_values(manifest_path, output, manifest))


def test_complete_manifest_state_and_caller_without_baf_are_accepted(tmp_path: Path) -> None:
    manifest_path, manifest = make_inputs(tmp_path)
    manifest["status"] = "PASS"
    manifest["reference"].pop("fai_bytes")
    manifest["reference"].pop("dictionary_bytes")
    write_json(manifest_path, manifest)
    output = tmp_path / "manta-gate.json"
    report = run_input_gate(
        manifest_path,
        SYNTHETIC_DATASET,
        "manta",
        output,
        policy=InputGatePolicy(minimum_free_bytes=1),
    )
    assert report["status"] == "PASS"
    assert report["input_manifest_state"] == "PASS"
    assert report["baf_vcf"] is None


def test_declared_technical_source_vcf_alias_is_verified_without_persistence(
    tmp_path: Path,
) -> None:
    technical_alias = "SYNTHETIC_TECHNICAL_SOURCE"
    manifest_path, manifest = make_inputs(tmp_path, baf_sample=technical_alias)
    output = tmp_path / "aliased-baf-gate.json"
    report = run_input_gate(
        manifest_path,
        SYNTHETIC_DATASET,
        "cnvpytor",
        output,
        policy=InputGatePolicy(minimum_free_bytes=1),
    )
    assert report["status"] == "PASS"
    assert report["baf_vcf"]["sample_identity_match"] is True
    payload = output.read_text(encoding="utf-8")
    assert technical_alias not in payload
    assert not any(value in payload for value in forbidden_values(manifest_path, output, manifest))


@pytest.mark.parametrize(
    ("mutation", "error_code"),
    [
        ("wrong_baf_sample", "baf_vcf_sample_identity_mismatch"),
        ("wrong_alignment_digest", "alignment_digest_changed"),
        ("wrong_dictionary_digest", "reference_dictionary_digest_changed"),
        ("invalid_manifest_state", "input_manifest_state_invalid"),
    ],
)
def test_failures_are_retained_without_private_bindings(
    tmp_path: Path, mutation: str, error_code: str
) -> None:
    manifest_path, manifest = make_inputs(tmp_path)
    dataset = manifest["datasets"][0]
    assert isinstance(dataset, dict)
    if mutation == "wrong_baf_sample":
        dataset["baf_vcf"]["sample"] = "DIFFERENT_SYNTHETIC_SUBJECT"
    elif mutation == "wrong_alignment_digest":
        dataset["original_markdup_bam"]["sha256"] = "0" * 64
    elif mutation == "wrong_dictionary_digest":
        manifest["reference"]["dictionary_sha256"] = "0" * 64
    elif mutation == "invalid_manifest_state":
        manifest["status"] = "RUNNING"
    write_json(manifest_path, manifest)
    output = tmp_path / f"failed-{mutation}.json"
    with pytest.raises(InputGateFailure, match=error_code):
        run_input_gate(
            manifest_path,
            SYNTHETIC_DATASET,
            "cnvpytor",
            output,
            policy=InputGatePolicy(minimum_free_bytes=1),
        )
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["status"] == "FAIL_RETAINED"
    assert report["error_code"] == error_code
    assert report["error"] == "input gate failed; inspect the authorized private process log"
    payload = output.read_text(encoding="utf-8")
    assert not any(value in payload for value in forbidden_values(manifest_path, output, manifest))


def test_policy_and_output_replacement_fail_before_writing(tmp_path: Path) -> None:
    manifest_path, _ = make_inputs(tmp_path)
    invalid_output = tmp_path / "invalid-policy.json"
    with pytest.raises(InputGateFailure, match="input_gate_policy_invalid"):
        run_input_gate(
            manifest_path,
            SYNTHETIC_DATASET,
            "manta",
            invalid_output,
            policy=InputGatePolicy(minimum_free_bytes=-1),
        )
    assert not invalid_output.exists()

    existing = tmp_path / "existing.json"
    existing.write_text("do not replace\n", encoding="ascii")
    with pytest.raises(InputGateFailure, match="output_already_exists"):
        run_input_gate(manifest_path, SYNTHETIC_DATASET, "manta", existing)
    assert existing.read_text(encoding="ascii") == "do not replace\n"


def test_short_private_binding_is_token_aware() -> None:
    _assert_portable({"status": "PASS", "message": "PASSING"}, {"A"})
    with pytest.raises(InputGateFailure, match="portable_report_privacy_violation"):
        _assert_portable({"status": "PASS", "artifact": "artifact-A"}, {"A"})


def test_dictionary_semantic_mismatch_fails_even_with_matching_digest(tmp_path: Path) -> None:
    manifest_path, manifest = make_inputs(tmp_path)
    reference = manifest["reference"]
    assert isinstance(reference, dict)
    dictionary = Path(str(reference["dictionary_path"]))
    dictionary.write_text(
        dictionary.read_text(encoding="utf-8").replace("LN:200", "LN:199"),
        encoding="utf-8",
    )
    reference["dictionary_bytes"] = dictionary.stat().st_size
    reference["dictionary_sha256"] = sha256_file(dictionary)
    write_json(manifest_path, manifest)
    output = tmp_path / "dictionary-mismatch.json"
    with pytest.raises(InputGateFailure, match="reference_index_dictionary_boundary_mismatch"):
        run_input_gate(
            manifest_path,
            SYNTHETIC_DATASET,
            "manta",
            output,
            policy=InputGatePolicy(minimum_free_bytes=1),
        )
    assert json.loads(output.read_text(encoding="utf-8"))["error_code"] == (
        "reference_index_dictionary_boundary_mismatch"
    )
