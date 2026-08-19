from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import mt_workflow
from genome_annotation_builder.release_c.mt_workflow import (
    MitochondrialWorkflowError,
    load_policy,
    validate_resource_root,
)


def write_policy(path: Path, **updates: object) -> None:
    payload: dict[str, object] = {
        "schema_version": 1,
        "source_contig": "MT",
        "threads": 2,
        "max_reads_per_alignment_start": 75,
        "max_alt_allele_count": 4,
        "min_allele_fraction": 0.0,
        "contamination_estimate": None,
    }
    payload.update(updates)
    path.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")


def test_policy_is_strict_finite_and_identifier_free(tmp_path: Path) -> None:
    path = tmp_path / "policy.json"
    write_policy(path, source_contig="chrM", contamination_estimate=0.001)
    assert load_policy(path) == {
        "schema_version": 1,
        "source_contig": "chrM",
        "threads": 2,
        "max_reads_per_alignment_start": 75,
        "max_alt_allele_count": 4,
        "min_allele_fraction": 0.0,
        "contamination_estimate": 0.001,
    }
    write_policy(path, contamination_estimate=float("nan"))
    with pytest.raises(MitochondrialWorkflowError, match="policy_contamination_invalid"):
        load_policy(path)
    write_policy(path, unexpected_subject="not permitted")
    with pytest.raises(MitochondrialWorkflowError, match="policy_schema_invalid"):
        load_policy(path)


def test_resource_boundary_checks_exact_rotation_and_required_members(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    normal = "ACGTACGTACGT"
    shift = 5
    shifted = normal[shift:] + normal[:shift]
    monkeypatch.setattr(mt_workflow, "MT_LENGTH", len(normal))
    monkeypatch.setattr(mt_workflow, "SHIFT_BASES", shift)
    monkeypatch.setattr(
        mt_workflow,
        "MT_MD5",
        hashlib.md5(normal.encode(), usedforsecurity=False).hexdigest(),
    )
    monkeypatch.setattr(
        mt_workflow,
        "SHIFTED_MT_MD5",
        hashlib.md5(shifted.encode(), usedforsecurity=False).hexdigest(),
    )
    for name in mt_workflow.RESOURCE_FILES.values():
        (tmp_path / name).write_text("resource\n", encoding="ascii")
    (tmp_path / "mt.fa").write_text(f">chrM\n{normal}\n", encoding="ascii")
    (tmp_path / "mt.shifted.fa").write_text(f">chrM\n{shifted}\n", encoding="ascii")
    validated = validate_resource_root(tmp_path)
    assert set(validated) == set(mt_workflow.RESOURCE_FILES)

    (tmp_path / "mt.shifted.fa").write_text(f">chrM\n{normal}\n", encoding="ascii")
    with pytest.raises(MitochondrialWorkflowError, match="shifted_mt_reference_md5_mismatch"):
        validate_resource_root(tmp_path)


def test_subprocess_boundary_uses_argv_and_retains_private_log(tmp_path: Path) -> None:
    logs = tmp_path / "logs"
    logs.mkdir()
    output = tmp_path / "output.txt"
    mt_workflow._run(
        "synthetic",
        [sys.executable, "-c", "print('ok')"],
        cwd=tmp_path,
        logs=logs,
        stdout_path=output,
    )
    assert output.read_text(encoding="utf-8") == "ok\n"
    assert (logs / "synthetic.log").is_file()

    with pytest.raises(MitochondrialWorkflowError, match="step_failed:failure:7"):
        mt_workflow._run(
            "failure",
            [sys.executable, "-c", "raise SystemExit(7)"],
            cwd=tmp_path,
            logs=logs,
        )


def test_pinned_gatk_command_is_no_shell_and_uses_private_tmpfs() -> None:
    assert mt_workflow._gatk("8g", "Mutect2", "--mitochondria-mode") == [
        "gatk",
        "--java-options",
        "-Xmx8g -Djava.io.tmpdir=/work/tmp",
        "Mutect2",
        "--mitochondria-mode",
    ]


def test_vcf_scientific_record_stream_excludes_volatile_header(tmp_path: Path) -> None:
    one = tmp_path / "one.vcf"
    two = tmp_path / "two.vcf"
    record = "chrM\t20\t.\tT\tA\t.\tPASS\t.\tGT\t1\n"
    one.write_text("##generated=first\n#CHROM\n" + record, encoding="ascii")
    two.write_text("##generated=second\n#CHROM\n" + record, encoding="ascii")
    assert mt_workflow._vcf_record_stream(one) == mt_workflow._vcf_record_stream(two)
    assert mt_workflow._vcf_record_stream(one)[1] == 1
