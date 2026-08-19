from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c import cnvpytor_combined_hdf5_probe
from genome_annotation_builder.release_c.cnvpytor_combined_adapter import (
    CnvpytorCombinedProvenance,
    materialize_cnvpytor_combined_calls,
    parse_cnvpytor_combined_calls,
)
from genome_annotation_builder.release_c.cnvpytor_combined_hdf5_probe import (
    CombinedProbeError,
    _write_json_atomic,
    required_combined_signal_names,
    validate_combined_row,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

REFERENCE = "1" * 64


def combined_row(
    *,
    heterozygous_sites: int = 20,
    strict_mask_fraction: float = 0.9,
    top_likelihood: float = 0.9,
    reference_gap_fraction: float = 0.0,
) -> list[float]:
    base = [
        20.0,
        -1.0,
        100001.0,
        160000.0,
        60000.0,
        0.5,
        1e-8,
        0.8,
        0.1,
        0.1,
        0.1,
        60.0,
        0.2,
        1e-8,
        1e-6,
        0.0,
        float(heterozygous_sites),
        5.0,
        reference_gap_fraction,
        0.0,
        strict_mask_fraction,
    ]
    models = [
        1.0,
        1.0,
        0.0,
        top_likelihood,
        1.0,
        2.0,
        1.0,
        1.0,
        0.1,
        1.0,
        *([0.0] * 40),
    ]
    row = base + models
    assert len(row) == 71
    return row


def write_artifacts(
    root: Path,
    *,
    heterozygous_sites: int = 20,
    strict_mask_fraction: float = 0.9,
    hdf5_top_likelihood: float = 0.9,
    hdf5_reference_gap_fraction: float = 0.0,
    tsv_reference_gap_text: str = "0.0000",
    reported_rd_text: str = "1.000000e-08",
) -> tuple[Path, Path, Path, Path, CnvpytorCombinedProvenance]:
    root.mkdir(parents=True, exist_ok=True)
    pytor = root / "combined.pytor"
    pytor.write_bytes(b"synthetic combined HDF5 stand-in")
    pytor_sha256 = hashlib.sha256(pytor.read_bytes()).hexdigest()
    baf_vcf = root / "baf.vcf.gz"
    baf_vcf.write_bytes(b"synthetic accepted heterozygous-SNV source")
    baf_sha256 = hashlib.sha256(baf_vcf.read_bytes()).hexdigest()
    calls = root / "combined.tsv"
    fields = [
        "deletion",
        "chr1:100001-160000",
        "60000",
        "0.5000",
        "1.000000e-08",
        "8.000000e-01",
        "1.000000e-01",
        "1.000000e-01",
        "0.1000",
        tsv_reference_gap_text,
        "0.0000",
        f"{strict_mask_fraction:.4f}",
        "1000",
        "60",
        "0.2000",
        reported_rd_text,
        "1.000000e-06",
        str(heterozygous_sites),
        "5",
        "1",
        "CN1/CN0",
        "9.000000e-01",
        "1.0000",
        "2",
        "CN1/CN1",
        "1.000000e-01",
        "1.0000",
    ]
    assert len(fields) == 27
    calls.write_text("\t".join(fields) + "\n", encoding="ascii")
    audit = root / "combined-audit.json"
    dataset_name = "calls_2d_1_1000_mask_GC_mask"
    audit.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": "PASS",
                "caller_model": "combined_rd_baf",
                "cnvpytor_version": "1.3.2",
                "pytor_sha256": pytor_sha256,
                "bin_size": 1000,
                "rd_flags_suffix": "_GC_mask",
                "snp_flags_suffix": "_mask",
                "volatile_attribute_values_persisted": False,
                "private_paths_persisted": False,
                "call_datasets": {
                    dataset_name: {
                        "contig": "1",
                        "bin_size": 1000,
                        "rd_flags_suffix": "_GC_mask",
                        "snp_flags_suffix": "_mask",
                        "rows": [
                            combined_row(
                                heterozygous_sites=heterozygous_sites,
                                strict_mask_fraction=strict_mask_fraction,
                                top_likelihood=hdf5_top_likelihood,
                                reference_gap_fraction=hdf5_reference_gap_fraction,
                            )
                        ],
                    }
                },
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    provenance = CnvpytorCombinedProvenance(
        reference_sha256=REFERENCE,
        technical_source_role="primary",
        caller_run_identity_sha256="2" * 64,
        resource_manifest_sha256="3" * 64,
        source_alignment_sha256="4" * 64,
        raw_artifact_id="cnvpytor-combined-artifacts",
        raw_artifact_sha256="6" * 64,
        source_baf_vcf_sha256=baf_sha256,
        bin_size=1000,
        rd_flags_suffix="_GC_mask",
        snp_flags_suffix="_mask",
        rd_mask_applied=True,
        baf_mask_applied=True,
    )
    return calls, pytor, audit, baf_vcf, provenance


def test_combined_rd_baf_reconciliation_is_typed_deterministic_and_explicit(
    tmp_path: Path,
) -> None:
    calls, pytor, audit, baf_vcf, provenance = write_artifacts(tmp_path)
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    result_one = materialize_cnvpytor_combined_calls(
        calls, pytor, audit, baf_vcf, first, provenance
    )
    result_two = materialize_cnvpytor_combined_calls(
        calls, pytor, audit, baf_vcf, second, provenance
    )
    assert result_one["parsing"]["typed_rows"] == 1
    assert result_one["parsing"]["supported_rows"] == 1
    assert result_one["parsing"]["reported_rd_p_value_duplicates_combined"] is True
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(first)
    assert table.schema == schema_for("copy_number")
    row = table.to_pylist()[0]
    assert row["caller_model"] == "combined_rd_baf"
    assert row["call_type"] == "deletion"
    assert row["copy_number"] == pytest.approx(1)
    assert row["read_depth_copy_number"] == pytest.approx(1)
    assert row["major_copy_number"] == pytest.approx(1)
    assert row["minor_copy_number"] == pytest.approx(0)
    assert row["baf"] is None
    assert row["baf_deviation_from_half"] == pytest.approx(0.2)
    assert row["caller_reported_rd_p_value"] == pytest.approx(row["combined_p_value"])
    assert len(row["copy_number_models"]) == 10
    assert row["copy_number_models"][0]["reported_in_tsv"] is True
    assert row["copy_number_models"][2]["padding"] is True
    assert row["source_baf_vcf_sha256"] == provenance.source_baf_vcf_sha256
    assert row["raw_artifact_sha256"] == provenance.raw_artifact_sha256
    assert result_one["parsing"]["caller_raw_artifact_sha256"] == (provenance.raw_artifact_sha256)
    assert result_one["parsing"]["adapter_input_set_sha256"] != (provenance.raw_artifact_sha256)


def test_low_baf_support_and_mask_coverage_remain_review(tmp_path: Path) -> None:
    calls, pytor, audit, baf_vcf, provenance = write_artifacts(
        tmp_path, heterozygous_sites=1, strict_mask_fraction=0.1
    )
    rows, summary = parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, provenance)
    assert summary["review_rows"] == 1
    assert rows[0]["confidence_tier"] == "REVIEW"
    assert "HET_SITES_LT_5" in rows[0]["filter_reasons"]
    assert "STRICT_MASK_FRACTION_LT_0.5" in rows[0]["filter_reasons"]


def test_four_decimal_half_unit_reconciliation_boundary_is_inclusive(
    tmp_path: Path,
) -> None:
    calls, pytor, audit, baf_vcf, provenance = write_artifacts(
        tmp_path / "inclusive",
        hdf5_reference_gap_fraction=0.03125,
        tsv_reference_gap_text="0.0312",
    )
    rows, summary = parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, provenance)
    assert len(rows) == 1
    assert summary["record_reconciliation_delta"] == 0

    calls, pytor, audit, baf_vcf, provenance = write_artifacts(
        tmp_path / "outside",
        hdf5_reference_gap_fraction=0.0312501,
        tsv_reference_gap_text="0.0312",
    )
    with pytest.raises(PortableEvidenceError, match="does not match"):
        parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, provenance)


def test_pvalue_hdf5_baf_digest_and_mask_boundaries_fail_closed(tmp_path: Path) -> None:
    calls, pytor, audit, baf_vcf, provenance = write_artifacts(
        tmp_path, reported_rd_text="0.00000001"
    )
    with pytest.raises(PortableEvidenceError, match="reported RD p-value"):
        parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, provenance)

    calls, pytor, audit, baf_vcf, provenance = write_artifacts(
        tmp_path / "mismatch", hdf5_top_likelihood=0.8
    )
    with pytest.raises(PortableEvidenceError, match="does not match"):
        parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, provenance)

    wrong_baf = CnvpytorCombinedProvenance(
        **{**provenance.__dict__, "source_baf_vcf_sha256": "0" * 64}
    )
    with pytest.raises(PortableEvidenceError, match="BAF VCF differs"):
        parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, wrong_baf)

    wrong_artifact = CnvpytorCombinedProvenance(
        **{**provenance.__dict__, "raw_artifact_sha256": "invalid"}
    )
    with pytest.raises(PortableEvidenceError, match="artifact SHA-256 is invalid"):
        parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, wrong_artifact)

    no_mask = CnvpytorCombinedProvenance(
        **{
            **provenance.__dict__,
            "snp_flags_suffix": "",
            "baf_mask_applied": False,
        }
    )
    with pytest.raises(PortableEvidenceError, match="strict mask is mandatory"):
        parse_cnvpytor_combined_calls(calls, pytor, audit, baf_vcf, no_mask)


def test_hdf5_row_probe_rejects_type_and_model_corruption() -> None:
    row = combined_row()
    assert validate_combined_row(row) == row

    wrong_type = list(row)
    wrong_type[1] = 1.0
    with pytest.raises(CombinedProbeError, match="type differs"):
        validate_combined_row(wrong_type)

    wrong_model = list(row)
    wrong_model[22] = 2.0
    with pytest.raises(CombinedProbeError, match="genotype is inconsistent"):
        validate_combined_row(wrong_model)


def test_hdf5_probe_matches_cnvpytor_132_combined_signal_contract() -> None:
    names, call_name = required_combined_signal_names(
        "1",
        bin_size=1000,
        rd_flags_suffix="_GC_mask",
        snp_flags_suffix="_mask",
    )
    assert call_name == "calls_2d_1_1000_mask_GC_mask"
    assert names == {
        "1_rd_p",
        "1_rd_u",
        "his_rd_p_1_1000_GC_mask",
        "his_rd_p_1_1000",
        "his_rd_u_1_1000",
        "snp_bafc_01_1_1000_mask",
        "snp_bafc_10_1_1000_mask",
        "snp_bafc_11_1_1000_mask",
        "snp_baf_1_1000_mask",
        "calls_2d_1_1000_mask_GC_mask",
    }
    assert "his_rd_u_1_1000_GC_mask" not in names


def test_hdf5_audit_publish_is_no_replace_and_cli_failure_is_sanitized(
    tmp_path: Path,
) -> None:
    output = tmp_path / "audit.json"
    _write_json_atomic(output, {"schema_version": 1, "status": "PASS"})
    original = output.read_bytes()
    with pytest.raises(CombinedProbeError, match="refusing to replace"):
        _write_json_atomic(output, {"schema_version": 1, "status": "CHANGED"})
    assert output.read_bytes() == original

    private_token = "privateid"
    private_root = tmp_path / private_token
    private_root.mkdir()
    failed = subprocess.run(
        [
            sys.executable,
            str(Path(cnvpytor_combined_hdf5_probe.__file__)),
            "--input",
            str(private_root / "missing.pytor"),
            "--output",
            str(private_root / "audit.json"),
            "--bin-size",
            "10000",
            "--rd-flags-suffix",
            "_GC_mask",
            "--snp-flags-suffix",
            "_mask",
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
