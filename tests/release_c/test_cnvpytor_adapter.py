from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.cnvpytor_adapter import (
    CnvpytorProvenance,
    materialize_cnvpytor_mean_shift_calls,
    parse_cnvpytor_mean_shift_calls,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

REFERENCE = "1" * 64
PROVENANCE = CnvpytorProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="cnvpytor-mean-shift-artifacts",
    bin_size=1000,
    flags_suffix="_GC_mask",
    rd_mask_applied=True,
    x_ploidy=2,
    y_ploidy=1,
)


def write_artifacts(
    root: Path,
    *,
    hdf5_depth: float = 0.5,
    tsv_depth: float = 0.5,
    flags_suffix: str = "_GC_mask",
) -> tuple[Path, Path, Path]:
    pytor = root / "calls.pytor"
    pytor.write_bytes(b"synthetic HDF5 stand-in bound only by the audit test")
    pytor_sha256 = hashlib.sha256(pytor.read_bytes()).hexdigest()
    calls = root / "calls.tsv"
    calls.write_text(
        f"deletion\tchr1:100001-160000\t60000\t{tsv_depth:.4f}\t1.000000e-08\t"
        "1.000000e-09\t1.000000e-07\t1.000000e-08\t0.1000\t0.0000\t200000\n",
        encoding="ascii",
    )
    audit = root / "audit.json"
    payload = {
        "schema_version": 1,
        "status": "PASS",
        "cnvpytor_version": "1.3.2",
        "pytor_sha256": pytor_sha256,
        "bin_size": 1000,
        "flags_suffix": flags_suffix,
        "volatile_attribute_values_persisted": False,
        "private_paths_persisted": False,
        "call_datasets": {
            f"calls_1_1000{flags_suffix}": {
                "contig": "1",
                "bin_size": 1000,
                "flags_suffix": flags_suffix,
                "rows": [
                    [
                        12.0,
                        -1.0,
                        100001.0,
                        160000.0,
                        60000.0,
                        hdf5_depth,
                        1e-8,
                        1e-9,
                        1e-7,
                        1e-8,
                        0.1,
                        0.0,
                        200000.0,
                    ]
                ],
            }
        },
    }
    audit.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    return calls, pytor, audit


def test_tsv_hdf5_reconciliation_and_copy_number_basis(tmp_path: Path) -> None:
    calls, pytor, audit = write_artifacts(tmp_path)
    rows, summary = parse_cnvpytor_mean_shift_calls(calls, pytor, audit, PROVENANCE)
    assert summary["tsv_calls"] == summary["hdf5_calls"] == summary["typed_rows"] == 1
    assert summary["record_reconciliation_delta"] == 0
    assert summary["supported_rows"] == 1
    row = rows[0]
    assert row["call_type"] == "deletion"
    assert row["copy_number"] == pytest.approx(1.0)
    assert row["copy_number_basis"] == "normalized_read_depth_times_declared_ploidy"
    assert row["log2_ratio"] == pytest.approx(-1.0)
    assert row["bins_observed"] == 60
    assert row["rd_mask_applied"] is True


def test_audited_calls_control_materializes_deterministically(
    tmp_path: Path,
) -> None:
    calls, pytor, audit = write_artifacts(
        tmp_path,
        hdf5_depth=0.5244,
        tsv_depth=0.5244,
        flags_suffix="_GC",
    )
    provenance = CnvpytorProvenance(
        **{
            **PROVENANCE.__dict__,
            "bin_size": 1000,
            "flags_suffix": "_GC",
            "rd_mask_applied": False,
        }
    )
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    result_one = materialize_cnvpytor_mean_shift_calls(calls, pytor, audit, first, provenance)
    result_two = materialize_cnvpytor_mean_shift_calls(calls, pytor, audit, second, provenance)
    assert result_one["parsing"]["typed_rows"] == 1
    assert result_one["parsing"]["supported_rows"] == 1
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(first)
    assert table.schema == schema_for("copy_number")
    row = table.to_pylist()[0]
    assert row["chrom"] == "chr1"
    assert row["normalized_read_depth"] == pytest.approx(0.5244)
    assert row["copy_number"] == pytest.approx(1.0488)
    assert row["filter_reasons"] == []


def test_metric_and_artifact_mismatches_fail_closed(tmp_path: Path) -> None:
    calls, pytor, audit = write_artifacts(tmp_path, hdf5_depth=0.7)
    with pytest.raises(PortableEvidenceError, match="does not match"):
        parse_cnvpytor_mean_shift_calls(calls, pytor, audit, PROVENANCE)

    payload = json.loads(audit.read_text(encoding="utf-8"))
    payload["pytor_sha256"] = "0" * 64
    audit.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="does not bind"):
        parse_cnvpytor_mean_shift_calls(calls, pytor, audit, PROVENANCE)


def test_logged_error_or_sample_column_cannot_parse_as_calls(tmp_path: Path) -> None:
    calls, pytor, audit = write_artifacts(tmp_path)
    calls.write_text("ERROR\tinput unavailable\n", encoding="ascii")
    with pytest.raises(PortableEvidenceError, match="does not have 11 columns"):
        parse_cnvpytor_mean_shift_calls(calls, pytor, audit, PROVENANCE)
