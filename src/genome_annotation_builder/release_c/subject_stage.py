"""Assemble one subject's completed caller matrix into a portable extension stage."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from typing import Any

import pyarrow.parquet as pq

from .cnv_concordance import annotate_cnv_concordance
from .extension_stage import (
    stage_release_c_extensions,
    subject_binding_from_private_manifest,
    validate_release_c_extension_stage,
)
from .hla_concordance import annotate_hla_concordance
from .mt_concordance import annotate_mitochondrial_replicate_concordance
from .pgx_concordance import annotate_pgx_concordance
from .portable_model import (
    DOMAINS,
    PortableEvidenceError,
    canonical_row_multiset,
    schema_for,
    sha256_file,
    validate_rows,
    write_extension,
)
from .readback_concordance import annotate_readback_replicate_concordance
from .repeat_concordance import annotate_repeat_replicate_concordance
from .sv_concordance import annotate_sv_concordance

PLAN_NAME = "subject-stage-plan.private.json"
COMPLETE_NAME = "subject-stage-complete.private.json"
SUMMARY_NAME = "subject-stage-summary.json"


class SubjectStageError(RuntimeError):
    """The subject stage cannot progress without violating an immutable boundary."""


def _fail(message: str) -> None:
    raise SubjectStageError(message)


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _load_object(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        _fail(f"{label}_missing_or_symlinked")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise SubjectStageError(f"{label}_unreadable") from exc
    if not isinstance(value, dict):
        _fail(f"{label}_not_an_object")
    return value


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def _atomic_no_replace(path: Path, value: object) -> None:
    if path.exists():
        _fail("immutable_subject_stage_checkpoint_already_exists")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, raw_name = tempfile.mkstemp(prefix=f".{path.name}.partial.", dir=path.parent)
    temporary = Path(raw_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(_json_bytes(value))
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as exc:
            raise SubjectStageError("immutable_subject_stage_checkpoint_already_exists") from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _disk_gate(path: Path, minimum_free_bytes: int) -> None:
    if minimum_free_bytes < 1:
        _fail("subject_stage_minimum_free_bytes_invalid")
    candidate = path
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    if not candidate.exists() or shutil.disk_usage(candidate).free < minimum_free_bytes:
        _fail("subject_stage_free_space_below_required_floor")


def _validate_canonical_identity(value: dict[str, Any], field: str, label: str) -> str:
    claimed = value.get(field)
    if not isinstance(claimed, str):
        _fail(f"{label}_identity_missing")
    unsigned = dict(value)
    unsigned.pop(field, None)
    if claimed != _canonical_sha256(unsigned):
        _fail(f"{label}_identity_changed")
    return claimed


def _attempt_paths(root: Path, ordinal: int, attempt: int) -> dict[str, Path]:
    attempt_root = root / f"job-{ordinal:04d}" / f"attempt-{attempt:04d}"
    return {
        "bindings": attempt_root / "bindings.private.json",
        "lifecycle": attempt_root / "lifecycle",
        "typed": attempt_root / "typed",
    }


def _load_completed_matrix(
    matrix_root: Path, manifest_path: Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = _load_object(matrix_root / "matrix-plan.private.json", "matrix_plan")
    runtime = _load_object(matrix_root / "matrix-runtime.private.json", "matrix_runtime")
    summary = _load_object(matrix_root / "matrix-summary.json", "matrix_summary")
    plan_identity = _validate_canonical_identity(plan, "plan_identity_sha256", "matrix_plan")
    _validate_canonical_identity(summary, "matrix_summary_identity_sha256", "matrix_summary")
    plan_jobs = plan.get("jobs")
    runtime_jobs = runtime.get("jobs")
    if (
        plan.get("schema_version") != 1
        or plan.get("status") != "READY"
        or not isinstance(plan_jobs, list)
        or runtime.get("schema_version") != 1
        or runtime.get("state") != "complete"
        or runtime.get("plan_identity_sha256") != plan_identity
        or not isinstance(runtime_jobs, list)
        or len(runtime_jobs) != len(plan_jobs)
        or runtime.get("completed_jobs") != len(plan_jobs)
        or summary.get("status") != "PASS"
        or summary.get("completed_jobs") != len(plan_jobs)
        or summary.get("job_count") != len(plan_jobs)
    ):
        _fail("matrix_not_complete_or_boundary_changed")
    identities = plan.get("input_identities")
    if not isinstance(identities, dict) or identities.get("private_manifest_sha256") != sha256_file(
        manifest_path
    ):
        _fail("matrix_private_manifest_boundary_changed")
    for plan_job, runtime_job in zip(plan_jobs, runtime_jobs, strict=True):
        if (
            not isinstance(plan_job, dict)
            or not isinstance(runtime_job, dict)
            or runtime_job.get("ordinal") != plan_job.get("ordinal")
            or runtime_job.get("dataset_id") != plan_job.get("dataset_id")
            or runtime_job.get("profile_id") != plan_job.get("profile_id")
            or runtime_job.get("state") != "complete"
            or not isinstance(runtime_job.get("attempt"), int)
            or runtime_job["attempt"] < 1
        ):
            _fail("matrix_runtime_job_boundary_changed")
    return plan, runtime


def _safe_adapter_file(typed_root: Path, raw_path: Any) -> Path:
    if not isinstance(raw_path, str):
        _fail("adapter_file_path_invalid")
    portable = PurePosixPath(raw_path)
    if portable.is_absolute() or len(portable.parts) != 1 or portable.name != raw_path:
        _fail("adapter_file_path_not_canonical")
    path = typed_root / raw_path
    if path.is_symlink() or not path.is_file():
        _fail("adapter_file_missing_or_symlinked")
    return path


def _configuration_count(domain: str, rows: list[dict[str, Any]]) -> int:
    if domain == "copy_number":
        return len(
            {
                (
                    row["caller_name"],
                    row["caller_model"],
                    row["technical_source_role"],
                    row["source_alignment_sha256"],
                    int(row["bin_size"] or 0),
                    row["caller_run_identity_sha256"],
                )
                for row in rows
            }
        )
    if domain in {"structural_variants", "hla_typing", "complex_pharmacogenes"}:
        return len(
            {
                (
                    row["caller_name"],
                    row["technical_source_role"],
                    row["source_alignment_sha256"],
                )
                for row in rows
            }
        )
    return len({row["source_alignment_sha256"] for row in rows})


CONCORDANCE: dict[
    str,
    Callable[..., tuple[list[dict[str, Any]], dict[str, Any]]],
] = {
    "structural_variants": annotate_sv_concordance,
    "copy_number": annotate_cnv_concordance,
    "repeat_expansions": annotate_repeat_replicate_concordance,
    "hla_typing": annotate_hla_concordance,
    "complex_pharmacogenes": annotate_pgx_concordance,
    "mitochondrial_heteroplasmy": annotate_mitochondrial_replicate_concordance,
    "read_backed_evidence": annotate_readback_replicate_concordance,
}


def _apply_concordance(
    domain: str,
    rows: list[dict[str, Any]],
    *,
    reference_sha256: str,
    forbidden_bindings: tuple[str, ...],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    configurations = _configuration_count(domain, rows)
    if configurations < 2:
        validated = validate_rows(
            domain,
            rows,
            expected_reference_sha256=reference_sha256,
            forbidden_bindings=forbidden_bindings,
        )
        return validated, {
            "schema_version": 1,
            "status": "NOT_APPLICABLE_SINGLE_SOURCE",
            "input_rows": len(validated),
            "output_rows": len(validated),
            "independent_configurations": configurations,
            "source_rows_removed": 0,
            "majority_vote_applied": False,
            "canonical_row_multiset": canonical_row_multiset(validated),
            "private_paths_persisted": False,
        }
    annotated, summary = CONCORDANCE[domain](
        rows,
        expected_reference_sha256=reference_sha256,
        forbidden_bindings=forbidden_bindings,
    )
    if summary.get("source_rows_removed", 0) != 0 or len(annotated) != len(rows):
        _fail("subject_stage_concordance_removed_source_rows")
    normalized_summary = dict(summary)
    normalized_summary.setdefault("source_rows_removed", 0)
    normalized_summary.setdefault("majority_vote_applied", False)
    return annotated, normalized_summary


def _collect_subject_inputs(
    manifest_path: Path,
    matrix_root: Path,
    dataset_ids: list[str],
) -> tuple[
    dict[str, list[dict[str, Any]]],
    list[Path],
    list[Path],
    list[Path],
    dict[str, Any],
    tuple[str, ...],
]:
    if not dataset_ids or len(dataset_ids) != len(set(dataset_ids)):
        _fail("subject_stage_dataset_ids_empty_or_duplicated")
    subject = subject_binding_from_private_manifest(manifest_path, dataset_ids)
    plan, runtime = _load_completed_matrix(matrix_root, manifest_path)
    selected = [
        (plan_job, runtime_job)
        for plan_job, runtime_job in zip(plan["jobs"], runtime["jobs"], strict=True)
        if plan_job.get("dataset_id") in set(dataset_ids)
    ]
    if (
        not selected
        or {str(plan_job.get("dataset_id")) for plan_job, _ in selected} != set(dataset_ids)
        or {str(plan_job.get("logical_subject")) for plan_job, _ in selected}
        != {subject.logical_subject}
    ):
        _fail("subject_stage_matrix_dataset_binding_mismatch")
    profile_sets: dict[str, set[str]] = {dataset_id: set() for dataset_id in dataset_ids}
    for plan_job, _ in selected:
        profile_sets[str(plan_job["dataset_id"])].add(str(plan_job["profile_id"]))
    if len({frozenset(value) for value in profile_sets.values()}) != 1:
        _fail("technical_replicates_use_different_profile_sets")

    forbidden = tuple(
        sorted(
            {
                *subject.forbidden_bindings,
                str(matrix_root),
                str(manifest_path),
                *dataset_ids,
            }
        )
    )
    rows_by_domain: dict[str, list[dict[str, Any]]] = {}
    gates: list[Path] = []
    receipts: list[Path] = []
    adapters: list[Path] = []
    job_boundaries: list[dict[str, Any]] = []
    for plan_job, runtime_job in selected:
        paths = _attempt_paths(matrix_root, int(plan_job["ordinal"]), int(runtime_job["attempt"]))
        adapter_path = paths["typed"] / "adapter-summary.json"
        adapter = _load_object(adapter_path, "adapter_summary")
        _validate_canonical_identity(adapter, "adapter_identity_sha256", "adapter_summary")
        if (
            adapter.get("schema_version") != 1
            or adapter.get("status") != "PASS"
            or adapter.get("profile_id") != plan_job["profile_id"]
            or adapter.get("technical_source_role") != plan_job["technical_source_role"]
            or adapter.get("private_paths_persisted") is not False
            or adapter.get("subject_identifier_persisted") is not False
            or adapter.get("source_rows_removed") != 0
            or not isinstance(adapter.get("files"), list)
        ):
            _fail("adapter_summary_boundary_invalid")
        job_has_rows = False
        for metadata in adapter["files"]:
            if not isinstance(metadata, dict) or metadata.get("domain") not in DOMAINS:
                _fail("adapter_file_inventory_invalid")
            domain = str(metadata["domain"])
            path = _safe_adapter_file(paths["typed"], metadata.get("path"))
            if path.stat().st_size != metadata.get("bytes") or sha256_file(path) != metadata.get(
                "sha256"
            ):
                _fail("adapter_file_digest_or_size_changed")
            table = pq.read_table(path)
            if table.schema != schema_for(domain):
                _fail("adapter_file_physical_schema_changed")
            rows = validate_rows(
                domain,
                table.to_pylist(),
                expected_reference_sha256=subject.reference_sha256,
                forbidden_bindings=forbidden,
            )
            if len(rows) != metadata.get("rows") or canonical_row_multiset(rows) != metadata.get(
                "canonical_row_multiset"
            ):
                _fail("adapter_file_row_boundary_changed")
            if rows:
                rows_by_domain.setdefault(domain, []).extend(rows)
                job_has_rows = True
        gate = paths["lifecycle"] / "input-gate.json"
        receipt = paths["lifecycle"] / "caller-receipt.json"
        if gate.is_symlink() or not gate.is_file() or receipt.is_symlink() or not receipt.is_file():
            _fail("matrix_job_gate_or_receipt_missing")
        gates.append(gate)
        receipts.append(receipt)
        adapters.append(adapter_path)
        job_boundaries.append(
            {
                "ordinal": int(plan_job["ordinal"]),
                "attempt": int(runtime_job["attempt"]),
                "profile_id": str(plan_job["profile_id"]),
                "adapter_summary_sha256": sha256_file(adapter_path),
                "input_gate_sha256": sha256_file(gate),
                "caller_receipt_sha256": sha256_file(receipt),
                "has_typed_rows": job_has_rows,
            }
        )
    if not rows_by_domain or not gates or not receipts or not adapters:
        _fail("subject_stage_has_no_receipt_bound_rows")
    private_boundary = {
        "schema_version": 1,
        "input_manifest_sha256": sha256_file(manifest_path),
        "matrix_plan_identity_sha256": plan["plan_identity_sha256"],
        "dataset_ids": sorted(dataset_ids),
        "reference_sha256": subject.reference_sha256,
        "jobs": sorted(job_boundaries, key=lambda row: row["ordinal"]),
    }
    return rows_by_domain, gates, receipts, adapters, private_boundary, forbidden


def _stage_matches_plan(stage: Path, plan: dict[str, Any]) -> dict[str, Any]:
    validation = validate_release_c_extension_stage(
        stage, expected_reference_sha256=str(plan["reference_sha256"])
    )
    manifest = _load_object(stage / "release_c_extensions.json", "extension_manifest")
    expected_receipts = {row["caller_receipt_sha256"] for row in plan["jobs"]}
    observed_receipts = {row.get("receipt_file_sha256") for row in manifest.get("caller_runs", [])}
    expected_adapters = {row["adapter_summary_sha256"] for row in plan["jobs"]}
    observed_adapters = {
        row.get("adapter_summary_file_sha256") for row in manifest.get("adapter_runs", [])
    }
    if (
        manifest.get("source_dataset_count") != len(plan["dataset_ids"])
        or observed_receipts != expected_receipts
        or observed_adapters != expected_adapters
    ):
        _fail("published_subject_stage_differs_from_private_plan")
    return validation


def build_subject_stage(
    manifest_path: str | Path,
    matrix_work_root: str | Path,
    dataset_ids: list[str],
    work_directory: str | Path,
    output_directory: str | Path,
    *,
    minimum_free_bytes: int,
) -> dict[str, Any]:
    """Build or validate one no-overwrite, one-subject portable extension stage."""

    manifest = Path(manifest_path).expanduser().resolve()
    matrix_root = Path(matrix_work_root).expanduser().resolve()
    work = Path(work_directory).expanduser().resolve()
    output = Path(output_directory).expanduser().resolve()
    if any(path.is_symlink() for path in (manifest, matrix_root, work, output)):
        _fail("subject_stage_top_level_path_must_not_be_symlink")
    _disk_gate(work, minimum_free_bytes)
    rows_by_domain, gates, receipts, adapters, private_plan, forbidden = _collect_subject_inputs(
        manifest, matrix_root, dataset_ids
    )
    private_plan.update(
        {
            "work_directory": str(work),
            "output_directory": str(output),
            "minimum_free_bytes": minimum_free_bytes,
        }
    )
    private_plan["plan_identity_sha256"] = _canonical_sha256(private_plan)

    plan_path = work / PLAN_NAME
    complete_path = work / COMPLETE_NAME
    if work.exists():
        existing = _load_object(plan_path, "subject_stage_plan")
        if existing != private_plan:
            _fail("subject_stage_private_plan_boundary_changed")
        if output.is_dir():
            validation = _stage_matches_plan(output, existing)
            completion = {
                "schema_version": 1,
                "status": "PASS",
                "plan_identity_sha256": existing["plan_identity_sha256"],
                "extension_set_identity_sha256": validation["extension_set_identity_sha256"],
            }
            if complete_path.exists():
                if _load_object(complete_path, "subject_stage_completion") != completion:
                    _fail("subject_stage_completion_boundary_changed")
                resume = "reused"
            else:
                _atomic_no_replace(complete_path, completion)
                resume = "finalized_after_interruption"
            return {
                "schema_version": 1,
                "status": "PASS",
                "resume": resume,
                "domains": validation["domains"],
                "rows": validation["rows"],
                "source_dataset_count": len(dataset_ids),
                "source_rows_removed": 0,
                "private_paths_persisted": False,
                "subject_identifier_persisted": False,
            }
        _fail("retained_incomplete_subject_stage_requires_new_work_directory")
    if output.exists():
        _fail("refusing_to_replace_existing_subject_stage")

    work.mkdir(parents=True)
    _atomic_no_replace(plan_path, private_plan)
    consolidated = work / "consolidated"
    consolidated.mkdir()
    concordance_summaries: dict[str, Any] = {}
    domain_paths: dict[str, Path] = {}
    source_rows = 0
    output_rows = 0
    for domain in sorted(rows_by_domain):
        _disk_gate(consolidated, minimum_free_bytes)
        rows = rows_by_domain[domain]
        source_rows += len(rows)
        annotated, concordance = _apply_concordance(
            domain,
            rows,
            reference_sha256=str(private_plan["reference_sha256"]),
            forbidden_bindings=(*forbidden, str(work), str(output)),
        )
        path = consolidated / f"{domain}.parquet"
        write_extension(
            domain,
            annotated,
            path,
            expected_reference_sha256=str(private_plan["reference_sha256"]),
            forbidden_bindings=(*forbidden, str(work), str(output)),
        )
        output_rows += len(annotated)
        domain_paths[domain] = path
        concordance_summaries[domain] = concordance
    if output_rows != source_rows:
        _fail("subject_stage_row_count_changed")
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "domains": len(domain_paths),
        "input_rows": source_rows,
        "output_rows": output_rows,
        "source_rows_removed": 0,
        "majority_vote_applied": False,
        "concordance": concordance_summaries,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
    }
    summary["summary_identity_sha256"] = _canonical_sha256(summary)
    _atomic_no_replace(work / SUMMARY_NAME, summary)

    _disk_gate(output, minimum_free_bytes)
    subject_binding = subject_binding_from_private_manifest(manifest, dataset_ids)
    staged = stage_release_c_extensions(
        output,
        domain_paths,
        gates,
        receipts,
        adapters,
        subject_binding,
    )
    validation = staged["validation"]
    if validation["rows"] != source_rows or validation["source_rows_removed"] != 0:
        _fail("published_subject_stage_row_boundary_changed")
    _stage_matches_plan(output, private_plan)
    completion = {
        "schema_version": 1,
        "status": "PASS",
        "plan_identity_sha256": private_plan["plan_identity_sha256"],
        "extension_set_identity_sha256": validation["extension_set_identity_sha256"],
    }
    _atomic_no_replace(complete_path, completion)
    return {
        "schema_version": 1,
        "status": "PASS",
        "resume": "created",
        "domains": validation["domains"],
        "rows": validation["rows"],
        "source_dataset_count": len(dataset_ids),
        "source_rows_removed": 0,
        "private_paths_persisted": False,
        "subject_identifier_persisted": False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="operation", required=True)
    build = subparsers.add_parser("build")
    build.add_argument("--manifest", type=Path, required=True)
    build.add_argument("--matrix-work-root", type=Path, required=True)
    build.add_argument("--dataset-id", action="append", default=[], required=True)
    build.add_argument("--work-directory", type=Path, required=True)
    build.add_argument("--output-directory", type=Path, required=True)
    build.add_argument("--minimum-free-bytes", type=int, required=True)
    validate = subparsers.add_parser("validate")
    validate.add_argument("--stage", type=Path, required=True)
    validate.add_argument("--reference-sha256")
    args = parser.parse_args(argv)
    try:
        if args.operation == "build":
            result = build_subject_stage(
                args.manifest,
                args.matrix_work_root,
                args.dataset_id,
                args.work_directory,
                args.output_directory,
                minimum_free_bytes=args.minimum_free_bytes,
            )
        else:
            result = validate_release_c_extension_stage(
                args.stage, expected_reference_sha256=args.reference_sha256
            )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, PortableEvidenceError, SubjectStageError):
        print("error: release_c_subject_stage_failed", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
