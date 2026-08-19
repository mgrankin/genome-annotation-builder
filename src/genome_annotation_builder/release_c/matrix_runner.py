"""Sequential, restart-safe execution of a generic Release-C caller matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .binding_builder import build_caller_bindings, validate_caller_bindings
from .caller_lifecycle import run_caller_lifecycle
from .extension_stage import ACCEPTED_PRIVATE_MANIFEST_STATES
from .portable_model import (
    FORBIDDEN_PORTABLE_TEXT,
    PortableEvidenceError,
    contains_forbidden_binding,
    sha256_file,
)
from .profile_adapter import (
    SUPPORTED_PROFILES,
    materialize_profile_adapter,
    validate_profile_adapter_output,
)
from .profile_catalog import build_profile_catalog, resolve_shipped_profile
from .resource_manifest import validate_resource_manifest_payload

MATRIX_SCHEMA_VERSION = 1
PLAN_NAME = "matrix-plan.private.json"
RUNTIME_NAME = "matrix-runtime.private.json"
SUMMARY_NAME = "matrix-summary.json"
DEFAULT_FIXED_PROFILES = (
    "manta-germline-1.6.0",
    "delly-short-read-2.5.1",
    "cnvpytor-combined-rd-baf-1.3.2",
    "hla-la-short-read-1.0.4",
    "spechla-short-read-bam-1.0.12",
    "gatk-mitochondria-normal-shifted-4.6.2.0",
    "readback-bounded-pysam-0.24.0",
)
DEFAULT_KARYOTYPE_PROFILES = {
    "XX": "expansionhunter-grch38-female-5.0.0",
    "XY": "expansionhunter-grch38-male-5.0.0",
}
HLA_PROFILES = frozenset(
    {
        "hla-la-short-read-1.0.4",
        "spechla-short-read-bam-1.0.12",
    }
)


class MatrixRunnerError(RuntimeError):
    """The private matrix cannot progress without violating a checkpoint."""


def _fail(message: str) -> None:
    raise MatrixRunnerError(message)


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
        _fail(f"{label}_missing")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise MatrixRunnerError(f"{label}_unreadable") from exc
    if not isinstance(value, dict):
        _fail(f"{label}_invalid")
    return value


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def _atomic_no_replace(path: Path, value: object) -> None:
    if path.exists():
        _fail("immutable_matrix_checkpoint_already_exists")
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _json_bytes(value)
    descriptor, raw_name = tempfile.mkstemp(prefix=f".{path.name}.partial.", dir=path.parent)
    temporary = Path(raw_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as exc:
            raise MatrixRunnerError("immutable_matrix_checkpoint_already_exists") from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _atomic_replace_private(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _json_bytes(value)
    descriptor, raw_name = tempfile.mkstemp(prefix=f".{path.name}.partial.", dir=path.parent)
    temporary = Path(raw_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _timestamp() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _load_matrix_config(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {
            "schema_version": MATRIX_SCHEMA_VERSION,
            "fixed_profiles": list(DEFAULT_FIXED_PROFILES),
            "karyotype_profiles": dict(DEFAULT_KARYOTYPE_PROFILES),
        }
    config = _load_object(path, "matrix_config")
    if set(config) != {"schema_version", "fixed_profiles", "karyotype_profiles"}:
        _fail("matrix_config_fields_invalid")
    return config


def _validate_matrix_config(config: dict[str, Any], catalog_ids: set[str]) -> None:
    fixed = config.get("fixed_profiles")
    karyotype = config.get("karyotype_profiles")
    if (
        config.get("schema_version") != MATRIX_SCHEMA_VERSION
        or not isinstance(fixed, list)
        or not fixed
        or any(not isinstance(value, str) for value in fixed)
        or len(fixed) != len(set(fixed))
        or not isinstance(karyotype, dict)
        or set(karyotype) != {"XX", "XY"}
        or any(not isinstance(value, str) for value in karyotype.values())
    ):
        _fail("matrix_config_invalid")
    selected = set(fixed) | set(karyotype.values())
    if not selected.issubset(catalog_ids) or not selected.issubset(SUPPORTED_PROFILES):
        _fail("matrix_config_profile_not_shipped_or_unsupported")
    if set(fixed) & set(karyotype.values()) or len(set(karyotype.values())) != 2:
        _fail("matrix_config_profile_roles_overlap")


def _validate_manifest_datasets(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    if (
        manifest.get("schema_version") != 1
        or manifest.get("status") not in ACCEPTED_PRIVATE_MANIFEST_STATES
    ):
        _fail("private_manifest_schema_or_status_invalid")
    datasets = manifest.get("datasets")
    if (
        not isinstance(datasets, list)
        or not datasets
        or any(not isinstance(row, dict) for row in datasets)
    ):
        _fail("private_manifest_datasets_invalid")
    dataset_ids: set[str] = set()
    by_subject: dict[str, list[dict[str, Any]]] = {}
    for row in datasets:
        dataset_id = row.get("dataset_id")
        subject = row.get("logical_subject")
        role = row.get("technical_role")
        karyotype = row.get("expected_karyotype")
        if (
            not isinstance(dataset_id, str)
            or not dataset_id
            or dataset_id in dataset_ids
            or not isinstance(subject, str)
            or not subject
            or role not in {"primary", "technical_replicate"}
            or karyotype not in {"XX", "XY"}
        ):
            _fail("private_manifest_dataset_identity_or_role_invalid")
        dataset_ids.add(dataset_id)
        by_subject.setdefault(subject, []).append(row)
    for rows in by_subject.values():
        if sum(row["technical_role"] == "primary" for row in rows) != 1:
            _fail("logical_subject_requires_exactly_one_primary_dataset")
        if len({row["expected_karyotype"] for row in rows}) != 1:
            _fail("technical_replicates_disagree_on_expected_karyotype")
    return datasets


def _validate_resource_manifest_boundary(
    resource_manifest: dict[str, Any], catalog: dict[str, Any]
) -> None:
    if resource_manifest.get("schema_version") != 1 or resource_manifest.get("status") not in {
        "PASS",
        "PASS_WITH_DOCUMENTED_LIMITATIONS",
    }:
        _fail("resource_manifest_schema_or_status_invalid")
    declared = resource_manifest.get("profile_catalog")
    if not isinstance(declared, dict) or declared != catalog:
        _fail("resource_manifest_profile_catalog_boundary_mismatch")
    embedded = dict(declared)
    claimed_identity = embedded.pop("catalog_identity_sha256", None)
    if claimed_identity != _canonical_sha256(embedded):
        _fail("resource_manifest_profile_catalog_identity_invalid")
    try:
        validate_resource_manifest_payload(resource_manifest)
    except PortableEvidenceError as exc:
        raise MatrixRunnerError("resource_manifest_release_boundary_invalid") from exc


def build_matrix_plan(
    manifest_path: str | Path,
    resource_bindings_path: str | Path,
    resource_manifest_path: str | Path,
    *,
    adapter_options_path: str | Path | None = None,
    matrix_config_path: str | Path | None = None,
    runtime: str,
    user: str,
    minimum_free_bytes: int,
    indexed_probe_bases: int,
    selected_dataset_ids: list[str] | None = None,
) -> dict[str, Any]:
    """Compile a deterministic private matrix without executing a caller."""

    manifest_file = Path(manifest_path)
    resources_file = Path(resource_bindings_path)
    resource_manifest_file = Path(resource_manifest_path)
    options_file = Path(adapter_options_path) if adapter_options_path else None
    matrix_file = Path(matrix_config_path) if matrix_config_path else None
    manifest = _load_object(manifest_file, "private_manifest")
    resources = _load_object(resources_file, "resource_bindings")
    resource_manifest = _load_object(resource_manifest_file, "resource_manifest")
    options = _load_object(options_file, "adapter_options") if options_file else None
    config = _load_matrix_config(matrix_file)
    datasets = _validate_manifest_datasets(manifest)
    catalog = build_profile_catalog()
    catalog_ids = {str(row["profile_id"]) for row in catalog["profiles"]}
    _validate_matrix_config(config, catalog_ids)
    if runtime not in {"docker", "podman"} or re.fullmatch(r"[0-9]+:[0-9]+", user) is None:
        _fail("matrix_runtime_or_user_invalid")
    if minimum_free_bytes < 1 or indexed_probe_bases < 1:
        _fail("matrix_safety_policy_invalid")
    if resources.get("schema_version") != 1 or not isinstance(resources.get("profiles"), dict):
        _fail("resource_bindings_schema_invalid")
    _validate_resource_manifest_boundary(resource_manifest, catalog)
    selected_profiles = set(config["fixed_profiles"]) | set(config["karyotype_profiles"].values())
    if not selected_profiles.issubset(resources["profiles"]):
        _fail("resource_bindings_profile_coverage_incomplete")
    if selected_profiles & HLA_PROFILES and (
        not isinstance(options, dict)
        or options.get("schema_version") != 1
        or not HLA_PROFILES.issubset(set(options.get("profiles", {})))
    ):
        _fail("matrix_hla_adapter_options_incomplete")
    requested = selected_dataset_ids or [str(row["dataset_id"]) for row in datasets]
    if not requested or len(requested) != len(set(requested)):
        _fail("selected_dataset_ids_invalid")
    dataset_by_id = {str(row["dataset_id"]): row for row in datasets}
    if set(requested) - set(dataset_by_id):
        _fail("selected_dataset_id_not_in_manifest")
    jobs: list[dict[str, Any]] = []
    ordinal = 0
    for dataset_id in requested:
        dataset = dataset_by_id[dataset_id]
        profiles = [
            *config["fixed_profiles"],
            config["karyotype_profiles"][dataset["expected_karyotype"]],
        ]
        for profile_id in profiles:
            ordinal += 1
            jobs.append(
                {
                    "ordinal": ordinal,
                    "dataset_id": dataset_id,
                    "logical_subject": dataset["logical_subject"],
                    "technical_source_role": dataset["technical_role"],
                    "expected_karyotype": dataset["expected_karyotype"],
                    "profile_id": profile_id,
                }
            )
    plan: dict[str, Any] = {
        "schema_version": MATRIX_SCHEMA_VERSION,
        "status": "READY",
        "input_identities": {
            "private_manifest_sha256": sha256_file(manifest_file),
            "resource_bindings_sha256": sha256_file(resources_file),
            "resource_manifest_sha256": sha256_file(resource_manifest_file),
            "adapter_options_sha256": sha256_file(options_file) if options_file else None,
            "matrix_config_identity_sha256": _canonical_sha256(config),
            "profile_catalog_identity_sha256": catalog["catalog_identity_sha256"],
            "resource_manifest_status": resource_manifest.get("status"),
        },
        "execution": {
            "runtime": runtime,
            "user": user,
            "minimum_free_bytes": minimum_free_bytes,
            "indexed_probe_bases": indexed_probe_bases,
            "sequential_jobs": True,
        },
        "dataset_count": len(requested),
        "logical_subject_count": len(
            {dataset_by_id[value]["logical_subject"] for value in requested}
        ),
        "job_count": len(jobs),
        "jobs": jobs,
    }
    plan["plan_identity_sha256"] = _canonical_sha256(plan)
    return plan


def _initial_runtime(plan: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": MATRIX_SCHEMA_VERSION,
        "state": "ready",
        "plan_identity_sha256": plan["plan_identity_sha256"],
        "active_job_ordinal": None,
        "completed_jobs": 0,
        "total_jobs": plan["job_count"],
        "updated_at": _timestamp(),
        "jobs": [
            {
                "ordinal": job["ordinal"],
                "dataset_id": job["dataset_id"],
                "profile_id": job["profile_id"],
                "state": "pending",
                "attempt": 1,
                "last_completed_stage": None,
            }
            for job in plan["jobs"]
        ],
    }


def _validate_runtime(runtime: dict[str, Any], plan: dict[str, Any]) -> None:
    jobs = runtime.get("jobs")
    if (
        runtime.get("schema_version") != MATRIX_SCHEMA_VERSION
        or runtime.get("plan_identity_sha256") != plan["plan_identity_sha256"]
        or runtime.get("total_jobs") != plan["job_count"]
        or not isinstance(jobs, list)
        or len(jobs) != plan["job_count"]
    ):
        _fail("matrix_runtime_plan_boundary_changed")
    for expected, observed in zip(plan["jobs"], jobs, strict=True):
        if (
            not isinstance(observed, dict)
            or observed.get("ordinal") != expected["ordinal"]
            or observed.get("dataset_id") != expected["dataset_id"]
            or observed.get("profile_id") != expected["profile_id"]
            or observed.get("state") not in {"pending", "running", "failed", "complete"}
            or not isinstance(observed.get("attempt"), int)
            or observed["attempt"] < 1
        ):
            _fail("matrix_runtime_job_boundary_changed")


def _save_runtime(path: Path, runtime: dict[str, Any]) -> None:
    runtime["updated_at"] = _timestamp()
    runtime["completed_jobs"] = sum(row["state"] == "complete" for row in runtime["jobs"])
    _atomic_replace_private(path, runtime)


def _free_space_gate(root: Path, minimum_free_bytes: int) -> int:
    free = shutil.disk_usage(root).free
    if free < minimum_free_bytes:
        _fail("matrix_free_space_floor_reached")
    return free


def _attempt_paths(work_root: Path, ordinal: int, attempt: int) -> dict[str, Path]:
    root = work_root / f"job-{ordinal:04d}" / f"attempt-{attempt:04d}"
    return {
        "root": root,
        "bindings": root / "bindings.private.json",
        "caller_workspace": root / "caller-workspace",
        "lifecycle": root / "lifecycle",
        "typed": root / "typed",
    }


def _checkpoint_complete(paths: dict[str, Path]) -> bool:
    lifecycle = paths["lifecycle"]
    expected = {
        lifecycle / "input-gate.json",
        lifecycle / "plan",
        lifecycle / "execution-report.json",
        lifecycle / "caller-receipt.json",
        lifecycle / "lifecycle-summary.json",
    }
    return (
        paths["bindings"].is_file()
        and all(path.exists() for path in expected)
        and (paths["typed"] / "adapter-summary.json").is_file()
    )


def _portable_summary(
    plan: dict[str, Any], work_root: Path, runtime: dict[str, Any]
) -> dict[str, Any]:
    profile_counts: Counter[str] = Counter()
    domain_files: Counter[str] = Counter()
    domain_rows: Counter[str] = Counter()
    caller_receipts: set[str] = set()
    adapter_identities: set[str] = set()
    for plan_job, runtime_job in zip(plan["jobs"], runtime["jobs"], strict=True):
        paths = _attempt_paths(work_root, plan_job["ordinal"], runtime_job["attempt"])
        adapter = _load_object(paths["typed"] / "adapter-summary.json", "adapter_summary")
        lifecycle = _load_object(paths["lifecycle"] / "lifecycle-summary.json", "lifecycle_summary")
        profile_counts[str(plan_job["profile_id"])] += 1
        caller_receipts.add(str(lifecycle["receipt_identity_sha256"]))
        adapter_identities.add(str(adapter["adapter_identity_sha256"]))
        for row in adapter.get("files", []):
            if not isinstance(row, dict):
                _fail("adapter_summary_file_inventory_invalid")
            domain_files[str(row["domain"])] += 1
            domain_rows[str(row["domain"])] += int(row["rows"])
    summary: dict[str, Any] = {
        "schema_version": MATRIX_SCHEMA_VERSION,
        "status": "PASS",
        "execution_model": "strictly_sequential",
        "dataset_count": plan["dataset_count"],
        "logical_subject_count": plan["logical_subject_count"],
        "job_count": plan["job_count"],
        "completed_jobs": plan["job_count"],
        "profile_run_counts": dict(sorted(profile_counts.items())),
        "domain_file_counts": dict(sorted(domain_files.items())),
        "domain_row_counts": dict(sorted(domain_rows.items())),
        "unique_caller_receipts": len(caller_receipts),
        "unique_adapter_identities": len(adapter_identities),
        "source_rows_removed": 0,
        "caller_or_replicate_rows_collapsed": 0,
        "logical_subject_identifiers_persisted": False,
        "dataset_identifiers_persisted": False,
        "private_paths_persisted": False,
        "private_manifest_digest_persisted": False,
    }
    summary["matrix_summary_identity_sha256"] = _canonical_sha256(summary)
    return summary


def _forbidden_tokens(plan: dict[str, Any], paths: list[Path]) -> tuple[str, ...]:
    values = {str(path) for path in paths}
    for job in plan["jobs"]:
        values.update((str(job["dataset_id"]), str(job["logical_subject"])))
    return tuple(sorted(value for value in values if value))


def _assert_portable_summary(summary: dict[str, Any], forbidden: tuple[str, ...]) -> None:
    rendered = json.dumps(summary, sort_keys=True, ensure_ascii=True)
    if FORBIDDEN_PORTABLE_TEXT.search(rendered) or any(
        contains_forbidden_binding(rendered, value) for value in forbidden
    ):
        _fail("portable_matrix_summary_contains_private_binding")


def run_release_c_matrix(
    manifest_path: str | Path,
    resource_bindings_path: str | Path,
    resource_manifest_path: str | Path,
    work_root: str | Path,
    *,
    adapter_options_path: str | Path | None = None,
    matrix_config_path: str | Path | None = None,
    forbidden_binding_file: str | Path | None = None,
    runtime: str = "docker",
    user: str,
    minimum_free_bytes: int,
    indexed_probe_bases: int,
    selected_dataset_ids: list[str] | None = None,
    retry_failed: bool = False,
    validate_only: bool = False,
) -> dict[str, Any]:
    """Run or validate every selected dataset/profile pair in strict sequence."""

    manifest = Path(manifest_path).expanduser().resolve()
    resource_bindings = Path(resource_bindings_path).expanduser().resolve()
    resource_manifest = Path(resource_manifest_path).expanduser().resolve()
    options = Path(adapter_options_path).expanduser().resolve() if adapter_options_path else None
    matrix_config = Path(matrix_config_path).expanduser().resolve() if matrix_config_path else None
    forbidden_file = (
        Path(forbidden_binding_file).expanduser().resolve() if forbidden_binding_file else None
    )
    root = Path(work_root).expanduser().resolve()
    if root.is_symlink():
        _fail("matrix_work_root_must_not_be_symlink")
    root.mkdir(parents=True, exist_ok=True)
    plan = build_matrix_plan(
        manifest,
        resource_bindings,
        resource_manifest,
        adapter_options_path=options,
        matrix_config_path=matrix_config,
        runtime=runtime,
        user=user,
        minimum_free_bytes=minimum_free_bytes,
        indexed_probe_bases=indexed_probe_bases,
        selected_dataset_ids=selected_dataset_ids,
    )
    plan_path = root / PLAN_NAME
    runtime_path = root / RUNTIME_NAME
    summary_path = root / SUMMARY_NAME
    if plan_path.exists():
        if _load_object(plan_path, "matrix_plan") != plan:
            _fail("existing_matrix_plan_differs")
    elif validate_only:
        _fail("validate_only_requires_existing_matrix_plan")
    else:
        _atomic_no_replace(plan_path, plan)
    if runtime_path.exists():
        runtime_state = _load_object(runtime_path, "matrix_runtime")
        _validate_runtime(runtime_state, plan)
    elif validate_only:
        _fail("validate_only_requires_existing_matrix_runtime")
    else:
        runtime_state = _initial_runtime(plan)
        _atomic_no_replace(runtime_path, runtime_state)

    private_paths = [
        manifest,
        resource_bindings,
        resource_manifest,
        root,
        *([options] if options else []),
        *([matrix_config] if matrix_config else []),
        *([forbidden_file] if forbidden_file else []),
    ]
    for state in runtime_state["jobs"]:
        if state["state"] != "failed":
            continue
        if validate_only or not retry_failed:
            _fail("matrix_has_retained_failed_attempt; use --retry-failed after review")
        adapter_only_retry = state.get("last_completed_stage") == "caller_lifecycle"
        state.update(
            {
                "state": "pending",
                "attempt": state["attempt"] + (0 if adapter_only_retry else 1),
                "last_completed_stage": ("caller_lifecycle" if adapter_only_retry else None),
                "error_type": None,
                "error_message": None,
            }
        )
    _save_runtime(runtime_path, runtime_state)

    # Compile and verify every binding before the first expensive caller starts.
    # This catches a missing late-matrix resource without wasting earlier calls.
    for plan_job, state in zip(plan["jobs"], runtime_state["jobs"], strict=True):
        paths = _attempt_paths(root, plan_job["ordinal"], state["attempt"])
        if validate_only and not _checkpoint_complete(paths):
            _fail("validate_only_found_incomplete_job_checkpoint")
        runtime_state["state"] = "preflight"
        runtime_state["active_job_ordinal"] = plan_job["ordinal"]
        _save_runtime(runtime_path, runtime_state)
        try:
            _free_space_gate(root, minimum_free_bytes)
            profile_path = resolve_shipped_profile(str(plan_job["profile_id"]))
            if paths["bindings"].exists():
                validate_caller_bindings(
                    manifest,
                    str(plan_job["dataset_id"]),
                    profile_path,
                    resource_bindings,
                    resource_manifest,
                    paths["caller_workspace"],
                    paths["bindings"],
                    runtime=runtime,
                    user=user,
                    minimum_free_bytes=minimum_free_bytes,
                )
            elif validate_only:
                _fail("validate_only_missing_binding_checkpoint")
            else:
                build_caller_bindings(
                    manifest,
                    str(plan_job["dataset_id"]),
                    profile_path,
                    resource_bindings,
                    resource_manifest,
                    paths["caller_workspace"],
                    paths["bindings"],
                    runtime=runtime,
                    user=user,
                    minimum_free_bytes=minimum_free_bytes,
                )
            if state["state"] != "complete":
                state["state"] = "pending"
                state["last_completed_stage"] = "bindings"
            _save_runtime(runtime_path, runtime_state)
        except BaseException as exc:
            state.update(
                {
                    "state": "failed",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                }
            )
            runtime_state["state"] = "failed"
            _save_runtime(runtime_path, runtime_state)
            raise MatrixRunnerError(
                f"matrix_preflight_job_{plan_job['ordinal']:04d}_failed; inspect private runtime"
            ) from exc

    for plan_job, state in zip(plan["jobs"], runtime_state["jobs"], strict=True):
        paths = _attempt_paths(root, plan_job["ordinal"], state["attempt"])
        state["state"] = "running"
        runtime_state["state"] = "running"
        runtime_state["active_job_ordinal"] = plan_job["ordinal"]
        _save_runtime(runtime_path, runtime_state)
        try:
            _free_space_gate(root, minimum_free_bytes)
            lifecycle_report = run_caller_lifecycle(
                manifest,
                str(plan_job["dataset_id"]),
                paths["bindings"],
                resource_manifest,
                paths["lifecycle"],
                profile_id=str(plan_job["profile_id"]),
                minimum_free_bytes=minimum_free_bytes,
                indexed_probe_bases=indexed_probe_bases,
                forbidden_binding_file=forbidden_file,
            )
            state["last_completed_stage"] = "caller_lifecycle"
            state["lifecycle_identity_sha256"] = lifecycle_report["lifecycle_identity_sha256"]
            _save_runtime(runtime_path, runtime_state)

            _free_space_gate(root, minimum_free_bytes)
            adapter_keywords = {
                "profile_id": str(plan_job["profile_id"]),
                "adapter_options_path": options,
                "forbidden_binding_file": forbidden_file,
            }
            if validate_only:
                adapter_report = validate_profile_adapter_output(
                    paths["typed"],
                    paths["lifecycle"],
                    paths["bindings"],
                    str(plan_job["technical_source_role"]),
                    **adapter_keywords,
                )
            else:
                adapter_report = materialize_profile_adapter(
                    paths["lifecycle"],
                    paths["bindings"],
                    str(plan_job["technical_source_role"]),
                    paths["typed"],
                    **adapter_keywords,
                )
            state.update(
                {
                    "state": "complete",
                    "last_completed_stage": "profile_adapter",
                    "adapter_identity_sha256": adapter_report["adapter_identity_sha256"],
                    "typed_files": adapter_report["typed_files"],
                    "typed_rows": adapter_report["typed_rows"],
                    "error_type": None,
                    "error_message": None,
                }
            )
            _save_runtime(runtime_path, runtime_state)
        except BaseException as exc:
            state.update(
                {
                    "state": "failed",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                }
            )
            runtime_state["state"] = "failed"
            runtime_state["active_job_ordinal"] = plan_job["ordinal"]
            _save_runtime(runtime_path, runtime_state)
            raise MatrixRunnerError(
                f"matrix_job_{plan_job['ordinal']:04d}_failed; inspect private runtime"
            ) from exc

    runtime_state["state"] = "complete"
    runtime_state["active_job_ordinal"] = None
    _save_runtime(runtime_path, runtime_state)
    summary = _portable_summary(plan, root, runtime_state)
    forbidden = _forbidden_tokens(plan, private_paths)
    _assert_portable_summary(summary, forbidden)
    if summary_path.exists():
        if _load_object(summary_path, "matrix_summary") != summary:
            _fail("existing_matrix_summary_changed")
    elif validate_only:
        _fail("validate_only_requires_existing_matrix_summary")
    else:
        _atomic_no_replace(summary_path, summary)
    return {
        **summary,
        "resume": {
            "validated_only": validate_only,
            "retained_attempts": sum(row["attempt"] for row in runtime_state["jobs"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Sequential generic Release-C caller matrix")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--resource-bindings", type=Path, required=True)
    parser.add_argument("--resource-manifest", type=Path, required=True)
    parser.add_argument("--adapter-options", type=Path)
    parser.add_argument("--matrix-config", type=Path)
    parser.add_argument("--forbidden-bindings-file", type=Path)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--runtime", choices=("docker", "podman"), default="docker")
    parser.add_argument("--user", required=True)
    parser.add_argument("--minimum-free-bytes", type=int, required=True)
    parser.add_argument("--indexed-probe-bases", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--dataset-id", action="append")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    try:
        result = run_release_c_matrix(
            args.manifest,
            args.resource_bindings,
            args.resource_manifest,
            args.work_root,
            adapter_options_path=args.adapter_options,
            matrix_config_path=args.matrix_config,
            forbidden_binding_file=args.forbidden_bindings_file,
            runtime=args.runtime,
            user=args.user,
            minimum_free_bytes=args.minimum_free_bytes,
            indexed_probe_bases=args.indexed_probe_bases,
            selected_dataset_ids=args.dataset_id,
            retry_failed=args.retry_failed,
            validate_only=args.validate_only,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, PortableEvidenceError, MatrixRunnerError):
        print(
            "Release-C matrix failed; inspect the authorized private runtime",
            file=os.sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
