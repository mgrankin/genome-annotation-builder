from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import matrix_runner
from genome_annotation_builder.release_c.matrix_runner import (
    MatrixRunnerError,
    build_matrix_plan,
    run_release_c_matrix,
)
from genome_annotation_builder.release_c.resource_manifest import build_resource_manifest


def _write(path: Path, value: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _inputs(tmp_path: Path, *, two_datasets: bool = True):
    datasets = [
        {
            "dataset_id": "PRIVATE_DATASET_ALPHA",
            "logical_subject": "PRIVATE_SUBJECT_ALPHA",
            "technical_role": "primary",
            "expected_karyotype": "XX",
        }
    ]
    if two_datasets:
        datasets.append(
            {
                "dataset_id": "PRIVATE_DATASET_BETA",
                "logical_subject": "PRIVATE_SUBJECT_BETA",
                "technical_role": "primary",
                "expected_karyotype": "XY",
            }
        )
    manifest = _write(
        tmp_path / "manifest.json",
        {"schema_version": 1, "status": "PASS", "datasets": datasets},
    )
    profile_ids = {
        *matrix_runner.DEFAULT_FIXED_PROFILES,
        *matrix_runner.DEFAULT_KARYOTYPE_PROFILES.values(),
    }
    resources = _write(
        tmp_path / "resource-bindings.json",
        {
            "schema_version": 1,
            "profiles": {profile_id: {"roles": {}} for profile_id in profile_ids},
            "dataset_profiles": {},
        },
    )
    freeze = _write(
        tmp_path / "resource-manifest.json",
        build_resource_manifest(
            {
                "schema_version": 1,
                "status": "PASS",
                "scope": "synthetic matrix-runner control",
                "production_dataset_execution_claimed": False,
                "tools": {"synthetic": {"status": "PASS"}},
                "resources": {"synthetic": {"status": "PASS"}},
                "runtime_resource_trees": {},
                "validation_evidence": {"synthetic": {"status": "PASS"}},
                "warnings": [],
            }
        ),
    )
    options = _write(
        tmp_path / "adapter-options.json",
        {
            "schema_version": 1,
            "profiles": {
                profile_id: {
                    "caller_database_name": "synthetic",
                    "caller_database_version": "1",
                    "caller_database_identity_sha256": "1" * 64,
                    "nomenclature_database_name": "synthetic",
                    "nomenclature_version": "1",
                    "nomenclature_identity_sha256": "2" * 64,
                    "database_age_limitation": False,
                }
                for profile_id in matrix_runner.HLA_PROFILES
            },
        },
    )
    return manifest, resources, freeze, options


def _small_matrix(tmp_path: Path) -> Path:
    return _write(
        tmp_path / "small-matrix.json",
        {
            "schema_version": 1,
            "fixed_profiles": ["manta-germline-1.6.0"],
            "karyotype_profiles": {
                "XX": "expansionhunter-grch38-female-5.0.0",
                "XY": "expansionhunter-grch38-male-5.0.0",
            },
        },
    )


def test_matrix_plan_selects_exact_karyotype_profile(tmp_path: Path) -> None:
    manifest, resources, freeze, options = _inputs(tmp_path)
    plan = build_matrix_plan(
        manifest,
        resources,
        freeze,
        adapter_options_path=options,
        runtime="docker",
        user="1000:1000",
        minimum_free_bytes=1,
        indexed_probe_bases=1,
    )

    assert plan["dataset_count"] == 2
    assert plan["logical_subject_count"] == 2
    assert plan["job_count"] == 16
    alpha = [row["profile_id"] for row in plan["jobs"] if row["dataset_id"].endswith("ALPHA")]
    beta = [row["profile_id"] for row in plan["jobs"] if row["dataset_id"].endswith("BETA")]
    assert "expansionhunter-grch38-female-5.0.0" in alpha
    assert "expansionhunter-grch38-male-5.0.0" not in alpha
    assert "expansionhunter-grch38-male-5.0.0" in beta
    assert "expansionhunter-grch38-female-5.0.0" not in beta


def test_matrix_plan_rejects_resource_manifest_for_another_profile_release(
    tmp_path: Path,
) -> None:
    manifest, resources, freeze, options = _inputs(tmp_path)
    payload = json.loads(freeze.read_text(encoding="utf-8"))
    payload["profile_catalog"]["profiles"][0]["profile_identity_sha256"] = "0" * 64
    _write(freeze, payload)
    with pytest.raises(
        MatrixRunnerError, match="resource_manifest_profile_catalog_boundary_mismatch"
    ):
        build_matrix_plan(
            manifest,
            resources,
            freeze,
            adapter_options_path=options,
            runtime="docker",
            user="1000:1000",
            minimum_free_bytes=1,
            indexed_probe_bases=1,
        )


def test_matrix_plan_rejects_legacy_resource_manifest_without_identity(
    tmp_path: Path,
) -> None:
    manifest, resources, freeze, options = _inputs(tmp_path)
    payload = json.loads(freeze.read_text(encoding="utf-8"))
    payload.pop("manifest_identity_sha256")
    payload.pop("private_paths_persisted")
    _write(freeze, payload)
    with pytest.raises(MatrixRunnerError, match="resource_manifest_release_boundary_invalid"):
        build_matrix_plan(
            manifest,
            resources,
            freeze,
            adapter_options_path=options,
            runtime="docker",
            user="1000:1000",
            minimum_free_bytes=1,
            indexed_probe_bases=1,
        )


def _install_fake_execution(
    monkeypatch: pytest.MonkeyPatch,
    *,
    fail_once: bool = False,
    fail_adapter_once: bool = False,
    fail_binding_ordinal: int | None = None,
):
    calls = {"build": 0, "lifecycle": 0, "materialize": 0, "validate": 0}
    failed = {"value": False}
    adapter_failed = {"value": False}

    def fake_profile(profile_id: str) -> Path:
        return Path(f"/synthetic/{profile_id}.json")

    def fake_build(*args, **kwargs):
        calls["build"] += 1
        if fail_binding_ordinal == calls["build"]:
            raise RuntimeError("synthetic missing late-matrix resource")
        output = Path(args[6])
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text('{"schema_version":1}\n', encoding="utf-8")
        Path(args[5]).mkdir(parents=True, exist_ok=True)
        return {"status": "PASS"}

    def fake_validate_bindings(*args, **kwargs):
        assert Path(args[6]).is_file()
        return {"status": "PASS"}

    def fake_lifecycle(*args, **kwargs):
        calls["lifecycle"] += 1
        if fail_once and not failed["value"]:
            failed["value"] = True
            raise RuntimeError("synthetic retained failure")
        root = Path(args[4])
        (root / "plan").mkdir(parents=True, exist_ok=True)
        for name in (
            "input-gate.json",
            "execution-report.json",
            "caller-receipt.json",
        ):
            (root / name).write_text("{}\n", encoding="utf-8")
        profile = str(kwargs["profile_id"])
        identity = ("a" if profile.startswith("manta") else "b") * 64
        _write(
            root / "lifecycle-summary.json",
            {
                "lifecycle_identity_sha256": identity,
                "receipt_identity_sha256": identity,
            },
        )
        return {"lifecycle_identity_sha256": identity}

    def fake_materialize(*args, **kwargs):
        calls["materialize"] += 1
        if fail_adapter_once and not adapter_failed["value"]:
            adapter_failed["value"] = True
            raise RuntimeError("synthetic adapter failure")
        output = Path(args[3])
        output.mkdir(parents=True, exist_ok=True)
        profile = str(kwargs["profile_id"])
        domain = "structural_variants" if profile.startswith("manta") else "repeat_expansions"
        identity = ("c" if profile.startswith("manta") else "d") * 64
        _write(
            output / "adapter-summary.json",
            {
                "adapter_identity_sha256": identity,
                "files": [{"domain": domain, "path": "typed.parquet", "rows": 3}],
            },
        )
        return {
            "adapter_identity_sha256": identity,
            "typed_files": 1,
            "typed_rows": 3,
        }

    def fake_validate_adapter(*args, **kwargs):
        calls["validate"] += 1
        summary = json.loads((Path(args[0]) / "adapter-summary.json").read_text())
        return {
            "adapter_identity_sha256": summary["adapter_identity_sha256"],
            "typed_files": len(summary["files"]),
            "typed_rows": sum(row["rows"] for row in summary["files"]),
        }

    monkeypatch.setattr(matrix_runner, "resolve_shipped_profile", fake_profile)
    monkeypatch.setattr(matrix_runner, "build_caller_bindings", fake_build)
    monkeypatch.setattr(matrix_runner, "validate_caller_bindings", fake_validate_bindings)
    monkeypatch.setattr(matrix_runner, "run_caller_lifecycle", fake_lifecycle)
    monkeypatch.setattr(matrix_runner, "materialize_profile_adapter", fake_materialize)
    monkeypatch.setattr(matrix_runner, "validate_profile_adapter_output", fake_validate_adapter)
    return calls


def test_matrix_run_and_validate_only_replay_are_identifier_free(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, resources, freeze, _options = _inputs(tmp_path, two_datasets=False)
    matrix = _small_matrix(tmp_path)
    work = tmp_path / "private-work"
    calls = _install_fake_execution(monkeypatch)
    keywords = {
        "matrix_config_path": matrix,
        "runtime": "docker",
        "user": "1000:1000",
        "minimum_free_bytes": 1,
        "indexed_probe_bases": 1,
    }

    first = run_release_c_matrix(manifest, resources, freeze, work, **keywords)
    assert first["status"] == "PASS"
    assert first["job_count"] == 2
    assert first["completed_jobs"] == 2
    assert first["source_rows_removed"] == 0
    assert first["caller_or_replicate_rows_collapsed"] == 0
    assert calls["build"] == 2
    assert calls["materialize"] == 2
    portable = (work / matrix_runner.SUMMARY_NAME).read_text(encoding="utf-8")
    assert "PRIVATE_DATASET_ALPHA" not in portable
    assert "PRIVATE_SUBJECT_ALPHA" not in portable
    assert str(tmp_path) not in portable

    replay = run_release_c_matrix(
        manifest,
        resources,
        freeze,
        work,
        validate_only=True,
        **keywords,
    )
    assert replay["status"] == "PASS"
    assert replay["resume"]["validated_only"] is True
    assert calls["validate"] == 2


def test_failed_attempt_is_retained_and_retry_uses_new_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, resources, freeze, _options = _inputs(tmp_path, two_datasets=False)
    matrix = _small_matrix(tmp_path)
    work = tmp_path / "private-work"
    _install_fake_execution(monkeypatch, fail_once=True)
    keywords = {
        "matrix_config_path": matrix,
        "runtime": "docker",
        "user": "1000:1000",
        "minimum_free_bytes": 1,
        "indexed_probe_bases": 1,
    }

    with pytest.raises(MatrixRunnerError, match="matrix_job_0001_failed"):
        run_release_c_matrix(manifest, resources, freeze, work, **keywords)
    first_attempt = work / "job-0001" / "attempt-0001"
    assert first_attempt.is_dir()
    with pytest.raises(MatrixRunnerError, match="retained_failed_attempt"):
        run_release_c_matrix(manifest, resources, freeze, work, **keywords)

    result = run_release_c_matrix(
        manifest,
        resources,
        freeze,
        work,
        retry_failed=True,
        **keywords,
    )
    assert result["status"] == "PASS"
    assert first_attempt.is_dir()
    assert (work / "job-0001" / "attempt-0002").is_dir()


def test_adapter_failure_reuses_completed_caller_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, resources, freeze, _options = _inputs(tmp_path, two_datasets=False)
    matrix = _small_matrix(tmp_path)
    work = tmp_path / "private-work"
    calls = _install_fake_execution(monkeypatch, fail_adapter_once=True)
    keywords = {
        "matrix_config_path": matrix,
        "runtime": "docker",
        "user": "1000:1000",
        "minimum_free_bytes": 1,
        "indexed_probe_bases": 1,
    }

    with pytest.raises(MatrixRunnerError, match="matrix_job_0001_failed"):
        run_release_c_matrix(manifest, resources, freeze, work, **keywords)
    runtime = json.loads((work / matrix_runner.RUNTIME_NAME).read_text(encoding="utf-8"))
    assert runtime["jobs"][0]["last_completed_stage"] == "caller_lifecycle"

    result = run_release_c_matrix(
        manifest,
        resources,
        freeze,
        work,
        retry_failed=True,
        **keywords,
    )
    assert result["status"] == "PASS"
    assert (work / "job-0001" / "attempt-0001" / "typed").is_dir()
    assert not (work / "job-0001" / "attempt-0002").exists()
    assert calls["materialize"] == 3


def test_all_bindings_preflight_before_first_caller(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, resources, freeze, _options = _inputs(tmp_path, two_datasets=False)
    matrix = _small_matrix(tmp_path)
    work = tmp_path / "private-work"
    calls = _install_fake_execution(monkeypatch, fail_binding_ordinal=2)

    with pytest.raises(MatrixRunnerError, match="matrix_preflight_job_0002_failed"):
        run_release_c_matrix(
            manifest,
            resources,
            freeze,
            work,
            matrix_config_path=matrix,
            runtime="docker",
            user="1000:1000",
            minimum_free_bytes=1,
            indexed_probe_bases=1,
        )
    assert calls["build"] == 2
    assert calls["lifecycle"] == 0


def test_default_matrix_uses_the_resource_catalog_readback_version(tmp_path: Path) -> None:
    manifest, resources, freeze, options = _inputs(tmp_path)
    current_freeze = json.loads(freeze.read_bytes())
    metadata = {
        k: v
        for k, v in current_freeze.items()
        if k not in {"profile_catalog", "private_paths_persisted", "manifest_identity_sha256"}
    }
    bindings = json.loads(resources.read_bytes())
    legacy = matrix_runner.LEGACY_READBACK
    current = matrix_runner.CURRENT_READBACK
    bindings["profiles"][legacy] = bindings["profiles"].pop(current)
    _write(resources, bindings)
    _write(freeze, build_resource_manifest(metadata, catalog_version=1))
    plan = build_matrix_plan(
        manifest,
        resources,
        freeze,
        adapter_options_path=options,
        runtime="docker",
        user="1000:1000",
        minimum_free_bytes=1,
        indexed_probe_bases=1,
    )
    ids = [job["profile_id"] for job in plan["jobs"]]
    assert len(ids) == 16 and ids.count(legacy) == 2 and current not in ids
    _write(freeze, current_freeze)
    with pytest.raises(MatrixRunnerError, match="resource_bindings_profile_coverage_incomplete"):
        build_matrix_plan(
            manifest,
            resources,
            freeze,
            adapter_options_path=options,
            runtime="docker",
            user="1000:1000",
            minimum_free_bytes=1,
            indexed_probe_bases=1,
        )
