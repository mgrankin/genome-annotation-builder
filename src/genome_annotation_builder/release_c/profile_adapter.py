"""Restart-safe materialization of one shipped caller profile into typed evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq
from pyarrow import ArrowInvalid

from .caller_profile import load_caller_profile
from .caller_receipt import load_caller_run_receipt
from .cnvpytor_combined_adapter import (
    CnvpytorCombinedProvenance,
    materialize_cnvpytor_combined_calls,
)
from .delly_adapter import DellyProvenance, materialize_delly_calls
from .expansionhunter_adapter import (
    ExpansionHunterProvenance,
    materialize_expansionhunter_calls,
)
from .hla_la_adapter import HlaLaProvenance, materialize_hla_la_bestguess
from .manta_adapter import MantaProvenance, materialize_manta_germline_calls
from .mt_adapter import MitochondrialProvenance, materialize_filtered_mutect_vcf
from .portable_model import (
    SOURCE_ROLES,
    canonical_row_multiset,
    contains_forbidden_binding,
    schema_for,
    sha256_file,
    validate_rows,
)
from .profile_catalog import CURRENT_READBACK, LEGACY_READBACK, resolve_shipped_profile
from .readback_adapter import ReadbackProvenance, materialize_readback_capture
from .spechla_adapter import SpecHLAProvenance, materialize_spechla_results

ADAPTER_SCHEMA_VERSION = 1
SUMMARY_NAME = "adapter-summary.json"
SUPPORTED_PROFILES = frozenset(
    {
        "manta-germline-1.6.0",
        "delly-short-read-2.5.1",
        "cnvpytor-combined-rd-baf-1.3.2",
        "expansionhunter-grch38-female-5.0.0",
        "expansionhunter-grch38-male-5.0.0",
        "hla-la-short-read-1.0.4",
        "spechla-short-read-bam-1.0.12",
        "gatk-mitochondria-normal-shifted-4.6.2.0",
        LEGACY_READBACK,
        CURRENT_READBACK,
    }
)
HLA_OPTION_FIELDS = frozenset(
    {
        "caller_database_name",
        "caller_database_version",
        "caller_database_identity_sha256",
        "nomenclature_database_name",
        "nomenclature_version",
        "nomenclature_identity_sha256",
        "database_age_limitation",
    }
)
PRIVATE_PATH_KEYS = frozenset({"host_path", "path", "tree_manifest_path"})


class ProfileAdapterError(RuntimeError):
    """A caller output could not be bound to portable typed evidence."""


def _fail(message: str) -> None:
    raise ProfileAdapterError(message)


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _load_object(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        _fail(f"{label}_missing")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ProfileAdapterError(f"{label}_unreadable") from exc
    if not isinstance(value, dict):
        _fail(f"{label}_invalid")
    return value


def _path_values(value: Any) -> set[str]:
    paths: set[str] = set()
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if key in PRIVATE_PATH_KEYS and isinstance(nested, str) and nested:
                paths.add(nested)
            paths.update(_path_values(nested))
    elif isinstance(value, list):
        for nested in value:
            paths.update(_path_values(nested))
    return paths


def _forbidden_values(
    binding_spec: Mapping[str, Any],
    lifecycle: Path,
    output: Path,
    options_path: Path | None,
    forbidden_path: Path | None,
) -> tuple[str, ...]:
    values = _path_values(binding_spec) | {str(lifecycle), str(output)}
    if options_path is not None:
        values.add(str(options_path))
    if forbidden_path is not None:
        values.add(str(forbidden_path))
        try:
            values.update(
                line.strip()
                for line in forbidden_path.read_text(encoding="utf-8").splitlines()
                if line.strip() and not line.lstrip().startswith("#")
            )
        except (OSError, UnicodeError) as exc:
            raise ProfileAdapterError("forbidden_binding_file_unreadable") from exc
    return tuple(sorted(value for value in values if value))


def _assert_portable(value: object, forbidden: Sequence[str]) -> None:
    rendered = json.dumps(value, sort_keys=True, ensure_ascii=True, allow_nan=False)
    if any(contains_forbidden_binding(rendered, binding) for binding in forbidden):
        _fail("portable_adapter_summary_contains_private_binding")


def _load_options(path: Path | None, profile_id: str) -> tuple[dict[str, Any], str]:
    if path is None:
        selected: dict[str, Any] = {}
    else:
        payload = _load_object(path, "adapter_options")
        if set(payload) != {"schema_version", "profiles"} or payload.get("schema_version") != 1:
            _fail("adapter_options_schema_invalid")
        profiles = payload.get("profiles")
        if not isinstance(profiles, dict):
            _fail("adapter_options_profiles_invalid")
        raw = profiles.get(profile_id, {})
        if not isinstance(raw, dict):
            _fail("selected_adapter_options_invalid")
        selected = dict(raw)
    hla_profile = profile_id in {
        "hla-la-short-read-1.0.4",
        "spechla-short-read-bam-1.0.12",
    }
    if hla_profile:
        if set(selected) != HLA_OPTION_FIELDS:
            _fail("hla_adapter_options_incomplete")
        if not all(
            isinstance(selected[field], str) and selected[field]
            for field in HLA_OPTION_FIELDS - {"database_age_limitation"}
        ) or not isinstance(selected["database_age_limitation"], bool):
            _fail("hla_adapter_options_invalid")
    elif selected:
        _fail("profile_does_not_accept_adapter_options")
    return selected, _canonical_sha256(selected)


def _binding_mount(binding_spec: Mapping[str, Any], role: str) -> dict[str, Any]:
    mounts = binding_spec.get("mounts")
    if not isinstance(mounts, dict) or not isinstance(mounts.get(role), dict):
        _fail(f"required_binding_role_missing_{role}")
    return dict(mounts[role])


def _verify_context(
    profile_id: str,
    lifecycle: Path,
    binding_spec: Mapping[str, Any],
    forbidden: Sequence[str],
) -> dict[str, Any]:
    if profile_id not in SUPPORTED_PROFILES:
        _fail("profile_has_no_strict_adapter")
    profile = load_caller_profile(resolve_shipped_profile(profile_id), forbidden_bindings=forbidden)
    lifecycle_summary_path = lifecycle / "lifecycle-summary.json"
    lifecycle_summary = _load_object(lifecycle_summary_path, "lifecycle_summary")
    receipt_path = lifecycle / "caller-receipt.json"
    receipt, receipt_file_sha256 = load_caller_run_receipt(
        receipt_path, forbidden_bindings=forbidden
    )
    artifact_spec = _load_object(lifecycle / "plan/artifact-spec.json", "artifact_spec")
    boundary = receipt.get("run_boundary")
    artifact = receipt.get("raw_artifact")
    if not isinstance(boundary, dict) or not isinstance(artifact, dict):
        _fail("caller_receipt_boundary_invalid")
    if (
        lifecycle_summary.get("status") != "PASS"
        or lifecycle_summary.get("profile_id") != profile_id
        or lifecycle_summary.get("profile_identity_sha256") != profile.identity_sha256
        or lifecycle_summary.get("caller_run_identity_sha256")
        != boundary.get("caller_run_identity_sha256")
        or lifecycle_summary.get("receipt_identity_sha256")
        != receipt.get("receipt_identity_sha256")
        or boundary.get("caller_name") != profile.payload["caller_name"]
        or boundary.get("caller_version") != profile.payload["caller_version"]
    ):
        _fail("lifecycle_profile_receipt_identity_mismatch")
    if (
        artifact_spec.get("schema_version") != 1
        or artifact_spec.get("raw_artifact_id") != artifact.get("raw_artifact_id")
        or artifact_spec.get("identity_scheme") != artifact.get("identity_scheme")
        or artifact_spec.get("derived_payload") != artifact.get("derived_payload")
    ):
        _fail("artifact_spec_receipt_mismatch")
    component_paths = artifact_spec.get("components")
    receipt_components = artifact.get("components")
    if not isinstance(component_paths, dict) or not isinstance(receipt_components, list):
        _fail("artifact_components_invalid")
    receipt_by_role = {row.get("role"): row for row in receipt_components if isinstance(row, dict)}
    if set(component_paths) != set(receipt_by_role) or len(receipt_by_role) != len(
        receipt_components
    ):
        _fail("artifact_component_roles_mismatch")
    paths: dict[str, Path] = {}
    for role, raw_path in component_paths.items():
        if not isinstance(role, str) or not isinstance(raw_path, str):
            _fail("artifact_component_path_invalid")
        path = Path(raw_path)
        expected = receipt_by_role[role]
        if path.is_symlink() or not path.is_file():
            _fail("artifact_component_missing")
        if sha256_file(path) != expected.get("sha256") or path.stat().st_size != expected.get(
            "bytes"
        ):
            _fail("artifact_component_identity_changed")
        paths[role] = path
    mounts = binding_spec.get("mounts")
    expected_roles = {row["role"] for row in profile.payload["mounts"]}
    if (
        binding_spec.get("schema_version") != 1
        or not isinstance(mounts, dict)
        or set(mounts) != expected_roles
    ):
        _fail("binding_spec_profile_roles_mismatch")
    return {
        "profile": profile,
        "lifecycle_summary_path": lifecycle_summary_path,
        "lifecycle_summary": lifecycle_summary,
        "receipt_path": receipt_path,
        "receipt": receipt,
        "receipt_file_sha256": receipt_file_sha256,
        "boundary": boundary,
        "artifact": artifact,
        "components": paths,
    }


def _common_provenance(context: Mapping[str, Any], technical_source_role: str) -> dict[str, str]:
    boundary = context["boundary"]
    artifact = context["artifact"]
    return {
        "reference_sha256": str(boundary["reference_sha256"]),
        "technical_source_role": technical_source_role,
        "caller_run_identity_sha256": str(boundary["caller_run_identity_sha256"]),
        "resource_manifest_sha256": str(boundary["resource_manifest_sha256"]),
        "source_alignment_sha256": str(boundary["source_alignment_sha256"]),
        "raw_artifact_id": str(artifact["raw_artifact_id"]),
    }


def _materialize(
    profile_id: str,
    staging: Path,
    context: Mapping[str, Any],
    binding_spec: Mapping[str, Any],
    technical_source_role: str,
    options: Mapping[str, Any],
    forbidden: tuple[str, ...],
) -> list[tuple[str, Path, dict[str, Any]]]:
    profile = context["profile"]
    version = str(profile.payload["caller_version"])
    components = context["components"]
    common = _common_provenance(context, technical_source_role)
    results: list[tuple[str, Path, dict[str, Any]]] = []
    if profile_id == "manta-germline-1.6.0":
        output = staging / "manta.parquet"
        provenance = MantaProvenance(**common, caller_version=version)
        result = materialize_manta_germline_calls(
            components["primary_vcf"], output, provenance, forbidden_bindings=forbidden
        )
        results.append(("structural_variants", output, result))
    elif profile_id == "delly-short-read-2.5.1":
        output = staging / "delly.parquet"
        provenance = DellyProvenance(**common, caller_version=version)
        result = materialize_delly_calls(
            components["primary_bcf"], output, provenance, forbidden_bindings=forbidden
        )
        results.append(("structural_variants", output, result))
    elif profile_id == "cnvpytor-combined-rd-baf-1.3.2":
        baf = _binding_mount(binding_spec, "source_baf_vcf")
        baf_path = Path(str(baf.get("host_path")))
        baf_sha256 = str(baf.get("sha256"))
        if baf_path.is_symlink() or not baf_path.is_file() or sha256_file(baf_path) != baf_sha256:
            _fail("source_baf_binding_identity_changed")
        for bin_size in (10_000, 100_000):
            output = staging / f"cnvpytor-{bin_size}.parquet"
            provenance = CnvpytorCombinedProvenance(
                **common,
                raw_artifact_sha256=str(context["artifact"]["raw_artifact_sha256"]),
                source_baf_vcf_sha256=baf_sha256,
                bin_size=bin_size,
                rd_flags_suffix="_GC_mask",
                snp_flags_suffix="_mask",
                rd_mask_applied=True,
                baf_mask_applied=True,
                caller_version=version,
            )
            result = materialize_cnvpytor_combined_calls(
                components[f"calls_{bin_size}"],
                components["pytor"],
                components[f"hdf5_audit_{bin_size}"],
                baf_path,
                output,
                provenance,
                forbidden_bindings=forbidden,
            )
            results.append(("copy_number", output, result))
    elif profile_id.startswith("expansionhunter-grch38-"):
        output = staging / "expansionhunter.parquet"
        provenance = ExpansionHunterProvenance(**common, caller_version=version)
        result = materialize_expansionhunter_calls(
            components["expansionhunter_json_sha256"],
            components["expansionhunter_vcf_sha256"],
            output,
            provenance,
            forbidden_bindings=forbidden,
        )
        results.append(("repeat_expansions", output, result))
    elif profile_id == "hla-la-short-read-1.0.4":
        output = staging / "hla-la.parquet"
        provenance = HlaLaProvenance(
            **common,
            **options,
            caller_version=version,
        )
        result = materialize_hla_la_bestguess(
            components["bestguess_g_group"],
            output,
            provenance,
            forbidden_bindings=forbidden,
        )
        results.append(("hla_typing", output, result))
    elif profile_id == "spechla-short-read-bam-1.0.12":
        output = staging / "spechla.parquet"
        provenance = SpecHLAProvenance(
            **common,
            **options,
            caller_version=version,
        )
        result = materialize_spechla_results(
            components["full_result"],
            components["g_group_result"],
            components["details_result"],
            output,
            provenance,
            forbidden_bindings=forbidden,
        )
        results.append(("hla_typing", output, result))
    elif profile_id == "gatk-mitochondria-normal-shifted-4.6.2.0":
        output = staging / "mitochondrial.parquet"
        provenance = MitochondrialProvenance(
            **common, caller_version=version, circular_boundary_audited=True
        )
        result = materialize_filtered_mutect_vcf(
            components["filtered_vcf"],
            output,
            provenance,
            forbidden_bindings=forbidden,
        )
        results.append(("mitochondrial_heteroplasmy", output, result))
    elif profile_id in {LEGACY_READBACK, CURRENT_READBACK}:
        output = staging / "readback.parquet"
        provenance = ReadbackProvenance(**common, caller_version=version)
        result = materialize_readback_capture(
            components["raw_capture_json"],
            output,
            provenance,
            allow_empty=profile_id == CURRENT_READBACK,
            forbidden_bindings=forbidden,
        )
        results.append(("read_backed_evidence", output, result))
    else:
        raise AssertionError("unreachable profile adapter")
    return results


def _file_metadata(
    domain: str,
    path: Path,
    context: Mapping[str, Any],
    technical_source_role: str,
    forbidden: tuple[str, ...],
) -> dict[str, Any]:
    try:
        table = pq.read_table(path)
    except (OSError, ArrowInvalid) as exc:
        raise ProfileAdapterError("adapter_output_parquet_unreadable") from exc
    if table.schema != schema_for(domain):
        _fail("adapter_output_schema_invalid")
    boundary = context["boundary"]
    artifact = context["artifact"]
    rows = validate_rows(
        domain,
        table.to_pylist(),
        expected_reference_sha256=str(boundary["reference_sha256"]),
        forbidden_bindings=forbidden,
    )
    for row in rows:
        if (
            row["technical_source_role"] != technical_source_role
            or row["caller_name"] != boundary["caller_name"]
            or row["caller_version"] != boundary["caller_version"]
            or row["caller_run_identity_sha256"] != boundary["caller_run_identity_sha256"]
            or row["resource_manifest_sha256"] != boundary["resource_manifest_sha256"]
            or row["source_alignment_sha256"] != boundary["source_alignment_sha256"]
            or row["raw_artifact_id"] != artifact["raw_artifact_id"]
            or row["raw_artifact_sha256"] != artifact["raw_artifact_sha256"]
        ):
            _fail("typed_row_receipt_boundary_mismatch")
    return {
        "domain": domain,
        "path": path.name,
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "rows": len(rows),
        "canonical_row_multiset": canonical_row_multiset(rows),
    }


def _summary_without_identity(summary: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in summary.items() if key != "adapter_identity_sha256"}


def validate_profile_adapter_output(
    output_directory: str | Path,
    lifecycle_directory: str | Path,
    binding_spec_path: str | Path,
    technical_source_role: str,
    *,
    profile_id: str,
    adapter_options_path: str | Path | None = None,
    forbidden_binding_file: str | Path | None = None,
) -> dict[str, Any]:
    output = Path(output_directory)
    lifecycle = Path(lifecycle_directory)
    binding_path = Path(binding_spec_path)
    options_path = Path(adapter_options_path) if adapter_options_path else None
    forbidden_path = Path(forbidden_binding_file) if forbidden_binding_file else None
    if output.is_symlink() or not output.is_dir():
        _fail("adapter_output_directory_invalid")
    binding_spec = _load_object(binding_path, "binding_spec")
    forbidden = _forbidden_values(binding_spec, lifecycle, output, options_path, forbidden_path)
    options, options_identity = _load_options(options_path, profile_id)
    context = _verify_context(profile_id, lifecycle, binding_spec, forbidden)
    summary = _load_object(output / SUMMARY_NAME, "adapter_summary")
    identity = summary.get("adapter_identity_sha256")
    if not isinstance(identity, str) or identity != _canonical_sha256(
        _summary_without_identity(summary)
    ):
        _fail("adapter_summary_identity_changed")
    if (
        summary.get("schema_version") != ADAPTER_SCHEMA_VERSION
        or summary.get("status") != "PASS"
        or summary.get("profile_id") != profile_id
        or summary.get("profile_identity_sha256") != context["profile"].identity_sha256
        or summary.get("adapter_options_identity_sha256") != options_identity
        or summary.get("technical_source_role") != technical_source_role
        or summary.get("lifecycle_summary_sha256") != sha256_file(context["lifecycle_summary_path"])
        or summary.get("caller_receipt_file_sha256") != context["receipt_file_sha256"]
        or summary.get("caller_receipt_identity_sha256")
        != context["receipt"]["receipt_identity_sha256"]
    ):
        _fail("adapter_summary_boundary_changed")
    declared = summary.get("files")
    if not isinstance(declared, list) or any(not isinstance(row, dict) for row in declared):
        _fail("adapter_summary_files_invalid")
    actual_files = {path.name for path in output.iterdir() if path.is_file()}
    expected_files = {SUMMARY_NAME, *(str(row.get("path")) for row in declared)}
    if actual_files != expected_files:
        _fail("adapter_output_file_inventory_changed")
    for row in declared:
        raw_name = row.get("path")
        if (
            not isinstance(raw_name, str)
            or Path(raw_name).name != raw_name
            or "/" in raw_name
            or "\\" in raw_name
        ):
            _fail("adapter_output_file_name_invalid")
        path = output / raw_name
        if path.is_symlink() or not path.is_file():
            _fail("adapter_output_file_missing")
        if sha256_file(path) != row.get("sha256") or path.stat().st_size != row.get("bytes"):
            _fail("adapter_output_file_identity_changed")
    observed = [
        _file_metadata(
            str(row.get("domain")),
            output / str(row.get("path")),
            context,
            technical_source_role,
            forbidden,
        )
        for row in declared
    ]
    if observed != declared:
        _fail("adapter_output_file_identity_changed")
    if bool(declared) != bool(summary.get("extension_required")):
        _fail("adapter_extension_required_boundary_changed")
    _assert_portable(summary, forbidden)
    return {
        "schema_version": ADAPTER_SCHEMA_VERSION,
        "status": "PASS",
        "profile_id": profile_id,
        "typed_files": len(declared),
        "typed_rows": sum(int(row["rows"]) for row in declared),
        "adapter_identity_sha256": identity,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "options_revalidated": options == _load_options(options_path, profile_id)[0],
    }


def materialize_profile_adapter(
    lifecycle_directory: str | Path,
    binding_spec_path: str | Path,
    technical_source_role: str,
    output_directory: str | Path,
    *,
    profile_id: str,
    adapter_options_path: str | Path | None = None,
    forbidden_binding_file: str | Path | None = None,
) -> dict[str, Any]:
    if technical_source_role not in SOURCE_ROLES:
        _fail("technical_source_role_invalid")
    lifecycle = Path(lifecycle_directory)
    binding_path = Path(binding_spec_path)
    output = Path(output_directory)
    options_path = Path(adapter_options_path) if adapter_options_path else None
    forbidden_path = Path(forbidden_binding_file) if forbidden_binding_file else None
    if output.exists():
        validated = validate_profile_adapter_output(
            output,
            lifecycle,
            binding_path,
            technical_source_role,
            profile_id=profile_id,
            adapter_options_path=options_path,
            forbidden_binding_file=forbidden_path,
        )
        return {**validated, "resume": "reused"}
    binding_spec = _load_object(binding_path, "binding_spec")
    forbidden = _forbidden_values(binding_spec, lifecycle, output, options_path, forbidden_path)
    options, options_identity = _load_options(options_path, profile_id)
    context = _verify_context(profile_id, lifecycle, binding_spec, forbidden)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent))
    results = _materialize(
        profile_id,
        staging,
        context,
        binding_spec,
        technical_source_role,
        options,
        forbidden,
    )
    files: list[dict[str, Any]] = []
    parsing: list[dict[str, Any]] = []
    for domain, path, result in results:
        parse_summary = result.get("parsing")
        if not isinstance(parse_summary, dict) or result.get("status") != "PASS":
            _fail("strict_adapter_did_not_return_pass")
        parsing.append({"domain": domain, "summary": parse_summary})
        if result.get("output") is not None:
            files.append(_file_metadata(domain, path, context, technical_source_role, forbidden))
        elif path.exists():
            _fail("adapter_created_undeclared_empty_output")
    files.sort(key=lambda row: (row["domain"], row["path"]))
    parsing.sort(key=lambda row: (row["domain"], _canonical_sha256(row["summary"])))
    boundary = context["boundary"]
    artifact = context["artifact"]
    summary: dict[str, Any] = {
        "schema_version": ADAPTER_SCHEMA_VERSION,
        "status": "PASS",
        "profile_id": profile_id,
        "profile_identity_sha256": context["profile"].identity_sha256,
        "adapter_options_identity_sha256": options_identity,
        "technical_source_role": technical_source_role,
        "reference_sha256": boundary["reference_sha256"],
        "source_alignment_sha256": boundary["source_alignment_sha256"],
        "resource_manifest_sha256": boundary["resource_manifest_sha256"],
        "caller_run_identity_sha256": boundary["caller_run_identity_sha256"],
        "caller_receipt_identity_sha256": context["receipt"]["receipt_identity_sha256"],
        "caller_receipt_file_sha256": context["receipt_file_sha256"],
        "lifecycle_summary_sha256": sha256_file(context["lifecycle_summary_path"]),
        "raw_artifact_id": artifact["raw_artifact_id"],
        "raw_artifact_sha256": artifact["raw_artifact_sha256"],
        "raw_component_roles": sorted(context["components"]),
        "files": files,
        "parsing": parsing,
        "extension_required": bool(files),
        "source_rows_removed": 0,
        "logical_subject_count": 1,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "raw_artifacts_embedded": False,
    }
    _assert_portable(summary, forbidden)
    summary["adapter_identity_sha256"] = _canonical_sha256(summary)
    with (staging / SUMMARY_NAME).open("xb") as handle:
        handle.write((json.dumps(summary, indent=2, sort_keys=True) + "\n").encode())
        handle.flush()
        os.fsync(handle.fileno())
    if output.exists():
        _fail("adapter_output_directory_already_exists")
    os.rename(staging, output)
    validated = validate_profile_adapter_output(
        output,
        lifecycle,
        binding_path,
        technical_source_role,
        profile_id=profile_id,
        adapter_options_path=options_path,
        forbidden_binding_file=forbidden_path,
    )
    return {**validated, "resume": "created"}


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("materialize", "validate"):
        command = subparsers.add_parser(name)
        command.add_argument("--profile-id", required=True)
        command.add_argument("--lifecycle-directory", type=Path, required=True)
        command.add_argument("--bindings", type=Path, required=True)
        command.add_argument("--technical-source-role", required=True)
        command.add_argument("--output-directory", type=Path, required=True)
        command.add_argument("--adapter-options", type=Path)
        command.add_argument("--forbidden-bindings-file", type=Path)
    args = parser.parse_args()
    try:
        keywords = {
            "profile_id": args.profile_id,
            "adapter_options_path": args.adapter_options,
            "forbidden_binding_file": args.forbidden_bindings_file,
        }
        if args.command == "materialize":
            result = materialize_profile_adapter(
                args.lifecycle_directory,
                args.bindings,
                args.technical_source_role,
                args.output_directory,
                **keywords,
            )
        else:
            result = validate_profile_adapter_output(
                args.output_directory,
                args.lifecycle_directory,
                args.bindings,
                args.technical_source_role,
                **keywords,
            )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except Exception:
        print(
            "profile adapter failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
