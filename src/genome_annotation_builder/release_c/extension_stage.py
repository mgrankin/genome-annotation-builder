"""Atomic, subject-bound staging for portable Release-C extension tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from .caller_profile import load_caller_profile
from .caller_receipt import (
    load_caller_run_receipt,
    validate_caller_run_receipt_payload,
)
from .portable_model import (
    DOMAINS,
    FORBIDDEN_PORTABLE_TEXT,
    SHA256,
    SOURCE_ROLES,
    PortableEvidenceError,
    canonical_row_multiset,
    contains_forbidden_binding,
    schema_for,
    sha256_file,
    validate_rows,
)
from .profile_catalog import resolve_shipped_profile

DOMAIN_FILENAMES = {domain: f"{domain}.parquet" for domain in DOMAINS}
ROW_CALLER_TO_GATE_CALLER = {
    "Manta": "manta",
    "DELLY": "delly",
    "CNVpytor": "cnvpytor",
    "ExpansionHunter": "expansionhunter",
    "SpecHLA": "spechla",
    "HLA*LA": "hla_la",
    "Cyrius": "cyrius",
    "Aldy": "aldy",
    "GATK Mutect2 mitochondrial mode": "gatk_mitochondria",
    "GAB bounded local pileup": "readback",
}
ACCEPTED_PRIVATE_MANIFEST_STATES = frozenset(
    {
        "PASS",
        "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH",
    }
)
MANIFEST_FILENAME = "release_c_extensions.json"
MANTA_PROFILE_ID = "manta-germline-1.6.0"
ADAPTER_SUMMARY_FIELDS = frozenset(
    {
        "schema_version",
        "status",
        "profile_id",
        "profile_identity_sha256",
        "adapter_options_identity_sha256",
        "technical_source_role",
        "reference_sha256",
        "source_alignment_sha256",
        "resource_manifest_sha256",
        "caller_run_identity_sha256",
        "caller_receipt_identity_sha256",
        "caller_receipt_file_sha256",
        "lifecycle_summary_sha256",
        "raw_artifact_id",
        "raw_artifact_sha256",
        "raw_component_roles",
        "files",
        "parsing",
        "extension_required",
        "source_rows_removed",
        "logical_subject_count",
        "subject_identifier_persisted",
        "private_paths_persisted",
        "raw_artifacts_embedded",
        "adapter_identity_sha256",
    }
)


@dataclass(frozen=True)
class SubjectBinding:
    """Private in-memory proof that supplied technical datasets belong together."""

    logical_subject: str = field(repr=False)
    reference_sha256: str
    allowed_alignment_sha256: frozenset[str]
    dataset_count: int
    forbidden_bindings: tuple[str, ...] = field(repr=False)


@dataclass(frozen=True)
class GateBinding:
    caller: str
    alignment_sha256: str
    reference_sha256: str
    run_identity_sha256: str
    report_sha256: str

    @property
    def key(self) -> tuple[str, str]:
        return self.caller, self.alignment_sha256


@dataclass(frozen=True)
class CallerReceiptBinding:
    caller_name: str
    caller_version: str
    gate_caller: str
    alignment_sha256: str
    reference_sha256: str
    resource_manifest_sha256: str
    caller_run_identity_sha256: str
    raw_artifact_id: str
    raw_artifact_sha256: str
    gate_run_identity_sha256: str
    gate_report_sha256: str
    receipt_identity_sha256: str
    receipt_file_sha256: str
    payload: dict[str, Any] = field(repr=False)

    @property
    def gate_key(self) -> tuple[str, str]:
        return self.gate_caller, self.alignment_sha256

    @property
    def row_key(self) -> tuple[str, ...]:
        return (
            self.caller_name,
            self.caller_version,
            self.caller_run_identity_sha256,
            self.resource_manifest_sha256,
            self.alignment_sha256,
            self.reference_sha256,
            self.raw_artifact_id,
            self.raw_artifact_sha256,
        )


@dataclass(frozen=True)
class AdapterBinding:
    """Portable strict-adapter checkpoint bound to exactly one caller receipt."""

    receipt_key: tuple[str, ...]
    gate_key: tuple[str, str]
    profile_id: str
    technical_source_role: str
    adapter_identity_sha256: str
    adapter_summary_file_sha256: str
    domain_rows: tuple[tuple[str, int], ...]
    payload: dict[str, Any] = field(repr=False)


def _sha256_json(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError(f"{label} is unreadable") from exc
    if not isinstance(value, dict):
        raise PortableEvidenceError(f"{label} must be a JSON object")
    return value


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise PortableEvidenceError(f"{label} is not a SHA-256")
    return value


def _nonnegative_integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise PortableEvidenceError(f"{label} is not a non-negative integer")
    return value


def _validate_canonical_multiset(
    value: Any,
    label: str,
    *,
    expected_count: int | None = None,
) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {
        "count",
        "sha256_xor",
        "sha256_sum_mod_2_256",
    }:
        raise PortableEvidenceError(f"{label} canonical row multiset is invalid")
    count = _nonnegative_integer(value.get("count"), f"{label} row count")
    _sha(value.get("sha256_xor"), f"{label} XOR digest")
    _sha(value.get("sha256_sum_mod_2_256"), f"{label} sum digest")
    if expected_count is not None and count != expected_count:
        raise PortableEvidenceError(f"{label} canonical row count differs")
    return value


def _domain_multiset_key(domain: str, multiset: dict[str, Any]) -> tuple[str, int, str, str]:
    return (
        domain,
        int(multiset["count"]),
        str(multiset["sha256_xor"]),
        str(multiset["sha256_sum_mod_2_256"]),
    )


@cache
def _profile_boundary(profile_id: str) -> tuple[str, str, str, str]:
    profile = load_caller_profile(resolve_shipped_profile(profile_id))
    payload = profile.payload
    return (
        profile.identity_sha256,
        str(payload["caller_name"]),
        str(payload["caller_version"]),
        str(payload["gate_caller"]),
    )


def _validate_manta_reconciliation(summary: dict[str, Any], typed_rows: int) -> None:
    integer_fields = (
        "source_records",
        "typed_rows",
        "record_reconciliation_delta",
        "omitted_noncanonical_records",
        "omitted_noncanonical_non_bnd_records",
        "omitted_noncanonical_bnd_records",
        "omitted_noncanonical_bnd_pairs",
        "core_to_noncanonical_bnd_pairs",
        "noncanonical_to_noncanonical_bnd_pairs",
        "source_bnd_records",
        "source_reciprocal_bnd_pairs",
        "bnd_records",
        "reciprocal_bnd_pairs",
    )
    counts = {
        field: _nonnegative_integer(summary.get(field), f"Manta {field}")
        for field in integer_fields
    }
    if summary.get("portable_contig_scope") != "grch38_core_chr1_22_x_y_m":
        raise PortableEvidenceError("Manta portable contig scope is invalid")
    if (
        counts["typed_rows"] != typed_rows
        or counts["source_records"] != counts["typed_rows"] + counts["omitted_noncanonical_records"]
        or counts["record_reconciliation_delta"] != counts["omitted_noncanonical_records"]
        or counts["omitted_noncanonical_records"]
        != counts["omitted_noncanonical_non_bnd_records"]
        + counts["omitted_noncanonical_bnd_records"]
        or counts["omitted_noncanonical_bnd_records"]
        != 2 * counts["omitted_noncanonical_bnd_pairs"]
        or counts["omitted_noncanonical_bnd_pairs"]
        != counts["core_to_noncanonical_bnd_pairs"]
        + counts["noncanonical_to_noncanonical_bnd_pairs"]
        or counts["bnd_records"] != 2 * counts["reciprocal_bnd_pairs"]
        or counts["source_bnd_records"] != 2 * counts["source_reciprocal_bnd_pairs"]
        or counts["source_bnd_records"]
        != counts["bnd_records"] + counts["omitted_noncanonical_bnd_records"]
    ):
        raise PortableEvidenceError("Manta raw/typed/scope-omission arithmetic differs")


def _validate_parsing_summary(
    profile_id: str,
    domain: str,
    summary: Any,
    *,
    raw_artifact_id: str,
    expected_rows: int,
) -> dict[str, Any]:
    if not isinstance(summary, dict):
        raise PortableEvidenceError("strict-adapter parsing summary is not an object")
    if (
        summary.get("schema_version") != 1
        or summary.get("status") != "PASS"
        or summary.get("private_paths_persisted") is not False
        or summary.get("raw_artifact_id") != raw_artifact_id
    ):
        raise PortableEvidenceError("strict-adapter parsing summary boundary is invalid")
    identifier_flags = [
        summary[key]
        for key in ("sample_identifier_persisted", "subject_identifier_persisted")
        if key in summary
    ]
    if not identifier_flags or any(value is not False for value in identifier_flags):
        raise PortableEvidenceError("strict-adapter parsing summary persists an identifier")
    _validate_canonical_multiset(
        summary.get("canonical_row_multiset"),
        "strict-adapter parsing summary",
        expected_count=expected_rows,
    )
    for metric in ("source_records", "typed_rows", "record_reconciliation_delta"):
        if metric in summary:
            _nonnegative_integer(summary[metric], f"strict-adapter {metric}")
    if all(
        field in summary
        for field in ("source_records", "typed_rows", "record_reconciliation_delta")
    ) and (
        summary["typed_rows"] != expected_rows
        or summary["source_records"] - summary["typed_rows"]
        != summary["record_reconciliation_delta"]
    ):
        raise PortableEvidenceError("strict-adapter record reconciliation arithmetic differs")
    if profile_id == MANTA_PROFILE_ID:
        if domain != "structural_variants":
            raise PortableEvidenceError("Manta parsing summary names the wrong domain")
        _validate_manta_reconciliation(summary, expected_rows)
    return summary


def _adapter_binding_from_payload(
    payload: dict[str, Any],
    adapter_summary_file_sha256: str,
    receipts: dict[tuple[str, ...], CallerReceiptBinding],
    *,
    forbidden_bindings: tuple[str, ...],
) -> AdapterBinding:
    _assert_portable_payload(payload, forbidden_bindings)
    if set(payload) != ADAPTER_SUMMARY_FIELDS:
        raise PortableEvidenceError("strict-adapter summary schema is invalid")
    identity = _sha(payload.get("adapter_identity_sha256"), "adapter identity")
    unsigned = dict(payload)
    unsigned.pop("adapter_identity_sha256")
    if identity != _sha256_json(unsigned):
        raise PortableEvidenceError("strict-adapter summary identity changed")
    if (
        payload.get("schema_version") != 1
        or payload.get("status") != "PASS"
        or payload.get("source_rows_removed") != 0
        or payload.get("logical_subject_count") != 1
        or payload.get("subject_identifier_persisted") is not False
        or payload.get("private_paths_persisted") is not False
        or payload.get("raw_artifacts_embedded") is not False
    ):
        raise PortableEvidenceError(
            "strict-adapter summary privacy/losslessness boundary is invalid"
        )
    profile_id = payload.get("profile_id")
    if not isinstance(profile_id, str):
        raise PortableEvidenceError("strict-adapter profile ID is invalid")
    profile_identity, caller_name, caller_version, gate_caller = _profile_boundary(profile_id)
    if payload.get("profile_identity_sha256") != profile_identity:
        raise PortableEvidenceError("strict-adapter shipped profile identity differs")
    for digest_field in (
        "adapter_options_identity_sha256",
        "reference_sha256",
        "source_alignment_sha256",
        "resource_manifest_sha256",
        "caller_run_identity_sha256",
        "caller_receipt_identity_sha256",
        "caller_receipt_file_sha256",
        "lifecycle_summary_sha256",
        "raw_artifact_sha256",
    ):
        _sha(payload.get(digest_field), f"strict-adapter {digest_field}")
    role = payload.get("technical_source_role")
    if role not in SOURCE_ROLES:
        raise PortableEvidenceError("strict-adapter technical source role is invalid")
    raw_artifact_id = payload.get("raw_artifact_id")
    if not isinstance(raw_artifact_id, str) or not raw_artifact_id:
        raise PortableEvidenceError("strict-adapter raw artifact ID is invalid")

    matching_receipts = [
        receipt
        for receipt in receipts.values()
        if receipt.caller_name == caller_name
        and receipt.caller_version == caller_version
        and receipt.gate_caller == gate_caller
        and receipt.reference_sha256 == payload["reference_sha256"]
        and receipt.alignment_sha256 == payload["source_alignment_sha256"]
        and receipt.resource_manifest_sha256 == payload["resource_manifest_sha256"]
        and receipt.caller_run_identity_sha256 == payload["caller_run_identity_sha256"]
        and receipt.receipt_identity_sha256 == payload["caller_receipt_identity_sha256"]
        and receipt.receipt_file_sha256 == payload["caller_receipt_file_sha256"]
        and receipt.raw_artifact_id == raw_artifact_id
        and receipt.raw_artifact_sha256 == payload["raw_artifact_sha256"]
    ]
    if len(matching_receipts) != 1:
        raise PortableEvidenceError("strict-adapter summary lacks one exact caller receipt")
    receipt = matching_receipts[0]
    raw_components = receipt.payload["raw_artifact"]["components"]
    expected_component_roles = sorted(str(component["role"]) for component in raw_components)
    raw_component_roles = payload.get("raw_component_roles")
    if raw_component_roles != expected_component_roles:
        raise PortableEvidenceError("strict-adapter raw component roles differ from its receipt")

    files = payload.get("files")
    if not isinstance(files, list) or any(not isinstance(row, dict) for row in files):
        raise PortableEvidenceError("strict-adapter file inventory is invalid")
    normalized_files: list[dict[str, Any]] = []
    domains_seen: set[str] = set()
    paths_seen: set[str] = set()
    domain_rows: dict[str, int] = {}
    unbound_file_multisets: Counter[tuple[str, int, str, str]] = Counter()
    for metadata in files:
        if set(metadata) != {
            "domain",
            "path",
            "sha256",
            "bytes",
            "rows",
            "canonical_row_multiset",
        }:
            raise PortableEvidenceError("strict-adapter file metadata schema is invalid")
        domain = metadata.get("domain")
        raw_path = metadata.get("path")
        if (
            domain not in DOMAINS
            or not isinstance(raw_path, str)
            or Path(raw_path).name != raw_path
            or "/" in raw_path
            or "\\" in raw_path
            or raw_path in paths_seen
        ):
            raise PortableEvidenceError("strict-adapter file metadata is not canonical")
        rows = _nonnegative_integer(metadata.get("rows"), "strict-adapter file rows")
        if rows < 1 or _nonnegative_integer(metadata.get("bytes"), "strict-adapter file bytes") < 1:
            raise PortableEvidenceError("strict-adapter declared an empty evidence file")
        _sha(metadata.get("sha256"), "strict-adapter file digest")
        multiset = _validate_canonical_multiset(
            metadata.get("canonical_row_multiset"),
            "strict-adapter file",
            expected_count=rows,
        )
        domain = str(domain)
        domains_seen.add(domain)
        paths_seen.add(raw_path)
        domain_rows[domain] = domain_rows.get(domain, 0) + rows
        unbound_file_multisets[_domain_multiset_key(domain, multiset)] += 1
        normalized_files.append(metadata)
    if normalized_files != sorted(
        normalized_files, key=lambda row: (str(row["domain"]), str(row["path"]))
    ) or bool(normalized_files) is not bool(payload.get("extension_required")):
        raise PortableEvidenceError("strict-adapter file ordering/extension boundary differs")

    parsing = payload.get("parsing")
    if not isinstance(parsing, list) or not parsing:
        raise PortableEvidenceError("strict-adapter parsing inventory is invalid")
    parsed_domains: set[str] = set()
    normalized_parsing: list[dict[str, Any]] = []
    for row in parsing:
        if not isinstance(row, dict) or set(row) != {"domain", "summary"}:
            raise PortableEvidenceError("strict-adapter parsing row is invalid")
        domain = row.get("domain")
        if domain not in DOMAINS:
            raise PortableEvidenceError("strict-adapter parsing domain is invalid")
        domain = str(domain)
        summary = row.get("summary")
        if not isinstance(summary, dict):
            raise PortableEvidenceError("strict-adapter parsing summary is not an object")
        multiset = _validate_canonical_multiset(
            summary.get("canonical_row_multiset"),
            "strict-adapter parsing summary",
        )
        expected_rows = int(multiset["count"])
        _validate_parsing_summary(
            profile_id,
            domain,
            summary,
            raw_artifact_id=raw_artifact_id,
            expected_rows=expected_rows,
        )
        if expected_rows:
            key = _domain_multiset_key(domain, multiset)
            if not unbound_file_multisets[key]:
                raise PortableEvidenceError(
                    "strict-adapter parsing summary lacks exact typed-file metadata"
                )
            unbound_file_multisets[key] -= 1
        parsed_domains.add(domain)
        normalized_parsing.append(row)
    if (
        domains_seen - parsed_domains
        or any(unbound_file_multisets.values())
        or normalized_parsing
        != sorted(
            normalized_parsing,
            key=lambda row: (str(row["domain"]), _sha256_json(row["summary"])),
        )
    ):
        raise PortableEvidenceError("strict-adapter parsing inventory/order differs")
    return AdapterBinding(
        receipt_key=receipt.row_key,
        gate_key=receipt.gate_key,
        profile_id=profile_id,
        technical_source_role=str(role),
        adapter_identity_sha256=identity,
        adapter_summary_file_sha256=_sha(
            adapter_summary_file_sha256, "adapter summary file digest"
        ),
        domain_rows=tuple(sorted(domain_rows.items())),
        payload=payload,
    )


def _load_adapter_summary(
    path: Path,
    subject_binding: SubjectBinding,
    receipts: dict[tuple[str, ...], CallerReceiptBinding],
) -> AdapterBinding:
    if path.is_symlink() or not path.is_file():
        raise PortableEvidenceError("strict-adapter summary must be a regular file")
    payload = _load_json(path, "strict-adapter summary")
    return _adapter_binding_from_payload(
        payload,
        sha256_file(path),
        receipts,
        forbidden_bindings=subject_binding.forbidden_bindings,
    )


def _walk_strings(value: Any):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from _walk_strings(key)
            yield from _walk_strings(item)
    elif isinstance(value, (list, tuple, set, frozenset)):
        for item in value:
            yield from _walk_strings(item)


def _assert_portable_payload(value: Any, forbidden: tuple[str, ...]) -> None:
    for text in _walk_strings(value):
        if FORBIDDEN_PORTABLE_TEXT.search(text):
            raise PortableEvidenceError("Release-C stage contains a private path or secret")
        if any(contains_forbidden_binding(text, binding) for binding in forbidden):
            raise PortableEvidenceError("Release-C stage contains a private binding")


def subject_binding_from_private_manifest(
    manifest_path: str | Path, dataset_ids: list[str]
) -> SubjectBinding:
    """Resolve one logical subject without persisting its identifiers or cohort digest."""

    path = Path(manifest_path)
    manifest = _load_json(path, "private input manifest")
    if (
        manifest.get("schema_version") != 1
        or manifest.get("status") not in ACCEPTED_PRIVATE_MANIFEST_STATES
    ):
        raise PortableEvidenceError("private input manifest schema/status is invalid")
    if not dataset_ids or len(dataset_ids) != len(set(dataset_ids)):
        raise PortableEvidenceError("private dataset binding list is empty or duplicated")
    datasets = manifest.get("datasets")
    if not isinstance(datasets, list) or any(not isinstance(row, dict) for row in datasets):
        raise PortableEvidenceError("private input manifest datasets are invalid")
    resolved: list[dict[str, Any]] = []
    for dataset_id in dataset_ids:
        matches = [row for row in datasets if row.get("dataset_id") == dataset_id]
        if len(matches) != 1:
            raise PortableEvidenceError("private dataset binding is not unique")
        resolved.append(matches[0])
    subjects = {row.get("logical_subject") for row in resolved}
    if len(subjects) != 1 or not isinstance(next(iter(subjects)), str) or not next(iter(subjects)):
        raise PortableEvidenceError("technical datasets do not resolve to one logical subject")
    logical_subject = str(next(iter(subjects)))
    reference = manifest.get("reference")
    if not isinstance(reference, dict) or not str(reference.get("assembly", "")).startswith(
        "GRCh38"
    ):
        raise PortableEvidenceError("private subject binding reference is not GRCh38")
    reference_sha256 = _sha(reference.get("fasta_sha256"), "private reference digest")

    alignment_hashes: set[str] = set()
    dataset_primary_alignment_hashes: list[str] = []
    identifiers: set[str] = {logical_subject, *dataset_ids}
    paths: set[str] = {str(path)}
    for key, value in reference.items():
        if key.endswith("_path") and isinstance(value, str) and value:
            paths.add(value)
    for row in resolved:
        for identity_key in ("laboratory",):
            value = row.get(identity_key)
            if isinstance(value, str) and value:
                identifiers.add(value)
        row_alignment_hashes: list[str] = []
        for alignment_key in ("original_markdup_bam", "bqsr_cram"):
            spec = row.get(alignment_key)
            if not isinstance(spec, dict):
                continue
            alignment_digest = _sha(spec.get("sha256"), "private alignment digest")
            alignment_hashes.add(alignment_digest)
            row_alignment_hashes.append(alignment_digest)
            for path_key in ("path", "index_path"):
                value = spec.get(path_key)
                if isinstance(value, str) and value:
                    paths.add(value)
        if not row_alignment_hashes:
            raise PortableEvidenceError("private dataset binding has no alignment digest")
        dataset_primary_alignment_hashes.append(row_alignment_hashes[0])
        baf = row.get("baf_vcf")
        if isinstance(baf, dict):
            sample = baf.get("sample")
            if isinstance(sample, str) and sample:
                identifiers.add(sample)
            for path_key in ("path", "index_path"):
                value = baf.get(path_key)
                if isinstance(value, str) and value:
                    paths.add(value)
    if not alignment_hashes:
        raise PortableEvidenceError("private subject binding has no alignment digests")
    if len(dataset_primary_alignment_hashes) != len(set(dataset_primary_alignment_hashes)):
        raise PortableEvidenceError("private subject binding duplicates a technical alignment")
    forbidden = identifiers | paths
    forbidden.update(hashlib.sha256(value.encode()).hexdigest() for value in identifiers)
    return SubjectBinding(
        logical_subject=logical_subject,
        reference_sha256=reference_sha256,
        allowed_alignment_sha256=frozenset(alignment_hashes),
        dataset_count=len(resolved),
        forbidden_bindings=tuple(sorted(forbidden)),
    )


def _load_gate(path: Path, subject_binding: SubjectBinding) -> GateBinding:
    report = _load_json(path, "input-gate report")
    _assert_portable_payload(report, subject_binding.forbidden_bindings)
    if report.get("schema_version") != 1 or report.get("status") != "PASS":
        raise PortableEvidenceError("input-gate report is not PASS")
    if (
        report.get("subject_identifier_persisted") is not False
        or report.get("private_paths_persisted") is not False
        or report.get("private_manifest_digest_persisted") is not False
        or report.get("logical_subject_count") != 1
    ):
        raise PortableEvidenceError("input-gate report violates its privacy/subject boundary")
    caller = report.get("caller")
    if not isinstance(caller, str) or caller not in set(ROW_CALLER_TO_GATE_CALLER.values()):
        raise PortableEvidenceError("input-gate caller is unsupported")
    alignment = report.get("selected_alignment")
    reference = report.get("reference")
    if not isinstance(alignment, dict) or not isinstance(reference, dict):
        raise PortableEvidenceError("input-gate artifact binding is incomplete")
    alignment_sha256 = _sha(alignment.get("sha256"), "input-gate alignment digest")
    reference_sha256 = _sha(reference.get("sha256"), "input-gate reference digest")
    if alignment_sha256 not in subject_binding.allowed_alignment_sha256:
        raise PortableEvidenceError("input-gate alignment is outside the private subject binding")
    if reference_sha256 != subject_binding.reference_sha256:
        raise PortableEvidenceError("input-gate reference differs from the private subject binding")
    return GateBinding(
        caller=caller,
        alignment_sha256=alignment_sha256,
        reference_sha256=reference_sha256,
        run_identity_sha256=_sha(report.get("run_identity_sha256"), "input-gate run identity"),
        report_sha256=sha256_file(path),
    )


def _receipt_binding_from_payload(
    payload: dict[str, Any], receipt_file_sha256: str
) -> CallerReceiptBinding:
    boundary = payload["run_boundary"]
    input_gate = boundary["input_gate"]
    artifact = payload["raw_artifact"]
    return CallerReceiptBinding(
        caller_name=str(boundary["caller_name"]),
        caller_version=str(boundary["caller_version"]),
        gate_caller=str(boundary["gate_caller"]),
        alignment_sha256=str(boundary["source_alignment_sha256"]),
        reference_sha256=str(boundary["reference_sha256"]),
        resource_manifest_sha256=str(boundary["resource_manifest_sha256"]),
        caller_run_identity_sha256=str(boundary["caller_run_identity_sha256"]),
        raw_artifact_id=str(artifact["raw_artifact_id"]),
        raw_artifact_sha256=str(artifact["raw_artifact_sha256"]),
        gate_run_identity_sha256=str(input_gate["run_identity_sha256"]),
        gate_report_sha256=str(input_gate["report_sha256"]),
        receipt_identity_sha256=str(payload["receipt_identity_sha256"]),
        receipt_file_sha256=_sha(receipt_file_sha256, "caller receipt file digest"),
        payload=payload,
    )


def _validate_receipt_boundary(
    receipt: CallerReceiptBinding,
    gates: dict[tuple[str, str], GateBinding],
    *,
    subject_binding: SubjectBinding | None = None,
) -> None:
    expected_gate_caller = ROW_CALLER_TO_GATE_CALLER.get(receipt.caller_name)
    if expected_gate_caller is None or receipt.gate_caller != expected_gate_caller:
        raise PortableEvidenceError("caller receipt has no exact caller/input-gate mapping")
    if subject_binding is not None:
        if receipt.alignment_sha256 not in subject_binding.allowed_alignment_sha256:
            raise PortableEvidenceError("caller receipt is outside the private subject binding")
        if receipt.reference_sha256 != subject_binding.reference_sha256:
            raise PortableEvidenceError("caller receipt uses another private reference")
    gate = gates.get(receipt.gate_key)
    if gate is None:
        raise PortableEvidenceError("caller receipt lacks its exact input-gate report")
    if (
        receipt.reference_sha256 != gate.reference_sha256
        or receipt.gate_run_identity_sha256 != gate.run_identity_sha256
        or receipt.gate_report_sha256 != gate.report_sha256
    ):
        raise PortableEvidenceError("caller receipt and input-gate boundary differ")


def _load_receipt(
    path: Path,
    subject_binding: SubjectBinding,
    gates: dict[tuple[str, str], GateBinding],
) -> CallerReceiptBinding:
    payload, receipt_file_sha256 = load_caller_run_receipt(
        path, forbidden_bindings=subject_binding.forbidden_bindings
    )
    receipt = _receipt_binding_from_payload(payload, receipt_file_sha256)
    _validate_receipt_boundary(receipt, gates, subject_binding=subject_binding)
    return receipt


def _read_domain(domain: str, path: Path, subject_binding: SubjectBinding) -> list[dict[str, Any]]:
    if path.is_symlink() or not path.is_file():
        raise PortableEvidenceError("Release-C extension input must be a regular Parquet file")
    table = pq.read_table(path)
    if table.schema != schema_for(domain):
        raise PortableEvidenceError(f"{domain} input has the wrong physical schema")
    return validate_rows(
        domain,
        table.to_pylist(),
        expected_reference_sha256=subject_binding.reference_sha256,
        forbidden_bindings=subject_binding.forbidden_bindings,
    )


def _bind_rows_to_gates(
    rows_by_domain: dict[str, list[dict[str, Any]]],
    gates: dict[tuple[str, str], GateBinding],
    subject_binding: SubjectBinding,
) -> set[tuple[str, str]]:
    used: set[tuple[str, str]] = set()
    for rows in rows_by_domain.values():
        for row in rows:
            gate_caller = ROW_CALLER_TO_GATE_CALLER.get(str(row["caller_name"]))
            if gate_caller is None:
                raise PortableEvidenceError("typed row has no declared input-gate caller mapping")
            alignment_sha256 = str(row["source_alignment_sha256"])
            if alignment_sha256 not in subject_binding.allowed_alignment_sha256:
                raise PortableEvidenceError("typed row is outside the private subject binding")
            key = gate_caller, alignment_sha256
            gate = gates.get(key)
            if gate is None:
                raise PortableEvidenceError("typed row lacks its exact caller/alignment input gate")
            if row["reference_sha256"] != gate.reference_sha256:
                raise PortableEvidenceError("typed row and input gate use different references")
            used.add(key)
    return used


def _row_receipt_key(row: dict[str, Any]) -> tuple[str, ...]:
    return (
        str(row["caller_name"]),
        str(row["caller_version"]),
        str(row["caller_run_identity_sha256"]),
        str(row["resource_manifest_sha256"]),
        str(row["source_alignment_sha256"]),
        str(row["reference_sha256"]),
        str(row["raw_artifact_id"]),
        str(row["raw_artifact_sha256"]),
    )


def _bind_rows_to_receipts(
    rows_by_domain: dict[str, list[dict[str, Any]]],
    receipts: dict[tuple[str, ...], CallerReceiptBinding],
) -> set[tuple[str, ...]]:
    used: set[tuple[str, ...]] = set()
    for rows in rows_by_domain.values():
        for row in rows:
            key = _row_receipt_key(row)
            if key not in receipts:
                raise PortableEvidenceError(
                    "typed row lacks its exact caller-run/raw-artifact receipt"
                )
            used.add(key)
    return used


def _copy_verified(source: Path, destination: Path, expected_sha256: str) -> None:
    with source.open("rb") as read_handle, destination.open("xb") as write_handle:
        while chunk := read_handle.read(16 * 1024 * 1024):
            write_handle.write(chunk)
        write_handle.flush()
        os.fsync(write_handle.fileno())
    if sha256_file(destination) != expected_sha256:
        raise PortableEvidenceError("copied Release-C extension digest changed")


def _canonical_gate_rows(gates: dict[tuple[str, str], GateBinding]) -> list[dict[str, Any]]:
    return [
        {
            "caller": gate.caller,
            "selected_alignment_sha256": gate.alignment_sha256,
            "reference_sha256": gate.reference_sha256,
            "gate_run_identity_sha256": gate.run_identity_sha256,
            "gate_report_sha256": gate.report_sha256,
        }
        for _, gate in sorted(gates.items())
    ]


def _canonical_receipt_rows(
    receipts: dict[tuple[str, ...], CallerReceiptBinding],
) -> list[dict[str, Any]]:
    return [
        {
            "receipt_file_sha256": receipt.receipt_file_sha256,
            "receipt": receipt.payload,
        }
        for _, receipt in sorted(receipts.items())
    ]


def _canonical_adapter_rows(
    adapters: dict[tuple[str, ...], AdapterBinding],
) -> list[dict[str, Any]]:
    return [
        {
            "adapter_summary_file_sha256": adapter.adapter_summary_file_sha256,
            "adapter_summary": adapter.payload,
        }
        for _, adapter in sorted(adapters.items())
    ]


def _adapter_domain_counts(
    adapters: dict[tuple[str, ...], AdapterBinding],
) -> dict[tuple[tuple[str, ...], str], int]:
    return {
        (receipt_key, domain): rows
        for receipt_key, adapter in adapters.items()
        for domain, rows in adapter.domain_rows
    }


def _row_domain_counts(
    rows_by_domain: dict[str, list[dict[str, Any]]],
) -> dict[tuple[tuple[str, ...], str], int]:
    counts: dict[tuple[tuple[str, ...], str], int] = {}
    for domain, rows in rows_by_domain.items():
        for row in rows:
            key = (_row_receipt_key(row), domain)
            counts[key] = counts.get(key, 0) + 1
    return counts


def _adapter_reconciliation(
    adapters: dict[tuple[str, ...], AdapterBinding],
) -> dict[str, Any]:
    parsing_summaries = 0
    typed_rows = 0
    exact_summaries = 0
    exact_source_records = 0
    exact_typed_rows = 0
    exact_delta = 0
    scope = {
        "policy": "grch38_core_chr1_22_x_y_m",
        "adapter_runs": 0,
        "source_records": 0,
        "typed_rows": 0,
        "omitted_records": 0,
        "omitted_non_bnd_records": 0,
        "omitted_bnd_records": 0,
        "omitted_bnd_pairs": 0,
        "core_to_noncanonical_bnd_pairs": 0,
        "noncanonical_to_noncanonical_bnd_pairs": 0,
    }
    for adapter in adapters.values():
        typed_rows += sum(rows for _, rows in adapter.domain_rows)
        parsing = adapter.payload["parsing"]
        parsing_summaries += len(parsing)
        for row in parsing:
            summary = row["summary"]
            if all(
                field in summary
                for field in ("source_records", "typed_rows", "record_reconciliation_delta")
            ):
                exact_summaries += 1
                exact_source_records += int(summary["source_records"])
                exact_typed_rows += int(summary["typed_rows"])
                exact_delta += int(summary["record_reconciliation_delta"])
            if adapter.profile_id == MANTA_PROFILE_ID:
                scope["adapter_runs"] += 1
                scope["source_records"] += int(summary["source_records"])
                scope["typed_rows"] += int(summary["typed_rows"])
                scope["omitted_records"] += int(summary["omitted_noncanonical_records"])
                scope["omitted_non_bnd_records"] += int(
                    summary["omitted_noncanonical_non_bnd_records"]
                )
                scope["omitted_bnd_records"] += int(summary["omitted_noncanonical_bnd_records"])
                scope["omitted_bnd_pairs"] += int(summary["omitted_noncanonical_bnd_pairs"])
                scope["core_to_noncanonical_bnd_pairs"] += int(
                    summary["core_to_noncanonical_bnd_pairs"]
                )
                scope["noncanonical_to_noncanonical_bnd_pairs"] += int(
                    summary["noncanonical_to_noncanonical_bnd_pairs"]
                )
    return {
        "adapter_runs": len(adapters),
        "parsing_summaries": parsing_summaries,
        "typed_rows": typed_rows,
        "exact_record_reconciliations": {
            "summaries": exact_summaries,
            "source_records": exact_source_records,
            "typed_rows": exact_typed_rows,
            "record_reconciliation_delta": exact_delta,
        },
        "portable_core_scope": scope,
    }


def stage_release_c_extensions(
    output_dir: str | Path,
    domain_parquets: dict[str, str | Path],
    input_gate_reports: list[str | Path],
    caller_run_receipts: list[str | Path],
    adapter_summaries: list[str | Path],
    subject_binding: SubjectBinding,
) -> dict[str, Any]:
    """Build a deterministic, non-empty, portable extension-set directory."""

    output = Path(output_dir)
    if output.exists():
        raise PortableEvidenceError("refusing to replace an existing Release-C extension stage")
    if not domain_parquets:
        raise PortableEvidenceError("Release-C extension stage has no evidence domains")
    unknown = set(domain_parquets) - set(DOMAINS)
    if unknown:
        raise PortableEvidenceError("Release-C extension stage names an unknown domain")
    if not input_gate_reports:
        raise PortableEvidenceError("Release-C extension stage has no input-gate reports")
    if not caller_run_receipts:
        raise PortableEvidenceError("Release-C extension stage has no caller-run receipts")
    if not adapter_summaries:
        raise PortableEvidenceError("Release-C extension stage has no strict-adapter summaries")

    forbidden = set(subject_binding.forbidden_bindings)
    forbidden.update(str(Path(path)) for path in domain_parquets.values())
    forbidden.update(str(Path(path)) for path in input_gate_reports)
    forbidden.update(str(Path(path)) for path in caller_run_receipts)
    forbidden.update(str(Path(path)) for path in adapter_summaries)
    forbidden.add(str(output))
    portable_forbidden = tuple(sorted(value for value in forbidden if value))

    rows_by_domain = {
        domain: _read_domain(domain, Path(path), subject_binding)
        for domain, path in sorted(domain_parquets.items())
    }
    gates: dict[tuple[str, str], GateBinding] = {}
    for raw_path in input_gate_reports:
        path = Path(raw_path)
        if path.is_symlink() or not path.is_file():
            raise PortableEvidenceError("input-gate report must be a regular file")
        gate = _load_gate(path, subject_binding)
        if gate.key in gates:
            raise PortableEvidenceError("duplicate caller/alignment input-gate report")
        gates[gate.key] = gate
    used_gates = _bind_rows_to_gates(rows_by_domain, gates, subject_binding)
    receipts: dict[tuple[str, ...], CallerReceiptBinding] = {}
    for raw_path in caller_run_receipts:
        path = Path(raw_path)
        if path.is_symlink() or not path.is_file():
            raise PortableEvidenceError("caller-run receipt must be a regular file")
        receipt = _load_receipt(path, subject_binding, gates)
        if receipt.row_key in receipts:
            raise PortableEvidenceError("duplicate caller-run/raw-artifact receipt")
        receipts[receipt.row_key] = receipt
    used_receipts = _bind_rows_to_receipts(rows_by_domain, receipts)
    adapters: dict[tuple[str, ...], AdapterBinding] = {}
    adapter_identities: set[str] = set()
    for raw_path in adapter_summaries:
        adapter = _load_adapter_summary(Path(raw_path), subject_binding, receipts)
        if adapter.receipt_key in adapters:
            raise PortableEvidenceError("duplicate strict-adapter summary for one caller receipt")
        if adapter.adapter_identity_sha256 in adapter_identities:
            raise PortableEvidenceError("strict-adapter identity is duplicated")
        adapters[adapter.receipt_key] = adapter
        adapter_identities.add(adapter.adapter_identity_sha256)
    if set(adapters) != set(receipts):
        raise PortableEvidenceError("caller receipts and strict-adapter summaries differ")
    adapter_gate_keys = {adapter.gate_key for adapter in adapters.values()}
    if adapter_gate_keys != set(gates):
        raise PortableEvidenceError("input gates and strict-adapter caller runs differ")
    if not used_gates.issubset(adapter_gate_keys) or not used_receipts.issubset(set(adapters)):
        raise PortableEvidenceError("typed evidence is outside its strict-adapter boundary")
    if _row_domain_counts(rows_by_domain) != _adapter_domain_counts(adapters):
        raise PortableEvidenceError("typed evidence rows differ from strict-adapter accounting")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent))
    try:
        domain_manifest: dict[str, dict[str, Any]] = {}
        for domain, rows in rows_by_domain.items():
            source = Path(domain_parquets[domain])
            destination = staging / DOMAIN_FILENAMES[domain]
            digest = sha256_file(source)
            _copy_verified(source, destination, digest)
            caller_boundaries = sorted(
                {
                    (
                        str(row["caller_name"]),
                        str(row["caller_version"]),
                        str(row["resource_manifest_sha256"]),
                        str(row["caller_run_identity_sha256"]),
                    )
                    for row in rows
                }
            )
            domain_manifest[domain] = {
                "path": DOMAIN_FILENAMES[domain],
                "sha256": digest,
                "bytes": destination.stat().st_size,
                "rows": len(rows),
                "canonical_row_multiset": canonical_row_multiset(rows),
                "caller_boundaries": [
                    {
                        "caller_name": caller_name,
                        "caller_version": caller_version,
                        "resource_manifest_sha256": resource_sha256,
                        "caller_run_identity_sha256": run_sha256,
                    }
                    for caller_name, caller_version, resource_sha256, run_sha256 in caller_boundaries
                ],
            }
        manifest: dict[str, Any] = {
            "schema_version": 1,
            "extension_schema_version": 1,
            "status": "PASS",
            "assembly": "GRCh38",
            "reference_sha256": subject_binding.reference_sha256,
            "subject_scope": "one_logical_subject",
            "logical_subject_count": 1,
            "source_dataset_count": subject_binding.dataset_count,
            "subject_binding_verified": True,
            "subject_identifier_persisted": False,
            "private_paths_persisted": False,
            "private_manifest_digest_persisted": False,
            "source_artifacts_embedded": False,
            "source_rows_removed": 0,
            "empty_optional_files": 0,
            "domains": domain_manifest,
            "input_gates": _canonical_gate_rows(gates),
            "caller_runs": _canonical_receipt_rows(receipts),
            "adapter_runs": _canonical_adapter_rows(adapters),
            "adapter_reconciliation": _adapter_reconciliation(adapters),
        }
        manifest["extension_set_identity_sha256"] = _sha256_json(manifest)
        _assert_portable_payload(manifest, portable_forbidden)
        manifest_path = staging / MANIFEST_FILENAME
        with manifest_path.open("xb") as handle:
            handle.write((json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode())
            handle.flush()
            os.fsync(handle.fileno())
        validation = validate_release_c_extension_stage(
            staging,
            expected_reference_sha256=subject_binding.reference_sha256,
            forbidden_bindings=portable_forbidden,
        )
        if output.exists():
            raise PortableEvidenceError("refusing to replace an existing Release-C extension stage")
        os.rename(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return {"status": "PASS", "manifest": manifest, "validation": validation}


def validate_release_c_extension_stage(
    stage_dir: str | Path,
    *,
    expected_reference_sha256: str | None = None,
    forbidden_bindings: tuple[str, ...] = (),
) -> dict[str, Any]:
    root = Path(stage_dir)
    if root.is_symlink() or not root.is_dir():
        raise PortableEvidenceError("Release-C extension stage is not a regular directory")
    manifest = _load_json(root / MANIFEST_FILENAME, "Release-C extension manifest")
    _assert_portable_payload(manifest, forbidden_bindings)
    if manifest.get("schema_version") != 1 or manifest.get("status") != "PASS":
        raise PortableEvidenceError("Release-C extension manifest schema/status is invalid")
    if manifest.get("assembly") != "GRCh38":
        raise PortableEvidenceError("Release-C extension stage is not GRCh38")
    reference_sha256 = _sha(manifest.get("reference_sha256"), "extension reference digest")
    if expected_reference_sha256 and reference_sha256 != expected_reference_sha256:
        raise PortableEvidenceError("Release-C extension stage uses the wrong reference")
    if (
        manifest.get("subject_scope") != "one_logical_subject"
        or manifest.get("logical_subject_count") != 1
        or manifest.get("subject_binding_verified") is not True
        or manifest.get("subject_identifier_persisted") is not False
        or manifest.get("private_paths_persisted") is not False
        or manifest.get("private_manifest_digest_persisted") is not False
        or manifest.get("source_artifacts_embedded") is not False
        or manifest.get("source_rows_removed") != 0
        or manifest.get("empty_optional_files") != 0
        or not isinstance(manifest.get("source_dataset_count"), int)
        or manifest.get("source_dataset_count") < 1
    ):
        raise PortableEvidenceError("Release-C extension privacy/losslessness boundary is invalid")
    identity = manifest.pop("extension_set_identity_sha256", None)
    try:
        if identity != _sha256_json(manifest):
            raise PortableEvidenceError("Release-C extension-set identity changed")
    finally:
        manifest["extension_set_identity_sha256"] = identity

    raw_gates = manifest.get("input_gates")
    if not isinstance(raw_gates, list) or not raw_gates:
        raise PortableEvidenceError("Release-C extension manifest has no input gates")
    gate_keys: set[tuple[str, str]] = set()
    manifest_gates: dict[tuple[str, str], GateBinding] = {}
    for gate in raw_gates:
        if not isinstance(gate, dict) or set(gate) != {
            "caller",
            "selected_alignment_sha256",
            "reference_sha256",
            "gate_run_identity_sha256",
            "gate_report_sha256",
        }:
            raise PortableEvidenceError("Release-C extension input-gate row is invalid")
        caller = gate.get("caller")
        alignment = _sha(gate.get("selected_alignment_sha256"), "gate alignment digest")
        if caller not in set(ROW_CALLER_TO_GATE_CALLER.values()):
            raise PortableEvidenceError("Release-C extension input-gate caller is invalid")
        if _sha(gate.get("reference_sha256"), "gate reference digest") != reference_sha256:
            raise PortableEvidenceError("Release-C input gate uses another reference")
        gate_run_identity = _sha(gate.get("gate_run_identity_sha256"), "gate run identity")
        gate_report_sha256 = _sha(gate.get("gate_report_sha256"), "gate report digest")
        key = str(caller), alignment
        if key in gate_keys:
            raise PortableEvidenceError("Release-C extension input gate is duplicated")
        gate_keys.add(key)
        manifest_gates[key] = GateBinding(
            caller=str(caller),
            alignment_sha256=alignment,
            reference_sha256=reference_sha256,
            run_identity_sha256=gate_run_identity,
            report_sha256=gate_report_sha256,
        )

    raw_receipts = manifest.get("caller_runs")
    if not isinstance(raw_receipts, list) or not raw_receipts:
        raise PortableEvidenceError("Release-C extension manifest has no caller-run receipts")
    receipt_keys: set[tuple[str, ...]] = set()
    receipt_identities: set[str] = set()
    manifest_receipts: dict[tuple[str, ...], CallerReceiptBinding] = {}
    for wrapper in raw_receipts:
        if not isinstance(wrapper, dict) or set(wrapper) != {"receipt_file_sha256", "receipt"}:
            raise PortableEvidenceError("Release-C caller-run receipt wrapper is invalid")
        receipt_file_sha256 = _sha(wrapper.get("receipt_file_sha256"), "caller receipt file digest")
        payload = wrapper.get("receipt")
        if not isinstance(payload, dict):
            raise PortableEvidenceError("Release-C caller-run receipt payload is invalid")
        normalized = validate_caller_run_receipt_payload(
            payload, forbidden_bindings=forbidden_bindings
        )
        if normalized != payload:
            raise PortableEvidenceError("Release-C caller-run receipt is not canonical")
        receipt = _receipt_binding_from_payload(normalized, receipt_file_sha256)
        _validate_receipt_boundary(receipt, manifest_gates)
        if receipt.reference_sha256 != reference_sha256:
            raise PortableEvidenceError("Release-C caller receipt uses another reference")
        if receipt.row_key in receipt_keys:
            raise PortableEvidenceError("Release-C caller-run receipt is duplicated")
        if receipt.receipt_identity_sha256 in receipt_identities:
            raise PortableEvidenceError("Release-C caller receipt identity is duplicated")
        receipt_keys.add(receipt.row_key)
        receipt_identities.add(receipt.receipt_identity_sha256)
        manifest_receipts[receipt.row_key] = receipt

    raw_adapters = manifest.get("adapter_runs")
    if not isinstance(raw_adapters, list) or not raw_adapters:
        raise PortableEvidenceError("Release-C extension manifest has no strict-adapter runs")
    manifest_adapters: dict[tuple[str, ...], AdapterBinding] = {}
    adapter_identities: set[str] = set()
    for wrapper in raw_adapters:
        if not isinstance(wrapper, dict) or set(wrapper) != {
            "adapter_summary_file_sha256",
            "adapter_summary",
        }:
            raise PortableEvidenceError("Release-C strict-adapter wrapper is invalid")
        payload = wrapper.get("adapter_summary")
        if not isinstance(payload, dict):
            raise PortableEvidenceError("Release-C strict-adapter payload is invalid")
        adapter = _adapter_binding_from_payload(
            payload,
            _sha(wrapper.get("adapter_summary_file_sha256"), "adapter summary file digest"),
            manifest_receipts,
            forbidden_bindings=forbidden_bindings,
        )
        if adapter.receipt_key in manifest_adapters:
            raise PortableEvidenceError("Release-C strict-adapter receipt is duplicated")
        if adapter.adapter_identity_sha256 in adapter_identities:
            raise PortableEvidenceError("Release-C strict-adapter identity is duplicated")
        manifest_adapters[adapter.receipt_key] = adapter
        adapter_identities.add(adapter.adapter_identity_sha256)
    if set(manifest_adapters) != receipt_keys:
        raise PortableEvidenceError("Release-C caller receipts and strict-adapter runs differ")
    declared_reconciliation = manifest.get("adapter_reconciliation")
    observed_reconciliation = _adapter_reconciliation(manifest_adapters)
    if declared_reconciliation != observed_reconciliation:
        raise PortableEvidenceError("Release-C strict-adapter reconciliation changed")

    domains = manifest.get("domains")
    if not isinstance(domains, dict) or not domains or set(domains) - set(DOMAINS):
        raise PortableEvidenceError("Release-C extension domain inventory is invalid")
    expected_files = {MANIFEST_FILENAME}
    total_rows = 0
    concordance: dict[str, dict[str, int]] = {}
    technical_source_roles: dict[str, dict[str, int]] = {}
    used_gates: set[tuple[str, str]] = set()
    used_receipts: set[tuple[str, ...]] = set()
    rows_by_domain: dict[str, list[dict[str, Any]]] = {}
    for domain, metadata in sorted(domains.items()):
        if not isinstance(metadata, dict) or metadata.get("path") != DOMAIN_FILENAMES[domain]:
            raise PortableEvidenceError("Release-C extension path is not canonical")
        expected_files.add(DOMAIN_FILENAMES[domain])
        path = root / DOMAIN_FILENAMES[domain]
        if path.is_symlink() or not path.is_file():
            raise PortableEvidenceError("Release-C extension table is missing")
        if sha256_file(path) != _sha(metadata.get("sha256"), "extension table digest"):
            raise PortableEvidenceError("Release-C extension table digest changed")
        if path.stat().st_size != metadata.get("bytes"):
            raise PortableEvidenceError("Release-C extension table byte count changed")
        table = pq.read_table(path)
        if table.schema != schema_for(domain):
            raise PortableEvidenceError("Release-C extension physical schema changed")
        rows = validate_rows(
            domain,
            table.to_pylist(),
            expected_reference_sha256=reference_sha256,
            forbidden_bindings=forbidden_bindings,
        )
        if len(rows) != metadata.get("rows") or not rows:
            raise PortableEvidenceError("Release-C extension row count changed or is empty")
        if canonical_row_multiset(rows) != metadata.get("canonical_row_multiset"):
            raise PortableEvidenceError("Release-C extension row multiset changed")
        rows_by_domain[domain] = rows
        observed_boundaries = sorted(
            {
                (
                    str(row["caller_name"]),
                    str(row["caller_version"]),
                    str(row["resource_manifest_sha256"]),
                    str(row["caller_run_identity_sha256"]),
                )
                for row in rows
            }
        )
        declared_boundaries = [
            (
                item.get("caller_name"),
                item.get("caller_version"),
                item.get("resource_manifest_sha256"),
                item.get("caller_run_identity_sha256"),
            )
            for item in metadata.get("caller_boundaries", [])
            if isinstance(item, dict)
        ]
        if observed_boundaries != declared_boundaries:
            raise PortableEvidenceError("Release-C extension caller boundary changed")
        concordance_counts: dict[str, int] = {}
        source_role_counts: dict[str, int] = {}
        for row in rows:
            gate_caller = ROW_CALLER_TO_GATE_CALLER.get(str(row["caller_name"]))
            key = str(gate_caller), str(row["source_alignment_sha256"])
            if gate_caller is None or key not in gate_keys:
                raise PortableEvidenceError("Release-C extension row lacks its input gate")
            used_gates.add(key)
            receipt_key = _row_receipt_key(row)
            if receipt_key not in receipt_keys:
                raise PortableEvidenceError("Release-C extension row lacks its caller receipt")
            used_receipts.add(receipt_key)
            status = str(row["concordance_status"])
            role = str(row["technical_source_role"])
            concordance_counts[status] = concordance_counts.get(status, 0) + 1
            source_role_counts[role] = source_role_counts.get(role, 0) + 1
        concordance[domain] = dict(sorted(concordance_counts.items()))
        technical_source_roles[domain] = dict(sorted(source_role_counts.items()))
        total_rows += len(rows)
    actual_files = {path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file()}
    if actual_files != expected_files:
        raise PortableEvidenceError("Release-C extension stage has undeclared files")
    adapter_gate_keys = {adapter.gate_key for adapter in manifest_adapters.values()}
    if adapter_gate_keys != gate_keys:
        raise PortableEvidenceError("Release-C input gates and strict-adapter runs differ")
    if not used_gates.issubset(gate_keys) or not used_receipts.issubset(receipt_keys):
        raise PortableEvidenceError("Release-C typed evidence escaped its caller boundary")
    if _row_domain_counts(rows_by_domain) != _adapter_domain_counts(manifest_adapters):
        raise PortableEvidenceError("Release-C typed rows differ from adapter accounting")
    return {
        "schema_version": 1,
        "status": "PASS",
        "domains": len(domains),
        "rows": total_rows,
        "input_gates": len(gate_keys),
        "caller_runs": len(receipt_keys),
        "adapter_runs": len(manifest_adapters),
        "adapter_reconciliation": observed_reconciliation,
        "concordance": concordance,
        "technical_source_roles": technical_source_roles,
        "extension_set_identity_sha256": identity,
        "subject_identifier_persisted": False,
        "private_paths_persisted": False,
        "source_rows_removed": 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    stage_parser = subparsers.add_parser("stage")
    stage_parser.add_argument("--private-manifest", type=Path, required=True)
    stage_parser.add_argument("--dataset-id", action="append", required=True)
    stage_parser.add_argument(
        "--domain",
        action="append",
        required=True,
        metavar="NAME=PARQUET",
    )
    stage_parser.add_argument("--input-gate", action="append", type=Path, required=True)
    stage_parser.add_argument("--caller-receipt", action="append", type=Path, required=True)
    stage_parser.add_argument("--adapter-summary", action="append", type=Path, required=True)
    stage_parser.add_argument("--output", type=Path, required=True)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--stage", type=Path, required=True)
    validate_parser.add_argument("--reference-sha256")
    args = parser.parse_args()
    if args.command == "stage":
        domain_parquets: dict[str, Path] = {}
        for raw in args.domain:
            if "=" not in raw:
                parser.error("--domain must use NAME=PARQUET")
            name, path = raw.split("=", 1)
            if not name or not path or name in domain_parquets:
                parser.error("--domain names and paths must be non-empty and unique")
            domain_parquets[name] = Path(path)
        binding = subject_binding_from_private_manifest(
            args.private_manifest,
            args.dataset_id,
        )
        result = stage_release_c_extensions(
            args.output,
            domain_parquets,
            args.input_gate,
            args.caller_receipt,
            args.adapter_summary,
            binding,
        )
        print(json.dumps(result["validation"], indent=2, sort_keys=True))
        return 0
    if args.command == "validate":
        report = validate_release_c_extension_stage(
            args.stage, expected_reference_sha256=args.reference_sha256
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
