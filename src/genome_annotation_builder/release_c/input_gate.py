"""Privacy-minimized, fail-closed input gate for read-backed callers."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

ACCEPTED_MANIFEST_STATES = frozenset(
    {
        "PASS",
        "PARTIAL_PASS_PENDING_IMMEDIATE_PRECALL_LARGE_FILE_REHASH",
    }
)
CALLER_NAME = re.compile(r"^[a-z0-9][a-z0-9_.-]{0,63}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class InputGatePolicy:
    minimum_free_bytes: int = 50 * 1024**3
    indexed_probe_bases: int = 100_000


DEFAULT_INPUT_GATE_POLICY = InputGatePolicy()


class InputGateFailure(RuntimeError):
    """A sanitized failure whose details belong only in an authorized private log."""

    def __init__(self, code: str):
        super().__init__(code)
        self.code = code


def _fail(code: str) -> None:
    raise InputGateFailure(code)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_json(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def _atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
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
            raise InputGateFailure("output_already_exists") from exc
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _run(command: list[str]) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(command, check=False, capture_output=True, text=True)
    except OSError as exc:
        raise InputGateFailure("required_tool_unavailable") from exc


def _safe_int(value: Any, code: str) -> int:
    if isinstance(value, bool):
        _fail(code)
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise InputGateFailure(code) from exc
    if parsed < 0:
        _fail(code)
    return parsed


def _digest(value: Any, code: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        _fail(code)
    return value


def _require_bound_file(spec: dict[str, Any], *, prefix: str) -> tuple[Path, Path]:
    try:
        path = Path(spec["path"])
        index = Path(spec["index_path"])
    except (KeyError, TypeError) as exc:
        raise InputGateFailure(f"{prefix}_binding_invalid") from exc
    expected_bytes = _safe_int(spec.get("bytes"), f"{prefix}_binding_invalid")
    expected_index_bytes = _safe_int(spec.get("index_bytes"), f"{prefix}_binding_invalid")
    if not path.is_file() or path.stat().st_size != expected_bytes:
        _fail(f"{prefix}_missing_or_size_changed")
    if not index.is_file() or index.stat().st_size != expected_index_bytes:
        _fail(f"{prefix}_index_missing_or_size_changed")
    return path, index


def _parse_fai(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise InputGateFailure("reference_fai_invalid") from exc
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 5:
            _fail("reference_fai_invalid")
        length = _safe_int(fields[1], "reference_fai_invalid")
        if not fields[0] or length < 1:
            _fail("reference_fai_invalid")
        rows.append((fields[0], length))
    if not rows or len(rows) != len({name for name, _ in rows}):
        _fail("reference_fai_invalid")
    return rows


def _parse_dictionary(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise InputGateFailure("reference_dictionary_invalid") from exc
    for line in lines:
        if not line.startswith("@SQ\t"):
            continue
        values = dict(
            field.split(":", maxsplit=1) for field in line.split("\t")[1:] if ":" in field
        )
        if "SN" not in values or "LN" not in values:
            _fail("reference_dictionary_invalid")
        rows.append((values["SN"], _safe_int(values["LN"], "reference_dictionary_invalid")))
    if not rows or len(rows) != len({name for name, _ in rows}):
        _fail("reference_dictionary_invalid")
    return rows


@contextmanager
def _indexed_path(data: Path, index: Path, name: str, index_suffix: str) -> Iterator[Path]:
    """Expose only the declared pair to tools that discover sidecars by filename."""
    with tempfile.TemporaryDirectory() as directory:
        local = Path(directory) / name
        local.symlink_to(data.resolve())
        Path(f"{local}{index_suffix}").symlink_to(index.resolve())
        if index_suffix == ".fai":
            # HTSlib also uses a block-offset index for BGZF references.
            gzi = Path(f"{data}.gzi")
            if gzi.is_file():
                Path(f"{local}.gzi").symlink_to(gzi.resolve())
        yield local


def _alignment_header(alignment: Path) -> tuple[str | None, list[str], list[tuple[str, int]]]:
    completed = _run(["samtools", "view", "-H", str(alignment)])
    if completed.returncode:
        _fail("alignment_header_unreadable")
    sort_order = None
    samples: set[str] = set()
    contigs: list[tuple[str, int]] = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if not fields:
            continue
        if fields[0] == "@HD":
            sort_order = next((field[3:] for field in fields[1:] if field.startswith("SO:")), None)
        elif fields[0] == "@RG":
            sample = next((field[3:] for field in fields[1:] if field.startswith("SM:")), None)
            if sample:
                samples.add(sample)
        elif fields[0] == "@SQ":
            name = next((field[3:] for field in fields[1:] if field.startswith("SN:")), None)
            length = next((field[3:] for field in fields[1:] if field.startswith("LN:")), None)
            if name is None or length is None:
                _fail("alignment_header_invalid")
            contigs.append((name, _safe_int(length, "alignment_header_invalid")))
    return sort_order, sorted(samples), contigs


def _existing_parent(path: Path) -> Path:
    candidate = path.parent
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    if not candidate.exists():
        _fail("output_parent_unavailable")
    return candidate


def _private_bindings(
    manifest_path: Path,
    output_path: Path,
    dataset_id: str,
    dataset: dict[str, Any] | None,
    paths: list[Path],
) -> set[str]:
    identifiers = {dataset_id}
    if dataset is not None:
        identifiers.update(
            {
                str(dataset.get("logical_subject", "")),
                str(dataset.get("laboratory", "")),
            }
        )
        baf = dataset.get("baf_vcf")
        if isinstance(baf, dict):
            identifiers.add(str(baf.get("sample", "")))
    identifiers.discard("")
    bindings = {str(manifest_path), str(output_path), *(str(path) for path in paths)}
    bindings.update(identifiers)
    bindings.update(hashlib.sha256(value.encode()).hexdigest() for value in identifiers)
    bindings.discard("")
    return bindings


def _assert_portable(report: dict[str, Any], forbidden: set[str]) -> None:
    def strings(value: Any):
        if isinstance(value, str):
            yield value
        elif isinstance(value, dict):
            for key, item in value.items():
                yield from strings(key)
                yield from strings(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                yield from strings(item)

    def contains(text: str, binding: str) -> bool:
        if not binding:
            return False
        if len(binding) >= 16 or "/" in binding or "\\" in binding or SHA256.fullmatch(binding):
            return binding in text
        start = 0
        while True:
            index = text.find(binding, start)
            if index < 0:
                return False
            before = text[index - 1] if index else ""
            after_index = index + len(binding)
            after = text[after_index] if after_index < len(text) else ""
            if (not before or not before.isascii() or not before.isalnum()) and (
                not after or not after.isascii() or not after.isalnum()
            ):
                return True
            start = index + 1

    if any(contains(text, binding) for text in strings(report) for binding in forbidden):
        _fail("portable_report_privacy_violation")


def run_input_gate(
    manifest_path: str | Path,
    dataset_id: str,
    caller: str,
    output_path: str | Path,
    *,
    policy: InputGatePolicy = DEFAULT_INPUT_GATE_POLICY,
) -> dict[str, Any]:
    """Validate one caller input and atomically write a path/identifier-free report."""

    manifest_path = Path(manifest_path)
    output_path = Path(output_path)
    if output_path.exists():
        _fail("output_already_exists")
    if not isinstance(dataset_id, str) or not dataset_id:
        _fail("dataset_id_invalid")
    if not isinstance(caller, str) or CALLER_NAME.fullmatch(caller) is None:
        _fail("caller_name_invalid")
    if policy.minimum_free_bytes < 0 or policy.indexed_probe_bases < 1:
        _fail("input_gate_policy_invalid")

    started_at = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    report: dict[str, Any] = {
        "schema_version": 1,
        "status": "RUNNING",
        "started_at": started_at,
        "caller": caller,
    }
    dataset: dict[str, Any] | None = None
    private_paths: list[Path] = []
    forbidden = _private_bindings(manifest_path, output_path, dataset_id, dataset, private_paths)
    try:
        free_bytes = shutil.disk_usage(_existing_parent(output_path)).free
        report["free_bytes_before"] = free_bytes
        report["minimum_free_bytes"] = policy.minimum_free_bytes
        if free_bytes < policy.minimum_free_bytes:
            _fail("free_space_below_hard_stop")

        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise InputGateFailure("input_manifest_unreadable") from exc
        if not isinstance(manifest, dict):
            _fail("input_manifest_invalid")
        manifest_state = manifest.get("status")
        if manifest.get("schema_version") != 1 or manifest_state not in ACCEPTED_MANIFEST_STATES:
            _fail("input_manifest_state_invalid")
        datasets = manifest.get("datasets")
        if not isinstance(datasets, list) or any(not isinstance(row, dict) for row in datasets):
            _fail("input_manifest_datasets_invalid")
        matches = [row for row in datasets if row.get("dataset_id") == dataset_id]
        if len(matches) != 1:
            _fail("dataset_binding_not_unique")
        dataset = matches[0]
        logical_subject = dataset.get("logical_subject")
        if not isinstance(logical_subject, str) or not logical_subject:
            _fail("logical_subject_binding_invalid")

        caller_policy = manifest.get("caller_alignment_policy")
        if not isinstance(caller_policy, dict):
            _fail("caller_alignment_policy_invalid")
        alignment_key = caller_policy.get(caller)
        if alignment_key not in {"original_markdup_bam", "bqsr_cram"}:
            _fail("caller_alignment_policy_invalid")

        reference_spec = manifest.get("reference")
        if not isinstance(reference_spec, dict):
            _fail("reference_binding_invalid")
        assembly = reference_spec.get("assembly")
        if not isinstance(assembly, str) or not assembly.startswith("GRCh38"):
            _fail("reference_assembly_invalid")
        try:
            reference = Path(reference_spec["fasta_path"])
            fai = Path(reference_spec["fai_path"])
            dictionary = Path(reference_spec["dictionary_path"])
        except (KeyError, TypeError) as exc:
            raise InputGateFailure("reference_binding_invalid") from exc
        private_paths.extend((reference, fai, dictionary))
        for path, size_key, digest_key, error_prefix in (
            (reference, "fasta_bytes", "fasta_sha256", "reference_fasta"),
            (fai, "fai_bytes", "fai_sha256", "reference_fai"),
            (dictionary, "dictionary_bytes", "dictionary_sha256", "reference_dictionary"),
        ):
            raw_expected_size = reference_spec.get(size_key)
            expected_size = (
                _safe_int(raw_expected_size, f"{error_prefix}_binding_invalid")
                if raw_expected_size is not None
                else None
            )
            if error_prefix == "reference_fasta" and expected_size is None:
                _fail("reference_fasta_binding_invalid")
            expected_digest = _digest(
                reference_spec.get(digest_key), f"{error_prefix}_binding_invalid"
            )
            if not path.is_file() or (
                expected_size is not None and path.stat().st_size != expected_size
            ):
                _fail(f"{error_prefix}_missing_or_size_changed")
            if _sha256_file(path) != expected_digest:
                _fail(f"{error_prefix}_digest_changed")
        fai_contigs = _parse_fai(fai)
        dictionary_contigs = _parse_dictionary(dictionary)
        if fai_contigs != dictionary_contigs:
            _fail("reference_index_dictionary_boundary_mismatch")

        alignment_spec = dataset.get(alignment_key)
        if not isinstance(alignment_spec, dict):
            _fail("alignment_binding_invalid")
        alignment, alignment_index = _require_bound_file(alignment_spec, prefix="alignment")
        private_paths.extend((alignment, alignment_index))
        alignment_sha256 = _sha256_file(alignment)
        alignment_index_sha256 = _sha256_file(alignment_index)
        if alignment_sha256 != _digest(alignment_spec.get("sha256"), "alignment_binding_invalid"):
            _fail("alignment_digest_changed")
        if alignment_index_sha256 != _digest(
            alignment_spec.get("index_sha256"), "alignment_binding_invalid"
        ):
            _fail("alignment_index_digest_changed")
        if _run(["samtools", "quickcheck", "-vv", str(alignment)]).returncode:
            _fail("alignment_quickcheck_failed")
        sort_order, samples, alignment_contigs = _alignment_header(alignment)
        if sort_order != "coordinate":
            _fail("alignment_not_coordinate_sorted")
        if samples != [logical_subject]:
            _fail("alignment_sample_identity_mismatch")
        if alignment_contigs != fai_contigs:
            _fail("alignment_reference_boundary_mismatch")

        probe_contig, probe_length = next(
            (row for preferred in ("1", "chr1") for row in fai_contigs if row[0] == preferred),
            fai_contigs[0],
        )
        probe_region = f"{probe_contig}:1-{min(probe_length, policy.indexed_probe_bases)}"
        with _indexed_path(reference, fai, "reference.fa", ".fai") as indexed_reference:
            reference_probe = _run(["samtools", "faidx", str(indexed_reference), probe_region])
            reference_probe_lines = reference_probe.stdout.splitlines()
            expected_probe_bases = min(probe_length, policy.indexed_probe_bases)
            if (
                reference_probe.returncode
                or not reference_probe_lines
                or not reference_probe_lines[0].startswith(">")
                or len("".join(reference_probe_lines[1:])) != expected_probe_bases
            ):
                _fail("reference_index_query_failed")
            alignment_probe_command = ["samtools", "view", "-c"]
            if alignment.suffix.lower() == ".cram":
                alignment_probe_command.extend(["-T", str(indexed_reference)])
            alignment_probe_command.extend(
                ["-X", str(alignment), str(alignment_index), probe_region]
            )
            alignment_probe = _run(alignment_probe_command)
        try:
            indexed_probe_records = int(alignment_probe.stdout.strip())
        except ValueError as exc:
            raise InputGateFailure("alignment_index_query_failed") from exc
        if alignment_probe.returncode or indexed_probe_records < 0:
            _fail("alignment_index_query_failed")

        baf_identity = None
        baf_spec = dataset.get("baf_vcf")
        if caller == "cnvpytor":
            if not isinstance(baf_spec, dict):
                _fail("baf_vcf_binding_invalid")
            baf, baf_index = _require_bound_file(baf_spec, prefix="baf_vcf")
            private_paths.extend((baf, baf_index))
            baf_sha256 = _sha256_file(baf)
            baf_index_sha256 = _sha256_file(baf_index)
            if baf_sha256 != _digest(baf_spec.get("sha256"), "baf_vcf_binding_invalid"):
                _fail("baf_vcf_digest_changed")
            if baf_index_sha256 != _digest(baf_spec.get("index_sha256"), "baf_vcf_binding_invalid"):
                _fail("baf_vcf_index_digest_changed")
            declared_baf_sample = baf_spec.get("sample")
            if not isinstance(declared_baf_sample, str) or not declared_baf_sample:
                _fail("baf_vcf_sample_binding_invalid")
            sample_query = _run(["bcftools", "query", "-l", str(baf)])
            observed_samples = sample_query.stdout.splitlines()
            if sample_query.returncode or observed_samples != [declared_baf_sample]:
                _fail("baf_vcf_sample_identity_mismatch")
            # HTSlib detects CSI and TBI by content at the conventional index path.
            with _indexed_path(baf, baf_index, "baf.vcf.gz", ".csi") as indexed_baf:
                count_query = _run(["bcftools", "index", "-n", str(indexed_baf)])
            try:
                baf_records = int(count_query.stdout.strip())
            except ValueError as exc:
                raise InputGateFailure("baf_vcf_index_query_failed") from exc
            if count_query.returncode or baf_records < 1:
                _fail("baf_vcf_index_query_failed")
            baf_identity = {
                "sha256": baf_sha256,
                "index_sha256": baf_index_sha256,
                "sample_count": 1,
                "sample_identity_match": True,
                "records": baf_records,
            }

        identity = {
            "schema_version": 1,
            "caller": caller,
            "alignment_policy": alignment_key,
            "alignment_sha256": alignment_sha256,
            "alignment_index_sha256": alignment_index_sha256,
            "reference_sha256": reference_spec["fasta_sha256"],
            "reference_fai_sha256": reference_spec["fai_sha256"],
            "reference_dictionary_sha256": reference_spec["dictionary_sha256"],
            "baf": baf_identity,
        }
        report.update(
            {
                "status": "PASS",
                "finished_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
                "input_manifest_state": manifest_state,
                "alignment_policy": alignment_key,
                "selected_alignment": {
                    "sha256": alignment_sha256,
                    "index_sha256": alignment_index_sha256,
                    "bytes": alignment.stat().st_size,
                    "sort_order": sort_order,
                    "contigs": len(alignment_contigs),
                    "indexed_probe_records": indexed_probe_records,
                },
                "reference": {
                    "sha256": reference_spec["fasta_sha256"],
                    "fai_sha256": reference_spec["fai_sha256"],
                    "dictionary_sha256": reference_spec["dictionary_sha256"],
                    "contigs": len(fai_contigs),
                    "indexed_probe_succeeded": True,
                },
                "baf_vcf": baf_identity,
                "run_identity_sha256": _sha256_json(identity),
                "logical_subject_count": 1,
                "private_paths_persisted": False,
                "subject_identifier_persisted": False,
                "private_manifest_digest_persisted": False,
            }
        )
        forbidden = _private_bindings(
            manifest_path, output_path, dataset_id, dataset, private_paths
        )
        _assert_portable(report, forbidden)
        report["privacy_scan"] = {
            "status": "PASS",
            "forbidden_private_bindings_checked": len(forbidden),
            "raw_identifier_hashes_persisted": False,
        }
        _assert_portable(report, forbidden)
    except Exception as exc:
        failure = (
            exc if isinstance(exc, InputGateFailure) else InputGateFailure("unexpected_failure")
        )
        report = {
            "schema_version": 1,
            "status": "FAIL_RETAINED",
            "started_at": started_at,
            "finished_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
            "caller": caller,
            "error_code": failure.code,
            "error": "input gate failed; inspect the authorized private process log",
            "private_paths_persisted": False,
            "subject_identifier_persisted": False,
            "private_manifest_digest_persisted": False,
        }
        forbidden = _private_bindings(
            manifest_path, output_path, dataset_id, dataset, private_paths
        )
        try:
            _assert_portable(report, forbidden)
        except InputGateFailure:
            report["error_code"] = "portable_report_privacy_violation"
        _atomic_json(output_path, report)
        if isinstance(exc, InputGateFailure):
            raise
        raise failure from exc
    _atomic_json(output_path, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--caller", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--minimum-free-bytes",
        type=int,
        default=DEFAULT_INPUT_GATE_POLICY.minimum_free_bytes,
    )
    parser.add_argument(
        "--indexed-probe-bases",
        type=int,
        default=DEFAULT_INPUT_GATE_POLICY.indexed_probe_bases,
    )
    args = parser.parse_args()
    try:
        report = run_input_gate(
            args.manifest,
            args.dataset_id,
            args.caller,
            args.output,
            policy=InputGatePolicy(
                minimum_free_bytes=args.minimum_free_bytes,
                indexed_probe_bases=args.indexed_probe_bases,
            ),
        )
    except InputGateFailure as exc:
        print(json.dumps({"status": "FAIL_RETAINED", "error_code": exc.code}, sort_keys=True))
        return 2
    print(
        json.dumps(
            {
                "status": report["status"],
                "caller": report["caller"],
                "alignment_policy": report["alignment_policy"],
                "run_identity_sha256": report["run_identity_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
