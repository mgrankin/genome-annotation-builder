"""Deterministic exact-member manifests for read-only caller resource trees."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

from .caller_receipt import canonical_json_sha256
from .portable_model import (
    FORBIDDEN_PORTABLE_TEXT,
    SHA256,
    PortableEvidenceError,
    contains_forbidden_binding,
    sha256_file,
)

TREE_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class ValidatedResourceTree:
    manifest_sha256: str
    manifest_bytes: int
    tree_identity_sha256: str
    files: int
    directories: int
    logical_bytes: int
    normalized_manifest: dict[str, Any] = field(repr=False)
    signatures: dict[Path, tuple[int, int, int, int, int]] = field(repr=False)

    def portable_summary(self) -> dict[str, Any]:
        return {
            "schema_version": TREE_SCHEMA_VERSION,
            "status": "PASS",
            "manifest_sha256": self.manifest_sha256,
            "manifest_bytes": self.manifest_bytes,
            "tree_identity_sha256": self.tree_identity_sha256,
            "files": self.files,
            "directories": self.directories,
            "logical_bytes": self.logical_bytes,
            "private_paths_persisted": False,
        }


def _fail(message: str) -> None:
    raise PortableEvidenceError(message)


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        _fail(label)
    return value


def _integer(value: Any, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _fail(label)
    return value


def _regular_root(path: Path) -> None:
    if not path.is_absolute() or path.is_symlink() or not path.is_dir():
        _fail("resource tree root must be an absolute regular directory")
    try:
        if path.resolve(strict=True) != path:
            _fail("resource tree root contains a symbolic path component")
    except OSError as exc:
        raise PortableEvidenceError("resource tree root is unreadable") from exc


def _regular_manifest(path: Path) -> None:
    if not path.is_absolute() or path.is_symlink() or not path.is_file():
        _fail("resource tree manifest must be an absolute regular file")
    try:
        if path.resolve(strict=True) != path:
            _fail("resource tree manifest contains a symbolic path component")
    except OSError as exc:
        raise PortableEvidenceError("resource tree manifest is unreadable") from exc


def _relative_member(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        _fail(label)
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or str(path) != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        _fail(label)
    return value


def _stat_signature(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def _inventory(
    root: Path,
) -> tuple[dict[str, Path], dict[str, Path], dict[Path, tuple[int, int, int, int, int]]]:
    files: dict[str, Path] = {}
    directories: dict[str, Path] = {}
    signatures = {root: _stat_signature(root)}
    try:
        for path in root.rglob("*"):
            relative = path.relative_to(root).as_posix()
            if path.is_symlink():
                _fail("resource tree contains a symbolic link")
            if path.is_file():
                files[relative] = path
            elif path.is_dir():
                directories[relative] = path
            else:
                _fail("resource tree contains a special file")
            signatures[path] = _stat_signature(path)
    except OSError as exc:
        raise PortableEvidenceError("resource tree inventory is unreadable") from exc
    if not files:
        _fail("resource tree contains no regular files")
    return files, directories, signatures


def _assert_private_free(value: Any, forbidden_bindings: Sequence[str]) -> None:
    forbidden = tuple(binding for binding in forbidden_bindings if binding)

    def strings(item: Any):
        if isinstance(item, str):
            yield item
        elif isinstance(item, dict):
            for key, child in item.items():
                yield from strings(key)
                yield from strings(child)
        elif isinstance(item, (list, tuple)):
            for child in item:
                yield from strings(child)

    for text in strings(value):
        if FORBIDDEN_PORTABLE_TEXT.search(text):
            _fail("resource tree manifest contains a private path or secret")
        if any(contains_forbidden_binding(text, binding) for binding in forbidden):
            _fail("resource tree manifest contains a forbidden private binding")


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    if path.exists():
        _fail("refusing to replace an existing resource tree manifest")
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
            raise PortableEvidenceError(
                "refusing to replace an existing resource tree manifest"
            ) from exc
        temporary.unlink()
        _fsync_directory(path.parent)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _fsync_directory(path: Path) -> None:
    """Persist directory-entry changes after an atomic no-replace publish."""

    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _output_outside_root(root: Path, output: Path) -> None:
    absolute_output = output.absolute()
    if absolute_output == root or root in absolute_output.parents:
        _fail("resource tree manifest must be written outside the resource root")


def create_resource_tree_manifest(
    root: str | Path,
    output_path: str | Path,
    *,
    forbidden_bindings: Sequence[str] = (),
) -> dict[str, Any]:
    """Hash one stable tree and atomically publish a path-free exact inventory."""

    root_path = Path(root)
    output = Path(output_path)
    _regular_root(root_path)
    if not output.is_absolute():
        _fail("resource tree manifest output must be an absolute path")
    if output.exists():
        _fail("refusing to replace an existing resource tree manifest")
    _output_outside_root(root_path, output)
    files, directories, initial_signatures = _inventory(root_path)
    file_rows: list[dict[str, Any]] = []
    logical_bytes = 0
    for relative, path in sorted(files.items()):
        before = _stat_signature(path)
        digest = sha256_file(path)
        after = _stat_signature(path)
        if before != after or initial_signatures[path] != after:
            _fail("resource tree changed while it was being hashed")
        byte_count = path.stat().st_size
        logical_bytes += byte_count
        file_rows.append({"path": relative, "bytes": byte_count, "sha256": digest})
    directory_rows = sorted(directories)
    final_files, final_directories, final_signatures = _inventory(root_path)
    if set(final_files) != set(files) or set(final_directories) != set(directories):
        _fail("resource tree inventory changed while it was being hashed")
    if any(final_signatures[path] != signature for path, signature in initial_signatures.items()):
        _fail("resource tree metadata changed while it was being hashed")
    identity_payload = {
        "schema_version": TREE_SCHEMA_VERSION,
        "directories": directory_rows,
        "files": file_rows,
    }
    manifest = {
        "schema_version": TREE_SCHEMA_VERSION,
        "status": "PASS",
        "directories": directory_rows,
        "files": file_rows,
        "tree_identity_sha256": canonical_json_sha256(identity_payload),
    }
    _assert_private_free(
        manifest,
        (*forbidden_bindings, str(root_path), str(output)),
    )
    _atomic_json(output, manifest)
    return {
        "schema_version": TREE_SCHEMA_VERSION,
        "status": "PASS",
        "manifest_sha256": sha256_file(output),
        "manifest_bytes": output.stat().st_size,
        "tree_identity_sha256": manifest["tree_identity_sha256"],
        "files": len(file_rows),
        "directories": len(directory_rows),
        "logical_bytes": logical_bytes,
        "private_paths_persisted": False,
    }


def validate_resource_tree_manifest(
    root: str | Path,
    manifest_path: str | Path,
    *,
    expected_manifest_sha256: str | None = None,
    expected_manifest_bytes: int | None = None,
    forbidden_bindings: Sequence[str] = (),
) -> ValidatedResourceTree:
    """Independently re-hash every declared member and reject any undeclared entry."""

    root_path = Path(root)
    path = Path(manifest_path)
    _regular_root(root_path)
    _regular_manifest(path)
    _output_outside_root(root_path, path)
    observed_manifest_sha256 = sha256_file(path)
    observed_manifest_bytes = path.stat().st_size
    if expected_manifest_sha256 is not None and observed_manifest_sha256 != _sha(
        expected_manifest_sha256, "expected resource tree manifest digest is invalid"
    ):
        _fail("resource tree manifest digest changed")
    if expected_manifest_bytes is not None and observed_manifest_bytes != _integer(
        expected_manifest_bytes, "expected resource tree manifest byte count is invalid", minimum=1
    ):
        _fail("resource tree manifest byte count changed")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError("resource tree manifest is unreadable") from exc
    if not isinstance(manifest, dict) or set(manifest) != {
        "schema_version",
        "status",
        "directories",
        "files",
        "tree_identity_sha256",
    }:
        _fail("resource tree manifest schema is invalid")
    if manifest.get("schema_version") != TREE_SCHEMA_VERSION or manifest.get("status") != "PASS":
        _fail("resource tree manifest state is invalid")
    raw_directories = manifest.get("directories")
    raw_files = manifest.get("files")
    if not isinstance(raw_directories, list) or any(
        not isinstance(value, str) for value in raw_directories
    ):
        _fail("resource tree directory inventory is invalid")
    directories = [
        _relative_member(value, "resource tree directory member is invalid")
        for value in raw_directories
    ]
    if directories != sorted(set(directories)):
        _fail("resource tree directory inventory is not canonical")
    if not isinstance(raw_files, list) or not raw_files:
        _fail("resource tree file inventory is invalid")
    files: list[dict[str, Any]] = []
    seen_files: set[str] = set()
    for item in raw_files:
        if not isinstance(item, dict) or set(item) != {"path", "bytes", "sha256"}:
            _fail("resource tree file member is invalid")
        relative = _relative_member(item.get("path"), "resource tree file member is invalid")
        if relative in seen_files:
            _fail("resource tree file member is duplicated")
        seen_files.add(relative)
        files.append(
            {
                "path": relative,
                "bytes": _integer(item.get("bytes"), "resource tree file member is invalid"),
                "sha256": _sha(item.get("sha256"), "resource tree file member is invalid"),
            }
        )
    if files != sorted(files, key=lambda item: item["path"]):
        _fail("resource tree file inventory is not canonical")
    identity_payload = {
        "schema_version": TREE_SCHEMA_VERSION,
        "directories": directories,
        "files": files,
    }
    expected_tree_identity = canonical_json_sha256(identity_payload)
    if manifest.get("tree_identity_sha256") != expected_tree_identity:
        _fail("resource tree identity changed")
    _assert_private_free(
        manifest,
        (*forbidden_bindings, str(root_path), str(path)),
    )
    actual_files, actual_directories, signatures = _inventory(root_path)
    if set(actual_files) != seen_files or set(actual_directories) != set(directories):
        _fail("resource tree exact inventory changed")
    logical_bytes = 0
    for item in files:
        member = actual_files[item["path"]]
        before = _stat_signature(member)
        digest = sha256_file(member)
        after = _stat_signature(member)
        if before != after:
            _fail("resource tree changed while it was being validated")
        if after[2] != item["bytes"] or digest != item["sha256"]:
            _fail("resource tree member digest or byte count changed")
        signatures[member] = after
        logical_bytes += int(item["bytes"])
    signatures[path] = _stat_signature(path)
    return ValidatedResourceTree(
        manifest_sha256=observed_manifest_sha256,
        manifest_bytes=observed_manifest_bytes,
        tree_identity_sha256=expected_tree_identity,
        files=len(files),
        directories=len(directories),
        logical_bytes=logical_bytes,
        normalized_manifest={
            "schema_version": TREE_SCHEMA_VERSION,
            "status": "PASS",
            "directories": directories,
            "files": files,
            "tree_identity_sha256": expected_tree_identity,
        },
        signatures=signatures,
    )


def _read_forbidden_bindings(path: Path | None) -> tuple[str, ...]:
    if path is None:
        return ()
    if path.is_symlink() or not path.is_file():
        _fail("forbidden-binding file is unreadable")
    try:
        return tuple(
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    except (OSError, UnicodeError) as exc:
        raise PortableEvidenceError("forbidden-binding file is unreadable") from exc


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    create_parser = subparsers.add_parser("create")
    create_parser.add_argument("--root", type=Path, required=True)
    create_parser.add_argument("--output", type=Path, required=True)
    create_parser.add_argument("--forbidden-bindings-file", type=Path)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--root", type=Path, required=True)
    validate_parser.add_argument("--manifest", type=Path, required=True)
    validate_parser.add_argument("--expected-manifest-sha256")
    validate_parser.add_argument("--expected-manifest-bytes", type=int)
    validate_parser.add_argument("--forbidden-bindings-file", type=Path)
    args = parser.parse_args()
    try:
        forbidden = _read_forbidden_bindings(args.forbidden_bindings_file)
        if args.command == "create":
            report = create_resource_tree_manifest(
                args.root,
                args.output,
                forbidden_bindings=(
                    *forbidden,
                    *((str(args.forbidden_bindings_file),) if args.forbidden_bindings_file else ()),
                ),
            )
        else:
            validated = validate_resource_tree_manifest(
                args.root,
                args.manifest,
                expected_manifest_sha256=args.expected_manifest_sha256,
                expected_manifest_bytes=args.expected_manifest_bytes,
                forbidden_bindings=(
                    *forbidden,
                    *((str(args.forbidden_bindings_file),) if args.forbidden_bindings_file else ()),
                ),
            )
            report = validated.portable_summary()
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    except Exception:
        print(
            "resource tree operation failed; inspect the authorized private process log",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
