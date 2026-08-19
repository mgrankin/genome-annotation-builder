"""Small deterministic filesystem and hashing primitives."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any, BinaryIO

CHUNK_SIZE = 8 * 1024 * 1024
UINT256_MODULUS = 1 << 256


class MultisetDigest:
    """Order-independent cryptographic accounting for a multiset of text IDs."""

    def __init__(self) -> None:
        self.count = 0
        self._xor = 0
        self._sum = 0

    def add(self, value: str) -> None:
        number = int.from_bytes(hashlib.sha256(value.encode("utf-8")).digest(), "big")
        self.count += 1
        self._xor ^= number
        self._sum = (self._sum + number) % UINT256_MODULUS

    def as_dict(self) -> dict[str, int | str]:
        return {
            "count": self.count,
            "sha256_xor": f"{self._xor:064x}",
            "sha256_sum_mod_2_256": f"{self._sum:064x}",
        }


def sha256_stream(handle: BinaryIO) -> str:
    digest = hashlib.sha256()
    while chunk := handle.read(CHUNK_SIZE):
        digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: str | Path) -> str:
    with Path(path).open("rb") as handle:
        return sha256_stream(handle)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json_bytes(value))


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"
    ).encode()


def directory_digest(path: str | Path) -> str:
    """Hash relative names, NUL separators, and bytes in lexical order."""
    root = Path(path)
    digest = hashlib.sha256()
    for item in sorted(p for p in root.rglob("*") if p.is_file()):
        digest.update(item.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        with item.open("rb") as handle:
            while chunk := handle.read(CHUNK_SIZE):
                digest.update(chunk)
    return digest.hexdigest()


def source_tree_digest(path: str | Path) -> str:
    """Hash shipped implementation bytes while excluding runtime bytecode caches."""

    root = Path(path)
    digest = hashlib.sha256()
    files = sorted(
        item
        for item in root.rglob("*")
        if item.is_file()
        and "__pycache__" not in item.relative_to(root).parts
        and item.suffix not in {".pyc", ".pyo"}
    )
    for item in files:
        if item.is_symlink():
            raise ValueError(f"implementation tree contains a symbolic link: {item}")
        digest.update(item.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        with item.open("rb") as handle:
            while chunk := handle.read(CHUNK_SIZE):
                digest.update(chunk)
    return digest.hexdigest()


def fsync_path(path: str | Path) -> None:
    descriptor = os.open(Path(path), os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write_bytes(path: str | Path, payload: bytes, mode: int = 0o644) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.partial.", dir=destination.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, destination)
        fsync_path(destination.parent)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_json(path: str | Path, value: Any) -> None:
    atomic_write_bytes(path, canonical_json_bytes(value))


@contextmanager
def atomic_output_directory(destination: str | Path, overwrite: bool = False) -> Iterator[Path]:
    """Yield a sibling staging directory and promote it atomically on success."""
    final = Path(destination)
    final.parent.mkdir(parents=True, exist_ok=True)
    if final.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing output: {final}")
    staging = Path(tempfile.mkdtemp(prefix=f".{final.name}.partial.", dir=final.parent))
    try:
        yield staging
        fsync_path(staging)
        if final.exists():
            # The caller opted in to overwrite. Rename the old output aside so
            # promotion is recoverable; cleanup is a separate explicit action.
            previous = final.with_name(f".{final.name}.previous")
            if previous.exists():
                raise FileExistsError(f"previous-output guard already exists: {previous}")
            os.replace(final, previous)
        os.replace(staging, final)
        fsync_path(final.parent)
    except BaseException:
        # Preserve a populated staging directory for forensic recovery.
        if staging.exists() and not any(staging.iterdir()):
            staging.rmdir()
        raise
