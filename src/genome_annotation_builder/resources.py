"""Checksum-pinned synchronization of normative upstream resources."""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .config import BuildConfig
from .constants import (
    GENOME_SPEC_COMMIT,
    GENOME_SPEC_RAW_ROOT,
    GENOME_SPEC_REPOSITORY,
    PINNED_UPSTREAM_FILES,
)
from .errors import ResourceError
from .util import atomic_write_bytes, atomic_write_json, sha256_bytes, sha256_file


def genome_spec_root(config: BuildConfig) -> Path:
    return config.resource_cache / "genome-spec" / GENOME_SPEC_COMMIT


def resource_path(config: BuildConfig, relative_path: str) -> Path:
    return genome_spec_root(config) / relative_path


def _vendored_path(relative_path: str) -> Path:
    return Path(__file__).with_name("vendor") / "genome_spec" / relative_path


def _download(relative_path: str, timeout_seconds: int = 60) -> bytes:
    url = f"{GENOME_SPEC_RAW_ROOT}/{relative_path}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "genome-annotation-builder-resource-sync/0.1"},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            if response.status != 200:
                raise ResourceError(f"unexpected HTTP status {response.status} for {url}")
            return response.read()
    except (urllib.error.URLError, TimeoutError) as exc:
        raise ResourceError(f"failed to download pinned resource {url}: {exc}") from exc


def verify_pinned_resources(config: BuildConfig) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    errors: list[str] = []
    for relative_path, expected in PINNED_UPSTREAM_FILES.items():
        path = resource_path(config, relative_path)
        if not path.is_file():
            results[relative_path] = {"status": "missing", "expected_sha256": expected}
            errors.append(f"missing: {relative_path}")
            continue
        actual = sha256_file(path)
        status = "verified" if actual == expected else "digest_mismatch"
        results[relative_path] = {
            "status": status,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "bytes": path.stat().st_size,
        }
        if status != "verified":
            errors.append(f"digest mismatch: {relative_path}")
    if errors:
        raise ResourceError(
            "pinned genome-spec resources failed verification: " + "; ".join(errors)
        )
    return results


def sync_pinned_resources(config: BuildConfig) -> dict[str, Any]:
    root = genome_spec_root(config)
    root.mkdir(parents=True, exist_ok=True)
    changed: list[str] = []
    reused: list[str] = []
    for relative_path, expected in PINNED_UPSTREAM_FILES.items():
        path = root / relative_path
        if path.is_file() and sha256_file(path) == expected:
            reused.append(relative_path)
            continue
        vendored = _vendored_path(relative_path)
        if vendored.is_file() and sha256_file(vendored) == expected:
            atomic_write_bytes(path, vendored.read_bytes())
            changed.append(relative_path)
            continue
        if not config.allow_download:
            raise ResourceError(
                f"resource {relative_path} is absent or invalid and "
                "resources.allow_download is false"
            )
        payload = _download(relative_path)
        actual = sha256_bytes(payload)
        if actual != expected:
            raise ResourceError(
                f"downloaded resource digest mismatch for {relative_path}: "
                f"expected {expected}, got {actual}"
            )
        atomic_write_bytes(path, payload)
        changed.append(relative_path)

    verified = verify_pinned_resources(config)
    manifest = {
        "resource_set": "genome-spec",
        "repository": GENOME_SPEC_REPOSITORY,
        "commit": GENOME_SPEC_COMMIT,
        "verified_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "changed": changed,
        "reused": reused,
        "files": verified,
    }
    atomic_write_json(root / "resource_manifest.json", manifest)
    return manifest


def load_normative_schema(config: BuildConfig) -> dict[str, Any]:
    schema_path = resource_path(config, "schema.json")
    expected = PINNED_UPSTREAM_FILES["schema.json"]
    if not schema_path.is_file() or sha256_file(schema_path) != expected:
        raise ResourceError("normative schema is missing or invalid; run `gab resources sync`")
    try:
        return json.loads(schema_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ResourceError(f"normative schema is invalid JSON: {exc}") from exc
