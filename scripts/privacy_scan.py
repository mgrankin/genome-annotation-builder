#!/usr/bin/env python3
"""Fail closed on accidental genomic data, secrets, or local paths in Git."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

FORBIDDEN_SUFFIXES = (
    ".bam",
    ".bai",
    ".cram",
    ".crai",
    ".sam",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".gvcf",
    ".gvcf.gz",
    ".bcf",
    ".vcf.gz",
    ".parquet",
    ".duckdb",
    ".sqlite",
    ".db",
    ".bw",
    ".bigwig",
    ".zip",
    ".tar",
    ".tar.gz",
    ".xlsx",
    ".xls",
    ".genome.zip",
)
SYNTHETIC_VCF = re.compile(r"^tests/fixtures/synthetic[^/]*\.vcf$")
SUSPICIOUS_TEXT = {
    "absolute Linux data path": re.compile(rb"/(?:raid|mnt|media|home)/[^\s\"']+"),
    "absolute macOS user path": re.compile(rb"/" + rb"Users/[^\s\"']+"),
    "private key": re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "GitHub token": re.compile(rb"gh[pousr]_[A-Za-z0-9_]{20,}"),
    "AWS access key": re.compile(rb"AKIA[0-9A-Z]{16}"),
    "credential assignment": re.compile(
        rb"(?i)(?:password|passwd|api[_-]?key|access[_-]?token)\s*[:=]\s*[^\s$<{][^\s]{5,}"
    ),
    "share URL with token": re.compile(rb"https?://[^\s]+(?:ssid|token|key)=[^\s&]+"),
}
MAX_TRACKED_BYTES = 5 * 1024 * 1024


@dataclass(frozen=True)
class GitIndexEntry:
    relative: str
    mode: str
    object_id: str


def tracked_entries(root: Path) -> list[GitIndexEntry]:
    completed = subprocess.run(
        ["git", "ls-files", "--stage", "-z"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    entries: list[GitIndexEntry] = []
    for item in completed.stdout.split(b"\0"):
        if not item:
            continue
        metadata, raw_path = item.split(b"\t", maxsplit=1)
        mode, object_id, stage = metadata.decode("ascii").split()
        if stage != "0":
            raise RuntimeError("privacy scan refuses an unmerged Git index")
        entries.append(
            GitIndexEntry(
                relative=raw_path.decode("utf-8"),
                mode=mode,
                object_id=object_id,
            )
        )
    return entries


def deny_values(path: Path | None) -> list[bytes]:
    if path is None or not path.exists():
        return []
    return [
        line.strip().encode("utf-8")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def scan_payload(*, relative: str, payload: bytes, size: int, denied: list[bytes]) -> list[str]:
    errors: list[str] = []
    lower = relative.lower()
    if lower.endswith(FORBIDDEN_SUFFIXES) and not SYNTHETIC_VCF.fullmatch(relative):
        errors.append(f"forbidden genomic/bundle extension: {relative}")
    if size > MAX_TRACKED_BYTES:
        errors.append(f"tracked file exceeds {MAX_TRACKED_BYTES} bytes: {relative}")
        return errors
    if b"\0" in payload:
        errors.append(f"binary payload is not allowed in public source: {relative}")
        return errors
    for label, pattern in SUSPICIOUS_TEXT.items():
        if pattern.search(payload):
            errors.append(f"{label} found in {relative}")
    for index, value in enumerate(denied, start=1):
        if value in payload:
            errors.append(f"local deny-list entry #{index} found in {relative}")
    return errors


def scan_paths(paths: list[Path], root: Path, denied: list[bytes]) -> list[str]:
    errors: list[str] = []
    for path in paths:
        try:
            relative = path.relative_to(root).as_posix()
        except ValueError:
            errors.append("scan path is outside the repository")
            continue
        if path.is_symlink():
            errors.append(f"symbolic link supplied for scan: {relative}")
        elif not path.is_file():
            errors.append(f"non-regular scan path: {relative}")
        else:
            payload = path.read_bytes()
            errors.extend(
                scan_payload(
                    relative=relative,
                    payload=payload,
                    size=len(payload),
                    denied=denied,
                )
            )
    return errors


def scan_index(root: Path, denied: list[bytes]) -> tuple[list[str], int]:
    errors: list[str] = []
    entries = tracked_entries(root)
    for entry in entries:
        if entry.mode == "120000":
            errors.append(f"symbolic link is tracked: {entry.relative}")
            continue
        if entry.mode not in {"100644", "100755"}:
            errors.append(f"unsupported Git mode {entry.mode}: {entry.relative}")
            continue
        size = int(
            subprocess.run(
                ["git", "cat-file", "-s", entry.object_id],
                cwd=root,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
        payload = b""
        if size <= MAX_TRACKED_BYTES:
            payload = subprocess.run(
                ["git", "cat-file", "blob", entry.object_id],
                cwd=root,
                check=True,
                capture_output=True,
            ).stdout
        errors.extend(
            scan_payload(
                relative=entry.relative,
                payload=payload,
                size=size,
                denied=denied,
            )
        )
    return errors, len(entries)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tracked", action="store_true", help="scan the Git index")
    parser.add_argument("--deny-file", type=Path)
    parser.add_argument("paths", nargs="*", type=Path)
    arguments = parser.parse_args(argv)
    root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    )
    if arguments.tracked:
        errors, count = scan_index(root, deny_values(arguments.deny_file))
    else:
        paths = [(root / path).resolve() for path in arguments.paths]
        errors = scan_paths(paths, root, deny_values(arguments.deny_file))
        count = len(paths)
    if errors:
        print("privacy scan failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1
    scope = "Git-index" if arguments.tracked else "working-tree"
    print(f"privacy scan PASS: {count} {scope} file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
