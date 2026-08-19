#!/usr/bin/env python3
"""Fail closed when release archives contain private or generated material."""

from __future__ import annotations

import argparse
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath

if __package__:
    from .privacy_scan import deny_values, scan_payload
else:
    from privacy_scan import deny_values, scan_payload

REQUIRED_SUFFIXES = (
    "vendor/clingen/clingen_gene_level_20260815.json",
    "vendor_patches/loftee_biodbhts_faidx_v1.patch",
    "release_c/profiles/cnvpytor-combined-rd-baf-1.3.2.json",
    "release_c/profiles/delly-short-read-2.5.1.json",
    "release_c/profiles/expansionhunter-grch38-female-5.0.0.json",
    "release_c/profiles/expansionhunter-grch38-male-5.0.0.json",
    "release_c/profiles/gatk-mitochondria-normal-shifted-4.6.2.0.json",
    "release_c/profiles/hla-la-short-read-1.0.4.json",
    "release_c/profiles/manta-germline-1.6.0.json",
    "release_c/profiles/readback-bounded-pysam-0.24.0.json",
    "release_c/profiles/spechla-short-read-bam-1.0.12.json",
)
FORBIDDEN_PARTS = {
    "__pycache__",
    "private",
    "inputs",
    "outputs",
    "resources",
    "results",
}
FORBIDDEN_SUFFIXES = (
    ".bam",
    ".bai",
    ".cram",
    ".crai",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".gvcf",
    ".gvcf.gz",
    ".parquet",
    ".pyc",
    ".pyo",
    ".sqlite",
    ".db",
    ".bw",
    ".bigwig",
    ".zip",
    ".tar",
    ".tar.gz",
    ".xlsx",
    ".xls",
    ".vcf",
    ".vcf.gz",
)


def _logical_name(name: str) -> str:
    """Remove the conventional sdist root without weakening path checks."""

    parts = PurePosixPath(name).parts
    if len(parts) > 1 and parts[0].startswith("genome_annotation_builder-"):
        return PurePosixPath(*parts[1:]).as_posix()
    return name


def _synthetic_fixture(name: str) -> bool:
    path = PurePosixPath(name)
    return (
        path.parent == PurePosixPath("tests/fixtures")
        and path.name.startswith("synthetic")
        and path.suffix == ".vcf"
    )


def _members(path: Path) -> list[tuple[str, bytes]]:
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            members: list[tuple[str, bytes]] = []
            for info in archive.infolist():
                mode = (info.external_attr >> 16) & 0o170000
                if mode == stat.S_IFLNK:
                    raise SystemExit(f"symbolic link in wheel {path}: {info.filename}")
                if not info.is_dir():
                    members.append((info.filename, archive.read(info)))
            return members
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path, "r:gz") as archive:
            members = []
            for info in archive.getmembers():
                if info.isdir():
                    continue
                if not info.isfile():
                    raise SystemExit(
                        f"non-regular member in source distribution {path}: {info.name}"
                    )
                extracted = archive.extractfile(info)
                if extracted is None:
                    raise SystemExit(f"cannot read source-distribution member: {info.name}")
                members.append((info.name, extracted.read()))
            return members
    raise ValueError(f"unsupported distribution: {path}")


def _validate(path: Path, denied: list[bytes]) -> None:
    members = _members(path)
    if not members:
        raise SystemExit(f"empty distribution: {path}")

    bad: list[str] = []
    privacy_errors: list[str] = []
    for name, payload in members:
        member = PurePosixPath(name)
        logical_name = _logical_name(name)
        lowered = name.lower()
        if member.is_absolute() or ".." in member.parts:
            bad.append(name)
            continue
        if set(part.lower() for part in member.parts) & FORBIDDEN_PARTS:
            bad.append(name)
            continue
        if (
            lowered.endswith(FORBIDDEN_SUFFIXES) and not _synthetic_fixture(logical_name)
        ) or member.name == "TASK.md":
            bad.append(name)
        privacy_errors.extend(
            scan_payload(
                relative=logical_name,
                payload=payload,
                size=len(payload),
                denied=denied,
            )
        )

    names = [name for name, _payload in members]
    missing = [
        suffix for suffix in REQUIRED_SUFFIXES if not any(name.endswith(suffix) for name in names)
    ]
    if bad or missing or privacy_errors:
        details = []
        if bad:
            details.append("forbidden members: " + ", ".join(sorted(bad)))
        if missing:
            details.append("missing required members: " + ", ".join(missing))
        if privacy_errors:
            details.append("privacy errors: " + ", ".join(sorted(privacy_errors)))
        raise SystemExit(f"distribution validation failed for {path}: " + "; ".join(details))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", nargs="?", default="dist", type=Path)
    parser.add_argument("--deny-file", type=Path)
    args = parser.parse_args()
    denied = deny_values(args.deny_file)

    wheels = sorted(args.directory.glob("*.whl"))
    sdists = sorted(args.directory.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise SystemExit(
            f"expected exactly one wheel and one sdist in {args.directory}; "
            f"found {len(wheels)} wheel(s), {len(sdists)} sdist(s)"
        )
    for path in (*wheels, *sdists):
        _validate(path, denied)
        print(f"PASS\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
