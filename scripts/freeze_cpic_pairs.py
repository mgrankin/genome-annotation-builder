#!/usr/bin/env python3
"""Freeze the public CPIC pair_view into a deterministic, checksum-ready JSON file."""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

from genome_annotation_builder.cpic_freeze import freeze_cpic_pairs
from genome_annotation_builder.util import atomic_write_bytes, canonical_json_bytes, sha256_bytes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--snapshot-date", required=True, type=date.fromisoformat)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    document = freeze_cpic_pairs(arguments.source_version, arguments.snapshot_date)
    payload = canonical_json_bytes(document)
    atomic_write_bytes(arguments.output, payload)
    print(f"{sha256_bytes(payload)}  {arguments.output}")
    print(f"rows={document['row_count']} source_version={arguments.source_version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
