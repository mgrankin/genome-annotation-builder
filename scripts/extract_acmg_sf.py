#!/usr/bin/env python3
"""Regenerate the deterministic vendored ACMG/ClinGen snapshot."""

from __future__ import annotations

import argparse
from pathlib import Path

from genome_annotation_builder.acmg_resources import extract_acmg_workbook
from genome_annotation_builder.util import atomic_write_json, sha256_file


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    snapshot = extract_acmg_workbook(arguments.workbook)
    atomic_write_json(arguments.output, snapshot)
    print(f"{sha256_file(arguments.output)}  {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
