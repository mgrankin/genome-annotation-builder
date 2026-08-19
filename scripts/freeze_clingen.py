#!/usr/bin/env python3
"""Freeze official ClinGen gene-level sources into a canonical public snapshot."""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

from genome_annotation_builder.clingen_freeze import freeze_clingen_snapshot, snapshot_identity
from genome_annotation_builder.util import atomic_write_bytes, canonical_json_bytes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-disease-csv", required=True, type=Path)
    parser.add_argument("--dosage-tsv", required=True, type=Path)
    parser.add_argument("--gene-disease-snapshot", required=True, type=date.fromisoformat)
    parser.add_argument("--dosage-snapshot", required=True, type=date.fromisoformat)
    parser.add_argument("--accessed-at", required=True, type=date.fromisoformat)
    parser.add_argument("--gene-disease-url", required=True)
    parser.add_argument("--dosage-url", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    document = freeze_clingen_snapshot(
        arguments.gene_disease_csv,
        arguments.dosage_tsv,
        gene_disease_snapshot=arguments.gene_disease_snapshot,
        dosage_snapshot=arguments.dosage_snapshot,
        accessed_at=arguments.accessed_at,
        gene_disease_url=arguments.gene_disease_url,
        dosage_url=arguments.dosage_url,
    )
    payload = canonical_json_bytes(document)
    atomic_write_bytes(arguments.output, payload)
    identity = snapshot_identity(document)
    print(f"{identity['canonical_sha256']}  {arguments.output}")
    print(
        f"gene_disease_rows={identity['gene_disease_rows']} "
        f"dosage_rows={identity['dosage_rows']} snapshot={identity['snapshot']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
