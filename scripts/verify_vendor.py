#!/usr/bin/env python3
"""Verify byte identity and structure of every vendored public resource."""

from pathlib import Path

from genome_annotation_builder.acmg_resources import verify_vendor_acmg_snapshot
from genome_annotation_builder.clingen_resources import verify_vendor_clingen_snapshot
from genome_annotation_builder.constants import (
    LOFTEE_PATCH_FILE,
    LOFTEE_PATCH_SHA256,
    PINNED_UPSTREAM_FILES,
)
from genome_annotation_builder.util import sha256_file


def main() -> int:
    root = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "genome_annotation_builder"
        / "vendor"
        / "genome_spec"
    )
    for relative, expected in PINNED_UPSTREAM_FILES.items():
        path = root / relative
        actual = sha256_file(path)
        if actual != expected:
            raise SystemExit(f"digest mismatch: {relative}: {actual}")
    patch = root.parent.parent / "vendor_patches" / LOFTEE_PATCH_FILE
    patch_digest = sha256_file(patch)
    if patch_digest != LOFTEE_PATCH_SHA256:
        raise SystemExit(f"digest mismatch: vendor_patches/{LOFTEE_PATCH_FILE}: {patch_digest}")
    acmg = verify_vendor_acmg_snapshot()
    clingen = verify_vendor_clingen_snapshot()
    print(
        "vendored resource verification PASS: "
        f"{len(PINNED_UPSTREAM_FILES)} genome-spec file(s), "
        f"ACMG SF v{acmg['version']} ({acmg['counts']['genes']} genes), "
        f"ClinGen {clingen['snapshot']} "
        f"({clingen['gene_disease_rows']} gene-disease + "
        f"{clingen['dosage_rows']} dosage rows), "
        "1 LOFTEE patch"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
