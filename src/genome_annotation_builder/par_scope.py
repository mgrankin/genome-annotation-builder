"""Import the family VCF producer's negative-query limits into portable coordinates."""

from __future__ import annotations

import re
from typing import Any

import pysam

from .errors import InputValidationError

TAG = "family_reference_confidence_par_query"
MASKED = "GRCh38-Ensembl-Y-PAR-hardmask"
UNMASKED = "GRCh38-Ensembl-unmasked-XY"
DECLARED = "DECLARED_PAR_MAPPING_UNASSESSED"
GRCH38 = {
    "X": (156040895, ((10000, 2781479), (155701382, 156030895))),
    "Y": (57227415, ((10000, 2781479), (56887902, 57217415))),
}
PINS = ("CanonicalSHA256", "MappingSHA256", "PARBedSHA256")


def _compatible_family_intervals(header: pysam.VariantHeader) -> list[dict[str, Any]]:
    # Length agreement supports a coordinate-only restriction. It does not identify a
    # mapping procedure, a person's sex or the ploidy of an individual genotype.
    for prefix in ("", "chr"):
        if all(
            prefix + c in header.contigs and header.contigs[prefix + c].length == length
            for c, (length, _) in GRCH38.items()
        ):
            return [
                {
                    "chrom": "chr" + c,
                    "start0": lo,
                    "end0": hi,
                    "profile": "legacy_family_GRCh38_compatible",
                    "scope": "negative",
                    "reason": "GRCh38_compatible_PAR_coordinates_mapping_profile_unverified",
                }
                for c, (_, intervals) in GRCH38.items()
                for lo, hi in intervals
            ]
    return []


def source_par_scope(
    header: pysam.VariantHeader, *, reference_sha256: str | None = None
) -> list[dict[str, Any]]:
    rows = [r for r in header.records if r.key == TAG]
    if not rows:
        if any(
            r.key
            in {
                "family_reference_confidence_overlay",
                "family_confidence_policy",
                "family_ploidy_normalization",
            }
            for r in header.records
        ):
            return _compatible_family_intervals(header)
        return []
    out = []
    seen = set()
    profile = rows[0].get("Profile")
    pins = {k: rows[0].get(k) for k in PINS} if profile == MASKED else {}
    if profile not in {MASKED, UNMASKED, DECLARED}:
        raise InputValidationError("unrecognized family PAR query profile")
    if pins and (
        not all(re.fullmatch("[0-9a-f]{64}", str(v)) for v in pins.values())
        or pins[PINS[0]] == pins[PINS[1]]
    ):
        raise InputValidationError("invalid family PAR mapping SHA-256 declarations")
    if pins and reference_sha256 and pins["CanonicalSHA256"] != reference_sha256:
        raise InputValidationError("family PAR canonical reference differs from the pinned FASTA")
    for r in rows:
        contig = str(r.get("Contig"))
        chrom = contig.removeprefix("chr")
        if (
            chrom not in GRCH38
            or contig not in header.contigs
            or r.get("Profile") != profile
            or any(r.get(k) != v for k, v in pins.items())
        ):
            raise InputValidationError("inconsistent family PAR query declaration")
        try:
            lo, hi = int(r["Start0"]), int(r["End0"])
        except (KeyError, ValueError, TypeError) as exc:
            raise InputValidationError("invalid family PAR query interval") from exc
        length = header.contigs[contig].length
        if (
            not 0 <= lo < hi
            or (length is not None and hi > length)
            or any(c == chrom and a < hi and b > lo for c, a, b in seen)
        ):
            raise InputValidationError("invalid or overlapping family PAR query interval")
        seen.add((chrom, lo, hi))
        mask_y = profile == MASKED and chrom == "Y"
        if profile == MASKED:
            reason = (
                "Y_PAR_unassessed_homologs_represented_together_on_X"
                if mask_y
                else "diploid_X_PAR_combined_homologs_negative_sensitivity_unqualified"
            )
        elif profile == UNMASKED:
            reason = "unmasked_GRCh38_homologous_PAR_negative_unassessed"
        else:
            reason = "declared_PAR_mapping_and_negative_sensitivity_unassessed"
        out.append(
            {
                "chrom": "chr" + chrom,
                "start0": lo,
                "end0": hi,
                "profile": profile,
                "scope": "all" if mask_y else "negative",
                "reason": reason,
                **pins,
            }
        )
    if profile in {MASKED, UNMASKED} and seen != {
        (c, lo, hi) for c, (_, intervals) in GRCH38.items() for lo, hi in intervals
    }:
        raise InputValidationError("incomplete GRCh38 family PAR declarations")
    if profile == DECLARED:
        # The producer's declared intervals replace the fallback on that contig only.
        declared = {r["chrom"] for r in out}
        out.extend(r for r in _compatible_family_intervals(header) if r["chrom"] not in declared)
    return out
