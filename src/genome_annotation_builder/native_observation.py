"""Original caller ALT evidence retained when a family projection sets haploid GT to missing."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any

import pysam

from .errors import InputValidationError

REASON = "native_carried_alt_with_missing_haploid_GT"
PROJECTION_FLAG = "GAB_NATIVE_NO_CALL"
SCHEMA = "family-native-projection/2.2"


def normalized_native_header(header: pysam.VariantHeader) -> bool:
    schemas = [r.value for r in header.records if r.key == "native_projection_schema"]
    if SCHEMA not in schemas:
        return False
    if schemas != [SCHEMA]:
        raise InputValidationError("ambiguous native projection schema")
    policies = [r.get("ID") for r in header.records if r.key == "family_ploidy_normalization"]
    if not policies:
        return False
    if policies != ["Confirmed_XX_XY_MT"]:
        raise InputValidationError("unsupported native ploidy normalization declaration")
    for fields, declarations in (
        (header.info, (("NID", 1, "String"), ("NAI", "A", "Integer"))),
        (
            header.formats,
            (
                ("NGT", 1, "String"),
                ("BGT", 1, "String"),
                *((tag, 1, "String") for tag in ("NAD", "NPL", "NDP", "NGQ", "NFT")),
                ("NSTATE", 1, "String"),
                ("NDOS", 1, "Integer"),
                ("NOTH", 1, "Integer"),
                ("NUNK", 1, "Integer"),
            ),
        ),
    ):
        for tag, number, kind in declarations:
            if tag not in fields or (fields[tag].number, fields[tag].type) != (number, kind):
                raise InputValidationError(f"invalid native projection header for {tag}")
    return True


def validate_native_no_call(
    info: Mapping[str, Any],
    sample: Mapping[str, Any],
    gt: Sequence[int | None],
    alts: Sequence[str],
    chrom: str,
    pos: int,
    end: int,
    par_scope: list[dict[str, Any]],
    *,
    sex: int = 0,
) -> bool:
    """Check the serialized native identity, source slots, and expected haploid collapse.

    Native ALT slots describe the original caller observation.
    The selected ALT index is one on the delivered biallelic row, while NAI
    refers to the unsplit native GT and can be larger.
    """
    indices = info.get("NAI")
    nid = info.get("NID")
    if (
        not isinstance(nid, str)
        or not re.fullmatch(r"[a-f0-9]{64}", nid)
        or len(alts) != 1
        or not isinstance(indices, (list, tuple))
        or len(indices) != 1
        or type(indices[0]) is not int
        or indices[0] < 1
    ):
        raise InputValidationError(
            "invalid native allele identity on a genotype without target ALT"
        )
    slots = sample.get("NGT")
    if not isinstance(slots, str) or not re.fullmatch(
        r"(?:[0-9]+|\.)(?:[|/](?:[0-9]+|\.))?", slots
    ):
        raise InputValidationError("invalid native genotype on projected no-call")
    index = indices[0]
    alleles = [None if a == "." else int(a) for a in re.split(r"[/|]", slots)]
    dosage = alleles.count(index)
    other = sum(a not in (None, 0, index) for a in alleles)
    unknown = alleles.count(None)
    observed = [sample.get(tag) for tag in ("NDOS", "NOTH", "NUNK")]
    if not all(type(n) is int for n in observed) or observed != [dosage, other, unknown]:
        raise InputValidationError("native allele counts disagree with NGT and NAI")
    if not dosage:
        if other:
            state = "OTHER_ALT_AND_UNKNOWN" if unknown else "OTHER_KNOWN_ALT"
        else:
            state = "UNKNOWN" if unknown else "REFERENCE_GT_UNQUALIFIED"
        if sample.get("NSTATE") != state:
            raise InputValidationError("native allele state disagrees with NGT and NAI")
        return False
    projected = re.sub(r"[0-9]+", lambda m: "1" if int(m[0]) == index else "0", slots)
    if (
        tuple(gt) != (None,)
        or (dosage, other, unknown) != (1, 0, 1)
        or sample.get("NSTATE") != "CARRIED_PARTIAL_NATIVE"
        or sample.get("BGT") != projected
    ):
        raise InputValidationError("native target ALT disagrees with projected missing haploid GT")
    if chrom not in {"chrM", "chrX", "chrY"} or (
        chrom in {"chrX", "chrY"}
        and (
            sex == 2
            or any(
                row["chrom"] == chrom and pos <= row["end0"] and row["start0"] < end
                for row in par_scope
            )
        )
    ):
        raise InputValidationError("native carried observation contradicts the haploid scope")
    return True
