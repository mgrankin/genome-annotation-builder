"""Shared, caller-neutral HLA parsing and identity helpers."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from pathlib import Path

from .portable_model import PortableEvidenceError, sha256_file

PRIMARY_HLA_LOCI = ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
UNRESOLVED_ALLELES = frozenset({"", ".", "-", "NA", "N/A", "None"})
HLA_ALLELE = re.compile(r"^(?P<gene>[A-Z0-9]+)\*(?P<fields>[0-9]+(?::[0-9]+)*)(?P<suffix>[A-Z]?)$")


def parse_hla_allele(value: str, *, expected_gene: str) -> str | None:
    """Return an exact allele or ``None`` for an explicit unresolved token."""
    allele = value.strip()
    if allele in UNRESOLVED_ALLELES:
        return None
    match = HLA_ALLELE.fullmatch(allele)
    if match is None or match.group("gene") != expected_gene:
        raise PortableEvidenceError("invalid or locus-mismatched HLA allele")
    return allele


def two_field_allele(allele: str | None) -> str | None:
    """Project an exact HLA allele to the common two-field comparison boundary."""
    if allele is None:
        return None
    match = HLA_ALLELE.fullmatch(allele)
    if match is None:
        raise PortableEvidenceError("cannot project an invalid HLA allele")
    fields = match.group("fields").split(":")
    if len(fields) < 2:
        return None
    return f"{match.group('gene')}*{fields[0]}:{fields[1]}"


def require_g_group(allele: str | None) -> str | None:
    if allele is not None and not allele.endswith("G"):
        raise PortableEvidenceError("HLA G-group result lacks the G suffix")
    return allele


def finite_float(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise PortableEvidenceError(f"invalid {label}") from exc
    if parsed != parsed or parsed in {float("inf"), float("-inf")}:
        raise PortableEvidenceError(f"non-finite {label}")
    return parsed


def artifact_set_identity(paths: Mapping[str, str | Path]) -> tuple[str, dict[str, str]]:
    """Digest a named artifact set without persisting local paths or basenames."""
    components = {role: sha256_file(Path(path)) for role, path in sorted(paths.items())}
    payload = json.dumps(components, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest(), components


def semantic_record_digest(payload: Mapping[str, object]) -> str:
    """Digest a privacy-minimized semantic source projection."""
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()
