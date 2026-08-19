"""Loss-minimized DELLY cohort-site planning and reconciliation.

The module does not execute a shell and never serializes sample identifiers into
portable summaries.  It builds argv vectors for an external restart-safe runner
and validates the resulting BCF/VCF artifacts.  Discovery calls remain separate
evidence; a force-genotyped ``0/0`` only means reference at a supplied candidate
site under DELLY's diploid model.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from .portable_model import PortableEvidenceError, sha256_file

BND_PARTNER = re.compile(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]")
SV_TYPES = frozenset({"DEL", "DUP", "INV", "INS", "BND"})
GENOTYPE_FIELDS = ("GT", "GQ", "PL", "DR", "DV", "RR", "RV", "RC", "RDCN", "FT")


@dataclass(frozen=True)
class DellyCohortPolicy:
    """Explicit policy for a candidate-complete, non-refining DELLY workflow."""

    caller_version: str = "2.5.1"
    merge_min_quality: int = 0
    merge_recurrent: int = 1
    merge_min_vaf: float = 0.0
    merge_min_coverage: int = 0
    merge_min_size: int = 0
    merge_max_size: int = 2_000_000_000
    merge_representative_min_af: float = 0.0
    genotype_min_mapq: int = 5
    genotype_max_reads: int = 250
    threads: int = 4
    filter_min_quality: int = 300
    filter_min_size: int = 0
    filter_max_size: int = 500_000_000

    def __post_init__(self) -> None:
        if self.caller_version != "2.5.1":
            raise PortableEvidenceError("DELLY cohort policy requires the audited 2.5.1 boundary")
        if self.merge_min_quality != 0:
            raise PortableEvidenceError(
                "candidate-complete DELLY merge requires quality threshold zero"
            )
        if self.merge_recurrent != 1:
            raise PortableEvidenceError(
                "candidate-complete DELLY merge requires recurrent threshold one"
            )
        if self.merge_min_vaf != 0 or self.merge_min_coverage != 0:
            raise PortableEvidenceError(
                "candidate-complete DELLY merge requires zero VAF/coverage gates"
            )
        if self.merge_min_size != 0 or self.merge_representative_min_af != 0:
            raise PortableEvidenceError(
                "candidate-complete DELLY merge cannot drop small or rare alleles"
            )
        if self.merge_max_size < self.filter_max_size or self.merge_max_size <= 0:
            raise PortableEvidenceError("invalid DELLY maximum-size policy")
        if self.genotype_min_mapq < 0 or self.genotype_max_reads < 1 or self.threads < 1:
            raise PortableEvidenceError("invalid DELLY genotyping policy")
        if self.filter_min_quality < 0 or self.filter_min_size < 0:
            raise PortableEvidenceError("invalid DELLY tag-only filter policy")


@dataclass(frozen=True)
class DellyTechnicalInput:
    """Private runtime paths for one technical source; never a portable record."""

    discovery_vcf: Path
    sample_free_vcf: Path
    alignment: Path
    sample_name: str
    genotyped_vcf: Path
    tagged_single_sample_vcf: Path

    def __post_init__(self) -> None:
        if not self.sample_name or any(char in self.sample_name for char in "\t\r\n"):
            raise PortableEvidenceError("invalid DELLY runtime sample name")


def _argv_path(value: str | Path) -> str:
    rendered = str(value)
    if not rendered or any(char in rendered for char in "\r\n\0"):
        raise PortableEvidenceError("invalid path in DELLY command plan")
    return rendered


def build_delly_cohort_command_plan(
    *,
    delly: str | Path,
    bcftools: str | Path,
    reference_fasta: str | Path,
    exclude_regions: str | Path,
    technical_inputs: Sequence[DellyTechnicalInput],
    candidate_sites_vcf: str | Path,
    raw_multisample_vcf: str | Path,
    tagged_multisample_vcf: str | Path,
    policy: DellyCohortPolicy | None = None,
) -> list[dict[str, Any]]:
    """Return shell-free argv vectors for discovery projection through tagged split."""

    if policy is None:
        policy = DellyCohortPolicy()
    if not technical_inputs:
        raise PortableEvidenceError("DELLY cohort workflow requires at least one technical input")
    sample_names = [item.sample_name for item in technical_inputs]
    if len(sample_names) != len(set(sample_names)):
        raise PortableEvidenceError("DELLY runtime sample names must be unique")
    delly_bin = _argv_path(delly)
    bcftools_bin = _argv_path(bcftools)
    reference = _argv_path(reference_fasta)
    exclusion = _argv_path(exclude_regions)
    sites = _argv_path(candidate_sites_vcf)
    raw_multi = _argv_path(raw_multisample_vcf)
    tagged_multi = _argv_path(tagged_multisample_vcf)
    plan: list[dict[str, Any]] = []
    for ordinal, item in enumerate(technical_inputs, start=1):
        projected = _argv_path(item.sample_free_vcf)
        plan.append(
            {
                "stage": "project_discovery_to_sample_free_site_input",
                "source_ordinal": ordinal,
                "argv": [
                    bcftools_bin,
                    "view",
                    "--no-version",
                    "-G",
                    "-Ob",
                    "-o",
                    projected,
                    _argv_path(item.discovery_vcf),
                ],
            }
        )
        plan.append(
            {
                "stage": "index_sample_free_site_input",
                "source_ordinal": ordinal,
                "argv": [bcftools_bin, "index", "-f", projected],
            }
        )
    plan.append(
        {
            "stage": "merge_candidate_sites",
            "argv": [
                delly_bin,
                "merge",
                "-g",
                reference,
                "-o",
                sites,
                "-y",
                str(policy.merge_min_quality),
                "-q",
                str(policy.merge_recurrent),
                "-a",
                str(policy.merge_min_vaf),
                "-v",
                str(policy.merge_min_coverage),
                "-m",
                str(policy.merge_min_size),
                "-n",
                str(policy.merge_max_size),
                "--rep-min-af",
                str(policy.merge_representative_min_af),
                *[_argv_path(item.sample_free_vcf) for item in technical_inputs],
            ],
        }
    )
    plan.append({"stage": "index_candidate_sites", "argv": [bcftools_bin, "index", "-f", sites]})
    for ordinal, item in enumerate(technical_inputs, start=1):
        genotyped = _argv_path(item.genotyped_vcf)
        plan.append(
            {
                "stage": "force_genotype_candidate_sites",
                "source_ordinal": ordinal,
                "argv": [
                    delly_bin,
                    "sr",
                    "-g",
                    reference,
                    "-x",
                    exclusion,
                    "-v",
                    sites,
                    "-u",
                    str(policy.genotype_min_mapq),
                    "-a",
                    str(policy.genotype_max_reads),
                    "-h",
                    str(policy.threads),
                    "-o",
                    genotyped,
                    _argv_path(item.alignment),
                ],
            }
        )
        plan.append(
            {
                "stage": "index_force_genotype_output",
                "source_ordinal": ordinal,
                "argv": [bcftools_bin, "index", "-f", genotyped],
            }
        )
    plan.append(
        {
            "stage": "merge_force_genotypes_without_site_loss",
            "argv": [
                bcftools_bin,
                "merge",
                "--no-version",
                "-m",
                "id",
                "-Ob",
                "-o",
                raw_multi,
                *[_argv_path(item.genotyped_vcf) for item in technical_inputs],
            ],
        }
    )
    plan.append(
        {"stage": "index_raw_multisample", "argv": [bcftools_bin, "index", "-f", raw_multi]}
    )
    plan.append(
        {
            "stage": "tag_germline_filters_without_refinement_or_collapse",
            "argv": [
                delly_bin,
                "filter",
                "-f",
                "germline",
                "-t",
                "--no-collapse",
                "--no-refine",
                "-y",
                str(policy.filter_min_quality),
                "-m",
                str(policy.filter_min_size),
                "-n",
                str(policy.filter_max_size),
                "-o",
                tagged_multi,
                raw_multi,
            ],
        }
    )
    plan.append(
        {
            "stage": "index_tagged_multisample",
            "argv": [bcftools_bin, "index", "-f", tagged_multi],
        }
    )
    for ordinal, item in enumerate(technical_inputs, start=1):
        split = _argv_path(item.tagged_single_sample_vcf)
        plan.append(
            {
                "stage": "split_tagged_single_sample",
                "source_ordinal": ordinal,
                "argv": [
                    bcftools_bin,
                    "view",
                    "--no-version",
                    "-s",
                    item.sample_name,
                    "-Ob",
                    "-o",
                    split,
                    tagged_multi,
                ],
            }
        )
        plan.append(
            {
                "stage": "index_tagged_single_sample",
                "source_ordinal": ordinal,
                "argv": [bcftools_bin, "index", "-f", split],
            }
        )
    return plan


def _canonical_chrom(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise PortableEvidenceError("invalid DELLY chromosome")
    core = value[3:] if value.startswith("chr") else value
    if core in {str(index) for index in range(1, 23)} | {"X", "Y"}:
        return f"chr{core}"
    raise PortableEvidenceError("DELLY cohort site uses a non-canonical nuclear chromosome")


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in sorted(value.items())}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    if isinstance(value, bytes):
        return value.decode()
    return value


def _digest(value: Any) -> str:
    rendered = json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode()
    return hashlib.sha256(rendered).hexdigest()


def _info(record: pysam.VariantRecord, key: str) -> Any:
    return record.info.get(key) if key in record.header.info else None


def _ci(record: pysam.VariantRecord, key: str) -> tuple[int, int]:
    value = _info(record, key)
    if value is None:
        if bool(_info(record, "PRECISE")):
            return (0, 0)
        raise PortableEvidenceError(f"imprecise DELLY cohort site lacks INFO/{key}")
    values = tuple(value) if isinstance(value, tuple) else (value,)
    if len(values) != 2 or any(item is None for item in values):
        raise PortableEvidenceError(f"invalid DELLY INFO/{key}")
    result = (int(values[0]), int(values[1]))
    if result[0] > result[1]:
        raise PortableEvidenceError(f"reversed DELLY INFO/{key}")
    return result


def _site_definition(record: pysam.VariantRecord, *, require_merge_fields: bool) -> dict[str, Any]:
    record_id = record.id
    if not record_id or record_id == ".":
        raise PortableEvidenceError("DELLY cohort site lacks a caller ID")
    sv_type = _info(record, "SVTYPE")
    if sv_type not in SV_TYPES:
        raise PortableEvidenceError("unsupported DELLY cohort SVTYPE")
    orientation = _info(record, "CT")
    if not isinstance(orientation, str) or not orientation:
        raise PortableEvidenceError("DELLY cohort site lacks INFO/CT")
    precise = bool(_info(record, "PRECISE"))
    imprecise = bool(_info(record, "IMPRECISE"))
    if precise == imprecise:
        raise PortableEvidenceError("DELLY cohort precision flags are ambiguous")
    chrom1 = _canonical_chrom(record.contig)
    pos1 = int(record.pos)
    alts = tuple(str(value) for value in (record.alts or ()))
    if pos1 < 1 or len(alts) != 1:
        raise PortableEvidenceError("invalid DELLY cohort site alleles")
    if sv_type == "BND":
        chrom2_raw = _info(record, "CHR2")
        pos2_raw = _info(record, "POS2")
        partner = BND_PARTNER.search(alts[0])
        if chrom2_raw is None or pos2_raw is None or partner is None:
            raise PortableEvidenceError("DELLY BND site lacks reciprocal partner coordinates")
        chrom2 = _canonical_chrom(str(chrom2_raw))
        pos2 = int(pos2_raw)
        if (_canonical_chrom(partner.group(1)), int(partner.group(2))) != (chrom2, pos2):
            raise PortableEvidenceError("DELLY BND ALT and INFO partner coordinates disagree")
    else:
        chrom2 = chrom1
        pos2 = int(record.stop)
    if pos2 < 1 or (sv_type != "BND" and pos2 < pos1):
        raise PortableEvidenceError("invalid DELLY cohort secondary breakpoint")
    allele_id = _info(record, "ALLELEID")
    allele_count = _info(record, "NALLELE")
    support = _info(record, "SUPP")
    if require_merge_fields:
        if allele_id is None or int(allele_id) < 0:
            raise PortableEvidenceError("merged DELLY site lacks INFO/ALLELEID")
        if allele_count is None or int(allele_count) < 1:
            raise PortableEvidenceError("merged DELLY site lacks INFO/NALLELE")
        if support is None or int(support) < 1:
            raise PortableEvidenceError("merged DELLY site lacks INFO/SUPP")
    definition = {
        "record_id": str(record_id),
        "chrom1": chrom1,
        "pos1": pos1,
        "chrom2": chrom2,
        "pos2": pos2,
        "sv_type": str(sv_type),
        "orientation": orientation,
        "ref": str(record.ref),
        "alt": alts[0],
        "sv_length": int(_info(record, "SVLEN")) if _info(record, "SVLEN") is not None else None,
        "insertion_length": (
            int(_info(record, "INSLEN")) if _info(record, "INSLEN") is not None else None
        ),
        "precise": precise,
        "imprecise": imprecise,
        "cipos": _ci(record, "CIPOS"),
        "ciend": _ci(record, "CIEND"),
        "allele_id": int(allele_id) if allele_id is not None else None,
        "allele_count": int(allele_count) if allele_count is not None else None,
        "support": int(support) if support is not None else None,
    }
    definition["structural_digest"] = _digest(
        {
            key: value
            for key, value in definition.items()
            if key
            not in {
                "record_id",
                "allele_id",
                "allele_count",
                "support",
                "sv_length",
            }
        }
    )
    definition["definition_digest"] = _digest(
        {key: value for key, value in definition.items() if key != "structural_digest"}
    )
    return definition


def _read_definitions(
    path: str | Path,
    *,
    expected_sample_count: int | None,
    require_merge_fields: bool,
) -> tuple[list[dict[str, Any]], tuple[str, ...]]:
    try:
        source = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise PortableEvidenceError("invalid DELLY cohort BCF/VCF artifact") from exc
    definitions: list[dict[str, Any]] = []
    try:
        samples = tuple(source.header.samples)
        if expected_sample_count is not None and len(samples) != expected_sample_count:
            raise PortableEvidenceError("DELLY cohort artifact has the wrong sample cardinality")
        seen: set[str] = set()
        for record in source:
            definition = _site_definition(record, require_merge_fields=require_merge_fields)
            record_id = definition["record_id"]
            if record_id in seen:
                raise PortableEvidenceError("duplicate DELLY cohort caller record ID")
            seen.add(record_id)
            definitions.append(definition)
    finally:
        source.close()
    if not definitions:
        raise PortableEvidenceError("empty DELLY cohort artifact")
    return definitions, samples


def validate_sample_free_projection(
    discovery_vcf: str | Path, sample_free_vcf: str | Path
) -> dict[str, Any]:
    """Prove that sample stripping changed no site-level VCF columns."""

    def site_lines(path: str | Path, expected_samples: int) -> tuple[list[str], tuple[str, ...]]:
        try:
            source = pysam.VariantFile(str(path))
        except (OSError, ValueError) as exc:
            raise PortableEvidenceError("invalid DELLY discovery projection artifact") from exc
        try:
            samples = tuple(source.header.samples)
            if len(samples) != expected_samples:
                raise PortableEvidenceError(
                    "DELLY discovery projection has wrong sample cardinality"
                )
            lines = ["\t".join(str(record).rstrip("\n").split("\t")[:8]) for record in source]
        finally:
            source.close()
        if not lines:
            raise PortableEvidenceError("empty DELLY discovery projection")
        return lines, samples

    original, private_samples = site_lines(discovery_vcf, 1)
    projected, projected_samples = site_lines(sample_free_vcf, 0)
    if original != projected:
        raise PortableEvidenceError("sample-free DELLY projection changed site-level records")
    assert private_samples and not projected_samples
    return {
        "status": "PASS",
        "records": len(original),
        "site_stream_sha256": _digest(original),
        "discovery_artifact_sha256": sha256_file(Path(discovery_vcf)),
        "sample_free_artifact_sha256": sha256_file(Path(sample_free_vcf)),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }


def _candidate_match(source: dict[str, Any], site: dict[str, Any]) -> bool:
    if any(source[key] != site[key] for key in ("chrom1", "chrom2", "sv_type", "orientation")):
        return False
    pos1_low = site["pos1"] + site["cipos"][0]
    pos1_high = site["pos1"] + site["cipos"][1]
    pos2_low = site["pos2"] + site["ciend"][0]
    pos2_high = site["pos2"] + site["ciend"][1]
    return pos1_low <= source["pos1"] <= pos1_high and pos2_low <= source["pos2"] <= pos2_high


def validate_delly_candidate_union(
    sample_free_site_inputs: Sequence[str | Path], candidate_sites_vcf: str | Path
) -> dict[str, Any]:
    """Reconcile every sample-free discovery record to the native merged site list."""

    if not sample_free_site_inputs:
        raise PortableEvidenceError("DELLY candidate union has no discovery inputs")
    source_rows: list[dict[str, Any]] = []
    for source_ordinal, path in enumerate(sample_free_site_inputs, start=1):
        definitions, _ = _read_definitions(
            path, expected_sample_count=0, require_merge_fields=False
        )
        for definition in definitions:
            source_rows.append({"source_ordinal": source_ordinal, **definition})
    sites, _ = _read_definitions(
        candidate_sites_vcf, expected_sample_count=0, require_merge_fields=True
    )
    source_memberships: list[dict[str, Any]] = []
    covered_sites: set[str] = set()
    for source in source_rows:
        matches = sorted(
            site["definition_digest"] for site in sites if _candidate_match(source, site)
        )
        if not matches:
            raise PortableEvidenceError(
                "a DELLY discovery record is absent from the candidate site list"
            )
        covered_sites.update(matches)
        source_memberships.append(
            {
                "source_ordinal": source["source_ordinal"],
                "source_definition_digest": source["definition_digest"],
                "candidate_site_definition_digests": matches,
            }
        )
    missing_sites = [site for site in sites if site["definition_digest"] not in covered_sites]
    if missing_sites:
        raise PortableEvidenceError("a merged DELLY candidate site has no discovery provenance")
    ambiguous = sum(
        len(item["candidate_site_definition_digests"]) > 1 for item in source_memberships
    )
    return {
        "status": "PASS",
        "technical_sources": len(sample_free_site_inputs),
        "projected_discovery_records": len(source_rows),
        "candidate_sites": len(sites),
        "matched_discovery_records": len(source_memberships),
        "unmatched_discovery_records": 0,
        "ambiguous_discovery_memberships": ambiguous,
        "candidate_sites_without_discovery_provenance": 0,
        "candidate_site_artifact_sha256": sha256_file(Path(candidate_sites_vcf)),
        "candidate_site_definition_multiset_sha256": _digest(
            sorted(site["definition_digest"] for site in sites)
        ),
        "membership_multiset_sha256": _digest(
            sorted(source_memberships, key=lambda item: _digest(item))
        ),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }


def validate_force_genotype_outputs(
    candidate_sites_vcf: str | Path, genotyped_vcfs: Sequence[str | Path]
) -> dict[str, Any]:
    """Require every technical source to contain exactly the candidate site universe."""

    if not genotyped_vcfs:
        raise PortableEvidenceError("DELLY force-genotype validation has no outputs")
    sites, _ = _read_definitions(
        candidate_sites_vcf, expected_sample_count=0, require_merge_fields=True
    )
    expected = [site["structural_digest"] for site in sites]
    seen_samples: set[str] = set()
    artifact_digests: list[str] = []
    derived_length_disagreements = 0
    caller_id_mappings: list[list[tuple[str, str]]] = []
    for path in genotyped_vcfs:
        definitions, samples = _read_definitions(
            path, expected_sample_count=1, require_merge_fields=False
        )
        if samples[0] in seen_samples:
            raise PortableEvidenceError("duplicate sample in DELLY force-genotype outputs")
        seen_samples.add(samples[0])
        observed = [definition["structural_digest"] for definition in definitions]
        if observed != expected:
            raise PortableEvidenceError("DELLY force-genotype output differs from candidate sites")
        caller_id_mappings.append(
            [
                (site["record_id"], definition["record_id"])
                for site, definition in zip(sites, definitions, strict=True)
            ]
        )
        derived_length_disagreements += sum(
            site["sv_length"] != definition["sv_length"]
            for site, definition in zip(sites, definitions, strict=True)
        )
        artifact_digests.append(sha256_file(Path(path)))
    if any(mapping != caller_id_mappings[0] for mapping in caller_id_mappings[1:]):
        raise PortableEvidenceError("DELLY force-genotype outputs disagree on caller ID mapping")
    caller_id_reassignments = sum(first != second for first, second in caller_id_mappings[0])
    return {
        "status": "PASS",
        "technical_sources": len(genotyped_vcfs),
        "candidate_sites": len(sites),
        "expected_genotype_rows": len(sites) * len(genotyped_vcfs),
        "validated_genotype_rows": len(sites) * len(genotyped_vcfs),
        "missing_or_extra_rows": 0,
        "derived_sv_length_disagreements": derived_length_disagreements,
        "caller_id_reassignments": caller_id_reassignments,
        "caller_id_mapping_sha256": _digest(caller_id_mappings[0]),
        "candidate_site_definition_multiset_sha256": _digest(expected),
        "force_genotype_artifact_multiset_sha256": _digest(sorted(artifact_digests)),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }


def _sample_payload(record: pysam.VariantRecord, sample: str) -> dict[str, Any]:
    values = record.samples[sample]
    return {
        key: _jsonable(values.get(key)) if key in record.header.formats else None
        for key in GENOTYPE_FIELDS
    }


def _read_multisample_payload(path: str | Path) -> tuple[list[dict[str, Any]], tuple[str, ...]]:
    try:
        source = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise PortableEvidenceError("invalid DELLY multisample artifact") from exc
    rows: list[dict[str, Any]] = []
    try:
        samples = tuple(source.header.samples)
        if not samples:
            raise PortableEvidenceError("DELLY multisample artifact has no samples")
        for record in source:
            definition = _site_definition(record, require_merge_fields=False)
            rows.append(
                {
                    "definition_digest": definition["definition_digest"],
                    "filters": sorted(str(value) for value in record.filter),
                    "samples": [_sample_payload(record, sample) for sample in samples],
                }
            )
    finally:
        source.close()
    if not rows:
        raise PortableEvidenceError("empty DELLY multisample artifact")
    return rows, samples


def validate_tag_only_germline_filter(
    raw_multisample_vcf: str | Path, tagged_multisample_vcf: str | Path
) -> dict[str, Any]:
    """Prove that filtering retained every site and changed no GT/evidence payload."""

    raw_rows, raw_samples = _read_multisample_payload(raw_multisample_vcf)
    tagged_rows, tagged_samples = _read_multisample_payload(tagged_multisample_vcf)
    if raw_samples != tagged_samples:
        raise PortableEvidenceError("DELLY tag-only filter changed sample order or identity")
    raw_biology = [
        {"definition_digest": row["definition_digest"], "samples": row["samples"]}
        for row in raw_rows
    ]
    tagged_biology = [
        {"definition_digest": row["definition_digest"], "samples": row["samples"]}
        for row in tagged_rows
    ]
    if raw_biology != tagged_biology:
        raise PortableEvidenceError("DELLY tag-only filter changed sites, GT, or read evidence")
    filter_counts: dict[str, int] = {}
    for row in tagged_rows:
        label = ";".join(row["filters"]) if row["filters"] else "."
        filter_counts[label] = filter_counts.get(label, 0) + 1
    biological_digest = _digest(raw_biology)
    return {
        "status": "PASS",
        "technical_sources": len(raw_samples),
        "raw_sites": len(raw_rows),
        "tagged_sites": len(tagged_rows),
        "removed_sites": 0,
        "changed_site_definitions": 0,
        "changed_genotypes_or_evidence": 0,
        "biological_payload_sha256_before": biological_digest,
        "biological_payload_sha256_after": biological_digest,
        "tagged_filter_counts": dict(sorted(filter_counts.items())),
        "raw_multisample_artifact_sha256": sha256_file(Path(raw_multisample_vcf)),
        "tagged_multisample_artifact_sha256": sha256_file(Path(tagged_multisample_vcf)),
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }


def command_plan_semantic_summary(plan: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Summarize stages without serializing private argv tokens."""

    stages: dict[str, int] = {}
    commands = list(plan)
    for command in commands:
        stage = command.get("stage")
        argv = command.get("argv")
        if not isinstance(stage, str) or not stage or not isinstance(argv, list) or not argv:
            raise PortableEvidenceError("invalid DELLY command plan")
        stages[stage] = stages.get(stage, 0) + 1
    return {
        "status": "PASS",
        "commands": len(commands),
        "stage_counts": dict(sorted(stages.items())),
        "plan_shape_sha256": _digest([(item["stage"], len(item["argv"])) for item in commands]),
        "argv_persisted": False,
        "sample_identifier_persisted": False,
        "private_paths_persisted": False,
    }
