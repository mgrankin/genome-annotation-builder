"""Run fixed public semantic controls through the real offline VEP plugins."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pysam

from .config import BuildConfig
from .constants import VEP_VERSION
from .errors import InputValidationError, ResourceError
from .util import (
    atomic_output_directory,
    atomic_write_bytes,
    atomic_write_json,
    canonical_json_bytes,
    sha256_file,
    sha256_json,
)
from .vep_parser import VepParser
from .vep_runner import run_vep_annotation

DOCTOR_SCHEMA_VERSION = 1
CONTROL_SET = "gab-public-predictor-controls-grch38-v1"
SEMANTIC_BOUNDARY = (
    "Fixed public positive controls verify resource lookup and typed score extraction; "
    "they do not establish genome-wide sensitivity, specificity, or clinical validity."
)

CONTROL_VARIANTS: tuple[dict[str, Any], ...] = (
    {
        "control_id": "or4f5_missense",
        "chrom": "1",
        "pos": 69094,
        "ref": "G",
        "alt": "T",
        "gene": "OR4F5",
        "consequence": "missense_variant",
        "checks": {"cadd": 21.0, "spliceai": 0.03},
    },
    {
        "control_id": "or4f5_stop_gained",
        "chrom": "1",
        "pos": 69130,
        "ref": "C",
        "alt": "T",
        "gene": "OR4F5",
        "consequence": "stop_gained",
        "checks": {"loftee": "HC", "cadd": 25.4},
    },
    {
        "control_id": "tp53_missense",
        "chrom": "17",
        "pos": 7669613,
        "ref": "T",
        "alt": "G",
        "gene": "TP53",
        "consequence": "missense_variant",
        "checks": {"cadd": 23.4, "revel": 0.675, "alphamissense": 0.5749},
    },
)
CONTROL_DEFINITION_SHA256 = hashlib.sha256(canonical_json_bytes(CONTROL_VARIANTS)).hexdigest()


def _contig(config: BuildConfig, value: str) -> str:
    return "chr" + value if config.contig_style == "chr" else value


def _control_vcf_bytes(config: BuildConfig) -> bytes:
    contigs: dict[str, int] = {}
    with pysam.FastaFile(str(config.reference_fasta)) as reference:
        for control in CONTROL_VARIANTS:
            contig = _contig(config, control["chrom"])
            if contig not in reference.references:
                raise ResourceError(f"predictor control contig is absent from reference: {contig}")
            observed = reference.fetch(contig, control["pos"] - 1, control["pos"]).upper()
            if observed != control["ref"]:
                raise ResourceError(
                    f"predictor control REF mismatch at {contig}:{control['pos']}: "
                    f"{control['ref']} != {observed}"
                )
            contigs[contig] = reference.get_reference_length(contig)
    lines = ["##fileformat=VCFv4.2", "##reference=GRCh38"]
    lines.extend(f"##contig=<ID={name},length={length}>" for name, length in contigs.items())
    lines.extend(
        [
            '##INFO=<ID=GAB_SOURCE_ID,Number=1,Type=String,Description="Public control ID">',
            '##INFO=<ID=GAB_SOURCE_DIGEST,Number=1,Type=String,Description="Public control digest">',
            '##INFO=<ID=GAB_SOURCE_ORDINAL,Number=1,Type=Integer,Description="Public control ordinal">',
            '##INFO=<ID=GAB_ALT_INDEX,Number=1,Type=Integer,Description="Public control ALT index">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Synthetic depth">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Synthetic genotype quality">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Synthetic allelic depths">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
        ]
    )
    for ordinal, control in enumerate(CONTROL_VARIANTS, start=1):
        digest = hashlib.sha256(control["control_id"].encode()).hexdigest()
        info = (
            f"GAB_SOURCE_ID={digest};GAB_SOURCE_DIGEST={digest};"
            f"GAB_SOURCE_ORDINAL={ordinal};GAB_ALT_INDEX=1"
        )
        lines.append(
            "\t".join(
                [
                    _contig(config, control["chrom"]),
                    str(control["pos"]),
                    control["control_id"],
                    control["ref"],
                    control["alt"],
                    "60",
                    "PASS",
                    info,
                    "GT:DP:GQ:AD",
                    "0/1:30:60:15,15",
                ]
            )
        )
    return ("\n".join(lines) + "\n").encode()


def _write_control_vcf(config: BuildConfig, destination: Path) -> None:
    plain = destination.with_suffix("")
    atomic_write_bytes(plain, _control_vcf_bytes(config))
    pysam.tabix_compress(str(plain), str(destination), force=True)
    pysam.tabix_index(str(destination), preset="vcf", force=True)


def _identity(path: Path) -> tuple[list[tuple[Any, ...]], str]:
    rows: list[tuple[Any, ...]] = []
    digest = hashlib.sha256()
    with pysam.VariantFile(str(path)) as handle:
        if list(handle.header.samples) != ["SAMPLE"]:
            raise InputValidationError("predictor doctor VCF must contain exactly SAMPLE")
        for record in handle:
            sample = record.samples["SAMPLE"]
            rows.append(
                (
                    record.contig,
                    record.pos,
                    record.id,
                    record.ref,
                    tuple(record.alts or ()),
                    tuple(sample.get("GT") or ()),
                    sample.phased,
                )
            )
            gt = sample.get("GT")
            identity_fields = (
                record.contig,
                str(record.pos),
                record.ref,
                ",".join(record.alts or ()),
                "/".join("." if allele is None else str(allele) for allele in (gt or ())),
                "1" if sample.phased else "0",
                str(record.info.get("GAB_SOURCE_ID", "")),
                str(record.info.get("GAB_SOURCE_DIGEST", "")),
                str(record.info.get("GAB_SOURCE_ORDINAL", "")),
                str(record.info.get("GAB_ALT_INDEX", "")),
            )
            digest.update("\0".join(identity_fields).encode("utf-8"))
            digest.update(b"\n")
    return rows, digest.hexdigest()


def _score_check(plugin: str, expected: float, observed: Any, control_id: str) -> float:
    if not isinstance(observed, (int, float)) or not math.isfinite(observed):
        raise InputValidationError(
            f"predictor doctor {plugin} control {control_id} has no finite score"
        )
    value = float(observed)
    if not math.isclose(value, expected, rel_tol=0.0, abs_tol=1e-6):
        raise InputValidationError(
            f"predictor doctor {plugin} control {control_id} mismatch: "
            f"expected {expected}, observed {value}"
        )
    return value


def evaluate_predictor_controls(
    annotated_vcf: str | Path,
    enabled_plugins: tuple[str, ...],
    *,
    contig_style: str = "bare",
) -> list[dict[str, Any]]:
    """Validate exact public controls and return deterministic evidence rows."""

    annotated_vcf = Path(annotated_vcf)
    covered = {
        plugin
        for control in CONTROL_VARIANTS
        for plugin in control["checks"]
        if plugin in enabled_plugins
    }
    if covered != set(enabled_plugins):
        raise InputValidationError(
            "predictor doctor control coverage mismatch: "
            f"enabled={sorted(enabled_plugins)}, covered={sorted(covered)}"
        )
    expected_by_locus = {
        (
            ("chr" + control["chrom"] if contig_style == "chr" else control["chrom"]),
            control["pos"],
            control["ref"],
            control["alt"],
        ): control
        for control in CONTROL_VARIANTS
    }
    results: list[dict[str, Any]] = []
    seen: set[str] = set()
    with pysam.VariantFile(str(annotated_vcf)) as handle:
        parser = VepParser.from_header(handle.header)
        if not parser.enabled:
            raise InputValidationError("predictor doctor output lacks VEP CSQ")
        parser.require_pipeline_fields()
        parser.require_predictor_fields(enabled_plugins)
        for record in handle:
            if len(record.alts or ()) != 1:
                raise InputValidationError("predictor doctor output is not biallelic")
            locus = (record.contig, record.pos, record.ref, record.alts[0])
            control = expected_by_locus.get(locus)
            if control is None:
                raise InputValidationError(
                    f"predictor doctor output has unexpected allele: {locus}"
                )
            annotation = parser.parse(record)
            if annotation is None:
                raise InputValidationError(
                    f"predictor doctor control {control['control_id']} lacks annotation"
                )
            if annotation["gene"]["symbol"] != control["gene"]:
                raise InputValidationError(
                    f"predictor doctor control {control['control_id']} gene mismatch"
                )
            if annotation["consequence"]["most_severe"] != control["consequence"]:
                raise InputValidationError(
                    f"predictor doctor control {control['control_id']} consequence mismatch"
                )
            checks: dict[str, dict[str, Any]] = {}
            for plugin, expected in control["checks"].items():
                if plugin not in enabled_plugins:
                    continue
                if plugin == "loftee":
                    observed_values = sorted(
                        {
                            value["lof"]
                            for value in annotation["consequence"]["all_annotations"]
                            if value["symbol"] == control["gene"]
                            and control["consequence"] in value["consequence_terms"]
                            and value["lof"] is not None
                        }
                    )
                    if observed_values != [expected]:
                        raise InputValidationError(
                            f"predictor doctor LOFTEE control {control['control_id']} mismatch: "
                            f"expected {[expected]}, observed {observed_values}"
                        )
                    observed: Any = expected
                else:
                    field = {
                        "cadd": "cadd_phred",
                        "revel": "revel",
                        "alphamissense": "alphamissense",
                        "spliceai": "spliceai_max",
                    }[plugin]
                    observed = _score_check(
                        plugin,
                        float(expected),
                        annotation["predictions"][field],
                        control["control_id"],
                    )
                checks[plugin] = {"expected": expected, "observed": observed, "status": "PASS"}
            control_id = control["control_id"]
            if control_id in seen:
                raise InputValidationError(f"predictor doctor duplicated control: {control_id}")
            seen.add(control_id)
            results.append(
                {
                    "control_id": control_id,
                    "variant": f"{record.contig}:{record.pos}:{record.ref}:{record.alts[0]}",
                    "gene": control["gene"],
                    "consequence": control["consequence"],
                    "checks": checks,
                }
            )
    expected_ids = {control["control_id"] for control in CONTROL_VARIANTS}
    if seen != expected_ids:
        raise InputValidationError(
            f"predictor doctor control set is incomplete: missing={sorted(expected_ids - seen)}"
        )
    order = {control["control_id"]: index for index, control in enumerate(CONTROL_VARIANTS)}
    return sorted(results, key=lambda value: order[value["control_id"]])


def _doctor_boundary(config: BuildConfig) -> dict[str, Any]:
    if not config.predictors.manifest.is_file():
        raise ResourceError(f"predictor resource manifest not found: {config.predictors.manifest}")
    if not config.vep_cache_manifest.is_file():
        raise ResourceError(f"VEP cache manifest not found: {config.vep_cache_manifest}")
    reference_sha256 = config.reference_expected_sha256 or sha256_file(config.reference_fasta)
    boundary = {
        "doctor_schema_version": DOCTOR_SCHEMA_VERSION,
        "control_definition_sha256": CONTROL_DEFINITION_SHA256,
        "predictor_manifest_sha256": sha256_file(config.predictors.manifest),
        "reference_fasta_sha256": reference_sha256,
        "contig_style": config.contig_style,
        "vep_version": VEP_VERSION,
        "vep_cache_version": config.vep_cache_version,
        "vep_cache_manifest_sha256": sha256_file(config.vep_cache_manifest),
        "container_image": config.vep_container_image,
        "enabled_plugins": list(config.predictors.enabled_plugins),
    }
    return {**boundary, "doctor_run_key": sha256_json(boundary)}


def predictor_doctor_directory(config: BuildConfig) -> Path:
    boundary = _doctor_boundary(config)
    return config.predictors.root / "semantic-controls" / boundary["doctor_run_key"][:20]


def _expected_report(
    config: BuildConfig,
    vep_summary: dict[str, Any],
    control_results: list[dict[str, Any]],
    input_identity_sha256: str,
) -> dict[str, Any]:
    boundary = _doctor_boundary(config)
    return {
        "schema_version": DOCTOR_SCHEMA_VERSION,
        "status": "PASS",
        "control_set": CONTROL_SET,
        "control_definition_sha256": CONTROL_DEFINITION_SHA256,
        "semantic_boundary": SEMANTIC_BOUNDARY,
        "assembly": "GRCh38",
        **boundary,
        "vep_cache_tree_sha256": vep_summary["cache_tree_sha256"],
        "record_count": len(CONTROL_VARIANTS),
        "record_identity_sha256": input_identity_sha256,
        "variant_and_genotype_identity_unchanged": True,
        "control_results": control_results,
    }


def run_predictor_doctor(config: BuildConfig) -> dict[str, Any]:
    if not config.predictors.enabled:
        return {"status": "not_enabled"}
    destination = predictor_doctor_directory(config)
    if destination.exists():
        return {**verify_predictor_doctor(config), "reused": True}
    with atomic_output_directory(destination) as staging:
        controls = staging / "controls.vcf.gz"
        _write_control_vcf(config, controls)
        input_rows, input_identity_sha256 = _identity(controls)
        vep_summary = run_vep_annotation(config, controls, staging / "vep")
        annotated = staging / "vep" / "annotated.vcf.gz"
        output_rows, output_identity_sha256 = _identity(annotated)
        if output_rows != input_rows or output_identity_sha256 != input_identity_sha256:
            raise InputValidationError(
                "predictor doctor VEP output changed allele or genotype identity"
            )
        results = evaluate_predictor_controls(
            annotated,
            config.predictors.enabled_plugins,
            contig_style=config.contig_style,
        )
        report = _expected_report(config, vep_summary, results, input_identity_sha256)
        atomic_write_json(staging / "predictor_doctor_report.json", report)
    return {**verify_predictor_doctor(config), "reused": False}


def verify_predictor_doctor(config: BuildConfig) -> dict[str, Any]:
    if not config.predictors.enabled:
        return {"status": "not_enabled"}
    destination = predictor_doctor_directory(config)
    controls = destination / "controls.vcf.gz"
    annotated = destination / "vep" / "annotated.vcf.gz"
    summary_path = destination / "vep" / "vep_summary.json"
    report_path = destination / "predictor_doctor_report.json"
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
        vep_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ResourceError(
            "predictor semantic controls are absent; run resources predictors-doctor"
        ) from exc
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ResourceError(f"invalid predictor doctor report: {exc}") from exc
    if not isinstance(report, dict) or not isinstance(vep_summary, dict):
        raise ResourceError("predictor doctor reports must be JSON objects")
    input_rows, input_identity_sha256 = _identity(controls)
    output_rows, output_identity_sha256 = _identity(annotated)
    if input_rows != output_rows or input_identity_sha256 != output_identity_sha256:
        raise ResourceError("predictor doctor allele/genotype identity changed")
    if vep_summary.get("status") != "PASS":
        raise ResourceError("predictor doctor VEP checkpoint is not PASS")
    if vep_summary.get("record_identity_sha256") != input_identity_sha256:
        raise ResourceError("predictor doctor VEP identity digest mismatch")
    if vep_summary.get("vep_version") != VEP_VERSION:
        raise ResourceError("predictor doctor VEP version mismatch")
    if vep_summary.get("cache_version") != config.vep_cache_version:
        raise ResourceError("predictor doctor VEP cache version mismatch")
    if vep_summary.get("container_image") != config.vep_container_image:
        raise ResourceError("predictor doctor VEP container mismatch")
    predictors = vep_summary.get("predictors")
    if not isinstance(predictors, dict) or predictors.get("manifest_sha256") != sha256_file(
        config.predictors.manifest
    ):
        raise ResourceError("predictor doctor resource manifest mismatch")
    results = evaluate_predictor_controls(
        annotated,
        config.predictors.enabled_plugins,
        contig_style=config.contig_style,
    )
    expected = _expected_report(config, vep_summary, results, input_identity_sha256)
    if report != expected:
        raise ResourceError("predictor doctor report differs from recomputed controls")
    return {
        "status": "verified",
        "control_set": CONTROL_SET,
        "control_definition_sha256": CONTROL_DEFINITION_SHA256,
        "doctor_run_key": report["doctor_run_key"],
        "predictor_manifest_sha256": report["predictor_manifest_sha256"],
        "reference_fasta_sha256": report["reference_fasta_sha256"],
        "vep_cache_manifest_sha256": report["vep_cache_manifest_sha256"],
        "enabled_plugins": report["enabled_plugins"],
        "record_count": report["record_count"],
        "record_identity_sha256": report["record_identity_sha256"],
        "control_results": report["control_results"],
        "semantic_boundary": SEMANTIC_BOUNDARY,
        "report_sha256": sha256_file(report_path),
    }
