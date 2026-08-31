"""Typed Ensembl VEP CSQ parsing with deterministic transcript selection."""

from __future__ import annotations

import math
import re
from typing import Any
from urllib.parse import unquote

import pysam

from .errors import InputValidationError

# Ordered from most to least severe. The list follows Ensembl/VEP Sequence
# Ontology consequence ordering and includes current splice subtypes.
CONSEQUENCE_SEVERITY = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
)
SEVERITY_RANK = {term: index for index, term in enumerate(CONSEQUENCE_SEVERITY)}
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
SCORE_PATTERN = re.compile(r"^([^()]*)\(([-+0-9.eE]+)\)$")
PIPELINE_REQUIRED_FIELDS = frozenset(
    {
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",
        "BIOTYPE",
        "HGVSc",
        "HGVSp",
        "Existing_variation",
        "CANONICAL",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "APPRIS",
        "TSL",
        "HGNC_ID",
        "SIFT",
        "PolyPhen",
        "DOMAINS",
        "ALLELE_NUM",
        "gnomADg_AF",
        "gnomADe_AF",
    }
)
PREDICTOR_REQUIRED_FIELDS: dict[str, frozenset[str]] = {
    "loftee": frozenset({"LoF", "LoF_filter", "LoF_flags", "LoF_info"}),
    "cadd": frozenset({"CADD_PHRED"}),
    "revel": frozenset({"REVEL"}),
    "alphamissense": frozenset({"am_pathogenicity"}),
    "spliceai": frozenset(
        {
            "SpliceAI_pred_DS_AG",
            "SpliceAI_pred_DS_AL",
            "SpliceAI_pred_DS_DG",
            "SpliceAI_pred_DS_DL",
        }
    ),
}


def _text(value: str | None) -> str | None:
    if value is None or value in {"", ".", "-"}:
        return None
    return unquote(value)


def _float(value: str | None) -> float | None:
    value = _text(value)
    if value is None:
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def _bounded_score(
    value: float | None,
    field: str,
    *,
    minimum: float = 0.0,
    maximum: float | None = 1.0,
) -> float | None:
    if value is None:
        return None
    if value < minimum or (maximum is not None and value > maximum):
        upper = "unbounded" if maximum is None else str(maximum)
        raise InputValidationError(f"VEP {field} score {value} is outside [{minimum}, {upper}]")
    return value


def _integer(value: str | None) -> int | None:
    value = _text(value)
    if value is None:
        return None
    match = re.match(r"^-?\d+", value)
    return int(match.group(0)) if match else None


def _split(value: str | None, separator: str = "&") -> list[str] | None:
    value = _text(value)
    if value is None:
        return None
    values = [unquote(item) for item in value.split(separator) if item]
    return values or None


def _prediction(value: str | None) -> tuple[str | None, float | None]:
    value = _text(value)
    if value is None:
        return None, None
    match = SCORE_PATTERN.fullmatch(value)
    if not match:
        return value, None
    return match.group(1).replace("_", " ") or None, _float(match.group(2))


def _csq_fields(header: pysam.VariantHeader) -> tuple[str, ...] | None:
    if "CSQ" not in header.info:
        return None
    description = header.info["CSQ"].description or ""
    match = re.search(r"(?:Format:\s*)([^\"]+)", description)
    if not match:
        raise InputValidationError("VEP CSQ header has no parseable Format declaration")
    fields = tuple(field.strip() for field in match.group(1).strip().split("|"))
    if not fields or "Consequence" not in fields:
        raise InputValidationError("VEP CSQ header is missing Consequence")
    if len(set(fields)) != len(fields):
        raise InputValidationError("VEP CSQ header contains duplicate field names")
    return fields


def _severity(terms: list[str] | None) -> int:
    return min((SEVERITY_RANK.get(term, len(SEVERITY_RANK)) for term in terms or []), default=999)


def _hgnc(value: str | None) -> str | None:
    value = _text(value)
    if value and value.isdigit():
        return f"HGNC:{value}"
    return value


def _first_float(entries: list[dict[str, str]], aliases: tuple[str, ...]) -> float | None:
    for entry in entries:
        for key in aliases:
            value = _float(entry.get(key))
            if value is not None:
                return value
    return None


def _appris_rank(value: str | None) -> int | None:
    value = (value or "").strip().lower()
    match = re.fullmatch(r"(?:p|principal)([1-5])", value)
    return int(match.group(1)) if match else None


class VepParser:
    def __init__(self, fields: tuple[str, ...] | None, rare_threshold: float = 0.01) -> None:
        self.fields = fields
        self.rare_threshold = rare_threshold

    @classmethod
    def from_header(cls, header: pysam.VariantHeader, rare_threshold: float = 0.01) -> VepParser:
        return cls(_csq_fields(header), rare_threshold=rare_threshold)

    @property
    def enabled(self) -> bool:
        return self.fields is not None

    def require_pipeline_fields(self) -> None:
        missing = sorted(PIPELINE_REQUIRED_FIELDS - set(self.fields or ()))
        if missing:
            raise InputValidationError(
                "VEP CSQ header lacks required pipeline fields: " + ", ".join(missing)
            )

    def require_predictor_fields(self, enabled_plugins: tuple[str, ...]) -> None:
        declared = set(self.fields or ())
        unknown = sorted(set(enabled_plugins) - set(PREDICTOR_REQUIRED_FIELDS))
        if unknown:
            raise InputValidationError(
                "unknown enabled VEP predictor plugins: " + ", ".join(unknown)
            )
        missing_by_plugin = {
            plugin: sorted(PREDICTOR_REQUIRED_FIELDS[plugin] - declared)
            for plugin in enabled_plugins
            if PREDICTOR_REQUIRED_FIELDS[plugin] - declared
        }
        if missing_by_plugin:
            detail = "; ".join(
                f"{plugin}: {', '.join(fields)}"
                for plugin, fields in sorted(missing_by_plugin.items())
            )
            raise InputValidationError("VEP CSQ header lacks predictor fields: " + detail)

    def _entries(self, record: pysam.VariantRecord) -> list[dict[str, str]]:
        if self.fields is None:
            return []
        raw = record.info.get("CSQ")
        if raw is None:
            return []
        values = raw if isinstance(raw, tuple) else (raw,)
        output: list[dict[str, str]] = []
        for value in values:
            parts = str(value).split("|")
            if len(parts) > len(self.fields):
                raise InputValidationError(
                    f"CSQ value has {len(parts)} fields, header declares {len(self.fields)}"
                )
            parts.extend([""] * (len(self.fields) - len(parts)))
            entry = dict(zip(self.fields, parts, strict=True))
            allele_number = _integer(entry.get("ALLELE_NUM"))
            if allele_number != 1:
                raise InputValidationError(
                    f"biallelic VEP CSQ entry has invalid ALLELE_NUM={allele_number!r}"
                )
            output.append(entry)
        return output

    def parse(self, record: pysam.VariantRecord) -> dict[str, Any] | None:
        entries = self._entries(record)
        if not entries:
            return None
        annotations: list[dict[str, Any]] = []
        for entry in entries:
            sift_prediction, sift_score = _prediction(entry.get("SIFT"))
            polyphen_prediction, polyphen_score = _prediction(entry.get("PolyPhen"))
            consequence_terms = _split(entry.get("Consequence")) or []
            canonical_raw = _text(entry.get("CANONICAL"))
            mane_select = _text(entry.get("MANE_SELECT"))
            feature = _text(entry.get("Feature"))
            refseq_transcript = next(
                (
                    value
                    for value in (mane_select, feature)
                    if value and re.match(r"^(?:NM|NR|XM|XR)_", value)
                ),
                None,
            )
            annotations.append(
                {
                    "allele": _text(entry.get("Allele")),
                    "consequence_terms": consequence_terms,
                    "impact": _text(entry.get("IMPACT")),
                    "symbol": _text(entry.get("SYMBOL")),
                    "gene_id": _text(entry.get("Gene")),
                    "feature_type": _text(entry.get("Feature_type")),
                    "transcript_id": feature,
                    "biotype": _text(entry.get("BIOTYPE")),
                    "exon": _text(entry.get("EXON")),
                    "intron": _text(entry.get("INTRON")),
                    "hgvsg": _text(entry.get("HGVSg")),
                    "hgvsc": _text(entry.get("HGVSc")),
                    "hgvsp": _text(entry.get("HGVSp")),
                    "cdna_position": _text(entry.get("cDNA_position")),
                    "cds_position": _text(entry.get("CDS_position")),
                    "protein_position": _text(entry.get("Protein_position")),
                    "amino_acids": _text(entry.get("Amino_acids")),
                    "codons": _text(entry.get("Codons")),
                    "existing_variation": _split(entry.get("Existing_variation")),
                    "distance": _integer(entry.get("DISTANCE")),
                    "strand": _integer(entry.get("STRAND")),
                    "canonical": canonical_raw in {"YES", "1"},
                    "mane_select": mane_select,
                    "mane_plus_clinical": _text(entry.get("MANE_PLUS_CLINICAL")),
                    "appris": _text(entry.get("APPRIS")),
                    "transcript_support_level": _integer(entry.get("TSL")),
                    "transcript_source": _text(entry.get("SOURCE")),
                    "refseq_transcript_id": refseq_transcript,
                    "refseq_match": _text(entry.get("REFSEQ_MATCH")),
                    "ccds": _text(entry.get("CCDS")),
                    "protein_id": _text(entry.get("ENSP")),
                    "hgnc_id": _hgnc(entry.get("HGNC_ID")),
                    "symbol_source": _text(entry.get("SYMBOL_SOURCE")),
                    "sift_prediction": sift_prediction,
                    "sift_score": sift_score,
                    "polyphen_prediction": polyphen_prediction,
                    "polyphen_score": polyphen_score,
                    "domains": _split(entry.get("DOMAINS")),
                    "lof": _text(entry.get("LoF")),
                    "lof_filter": _text(entry.get("LoF_filter")),
                    "lof_flags": _text(entry.get("LoF_flags")),
                    "lof_info": _text(entry.get("LoF_info")),
                    "flags": _split(entry.get("FLAGS")),
                    "variant_class": _text(entry.get("VARIANT_CLASS")),
                }
            )

        def transcript_selection_key(annotation: dict[str, Any]) -> tuple[Any, ...]:
            if annotation["mane_select"]:
                transcript_priority = 0
            elif annotation["mane_plus_clinical"]:
                transcript_priority = 1
            elif annotation["canonical"]:
                transcript_priority = 2
            elif _appris_rank(annotation["appris"]) is not None:
                transcript_priority = 3
            else:
                transcript_priority = 4
            return (
                transcript_priority,
                _appris_rank(annotation["appris"]) or 99,
                _severity(annotation["consequence_terms"]),
                IMPACT_RANK.get(annotation["impact"], 99),
                annotation["transcript_id"] or "~",
            )

        # Transcript quality is meaningful within one gene, but comparing it
        # globally can make a nearby gene's MANE/APPRIS transcript override a
        # coding consequence in the gene actually affected by the variant.
        # First choose one representative transcript per gene, then select the
        # primary gene by the consequence of those representatives.
        gene_annotations: dict[tuple[str, str], list[dict[str, Any]]] = {}
        for annotation in annotations:
            if annotation["gene_id"]:
                gene_key = ("gene_id", annotation["gene_id"])
            elif annotation["symbol"]:
                gene_key = ("symbol", annotation["symbol"])
            else:
                gene_key = ("feature", annotation["transcript_id"] or "~")
            gene_annotations.setdefault(gene_key, []).append(annotation)

        gene_representatives = [
            min(values, key=transcript_selection_key) for values in gene_annotations.values()
        ]
        selected = min(
            gene_representatives,
            key=lambda annotation: (
                _severity(annotation["consequence_terms"]),
                IMPACT_RANK.get(annotation["impact"], 99),
                transcript_selection_key(annotation),
                annotation["gene_id"] or "~",
                annotation["symbol"] or "~",
            ),
        )
        severe_annotation = min(
            annotations,
            key=lambda value: (
                _severity(value["consequence_terms"]),
                IMPACT_RANK.get(value["impact"], 99),
                value["transcript_id"] or "~",
            ),
        )
        all_terms = [term for annotation in annotations for term in annotation["consequence_terms"]]
        most_severe = min(
            all_terms, key=lambda term: (SEVERITY_RANK.get(term, 999), term), default=None
        )
        existing = sorted(
            {
                identifier
                for annotation in annotations
                for identifier in annotation["existing_variation"] or []
            }
        )
        rsid = next(
            (identifier for identifier in existing if re.fullmatch(r"rs\d+", identifier)), None
        )

        sift_values = [
            value["sift_score"] for value in annotations if value["sift_score"] is not None
        ]
        polyphen_values = [
            value["polyphen_score"] for value in annotations if value["polyphen_score"] is not None
        ]
        gnomadg_af = _first_float(entries, ("gnomADg_AF",))
        gnomade_af = _first_float(entries, ("gnomADe_AF", "gnomAD_AF"))
        gnomad_af = gnomadg_af if gnomadg_af is not None else gnomade_af
        population_aliases = {
            "gnomad_af_nfe": ("gnomADg_NFE_AF", "gnomADe_NFE_AF", "gnomAD_NFE_AF"),
            "gnomad_af_eas": ("gnomADg_EAS_AF", "gnomADe_EAS_AF", "gnomAD_EAS_AF"),
            "gnomad_af_afr": ("gnomADg_AFR_AF", "gnomADe_AFR_AF", "gnomAD_AFR_AF"),
            "gnomad_af_sas": ("gnomADg_SAS_AF", "gnomADe_SAS_AF", "gnomAD_SAS_AF"),
            "gnomad_af_amr": ("gnomADg_AMR_AF", "gnomADe_AMR_AF", "gnomAD_AMR_AF"),
        }
        populations = {
            target: _first_float(entries, aliases) for target, aliases in population_aliases.items()
        }
        all_population_aliases = (
            "gnomADg_AFR_AF",
            "gnomADg_AMI_AF",
            "gnomADg_AMR_AF",
            "gnomADg_ASJ_AF",
            "gnomADg_EAS_AF",
            "gnomADg_FIN_AF",
            "gnomADg_MID_AF",
            "gnomADg_NFE_AF",
            "gnomADg_REMAINING_AF",
            "gnomADg_SAS_AF",
        )
        if gnomadg_af is None:
            all_population_aliases = tuple(
                value.replace("gnomADg_", "gnomADe_") for value in all_population_aliases
            )
        population_values = [
            value
            for key in all_population_aliases
            if (value := _first_float(entries, (key,))) is not None
        ]
        popmax = max(population_values) if population_values else gnomad_af

        cadd = _first_float(entries, ("CADD_PHRED", "CADD_PHRED_score"))
        revel = _first_float(entries, ("REVEL", "REVEL_score"))
        alphamissense = _first_float(
            entries, ("am_pathogenicity", "AlphaMissense", "AlphaMissense_score")
        )
        spliceai_values = [
            _first_float(entries, (key,))
            for key in (
                "SpliceAI_pred_DS_AG",
                "SpliceAI_pred_DS_AL",
                "SpliceAI_pred_DS_DG",
                "SpliceAI_pred_DS_DL",
            )
        ]
        spliceai = max((value for value in spliceai_values if value is not None), default=None)
        cadd = _bounded_score(cadd, "CADD_PHRED", maximum=None)
        revel = _bounded_score(revel, "REVEL")
        alphamissense = _bounded_score(alphamissense, "AlphaMissense")
        spliceai = _bounded_score(spliceai, "SpliceAI")
        sift_values = [_bounded_score(value, "SIFT") for value in sift_values]
        polyphen_values = [_bounded_score(value, "PolyPhen") for value in polyphen_values]
        frequency_values = {
            "gnomad_af": gnomad_af,
            "gnomad_af_popmax": popmax,
            **populations,
        }
        for field, value in frequency_values.items():
            _bounded_score(value, field)
        damaging_votes = sum(
            (
                revel is not None and revel > 0.75,
                alphamissense is not None and alphamissense > 0.564,
                cadd is not None and cadd > 20.0,
            )
        )
        return {
            "rsid": rsid,
            "gene": {
                "symbol": selected["symbol"],
                "ensembl_id": selected["gene_id"],
                "hgnc_id": selected["hgnc_id"],
                "biotype": selected["biotype"],
                "is_canonical": selected["canonical"],
            },
            "consequence": {
                "most_severe": most_severe,
                "impact": severe_annotation["impact"],
                "hgvsc": selected["hgvsc"],
                "hgvsp": selected["hgvsp"],
                "all_annotations": annotations,
            },
            "population_frequency": {
                "gnomad_af": gnomad_af,
                "gnomad_af_popmax": popmax,
                **populations,
                "is_rare": popmax is None or popmax < self.rare_threshold,
            },
            "predictions": {
                "cadd_phred": cadd,
                "revel": revel,
                "alphamissense": alphamissense,
                "spliceai_max": spliceai,
                "polyphen": max(polyphen_values) if polyphen_values else None,
                "sift": min(sift_values) if sift_values else None,
                "consensus_damaging": most_severe == "missense_variant" and damaging_votes >= 2,
            },
        }


def apply_vep_annotation(row: dict[str, Any], annotation: dict[str, Any] | None) -> None:
    if annotation is None:
        return
    if row.get("rsid") is None and annotation.get("rsid"):
        row["rsid"] = annotation["rsid"]
    for key in ("gene", "consequence", "population_frequency", "predictions"):
        row[key].update(annotation[key])
