"""Deterministic, aggregate-only annotation QC for portable project releases.

The report intentionally contains counts, distributions, completeness fractions,
and immutable reconciliation totals only.  It never exports loci, alleles,
genotypes, subject identifiers, or host paths.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import duckdb

from .portable_model import DOMAINS, PortableEvidenceError

ANNOTATION_QC_SCHEMA_VERSION = 1
QC_IDENTITY_FIELD = "annotation_qc_identity_sha256"
INDEL_BUCKETS = ("1", "2", "3", "4-5", "6-10", "11-20", "21-50", "51-100", "101+")
LOF_CONSEQUENCES = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
)
PREDICTOR_DEFINITIONS = {
    "predictor.alphamissense": "missense_variant rows; populated AlphaMissense score",
    "predictor.cadd": "concrete SNV/INS/DEL rows; populated CADD PHRED score",
    "predictor.loftee": "high-impact loss-of-function consequence rows; populated transcript LOFTEE status",
    "predictor.polyphen": "missense_variant rows; populated PolyPhen score",
    "predictor.revel": "missense_variant rows; populated REVEL score",
    "predictor.sift": "missense_variant rows; populated SIFT score",
    "predictor.spliceai": "concrete SNV/INS/DEL rows; populated SpliceAI maximum score",
}


def _canonical_identity(value: Mapping[str, Any]) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _fraction(numerator: int | float, denominator: int | float) -> float | None:
    if denominator == 0:
        return None
    return round(float(numerator) / float(denominator), 12)


def _load_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise PortableEvidenceError(f"{label} is missing or symlinked")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PortableEvidenceError(f"{label} is unreadable") from exc
    if not isinstance(value, dict):
        raise PortableEvidenceError(f"{label} must be a JSON object")
    return value


def _parquet_members(path: Path, *, required: bool = False) -> list[str]:
    if not path.exists():
        if required:
            raise PortableEvidenceError(f"required QC table is missing: {path.name}")
        return []
    if path.is_symlink():
        raise PortableEvidenceError(f"QC table is symlinked: {path.name}")
    if path.is_file():
        members = [path]
    elif path.is_dir():
        members = sorted(
            member
            for member in path.rglob("*.parquet")
            if member.name not in {"_common_metadata", "_metadata"}
        )
    else:
        raise PortableEvidenceError(f"QC table is not a regular file or directory: {path.name}")
    if not members:
        if required:
            raise PortableEvidenceError(f"required QC table has no Parquet members: {path.name}")
        return []
    if any(member.is_symlink() or not member.is_file() for member in members):
        raise PortableEvidenceError(f"QC table contains an unsafe member: {path.name}")
    return [member.resolve().as_posix() for member in members]


def _register(
    connection: duckdb.DuckDBPyConnection,
    view: str,
    path: Path,
    *,
    required: bool = False,
) -> bool:
    members = _parquet_members(path, required=required)
    if not members:
        return False
    connection.read_parquet(members, hive_partitioning=False).create_view(view)
    return True


def _counter(connection: duckdb.DuckDBPyConnection, query: str) -> dict[str, int]:
    output: dict[str, int] = {}
    for raw_category, raw_count in connection.execute(query).fetchall():
        category = "unknown" if raw_category is None else str(raw_category)
        count = int(raw_count)
        if count < 0 or category in output:
            raise PortableEvidenceError("annotation QC distribution is invalid")
        output[category] = count
    return dict(sorted(output.items()))


def _named_row(connection: duckdb.DuckDBPyConnection, query: str) -> dict[str, int]:
    row = connection.execute(query).fetchone()
    if row is None:
        raise PortableEvidenceError("annotation QC aggregate query returned no row")
    names = [description[0] for description in connection.description]
    return {name: int(value or 0) for name, value in zip(names, row, strict=True)}


def _new_measurements() -> dict[str, Any]:
    return {
        "counts": {},
        "distributions": {},
        "ratios": {},
        "completeness": {},
        "definitions": {
            "small_variants.indel_length": (
                "absolute length difference |len(ALT)-len(REF)| in fixed buckets: "
                + ",".join(INDEL_BUCKETS)
            ),
            "small_variants.ti_tv": "biallelic concrete SNV transitions divided by transversions",
            "vep.gene": "gene symbol or Ensembl gene identifier populated",
            "vep.transcript": "at least one VEP transcript annotation with transcript identifier",
            "vep.hgvs": "top-level or transcript-level HGVSc/HGVSp/HGVSg populated",
            "gnomad.any_af": "at least one typed gnomAD total, popmax, or population AF populated",
            "callability.reference_observed": (
                "explicit site-level reference evidence; interval callability alone is not reference"
            ),
            "callability.unknown": "target resolution is unknown or otherwise lacks a resolved genotype",
            "release_c.adapter_reconciliation": (
                "strict-adapter raw/typed reconciliation counters embedded without loci or alleles"
            ),
            "release_c.portable_scope_omission": (
                "Manta records omitted as whole reciprocal pairs when either endpoint is outside "
                "the portable GRCh38 chr1-22/X/Y/M scope"
            ),
            **PREDICTOR_DEFINITIONS,
        },
    }


def _add_count(
    measurements: dict[str, Any], metric: str, value: int | float, unit: str = "rows"
) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise PortableEvidenceError(f"annotation QC count {metric} is invalid")
    if value < 0 or metric in measurements["counts"]:
        raise PortableEvidenceError(f"annotation QC count {metric} is duplicated or negative")
    measurements["counts"][metric] = {"unit": unit, "value": value}


def _add_distribution(measurements: dict[str, Any], metric: str, values: Mapping[str, int]) -> None:
    if metric in measurements["distributions"]:
        raise PortableEvidenceError(f"annotation QC distribution {metric} is duplicated")
    normalized: dict[str, int] = {}
    for category, count in values.items():
        if not isinstance(category, str) or not category or not isinstance(count, int) or count < 0:
            raise PortableEvidenceError(f"annotation QC distribution {metric} is invalid")
        normalized[category] = count
    measurements["distributions"][metric] = dict(sorted(normalized.items()))


def _add_ratio(
    measurements: dict[str, Any],
    metric: str,
    numerator: int | float,
    denominator: int | float,
    numerator_label: str,
    denominator_label: str,
) -> None:
    if metric in measurements["ratios"] or min(numerator, denominator) < 0:
        raise PortableEvidenceError(f"annotation QC ratio {metric} is invalid")
    measurements["ratios"][metric] = {
        "numerator": numerator,
        "denominator": denominator,
        "numerator_label": numerator_label,
        "denominator_label": denominator_label,
        "value": _fraction(numerator, denominator),
    }


def _add_completeness(
    measurements: dict[str, Any], metric: str, applicable: int, populated: int
) -> None:
    if metric in measurements["completeness"] or not 0 <= populated <= applicable:
        raise PortableEvidenceError(f"annotation QC completeness {metric} is invalid")
    measurements["completeness"][metric] = {
        "applicable_rows": applicable,
        "populated_rows": populated,
        "missing_rows": applicable - populated,
        "fraction": _fraction(populated, applicable),
    }


def _finalize(measurements: dict[str, Any], *, scope: str, bundle_count: int) -> dict[str, Any]:
    document: dict[str, Any] = {
        "schema_version": ANNOTATION_QC_SCHEMA_VERSION,
        "status": "PASS",
        "scope": scope,
        "bundle_count": bundle_count,
        "counts": [
            {"metric": metric, **row} for metric, row in sorted(measurements["counts"].items())
        ],
        "distributions": [
            {"metric": metric, "category": category, "count": count}
            for metric, values in sorted(measurements["distributions"].items())
            for category, count in sorted(values.items())
        ],
        "ratios": [
            {"metric": metric, **row} for metric, row in sorted(measurements["ratios"].items())
        ],
        "completeness": [
            {"metric": metric, **row}
            for metric, row in sorted(measurements["completeness"].items())
        ],
        "definitions": dict(sorted(measurements["definitions"].items())),
    }
    document[QC_IDENTITY_FIELD] = _canonical_identity(document)
    validate_annotation_qc(document)
    return document


def _metric_map(document: Mapping[str, Any], section: str, key: str) -> dict[str, dict[str, Any]]:
    rows = document.get(section)
    if not isinstance(rows, list):
        raise PortableEvidenceError(f"annotation QC {section} is invalid")
    output: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get(key), str) or row[key] in output:
            raise PortableEvidenceError(f"annotation QC {section} contains an invalid row")
        output[row[key]] = row
    return output


def validate_annotation_qc(document: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(document, Mapping):
        raise PortableEvidenceError("annotation QC document is not an object")
    identity = document.get(QC_IDENTITY_FIELD)
    unsigned = {key: value for key, value in document.items() if key != QC_IDENTITY_FIELD}
    if (
        document.get("schema_version") != ANNOTATION_QC_SCHEMA_VERSION
        or document.get("status") != "PASS"
        or document.get("scope") not in {"bundle", "project"}
        or not isinstance(document.get("bundle_count"), int)
        or int(document["bundle_count"]) < 1
        or identity != _canonical_identity(unsigned)
    ):
        raise PortableEvidenceError("annotation QC document identity differs")
    expected_bundle_count = 1 if document.get("scope") == "bundle" else document["bundle_count"]
    if document["scope"] == "bundle" and document["bundle_count"] != expected_bundle_count:
        raise PortableEvidenceError("bundle annotation QC has an invalid bundle count")

    counts = _metric_map(document, "counts", "metric")
    for row in counts.values():
        value = row.get("value")
        if (
            set(row) != {"metric", "unit", "value"}
            or not isinstance(row.get("unit"), str)
            or isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
        ):
            raise PortableEvidenceError("annotation QC count row is invalid")

    distributions = document.get("distributions")
    if not isinstance(distributions, list):
        raise PortableEvidenceError("annotation QC distributions are invalid")
    distribution_map: dict[str, dict[str, int]] = defaultdict(dict)
    for row in distributions:
        if (
            not isinstance(row, dict)
            or set(row) != {"metric", "category", "count"}
            or not isinstance(row.get("metric"), str)
            or not isinstance(row.get("category"), str)
            or not isinstance(row.get("count"), int)
            or row["count"] < 0
            or row["category"] in distribution_map[row["metric"]]
        ):
            raise PortableEvidenceError("annotation QC distribution row is invalid")
        distribution_map[row["metric"]][row["category"]] = row["count"]

    ratios = _metric_map(document, "ratios", "metric")
    for row in ratios.values():
        if set(row) != {
            "metric",
            "numerator",
            "denominator",
            "numerator_label",
            "denominator_label",
            "value",
        }:
            raise PortableEvidenceError("annotation QC ratio row is invalid")
        numerator, denominator = row["numerator"], row["denominator"]
        if (
            isinstance(numerator, bool)
            or isinstance(denominator, bool)
            or not isinstance(numerator, (int, float))
            or not isinstance(denominator, (int, float))
            or min(numerator, denominator) < 0
            or row["value"] != _fraction(numerator, denominator)
        ):
            raise PortableEvidenceError("annotation QC ratio arithmetic differs")

    completeness = _metric_map(document, "completeness", "metric")
    for row in completeness.values():
        applicable = row.get("applicable_rows")
        populated = row.get("populated_rows")
        missing = row.get("missing_rows")
        if (
            set(row) != {"metric", "applicable_rows", "populated_rows", "missing_rows", "fraction"}
            or not all(isinstance(value, int) for value in (applicable, populated, missing))
            or min(applicable, populated, missing) < 0
            or populated + missing != applicable
            or row["fraction"] != _fraction(populated, applicable)
        ):
            raise PortableEvidenceError("annotation QC completeness arithmetic differs")
    definitions = document.get("definitions")
    if not isinstance(definitions, dict) or not all(
        isinstance(key, str) and isinstance(value, str) for key, value in definitions.items()
    ):
        raise PortableEvidenceError("annotation QC definitions are invalid")

    required_counts = {
        "records.source_total",
        "records.accepted_source_alleles",
        "records.normalized_output_alleles",
        "records.output_variant_rows",
        "records.missing_to_reference_conversions",
        "small_variants.total",
        "small_variants.indel_total",
    }
    if not required_counts.issubset(counts):
        raise PortableEvidenceError("annotation QC required record counters are absent")
    total = int(counts["small_variants.total"]["value"])
    type_counts = distribution_map.get("small_variants.variant_type", {})
    zygosity_counts = distribution_map.get("small_variants.zygosity", {})
    if sum(type_counts.values()) != total or sum(zygosity_counts.values()) != total:
        raise PortableEvidenceError("annotation QC variant distributions do not reconcile")
    if counts["records.output_variant_rows"]["value"] != total:
        raise PortableEvidenceError("annotation QC output record count differs from variants")
    if counts["records.missing_to_reference_conversions"]["value"] != 0:
        raise PortableEvidenceError("annotation QC found missing-to-reference conversions")
    ti_tv = ratios.get("small_variants.ti_tv")
    if ti_tv is None or ti_tv["numerator"] + ti_tv["denominator"] != type_counts.get("SNV", 0):
        raise PortableEvidenceError("annotation QC Ti/Tv does not reconcile with SNV count")
    for metric in ("vep.gene", "vep.transcript", "vep.hgvs", "gnomad.any_af"):
        if metric not in completeness or completeness[metric]["applicable_rows"] != total:
            raise PortableEvidenceError(f"annotation QC completeness {metric} is absent")
    return dict(document)


def collect_bundle_annotation_qc(
    bundle: str | Path,
    *,
    validation: Mapping[str, Any],
) -> dict[str, Any]:
    root = Path(bundle).expanduser()
    if root.is_symlink() or not root.is_dir():
        raise PortableEvidenceError("annotation QC bundle root is missing or symlinked")
    root = root.resolve()
    measurements = _new_measurements()

    reconciliation = _load_json(
        root / "extensions" / "reconciliation.json", "variant reconciliation"
    )
    if reconciliation.get("status") != "PASS":
        raise PortableEvidenceError("variant reconciliation did not pass")
    source_counts = reconciliation.get("source_record_counts")
    accepted = reconciliation.get("accepted_source_alleles")
    normalized = reconciliation.get("normalized_output_alleles")
    if not all(isinstance(value, dict) for value in (source_counts, accepted, normalized)):
        raise PortableEvidenceError("variant reconciliation counters are incomplete")

    connection = duckdb.connect(":memory:")
    try:
        # Nested VEP transcript arrays are large enough that an unbounded,
        # parallel correlated UNNEST can consume tens of gigabytes.  Keep the
        # aggregate scan bounded and use list lambdas, which preserve one row
        # per variant instead of materializing transcript cross-products.
        connection.execute("SET threads=2")
        connection.execute("SET memory_limit='4GB'")
        connection.execute("SET preserve_insertion_order=false")
        _register(connection, "variants", root / "variants.parquet", required=True)
        base = _named_row(
            connection,
            """
            SELECT count(*) AS total,
              count(*) FILTER (WHERE variant_type = 'SNV') AS snv,
              count(*) FILTER (WHERE variant_type = 'INS') AS ins,
              count(*) FILTER (WHERE variant_type = 'DEL') AS del,
              count(*) FILTER (WHERE variant_type = 'MNV') AS mnv,
              count(*) FILTER (WHERE variant_type = 'SV') AS sv,
              count(*) FILTER (WHERE genotype.zygosity = 'het') AS het,
              count(*) FILTER (WHERE genotype.zygosity = 'hom_alt') AS hom_alt,
              count(*) FILTER (WHERE genotype.zygosity = 'hom_ref') AS hom_ref,
              count(*) FILTER (WHERE list_count(genotype.gt) = 1) AS haploid,
              count(*) FILTER (
                WHERE variant_type = 'SNV' AND ref || alt IN ('AG','GA','CT','TC')
              ) AS transitions,
              count(*) FILTER (
                WHERE variant_type = 'SNV' AND ref || alt NOT IN ('AG','GA','CT','TC')
              ) AS transversions,
              count(*) FILTER (
                WHERE gene.symbol IS NOT NULL OR gene.ensembl_id IS NOT NULL
              ) AS gene_populated,
              count(*) FILTER (
                WHERE list_count(list_filter(
                  consequence.all_annotations,
                  annotation -> annotation.transcript_id IS NOT NULL
                )) > 0
              ) AS transcript_populated,
              count(*) FILTER (
                WHERE consequence.hgvsc IS NOT NULL OR consequence.hgvsp IS NOT NULL OR
                  list_count(list_filter(
                    consequence.all_annotations,
                    annotation -> annotation.hgvsg IS NOT NULL OR
                                  annotation.hgvsc IS NOT NULL OR
                                  annotation.hgvsp IS NOT NULL
                  )) > 0
              ) AS hgvs_populated,
              count(*) FILTER (WHERE pathogenicity.clinvar_id IS NOT NULL) AS clinvar_matched,
              count(*) FILTER (
                WHERE pathogenicity.clinvar_id IS NOT NULL AND
                      pathogenicity.clinvar_has_conflicts
              ) AS clinvar_conflicts,
              count(*) FILTER (
                WHERE population_frequency.gnomad_af IS NOT NULL OR
                      population_frequency.gnomad_af_popmax IS NOT NULL OR
                      population_frequency.gnomad_af_nfe IS NOT NULL OR
                      population_frequency.gnomad_af_eas IS NOT NULL OR
                      population_frequency.gnomad_af_afr IS NOT NULL OR
                      population_frequency.gnomad_af_sas IS NOT NULL OR
                      population_frequency.gnomad_af_amr IS NOT NULL
              ) AS gnomad_any,
              count(population_frequency.gnomad_af) AS gnomad_af,
              count(population_frequency.gnomad_af_popmax) AS gnomad_af_popmax,
              count(population_frequency.gnomad_af_nfe) AS gnomad_af_nfe,
              count(population_frequency.gnomad_af_eas) AS gnomad_af_eas,
              count(population_frequency.gnomad_af_afr) AS gnomad_af_afr,
              count(population_frequency.gnomad_af_sas) AS gnomad_af_sas,
              count(population_frequency.gnomad_af_amr) AS gnomad_af_amr,
              count(*) FILTER (
                WHERE variant_type IN ('SNV','INS','DEL')
              ) AS cadd_applicable,
              count(*) FILTER (
                WHERE variant_type IN ('SNV','INS','DEL') AND predictions.cadd_phred IS NOT NULL
              ) AS cadd_populated,
              count(*) FILTER (
                WHERE consequence.most_severe = 'missense_variant'
              ) AS missense_applicable,
              count(*) FILTER (
                WHERE consequence.most_severe = 'missense_variant' AND
                      predictions.revel IS NOT NULL
              ) AS revel_populated,
              count(*) FILTER (
                WHERE consequence.most_severe = 'missense_variant' AND
                      predictions.alphamissense IS NOT NULL
              ) AS alphamissense_populated,
              count(*) FILTER (
                WHERE consequence.most_severe = 'missense_variant' AND
                      predictions.polyphen IS NOT NULL
              ) AS polyphen_populated,
              count(*) FILTER (
                WHERE consequence.most_severe = 'missense_variant' AND
                      predictions.sift IS NOT NULL
              ) AS sift_populated,
              count(*) FILTER (
                WHERE variant_type IN ('SNV','INS','DEL')
              ) AS spliceai_applicable,
              count(*) FILTER (
                WHERE variant_type IN ('SNV','INS','DEL') AND
                      predictions.spliceai_max IS NOT NULL
              ) AS spliceai_populated,
              count(*) FILTER (
                WHERE consequence.most_severe IN (
                  'transcript_ablation','splice_acceptor_variant','splice_donor_variant',
                  'stop_gained','frameshift_variant'
                )
              ) AS loftee_applicable,
              count(*) FILTER (
                WHERE consequence.most_severe IN (
                  'transcript_ablation','splice_acceptor_variant','splice_donor_variant',
                  'stop_gained','frameshift_variant'
                ) AND list_count(list_filter(
                  consequence.all_annotations,
                  annotation -> annotation.lof IS NOT NULL
                )) > 0
              ) AS loftee_populated
            FROM variants
            """,
        )
        total = base["total"]
        _add_count(measurements, "small_variants.total", total)
        _add_count(measurements, "small_variants.indel_total", base["ins"] + base["del"])
        _add_count(measurements, "small_variants.haploid", base["haploid"])
        _add_distribution(
            measurements,
            "small_variants.variant_type",
            {
                "SNV": base["snv"],
                "INS": base["ins"],
                "DEL": base["del"],
                "MNV": base["mnv"],
                "SV": base["sv"],
            },
        )
        _add_distribution(
            measurements,
            "small_variants.zygosity",
            {"het": base["het"], "hom_alt": base["hom_alt"], "hom_ref": base["hom_ref"]},
        )
        _add_ratio(
            measurements,
            "small_variants.ti_tv",
            base["transitions"],
            base["transversions"],
            "transitions",
            "transversions",
        )
        _add_completeness(measurements, "vep.gene", total, base["gene_populated"])
        _add_completeness(measurements, "vep.transcript", total, base["transcript_populated"])
        _add_completeness(measurements, "vep.hgvs", total, base["hgvs_populated"])

        consequence_counts = _counter(
            connection,
            "SELECT coalesce(consequence.most_severe, 'missing'), count(*) "
            "FROM variants GROUP BY 1 ORDER BY 1",
        )
        _add_distribution(measurements, "vep.most_severe_consequence", consequence_counts)
        indel_counts = _counter(
            connection,
            """
            SELECT variant_type || ':' || CASE
              WHEN abs(length(alt)-length(ref)) = 1 THEN '1'
              WHEN abs(length(alt)-length(ref)) = 2 THEN '2'
              WHEN abs(length(alt)-length(ref)) = 3 THEN '3'
              WHEN abs(length(alt)-length(ref)) <= 5 THEN '4-5'
              WHEN abs(length(alt)-length(ref)) <= 10 THEN '6-10'
              WHEN abs(length(alt)-length(ref)) <= 20 THEN '11-20'
              WHEN abs(length(alt)-length(ref)) <= 50 THEN '21-50'
              WHEN abs(length(alt)-length(ref)) <= 100 THEN '51-100'
              ELSE '101+' END AS category,
              count(*)
            FROM variants WHERE variant_type IN ('INS','DEL')
            GROUP BY 1 ORDER BY 1
            """,
        )
        for variant_type in ("INS", "DEL"):
            _add_distribution(
                measurements,
                f"small_variants.indel_length.{variant_type}",
                {
                    bucket: indel_counts.get(f"{variant_type}:{bucket}", 0)
                    for bucket in INDEL_BUCKETS
                },
            )

        matched = base["clinvar_matched"]
        _add_distribution(
            measurements,
            "clinvar.match",
            {"matched": matched, "unmatched": total - matched},
        )
        _add_distribution(
            measurements,
            "clinvar.conflict",
            {
                "conflict": base["clinvar_conflicts"],
                "no_conflict": matched - base["clinvar_conflicts"],
            },
        )
        _add_distribution(
            measurements,
            "clinvar.review_stars",
            _counter(
                connection,
                "SELECT coalesce(CAST(pathogenicity.clinvar_review_stars AS VARCHAR), 'unknown'), "
                "count(*) FROM variants WHERE pathogenicity.clinvar_id IS NOT NULL "
                "GROUP BY 1 ORDER BY 1",
            ),
        )
        _add_distribution(
            measurements,
            "clinvar.significance",
            _counter(
                connection,
                "SELECT coalesce(pathogenicity.clinvar_significance, 'unknown'), count(*) "
                "FROM variants WHERE pathogenicity.clinvar_id IS NOT NULL GROUP BY 1 ORDER BY 1",
            ),
        )

        _add_completeness(measurements, "gnomad.any_af", total, base["gnomad_any"])
        for field in (
            "gnomad_af",
            "gnomad_af_popmax",
            "gnomad_af_nfe",
            "gnomad_af_eas",
            "gnomad_af_afr",
            "gnomad_af_sas",
            "gnomad_af_amr",
        ):
            _add_completeness(measurements, f"gnomad.{field}", total, base[field])
        for metric, applicable_key, populated_key in (
            ("predictor.cadd", "cadd_applicable", "cadd_populated"),
            ("predictor.revel", "missense_applicable", "revel_populated"),
            ("predictor.alphamissense", "missense_applicable", "alphamissense_populated"),
            ("predictor.polyphen", "missense_applicable", "polyphen_populated"),
            ("predictor.sift", "missense_applicable", "sift_populated"),
            ("predictor.spliceai", "spliceai_applicable", "spliceai_populated"),
            ("predictor.loftee", "loftee_applicable", "loftee_populated"),
        ):
            _add_completeness(measurements, metric, base[applicable_key], base[populated_key])

        source_total = int(source_counts.get("records_total", -1))
        accepted_count = int(accepted.get("count", -1))
        normalized_count = int(normalized.get("count", -1))
        unique_output = int(reconciliation.get("unique_normalized_variants", -1))
        missing_to_reference = int(reconciliation.get("missing_to_reference_conversions", -1))
        if (
            min(source_total, accepted_count, normalized_count, unique_output, missing_to_reference)
            < 0
        ):
            raise PortableEvidenceError("variant reconciliation has a negative or missing counter")
        if accepted_count != normalized_count or unique_output != total:
            raise PortableEvidenceError("variant reconciliation differs from typed variant rows")
        for metric, value in (
            ("records.source_total", source_total),
            ("records.reference_blocks", int(source_counts.get("records_reference_block", 0))),
            (
                "records.explicit_reference_sites",
                int(source_counts.get("records_explicit_reference_site", 0)),
            ),
            ("records.no_calls", int(source_counts.get("records_no_call", 0))),
            ("records.accepted_source_alleles", accepted_count),
            ("records.normalized_output_alleles", normalized_count),
            ("records.output_variant_rows", total),
            (
                "records.observed_partial_variant_rows",
                int(reconciliation.get("observed_partial_variant_rows", 0)),
            ),
            (
                "records.one_to_many_extra_rows",
                total
                + int(reconciliation.get("observed_partial_variant_rows", 0))
                - normalized_count,
            ),
            ("records.non_core_variant_rows", int(reconciliation.get("non_core_variant_rows", 0))),
            ("records.missing_to_reference_conversions", missing_to_reference),
        ):
            _add_count(measurements, metric, value)

        if _register(
            connection,
            "target_callability",
            root / "extensions" / "target_callability_evidence.parquet",
        ):
            target_total = int(
                connection.execute("SELECT count(*) FROM target_callability").fetchone()[0]
            )
            _add_count(measurements, "callability.target_rows", target_total)
            _add_distribution(
                measurements,
                "callability.status",
                {"available": 1},
            )
            for metric, expression in (
                ("callability.resolution", "resolution"),
                ("callability.evidence_scope", "evidence_scope"),
                (
                    "callability.reference_observed",
                    "CASE WHEN reference_observed IS TRUE THEN 'true' "
                    "WHEN reference_observed IS FALSE THEN 'false' ELSE 'unknown' END",
                ),
                (
                    "callability.callable",
                    "CASE WHEN callable IS TRUE THEN 'true' "
                    "WHEN callable IS FALSE THEN 'false' ELSE 'unknown' END",
                ),
            ):
                _add_distribution(
                    measurements,
                    metric,
                    _counter(
                        connection,
                        f"SELECT {expression}, count(*) FROM target_callability GROUP BY 1 ORDER BY 1",
                    ),
                )
            resolved_reference = int(
                connection.execute(
                    "SELECT count(*) FROM target_callability WHERE reference_observed IS TRUE"
                ).fetchone()[0]
            )
            unknown = int(
                connection.execute(
                    "SELECT count(*) FROM target_callability WHERE resolution = 'unknown'"
                ).fetchone()[0]
            )
            _add_ratio(
                measurements,
                "callability.reference_observed_fraction",
                resolved_reference,
                target_total,
                "reference_observed_rows",
                "target_rows",
            )
            _add_ratio(
                measurements,
                "callability.unknown_fraction",
                unknown,
                target_total,
                "unknown_rows",
                "target_rows",
            )
        else:
            _add_count(measurements, "callability.target_rows", 0)
            _add_distribution(measurements, "callability.status", {"not_available": 1})
            _add_ratio(
                measurements,
                "callability.reference_observed_fraction",
                0,
                0,
                "reference_observed_rows",
                "target_rows",
            )
            _add_ratio(
                measurements,
                "callability.unknown_fraction",
                0,
                0,
                "unknown_rows",
                "target_rows",
            )

        if _register(connection, "callable_regions", root / "callable_regions.parquet"):
            region = connection.execute(
                "SELECT count(*), coalesce(sum(end_pos-start_pos+1),0), "
                "coalesce(sum(end_pos-start_pos+1) FILTER (WHERE callable),0) "
                "FROM callable_regions"
            ).fetchone()
            _add_count(measurements, "callability.region_rows", int(region[0]))
            _add_count(measurements, "callability.region_bases", int(region[1]), "bases")
            _add_count(measurements, "callability.callable_region_bases", int(region[2]), "bases")
        else:
            _add_count(measurements, "callability.region_rows", 0)
            _add_count(measurements, "callability.region_bases", 0, "bases")
            _add_count(measurements, "callability.callable_region_bases", 0, "bases")

        family_available = _register(
            connection, "family_evidence", root / "extensions" / "family_evidence.parquet"
        )
        _add_distribution(
            measurements,
            "family.status",
            {"available" if family_available else "not_available": 1},
        )
        if family_available:
            family_total = int(
                connection.execute("SELECT count(*) FROM family_evidence").fetchone()[0]
            )
            complex_rows = int(
                connection.execute(
                    "SELECT count(*) FROM family_evidence WHERE complex_region"
                ).fetchone()[0]
            )
            _add_count(measurements, "family.rows", family_total)
            _add_count(measurements, "difficult_regions.complex_rows", complex_rows)
            _add_ratio(
                measurements,
                "difficult_regions.complex_fraction",
                complex_rows,
                family_total,
                "complex_rows",
                "family_rows",
            )
            for field in (
                "inheritance_state",
                "parental_origin",
                "mendelian_status",
                "loco_status",
                "caller_concordance",
                "technical_replicate_concordance",
                "phase_confidence",
                "read_evidence_status",
            ):
                _add_distribution(
                    measurements,
                    f"family.{field}",
                    _counter(
                        connection,
                        f"SELECT coalesce(CAST({field} AS VARCHAR), 'unknown'), count(*) "
                        "FROM family_evidence GROUP BY 1 ORDER BY 1",
                    ),
                )
            _add_distribution(
                measurements,
                "family.candidate_flags",
                _named_row(
                    connection,
                    """
                    SELECT count(*) FILTER (WHERE de_novo_candidate) AS de_novo_candidate,
                      count(*) FILTER (WHERE parental_mosaicism_candidate) AS parental_mosaicism_candidate,
                      count(*) FILTER (WHERE parental_dropout_suspected) AS parental_dropout_suspected,
                      count(*) FILTER (WHERE child_dropout_suspected) AS child_dropout_suspected
                    FROM family_evidence
                    """,
                ),
            )
            flags = _counter(
                connection,
                "SELECT flag, count(*) FROM family_evidence, "
                "UNNEST(difficult_region_flags) AS x(flag) GROUP BY 1 ORDER BY 1",
            )
            _add_distribution(measurements, "difficult_regions.flags", flags)
            _add_distribution(
                measurements,
                "difficult_regions.complex_status",
                {"complex": complex_rows, "not_complex": family_total - complex_rows},
            )
            _add_distribution(
                measurements,
                "difficult_regions.variant_type_by_status",
                _counter(
                    connection,
                    "SELECT (CASE WHEN complex_region THEN 'complex:' ELSE 'not_complex:' END) || "
                    "CASE WHEN length(ref)=1 AND length(alt)=1 THEN 'SNV' "
                    "WHEN length(ref)<length(alt) THEN 'INS' "
                    "WHEN length(ref)>length(alt) THEN 'DEL' ELSE 'MNV' END, count(*) "
                    "FROM family_evidence "
                    "GROUP BY 1 ORDER BY 1",
                ),
            )
        else:
            _add_count(measurements, "family.rows", 0)
            _add_count(measurements, "difficult_regions.complex_rows", 0)
            _add_ratio(
                measurements,
                "difficult_regions.complex_fraction",
                0,
                0,
                "complex_rows",
                "family_rows",
            )
            _add_distribution(measurements, "difficult_regions.flags", {})
            _add_distribution(measurements, "difficult_regions.complex_status", {})
            _add_distribution(measurements, "difficult_regions.variant_type_by_status", {})

        pgx_gene = _register(
            connection, "pgx_gene", root / "extensions" / "pgx_gene_call_audit.parquet"
        )
        pgx_locus = _register(
            connection, "pgx_locus", root / "extensions" / "pgx_locus_audit.parquet"
        )
        _add_distribution(
            measurements,
            "pgx.status",
            {"available" if pgx_gene or pgx_locus else "not_run": 1},
        )
        _add_distribution(
            measurements,
            "pgx.gene_call_status",
            _counter(connection, "SELECT call_status, count(*) FROM pgx_gene GROUP BY 1 ORDER BY 1")
            if pgx_gene
            else {},
        )
        _add_distribution(
            measurements,
            "pgx.locus_status",
            _counter(connection, "SELECT status, count(*) FROM pgx_locus GROUP BY 1 ORDER BY 1")
            if pgx_locus
            else {},
        )
        _add_count(
            measurements,
            "pgx.gene_rows",
            int(connection.execute("SELECT count(*) FROM pgx_gene").fetchone()[0])
            if pgx_gene
            else 0,
        )
        _add_count(
            measurements,
            "pgx.locus_rows",
            int(connection.execute("SELECT count(*) FROM pgx_locus").fetchone()[0])
            if pgx_locus
            else 0,
        )

        pgs_models = _register(
            connection, "pgs_models", root / "extensions" / "pgs_model_audit.parquet"
        )
        pgs_loci = _register(
            connection, "pgs_loci", root / "extensions" / "pgs_locus_audit.parquet"
        )
        _add_distribution(
            measurements,
            "prs.status",
            {"available" if pgs_models and pgs_loci else "not_run": 1},
        )
        if pgs_models and pgs_loci:
            model = connection.execute(
                "SELECT count(*), coalesce(sum(expected_loci),0), coalesce(sum(supported_loci),0), "
                "coalesce(sum(strand_flipped_loci),0), coalesce(sum(missing_loci),0), "
                "coalesce(sum(ambiguous_loci),0), coalesce(sum(low_confidence_loci),0), "
                "coalesce(sum(unsupported_loci),0), coalesce(sum(expected_absolute_weight),0), "
                "coalesce(sum(used_absolute_weight),0), count(*) FILTER (WHERE score_valid), "
                "count(*) FILTER (WHERE calibrated) FROM pgs_models"
            ).fetchone()
            names = (
                "model_rows",
                "expected_loci",
                "supported_loci",
                "strand_flipped_loci",
                "missing_loci",
                "ambiguous_loci",
                "low_confidence_loci",
                "unsupported_loci",
                "expected_weight",
                "used_weight",
                "valid_models",
                "calibrated_models",
            )
            model_values = dict(zip(names, model, strict=True))
            for metric in names[:8] + names[10:]:
                _add_count(measurements, f"prs.{metric}", int(model_values[metric]))
            _add_count(
                measurements,
                "prs.expected_absolute_weight",
                float(model_values["expected_weight"]),
                "weight",
            )
            _add_count(
                measurements,
                "prs.used_absolute_weight",
                float(model_values["used_weight"]),
                "weight",
            )
            _add_ratio(
                measurements,
                "prs.weight_completeness",
                float(model_values["used_weight"]),
                float(model_values["expected_weight"]),
                "used_absolute_weight",
                "expected_absolute_weight",
            )
            _add_distribution(
                measurements,
                "prs.publication_status",
                _counter(
                    connection,
                    "SELECT publication_status, count(*) FROM pgs_models GROUP BY 1 ORDER BY 1",
                ),
            )
            _add_distribution(
                measurements,
                "prs.locus_status",
                _counter(connection, "SELECT status, count(*) FROM pgs_loci GROUP BY 1 ORDER BY 1"),
            )
        else:
            for metric in (
                "model_rows",
                "expected_loci",
                "supported_loci",
                "strand_flipped_loci",
                "missing_loci",
                "ambiguous_loci",
                "low_confidence_loci",
                "unsupported_loci",
                "valid_models",
                "calibrated_models",
            ):
                _add_count(measurements, f"prs.{metric}", 0)
            _add_count(measurements, "prs.expected_absolute_weight", 0.0, "weight")
            _add_count(measurements, "prs.used_absolute_weight", 0.0, "weight")
            _add_ratio(
                measurements,
                "prs.weight_completeness",
                0.0,
                0.0,
                "used_absolute_weight",
                "expected_absolute_weight",
            )
            _add_distribution(measurements, "prs.publication_status", {})
            _add_distribution(measurements, "prs.locus_status", {})

        release_c = validation.get("release_c")
        if not isinstance(release_c, Mapping) or release_c.get("status") != "PASS":
            raise PortableEvidenceError("Release-C validation is absent from annotation QC input")
        concordance = release_c.get("concordance")
        source_roles = release_c.get("technical_source_roles")
        if not isinstance(concordance, Mapping) or not isinstance(source_roles, Mapping):
            raise PortableEvidenceError("Release-C validation summaries are incomplete")
        domain_counts: dict[str, int] = {}
        for domain in DOMAINS:
            domain_concordance = concordance.get(domain, {})
            domain_roles = source_roles.get(domain, {})
            if not isinstance(domain_concordance, Mapping) or not isinstance(domain_roles, Mapping):
                raise PortableEvidenceError("Release-C validation domain summary is invalid")
            counts = {str(key): int(value) for key, value in domain_concordance.items()}
            roles = {str(key): int(value) for key, value in domain_roles.items()}
            if sum(counts.values()) != sum(roles.values()):
                raise PortableEvidenceError("Release-C domain summaries do not reconcile")
            if counts:
                domain_counts[domain] = sum(counts.values())
                _add_distribution(measurements, f"release_c.concordance.{domain}", counts)
                _add_distribution(measurements, f"release_c.technical_source_role.{domain}", roles)
                view = f"release_c_{domain}"
                path = root / "extensions" / "release_c" / f"{domain}.parquet"
                if not _register(connection, view, path, required=True):
                    raise PortableEvidenceError("Release-C domain table is absent")
                _add_distribution(
                    measurements,
                    f"release_c.confidence_tier.{domain}",
                    _counter(
                        connection,
                        f"SELECT confidence_tier, count(*) FROM {view} GROUP BY 1 ORDER BY 1",
                    ),
                )
        _add_distribution(measurements, "release_c.domain_rows", domain_counts)
        _add_count(measurements, "release_c.rows", sum(domain_counts.values()))
        adapter_reconciliation = release_c.get("adapter_reconciliation")
        if not isinstance(adapter_reconciliation, Mapping):
            raise PortableEvidenceError("Release-C adapter reconciliation is absent")
        exact = adapter_reconciliation.get("exact_record_reconciliations")
        portable_scope = adapter_reconciliation.get("portable_core_scope")
        if not isinstance(exact, Mapping) or not isinstance(portable_scope, Mapping):
            raise PortableEvidenceError("Release-C adapter reconciliation is incomplete")
        adapter_runs = int(adapter_reconciliation.get("adapter_runs", -1))
        parsing_summaries = int(adapter_reconciliation.get("parsing_summaries", -1))
        adapter_typed_rows = int(adapter_reconciliation.get("typed_rows", -1))
        exact_summaries = int(exact.get("summaries", -1))
        exact_source = int(exact.get("source_records", -1))
        exact_typed = int(exact.get("typed_rows", -1))
        exact_delta = int(exact.get("record_reconciliation_delta", -1))
        scope_runs = int(portable_scope.get("adapter_runs", -1))
        scope_source = int(portable_scope.get("source_records", -1))
        scope_typed = int(portable_scope.get("typed_rows", -1))
        scope_omitted = int(portable_scope.get("omitted_records", -1))
        scope_non_bnd = int(portable_scope.get("omitted_non_bnd_records", -1))
        scope_bnd = int(portable_scope.get("omitted_bnd_records", -1))
        scope_bnd_pairs = int(portable_scope.get("omitted_bnd_pairs", -1))
        if (
            min(
                adapter_runs,
                parsing_summaries,
                adapter_typed_rows,
                exact_summaries,
                exact_source,
                exact_typed,
                exact_delta,
                scope_runs,
                scope_source,
                scope_typed,
                scope_omitted,
                scope_non_bnd,
                scope_bnd,
                scope_bnd_pairs,
            )
            < 0
            or adapter_typed_rows != sum(domain_counts.values())
            or exact_source != exact_typed + exact_delta
            or scope_source != scope_typed + scope_omitted
            or scope_omitted != scope_non_bnd + scope_bnd
            or scope_bnd != 2 * scope_bnd_pairs
        ):
            raise PortableEvidenceError("Release-C adapter reconciliation arithmetic differs")
        for metric, value in (
            ("release_c.adapter_runs", adapter_runs),
            ("release_c.parsing_summaries", parsing_summaries),
            ("release_c.adapter_typed_rows", adapter_typed_rows),
            ("release_c.exact_reconciliation_summaries", exact_summaries),
            ("release_c.exact_reconciliation_source_records", exact_source),
            ("release_c.exact_reconciliation_typed_rows", exact_typed),
            ("release_c.exact_reconciliation_delta", exact_delta),
            ("release_c.portable_scope_adapter_runs", scope_runs),
            ("release_c.portable_scope_source_records", scope_source),
            ("release_c.portable_scope_typed_rows", scope_typed),
            ("release_c.portable_scope_omitted_records", scope_omitted),
            ("release_c.portable_scope_omitted_non_bnd_records", scope_non_bnd),
            ("release_c.portable_scope_omitted_bnd_records", scope_bnd),
            ("release_c.portable_scope_omitted_bnd_pairs", scope_bnd_pairs),
        ):
            _add_count(measurements, metric, value)

        readback_view = "release_c_read_backed_evidence"
        if "read_backed_evidence" in domain_counts:
            for metric, expression in (
                ("read_backed.candidate_category", "candidate_category"),
                ("read_backed.assessment_status", "assessment_status"),
                ("read_backed.variant_type", "variant_type"),
                ("read_backed.local_complexity_status", "local_complexity_status"),
                ("read_backed.local_haplotype_status", "local_haplotype_status"),
            ):
                _add_distribution(
                    measurements,
                    metric,
                    _counter(
                        connection,
                        f"SELECT {expression}, count(*) FROM {readback_view} GROUP BY 1 ORDER BY 1",
                    ),
                )
        else:
            for metric in (
                "read_backed.candidate_category",
                "read_backed.assessment_status",
                "read_backed.variant_type",
                "read_backed.local_complexity_status",
                "read_backed.local_haplotype_status",
            ):
                _add_distribution(measurements, metric, {})
    finally:
        connection.close()

    return _finalize(measurements, scope="bundle", bundle_count=1)


def aggregate_annotation_qc(documents: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if not documents:
        raise PortableEvidenceError("project annotation QC requires bundle reports")
    validated = [validate_annotation_qc(document) for document in documents]
    if any(document.get("scope") != "bundle" for document in validated):
        raise PortableEvidenceError("project annotation QC can aggregate only bundle reports")
    measurements = _new_measurements()
    measurements["definitions"] = {}

    count_values: dict[tuple[str, str], list[int | float]] = defaultdict(list)
    distribution_values: dict[tuple[str, str], int] = defaultdict(int)
    ratio_values: dict[str, dict[str, Any]] = {}
    completeness_values: dict[str, dict[str, int]] = defaultdict(
        lambda: {"applicable": 0, "populated": 0}
    )
    for document in validated:
        for metric, definition in document["definitions"].items():
            previous = measurements["definitions"].setdefault(metric, definition)
            if previous != definition:
                raise PortableEvidenceError("annotation QC definitions differ between bundles")
        for row in document["counts"]:
            count_values[(row["metric"], row["unit"])].append(row["value"])
        for row in document["distributions"]:
            distribution_values[(row["metric"], row["category"])] += row["count"]
        for row in document["ratios"]:
            target = ratio_values.setdefault(
                row["metric"],
                {
                    "numerator": [],
                    "denominator": [],
                    "numerator_label": row["numerator_label"],
                    "denominator_label": row["denominator_label"],
                },
            )
            if (
                target["numerator_label"] != row["numerator_label"]
                or target["denominator_label"] != row["denominator_label"]
            ):
                raise PortableEvidenceError("annotation QC ratio definitions differ")
            target["numerator"].append(row["numerator"])
            target["denominator"].append(row["denominator"])
        for row in document["completeness"]:
            target = completeness_values[row["metric"]]
            target["applicable"] += row["applicable_rows"]
            target["populated"] += row["populated_rows"]

    for (metric, unit), values in sorted(count_values.items()):
        total = math.fsum(float(value) for value in values)
        if all(isinstance(value, int) for value in values):
            total = int(total)
        _add_count(measurements, metric, total, unit)
    grouped_distributions: dict[str, dict[str, int]] = defaultdict(dict)
    for (metric, category), count in sorted(distribution_values.items()):
        grouped_distributions[metric][category] = count
    for metric, values in sorted(grouped_distributions.items()):
        _add_distribution(measurements, metric, values)
    for metric, row in sorted(ratio_values.items()):
        numerator = math.fsum(float(value) for value in row["numerator"])
        denominator = math.fsum(float(value) for value in row["denominator"])
        if all(isinstance(value, int) for value in row["numerator"]):
            numerator = int(numerator)
        if all(isinstance(value, int) for value in row["denominator"]):
            denominator = int(denominator)
        _add_ratio(
            measurements,
            metric,
            numerator,
            denominator,
            row["numerator_label"],
            row["denominator_label"],
        )
    for metric, row in sorted(completeness_values.items()):
        _add_completeness(measurements, metric, row["applicable"], row["populated"])
    return _finalize(measurements, scope="project", bundle_count=len(validated))
