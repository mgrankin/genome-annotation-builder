from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.arrow_schemas import transcript_annotation_type, variant_schema
from genome_annotation_builder.release_c.annotation_qc import (
    aggregate_annotation_qc,
    collect_bundle_annotation_qc,
    validate_annotation_qc,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError
from genome_annotation_builder.variant_build import _blank_variant_row


def _transcript(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {}
    for field in transcript_annotation_type():
        row[field.name] = [] if pa.types.is_list(field.type) else None
    row.update(
        {
            "allele": "G",
            "consequence_terms": ["missense_variant"],
            "feature_type": "Transcript",
            "transcript_id": "ENST000001",
            **overrides,
        }
    )
    return row


def _row(
    *, variant_id: str, pos: int, ref: str, alt: str, gt: tuple[int, ...]
) -> dict[str, object]:
    row = _blank_variant_row(
        variant_id=variant_id,
        chrom="chr1",
        pos=pos,
        ref=ref,
        alt=alt,
        rsid=None,
        call_confidence="high",
        dp=30,
        gq=50,
        gt=gt,
        phased=False,
        allele_balance=0.5 if len(gt) == 2 else 1.0,
    )
    row["consequence"]["all_annotations"] = []
    return row


def _bundle(tmp_path: Path) -> Path:
    root = tmp_path / "synthetic.genome"
    variant_root = root / "variants.parquet"
    variant_root.mkdir(parents=True)

    missense = _row(variant_id="chr1:10:A:G", pos=10, ref="A", alt="G", gt=(0, 1))
    missense["gene"].update(
        {"symbol": "GENE1", "ensembl_id": "ENSG000001", "biotype": "protein_coding"}
    )
    missense["consequence"].update(
        {
            "most_severe": "missense_variant",
            "impact": "MODERATE",
            "hgvsc": "ENST000001:c.1A>G",
            "all_annotations": [_transcript(hgvsc="ENST000001:c.1A>G")],
        }
    )
    missense["population_frequency"].update({"gnomad_af": 0.001, "gnomad_af_popmax": 0.002})
    missense["pathogenicity"].update(
        {
            "clinvar_significance": "Conflicting_classifications",
            "clinvar_has_conflicts": True,
            "clinvar_review_stars": 2,
            "clinvar_id": "123",
        }
    )
    missense["predictions"].update(
        {
            "cadd_phred": 22.0,
            "revel": 0.8,
            "alphamissense": 0.7,
            "polyphen": 0.9,
            "sift": 0.01,
            "spliceai_max": 0.02,
        }
    )

    transversion = _row(variant_id="chr1:20:C:A", pos=20, ref="C", alt="A", gt=(1, 1))
    transversion["consequence"]["most_severe"] = "intergenic_variant"

    insertion = _row(variant_id="chr1:30:A:AT", pos=30, ref="A", alt="AT", gt=(1,))
    insertion["consequence"].update(
        {
            "most_severe": "frameshift_variant",
            "impact": "HIGH",
            "all_annotations": [
                _transcript(
                    allele="AT",
                    consequence_terms=["frameshift_variant"],
                    hgvsg="chr1:g.30dup",
                    lof="HC",
                )
            ],
        }
    )
    insertion["predictions"].update({"cadd_phred": 25.0, "spliceai_max": 0.1})

    table = pa.Table.from_pylist([missense, transversion, insertion], schema=variant_schema())
    pq.write_table(table, variant_root / "part-00000.parquet")
    extension_root = root / "extensions"
    extension_root.mkdir()
    reconciliation = {
        "status": "PASS",
        "accepted_source_alleles": {"count": 3},
        "normalized_output_alleles": {"count": 3},
        "missing_to_reference_conversions": 0,
        "non_core_variant_rows": 0,
        "unique_normalized_variants": 3,
        "source_record_counts": {
            "records_total": 8,
            "records_reference_block": 2,
            "records_explicit_reference_site": 1,
            "records_no_call": 2,
        },
    }
    (extension_root / "reconciliation.json").write_text(
        json.dumps(reconciliation, sort_keys=True), encoding="utf-8"
    )
    return root


def _release_c_validation() -> dict[str, object]:
    return {
        "release_c": {
            "status": "PASS",
            "concordance": {},
            "technical_source_roles": {},
            "adapter_reconciliation": {
                "adapter_runs": 0,
                "parsing_summaries": 0,
                "typed_rows": 0,
                "exact_record_reconciliations": {
                    "summaries": 0,
                    "source_records": 0,
                    "typed_rows": 0,
                    "record_reconciliation_delta": 0,
                },
                "portable_core_scope": {
                    "policy": "grch38_core_chr1_22_x_y_m",
                    "adapter_runs": 0,
                    "source_records": 0,
                    "typed_rows": 0,
                    "omitted_records": 0,
                    "omitted_non_bnd_records": 0,
                    "omitted_bnd_records": 0,
                    "omitted_bnd_pairs": 0,
                    "core_to_noncanonical_bnd_pairs": 0,
                    "noncanonical_to_noncanonical_bnd_pairs": 0,
                },
            },
        }
    }


def _row_by_metric(report: dict[str, object], section: str, metric: str) -> dict[str, object]:
    return next(row for row in report[section] if row["metric"] == metric)


def _distribution(report: dict[str, object], metric: str) -> dict[str, int]:
    return {
        row["category"]: row["count"] for row in report["distributions"] if row["metric"] == metric
    }


def test_annotation_qc_is_aggregate_only_reconciled_and_aggregatable(tmp_path: Path) -> None:
    report = collect_bundle_annotation_qc(_bundle(tmp_path), validation=_release_c_validation())

    assert validate_annotation_qc(report)["status"] == "PASS"
    assert _distribution(report, "small_variants.variant_type") == {
        "DEL": 0,
        "INS": 1,
        "MNV": 0,
        "SNV": 2,
        "SV": 0,
    }
    assert _distribution(report, "small_variants.zygosity") == {
        "het": 1,
        "hom_alt": 2,
        "hom_ref": 0,
    }
    assert _row_by_metric(report, "ratios", "small_variants.ti_tv")["value"] == 1.0
    assert _row_by_metric(report, "completeness", "vep.transcript")["populated_rows"] == 2
    assert _row_by_metric(report, "completeness", "gnomad.any_af")["populated_rows"] == 1
    assert _row_by_metric(report, "completeness", "predictor.revel")["fraction"] == 1.0
    assert _distribution(report, "clinvar.conflict") == {"conflict": 1, "no_conflict": 0}
    assert _distribution(report, "family.status") == {"not_available": 1}
    assert _distribution(report, "prs.status") == {"not_run": 1}

    aggregate = aggregate_annotation_qc([report, report])
    assert aggregate["scope"] == "project"
    assert aggregate["bundle_count"] == 2
    assert _row_by_metric(aggregate, "counts", "small_variants.total")["value"] == 6
    assert _distribution(aggregate, "small_variants.variant_type")["SNV"] == 4
    assert _row_by_metric(aggregate, "ratios", "small_variants.ti_tv")["value"] == 1.0

    tampered = json.loads(json.dumps(report))
    _row_by_metric(tampered, "counts", "small_variants.total")["value"] = 4
    with pytest.raises(PortableEvidenceError, match="identity differs"):
        validate_annotation_qc(tampered)


def test_annotation_qc_rejects_missing_to_reference_conversion(tmp_path: Path) -> None:
    root = _bundle(tmp_path)
    reconciliation_path = root / "extensions" / "reconciliation.json"
    reconciliation = json.loads(reconciliation_path.read_text(encoding="utf-8"))
    reconciliation["missing_to_reference_conversions"] = 1
    reconciliation_path.write_text(json.dumps(reconciliation), encoding="utf-8")

    with pytest.raises(PortableEvidenceError, match="missing-to-reference"):
        collect_bundle_annotation_qc(root, validation=_release_c_validation())
