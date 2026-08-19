from __future__ import annotations

from typing import Any

import pyarrow as pa

from genome_annotation_builder.acmg_evaluator import (
    AcmgSecondaryFindingsEvaluator,
    materialize_acmg_clinical_artifacts,
)
from genome_annotation_builder.acmg_resources import verify_vendor_acmg_snapshot
from genome_annotation_builder.arrow_schemas import (
    clinical_evidence_schema,
    clinical_finding_schema,
    interpretation_schema,
)


def _evidence(gene: str, condition: str, *, stars: int = 3) -> dict[str, Any]:
    return {
        "evidence_id": f"evidence-{gene}-{condition}",
        "clinvar_accession": "VCV000000001",
        "snapshot": "20260808",
        "aggregate_significance": "Pathogenic",
        "review_stars": stars,
        "has_conflicts": False,
        "condition_identifiers": [condition],
        "gene_info": [f"{gene}:1"],
    }


def _observe(
    evaluator: AcmgSecondaryFindingsEvaluator,
    gene: str,
    condition: str,
    *,
    variant_id: str = "chr1:1:A:G",
    gt: tuple[int, ...] = (0, 1),
    stars: int = 3,
    terms: list[str] | None = None,
    hgvsp: str | None = None,
) -> None:
    evaluator.observe(
        variant_id=variant_id,
        gene_symbols={gene},
        gt=gt,
        phased=False,
        call_confidence="high",
        transcripts=[
            {
                "consequence_terms": terms or ["missense_variant"],
                "hgvsp": hgvsp,
                "refseq_transcript_id": None,
                "transcript_id": None,
                "exon": None,
            }
        ],
        clinvar_rows=[_evidence(gene, condition, stars=stars)],
    )


def test_gene_membership_alone_is_not_actionable() -> None:
    evaluator = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    evaluator.observe(
        variant_id="chr17:1:A:G",
        gene_symbols={"BRCA1"},
        gt=(0, 1),
        phased=False,
        call_confidence="high",
        transcripts=[],
        clinvar_rows=[],
    )
    result = evaluator.finalize()
    assert result["rows"] == []
    assert result["actionable_variant_ids"] == []


def test_disease_specific_brca1_pathogenic_variant_is_candidate() -> None:
    evaluator = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    _observe(evaluator, "BRCA1", "OMIM:604370")
    result = evaluator.finalize()
    assert result["actionable_variant_ids"] == ["chr1:1:A:G"]
    assert result["actionable_variant_gene_pairs"] == [
        {"variant_id": "chr1:1:A:G", "gene_symbol": "BRCA1"}
    ]
    assert result["rows"][0]["decision_status"] == "REPORTABLE_CANDIDATE"
    assert result["rows"][0]["clinical_grade"] is False
    assert result["rows"][0]["confirmation_required"] is True


def test_low_review_evidence_does_not_become_actionable() -> None:
    evaluator = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot(), min_review_stars=2)
    _observe(evaluator, "BRCA1", "MIM:604370", stars=1)
    result = evaluator.finalize()
    assert result["actionable_variant_ids"] == []
    assert result["rows"][0]["decision_status"] == "INSUFFICIENT_CLINVAR_REVIEW"


def test_hfe_requires_exact_variant_and_homozygous_genotype() -> None:
    heterozygous = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    _observe(
        heterozygous,
        "HFE",
        "MONDO:MONDO:0021001",
        gt=(0, 1),
        hgvsp="ENSP00000337876.4:p.Cys282Tyr",
    )
    assert heterozygous.finalize()["rows"][0]["decision_status"] == "GENOTYPE_NOT_REPORTABLE"

    homozygous = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    _observe(
        homozygous,
        "HFE",
        "MONDO:0021001",
        gt=(1, 1),
        hgvsp="p.Cys282Tyr",
    )
    assert homozygous.finalize()["actionable_variant_ids"] == ["chr1:1:A:G"]


def test_recessive_compound_heterozygotes_require_phase_review() -> None:
    evaluator = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    _observe(evaluator, "MUTYH", "MONDO:0012041", variant_id="chr1:1:A:G")
    _observe(evaluator, "MUTYH", "MONDO:0012041", variant_id="chr1:2:C:T")
    result = evaluator.finalize()
    assert result["actionable_variant_ids"] == []
    assert {row["decision_status"] for row in result["rows"]} == {"PHASE_REVIEW_REQUIRED"}


def test_recessive_homozygous_variant_can_pass_policy_gate() -> None:
    evaluator = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    _observe(evaluator, "MUTYH", "OMIM:608456", gt=(1, 1))
    assert evaluator.finalize()["actionable_variant_ids"] == ["chr1:1:A:G"]


def test_apob_requires_missense_variant() -> None:
    evaluator = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    _observe(evaluator, "APOB", "MIM:144010", terms=["stop_gained"])
    result = evaluator.finalize()
    assert result["actionable_variant_ids"] == []
    assert result["rows"][0]["decision_status"] == "VARIANT_CLASS_NOT_REPORTABLE"


def test_acmg_decision_materializes_normative_unreviewed_clinical_chain() -> None:
    evaluator = AcmgSecondaryFindingsEvaluator(verify_vendor_acmg_snapshot())
    _observe(evaluator, "BRCA1", "OMIM:604370")
    artifacts = materialize_acmg_clinical_artifacts(evaluator.finalize())

    findings = pa.Table.from_pylist(
        artifacts["clinical_findings"], schema=clinical_finding_schema()
    )
    policy_evidence = pa.Table.from_pylist(
        artifacts["policy_evidence"], schema=clinical_evidence_schema()
    )
    interpretations = pa.Table.from_pylist(
        artifacts["interpretations"], schema=interpretation_schema()
    )
    assert findings.num_rows == policy_evidence.num_rows == interpretations.num_rows == 1
    finding = artifacts["clinical_findings"][0]
    interpretation = artifacts["interpretations"][0]
    assert finding["clinical_grade"] is False
    assert policy_evidence["finding_id"][0].as_py() == finding["finding_id"]
    assert interpretation["finding_id"] == finding["finding_id"]
    assert interpretation["interpretation_status"] == "algorithmic"
    assert interpretation["confirmation_required"] is True
    assert interpretation["reviewed_by"] is None
    assert interpretation["reviewed_at"] is None
