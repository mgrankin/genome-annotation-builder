from __future__ import annotations

from datetime import UTC, datetime, timedelta, timezone
from pathlib import Path

import duckdb
import pyarrow as pa

from genome_annotation_builder.arrow_schemas import clinvar_assertion_schema
from genome_annotation_builder.clinvar_parser import SortedClinVarLookup, aggregate_clinvar
from genome_annotation_builder.query import _clinical, _json_value


def _clinvar_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "clinvar.vcf"
    path.write_text(
        """##fileformat=VCFv4.2
##fileDate=2026-08-08
##reference=GRCh38
##contig=<ID=1,length=100>
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="Allele ID">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Significance">
##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicts">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review">
##INFO=<ID=CLNSIGSCV,Number=.,Type=String,Description="Submissions">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="Conditions">
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Condition IDs">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Genes">
##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP">
##INFO=<ID=ORIGIN,Number=.,Type=Integer,Description="Origin">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="HGVS">
##INFO=<ID=CLNVC,Number=1,Type=String,Description="Type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10\t100\tA\tG\t.\t.\tALLELEID=1000;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNSIGSCV=SCV1;CLNDN=Condition_one;CLNDISDB=MONDO:1;GENEINFO=GENE1:1;RS=123;ORIGIN=1;CLNHGVS=NC_000001.11:g.10A>G;CLNVC=single_nucleotide_variant
1\t10\t101\tA\tG\t.\t.\tALLELEID=1001;CLNSIG=Likely_benign;CLNSIGCONF=Pathogenic|Likely_benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIGSCV=SCV2;CLNDN=Condition_two;CLNDISDB=MONDO:2;GENEINFO=GENE1:1;RS=123;ORIGIN=1
1\t10\t102\tA\tT\t.\t.\tALLELEID=1002;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIGSCV=SCV3
1\t20\t103\tC\tT\t.\t.\tALLELEID=1003;CLNSIG=Benign;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIGSCV=SCV4|SCV5
""",
        encoding="utf-8",
    )
    return path


def test_streaming_exact_allele_match_and_conflict_preservation(tmp_path: Path) -> None:
    with SortedClinVarLookup(_clinvar_vcf(tmp_path), "20260808") as lookup:
        records = lookup.lookup("chr1", 10, "A", "G")
        assert len(records) == 2
        pathogenicity, evidence, rsid = aggregate_clinvar(records, "chr1:10:A:G", "20260808")
        assert pathogenicity == {
            "clinvar_significance": "Conflicting_interpretations",
            "clinvar_has_conflicts": True,
            "clinvar_conflict_summary": "Likely benign; Pathogenic",
            "clinvar_review_stars": 3,
            "clinvar_submitters_count": 2,
            "clinvar_id": "VCV000000100,VCV000000101",
        }
        assert rsid == "rs123"
        assert len(evidence) == 2
        assert evidence[0]["condition_names"] == ["Condition one"]
        assert evidence[0]["source_record_digest"] != evidence[1]["source_record_digest"]
        table = pa.Table.from_pylist(evidence, schema=clinvar_assertion_schema())
        assert table.num_rows == 2

        vus = lookup.lookup("chr1", 10, "A", "T")
        assert (
            aggregate_clinvar(vus, "chr1:10:A:T", "20260808")[0]["clinvar_significance"]
            == "Uncertain_significance"
        )
        benign = lookup.lookup("chr1", 20, "C", "T")
        benign_summary = aggregate_clinvar(benign, "chr1:20:C:T", "20260808")[0]
        assert benign_summary["clinvar_significance"] == "Benign"
        assert benign_summary["clinvar_review_stars"] == 2
        assert benign_summary["clinvar_submitters_count"] == 2


def test_clinical_query_searches_aggregate_variants_and_raw_evidence() -> None:
    connection = duckdb.connect(":memory:")
    try:
        connection.execute(
            """
            CREATE TABLE variants AS SELECT
              'chr1:10:A:G'::VARCHAR AS variant_id,
              {'symbol': 'GENE1'} AS gene,
              {
                'clinvar_significance': 'Conflicting_interpretations',
                'clinvar_has_conflicts': true,
                'clinvar_conflict_summary': 'Likely benign; Pathogenic',
                'clinvar_review_stars': 3::TINYINT,
                'clinvar_id': 'VCV000000100,VCV000000101',
                'is_actionable': false
              } AS pathogenicity,
              false AS clinical_grade
            """
        )
        connection.execute(
            """
            CREATE TABLE acmg_secondary_findings_evidence AS SELECT
              'decision-1'::VARCHAR AS decision_id,
              'chr1:10:A:G'::VARCHAR AS variant_id,
              'GENE1'::VARCHAR AS gene_symbol,
              'Condition one'::VARCHAR AS disease,
              ['MONDO:1']::VARCHAR[] AS disease_identifiers,
              'AD'::VARCHAR AS inheritance,
              'All P and LP'::VARCHAR AS variants_to_report,
              'REPORTABLE_CANDIDATE'::VARCHAR AS decision_status,
              true AS is_actionable,
              false AS clinical_grade,
              true AS confirmation_required,
              'high'::VARCHAR AS analytical_confidence,
              'Pathogenic'::VARCHAR AS clinvar_significance,
              3::TINYINT AS clinvar_review_stars,
              false AS clinvar_has_conflicts,
              ['VCV000000100']::VARCHAR[] AS clinvar_accessions,
              ['MONDO:1']::VARCHAR[] AS matched_condition_identifiers,
              '3.3'::VARCHAR AS policy_version,
              '20260815'::VARCHAR AS policy_snapshot,
              NULL::VARCHAR AS guidance,
              ['test policy gate']::VARCHAR[] AS reasons
            """
        )
        connection.execute(
            """
            CREATE TABLE clinvar_assertions AS SELECT
              'evidence-1'::VARCHAR AS evidence_id,
              'chr1:10:A:G'::VARCHAR AS variant_id,
              '100'::VARCHAR AS clinvar_variation_id,
              'VCV000000100'::VARCHAR AS clinvar_accession,
              ['Pathogenic']::VARCHAR[] AS raw_significance,
              'Pathogenic'::VARCHAR AS aggregate_significance,
              ['reviewed by expert panel']::VARCHAR[] AS review_status,
              3::TINYINT AS review_stars,
              false AS has_conflicts,
              []::VARCHAR[] AS conflict_classifications,
              ['Condition one']::VARCHAR[] AS condition_names,
              ['MONDO:1']::VARCHAR[] AS condition_identifiers,
              ['GENE1:1']::VARCHAR[] AS gene_info,
              ['rs123']::VARCHAR[] AS rsids,
              ['NC_000001.11:g.10A>G']::VARCHAR[] AS hgvs,
              '20260808'::VARCHAR AS snapshot
            """
        )
        connection.execute(
            """
            CREATE TABLE clinical_evidence AS SELECT
              'evidence-1'::VARCHAR AS evidence_id,
              NULL::VARCHAR AS finding_id,
              'chr1:10:A:G'::VARCHAR AS variant_id,
              'NCBI ClinVar'::VARCHAR AS source,
              'VCV000000100'::VARCHAR AS source_record_id,
              '20260808'::VARCHAR AS source_version,
              'Pathogenic'::VARCHAR AS assertion,
              'reviewed by expert panel'::VARCHAR AS review_status,
              TIMESTAMPTZ '2026-08-08 00:00:00+00' AS retrieved_at
            """
        )
        rows, truncated, extra = _clinical(
            connection,
            {
                "variants",
                "clinical_evidence",
                "clinvar_assertions",
                "acmg_secondary_findings_evidence",
            },
            "Conflicting_interpretations",
            10,
        )
        assert truncated is False
        assert [row["_table"] for row in rows] == ["variants"]
        assert extra["tables_searched"] == [
            "extensions/acmg_secondary_findings_evidence",
            "variants",
            "clinical_evidence",
            "extensions/clinvar_assertions",
        ]

        condition_rows, _, condition_extra = _clinical(
            connection,
            {
                "variants",
                "clinical_evidence",
                "clinvar_assertions",
                "acmg_secondary_findings_evidence",
            },
            "Condition one",
            10,
        )
        assert [row["_table"] for row in condition_rows] == [
            "extensions/acmg_secondary_findings_evidence",
            "extensions/clinvar_assertions",
        ]
        assert condition_extra["tables_searched"] == [
            "extensions/acmg_secondary_findings_evidence",
            "variants",
            "clinical_evidence",
            "extensions/clinvar_assertions",
        ]
    finally:
        connection.close()


def test_query_datetime_is_json_safe_and_normalized_to_utc() -> None:
    local = datetime(2026, 8, 15, 1, 0, tzinfo=timezone(timedelta(hours=1)))
    assert _json_value(local) == "2026-08-15T00:00:00Z"
    assert _json_value(datetime(2026, 8, 15, tzinfo=UTC)) == "2026-08-15T00:00:00Z"
