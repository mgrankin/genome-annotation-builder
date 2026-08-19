from __future__ import annotations

import csv
from copy import deepcopy
from pathlib import Path

import duckdb
import pyarrow as pa
import pytest

from genome_annotation_builder import validator
from genome_annotation_builder.arrow_schemas import (
    pgx_gene_call_audit_schema,
    pgx_guidance_audit_schema,
    pgx_locus_audit_schema,
    pharmacogenomics_schema,
    target_callability_evidence_schema,
    variant_schema,
)
from genome_annotation_builder.errors import InputValidationError
from genome_annotation_builder.pgx_materialize import _resolve_pgx_locus
from genome_annotation_builder.pgx_project import (
    CALLS_HEADER,
    _parse_calls_tsv,
    _project_gene_calls,
    _project_guidance,
)
from genome_annotation_builder.pharmcat_runner import (
    DANGEROUS_PHARMCAT_ARGUMENTS,
    _pipeline_arguments,
)


def _rows(*, complex_gene: bool = False) -> list[dict[str, object]]:
    shared: dict[str, object] = {
        "locus_id": "pgx:pharmcat:3.4.0:record:000000001",
        "record_number": 1,
        "source_record_digest": "a" * 64,
        "positions_vcf_sha256": "b" * 64,
        "chrom": "chr1",
        "pos": 10,
        "source_id": "rs1",
        "ref": "C",
        "genes": ["GENE1"],
        "position_of_interest": False,
        "structurally_complex": complex_gene,
    }
    output: list[dict[str, object]] = []
    for allele_index, alt in enumerate(("A", "T"), start=1):
        output.append(
            {
                **shared,
                "allele_index": allele_index,
                "alt": alt,
                "target_id": f"target-{allele_index}",
                "target_variant_id": f"chr1:10:C:{alt}",
                "normalization_changed": False,
                "callability_id": f"callability-{allele_index}",
                "resolution": "unknown",
                "reference_observed": None,
                "callable": None,
                "call_confidence": None,
                "depth": None,
                "genotype_quality": None,
                "observed_ploidy": None,
                "observed_variant_ids": [],
                "variant_gt": None,
                "variant_phased": None,
            }
        )
    return output


def test_pgx_resolver_reconstructs_multiallelic_genotype_slots() -> None:
    rows = _rows()
    rows[0].update(
        resolution="non_reference_observed",
        callable=True,
        call_confidence="high",
        depth=31,
        genotype_quality=45.0,
        observed_ploidy=2,
        observed_variant_ids=["chr1:10:C:A", "chr1:10:C:T"],
        variant_gt=[1, 0],
        variant_phased=True,
    )
    rows[1].update(
        resolution="non_reference_observed",
        callable=True,
        call_confidence="medium",
        depth=29,
        genotype_quality=38.0,
        observed_ploidy=2,
        observed_variant_ids=["chr1:10:C:A", "chr1:10:C:T"],
        variant_gt=[0, 1],
        variant_phased=True,
    )
    result = _resolve_pgx_locus(rows)
    assert result["status"] == "variant_call"
    assert result["genotype"] == [1, 2]
    assert result["phased"] is True
    assert result["depth"] == 29
    assert result["genotype_quality"] == 38.0
    assert result["eligible_for_matcher"] is True


def test_pgx_resolver_uses_only_complete_passing_reference_evidence() -> None:
    rows = _rows()
    for row in rows:
        row.update(
            resolution="reference_observed",
            reference_observed=True,
            callable=True,
            call_confidence="high",
            depth=25,
            genotype_quality=40.0,
            observed_ploidy=2,
        )
    result = _resolve_pgx_locus(rows)
    assert result["status"] == "reference_call"
    assert result["genotype"] == [0, 0]
    assert result["eligible_for_matcher"] is True

    low = deepcopy(rows)
    low[1]["resolution"] = "low_confidence_reference"
    low[1]["callable"] = False
    low[1]["call_confidence"] = "low"
    unresolved = _resolve_pgx_locus(low)
    assert unresolved["status"] == "no_call_low_confidence_reference"
    assert unresolved["genotype"] is None
    assert unresolved["eligible_for_matcher"] is False


def test_pgx_resolver_never_coerces_other_alt_interval_or_unknown_to_reference() -> None:
    other = _rows()
    for row in other:
        row.update(
            resolution="overlapping_variant_observed",
            callable=True,
            call_confidence="high",
            observed_ploidy=2,
            observed_variant_ids=["chr1:10:C:G"],
        )
    result = _resolve_pgx_locus(other)
    assert result["status"] == "no_call_other_variant"
    assert result["genotype"] is None

    interval = _rows()
    for row in interval:
        row.update(
            resolution="callable_without_reference",
            callable=True,
            reference_observed=None,
        )
    result = _resolve_pgx_locus(interval)
    assert result["status"] == "no_call_interval_only"
    assert result["genotype"] is None

    unknown = _resolve_pgx_locus(_rows())
    assert unknown["status"] == "no_call_unknown"
    assert unknown["genotype"] is None


def test_pgx_resolver_suppresses_complex_and_low_confidence_calls() -> None:
    complex_rows = _rows(complex_gene=True)
    for row in complex_rows:
        row.update(
            resolution="reference_observed",
            reference_observed=True,
            callable=True,
            call_confidence="high",
            depth=30,
            genotype_quality=50.0,
            observed_ploidy=2,
        )
    result = _resolve_pgx_locus(complex_rows)
    assert result["status"] == "no_call_structurally_complex"
    assert result["genotype"] is None

    low = _rows()
    low[0].update(
        resolution="non_reference_observed",
        callable=False,
        call_confidence="low",
        depth=4,
        genotype_quality=5.0,
        observed_ploidy=2,
        observed_variant_ids=["chr1:10:C:A"],
        variant_gt=[0, 1],
        variant_phased=False,
    )
    low[1].update(
        resolution="overlapping_variant_observed",
        callable=False,
        call_confidence="low",
        observed_ploidy=2,
        observed_variant_ids=["chr1:10:C:A"],
    )
    result = _resolve_pgx_locus(low)
    assert result["status"] == "no_call_low_confidence_variant"
    assert result["genotype"] is None


def test_pgx_resolver_preserves_haploid_reference_and_variant_calls() -> None:
    reference = _rows()
    for row in reference:
        row.update(
            resolution="reference_observed",
            reference_observed=True,
            callable=True,
            call_confidence="medium",
            depth=18,
            genotype_quality=25.0,
            observed_ploidy=1,
        )
    assert _resolve_pgx_locus(reference)["genotype"] == [0]

    variant = _rows()
    variant[1].update(
        resolution="non_reference_observed",
        callable=True,
        call_confidence="high",
        depth=22,
        genotype_quality=35.0,
        observed_ploidy=1,
        observed_variant_ids=["chr1:10:C:T"],
        variant_gt=[1],
        variant_phased=False,
    )
    variant[0].update(
        resolution="overlapping_variant_observed",
        callable=True,
        call_confidence="high",
        observed_ploidy=1,
        observed_variant_ids=["chr1:10:C:T"],
    )
    assert _resolve_pgx_locus(variant)["genotype"] == [2]


def test_pharmcat_arguments_never_enable_missing_to_reference_or_research_mode() -> None:
    arguments = _pipeline_arguments(["CYP2C19", "TPMT"])
    assert not (set(arguments) & DANGEROUS_PHARMCAT_ARGUMENTS)
    assert "CYP2D6" not in arguments[arguments.index("--genes") + 1].split(",")
    assert "--network" not in arguments  # container hardening is outside tool arguments


def _call(
    gene: str,
    diplotype: str,
    phenotype: str,
    *,
    missing: str = "no",
    outside: str = "no",
) -> dict[str, str]:
    return {
        "Gene": gene,
        "Source Diplotype": diplotype,
        "Phenotype": phenotype,
        "Activity Score": "2.0" if phenotype else "",
        "Haplotype 1": "",
        "Haplotype 1 Function": "",
        "Haplotype 1 Activity Value": "",
        "Haplotype 2": "",
        "Haplotype 2 Function": "",
        "Haplotype 2 Activity Value": "",
        "Outside Call": outside,
        "Match Score": "10",
        "Missing positions": missing,
        "Recommendation Lookup Diplotype": diplotype,
        "Recommendation Lookup Phenotype": phenotype,
        "Recommendation Lookup Activity Score": "2.0" if phenotype else "",
    }


def _pair(pair_id: int, gene: str, drug: str, level: str) -> dict[str, object]:
    return {
        "pair_id": pair_id,
        "drug_name": drug,
        "gene_symbol": gene,
        "guideline_url": f"https://example.test/{pair_id}",
        "cpic_level": level,
        "used_for_recommendation": True,
        "provisional": False,
    }


def test_pgx_projection_publishes_only_unambiguous_complete_noncomplex_calls() -> None:
    calls = [
        _call("CYP2C19", "*1/*1", "Normal Metabolizer"),
        _call("TPMT", "*1/*1 OR *1/*3A", ""),
        _call("UGT1A1", "*1/*1", "Normal Metabolizer", missing="rs123"),
        _call("CYP2D6", "*1/*1", "Normal Metabolizer", outside="yes"),
    ]
    pairs = [
        _pair(1, "CYP2C19", "clopidogrel", "A"),
        _pair(2, "CYP2C19", "amitriptyline", "B"),
        _pair(3, "TPMT", "azathioprine", "A"),
        _pair(4, "UGT1A1", "irinotecan", "A"),
        _pair(5, "CYP2D6", "codeine", "A"),
    ]
    core, audit = _project_gene_calls(
        calls,
        pairs,
        calls_sha256="a" * 64,
        report_sha256="b" * 64,
        data_version="2026-07-13-11-40",
    )
    assert len(core) == 1
    assert core[0]["gene_symbol"] == "CYP2C19"
    assert core[0]["cpic_level"] == "A"
    assert core[0]["affected_drugs"] == ["clopidogrel"]
    statuses = {row["gene_symbol"]: row["call_status"] for row in audit}
    assert statuses == {
        "CYP2C19": "published",
        "CYP2D6": "complex_requires_validated_outside_call",
        "TPMT": "ambiguous_diplotype",
        "UGT1A1": "incomplete_positions",
    }

    report = {
        "drugs": {
            "CPIC Guideline Annotation": {
                "clopidogrel": {
                    "name": "clopidogrel",
                    "source": "CPIC_GUIDELINE",
                    "guidelines": [
                        {
                            "id": "guide-1",
                            "url": "https://example.test/guide-1",
                            "annotations": [
                                {
                                    "lookupKey": [{"CYP2C19": "Normal Metabolizer"}],
                                    "implications": ["normal function"],
                                    "drugRecommendation": "standard guidance",
                                    "classification": "Strong",
                                    "population": "general",
                                }
                            ],
                        }
                    ],
                }
            }
        }
    }
    guidance = _project_guidance(report, pairs, audit, "b" * 64)
    assert len(guidance) == 1
    assert guidance[0]["gene_symbol"] == "CYP2C19"
    assert guidance[0]["cpic_level"] == "A"
    assert guidance[0]["matched"] is True


def test_pgx_projection_audits_complex_genes_excluded_from_matcher() -> None:
    core, audit = _project_gene_calls(
        [],
        [_pair(1, "CYP2D6", "codeine", "A")],
        calls_sha256="a" * 64,
        report_sha256="b" * 64,
        data_version="2026-07-13-11-40",
        required_complex_genes={"CYP2D6"},
    )
    assert core == []
    assert len(audit) == 1
    assert audit[0]["gene_symbol"] == "CYP2D6"
    assert audit[0]["source_diplotypes"] == []
    assert audit[0]["call_status"] == "complex_requires_validated_outside_call"
    assert audit[0]["affected_drugs"] == ["codeine"]


def test_calls_only_parser_enforces_exact_version_and_16_column_boundary(
    tmp_path: Path,
) -> None:
    path = tmp_path / "calls.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["PharmCAT 3.4.0"])
        writer.writerow(CALLS_HEADER)
        row = _call("TPMT", "*1/*1", "Normal Metabolizer")
        writer.writerow([row[column] for column in CALLS_HEADER])
    parsed = _parse_calls_tsv(path)
    assert parsed[0]["Gene"] == "TPMT"

    path.write_text(path.read_text(encoding="utf-8").replace("3.4.0", "3.3.0", 1))
    with pytest.raises(InputValidationError, match="version line mismatch"):
        _parse_calls_tsv(path)


def test_pgx_validator_recomputes_locus_and_projection(monkeypatch: pytest.MonkeyPatch) -> None:
    proof = _rows()[:1]
    proof[0].update(
        resolution="reference_observed",
        reference_observed=True,
        callable=True,
        call_confidence="high",
        depth=30,
        genotype_quality=60.0,
        observed_ploidy=2,
    )
    locus = _resolve_pgx_locus(proof)
    target = {
        "callability_id": "callability-1",
        "target_id": "target-1",
        "target_source": "PharmCAT",
        "target_source_version": "3.4.0@test",
        "catalog_sha256": "c" * 64,
        "catalog_row_digest": "d" * 64,
        "catalog_row_number": 2,
        "chrom": "chr1",
        "pos": 10,
        "end_pos": 10,
        "ref": "C",
        "alt": "A",
        "variant_id": "chr1:10:C:A",
        "gene_symbol": "GENE1",
        "resolution": "reference_observed",
        "reference_observed": True,
        "callable": True,
        "call_confidence": "high",
        "depth": 30,
        "genotype_quality": 60.0,
        "observed_ploidy": 2,
        "assay_scope": "test",
        "evidence_scope": "target_panel",
        "observed_variant_ids": [],
        "observed_variant_source_digests": [],
        "observed_non_core_evidence_ids": [],
        "overlapping_no_call_evidence_ids": [],
        "selected_reference_evidence_ids": ["reference-1"],
        "callable_region_ids": ["region-1"],
        "callable_min_dp": 10,
        "callable_min_gq": 20,
        "precedence_rule": "variant > no-call > explicit-reference > interval > unknown",
        "reason": "synthetic explicit reference proof",
    }
    gene = {
        "gene_symbol": "GENE1",
        "source_diplotypes": ["*1/*1"],
        "selected_diplotype": "*1/*1",
        "phenotypes": ["Normal Metabolizer"],
        "selected_phenotype": "Normal Metabolizer",
        "activity_score_raw": "2.0",
        "activity_score": 2.0,
        "outside_call": False,
        "match_score": "10",
        "missing_positions": [],
        "structurally_complex": False,
        "call_status": "published",
        "cpic_levels": ["A"],
        "affected_drugs": ["example-drug"],
        "guideline_urls": ["https://example.test/guideline"],
        "caller": "PharmCAT",
        "caller_version": "3.4.0",
        "data_version": "test-data",
        "calls_tsv_sha256": "e" * 64,
        "report_json_sha256": "f" * 64,
        "reason": "synthetic complete call",
    }
    core = {
        "gene_symbol": "GENE1",
        "diplotype": "*1/*1",
        "phenotype": "Normal Metabolizer",
        "activity_score": 2.0,
        "copy_number": None,
        "cpic_level": "A",
        "affected_drugs": ["example-drug"],
        "guideline_url": "https://example.test/guideline",
        "caller": "PharmCAT",
        "caller_version": "3.4.0",
    }
    guidance = {
        "guidance_id": "pgx-guidance:" + "1" * 64,
        "source": "CPIC_GUIDELINE",
        "gene_symbol": "GENE1",
        "drug": "example-drug",
        "guideline_id": "guide-1",
        "guideline_url": "https://example.test/guideline",
        "cpic_level": "A",
        "recommendation": "synthetic guidance",
        "implications": ["synthetic implication"],
        "classification": "Strong",
        "population": "test",
        "matched": True,
        "report_json_sha256": "f" * 64,
    }

    monkeypatch.setattr(validator, "PHARMCAT_POSITIONS_RECORDS", 1)
    monkeypatch.setattr(validator, "PHARMCAT_POSITIONS_ALLELES", 1)
    connection = duckdb.connect(":memory:")
    tables = {
        "pgx_locus_audit": pa.Table.from_pylist([locus], schema=pgx_locus_audit_schema()),
        "target_callability_evidence": pa.Table.from_pylist(
            [target], schema=target_callability_evidence_schema()
        ),
        "variants": pa.Table.from_pylist([], schema=variant_schema()),
        "pgx_gene_call_audit": pa.Table.from_pylist([gene], schema=pgx_gene_call_audit_schema()),
        "pharmacogenomics": pa.Table.from_pylist([core], schema=pharmacogenomics_schema()),
        "pgx_guidance_audit": pa.Table.from_pylist([guidance], schema=pgx_guidance_audit_schema()),
    }
    try:
        for name, table in tables.items():
            connection.register(f"source_{name}", table)
            connection.execute(f"CREATE VIEW {name} AS SELECT * FROM source_{name}")
        manifest = {
            "resource_provenance": {
                "pharmacogenomics": {
                    "pharmcat_version": "3.4.0",
                    "pharmcat_release_commit": "1" * 40,
                    "pharmcat_positions_vcf_sha256": "b" * 64,
                    "pharmcat_positions_records": 1,
                    "pharmcat_positions_alleles": 1,
                    "pharmcat_normalization_changed_alleles": 0,
                    "pharmcat_license": "MPL-2.0",
                    "pharmcat_container_image": "example@sha256:" + "2" * 64,
                    "pharmcat_redistribution_allowed": False,
                    "cpic_source": "CPIC PostgREST pair_view",
                    "cpic_source_version": "1.2.3",
                    "cpic_snapshot_date": "2026-08-15",
                    "cpic_snapshot_sha256": "3" * 64,
                    "cpic_rowset_canonical_sha256": "4" * 64,
                    "cpic_source_release": {
                        "tag": "v1.2.3",
                        "commit": "5" * 40,
                        "github_release_id": 42,
                        "published_at": "2026-08-01T00:00:00Z",
                        "release_url": ("https://github.com/cpicpgx/cpic-data/releases/tag/v1.2.3"),
                    },
                    "cpic_api_identity": {
                        "openapi_canonical_sha256": "6" * 64,
                        "postgrest_version": "12.0.2",
                        "swagger_version": "2.0",
                        "base_path": "/v1/",
                    },
                    "cpic_license": "reviewed-test-terms",
                    "cpic_redistribution_allowed": False,
                    "cpic_normative_recommendation_rows": 1,
                    "canonical_report_sha256": "f" * 64,
                    "calls_tsv_sha256": "e" * 64,
                    "missing_to_reference_conversion": False,
                    "network_during_matching": "none",
                    "complex_matcher_genes_excluded": sorted(
                        validator.STRUCTURALLY_COMPLEX_PGX_GENES
                    ),
                }
            }
        }
        errors: list[str] = []
        metrics = validator._validate_pgx(connection, set(tables), manifest, errors)
        tampered_manifest = deepcopy(manifest)
        tampered_manifest["resource_provenance"]["pharmacogenomics"][
            "cpic_rowset_canonical_sha256"
        ] = "not-a-digest"
        tamper_errors: list[str] = []
        validator._validate_pgx(connection, set(tables), tampered_manifest, tamper_errors)
    finally:
        connection.close()
    assert errors == []
    assert metrics["recomputation_mismatches"] == 0
    assert metrics["normative_mismatches"] == 0
    assert any("manifest provenance mismatch" in error for error in tamper_errors)
