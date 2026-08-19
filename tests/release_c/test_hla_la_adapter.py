from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.hla_la_adapter import (
    HlaLaProvenance,
    materialize_hla_la_bestguess,
    parse_hla_la_bestguess,
)
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError, schema_for

REFERENCE = "1" * 64
PROVENANCE = HlaLaProvenance(
    reference_sha256=REFERENCE,
    technical_source_role="primary",
    caller_run_identity_sha256="2" * 64,
    resource_manifest_sha256="3" * 64,
    source_alignment_sha256="4" * 64,
    raw_artifact_id="hlala-bestguess-g-group",
    caller_database_name="HLA*LA PRG MHC graph",
    caller_database_version="publisher graph snapshot 2017",
    caller_database_identity_sha256="5" * 64,
    nomenclature_database_name="IMGT nomenclature embedded in HLA*LA graph",
    nomenclature_version="publisher graph snapshot 2017",
    nomenclature_identity_sha256="6" * 64,
    database_age_limitation=True,
)


HLA_LA_COLUMNS = (
    "Locus\tChromosome\tAllele\tQ1\tQ2\tAverageCoverage\tCoverageFirstDecile\t"
    "MinimumCoverage\tproportionkMersCovered\tLocusAvgColumnError\t"
    "NColumns_UnaccountedAllele_fGT0.2\tperfectG"
)


def control_table(tmp_path: Path) -> Path:
    """Create a fully synthetic 13-locus parser control at test time."""

    calls = {
        "A": ("A*11:01:01G", "A*01:01:01G"),
        "B": ("B*08:01:01G", "B*56:01:01G"),
        "C": ("C*01:02:01G", "C*07:01:01G"),
        "DPA1": ("DPA1*01:03:01G", "DPA1*02:01:02"),
        "DPB1": ("DPB1*04:01:01G", "DPB1*14:01:01"),
        "DQA1": ("DQA1*01:01:01G", "DQA1*05:01:01G"),
        "DQB1": ("DQB1*02:01:01G", "DQB1*05:01:01G"),
        "DRB1": ("DRB1*01:01:01G", "DRB1*03:01:01G"),
        "DRB3": ("DRB3*01:01:02G", "DRB3*01:01:02G"),
        "DRB4": ("DRB4*03:01N", "DRB4*03:01N"),
        "E": ("E*01:01:01G", "E*01:01:01G"),
        "F": ("F*01:01:01G", "F*01:01:01G"),
        "G": ("G*01:01:03G", "G*01:01:02G"),
    }
    rows = [HLA_LA_COLUMNS]
    for locus, alleles in calls.items():
        for copy, allele in enumerate(alleles, start=1):
            q1 = "0.98459" if (locus, copy) == ("DPA1", 2) else "1"
            perfect_g = "0" if (locus, copy) in {("DPA1", 2), ("DPB1", 2)} else "1"
            # HLA*LA documents this as an allele-specific metric.  A null
            # allele can therefore carry the not-assessed sentinel while the
            # other called allele at the same locus has complete k-mer
            # coverage.
            proportion = "-1" if (locus, copy) == ("DRB4", 2) else "1"
            if locus == "G" and copy == 2:
                q1 = "0.99"
            rows.append(
                "\t".join(
                    (
                        locus,
                        str(copy),
                        allele,
                        q1,
                        "-5",
                        "49.1975",
                        "42",
                        "37",
                        proportion,
                        "0.001",
                        "0",
                        perfect_g,
                    )
                )
            )
    path = tmp_path / "synthetic-hlala-bestguess.txt"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_synthetic_control_preserves_all_calls_metrics_and_age_limitation(
    tmp_path: Path,
) -> None:
    source = control_table(tmp_path)
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    result_one = materialize_hla_la_bestguess(source, first, PROVENANCE)
    result_two = materialize_hla_la_bestguess(source, second, PROVENANCE)
    parsing = result_one["parsing"]
    assert parsing["source_allele_rows"] == parsing["reconciled_allele_slots"] == 26
    assert parsing["typed_loci"] == 13
    assert parsing["primary_loci"] == 8
    assert parsing["additional_loci"] == 5
    assert parsing["supported_loci"] == 9
    assert parsing["review_loci"] == 4
    assert parsing["unresolved_loci"] == 0
    assert result_one["output"]["parquet_sha256"] == result_two["output"]["parquet_sha256"]
    table = pq.read_table(first)
    assert table.schema == schema_for("hla_typing")
    rows = table.to_pylist()
    assert all(row["database_age_limitation"] for row in rows)
    assert all(row["phase_resolved"] is False for row in rows)
    a_locus = next(row for row in rows if row["locus"] == "HLA-A")
    assert (a_locus["allele1"], a_locus["allele2"]) == (
        "A*11:01:01G",
        "A*01:01:01G",
    )
    assert a_locus["allele1_two_field"] == "A*11:01"
    assert a_locus["allele1_hlala_metrics"]["q1"] == pytest.approx(1.0)
    assert a_locus["allele1_hlala_metrics"]["average_coverage"] == pytest.approx(49.1975)
    assert a_locus["confidence_tier"] == "SUPPORTED"
    dpa1 = next(row for row in rows if row["locus"] == "HLA-DPA1")
    assert dpa1["confidence_tier"] == "REVIEW"
    assert "ALLELE2:Q1_BELOW_ONE" in dpa1["filter_reasons"]
    assert "ALLELE2:PERFECT_G_FALSE" in dpa1["filter_reasons"]
    drb4 = next(row for row in rows if row["locus"] == "HLA-DRB4")
    assert drb4["allele1_hlala_metrics"]["proportion_kmers_covered"] == 1
    assert drb4["allele2_hlala_metrics"]["proportion_kmers_covered"] == -1
    assert "ALLELE2:KMER_COVERAGE_NOT_ASSESSED" in drb4["filter_reasons"]
    rendered = json.dumps({"rows": rows, "summary": parsing}, sort_keys=True)
    assert "/raid/" not in rendered
    assert "NA12878" not in rendered


def test_duplicate_missing_and_invalid_perfect_g_fail_closed(tmp_path: Path) -> None:
    source = control_table(tmp_path)
    lines = source.read_text(encoding="utf-8").splitlines()
    duplicate = tmp_path / "duplicate.txt"
    duplicate.write_text("\n".join([*lines, lines[1]]) + "\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match=r"duplicate.*locus-copy"):
        parse_hla_la_bestguess(duplicate, PROVENANCE)

    missing = tmp_path / "missing.txt"
    missing.write_text("\n".join([lines[0], *lines[2:]]) + "\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="lacks one diploid copy"):
        parse_hla_la_bestguess(missing, PROVENANCE)

    invalid = tmp_path / "invalid-perfect-g.txt"
    fields = lines[1].split("\t")
    fields[-1] = "2"
    invalid.write_text("\n".join([lines[0], "\t".join(fields), *lines[2:]]) + "\n")
    with pytest.raises(PortableEvidenceError, match="perfectG must be zero or one"):
        parse_hla_la_bestguess(invalid, PROVENANCE)


def test_repeated_true_locus_metric_disagreement_fails_closed(tmp_path: Path) -> None:
    source = control_table(tmp_path)
    lines = source.read_text(encoding="utf-8").splitlines()
    fields = lines[2].split("\t")
    fields[5] = "48.0"
    mismatch = tmp_path / "mismatched-locus-depth.txt"
    mismatch.write_text("\n".join([lines[0], lines[1], "\t".join(fields), *lines[3:]]) + "\n")
    with pytest.raises(PortableEvidenceError, match="repeated locus metrics disagree"):
        parse_hla_la_bestguess(mismatch, PROVENANCE)


def test_age_boundary_and_private_binding_fail_closed(tmp_path: Path) -> None:
    source = control_table(tmp_path)
    no_age_limit = HlaLaProvenance(**{**PROVENANCE.__dict__, "database_age_limitation": False})
    with pytest.raises(PortableEvidenceError, match="database-age limitation"):
        parse_hla_la_bestguess(source, no_age_limit)

    private_identifier = "PRIVATE_SUBJECT"
    leaking = HlaLaProvenance(
        **{
            **PROVENANCE.__dict__,
            "raw_artifact_id": f"hlala-{private_identifier}",
        }
    )
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        parse_hla_la_bestguess(source, leaking, forbidden_bindings=(private_identifier,))
