from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import annotation_qc, bundle_extension, project_release
from genome_annotation_builder.release_c.bundle_extension import build_enriched_bundle
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError

from .test_bundle_extension import _fake_base, _stage


def _validation(rows: int = 2) -> dict[str, object]:
    metrics = {
        "variant_rows": 10,
        "callable_region_rows": 20,
        "target_callability_rows": 3,
        "gene_index_rows": 4,
        "clinical_evidence_rows": 5,
        "clinical_findings_rows": 1,
        "interpretation_rows": 1,
        "prs_rows": 0,
        "pharmacogenomics_rows": 2,
        "gwas_association_rows": 6,
        "reference_rows_checked": 10,
    }
    return {
        "status": "PASS",
        "core": {"metrics": metrics, "warnings": []},
        "release_c": {
            "status": "PASS",
            "domains": 2,
            "rows": rows,
            "input_gates": 2,
            "caller_runs": 2,
            "source_rows_removed": 0,
            "concordance": {
                "structural_variants": {"comparable_concordant": rows},
            },
            "technical_source_roles": {"structural_variants": {"primary": rows}},
        },
        "query": {"status": "PASS", "term": "chr1:100", "rows": 1, "truncated": False},
    }


def _annotation_qc(rows: int = 10, *, gnomad_populated: int | None = None) -> dict[str, object]:
    measurements = annotation_qc._new_measurements()
    for metric, value in (
        ("records.source_total", rows + 5),
        ("records.accepted_source_alleles", rows),
        ("records.normalized_output_alleles", rows),
        ("records.output_variant_rows", rows),
        ("records.missing_to_reference_conversions", 0),
        ("small_variants.total", rows),
        ("small_variants.indel_total", 2),
        ("small_variants.haploid", 0),
        ("pgx.gene_rows", 1),
        ("prs.valid_models", 0),
        ("release_c.rows", 2),
    ):
        annotation_qc._add_count(measurements, metric, value)
    annotation_qc._add_distribution(
        measurements,
        "small_variants.variant_type",
        {"SNV": rows - 2, "INS": 1, "DEL": 1, "MNV": 0, "SV": 0},
    )
    annotation_qc._add_distribution(
        measurements,
        "small_variants.zygosity",
        {"het": rows - 1, "hom_alt": 1, "hom_ref": 0},
    )
    annotation_qc._add_distribution(
        measurements, "clinvar.match", {"matched": 1, "unmatched": rows - 1}
    )
    annotation_qc._add_ratio(
        measurements,
        "small_variants.ti_tv",
        rows - 4,
        2,
        "transitions",
        "transversions",
    )
    annotation_qc._add_ratio(
        measurements,
        "difficult_regions.complex_fraction",
        1,
        rows,
        "complex_rows",
        "family_rows",
    )
    annotation_qc._add_ratio(
        measurements,
        "callability.reference_observed_fraction",
        8,
        10,
        "reference_observed_rows",
        "target_rows",
    )
    annotation_qc._add_ratio(
        measurements,
        "callability.unknown_fraction",
        1,
        10,
        "unknown_rows",
        "target_rows",
    )
    for metric, populated in (
        ("vep.gene", rows - 1),
        ("vep.transcript", rows - 2),
        ("vep.hgvs", rows - 3),
        ("gnomad.any_af", rows - 1 if gnomad_populated is None else gnomad_populated),
    ):
        annotation_qc._add_completeness(measurements, metric, rows, populated)
    return annotation_qc._finalize(measurements, scope="bundle", bundle_count=1)


def _enriched_bundle(
    tmp_path: Path,
    output: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> Path:
    base = _fake_base(tmp_path / "base.genome")
    stage_work = tmp_path / "stage-work"
    stage_work.mkdir(parents=True)
    stage = _stage(stage_work)
    monkeypatch.setattr(
        bundle_extension,
        "validate_bundle",
        lambda *args, **kwargs: {"valid": True, "schema_version": "1.1.0"},
    )
    monkeypatch.setattr(
        bundle_extension,
        "query_bundle",
        lambda *args, **kwargs: {"count": 1, "truncated": False},
    )
    build_enriched_bundle(
        base,
        stage,
        output,
        materialization="copy",
        generated_at="2026-02-03T04:05:06Z",
        validate_base=False,
    )
    normative_schema = (
        Path(project_release.__file__).resolve().parents[1]
        / "vendor"
        / "genome_spec"
        / "schema.json"
    )
    (output / "schema.json").write_bytes(normative_schema.read_bytes())
    return output


def test_project_release_builds_reports_replays_and_detects_tamper(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "release"
    output.mkdir()
    first = _enriched_bundle(
        tmp_path / "first",
        output / "subject-01.genome",
        monkeypatch,
    )
    second = _enriched_bundle(
        tmp_path / "second",
        output / "subject-02.genome",
        monkeypatch,
    )
    reference = tmp_path / "GRCh38.fa"
    reference.write_text(">chr1\nA\n", encoding="ascii")
    monkeypatch.setattr(project_release, "validate_enriched_bundle", lambda *a, **k: _validation())
    monkeypatch.setattr(
        project_release, "collect_bundle_annotation_qc", lambda *a, **k: _annotation_qc()
    )
    bundles = [("subject-01", first), ("subject-02", second)]

    built = project_release.build_project_release(
        bundles,
        output,
        reference_fasta=reference,
        generated_at="2026-02-03T04:05:06Z",
    )

    assert built["status"] == "PASS"
    assert built["bundle_count"] == 2
    assert built["publication"] == {
        "SHA256SUMS": "published",
        "logs": "published",
        "manifests": "published",
        "reports": "published",
    }
    summary = json.loads(
        (output / "reports" / "validation_summary.json").read_text(encoding="utf-8")
    )
    assert summary["status"] == "PASS"
    assert summary["bundle_count"] == 2
    assert {row["bundle_id"] for row in summary["bundles"]} == {
        "subject-01",
        "subject-02",
    }
    combined = b"".join(
        path.read_bytes()
        for root in (output / "reports", output / "manifests", output / "logs")
        for path in sorted(root.rglob("*"))
        if path.is_file()
    )
    assert str(tmp_path).encode() not in combined
    qc_html = (output / "reports" / "annotation_qc.html").read_text(encoding="utf-8")
    assert "Ti/Tv" in qc_html
    assert "gnomAD any AF" in qc_html
    assert "PROJECT" in qc_html
    checksums = (output / "SHA256SUMS").read_text(encoding="ascii")
    assert "logs/project_release.jsonl" in checksums
    assert "reports/bundle_validation/subject-01.json" in checksums
    individual_report = json.loads(
        (output / "reports" / "bundle_validation" / "subject-01.json").read_text(encoding="utf-8")
    )
    assert individual_report["status"] == "PASS"
    assert individual_report["bundle_id"] == "subject-01"
    assert individual_report["errors"] == []
    assert individual_report["annotation_qc"]["status"] == "PASS"
    assert summary["bundles"][0]["annotation_qc_identity_sha256"]

    replay = project_release.build_project_release(
        bundles,
        output,
        reference_fasta=reference,
        generated_at="2026-02-03T04:05:06Z",
    )
    assert replay["status"] == "PASS"
    assert replay["publication"] == {
        "SHA256SUMS": "reused",
        "logs": "reused",
        "manifests": "reused",
        "reports": "reused",
    }
    validated = project_release.validate_project_release(
        bundles,
        output,
        reference_fasta=reference,
    )
    assert validated["status"] == "PASS"
    assert validated["bundles_revalidated"] is True

    log_path = output / "logs" / "project_release.jsonl"
    original_log = log_path.read_bytes()
    log_path.write_bytes(original_log + b"{}\n")
    with pytest.raises(PortableEvidenceError, match="project log differs"):
        project_release.validate_project_release(
            bundles,
            output,
            reference_fasta=reference,
            revalidate_bundles=False,
        )
    log_path.write_bytes(original_log)

    bundle_report_path = output / "reports" / "bundle_validation" / "subject-01.json"
    original_bundle_report = bundle_report_path.read_bytes()
    tampered_bundle_report = json.loads(original_bundle_report)
    tampered_bundle_report["status"] = "FAIL"
    bundle_report_path.write_text(
        json.dumps(tampered_bundle_report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(PortableEvidenceError, match="bundle validation report identity"):
        project_release.validate_project_release(
            bundles,
            output,
            reference_fasta=reference,
            revalidate_bundles=False,
        )
    bundle_report_path.write_bytes(original_bundle_report)

    with (output / "reports" / "annotation_qc.tsv").open("ab") as handle:
        handle.write(b"tamper\n")
    with pytest.raises(PortableEvidenceError, match="annotation QC TSV differs"):
        project_release.validate_project_release(
            bundles,
            output,
            reference_fasta=reference,
            revalidate_bundles=False,
        )


def test_project_release_rejects_duplicate_or_nonportable_bundle_ids(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    with pytest.raises(PortableEvidenceError, match="unique"):
        project_release._safe_bundle_specs([("subject", bundle), ("subject", bundle)])
    with pytest.raises(PortableEvidenceError, match="not portable"):
        project_release._safe_bundle_specs([("subject/path", bundle)])


def test_project_release_requires_explanations_for_cross_bundle_outliers(
    tmp_path: Path,
) -> None:
    records = [
        {"bundle_id": "subject-01", "annotation_qc": _annotation_qc()},
        {"bundle_id": "subject-02", "annotation_qc": _annotation_qc()},
        {
            "bundle_id": "subject-03",
            "annotation_qc": _annotation_qc(gnomad_populated=1),
        },
    ]
    with pytest.raises(PortableEvidenceError, match="outlier requires explanation"):
        project_release._outlier_analysis(records)

    explanation_path = tmp_path / "qc-explanations.json"
    explanation_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "explanations": [
                    {
                        "bundle_id": "subject-03",
                        "metric": "gnomad.any_af_fraction",
                        "code": "EXPECTED_INPUT_CLASS",
                        "explanation": "Synthetic fixture deliberately omits population annotations.",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    explanations = project_release._load_qc_explanations(explanation_path)
    analysis = project_release._outlier_analysis(records, explanations)
    assert analysis["status"] == "PASS"
    assert analysis["flag_count"] == 1
    assert analysis["flags"][0]["explanation_code"] == "EXPECTED_INPUT_CLASS"
