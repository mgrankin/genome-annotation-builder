"""A full build may resume only stages published for that build identity."""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from test_vcf_core import _environment

from genome_annotation_builder import pipeline
from genome_annotation_builder.errors import GenomeBuilderError
from genome_annotation_builder.resources import sync_pinned_resources


def test_build_rejects_transplanted_stages_and_resumes_its_own(tmp_path, monkeypatch):
    config, _ = _environment(tmp_path, pgs=True)
    sync_pinned_resources(config)
    built = pipeline.build_bundle(config)
    run_key = built["run_key"]
    run = config.work_dir / f"run-{run_key}"
    retry = replace(config, bundle_name="resumed.genome")
    summaries = sorted(run.glob("[0-9][0-9]-*/*_summary.json"))
    assert {p.parent.name for p in summaries} == {
        "01-ingest",
        "02-normalization",
        "04-variant-datasets",
        "05-pgs-resources",
        "07-target-callability",
        "08-pgs-scoring",
    }
    for index, path in enumerate(summaries):
        content = path.read_bytes()
        receipt = json.loads(content)
        assert receipt.pop("build_run_key") == run_key
        if index % 2:
            receipt["build_run_key"] = "0" * 20
        path.write_text(json.dumps(receipt))
        with pytest.raises(GenomeBuilderError, match="checkpoint belongs to a different") as exc:
            pipeline.build_bundle(retry)
        assert str(path.parent) in str(exc.value) and not retry.bundle_path.exists()
        path.write_bytes(content)

    def unexpected(*args, **kwargs):
        pytest.fail("an already completed stage was rebuilt")

    for func in (
        "ingest_reference_confidence_vcf",
        "normalize_projected_vcf",
        "build_variant_datasets",
        "prepare_pgs_resources",
        "build_target_callability",
        "score_pgs",
    ):
        monkeypatch.setattr(pipeline, func, unexpected)
    resumed = pipeline.build_bundle(retry)
    assert resumed["validation"]["valid"] and Path(resumed["bundle"]) == retry.bundle_path
