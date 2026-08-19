from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import bundle_extension
from genome_annotation_builder.release_c.bundle_extension import (
    DOMAIN_OMISSIONS,
    GENERIC_READ_BACKED_OMISSION,
    PORTABLE_CORE_SCOPE_OMISSION,
    build_enriched_bundle,
    validate_enriched_bundle,
)
from genome_annotation_builder.release_c.extension_stage import stage_release_c_extensions
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError
from genome_annotation_builder.util import directory_digest

from .test_extension_stage import complete_stage_inputs


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")


def _fake_base(root: Path) -> Path:
    (root / "extensions").mkdir(parents=True)
    (root / "variants.parquet").mkdir()
    (root / "variants.parquet" / "synthetic.bin").write_bytes(b"synthetic variants\n")
    (root / "extensions" / "base-evidence.bin").write_bytes(b"base evidence\n")
    (root / "schema.json").write_text('{"version":"1.1.0"}\n', encoding="utf-8")
    (root / "readmygenome" / "scripts").mkdir(parents=True)
    (root / "readmygenome" / "reference").mkdir()
    (root / "readmygenome" / "evals").mkdir()
    (root / "readmygenome.md").write_text("# Synthetic bundle\n", encoding="utf-8")
    (root / "readmygenome" / "SKILL.md").write_text("# Synthetic skill\n", encoding="utf-8")
    (root / "readmygenome" / "scripts" / "query.py").write_text(
        "# old synthetic query\n", encoding="utf-8"
    )
    (root / "readmygenome" / "reference" / "query_catalog.md").write_text(
        "# Synthetic catalogue\n", encoding="utf-8"
    )
    _write_json(root / "readmygenome" / "reference" / "synonyms.json", {})
    _write_json(root / "readmygenome" / "evals" / "canonical_questions.json", {"cases": []})
    _write_json(
        root / "raw_artifacts.json",
        {"schema_version": "1.1.0", "artifacts": []},
    )
    _write_json(
        root / "transformations.json",
        {"schema_version": "1.1.0", "stages": []},
    )
    _write_json(
        root / "omissions.json",
        {
            "schema_version": "1.1.0",
            "omissions": [
                {
                    "scope": GENERIC_READ_BACKED_OMISSION,
                    "reason": "synthetic base has no aligned-read inputs",
                    "pipeline_stage": "read_backed_extensions",
                    "interpretation_consequence": "no claim",
                    "recovery": "run the extension",
                }
            ],
        },
    )
    _write_json(
        root / "relationships.json",
        {"schema_version": "1.1.0", "relationships": []},
    )
    _write_json(
        root / "field_glossary.json",
        {"schema_version": "1.1.0", "fields": {}},
    )
    metadata_paths = [
        "manifest.json",
        "schema.json",
        "variants.parquet",
        "extensions/base-evidence.bin",
        "field_glossary.json",
        "inventory.json",
        "raw_artifacts.json",
        "transformations.json",
        "omissions.json",
        "relationships.json",
        "readmygenome.md",
        "readmygenome/SKILL.md",
        "readmygenome/scripts/query.py",
        "readmygenome/reference/query_catalog.md",
        "readmygenome/reference/synonyms.json",
        "readmygenome/evals/canonical_questions.json",
    ]
    entries = [
        {
            "authority": "metadata",
            "bytes": None,
            "logical_role": "synthetic test input",
            "media_type": "application/octet-stream",
            "path": path,
            "required": path in {"manifest.json", "schema.json", "variants.parquet"},
            "rows": None,
            "sha256": None,
        }
        for path in metadata_paths
    ]
    _write_json(root / "inventory.json", {"schema_version": "1.1.0", "entries": entries})
    files = {
        path: {"bytes": 0, "rows": None, "sha256": "0" * 64}
        for path in (
            "extensions",
            "field_glossary.json",
            "inventory.json",
            "omissions.json",
            "raw_artifacts.json",
            "readmygenome",
            "readmygenome.md",
            "relationships.json",
            "transformations.json",
        )
    }
    files.update(
        {
            "schema.json": {"bytes": 20, "rows": None, "sha256": "0" * 64},
            "variants.parquet": {"bytes": 19, "rows": 1, "sha256": "0" * 64},
        }
    )
    _write_json(
        root / "manifest.json",
        {
            "schema_version": "1.1.0",
            "pipeline_version": "genome-annotation-builder-0.1.18",
            "source_bundle_sha256": "8" * 64,
            "genome_build": "GRCh38",
            "generated_at": "2026-01-01T00:00:00Z",
            "annotation_sources": {},
            "consensus_definitions": {},
            "normalization": "synthetic",
            "profiles": {"agent_exploration": "1.0.0"},
            "privacy": {"input_paths_embedded": False, "sample_identifier_embedded": False},
            "resource_provenance": {
                "reference_fasta_sha256": "1" * 64,
                "reference_fai_sha256": "2" * 64,
                "builder_implementation_sha256": "3" * 64,
                "genome_spec_commit": "7354ae288784234f854562b7dceb87f067fad550",
                "schema_sha256": "4" * 64,
            },
            "files": files,
        },
    )
    return root


def _stage(tmp_path: Path) -> Path:
    subject_binding, inputs, input_gates, receipts, adapters = complete_stage_inputs(tmp_path)
    stage = tmp_path / "release-c-stage"
    stage_release_c_extensions(stage, inputs, input_gates, receipts, adapters, subject_binding)
    return stage


def test_bundle_extension_is_atomic_complete_and_preserves_base(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base = _fake_base(tmp_path / "base.genome")
    stage = _stage(tmp_path)
    output = tmp_path / "enriched.genome"
    base_manifest_before = hashlib.sha256((base / "manifest.json").read_bytes()).hexdigest()
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

    report = build_enriched_bundle(
        base,
        stage,
        output,
        materialization="hardlink",
        generated_at="2026-02-03T04:05:06Z",
        validate_base=False,
    )

    assert report["status"] == "PASS"
    assert output.is_dir()
    assert hashlib.sha256((base / "manifest.json").read_bytes()).hexdigest() == (
        base_manifest_before
    )
    assert (base / "schema.json").stat().st_ino == (output / "schema.json").stat().st_ino
    assert (base / "manifest.json").stat().st_ino != (output / "manifest.json").stat().st_ino
    assert directory_digest(stage) == directory_digest(output / "extensions" / "release_c")

    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    release_c = manifest["resource_provenance"]["release_c"]
    assert manifest["pipeline_version"].endswith("+release-c.1")
    assert manifest["generated_at"] == "2026-02-03T04:05:06Z"
    assert release_c["source_rows_removed"] == 0
    assert release_c["caller_or_replicate_rows_collapsed"] == 0
    assert set(release_c["domains"]) == {
        "mitochondrial_heteroplasmy",
        "structural_variants",
    }
    assert manifest["files"]["extensions"]["sha256"] == directory_digest(output / "extensions")

    omissions = json.loads((output / "omissions.json").read_text(encoding="utf-8"))["omissions"]
    scopes = {row["scope"] for row in omissions}
    assert GENERIC_READ_BACKED_OMISSION not in scopes
    assert DOMAIN_OMISSIONS["structural_variants"] not in scopes
    assert DOMAIN_OMISSIONS["mitochondrial_heteroplasmy"] not in scopes
    assert DOMAIN_OMISSIONS["complex_pharmacogenes"] in scopes
    assert "genome-wide mosaic variant calling" in scopes
    assert PORTABLE_CORE_SCOPE_OMISSION in scopes
    scope_row = next(row for row in omissions if row["scope"] == PORTABLE_CORE_SCOPE_OMISSION)
    assert scope_row["source_records"] == 2
    assert scope_row["typed_rows"] == 1
    assert scope_row["omitted_records"] == 1

    raw = json.loads((output / "raw_artifacts.json").read_text(encoding="utf-8"))
    assert len(raw["artifacts"]) == 2
    assert all(
        row["artifact_type"].startswith("Release-C caller raw artifact:")
        for row in raw["artifacts"]
    )
    transformations = json.loads((output / "transformations.json").read_text(encoding="utf-8"))
    assert transformations["stages"][-1]["stage_id"] == "release_c_bundle_extension"
    assert transformations["stages"][-1]["caller_or_replicate_rows_collapsed"] == 0

    glossary = json.loads((output / "field_glossary.json").read_text(encoding="utf-8"))
    assert "structural_variants.evidence_id" in glossary["fields"]
    assert "mitochondrial_heteroplasmy.allele_fraction" in glossary["fields"]
    query_text = (output / "readmygenome" / "scripts" / "query.py").read_text(encoding="utf-8")
    assert 'elif operation == "read-backed"' in query_text

    inventory = json.loads((output / "inventory.json").read_text(encoding="utf-8"))
    inventory_paths = {row["path"] for row in inventory["entries"]}
    assert "extensions/release_c/release_c_extensions.json" in inventory_paths
    assert "extensions/release_c/structural_variants.parquet" in inventory_paths
    assert "extensions/release_c/mitochondrial_heteroplasmy.parquet" in inventory_paths

    validation = validate_enriched_bundle(output, run_read_backed_query=True)
    assert validation["status"] == "PASS"
    assert validation["release_c"]["rows"] == 2

    with pytest.raises(PortableEvidenceError, match="refusing to replace"):
        build_enriched_bundle(base, stage, output, validate_base=False)
