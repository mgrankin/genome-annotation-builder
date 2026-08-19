from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from genome_annotation_builder.query import _digest, query_bundle
from genome_annotation_builder.release_c.portable_model import (
    make_evidence_id,
    schema_for,
    validate_rows,
)


def _required_default(field: pa.Field):
    if field.nullable:
        return None
    if pa.types.is_string(field.type):
        return ""
    if pa.types.is_integer(field.type):
        return 0
    if pa.types.is_floating(field.type):
        return 0.0
    if pa.types.is_boolean(field.type):
        return False
    if pa.types.is_list(field.type):
        return []
    raise AssertionError(f"no synthetic default for {field}")


def _readback_row() -> dict[str, object]:
    domain = "read_backed_evidence"
    schema = schema_for(domain)
    row = {field.name: _required_default(field) for field in schema}
    row.update(
        {
            "target_id": "synthetic-target-1",
            "target_manifest_sha256": "1" * 64,
            "chrom": "chr1",
            "pos": 101,
            "ref": "A",
            "alt": "T",
            "variant_type": "SNV",
            "assessment_status": "assessed",
            "eligible_fragments": 2,
            "ref_reads": 1,
            "alt_reads": 1,
            "other_reads": 0,
            "forward_ref_reads": 1,
            "reverse_ref_reads": 0,
            "forward_alt_reads": 0,
            "reverse_alt_reads": 1,
            "f1r2_ref_reads": 1,
            "f2r1_ref_reads": 0,
            "unclassified_orientation_ref_reads": 0,
            "f1r2_alt_reads": 0,
            "f2r1_alt_reads": 1,
            "unclassified_orientation_alt_reads": 0,
            "duplicate_reads_included": False,
            "overlapping_mate_policy": "collapse_by_query_name",
            "local_complexity_status": "standard",
            "local_haplotype_status": "not_assessed",
            "candidate_category": "MOSAIC_CANDIDATE",
            "event_id": "synthetic-event-1",
            "assembly": "GRCh38",
            "reference_sha256": "2" * 64,
            "technical_source_role": "primary",
            "caller_name": "GAB bounded local pileup",
            "caller_version": "1.0.0",
            "caller_run_identity_sha256": "3" * 64,
            "resource_manifest_sha256": "4" * 64,
            "source_alignment_sha256": "5" * 64,
            "raw_artifact_id": "bounded-readback-json",
            "raw_artifact_sha256": "6" * 64,
            "source_record_digest": "7" * 64,
            "filter_reasons": [],
            "confidence_tier": "SUPPORTED",
            "limitations": ["synthetic query fixture; not a clinical assertion"],
            "confirmation_required": True,
            "concordance_group_id": None,
            "concordance_status": "not_assessed",
        }
    )
    row["evidence_id"] = make_evidence_id(domain, row)
    return validate_rows(domain, [row], expected_reference_sha256="2" * 64)[0]


def _minimal_bundle(root: Path, *, include_release_c: bool = True) -> Path:
    variants = root / "variants.parquet"
    variants.mkdir(parents=True)
    pq.write_table(pa.table({"placeholder": [1]}), variants / "part-0.parquet")
    release_c = root / "extensions" / "release_c"
    release_c.mkdir(parents=True)
    if include_release_c:
        pq.write_table(
            pa.Table.from_pylist([_readback_row()], schema=schema_for("read_backed_evidence")),
            release_c / "read_backed_evidence.parquet",
        )
    schema_path = root / "schema.json"
    schema_path.write_text('{"version":"1.1.0"}\n', encoding="utf-8")
    manifest = {
        "schema_version": "1.1.0",
        "generated_at": "2026-01-01T00:00:00Z",
        "annotation_sources": {},
        "consensus_definitions": {},
        "files": {
            "schema.json": {
                "sha256": _digest(schema_path),
                "bytes": schema_path.stat().st_size,
                "rows": None,
            },
            "variants.parquet": {
                "sha256": _digest(variants),
                "bytes": None,
                "rows": 1,
            },
            "extensions": {
                "sha256": _digest(root / "extensions"),
                "bytes": None,
                "rows": None,
            },
        },
    }
    (root / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return root


def test_read_backed_query_is_bounded_and_source_specific(tmp_path: Path) -> None:
    bundle = _minimal_bundle(tmp_path / "synthetic.genome")

    by_locus = query_bundle(bundle, "read-backed", term="chr1:101:A:T", limit=10)
    assert by_locus["count"] == 1
    assert by_locus["rows"][0]["target_id"] == "synthetic-target-1"
    assert by_locus["rows"][0]["_table"] == ("extensions/release_c/read_backed_evidence")
    assert "do not replace, average, or vote away" in by_locus["interpretation"]

    by_category = query_bundle(bundle, "read-backed", term="mosaic", limit=10)
    assert by_category["count"] == 1
    assert by_category["rows"][0]["candidate_category"] == "MOSAIC_CANDIDATE"


def test_read_backed_query_reports_absent_extension(tmp_path: Path) -> None:
    bundle = _minimal_bundle(tmp_path / "synthetic.genome", include_release_c=False)

    result = query_bundle(bundle, "read-backed", term="chr1:101", limit=10)
    assert result["count"] == 0
    assert result["unavailable"] == "read-backed Release-C extensions were not built"
