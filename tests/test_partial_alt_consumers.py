from __future__ import annotations

import gzip
import json
from dataclasses import replace
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
import pytest
from test_vcf_core import _environment

from genome_annotation_builder.errors import BundleValidationError, GenomeBuilderError
from genome_annotation_builder.pipeline import build_bundle
from genome_annotation_builder.query import query_bundle
from genome_annotation_builder.resources import sync_pinned_resources
from genome_annotation_builder.validator import validate_bundle
from genome_annotation_builder.vcf_ingest import _canonical_json_text


def test_vcf_provenance_is_one_json_object() -> None:
    text = _canonical_json_text(
        {
            "slots": (1, None),
            "nested": {"flag": True, "text": b"source", "depth": 3.5},
            "nonfinite": [float("nan"), float("inf"), -float("inf")],
        }
    )
    assert json.loads(text) == {
        "slots": [1, None],
        "nested": {"flag": True, "text": "source", "depth": 3.5},
        "nonfinite": ["nan", "inf", "-inf"],
    }
    assert _canonical_json_text({}) is None


def test_partial_alleles_are_searchable_without_a_complete_genotype_or_score(
    tmp_path: Path,
) -> None:
    config, original = _environment(tmp_path, target_catalog=True, pgs=True)
    with gzip.open(original, "rt") as handle:
        text = handle.read()
    text = text.replace(
        "#CHROM",
        '##INFO=<ID=PROBE,Number=1,Type=String,Description="Source label">\n'
        '##FORMAT=<ID=EXTRA,Number=.,Type=Integer,Description="Other values">\n#CHROM',
    ).replace(
        "T,<*>\t60\tPASS\t.\tGT:DP:GQ:AD:FT\t0/1:30:55:15,15,0:PASS",
        "T,<*>\t60\tPASS\tPROBE=test\tGT:DP:GQ:AD:FT:EXTRA\t1/.:30:55:15,15,0:PASS:7,.,9",
    )
    # The repeat deletion normalizes to chr1:1, outside its source no-call interval.
    text += "1\t58\trs900002\tCGTAC\tC\t60\tPASS\t.\tGT:DP:GQ:AD:FT\t.|1:30:55:15,15:PASS\n"
    path = tmp_path / "partial.vcf"
    path.write_text(text)
    config = replace(config, input_vcf=Path(pysam.tabix_index(str(path), preset="vcf")))
    sync_pinned_resources(config)
    result = build_bundle(config)
    assert result["validation"]["valid"]
    bundle = config.bundle_path

    for operation, term in [
        ("variant", "rs100001"),
        ("search", "rs100001"),
        ("callability", "rs100001"),
        ("evidence", "chr1:10:C:T"),
    ]:
        result = query_bundle(bundle, operation, term=term)
        partial = [
            r for r in result["rows"] if r["_table"] == "extensions/observed_partial_variants"
        ]
        assert len(partial) == 1
        assert partial[0]["gt"] == [1, None]
        assert partial[0]["phased"] is False
        if operation == "callability":
            assert result["conclusion"] == "partial_non_reference_observed"
        if operation == "evidence":
            assert json.loads(partial[0]["source_other_format_json"])["EXTRA"] == [7, None, 9]
            assert json.loads(partial[0]["source_info_json"])["PROBE"] == "test"
    overlap = query_bundle(bundle, "callability", term="chr1:3:G:A")
    assert overlap["conclusion"] == "unknown_no_call"
    assert any(row.get("rsid") == "rs900002" for row in overlap["rows"])
    deletion = query_bundle(bundle, "variant", term="rs900002")["rows"]
    assert len(deletion) == 1
    assert (deletion[0]["pos"], deletion[0]["gt"], deletion[0]["phased"]) == (1, [None, 1], True)

    with duckdb.connect() as db:
        extensions = bundle / "extensions"
        no_call = db.execute(
            "SELECT source_info_json, other_format_json FROM read_parquet(?) "
            "WHERE start_pos=10 AND state='no_call'",
            [str(extensions / "reference_confidence_regions.parquet" / "**/*.parquet")],
        ).fetchone()
        assert json.loads(no_call[0])["PROBE"] == "test"
        assert json.loads(no_call[1])["EXTRA"] == [7, None, 9]
        scores = db.execute(
            "SELECT dosage, contribution FROM read_parquet(?) WHERE scoring_row_number=1",
            [str(extensions / "pgs_locus_audit.parquet")],
        ).fetchall()
        assert scores == [(None, None)]
        target = db.execute(
            "SELECT resolution, callable, reference_observed, overlapping_no_call_evidence_ids "
            "FROM read_parquet(?) WHERE target_id='explicit_reference_block'",
            [str(extensions / "target_callability_evidence.parquet")],
        ).fetchone()
        assert target[:3] == ("no_call", False, None)
        assert target[3] == [deletion[0]["source_id"]]
        assert db.execute(
            "SELECT count(*) FROM read_parquet(?) WHERE rsid = ?",
            [str(bundle / "variants.parquet" / "**/*.parquet"), "rs100001"],
        ).fetchone() == (0,)


def test_partial_gt_accounts_for_excluded_alleles_in_a_hybrid_build(tmp_path: Path) -> None:
    config, original = _environment(tmp_path, target_catalog=True)
    with config.reference_fasta.open("a") as handle:
        handle.write(">unplaced\nACGT\n")
    pysam.faidx(str(config.reference_fasta))
    with gzip.open(original, "rt") as handle:
        source = handle.read()
    source = source.replace("#CHROM", "##contig=<ID=unplaced,length=4>\n#CHROM").replace(
        "0/1:22:35:11,11:PASS", ".|1:22:35:11,11:PASS"
    )
    source += (
        "unplaced\t1\t.\tA\tC,G\t60\tPASS\t.\tGT\t./2\n"
        "unplaced\t2\t.\tC\tT\t60\tPASS\t.\tGT\t0/.\n"
        "unplaced\t3\t.\tG\tA\t60\tPASS\t.\tGT\t.\n"
    )
    plain = tmp_path / "excluded.vcf"
    plain.write_text(source)
    config = replace(config, input_vcf=Path(pysam.tabix_index(str(plain), preset="vcf")))
    sync_pinned_resources(config)
    assert build_bundle(config)["validation"]["valid"]
    extensions = config.bundle_path / "extensions"
    with duckdb.connect() as db:
        exclusions = db.execute(
            "SELECT chrom, alt, category, source_alt_index, source_gt, phased "
            "FROM read_parquet(?) ORDER BY chrom",
            [str(extensions / "non_core_variants.parquet/**/*.parquet")],
        ).fetchall()
        assert exclusions == [
            ("chr1", "<DEL>", "symbolic_or_non_acgt_allele", 1, [None, 1], True),
            ("unplaced", "G", "non_core_contig", 2, [None, 2], False),
        ]
        target = db.execute(
            "SELECT resolution, callable, reference_observed, observed_ploidy, "
            "len(observed_non_core_evidence_ids), len(overlapping_no_call_evidence_ids) "
            "FROM read_parquet(?) WHERE target_id='symbolic_overlap'",
            [str(extensions / "target_callability_evidence.parquet")],
        ).fetchone()
        assert target == ("no_call", False, None, 2, 1, 1)

    parquet = next((extensions / "non_core_variants.parquet").glob("**/*.parquet"))
    rows = pq.read_table(parquet)
    for changed_gt in ([0, 2], [None, 1]):
        changed = [
            dict(row, source_gt=changed_gt) if row["chrom"] == "unplaced" else row
            for row in rows.to_pylist()
        ]
        pq.write_table(pa.Table.from_pylist(changed, schema=rows.schema), parquet)
        with pytest.raises(
            BundleValidationError, match="excluded alleles differ from source genotype"
        ):
            validate_bundle(config.bundle_path, reference_fasta=config.reference_fasta)
    pq.write_table(
        rows.filter([value != "unplaced" for value in rows["chrom"].to_pylist()]), parquet
    )
    with pytest.raises(BundleValidationError, match="differ from explicit input GT alleles"):
        validate_bundle(config.bundle_path, reference_fasta=config.reference_fasta)


def test_existing_bundle_requires_matching_build_identity(tmp_path: Path) -> None:
    config, _ = _environment(tmp_path)
    sync_pinned_resources(config)
    built = build_bundle(config)
    assert build_bundle(config)["run_key"] == built["run_key"]
    changed = replace(config, quality=replace(config.quality, high_min_dp=99))
    with pytest.raises(GenomeBuilderError, match="different or unknown build"):
        build_bundle(changed)
    manifest_path = config.bundle_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    del manifest["resource_provenance"]["build_run_key"]
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(GenomeBuilderError, match="different or unknown build"):
        build_bundle(config)
