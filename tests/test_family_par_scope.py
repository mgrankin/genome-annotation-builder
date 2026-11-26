from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import duckdb
import pysam
import pytest
from test_pgx import _rows
from test_vcf_core import _environment

from genome_annotation_builder.errors import InputValidationError
from genome_annotation_builder.par_scope import DECLARED, GRCH38, MASKED, TAG, source_par_scope
from genome_annotation_builder.pgs import _resolve_locus
from genome_annotation_builder.pgx_materialize import _resolve_pgx_locus
from genome_annotation_builder.pipeline import build_bundle
from genome_annotation_builder.query import _register_views, par_query_limit, query_bundle
from genome_annotation_builder.resources import sync_pinned_resources
from genome_annotation_builder.validator import _validate_target_callability


def test_masked_profile_trust_and_legacy_scope() -> None:
    h = pysam.VariantHeader()
    for chrom, (length, _) in GRCH38.items():
        h.contigs.add(chrom, length=length)
    assert source_par_scope(h) == []
    h.add_line('##family_reference_confidence_overlay=<Version="2">')
    legacy = source_par_scope(h)
    assert len(legacy) == 4 and all(r["scope"] == "negative" for r in legacy)
    for c, (_, intervals) in GRCH38.items():
        for i, (a, b) in enumerate(intervals, 1):
            h.add_line(
                f"##{TAG}=<ID={c}_{i},Contig={c},Start0={a},End0={b},Profile={MASKED},"
                f"CanonicalSHA256={'a' * 64},MappingSHA256={'b' * 64},PARBedSHA256={'c' * 64}>"
            )
            if c == "X" and i == 1:
                with pytest.raises(InputValidationError, match="incomplete"):
                    source_par_scope(h)
    scopes = source_par_scope(h)
    with pytest.raises(InputValidationError, match="pinned FASTA"):
        source_par_scope(h, reference_sha256="d" * 64)
    for c, (_, intervals) in GRCH38.items():
        for lo, hi in intervals:
            assert par_query_limit(scopes, "chr" + c, lo, lo) is None
            limit = par_query_limit(scopes, "chr" + c, lo, lo + 1)
            assert limit["scope"] == ("all" if c == "Y" else "negative")
            assert par_query_limit(scopes, "chr" + c, hi, hi + 1) == limit
            assert par_query_limit(scopes, "chr" + c, hi + 1, hi + 1) is None


def test_scoped_reference_blocks_preserve_observation_and_ref_span(tmp_path: Path) -> None:
    cfg, _ = _environment(tmp_path, target_catalog=True)
    # A producer-declared interval in a short synthetic chromosome, spanning a
    # single source block whose left and right sides remain valid controls.
    ref = tmp_path / "scoped.fa"
    ref.write_text(">X\n" + "ACGT" * 50 + "\n")
    pysam.faidx(str(ref))
    vcf = tmp_path / "scoped.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n##contig=<ID=X,length=200>\n"
        f"##{TAG}=<ID=X_1,Contig=X,Start0=11,End0=20,Profile={DECLARED}>\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID\n"
        "X\t9\t.\tA\t<*>\t.\tPASS\tEND=22\tGT:DP:GQ\t0/0:35:60\n"
        "X\t16\t.\tT\tC\t.\tPASS\t.\tGT:DP:GQ\t0/1:35:60\n"
        "X\t17\t.\tA\tC\t.\tPASS\t.\tGT:DP:GQ\t1/.:35:60\n"
    )
    catalog = tmp_path / "targets.tsv"
    catalog.write_text(
        "target_id\tchrom\tpos\tref\talt\tgene_symbol\ttarget_source\t"
        "target_source_version\tevidence_scope\n"
        + "".join(
            f"{k}\tX\t{pos}\t{ref}\t{alt}\t\tcontrol\t1\tsite\n"
            for k, pos, ref, alt in [
                ("before", 10, "C", "T"),
                ("cross", 11, "GT", "G"),
                ("inside", 12, "T", "G"),
                ("observed", 16, "T", "C"),
                ("after", 21, "A", "C"),
            ]
        )
    )
    cfg = replace(
        cfg,
        input_vcf=Path(pysam.tabix_index(str(vcf), preset="vcf")),
        reference_fasta=ref,
        reference_expected_sha256=hashlib.sha256(ref.read_bytes()).hexdigest(),
        target_catalog_path=catalog,
        target_catalog_expected_sha256=hashlib.sha256(catalog.read_bytes()).hexdigest(),
    )
    sync_pinned_resources(cfg)
    assert build_bundle(cfg)["validation"]["valid"]
    for term, status in [
        ("chrX:10:C:T", "reference_observed"),
        ("chrX:11", "reference_observed"),
        ("chrX:11:GT:G", "unknown"),
        ("chrX:12:T:G", "unknown"),
        ("chrX:15:G:T", "unknown"),
        ("chrX:16:T:C", "non_reference_observed"),
        ("chrX:16:T:A", "overlapping_variant_observed"),
        ("chrX:17:A:C", "partial_non_reference_observed"),
        ("chrX:21:A:C", "reference_observed"),
        ("chrX:22:CG:C", "unknown"),
    ]:
        result = query_bundle(cfg.bundle_path, "callability", term=term)
        assert result["conclusion"] == status
    result = query_bundle(cfg.bundle_path, "callability", term="chrX:11:GT:G")
    raw = next(
        r for r in result["rows"] if r["_table"] == "extensions/reference_confidence_regions"
    )
    assert (raw["gt"], raw["dp"], raw["gq"], raw["callable"]) == ([0, 0], 35, 60, True)
    target = next(
        r for r in result["rows"] if r["_table"] == "extensions/target_callability_evidence"
    )
    assert target["selected_reference_evidence_ids"] == [raw["evidence_id"]]
    assert (target["resolution"], target["reference_observed"], target["callable"]) == (
        "mapping_scope_unassessed",
        None,
        False,
    )
    for operation in ("variant", "evidence", "search"):
        proof = query_bundle(cfg.bundle_path, operation, term="chrX:16:T:C")
        assert proof["rows"] and all(
            r["family_mapping_scope"]["scope"] == "negative" for r in proof["rows"]
        )
    # Audit a lost scope declaration independently of the materializer.
    with duckdb.connect() as connection:
        views = _register_views(connection, cfg.bundle_path)
        errors = []
        manifest = json.loads((cfg.bundle_path / "manifest.json").read_bytes())
        metrics = _validate_target_callability(connection, views, manifest, errors, [])
        assert metrics["target_mapping_scope_mismatches"] == 2
        assert any("full target REF span" in message for message in errors)
    script = cfg.bundle_path / "readmygenome/scripts/query.py"
    run = subprocess.run(
        [
            sys.executable,
            "-I",
            str(script),
            "--bundle",
            str(cfg.bundle_path),
            "callability",
            "chrX:14:C:T",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(run.stdout)["conclusion"] == "unknown"


def test_scoped_target_never_enters_matcher_or_score() -> None:
    rows = _rows()
    for r in rows:
        r.update(resolution="mapping_scope_unassessed", callable=False, reference_observed=None)
    rows[0]["variant_gt"] = [
        1,
        1,
    ]  # A source-reported GT can coexist with withheld Y-PAR resolution.
    pgx = _resolve_pgx_locus(rows)
    assert pgx["status"] == "no_call_unknown" and not pgx["eligible_for_matcher"]
    pgs = _resolve_locus(
        {
            **rows[0],
            "pgs_id": "PGS1",
            "scoring_row_number": 1,
            "scoring_row_digest": "a" * 64,
            "effect_allele": "A",
            "source_effect_allele": "A",
            "effect_weight": 1.0,
            "supported": True,
        }
    )
    assert (pgs["status"], pgs["dosage"], pgs["contribution"]) == ("missing_unknown", None, None)
