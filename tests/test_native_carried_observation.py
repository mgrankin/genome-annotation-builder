import gzip
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pysam
import pytest
from test_vcf_core import _environment

from genome_annotation_builder import native_observation
from genome_annotation_builder.errors import BundleValidationError, InputValidationError
from genome_annotation_builder.pgx_materialize import _resolve_pgx_locus
from genome_annotation_builder.pipeline import build_bundle
from genome_annotation_builder.query import query_bundle
from genome_annotation_builder.resources import sync_pinned_resources
from genome_annotation_builder.validator import validate_bundle


def native_header():
    return (
        "##native_projection_schema=family-native-projection/2.2\n"
        "##family_ploidy_normalization=<ID=Confirmed_XX_XY_MT>\n"
        '##INFO=<ID=NID,Number=1,Type=String,Description="Native record hash">\n'
        '##INFO=<ID=NAI,Number=A,Type=Integer,Description="Original ALT index">\n'
        + "".join(
            f'##FORMAT=<ID={t},Number=1,Type=String,Description="{t}">\n'
            for t in ["NGT", "NAD", "NPL", "NDP", "NGQ", "NFT", "BGT", "NSTATE"]
        )
        + "".join(
            f'##FORMAT=<ID={t},Number=1,Type=Integer,Description="{t}">\n'
            for t in ["NDOS", "NOTH", "NUNK"]
        )
    )


def row(pos, ref, alt, gt=".", original="2/.", before="1/.", counts="1:0:1:CARRIED_PARTIAL_NATIVE"):
    digest = hashlib.sha256(f"test-native-record-{pos}".encode()).hexdigest()
    return (
        f"X\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\tNID={digest};NAI=2\t"
        f"GT:NGT:BGT:NDOS:NOTH:NUNK:NSTATE\t{gt}:{original}:{before}:{counts}\n"
    )


def fixture(tmp_path, *, registered=True):
    config, source = _environment(tmp_path, pgs=True, pgs_calibration=False)
    (model,) = config.pgs_models
    # Move one existing scoring row onto the tested allele.
    score = model.scoring_file.read_text().replace(
        "rs100001\tT\tC\t1.0\tGENE1\tENSEMBL\trs100001\t1\t10\t",
        "rs100001\tC\tA\t1.0\tGENE1\tENSEMBL\trs100001\tX\t7\t",
    )
    model.scoring_file.write_text(score)
    config = replace(
        config,
        pgs_models=(replace(model, expected_sha256=hashlib.sha256(score.encode()).hexdigest()),),
    )
    with config.reference_fasta.open("a") as fh:
        fh.write(">X\n" + "A" * 100 + "\n")
    pysam.faidx(str(config.reference_fasta))
    with gzip.open(source, "rt") as fh:
        text = fh.read()
    headers = native_header()
    if not registered:
        headers = headers.replace("##native_projection_schema=family-native-projection/2.2\n", "")
    text = text.replace("#CHROM", "##contig=<ID=X,length=100>\n" + headers + "#CHROM")
    text += row(7, "A", "C") + row(12, "AA", "A") + row(15, "A", "*")
    text += row(30, "A", "T", original="3/.", before="./.", counts="0:1:1:OTHER_ALT_AND_UNKNOWN")
    path = tmp_path / "native.vcf"
    path.write_text(text)
    return replace(config, input_vcf=Path(pysam.tabix_index(str(path), preset="vcf")))


@pytest.mark.parametrize("registered", [True, False])
def test_native_observation_in_portable_build(tmp_path: Path, registered: bool):
    config = fixture(tmp_path, registered=registered)
    sync_pinned_resources(config)
    assert build_bundle(config)["validation"]["valid"]
    bundle = config.bundle_path
    exact = query_bundle(bundle, "variant", term="chrX:7:A:C")["rows"]
    count = json.loads((bundle / "extensions/reconciliation.json").read_bytes())[
        "native_carried_no_call_alleles"
    ]["count"]
    if not registered:
        assert exact == [] and count == 0
        assert (
            query_bundle(bundle, "callability", term="chrX:7:A:C")["conclusion"]
            == "unknown_no_call"
        )
        return
    assert count == 3
    audit = pq.read_table(bundle / "extensions/pgs_locus_audit.parquet").to_pylist()
    (score,) = [r for r in audit if r["variant_id"] == "chrX:7:A:C"]
    assert (score["resolution"], score["dosage"], score["contribution"]) == ("no_call", None, None)
    targets = pq.read_table(bundle / "extensions/target_callability_evidence.parquet").to_pylist()
    (target,) = [r for r in targets if r["variant_id"] == "chrX:7:A:C"]
    assert target["callable"] is False and target["observed_variant_ids"] == []
    # The pharmacogenomic adapter uses this same target/complete-variant join.
    pgx = _resolve_pgx_locus(
        [
            {
                **target,
                "locus_id": "test",
                "record_number": 1,
                "allele_index": 1,
                "target_variant_id": target["variant_id"],
                "normalization_changed": False,
                "structurally_complex": False,
                "source_record_digest": "a" * 64,
                "positions_vcf_sha256": "b" * 64,
            }
        ]
    )
    assert pgx["status"] == "no_call_source_missing" and pgx["eligible_for_matcher"] is False
    assert [
        (r["_table"], r["gt"], r["phased"], r["source_gt"], r["source_alt_index"]) for r in exact
    ] == [("extensions/observed_partial_variants", [None], False, ".", 1)]
    assert (
        query_bundle(bundle, "callability", term="chrX:7:A:C")["conclusion"]
        == "partial_non_reference_observed"
    )
    assert query_bundle(bundle, "callability", term="chrX:7:A:G")["conclusion"] == "unknown_no_call"
    assert query_bundle(bundle, "variant", term="chrX:30:A:T")["count"] == 0
    (evidence,) = query_bundle(bundle, "evidence", term="chrX:1:AA:A")["rows"]
    assert (
        evidence["source_pos"] == 12
        and evidence["source_evidence_reason"] == native_observation.REASON
    )
    assert json.loads(evidence["source_info_json"])["NAI"] == [2]
    assert json.loads(evidence["source_other_format_json"])["NGT"] == "2/."
    exclusions = pq.read_table(bundle / "extensions/non_core_variants.parquet").to_pylist()
    (star,) = [r for r in exclusions if r["alt"] == "*"]
    assert star["source_gt"] == [None] and star["category"] == "symbolic_or_non_acgt_allele"
    ref_path = next((bundle / "extensions/reference_confidence_regions.parquet").rglob("*.parquet"))
    table = pq.ParquetFile(ref_path).read()
    changed = table.to_pylist()
    parent = next(r for r in changed if r["reason"] == native_observation.REASON)
    fmt = json.loads(parent["other_format_json"])
    parent["other_format_json"] = json.dumps({**fmt, "NGT": "3/."})
    pq.write_table(pa.Table.from_pylist(changed, schema=table.schema), ref_path)
    with pytest.raises(BundleValidationError, match="invalid native carried ALT"):
        validate_bundle(bundle, reference_fasta=config.reference_fasta)


@pytest.mark.parametrize(
    ("info", "sample", "gt", "chrom", "pos", "end", "sex"),
    [
        ({"NID": "wrong"}, {}, [None], "chrX", 9, 9, 1),
        ({"NAI": [1]}, {}, [None], "chrX", 9, 9, 1),
        ({}, {"NSTATE": "REFERENCE_GT_UNQUALIFIED"}, [None], "chrX", 9, 9, 1),
        ({}, {"NUNK": 0}, [None], "chrX", 9, 9, 1),
        ({}, {"BGT": "./1"}, [None], "chrX", 9, 9, 1),
        ({}, {}, [None, None], "chrX", 9, 9, 1),
        ({}, {}, [0], "chrX", 9, 9, 1),
        ({}, {}, [None], "chrX", 9, 9, 2),
        ({}, {}, [None], "chrY", 9, 9, 2),
        ({}, {}, [None], "chrX", 9, 11, 1),
        ({}, {}, [None], "chr1", 9, 9, 1),
        ({}, {}, [None], "unplaced", 9, 9, 1),
    ],
)
def test_contradictory_native_evidence_is_rejected(info, sample, gt, chrom, pos, end, sex):
    with pytest.raises(InputValidationError):
        native_observation.validate_native_no_call(
            {"NID": "a" * 64, "NAI": [2], **info},
            {
                "NGT": "2/.",
                "BGT": "1/.",
                "NDOS": 1,
                "NOTH": 0,
                "NUNK": 1,
                "NSTATE": "CARRIED_PARTIAL_NATIVE",
                **sample,
            },
            gt,
            ["C"],
            chrom,
            pos,
            end,
            [{"chrom": "chrX", "start0": 10, "end0": 20}],
            sex=sex,
        )
