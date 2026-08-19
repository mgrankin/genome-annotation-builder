from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from genome_annotation_builder.release_c.readback_targets import (
    ReadbackTargetError,
    create_readback_targets,
    validate_readback_targets,
)
from genome_annotation_builder.util import sha256_file


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def synthetic_bundle(tmp_path: Path, *, clinical_ids: list[str] | None = None) -> Path:
    bundle = tmp_path / "PRIVATE_SUBJECT.genome"
    family = bundle / "extensions/family_evidence.parquet/chrom=chr1"
    family.mkdir(parents=True)
    clinical = bundle / "clinical_findings.parquet"
    ids = clinical_ids or ["v-clinical"]
    pq.write_table(pa.table({"variant_id": ids}), clinical)
    rows = [
        {
            "variant_id": "v-clinical",
            "chrom": "chr1",
            "pos": 10,
            "ref": "A",
            "alt": "C",
            "de_novo_candidate": False,
            "parental_mosaicism_candidate": False,
            "parental_dropout_suspected": False,
            "child_dropout_suspected": False,
            "mendelian_status": "consistent",
            "loco_status": "not_individually_persisted",
            "caller_concordance": "exact_genotype_concordant",
            "technical_replicate_concordance": "not_assessed",
            "complex_region": False,
            "difficult_region_flags": [],
            "phase_confidence": "family_phased",
        },
        {
            "variant_id": "v-denovo",
            "chrom": "chr1",
            "pos": 20,
            "ref": "A",
            "alt": "G",
            "de_novo_candidate": False,
            "parental_mosaicism_candidate": False,
            "parental_dropout_suspected": False,
            "child_dropout_suspected": False,
            "mendelian_status": "conflict",
            "loco_status": "de_novo_candidate",
            "caller_concordance": "exact_genotype_concordant",
            "technical_replicate_concordance": "not_assessed",
            "complex_region": False,
            "difficult_region_flags": [],
            "phase_confidence": "unphased",
        },
        {
            "variant_id": "v-loco",
            "chrom": "chr1",
            "pos": 30,
            "ref": "A",
            "alt": "T",
            "de_novo_candidate": False,
            "parental_mosaicism_candidate": False,
            "parental_dropout_suspected": False,
            "child_dropout_suspected": False,
            "mendelian_status": "consistent",
            "loco_status": "source_genotype_discordant",
            "caller_concordance": "exact_genotype_concordant",
            "technical_replicate_concordance": "not_assessed",
            "complex_region": True,
            "difficult_region_flags": ["low_mappability"],
            "phase_confidence": "family_phased",
        },
        {
            "variant_id": "v-mendel",
            "chrom": "chr1",
            "pos": 40,
            "ref": "A",
            "alt": "AC",
            "de_novo_candidate": False,
            "parental_mosaicism_candidate": False,
            "parental_dropout_suspected": False,
            "child_dropout_suspected": False,
            "mendelian_status": "conflict",
            "loco_status": "not_individually_persisted",
            "caller_concordance": "exact_genotype_concordant",
            "technical_replicate_concordance": "not_assessed",
            "complex_region": False,
            "difficult_region_flags": [],
            "phase_confidence": "unphased",
        },
        {
            "variant_id": "v-caller",
            "chrom": "chr1",
            "pos": 50,
            "ref": "A",
            "alt": "C",
            "de_novo_candidate": False,
            "parental_mosaicism_candidate": False,
            "parental_dropout_suspected": False,
            "child_dropout_suspected": False,
            "mendelian_status": "consistent",
            "loco_status": "not_individually_persisted",
            "caller_concordance": "discordant",
            "technical_replicate_concordance": "not_assessed",
            "complex_region": False,
            "difficult_region_flags": [],
            "phase_confidence": "unphased",
        },
        {
            "variant_id": "v-routine",
            "chrom": "chr1",
            "pos": 60,
            "ref": "A",
            "alt": "G",
            "de_novo_candidate": False,
            "parental_mosaicism_candidate": False,
            "parental_dropout_suspected": False,
            "child_dropout_suspected": False,
            "mendelian_status": "consistent",
            "loco_status": "not_individually_persisted",
            "caller_concordance": "exact_genotype_concordant",
            "technical_replicate_concordance": "not_assessed",
            "complex_region": False,
            "difficult_region_flags": [],
            "phase_confidence": "unphased",
        },
    ]
    pq.write_table(pa.Table.from_pylist(rows), family / "part-00000.parquet")
    write_json(
        bundle / "manifest.json",
        {
            "schema_version": "1.1.0",
            "genome_build": "GRCh38",
            "pipeline_version": "synthetic-test",
            "files": {
                "clinical_findings.parquet": {
                    "bytes": clinical.stat().st_size,
                    "rows": len(ids),
                    "sha256": sha256_file(clinical),
                },
                "variants.parquet": {"bytes": 1, "rows": len(rows), "sha256": "1" * 64},
                "extensions": {"bytes": 1, "rows": None, "sha256": "2" * 64},
            },
        },
    )
    return bundle


def test_target_selection_is_bounded_deterministic_and_revalidatable(tmp_path: Path) -> None:
    bundle = synthetic_bundle(tmp_path)
    first_targets = tmp_path / "first.tsv"
    first_receipt = tmp_path / "first.receipt.json"
    first = create_readback_targets(bundle, first_targets, first_receipt, maximum_targets=3)
    assert first["candidate_rows"] == 5
    assert first["selected_rows"] == 3
    assert first["truncated_rows"] == 2
    assert first["required_clinical_rows"] == 1
    assert "PRIVATE_SUBJECT" not in first_receipt.read_text(encoding="utf-8")
    validated = validate_readback_targets(first_targets, first_receipt, bundle_path=bundle)
    assert validated["status"] == "PASS"
    assert validated["targets"] == 3
    assert validated["source_bundle_revalidated"] is True

    second_targets = tmp_path / "second.tsv"
    second_receipt = tmp_path / "second.receipt.json"
    create_readback_targets(bundle, second_targets, second_receipt, maximum_targets=3)
    assert first_targets.read_bytes() == second_targets.read_bytes()
    assert first_receipt.read_bytes() == second_receipt.read_bytes()


def test_target_validation_rejects_changed_tsv(tmp_path: Path) -> None:
    bundle = synthetic_bundle(tmp_path)
    targets = tmp_path / "targets.tsv"
    receipt = tmp_path / "targets.receipt.json"
    create_readback_targets(bundle, targets, receipt, maximum_targets=3)
    targets.write_bytes(targets.read_bytes() + b"changed\n")
    with pytest.raises(ReadbackTargetError, match="target_manifest_identity_changed"):
        validate_readback_targets(targets, receipt)


def test_target_selection_never_truncates_required_clinical_rows(tmp_path: Path) -> None:
    bundle = synthetic_bundle(tmp_path, clinical_ids=["v-clinical", "v-denovo"])
    with pytest.raises(
        ReadbackTargetError, match="maximum_targets_cannot_retain_all_clinical_findings"
    ):
        create_readback_targets(
            bundle,
            tmp_path / "targets.tsv",
            tmp_path / "targets.receipt.json",
            maximum_targets=1,
        )
