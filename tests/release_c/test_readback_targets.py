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


def legacy_bundle(tmp_path: Path, *, version: str = "0.2.7", fault: str = "") -> Path:
    bundle = synthetic_bundle(tmp_path)
    manifest_file = bundle / "manifest.json"
    manifest = json.loads(manifest_file.read_bytes())
    findings = "clinical_findings.parquet"
    (bundle / findings).unlink()
    if fault != "registered_missing":
        del manifest["files"][findings]
    digest = "a" * 64
    policy = {
        "version": "3.3",
        "snapshot": "20260815",
        "canonical_sha256": digest,
        "decision_rows": 0,
        "actionable_candidate_variants": 0,
        "decision_status_counts": {},
    }
    stage = {
        "stage_id": "acmg_secondary_findings_policy_evaluation",
        "software": "genome-annotation-builder",
        "software_version": version.partition("+")[0],
        "status": "completed",
        "decision_rows": 0,
        "actionable_candidate_variants": 0,
        "inputs": ["sha256:" + digest],
    }
    manifest["pipeline_version"] = "genome-annotation-builder-" + version
    manifest["resource_provenance"] = {
        "acmg_secondary_findings": policy,
        "ensembl_vep": {"cache_tree_sha256": "b" * 64},
        "clinvar": {"normalized_vcf_sha256": "c" * 64},
    }
    manifest["annotation_sources"] = {
        "actionability_annotation": f"ACMG SF v3.3; ClinGen snapshot 20260815; canonical-sha256:{digest}"
    }
    inventory = ["variants.parquet", "transformations.json", "extensions/family_evidence.parquet"]
    if fault == "inventory_declared":
        inventory.append(findings)
    if fault == "leftover_decisions":
        (bundle / "extensions/acmg_secondary_findings_evidence.parquet").touch()
    if fault == "disabled":
        del manifest["resource_provenance"]["acmg_secondary_findings"]
    if fault == "missing_clinvar":
        del manifest["resource_provenance"]["clinvar"]
    if fault == "not_run":
        manifest["annotation_sources"]["actionability_annotation"] = "not_run"
    if fault == "interrupted":
        stage["status"] = "pending"
    if fault == "nonactionable_decision":
        policy.update(decision_rows=1, decision_status_counts={"SECOND_ALLELE_REQUIRED": 1})
        stage["decision_rows"] = 1
    if fault == "missing_count":
        del stage["decision_rows"]
    if fault == "false_count":
        stage["decision_rows"] = False
    if fault == "dangling_symlink":
        (bundle / findings).symlink_to(bundle / "missing.parquet")
    docs = {
        "transformations.json": {"stages": [stage]},
        "inventory.json": {"entries": [{"path": name} for name in inventory]},
    }
    for name, document in docs.items():
        path = bundle / name
        write_json(path, document)
        manifest["files"][name] = {"sha256": sha256_file(path), "bytes": path.stat().st_size}
    if fault == "changed_evaluation":
        write_json(bundle / "transformations.json", {"stages": []})
    write_json(manifest_file, manifest)
    return bundle


@pytest.mark.parametrize("version", ["0.2.7", "0.2.6+release-c.1"])
def test_missing_legacy_clinical_table_requires_evaluated_zero(
    tmp_path: Path, version: str
) -> None:
    bundle = legacy_bundle(tmp_path, version=version)
    targets, receipt = tmp_path / "targets.tsv", tmp_path / "receipt.json"
    result = create_readback_targets(bundle, targets, receipt)
    assert result["selected_rows"] == 4
    assert result["required_clinical_rows"] == 0
    assert result["source_tables"]["clinical_findings_sha256"] is None
    assert result["source_tables"]["clinical_findings_state"] == "legacy_evaluated_empty"
    assert validate_readback_targets(targets, receipt, bundle_path=bundle)["status"] == "PASS"
    (bundle / "inventory.json").write_bytes(b"{}")
    with pytest.raises(ReadbackTargetError, match="inventory_identity_changed"):
        validate_readback_targets(targets, receipt, bundle_path=bundle)


@pytest.mark.parametrize(
    "fault",
    [
        "registered_missing",
        "inventory_declared",
        "disabled",
        "not_run",
        "missing_clinvar",
        "interrupted",
        "nonactionable_decision",
        "missing_count",
        "false_count",
        "leftover_decisions",
        "changed_evaluation",
        "dangling_symlink",
    ],
)
def test_selector_does_not_infer_zero_from_incomplete_sources(tmp_path: Path, fault: str) -> None:
    bundle = legacy_bundle(tmp_path, fault=fault)
    targets, receipt = tmp_path / "targets.tsv", tmp_path / "receipt.json"
    with pytest.raises(ReadbackTargetError):
        create_readback_targets(bundle, targets, receipt)
    assert not targets.exists() and not receipt.exists()


def test_current_producer_requires_materialized_findings(tmp_path: Path) -> None:
    bundle = legacy_bundle(tmp_path, version="0.2.8+release-c.1")
    with pytest.raises(ReadbackTargetError, match="clinical_findings_missing"):
        create_readback_targets(bundle, tmp_path / "targets.tsv", tmp_path / "receipt.json")


def test_registered_findings_from_other_sources_and_unlinked_findings(tmp_path: Path) -> None:
    bundle = synthetic_bundle(tmp_path)
    path = bundle / "clinical_findings.parquet"
    pq.write_table(pa.table({"variant_id": [None, "v-routine", "v-routine"]}), path)
    manifest = json.loads((bundle / "manifest.json").read_bytes())
    entry = {"bytes": path.stat().st_size, "rows": 3, "sha256": sha256_file(path)}
    manifest["files"][path.name] = entry
    write_json(bundle / "manifest.json", manifest)
    result = create_readback_targets(bundle, tmp_path / "targets.tsv", tmp_path / "receipt.json")
    assert result["required_clinical_rows"] == 1
    assert result["source_tables"]["clinical_findings_sha256"] == entry["sha256"]
    assert "clinical_findings_state" not in result["source_tables"]
    entry["rows"] = 1
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(ReadbackTargetError, match="clinical_findings_row_count_changed"):
        create_readback_targets(bundle, tmp_path / "new.tsv", tmp_path / "new.json")
