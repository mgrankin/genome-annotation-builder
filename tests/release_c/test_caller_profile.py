from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from genome_annotation_builder.release_c import (
    caller_profile,
    cnvpytor_call_capture,
    cnvpytor_combined_hdf5_probe,
    oci_runner,
)
from genome_annotation_builder.release_c.caller_profile import (
    load_caller_profile,
    load_prepared_caller_profile_run,
    prepare_caller_profile_run,
    validate_caller_profile_payload,
)
from genome_annotation_builder.release_c.caller_receipt import (
    load_prepared_caller_run_boundary,
    write_caller_run_receipt,
)
from genome_annotation_builder.release_c.oci_runner import run_oci_caller
from genome_annotation_builder.release_c.portable_model import PortableEvidenceError
from genome_annotation_builder.release_c.resource_tree import create_resource_tree_manifest

PRIVATE_BINDING = "privateid"
GATE_RUN_IDENTITY = "5" * 64
PROFILE_ROOT = Path(caller_profile.__file__).with_name("profiles")
PROFILE_PATH = PROFILE_ROOT / "manta-germline-1.6.0.json"
PROFILE_IDENTITY = "aaa137729145e44caac0befe96a272c26aeb94a964c5d20450d7c09c02f8cd33"
HLA_PROFILE_PATH = PROFILE_ROOT / "hla-la-short-read-1.0.4.json"
HLA_PROFILE_IDENTITY = "31d97da9f13758bdab1b68eadff1324c80ef46f0ef5ff46cf28478ff111f7815"
SPECHLA_PROFILE_PATH = PROFILE_ROOT / "spechla-short-read-bam-1.0.12.json"
SPECHLA_PROFILE_IDENTITY = "85ce0109fd6764337d833dd42b89aecc873c34ef795f15e78cedc4186350e2d1"
CNVPYTOR_PROFILE_PATH = PROFILE_ROOT / "cnvpytor-combined-rd-baf-1.3.2.json"
CNVPYTOR_PROFILE_IDENTITY = "63874b343264e8ae5c858f865e6c53aafd8d5ba0f5668492a8549c40d50e1e33"
DELLY_PROFILE_PATH = PROFILE_ROOT / "delly-short-read-2.5.1.json"
DELLY_PROFILE_IDENTITY = "d081743a9c59fb3a848fefbcf9287bd5825f484209d0c711bcaa2c64490d84c8"
EXPANSIONHUNTER_FEMALE_PROFILE_PATH = PROFILE_ROOT / "expansionhunter-grch38-female-5.0.0.json"
EXPANSIONHUNTER_FEMALE_PROFILE_IDENTITY = (
    "4e2a03346cbbb6c2c694f5df1967c5219ce5549d6d5e17698a6d5ce043e30a26"
)
EXPANSIONHUNTER_MALE_PROFILE_PATH = PROFILE_ROOT / "expansionhunter-grch38-male-5.0.0.json"
EXPANSIONHUNTER_MALE_PROFILE_IDENTITY = (
    "84e93fb1e24cbebaeaee24424642fbbe63090d8c6fdcbbb689cc615849dc0063"
)
MT_PROFILE_PATH = PROFILE_ROOT / "gatk-mitochondria-normal-shifted-4.6.2.0.json"
MT_PROFILE_IDENTITY = "69910e7de8dd0b7a7d806307455a9b717cbdd2e7faa5d740f47d5a6353a5ff72"
MT_POLICY_EXAMPLE = Path(__file__).resolve().parents[2] / "config" / "mt-call-policy.example.json"
READBACK_PROFILE_PATH = PROFILE_ROOT / "readback-bounded-pysam-0.24.0.json"
READBACK_PROFILE_IDENTITY = "4a1a8d2e478d1158edebbcc496e3a7dd98b706d2f53967e433cd8c4ec1f79080"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def file_binding(path: Path) -> dict[str, object]:
    return {
        "host_path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def specialize_profile_for_bindings(
    profile_path: Path, bindings_path: Path, output_path: Path
) -> Path:
    """Make a synthetic-only profile whose pins match synthetic fixture bytes.

    Shipped profiles stay pinned to production GRCh38 resources. Behavioural unit
    tests use this private copy so they never need the multi-gigabyte resources.
    """

    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    bindings = json.loads(bindings_path.read_text(encoding="utf-8"))["mounts"]
    for mount in profile["mounts"]:
        binding = bindings[mount["role"]]
        if "expected_sha256" in mount:
            mount["expected_sha256"] = binding["sha256"]
        if "expected_tree_manifest_sha256" in mount:
            mount["expected_tree_manifest_sha256"] = binding["tree_manifest_sha256"]
    write_json(output_path, profile)
    return output_path


def make_fixture(tmp_path: Path) -> dict[str, Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    alignment = tmp_path / "alignment.bam"
    alignment_index = tmp_path / "alignment.bam.bai"
    reference = tmp_path / "reference.fa"
    reference_fai = tmp_path / "reference.fa.fai"
    resources = tmp_path / "resources.json"
    alignment.write_bytes(b"synthetic-alignment\n")
    alignment_index.write_bytes(b"synthetic-index\n")
    reference.write_text(">1\nACGTACGT\n", encoding="ascii")
    reference_fai.write_text("1\t8\t3\t8\t9\n", encoding="ascii")
    write_json(resources, {"schema_version": 1, "status": "PASS"})
    gate = tmp_path / "input-gate.json"
    write_json(
        gate,
        {
            "schema_version": 1,
            "status": "PASS",
            "caller": "manta",
            "selected_alignment": {"sha256": sha256_file(alignment)},
            "reference": {"sha256": sha256_file(reference)},
            "run_identity_sha256": GATE_RUN_IDENTITY,
            "logical_subject_count": 1,
            "subject_identifier_persisted": False,
            "private_paths_persisted": False,
            "private_manifest_digest_persisted": False,
        },
    )
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    bindings = tmp_path / "bindings.json"
    write_json(
        bindings,
        {
            "schema_version": 1,
            "runtime": "docker",
            "user": "1000:1000",
            "minimum_free_bytes": 1,
            "mounts": {
                "selected_alignment": file_binding(alignment),
                "selected_alignment_index": file_binding(alignment_index),
                "reference_fasta": file_binding(reference),
                "reference_fai": file_binding(reference_fai),
                "resource_manifest": file_binding(resources),
                "workspace": {"host_path": str(workspace)},
            },
        },
    )
    test_profile = specialize_profile_for_bindings(
        PROFILE_PATH, bindings, tmp_path / "manta-test-profile.json"
    )
    return {
        "alignment": alignment,
        "alignment_index": alignment_index,
        "reference": reference,
        "reference_fai": reference_fai,
        "resources": resources,
        "gate": gate,
        "workspace": workspace,
        "bindings": bindings,
        "profile": test_profile,
    }


def make_hla_fixture(tmp_path: Path) -> dict[str, Path]:
    fixture = make_fixture(tmp_path)
    gate = json.loads(fixture["gate"].read_text(encoding="utf-8"))
    gate["caller"] = "hla_la"
    write_json(fixture["gate"], gate)
    graph_root = tmp_path / "graphs"
    graph = graph_root / "PRG_MHC_GRCh38_withIMGT"
    graph.mkdir(parents=True)
    (graph / "sequences.txt").write_text("synthetic graph sequences\n", encoding="utf-8")
    (graph / "serializedGRAPH").write_bytes(b"synthetic graph\n")
    graph_manifest = tmp_path / "graph-tree.json"
    graph_summary = create_resource_tree_manifest(graph_root, graph_manifest)
    bindings = json.loads(fixture["bindings"].read_text(encoding="utf-8"))
    bindings["mounts"].pop("reference_fai")
    bindings["mounts"]["hla_graph"] = {
        "host_path": str(graph_root),
        "tree_manifest_path": str(graph_manifest),
        "tree_manifest_bytes": graph_summary["manifest_bytes"],
        "tree_manifest_sha256": graph_summary["manifest_sha256"],
    }
    write_json(fixture["bindings"], bindings)
    fixture["profile"] = specialize_profile_for_bindings(
        HLA_PROFILE_PATH,
        fixture["bindings"],
        tmp_path / "hla-la-test-profile.json",
    )
    fixture["graph"] = graph_root
    fixture["named_graph"] = graph
    fixture["graph_manifest"] = graph_manifest
    return fixture


def make_spechla_fixture(tmp_path: Path) -> dict[str, Path]:
    fixture = make_fixture(tmp_path)
    gate = json.loads(fixture["gate"].read_text(encoding="utf-8"))
    gate["caller"] = "spechla"
    write_json(fixture["gate"], gate)
    database = tmp_path / "spechla-database"
    (database / "HLA").mkdir(parents=True)
    (database / "ref").mkdir()
    (database / "HLA" / "Allelelist.txt").write_text(
        "# synthetic allele catalogue\n", encoding="utf-8"
    )
    (database / "ref" / "hla.ref.extend.fa").write_text(">HLA-A*01:01\nACGT\n", encoding="ascii")
    database_manifest = tmp_path / "spechla-database-tree.json"
    database_summary = create_resource_tree_manifest(database, database_manifest)
    bindings = json.loads(fixture["bindings"].read_text(encoding="utf-8"))
    bindings["mounts"]["spechla_database"] = {
        "host_path": str(database),
        "tree_manifest_path": str(database_manifest),
        "tree_manifest_bytes": database_summary["manifest_bytes"],
        "tree_manifest_sha256": database_summary["manifest_sha256"],
    }
    write_json(fixture["bindings"], bindings)
    fixture["profile"] = specialize_profile_for_bindings(
        SPECHLA_PROFILE_PATH,
        fixture["bindings"],
        tmp_path / "spechla-test-profile.json",
    )
    fixture["database"] = database
    fixture["database_manifest"] = database_manifest
    return fixture


def make_cnvpytor_fixture(tmp_path: Path) -> dict[str, Path]:
    fixture = make_fixture(tmp_path)
    baf_vcf = tmp_path / "baf.vcf.gz"
    baf_index = tmp_path / "baf.vcf.gz.tbi"
    baf_vcf.write_bytes(b"synthetic one-sample BAF VCF\n")
    baf_index.write_bytes(b"synthetic BAF index\n")
    gate = json.loads(fixture["gate"].read_text(encoding="utf-8"))
    gate["caller"] = "cnvpytor"
    gate["baf_vcf"] = {
        "sha256": sha256_file(baf_vcf),
        "index_sha256": sha256_file(baf_index),
        "sample_count": 1,
        "sample_identity_match": True,
        "records": 1,
    }
    write_json(fixture["gate"], gate)

    package_data = tmp_path / "cnvpytor-package-data"
    package_data.mkdir()
    (package_data / "mask_hg38.pytor").write_bytes(b"synthetic strict mask\n")
    (package_data / "gc_hg38.pytor").write_bytes(b"synthetic package GC\n")
    package_manifest = tmp_path / "cnvpytor-package-data-tree.json"
    package_summary = create_resource_tree_manifest(package_data, package_manifest)
    gc_resource = tmp_path / "gc_exact_grch38.pytor"
    gc_resource.write_bytes(b"synthetic exact-reference GC\n")
    reference_config = tmp_path / "reference_genomes_conf.py"
    reference_config.write_text(
        "import_reference_genomes = {'hg38': {'gc_file': "
        "'/gc/gc_exact_grch38.pytor', 'mask_file': "
        "'/source-data/mask_hg38.pytor'}}\n",
        encoding="utf-8",
    )
    call_capture = Path(cnvpytor_call_capture.__file__)
    hdf5_probe = Path(cnvpytor_combined_hdf5_probe.__file__)
    bindings = json.loads(fixture["bindings"].read_text(encoding="utf-8"))
    bindings["mounts"].update(
        {
            "call_capture": file_binding(call_capture),
            "cnvpytor_package_data": {
                "host_path": str(package_data),
                "tree_manifest_path": str(package_manifest),
                "tree_manifest_bytes": package_summary["manifest_bytes"],
                "tree_manifest_sha256": package_summary["manifest_sha256"],
            },
            "cnvpytor_source_data": {
                "host_path": str(package_data),
                "tree_manifest_path": str(package_manifest),
                "tree_manifest_bytes": package_summary["manifest_bytes"],
                "tree_manifest_sha256": package_summary["manifest_sha256"],
            },
            "gc_resource": file_binding(gc_resource),
            "hdf5_probe": file_binding(hdf5_probe),
            "reference_config": file_binding(reference_config),
            "source_baf_vcf": file_binding(baf_vcf),
            "source_baf_vcf_index": file_binding(baf_index),
        }
    )
    write_json(fixture["bindings"], bindings)
    fixture["profile"] = specialize_profile_for_bindings(
        CNVPYTOR_PROFILE_PATH,
        fixture["bindings"],
        tmp_path / "cnvpytor-test-profile.json",
    )
    fixture.update(
        {
            "baf_vcf": baf_vcf,
            "baf_index": baf_index,
            "package_data": package_data,
            "package_manifest": package_manifest,
            "gc_resource": gc_resource,
            "reference_config": reference_config,
        }
    )
    return fixture


def test_shipped_manta_profile_is_canonical_path_free_and_pinned() -> None:
    profile = load_caller_profile(PROFILE_PATH, forbidden_bindings=(PRIVATE_BINDING,))
    assert profile.identity_sha256 == PROFILE_IDENTITY
    summary = profile.portable_summary()
    assert summary == {
        "schema_version": 1,
        "status": "PASS",
        "profile_id": "manta-germline-1.6.0",
        "profile_identity_sha256": PROFILE_IDENTITY,
        "caller_name": "Manta",
        "caller_version": "1.6.0",
        "gate_caller": "manta",
        "image_identity_sha256": (
            "021aabe5f7314e02e60722c559f6f3c4bb6b1118c2e2ec0eb3095a8c7e431f26"
        ),
        "mounts": 6,
        "steps": 3,
        "private_paths_persisted": False,
    }
    retained = json.dumps(profile.payload, sort_keys=True)
    assert "/raid/" not in retained
    assert "/home/" not in retained
    assert PRIVATE_BINDING not in retained


def test_shipped_hlala_profile_uses_generic_alias_and_separate_graph_mount() -> None:
    profile = load_caller_profile(HLA_PROFILE_PATH, forbidden_bindings=(PRIVATE_BINDING,))
    assert profile.identity_sha256 == HLA_PROFILE_IDENTITY
    assert (
        profile.payload["command_semantics"]["parameters"]["caller_visible_sample_alias"]
        == "SUBJECT"
    )
    call = profile.payload["command_semantics"]["steps"][1]
    assert call["entrypoint"] == "HLA-LA.pl"
    assert call["arguments"][call["arguments"].index("--sampleID") + 1] == "SUBJECT"
    graph_mount = next(mount for mount in profile.payload["mounts"] if mount["role"] == "hla_graph")
    assert graph_mount == {
        "role": "hla_graph",
        "kind": "directory",
        "mode": "ro",
        "container_path": "/graphs",
        "expected_tree_manifest_sha256": (
            "37939ec3ad442f8f06e8e158abd3cc4d002f96a3f2e405f51c9056d3ae8f31c8"
        ),
    }
    retained = json.dumps(profile.payload, sort_keys=True)
    assert PRIVATE_BINDING not in retained
    assert "/raid/" not in retained


def test_shipped_spechla_profile_is_bam_only_uses_generic_alias_and_exact_database() -> None:
    profile = load_caller_profile(SPECHLA_PROFILE_PATH, forbidden_bindings=(PRIVATE_BINDING,))
    assert profile.identity_sha256 == SPECHLA_PROFILE_IDENTITY
    assert profile.portable_summary() == {
        "schema_version": 1,
        "status": "PASS",
        "profile_id": "spechla-short-read-bam-1.0.12",
        "profile_identity_sha256": SPECHLA_PROFILE_IDENTITY,
        "caller_name": "SpecHLA",
        "caller_version": "1.0.12",
        "gate_caller": "spechla",
        "image_identity_sha256": (
            "c6b76272e17f4f3f91efebe96bdbad52dd8a965dae64a50d33edf432f1a96144"
        ),
        "mounts": 7,
        "steps": 3,
        "private_paths_persisted": False,
    }
    parameters = profile.payload["command_semantics"]["parameters"]
    assert parameters["alignment_contract"] == "coordinate-sorted indexed BAM"
    assert parameters["caller_visible_sample_alias"] == "SUBJECT"
    assert parameters["population_prior"] == "Unknown"
    extraction, call = profile.payload["command_semantics"]["steps"][1:]
    assert parameters["image_environment_bootstrap"] == "env-execute"
    assert extraction["entrypoint"] == "bash"
    assert extraction["arguments"][3] == "spechla-extract-hla-reads"
    assert extraction["arguments"][extraction["arguments"].index("-s") + 1] == "SUBJECT"
    assert "{reference_fasta}" in extraction["arguments"]
    assert call["entrypoint"] == "bash"
    assert call["arguments"][3] == "spechla"
    assert call["arguments"][call["arguments"].index("-n") + 1] == "SUBJECT"
    assert call["arguments"][call["arguments"].index("-i") + 1] == "{spechla_database}"
    assert profile.payload["artifact"]["identity_scheme"] == "named_sha256_set"
    assert set(profile.payload["artifact"]["components"]) == {
        "full_result",
        "g_group_result",
        "details_result",
    }
    retained = json.dumps(profile.payload, sort_keys=True)
    assert PRIVATE_BINDING not in retained
    assert "/raid/" not in retained
    database_mount = next(
        mount for mount in profile.payload["mounts"] if mount["role"] == "spechla_database"
    )
    assert database_mount["expected_tree_manifest_sha256"] == (
        "f7f593c543d38ca4e668bff0004c05f2594b00e0b03751dc472b0601a4bce64c"
    )


def test_shipped_cnvpytor_profile_pins_false_zero_guard_and_two_exact_bin_sizes() -> None:
    profile = load_caller_profile(CNVPYTOR_PROFILE_PATH, forbidden_bindings=(PRIVATE_BINDING,))
    assert profile.identity_sha256 == CNVPYTOR_PROFILE_IDENTITY
    assert profile.portable_summary() == {
        "schema_version": 1,
        "status": "PASS",
        "profile_id": "cnvpytor-combined-rd-baf-1.3.2",
        "profile_identity_sha256": CNVPYTOR_PROFILE_IDENTITY,
        "caller_name": "CNVpytor",
        "caller_version": "1.3.2",
        "gate_caller": "cnvpytor",
        "image_identity_sha256": (
            "614b22c8fc1fcef147bcda48bb4d53cb63c988bffc588e186acff7335013fcf1"
        ),
        "mounts": 14,
        "steps": 12,
        "private_paths_persisted": False,
    }
    parameters = profile.payload["command_semantics"]["parameters"]
    assert parameters["bin_sizes"] == [10000, 100000]
    assert parameters["false_zero_exit_guard"] is True
    assert parameters["rd_strict_mask"] is True
    assert parameters["baf_strict_mask"] is True
    assert parameters["baf_sample_policy"].startswith("exactly one VCF sample")
    steps = profile.payload["command_semantics"]["steps"]
    assert [step["name"] for step in steps[-5:]] == [
        "capture_combined_10000",
        "audit_hdf5_10000",
        "capture_combined_100000",
        "audit_hdf5_100000",
        "audit_final_hdf5_10000",
    ]
    assert all(step["entrypoint"] == "python3" for step in steps[-5:])
    assert profile.payload["artifact"]["components"]["hdf5_audit_10000"].endswith(
        "combined-10000.final.hdf5-audit.json"
    )
    assert profile.payload["artifact"]["identity_scheme"] == "named_sha256_set"
    assert set(profile.payload["artifact"]["components"]) == {
        "pytor",
        "calls_10000",
        "calls_100000",
        "hdf5_audit_10000",
        "hdf5_audit_100000",
    }
    retained = json.dumps(profile.payload, sort_keys=True)
    assert PRIVATE_BINDING not in retained
    assert "/raid/" not in retained


def test_shipped_delly_profile_binds_static_binary_and_exact_reference_mask() -> None:
    profile = load_caller_profile(DELLY_PROFILE_PATH, forbidden_bindings=(PRIVATE_BINDING,))
    assert profile.identity_sha256 == DELLY_PROFILE_IDENTITY
    assert profile.portable_summary() == {
        "schema_version": 1,
        "status": "PASS",
        "profile_id": "delly-short-read-2.5.1",
        "profile_identity_sha256": DELLY_PROFILE_IDENTITY,
        "caller_name": "DELLY",
        "caller_version": "2.5.1",
        "gate_caller": "delly",
        "image_identity_sha256": (
            "561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea"
        ),
        "mounts": 8,
        "steps": 2,
        "private_paths_persisted": False,
    }
    parameters = profile.payload["command_semantics"]["parameters"]
    assert parameters["binary_sha256"] == (
        "98bab8261663f6bead10ebf9c864e22b71f087244350f098c7d631625e65e648"
    )
    command = profile.payload["command_semantics"]["steps"][1]
    assert command["entrypoint"] == "delly"
    assert command["arguments"][0] == "sr"
    assert "{delly_exclude}" in command["arguments"]
    assert "{reference_fasta}" in command["arguments"]
    binary_mount = next(
        mount for mount in profile.payload["mounts"] if mount["role"] == "delly_binary"
    )
    assert binary_mount["container_path"] == "/usr/local/bin/delly"
    assert profile.payload["artifact"]["components"] == {"primary_bcf": "delly/SUBJECT.sr.bcf"}


@pytest.mark.parametrize(
    ("profile_path", "identity", "sex_model"),
    [
        (
            EXPANSIONHUNTER_FEMALE_PROFILE_PATH,
            EXPANSIONHUNTER_FEMALE_PROFILE_IDENTITY,
            "female",
        ),
        (
            EXPANSIONHUNTER_MALE_PROFILE_PATH,
            EXPANSIONHUNTER_MALE_PROFILE_IDENTITY,
            "male",
        ),
    ],
)
def test_shipped_expansionhunter_profiles_pin_binary_catalog_and_sex_model(
    profile_path: Path, identity: str, sex_model: str
) -> None:
    profile = load_caller_profile(profile_path, forbidden_bindings=(PRIVATE_BINDING,))
    assert profile.identity_sha256 == identity
    assert profile.payload["caller_name"] == "ExpansionHunter"
    assert profile.payload["caller_version"] == "5.0.0"
    assert profile.payload["gate_caller"] == "expansionhunter"
    parameters = profile.payload["command_semantics"]["parameters"]
    assert parameters["binary_sha256"] == (
        "3dc7a14aa869b93bc45be413cc65a0c17801f7be8c5097fa7d86e608fefc8762"
    )
    assert parameters["sex_model"] == sex_model
    command = profile.payload["command_semantics"]["steps"][1]
    assert command["entrypoint"] == "ExpansionHunter"
    assert command["arguments"][command["arguments"].index("--sex") + 1] == sex_model
    assert "{variant_catalog}" in command["arguments"]
    assert profile.payload["artifact"] == {
        "raw_artifact_id": "expansionhunter-call-pair",
        "identity_scheme": "named_sha256_set",
        "components": {
            "expansionhunter_json_sha256": "expansionhunter/SUBJECT.json",
            "expansionhunter_vcf_sha256": "expansionhunter/SUBJECT.vcf",
        },
        "derived_payload": None,
    }
    retained = json.dumps(profile.payload, sort_keys=True)
    assert PRIVATE_BINDING not in retained
    assert "/raid/" not in retained


def test_shipped_mitochondrial_profile_pins_normal_shifted_workflow() -> None:
    profile = load_caller_profile(MT_PROFILE_PATH, forbidden_bindings=(PRIVATE_BINDING,))
    assert profile.identity_sha256 == MT_PROFILE_IDENTITY
    assert profile.portable_summary() == {
        "schema_version": 1,
        "status": "PASS",
        "profile_id": "gatk-mitochondria-normal-shifted-4.6.2.0",
        "profile_identity_sha256": MT_PROFILE_IDENTITY,
        "caller_name": "GATK Mutect2 mitochondrial mode",
        "caller_version": "4.6.2.0",
        "gate_caller": "gatk_mitochondria",
        "image_identity_sha256": (
            "71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb"
        ),
        "mounts": 11,
        "steps": 1,
        "private_paths_persisted": False,
    }
    parameters = profile.payload["command_semantics"]["parameters"]
    assert parameters["mt_length"] == 16569
    assert parameters["shift_bases"] == 8000
    assert "chrM:576-16024" in parameters["circular_boundary_control"]
    assert "chrM:8025-9144" in parameters["circular_boundary_control"]
    step = profile.payload["command_semantics"]["steps"][0]
    assert step["entrypoint"] == "python3"
    assert step["arguments"][0] == "{mt_workflow}"
    assert "{mt_resource_tree}" in step["arguments"]
    resource_mount = next(
        mount for mount in profile.payload["mounts"] if mount["role"] == "mt_resource_tree"
    )
    assert resource_mount["expected_tree_manifest_sha256"] == (
        "aec2a4d3bc1919a7d5896dc768cc584a82eaf690fc02a4fce770b25c51a0aa36"
    )
    policy_mount = next(
        mount for mount in profile.payload["mounts"] if mount["role"] == "mt_call_policy"
    )
    assert sha256_file(MT_POLICY_EXAMPLE) == policy_mount["expected_sha256"]
    assert json.loads(MT_POLICY_EXAMPLE.read_text(encoding="utf-8"))["source_contig"] == "MT"
    assert profile.payload["artifact"] == {
        "raw_artifact_id": "gatk-mitochondrial-filtered-vcf",
        "identity_scheme": "file_sha256",
        "components": {"filtered_vcf": "mt/filtered.vcf"},
        "derived_payload": None,
    }
    retained = json.dumps(profile.payload, sort_keys=True)
    assert PRIVATE_BINDING not in retained
    assert "/raid/" not in retained


def test_shipped_readback_profile_pins_capture_code_runtime_and_policy() -> None:
    profile = load_caller_profile(READBACK_PROFILE_PATH)
    assert profile.identity_sha256 == READBACK_PROFILE_IDENTITY
    assert profile.portable_summary() == {
        "schema_version": 1,
        "status": "PASS",
        "profile_id": "readback-bounded-pysam-0.24.0",
        "profile_identity_sha256": READBACK_PROFILE_IDENTITY,
        "caller_name": "GAB bounded local pileup",
        "caller_version": "1",
        "gate_caller": "readback",
        "image_identity_sha256": (
            "3097afa2a5d7dbd9cc3e2b1d53adf3f0e9ef1a1a3bdb5076a88c0bf8547f290b"
        ),
        "mounts": 8,
        "steps": 2,
        "private_paths_persisted": False,
    }
    capture_mount = next(
        mount for mount in profile.payload["mounts"] if mount["role"] == "readback_capture"
    )
    assert capture_mount["expected_sha256"] == (
        "b4c87b153d91250f9fa623c3b74d55bbb6812e96fdcf0a0f846012ed866d1b49"
    )
    parameters = profile.payload["command_semantics"]["parameters"]
    assert parameters["minimum_mapping_quality"] == 20
    assert parameters["minimum_base_quality"] == 20
    assert parameters["minimum_alt_fragments"] == 3
    assert parameters["duplicate_reads_included"] is False
    assert profile.payload["artifact"] == {
        "raw_artifact_id": "bounded-readback-capture-json",
        "identity_scheme": "file_sha256",
        "components": {"raw_capture_json": "readback/capture.json"},
        "derived_payload": None,
    }


def test_prepare_compiles_exact_private_specs_and_restart_marker(tmp_path: Path) -> None:
    fixture = make_fixture(tmp_path)
    profile_path = fixture["profile"]
    profile_identity = load_caller_profile(profile_path).identity_sha256
    run_dir = tmp_path / "run-plan"
    summary = prepare_caller_profile_run(
        profile_path,
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    assert summary["status"] == "PREPARED"
    assert summary["profile_identity_sha256"] == profile_identity
    assert summary["run_plan_scope"] == "private_runtime"
    assert summary["private_paths_persisted_in_run_plan"] is True
    assert summary["private_paths_persisted_in_portable_summary"] is False
    assert load_prepared_caller_profile_run(run_dir, forbidden_bindings=(PRIVATE_BINDING,)) == {
        "schema_version": 1,
        "status": "PREPARED",
        "profile_id": "manta-germline-1.6.0",
        "profile_identity_sha256": profile_identity,
        "caller_run_identity_sha256": summary["caller_run_identity_sha256"],
        "run_files": 5,
        "subject_identifier_persisted": False,
        "run_plan_scope": "private_runtime",
        "private_paths_persisted_in_run_plan": True,
        "private_paths_persisted_in_portable_summary": False,
    }
    checkpoint, _ = load_prepared_caller_run_boundary(run_dir / "prepared-boundary.json")
    assert checkpoint["run_boundary"]["command_semantics"]["parameters"] == {
        "jobs": 8,
        "mode": "local",
        "workflow": "germline diploid SV",
        "profile_id": "manta-germline-1.6.0",
        "profile_identity_sha256": profile_identity,
    }
    execution = json.loads((run_dir / "execution-spec.json").read_text(encoding="utf-8"))
    assert execution["image"].endswith(
        "@sha256:021aabe5f7314e02e60722c559f6f3c4bb6b1118c2e2ec0eb3095a8c7e431f26"
    )
    assert {mount["role"] for mount in execution["mounts"]} == {
        "selected_alignment",
        "selected_alignment_index",
        "reference_fasta",
        "reference_fai",
        "resource_manifest",
        "workspace",
    }
    artifact = json.loads((run_dir / "artifact-spec.json").read_text(encoding="utf-8"))
    assert artifact["components"] == {
        "primary_vcf": str(fixture["workspace"] / "manta/results/variants/diploidSV.vcf.gz")
    }


def test_profile_plan_runs_through_oci_and_finalizes_exact_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_fixture(tmp_path)
    profile_path = fixture["profile"]
    run_dir = tmp_path / "run-plan"
    prepare_caller_profile_run(
        profile_path,
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
    )
    commands: list[list[str]] = []
    profile = load_caller_profile(profile_path)

    def capture(command: list[str]) -> subprocess.CompletedProcess[str]:
        if command[1:] == ["--version"]:
            return subprocess.CompletedProcess(command, 0, "Docker synthetic\n", "")
        return subprocess.CompletedProcess(
            command,
            0,
            json.dumps([profile.payload["image"]]) + "\n",
            "",
        )

    def execute(command: list[str], stdout_handle, stderr_handle) -> int:
        commands.append(command)
        stdout_handle.write(b"synthetic caller output\n")
        stderr_handle.write(b"")
        if len(commands) == 3:
            output = fixture["workspace"] / "manta/results/variants/diploidSV.vcf.gz"
            output.parent.mkdir(parents=True)
            output.write_bytes(b"synthetic diploid SV VCF\n")
        return 0

    monkeypatch.setattr(oci_runner, "_capture", capture)
    monkeypatch.setattr(oci_runner, "_execute", execute)
    report_path = tmp_path / "execution-report.json"
    report = run_oci_caller(
        run_dir / "prepared-boundary.json",
        run_dir / "execution-spec.json",
        report_path,
        tmp_path / "logs",
    )
    assert report["status"] == "PASS"
    assert len(commands) == 3
    assert [command[command.index("--entrypoint") + 1] for command in commands] == [
        "mkdir",
        "configManta.py",
        "python2",
    ]
    artifact = json.loads((run_dir / "artifact-spec.json").read_text(encoding="utf-8"))
    checkpoint, _ = load_prepared_caller_run_boundary(run_dir / "prepared-boundary.json")
    receipt = write_caller_run_receipt(
        tmp_path / "caller-receipt.json",
        checkpoint["run_boundary"],
        fixture["gate"],
        fixture["resources"],
        execution_report=report_path,
        raw_artifact_id=artifact["raw_artifact_id"],
        identity_scheme=artifact["identity_scheme"],
        component_paths=artifact["components"],
    )
    primary = Path(artifact["components"]["primary_vcf"])
    assert receipt["raw_artifact"]["raw_artifact_sha256"] == sha256_file(primary)
    assert (
        receipt["run_boundary"]["caller_run_identity_sha256"]
        == checkpoint["run_boundary"]["caller_run_identity_sha256"]
    )


def test_hlala_profile_tree_binding_runs_with_no_private_sample_alias(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_hla_fixture(tmp_path)
    profile_path = fixture["profile"]
    run_dir = tmp_path / "run-plan"
    prepare_caller_profile_run(
        profile_path,
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    profile = load_caller_profile(profile_path)
    commands: list[list[str]] = []

    def capture(command: list[str]) -> subprocess.CompletedProcess[str]:
        if command[1:] == ["--version"]:
            return subprocess.CompletedProcess(command, 0, "Docker synthetic\n", "")
        return subprocess.CompletedProcess(
            command,
            0,
            json.dumps([profile.payload["image"]]) + "\n",
            "",
        )

    def execute(command: list[str], stdout_handle, stderr_handle) -> int:
        commands.append(command)
        stdout_handle.write(b"synthetic caller output\n")
        stderr_handle.write(b"")
        if len(commands) == 2:
            output = fixture["workspace"] / "SUBJECT/hla/R1_bestguess_G.txt"
            output.parent.mkdir(parents=True)
            output.write_text("Locus\tAllele\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(oci_runner, "_capture", capture)
    monkeypatch.setattr(oci_runner, "_execute", execute)
    report_path = tmp_path / "execution-report.json"
    report = run_oci_caller(
        run_dir / "prepared-boundary.json",
        run_dir / "execution-spec.json",
        report_path,
        tmp_path / "logs",
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    assert report["status"] == "PASS"
    graph_evidence = next(
        mount for mount in report["input_verification"]["mounts"] if mount["role"] == "hla_graph"
    )
    assert graph_evidence["tree_files"] == 2
    assert graph_evidence["tree_directories"] == 1
    assert [command[command.index("--entrypoint") + 1] for command in commands] == [
        "mkdir",
        "HLA-LA.pl",
    ]
    assert "SUBJECT" in commands[1]
    payload = report_path.read_text(encoding="utf-8")
    assert PRIVATE_BINDING not in payload
    assert str(fixture["workspace"]) not in payload


def test_spechla_profile_runs_two_stage_contract_and_receipts_three_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_spechla_fixture(tmp_path / PRIVATE_BINDING)
    profile_path = fixture["profile"]
    run_dir = tmp_path / PRIVATE_BINDING / "run-plan"
    prepare_caller_profile_run(
        profile_path,
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    profile = load_caller_profile(profile_path)
    commands: list[list[str]] = []

    def capture(command: list[str]) -> subprocess.CompletedProcess[str]:
        if command[1:] == ["--version"]:
            return subprocess.CompletedProcess(command, 0, "Docker synthetic\n", "")
        return subprocess.CompletedProcess(
            command,
            0,
            json.dumps([profile.payload["image"]]) + "\n",
            "",
        )

    def execute(command: list[str], stdout_handle, stderr_handle) -> int:
        commands.append(command)
        stdout_handle.write(b"synthetic caller output\n")
        stderr_handle.write(b"")
        if len(commands) == 3:
            result_dir = fixture["workspace"] / "results/SUBJECT"
            result_dir.mkdir(parents=True)
            (result_dir / "hla.result.txt").write_text(
                "Sample\tLocus\tChromosome\tAllele\n", encoding="utf-8"
            )
            (result_dir / "hla.result.g.group.txt").write_text(
                "Sample\tLocus\tChromosome\tAllele\n", encoding="utf-8"
            )
            (result_dir / "hla.result.details.txt").write_text(
                "Sample\tLocus\tChromosome\tAllele\tScore\n", encoding="utf-8"
            )
        return 0

    monkeypatch.setattr(oci_runner, "_capture", capture)
    monkeypatch.setattr(oci_runner, "_execute", execute)
    report_path = tmp_path / "execution-report.json"
    report = run_oci_caller(
        run_dir / "prepared-boundary.json",
        run_dir / "execution-spec.json",
        report_path,
        tmp_path / "logs",
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    assert report["status"] == "PASS"
    assert [command[command.index("--entrypoint") + 1] for command in commands] == [
        "mkdir",
        "bash",
        "bash",
    ]
    assert "spechla-extract-hla-reads" in commands[1]
    assert "spechla" in commands[2]
    assert "SUBJECT" in commands[1]
    assert "SUBJECT" in commands[2]
    database_evidence = next(
        mount
        for mount in report["input_verification"]["mounts"]
        if mount["role"] == "spechla_database"
    )
    assert database_evidence["tree_files"] == 2
    assert database_evidence["tree_directories"] == 2

    artifact = json.loads((run_dir / "artifact-spec.json").read_text(encoding="utf-8"))
    checkpoint, _ = load_prepared_caller_run_boundary(run_dir / "prepared-boundary.json")
    receipt = write_caller_run_receipt(
        tmp_path / "caller-receipt.json",
        checkpoint["run_boundary"],
        fixture["gate"],
        fixture["resources"],
        execution_report=report_path,
        raw_artifact_id=artifact["raw_artifact_id"],
        identity_scheme=artifact["identity_scheme"],
        component_paths=artifact["components"],
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    assert receipt["raw_artifact"]["identity_scheme"] == "named_sha256_set"
    assert [component["role"] for component in receipt["raw_artifact"]["components"]] == [
        "details_result",
        "full_result",
        "g_group_result",
    ]
    portable_report = report_path.read_text(encoding="utf-8")
    portable_receipt = json.dumps(receipt, sort_keys=True)
    assert PRIVATE_BINDING not in portable_report
    assert PRIVATE_BINDING not in portable_receipt
    assert str(fixture["workspace"]) not in portable_report
    assert str(fixture["workspace"]) not in portable_receipt


def test_cnvpytor_profile_runs_exact_mount_contract_and_receipts_all_audited_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_cnvpytor_fixture(tmp_path / PRIVATE_BINDING)
    profile_path = fixture["profile"]
    run_dir = tmp_path / PRIVATE_BINDING / "run-plan"
    prepare_caller_profile_run(
        profile_path,
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    profile = load_caller_profile(profile_path)
    commands: list[list[str]] = []

    def capture(command: list[str]) -> subprocess.CompletedProcess[str]:
        if command[1:] == ["--version"]:
            return subprocess.CompletedProcess(command, 0, "Docker synthetic\n", "")
        return subprocess.CompletedProcess(
            command,
            0,
            json.dumps([profile.payload["image"]]) + "\n",
            "",
        )

    def execute(command: list[str], stdout_handle, stderr_handle) -> int:
        commands.append(command)
        stdout_handle.write(b"synthetic caller output\n")
        stderr_handle.write(b"")
        if len(commands) == 12:
            result_dir = fixture["workspace"] / "cnvpytor"
            result_dir.mkdir(parents=True, exist_ok=True)
            (result_dir / "SUBJECT.pytor").write_bytes(b"synthetic CNVpytor HDF5\n")
            for bin_size in (10000, 100000):
                (result_dir / f"combined-{bin_size}.calls.tsv").write_text(
                    "synthetic validated combined row\n", encoding="utf-8"
                )
                audit_name = (
                    "combined-10000.final.hdf5-audit.json"
                    if bin_size == 10000
                    else "combined-100000.hdf5-audit.json"
                )
                (result_dir / audit_name).write_text(
                    '{"schema_version": 1, "status": "PASS"}\n', encoding="utf-8"
                )
        return 0

    monkeypatch.setattr(oci_runner, "_capture", capture)
    monkeypatch.setattr(oci_runner, "_execute", execute)
    report_path = tmp_path / "execution-report.json"
    report = run_oci_caller(
        run_dir / "prepared-boundary.json",
        run_dir / "execution-spec.json",
        report_path,
        tmp_path / "logs",
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    assert report["status"] == "PASS"
    assert len(commands) == 12
    assert [command[command.index("--entrypoint") + 1] for command in commands] == [
        "mkdir",
        "cnvpytor",
        "cnvpytor",
        "cnvpytor",
        "cnvpytor",
        "cnvpytor",
        "cnvpytor",
        "python3",
        "python3",
        "python3",
        "python3",
        "python3",
    ]
    tree_evidence = {
        mount["role"]: mount
        for mount in report["input_verification"]["mounts"]
        if mount["kind"] == "directory" and mount["mode"] == "ro"
    }
    assert set(tree_evidence) == {"cnvpytor_package_data", "cnvpytor_source_data"}
    assert all(mount["tree_files"] == 2 for mount in tree_evidence.values())

    artifact = json.loads((run_dir / "artifact-spec.json").read_text(encoding="utf-8"))
    checkpoint, _ = load_prepared_caller_run_boundary(run_dir / "prepared-boundary.json")
    receipt = write_caller_run_receipt(
        tmp_path / "caller-receipt.json",
        checkpoint["run_boundary"],
        fixture["gate"],
        fixture["resources"],
        execution_report=report_path,
        raw_artifact_id=artifact["raw_artifact_id"],
        identity_scheme=artifact["identity_scheme"],
        component_paths=artifact["components"],
        forbidden_bindings=(PRIVATE_BINDING,),
    )
    assert receipt["raw_artifact"]["identity_scheme"] == "named_sha256_set"
    assert len(receipt["raw_artifact"]["components"]) == 5
    portable_report = report_path.read_text(encoding="utf-8")
    portable_receipt = json.dumps(receipt, sort_keys=True)
    assert PRIVATE_BINDING not in portable_report
    assert PRIVATE_BINDING not in portable_receipt
    assert str(fixture["workspace"]) not in portable_report
    assert str(fixture["workspace"]) not in portable_receipt


def test_profile_and_bindings_fail_closed_before_publication(tmp_path: Path) -> None:
    profile_payload = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    profile_payload["licence"] = f"contains-{PRIVATE_BINDING}-binding"
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        validate_caller_profile_payload(profile_payload, forbidden_bindings=(PRIVATE_BINDING,))
    profile_payload = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    profile_payload["placeholder_bindings"]["{unused}"] = "/work"
    with pytest.raises(PortableEvidenceError, match="missing or unused"):
        validate_caller_profile_payload(profile_payload)

    fixture = make_fixture(tmp_path)
    bindings = json.loads(fixture["bindings"].read_text(encoding="utf-8"))
    bindings["mounts"].pop("reference_fai")
    write_json(fixture["bindings"], bindings)
    output = tmp_path / "invalid-run"
    with pytest.raises(PortableEvidenceError, match="missing or unexpected"):
        prepare_caller_profile_run(
            fixture["profile"],
            fixture["bindings"],
            fixture["gate"],
            fixture["resources"],
            output,
        )
    assert not output.exists()


def test_profile_file_digest_constraint_is_enforced_before_run_publication(
    tmp_path: Path,
) -> None:
    fixture = make_fixture(tmp_path)
    profile_payload = json.loads(fixture["profile"].read_text(encoding="utf-8"))
    alignment_mount = next(
        mount for mount in profile_payload["mounts"] if mount["role"] == "selected_alignment"
    )
    alignment_mount["expected_sha256"] = "f" * 64
    constrained_profile = tmp_path / "constrained-profile.json"
    write_json(constrained_profile, profile_payload)
    output = tmp_path / "must-not-publish"
    with pytest.raises(PortableEvidenceError, match="profile digest constraint"):
        prepare_caller_profile_run(
            constrained_profile,
            fixture["bindings"],
            fixture["gate"],
            fixture["resources"],
            output,
        )
    assert not output.exists()


def test_profile_tree_manifest_constraint_is_enforced_before_run_publication(
    tmp_path: Path,
) -> None:
    fixture = make_hla_fixture(tmp_path)
    profile_payload = json.loads(fixture["profile"].read_text(encoding="utf-8"))
    graph_mount = next(mount for mount in profile_payload["mounts"] if mount["role"] == "hla_graph")
    graph_mount["expected_tree_manifest_sha256"] = "f" * 64
    constrained_profile = tmp_path / "tree-constrained-profile.json"
    write_json(constrained_profile, profile_payload)
    output = tmp_path / "must-not-publish-tree-run"
    with pytest.raises(PortableEvidenceError, match="profile digest constraint"):
        prepare_caller_profile_run(
            constrained_profile,
            fixture["bindings"],
            fixture["gate"],
            fixture["resources"],
            output,
        )
    assert not output.exists()


def test_boundary_mismatch_empty_workspace_and_no_replace_are_enforced(
    tmp_path: Path,
) -> None:
    fixture = make_fixture(tmp_path)
    bindings = json.loads(fixture["bindings"].read_text(encoding="utf-8"))
    bindings["mounts"]["selected_alignment"]["sha256"] = "0" * 64
    write_json(fixture["bindings"], bindings)
    with pytest.raises(PortableEvidenceError, match="differ from the prepared input"):
        prepare_caller_profile_run(
            fixture["profile"],
            fixture["bindings"],
            fixture["gate"],
            fixture["resources"],
            tmp_path / "mismatch-run",
        )

    fixture = make_fixture(tmp_path / "not-empty")
    (fixture["workspace"] / "stale-output.txt").write_text("stale\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="workspace must be empty"):
        prepare_caller_profile_run(
            fixture["profile"],
            fixture["bindings"],
            fixture["gate"],
            fixture["resources"],
            tmp_path / "not-empty-run",
        )

    fixture = make_fixture(tmp_path / "valid")
    run_dir = tmp_path / "valid-run"
    prepare_caller_profile_run(
        fixture["profile"],
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
    )
    original_marker = (run_dir / "PREPARED.json").read_bytes()
    with pytest.raises(PortableEvidenceError, match="refusing to replace"):
        prepare_caller_profile_run(
            fixture["profile"],
            fixture["bindings"],
            fixture["gate"],
            fixture["resources"],
            run_dir,
        )
    assert (run_dir / "PREPARED.json").read_bytes() == original_marker


def test_prepared_run_detects_any_spec_or_marker_tamper(tmp_path: Path) -> None:
    fixture = make_fixture(tmp_path)
    run_dir = tmp_path / "run-plan"
    prepare_caller_profile_run(
        fixture["profile"],
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
    )
    execution = run_dir / "execution-spec.json"
    execution.write_bytes(execution.read_bytes() + b" ")
    with pytest.raises(PortableEvidenceError, match="file changed"):
        load_prepared_caller_profile_run(run_dir)

    fixture = make_fixture(tmp_path / "marker")
    run_dir = tmp_path / "marker-run"
    prepare_caller_profile_run(
        fixture["profile"],
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
    )
    marker_path = run_dir / "PREPARED.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    marker["profile_identity_sha256"] = "0" * 64
    write_json(marker_path, marker)
    with pytest.raises(PortableEvidenceError):
        load_prepared_caller_profile_run(run_dir)


def test_boundary_change_during_prepare_never_publishes_completion_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_fixture(tmp_path)
    original_writer = caller_profile.write_prepared_caller_run_boundary

    def mutate_gate_then_write(*args, **kwargs):
        gate = json.loads(fixture["gate"].read_text(encoding="utf-8"))
        gate["run_identity_sha256"] = "6" * 64
        write_json(fixture["gate"], gate)
        return original_writer(*args, **kwargs)

    monkeypatch.setattr(
        caller_profile, "write_prepared_caller_run_boundary", mutate_gate_then_write
    )
    run_dir = tmp_path / "raced-run"
    with pytest.raises(PortableEvidenceError, match="boundary changed"):
        prepare_caller_profile_run(
            fixture["profile"],
            fixture["bindings"],
            fixture["gate"],
            fixture["resources"],
            run_dir,
        )
    assert run_dir.is_dir()
    assert not (run_dir / "PREPARED.json").exists()


def test_prepare_fsyncs_completion_directory_and_parent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_fixture(tmp_path)
    run_dir = tmp_path / "run-plan"
    synchronized: list[Path] = []

    monkeypatch.setattr(
        caller_profile,
        "_fsync_directory",
        lambda path: synchronized.append(Path(path)),
    )
    prepare_caller_profile_run(
        fixture["profile"],
        fixture["bindings"],
        fixture["gate"],
        fixture["resources"],
        run_dir,
    )

    assert synchronized == [run_dir, run_dir.parent]
    assert (run_dir / "PREPARED.json").is_file()


def test_cli_success_and_failure_do_not_echo_private_bindings(tmp_path: Path) -> None:
    fixture = make_fixture(tmp_path / PRIVATE_BINDING)
    run_dir = tmp_path / PRIVATE_BINDING / "run-plan"
    module = "genome_annotation_builder.release_c.caller_profile"
    prepared = subprocess.run(
        [
            sys.executable,
            "-m",
            module,
            "prepare",
            "--profile",
            str(fixture["profile"]),
            "--bindings",
            str(fixture["bindings"]),
            "--input-gate",
            str(fixture["gate"]),
            "--resource-manifest",
            str(fixture["resources"]),
            "--output-directory",
            str(run_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert prepared.returncode == 0
    assert json.loads(prepared.stdout)["status"] == "PREPARED"
    assert PRIVATE_BINDING not in prepared.stdout
    assert str(fixture["workspace"]) not in prepared.stdout

    failed = subprocess.run(
        [
            sys.executable,
            "-m",
            module,
            "prepare",
            "--profile",
            str(fixture["profile"]),
            "--bindings",
            str(fixture["bindings"]),
            "--input-gate",
            str(fixture["gate"]),
            "--resource-manifest",
            str(fixture["resources"]),
            "--output-directory",
            str(run_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert failed.returncode == 2
    assert failed.stdout == ""
    assert "Traceback" not in failed.stderr
    assert PRIVATE_BINDING not in failed.stderr
    assert str(run_dir) not in failed.stderr

    blocked_parent = tmp_path / f"{PRIVATE_BINDING}-blocked-parent"
    blocked_parent.write_text("not a directory\n", encoding="utf-8")
    unexpected = subprocess.run(
        [
            sys.executable,
            "-m",
            module,
            "prepare",
            "--profile",
            str(fixture["profile"]),
            "--bindings",
            str(fixture["bindings"]),
            "--input-gate",
            str(fixture["gate"]),
            "--resource-manifest",
            str(fixture["resources"]),
            "--output-directory",
            str(blocked_parent / "run-plan"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert unexpected.returncode == 2
    assert unexpected.stdout == ""
    assert "Traceback" not in unexpected.stderr
    assert PRIVATE_BINDING not in unexpected.stderr
    assert str(blocked_parent) not in unexpected.stderr
