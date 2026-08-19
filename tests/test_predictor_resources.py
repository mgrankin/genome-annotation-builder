from __future__ import annotations

import gzip
import hashlib
import io
import sqlite3
import tarfile
import zipfile
from pathlib import Path

import duckdb
import pyarrow as pa
import pysam
import pytest

import genome_annotation_builder.predictor_resources as resources
import genome_annotation_builder.vep_runner as vep_runner
from genome_annotation_builder.bundle import _predictor_semantic_control_stage
from genome_annotation_builder.config import load_config
from genome_annotation_builder.errors import InputValidationError, ResourceError
from genome_annotation_builder.predictor_doctor import (
    CONTROL_DEFINITION_SHA256,
    CONTROL_SET,
    CONTROL_VARIANTS,
    SEMANTIC_BOUNDARY,
    evaluate_predictor_controls,
)
from genome_annotation_builder.predictor_resources import RemoteArtifact
from genome_annotation_builder.util import canonical_json_bytes
from genome_annotation_builder.validator import (
    PREDICTOR_ARTIFACT_KEYS,
    _validate_predictor_doctor_metadata,
    _validate_predictor_metadata,
    _validate_predictor_values,
)
from genome_annotation_builder.vep_parser import (
    PIPELINE_REQUIRED_FIELDS,
    PREDICTOR_REQUIRED_FIELDS,
)
from genome_annotation_builder.vep_runner import _vep_arguments


def _annotated_control_vcf(path: Path, *, first_cadd: str = "21.0") -> None:
    fields = tuple(
        sorted(PIPELINE_REQUIRED_FIELDS | frozenset().union(*PREDICTOR_REQUIRED_FIELDS.values()))
    )
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=248956422>",
        "##contig=<ID=17,length=83257441>",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Synthetic VEP. Format: '
        + "|".join(fields)
        + '">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    for control in CONTROL_VARIANTS:
        values = {field: "" for field in fields}
        values.update(
            {
                "Allele": control["alt"],
                "Consequence": control["consequence"],
                "IMPACT": "HIGH" if control["consequence"] == "stop_gained" else "MODERATE",
                "SYMBOL": control["gene"],
                "Gene": "ENSG_SYNTHETIC",
                "Feature_type": "Transcript",
                "Feature": "ENST_SYNTHETIC",
                "BIOTYPE": "protein_coding",
                "CANONICAL": "YES",
                "MANE_SELECT": "NM_SYNTHETIC.1",
                "HGNC_ID": "1",
                "ALLELE_NUM": "1",
                "gnomADg_AF": "0.001",
                "gnomADe_AF": "0.001",
            }
        )
        checks = control["checks"]
        if "cadd" in checks:
            values["CADD_PHRED"] = (
                first_cadd if control["control_id"] == "or4f5_missense" else str(checks["cadd"])
            )
        if "loftee" in checks:
            values["LoF"] = str(checks["loftee"])
        if "revel" in checks:
            values["REVEL"] = str(checks["revel"])
        if "alphamissense" in checks:
            values["am_pathogenicity"] = str(checks["alphamissense"])
        if "spliceai" in checks:
            values["SpliceAI_pred_DS_AG"] = str(checks["spliceai"])
            values["SpliceAI_pred_DS_AL"] = "0"
            values["SpliceAI_pred_DS_DG"] = "0"
            values["SpliceAI_pred_DS_DL"] = "0"
        csq = "|".join(values[field] for field in fields)
        lines.append(
            "\t".join(
                [
                    control["chrom"],
                    str(control["pos"]),
                    control["control_id"],
                    control["ref"],
                    control["alt"],
                    "60",
                    "PASS",
                    "CSQ=" + csq,
                    "GT",
                    "0/1",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_predictor_doctor_checks_all_five_real_plugin_semantics(tmp_path: Path) -> None:
    annotated = tmp_path / "annotated.vcf"
    _annotated_control_vcf(annotated)
    enabled = ("loftee", "cadd", "revel", "alphamissense", "spliceai")
    results = evaluate_predictor_controls(annotated, enabled)
    observed = {
        (row["control_id"], plugin): details["observed"]
        for row in results
        for plugin, details in row["checks"].items()
    }
    assert observed[("or4f5_missense", "cadd")] == 21.0
    assert observed[("or4f5_missense", "spliceai")] == 0.03
    assert observed[("or4f5_stop_gained", "loftee")] == "HC"
    assert observed[("or4f5_stop_gained", "cadd")] == 25.4
    assert observed[("tp53_missense", "cadd")] == 23.4
    assert observed[("tp53_missense", "revel")] == 0.675
    assert observed[("tp53_missense", "alphamissense")] == 0.5749

    tampered = tmp_path / "tampered.vcf"
    _annotated_control_vcf(tampered, first_cadd="20.9")
    with pytest.raises(InputValidationError, match=r"expected 21\.0, observed 20\.9"):
        evaluate_predictor_controls(tampered, enabled)


def test_remote_artifact_requires_sha256_and_checks_both_hashes(tmp_path: Path) -> None:
    path = tmp_path / "resource.bin"
    payload = b"immutable predictor resource\n"
    path.write_bytes(payload)
    md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    spec = RemoteArtifact(
        key="synthetic",
        plugin="cadd",
        relative_path="resource.bin",
        url="https://example.invalid/resource.bin",
        bytes=len(payload),
        sha256=hashlib.sha256(payload).hexdigest(),
        md5=md5,
    )
    identity = resources._verify_remote(path, spec)
    assert identity["sha256"] == spec.sha256
    assert identity["upstream_md5"] == md5

    with pytest.raises(ResourceError, match="not yet been frozen"):
        resources._verify_remote(
            path,
            RemoteArtifact(
                key="unfrozen",
                plugin="cadd",
                relative_path="resource.bin",
                url=spec.url,
                bytes=len(payload),
                sha256=None,
            ),
        )

    path.write_bytes(payload + b"tampered")
    with pytest.raises(ResourceError, match="byte count mismatch"):
        resources._verify_remote(path, spec)


def test_all_published_predictor_resource_specs_are_frozen() -> None:
    """A release must never advertise a mutable or only partly identified object."""

    assert resources.REMOTE_SPECS
    assert set(resources.REMOTE_SPECS) == {spec.key for spec in resources.REMOTE_SPECS.values()}
    relative_paths: set[str] = set()
    for key, spec in resources.REMOTE_SPECS.items():
        assert key == spec.key
        assert spec.plugin in resources.PLUGIN_ORDER
        assert spec.url.startswith("https://")
        assert spec.bytes > 0
        assert spec.sha256 is not None
        assert resources.SHA256_PATTERN.fullmatch(spec.sha256)
        assert not Path(spec.relative_path).is_absolute()
        assert ".." not in Path(spec.relative_path).parts
        assert spec.relative_path not in relative_paths
        relative_paths.add(spec.relative_path)
        if spec.md5 is not None:
            assert resources.MD5_PATTERN.fullmatch(spec.md5)


def test_predictor_manifest_schema_is_byte_stable_without_install_time(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
project: {output_dir: output, work_dir: work, bundle_name: synthetic.genome}
input: {vcf: synthetic.vcf.gz, sample: SAMPLE, genome_build: GRCh38}
reference: {fasta: synthetic.fa, contig_style: chr}
resources:
  cache_dir: resources
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
annotation:
  vep: {enabled: true, mode: docker, plugin_dir: /plugins}
  predictors:
    root: private-predictors
    cadd:
      enabled: true
      license: "CADD v1.7 terms; non-commercial use unless separately licensed"
      redistribution_allowed: false
""",
        encoding="utf-8",
    )
    config = load_config(config_path)
    artifacts = {
        key: {"plugin": "cadd", "sha256": str(index) * 64}
        for index, key in enumerate(sorted(resources._expected_artifact_keys("cadd")), start=1)
    }
    first = resources._manifest_document(config, {"cadd"}, artifacts, {}, {"cadd": {}})
    second = resources._manifest_document(config, {"cadd"}, artifacts, {}, {"cadd": {}})

    assert first == second
    assert first["schema_version"] == 3
    assert "created_at" not in first
    assert canonical_json_bytes(first) == canonical_json_bytes(second)


def _loftee_tar(path: Path, *, unsafe: bool = False) -> None:
    root = f"loftee-{resources.LOFTEE_COMMIT}"
    required = (
        "LoF.pm",
        "LICENSE",
        "maxEntScan/score3.pl",
        "maxEntScan/splicemodels/me2x5",
        "splice_data/donor_model.txt",
        "splice_data/acceptor_model.txt",
    )
    with tarfile.open(path, "w:gz") as archive:
        for name in required:
            payload = f"synthetic {name}\n".encode()
            member = tarfile.TarInfo(f"{root}/{name}")
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
        if unsafe:
            member = tarfile.TarInfo(f"{root}/escape")
            member.type = tarfile.SYMTYPE
            member.linkname = "/etc/passwd"
            archive.addfile(member)


def test_loftee_source_extraction_rejects_links_and_is_digestible(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    synthetic_identity = {
        "source_sha256": "a" * 64,
        "source_archive_sha256": "a" * 64,
        "patch_id": "synthetic",
        "patch_file": "synthetic.patch",
        "patch_sha256": "b" * 64,
        "original_lof_pm_sha256": "c" * 64,
        "patched_lof_pm_sha256": "d" * 64,
        "patch_command_semantics": ["synthetic"],
    }
    monkeypatch.setattr(resources, "_apply_loftee_patch", lambda _path: synthetic_identity)
    monkeypatch.setattr(resources, "_loftee_patch_identity", lambda _path: synthetic_identity)
    valid = tmp_path / "valid.tar.gz"
    _loftee_tar(valid)
    relative, artifact = resources._extract_loftee_source(tmp_path / "valid-root", valid)
    assert relative.endswith(f"/{resources.LOFTEE_PATCHED_SOURCE_DIR}")
    assert artifact["kind"] == "source_tree"
    assert artifact["tree_files"] == 6
    assert artifact["patch_id"] == "synthetic"

    unsafe = tmp_path / "unsafe.tar.gz"
    _loftee_tar(unsafe, unsafe=True)
    with pytest.raises(ResourceError, match="unsupported LOFTEE archive member type"):
        resources._extract_loftee_source(tmp_path / "unsafe-root", unsafe)


def test_loftee_patch_is_exact_preimage_to_exact_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    source.mkdir()
    original = b"old lookup\n"
    patched = b"use Bio::DB::HTS::Faidx;\nmy $faidx = Bio::DB::HTS::Faidx->new($path);\n"
    (source / "LoF.pm").write_bytes(original)
    patch_path = tmp_path / "synthetic.patch"
    patch_path.write_text(
        "--- a/LoF.pm\n"
        "+++ b/LoF.pm\n"
        "@@ -1 +1,2 @@\n"
        "-old lookup\n"
        "+use Bio::DB::HTS::Faidx;\n"
        "+my $faidx = Bio::DB::HTS::Faidx->new($path);\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(resources, "LOFTEE_SOURCE_SHA256", "a" * 64)
    monkeypatch.setattr(
        resources, "LOFTEE_LOF_PM_ORIGINAL_SHA256", hashlib.sha256(original).hexdigest()
    )
    monkeypatch.setattr(
        resources, "LOFTEE_LOF_PM_PATCHED_SHA256", hashlib.sha256(patched).hexdigest()
    )
    monkeypatch.setattr(resources, "LOFTEE_PATCH_ID", "synthetic-faidx")
    monkeypatch.setattr(resources, "LOFTEE_PATCH_FILE", patch_path.name)
    monkeypatch.setattr(resources, "LOFTEE_PATCH_SHA256", resources.sha256_file(patch_path))
    monkeypatch.setattr(resources, "_loftee_patch_path", lambda: patch_path)

    identity = resources._apply_loftee_patch(source)
    assert (source / "LoF.pm").read_bytes() == patched
    assert identity["patch_id"] == "synthetic-faidx"
    assert identity["patched_lof_pm_sha256"] == hashlib.sha256(patched).hexdigest()
    assert not list(source.glob("*.rej"))
    assert not list(source.glob("*.orig"))


def test_shipped_loftee_patch_is_frozen_and_reviewable() -> None:
    patch_path = resources._loftee_patch_path()
    text = patch_path.read_text(encoding="utf-8")
    assert resources.sha256_file(patch_path) == resources.LOFTEE_PATCH_SHA256
    assert "Bio::DB::HTS::Faidx->new" in text
    assert "-    my $faidx = `samtools faidx" in text


def test_loftee_sqlite_derivative_is_validated_and_source_bound(tmp_path: Path) -> None:
    plain = tmp_path / "source.sqlite"
    connection = sqlite3.connect(plain)
    connection.execute(
        "CREATE TABLE phylocsf_summary "
        "(chrom TEXT, transcript TEXT, exon INTEGER, corresponding_orf_score REAL)"
    )
    connection.execute("INSERT INTO phylocsf_summary VALUES ('1', 'ENST000001', 1, 3.5)")
    connection.commit()
    connection.close()
    compressed = tmp_path / "source.sqlite.gz"
    with (
        plain.open("rb") as source,
        compressed.open("wb") as compressed_handle,
        gzip.GzipFile(fileobj=compressed_handle, mode="wb", mtime=0) as output,
    ):
        output.write(source.read())

    relative, identity = resources._prepare_loftee_sql(tmp_path / "root", compressed)
    prepared = tmp_path / "root" / relative
    assert identity["phylocsf_summary_rows"] == 1
    assert resources._sqlite_identity(prepared)["sha256"] == identity["sha256"]


def test_revel_preparation_filters_missing_grch38_and_sorts_for_tabix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    header = b",".join(resources.REVEL_HEADER) + b"\n"
    rows = (
        b"2,20,20,A,G,K,R,0.9,ENST2\n"
        b"1,30,30,C,T,P,L,0.2,ENST1\n"
        b"1,31,.,C,G,P,R,0.3,ENST1\n"
        b"1,10,10,G,A,A,T,0.1,ENST1\n"
        b"Y,40,40,G,T,A,S,0.8,ENSTY\n"
    )
    member_payload = header + rows
    archive_path = tmp_path / "revel.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(resources.REVEL_MEMBER, member_payload)
    source_sha256 = resources.sha256_file(archive_path)
    monkeypatch.setattr(resources, "REVEL_MEMBER_BYTES", len(member_payload))
    monkeypatch.setattr(resources, "REVEL_ARCHIVE_SHA256", source_sha256)

    paths, report = resources._prepare_revel(tmp_path / "root", archive_path)
    assert report["source_rows"] == 5
    assert report["retained_grch38_rows"] == 4
    assert report["omitted_missing_grch38_position"] == 1
    tabix = pysam.TabixFile(
        str(tmp_path / "root" / paths["scores"]),
        index=str(tmp_path / "root" / paths["index"]),
    )
    try:
        chr1 = list(tabix.fetch("1"))
    finally:
        tabix.close()
    assert [int(row.split("\t")[2]) for row in chr1] == [10, 30]
    tabix = pysam.TabixFile(
        str(tmp_path / "root" / paths["scores"]),
        index=str(tmp_path / "root" / paths["index"]),
    )
    try:
        chromosome_y = list(tabix.fetch("Y"))
    finally:
        tabix.close()
    assert [int(row.split("\t")[2]) for row in chromosome_y] == [40]


def test_spliceai_empty_indel_sentinel_is_explicit_and_reusable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(resources, "SPLICEAI_SNV_SHA256", "a" * 64)
    paths, report = resources._prepare_spliceai_sentinel(tmp_path)
    assert report["record_count"] == 0
    identity = resources._tabix_identity(
        tmp_path / paths["indel_sentinel"],
        tmp_path / paths["indel_sentinel_index"],
        required_header_tokens=("EMPTY_SENTINEL_NO_SPLICEAI_INDEL_SCORES",),
        allow_empty=True,
    )
    assert identity["first_record_present"] is False
    assert resources._prepare_spliceai_sentinel(tmp_path)[1] == report


def test_vep_plugin_arguments_are_explicit_and_do_not_hide_spliceai_scope(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
project:
  output_dir: output
  work_dir: work
  bundle_name: synthetic.genome
input:
  vcf: synthetic.vcf.gz
  sample: SAMPLE
  genome_build: GRCh38
reference:
  fasta: synthetic.fa
  contig_style: chr
resources:
  cache_dir: resources
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
annotation:
  vep:
    enabled: true
    mode: docker
    plugin_dir: /plugins
  predictors:
    root: private-predictors
    cadd:
      enabled: true
      license: "CADD v1.7 terms; non-commercial use unless separately licensed"
      redistribution_allowed: false
    spliceai:
      enabled: true
      license: "Ensembl terms plus upstream SpliceAI terms; redistribution review required"
      redistribution_allowed: false
""",
        encoding="utf-8",
    )
    config = load_config(config_path)
    arguments = _vep_arguments(
        config,
        "INPUT",
        "OUTPUT",
        "REFERENCE",
        "CACHE",
        {
            "cadd_snv": "PREDICTORS/cadd.snv.gz",
            "cadd_indel": "PREDICTORS/cadd.indel.gz",
            "spliceai_snv": "PREDICTORS/spliceai.snv.vcf.gz",
            "spliceai_indel_sentinel": "PREDICTORS/NO_INDEL_SCORES.vcf.gz",
        },
    )
    plugins = [
        arguments[index + 1] for index, value in enumerate(arguments[:-1]) if value == "--plugin"
    ]
    assert plugins == [
        "CADD,snv=PREDICTORS/cadd.snv.gz,indels=PREDICTORS/cadd.indel.gz",
        "SpliceAI,snv=PREDICTORS/spliceai.snv.vcf.gz,"
        "indel=PREDICTORS/NO_INDEL_SCORES.vcf.gz,split_output=1",
    ]


def test_predictor_runtime_uses_the_host_identity_for_private_resource_access(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
project: {output_dir: output, work_dir: work, bundle_name: synthetic.genome}
input: {vcf: synthetic.vcf.gz, sample: SAMPLE, genome_build: GRCh38}
reference: {fasta: synthetic.fa, contig_style: chr}
resources:
  cache_dir: resources
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
annotation:
  vep:
    enabled: true
    mode: docker
    plugin_dir: /plugins
  predictors:
    root: private-predictors
    cadd:
      enabled: true
      license: "CADD v1.7 terms; non-commercial use unless separately licensed"
      redistribution_allowed: false
""",
        encoding="utf-8",
    )
    config = load_config(config_path)
    captured: dict[str, object] = {}

    def completed(command: list[str], **_kwargs: object) -> object:
        captured["command"] = command
        return type("Completed", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(vep_runner, "resolve_predictor_paths", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(vep_runner.subprocess, "run", completed)
    monkeypatch.setattr(vep_runner.os, "getuid", lambda: 1234)
    monkeypatch.setattr(vep_runner.os, "getgid", lambda: 5678)
    report = vep_runner.verify_predictor_plugin_runtime(config)
    command = captured["command"]
    assert isinstance(command, list)
    assert command[command.index("--user") + 1] == "1234:5678"
    assert report["modules"] == ["CADD"]


def test_loftee_runtime_probes_ancestor_with_biodbhts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
project: {output_dir: output, work_dir: work, bundle_name: synthetic.genome}
input: {vcf: synthetic.vcf.gz, sample: SAMPLE, genome_build: GRCh38}
reference: {fasta: synthetic.fa, contig_style: chr}
resources:
  cache_dir: resources
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
annotation:
  vep:
    enabled: true
    mode: docker
    plugin_dir: /plugins
  predictors:
    root: private-predictors
    loftee:
      enabled: true
      redistribution_allowed: false
""",
        encoding="utf-8",
    )
    config = load_config(config_path)
    captured: dict[str, object] = {}
    loftee_path = config.predictors.root / "prepared" / "loftee"
    ancestor_path = config.predictors.root / "archives" / "ancestor.fa.gz"

    def completed(command: list[str], **_kwargs: object) -> object:
        captured["command"] = command
        return type("Completed", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(
        vep_runner,
        "resolve_predictor_paths",
        lambda *_args, **_kwargs: {
            "loftee_path": loftee_path,
            "loftee_human_ancestor": ancestor_path,
        },
    )
    monkeypatch.setattr(vep_runner.subprocess, "run", completed)
    report = vep_runner.verify_predictor_plugin_runtime(config)
    command = captured["command"]
    assert isinstance(command, list)
    assert "-MBio::DB::HTS::Faidx" in command
    assert "GAB_LOFTEE_ANCESTOR=/predictors/archives/ancestor.fa.gz" in command
    assert any("get_sequence('1:1-1')" in item for item in command)
    assert report["ancestral_lookup"] == {
        "status": "verified",
        "backend": "Bio::DB::HTS::Faidx",
        "probe": "1:1-1",
    }


def test_predictor_warning_gate_rejects_silent_loftee_faidx_failure(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
project: {output_dir: output, work_dir: work, bundle_name: synthetic.genome}
input: {vcf: synthetic.vcf.gz, sample: SAMPLE, genome_build: GRCh38}
reference: {fasta: synthetic.fa, contig_style: chr}
resources:
  cache_dir: resources
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
annotation:
  vep:
    enabled: true
    mode: docker
    plugin_dir: /plugins
  predictors:
    root: private-predictors
    loftee:
      enabled: true
      redistribution_allowed: false
""",
        encoding="utf-8",
    )
    config = load_config(config_path)
    warning = tmp_path / "vep.stderr.log"
    warning.write_text(
        "WARNING: 20 : Use of uninitialized value $faidx in split at LoF.pm line 540\n",
        encoding="utf-8",
    )
    failures = vep_runner._predictor_stderr_failures(config, [warning])
    assert len(failures) == 1
    assert "$faidx" in failures[0]


def _portable_predictor_documents() -> dict[str, dict[str, object]]:
    manifest_digest = "a" * 64

    def integrity(plugin: str) -> str:
        artifacts = ",".join(
            f"{key}=sha256:{'b' * 64}" for key in sorted(PREDICTOR_ARTIFACT_KEYS[plugin])
        )
        return f"manifest-sha256:{manifest_digest};artifacts:{artifacts}"

    cadd_coverage = {
        "variant_types": ["SNV", "small indel"],
        "assembly": "GRCh38",
        "pre_scored_only": True,
    }
    spliceai_coverage = {
        "variant_types": ["SNV"],
        "assembly": "GRCh38",
        "scope": "Ensembl MANE v1.4 masked SNV scores, GRCh38",
        "indel_scores_available": False,
        "indel_runtime_input": "explicit empty sentinel",
    }
    return {
        "raw_artifacts.json": {
            "artifacts": [
                {
                    "artifact_type": "offline VEP predictor resource: cadd",
                    "disposition": "referenced",
                    "integrity_identifier": integrity("cadd"),
                    "license": "synthetic CADD terms",
                    "redistribution_allowed": False,
                    "coverage": cadd_coverage,
                },
                {
                    "artifact_type": "offline VEP predictor resource: spliceai",
                    "disposition": "referenced",
                    "integrity_identifier": integrity("spliceai"),
                    "license": "synthetic SpliceAI terms",
                    "redistribution_allowed": False,
                    "coverage": spliceai_coverage,
                },
            ]
        },
        "transformations.json": {
            "stages": [
                {
                    "stage_id": "ensembl_vep_annotation",
                    "enabled_predictors": ["cadd", "spliceai"],
                    "predictor_coverage": {
                        "cadd": cadd_coverage,
                        "spliceai": spliceai_coverage,
                    },
                    "inputs": ["normalized alleles", f"sha256:{manifest_digest}"],
                    "command_semantics": [
                        "vep",
                        "--dir_plugins",
                        "/plugins",
                        "--plugin",
                        "CADD,snv=PREDICTORS/cadd_snv,indels=PREDICTORS/cadd_indel",
                        "--plugin",
                        "SpliceAI,snv=PREDICTORS/spliceai_snv,"
                        "indel=PREDICTORS/spliceai_indel_sentinel,split_output=1",
                    ],
                }
            ]
        },
        "omissions.json": {
            "omissions": [
                {"scope": "LOFTEE predictor annotation"},
                {"scope": "REVEL predictor annotation"},
                {"scope": "AlphaMissense predictor annotation"},
                {"scope": "SpliceAI indel predictions"},
            ]
        },
    }


def test_portable_predictor_metadata_has_one_cross_checked_boundary() -> None:
    errors: list[str] = []
    metrics = _validate_predictor_metadata(_portable_predictor_documents(), errors)
    assert errors == []
    assert metrics == {
        "status": "validated",
        "enabled_plugins": ["cadd", "spliceai"],
        "resource_artifacts": 2,
        "resource_manifest_count": 1,
    }


def test_portable_predictor_doctor_metadata_is_cross_checked() -> None:
    documents = _portable_predictor_documents()
    enabled = ["cadd", "spliceai"]
    results = []
    for control in CONTROL_VARIANTS:
        checks = {
            plugin: {"expected": value, "observed": value, "status": "PASS"}
            for plugin, value in control["checks"].items()
            if plugin in enabled
        }
        results.append({"control_id": control["control_id"], "checks": checks})
    controls = {
        "control_set": CONTROL_SET,
        "control_definition_sha256": CONTROL_DEFINITION_SHA256,
        "doctor_run_key": "e" * 64,
        "report_sha256": "c" * 64,
        "record_count": len(CONTROL_VARIANTS),
        "record_identity_sha256": "d" * 64,
        "reference_fasta_sha256": "f" * 64,
        "vep_cache_manifest_sha256": "1" * 64,
        "control_results": results,
        "semantic_boundary": SEMANTIC_BOUNDARY,
    }
    documents["raw_artifacts.json"]["artifacts"].append(
        {
            "artifact_type": "offline VEP predictor semantic control report",
            "integrity_identifier": "sha256:" + controls["report_sha256"],
            "control_set": CONTROL_SET,
            "control_definition_sha256": CONTROL_DEFINITION_SHA256,
            "doctor_run_key": controls["doctor_run_key"],
            "enabled_plugins": enabled,
            "record_identity_sha256": controls["record_identity_sha256"],
            "reference_fasta_sha256": controls["reference_fasta_sha256"],
            "vep_cache_manifest_sha256": controls["vep_cache_manifest_sha256"],
        }
    )
    documents["transformations.json"]["stages"].append(
        {
            "stage_id": "offline_predictor_semantic_controls",
            **controls,
            "enabled_plugins": enabled,
            "status": "completed",
        }
    )
    manifest = {
        "resource_provenance": {
            "reference_fasta_sha256": controls["reference_fasta_sha256"],
            "offline_predictors": {
                "manifest_sha256": "a" * 64,
                "enabled_plugins": enabled,
                "semantic_controls": controls,
            },
        }
    }
    errors: list[str] = []
    metrics = _validate_predictor_doctor_metadata(documents, manifest, enabled, errors)
    assert errors == []
    assert metrics["status"] == "validated"
    assert metrics["control_count"] == 3
    assert metrics["covered_plugins"] == ["cadd", "spliceai"]

    controls["control_results"][0]["checks"]["cadd"]["observed"] = 20.9
    errors = []
    _validate_predictor_doctor_metadata(documents, manifest, enabled, errors)
    assert any("value mismatch" in error for error in errors)


def test_predictor_transformation_stage_retains_control_definition_digest() -> None:
    controls = {
        "control_set": CONTROL_SET,
        "control_definition_sha256": CONTROL_DEFINITION_SHA256,
        "doctor_run_key": "e" * 64,
        "report_sha256": "c" * 64,
        "enabled_plugins": ["cadd"],
        "record_count": len(CONTROL_VARIANTS),
        "record_identity_sha256": "d" * 64,
        "reference_fasta_sha256": "f" * 64,
        "vep_cache_manifest_sha256": "1" * 64,
        "control_results": [],
        "semantic_boundary": SEMANTIC_BOUNDARY,
    }
    stage = _predictor_semantic_control_stage(
        {"manifest_sha256": "a" * 64, "semantic_controls": controls},
        "b" * 64,
    )
    assert stage["control_definition_sha256"] == CONTROL_DEFINITION_SHA256
    assert stage["inputs"] == [
        "sha256:" + "a" * 64,
        "sha256:" + CONTROL_DEFINITION_SHA256,
    ]


def test_portable_predictor_metadata_rejects_tampered_scope_command_and_digest() -> None:
    documents = _portable_predictor_documents()
    raw_rows = documents["raw_artifacts.json"]["artifacts"]
    assert isinstance(raw_rows, list)
    raw_rows[0]["redistribution_allowed"] = True
    raw_rows[0]["integrity_identifier"] = str(raw_rows[0]["integrity_identifier"]).replace(
        "cadd_snv=sha256:", "missing_cadd_snv=sha256:"
    )
    raw_rows[1]["coverage"]["indel_scores_available"] = True
    stage = documents["transformations.json"]["stages"][0]
    stage["command_semantics"][-1] = str(stage["command_semantics"][-1]).replace(
        "spliceai_indel_sentinel", "unfrozen_indel_scores"
    )
    omissions = documents["omissions.json"]["omissions"]
    omissions.pop()

    errors: list[str] = []
    _validate_predictor_metadata(documents, errors)
    error_text = "\n".join(errors)
    assert "must not assert redistribution rights" in error_text
    assert "checksum coverage mismatch" in error_text
    assert "SNV-only boundary" in error_text
    assert "plugin arguments do not preserve resource semantics" in error_text
    assert "SpliceAI indel unavailability omission count" in error_text


def test_predictor_values_require_enabled_adapters_and_valid_score_domains() -> None:
    table = pa.Table.from_pylist(
        [
            {
                "predictions": {
                    "cadd_phred": -1.0,
                    "revel": 1.1,
                    "alphamissense": 0.5,
                    "spliceai_max": 0.2,
                },
                "consequence": {"all_annotations": [{"lof": "BAD"}]},
            }
        ]
    )
    connection = duckdb.connect(":memory:")
    try:
        connection.register("variants", table)
        errors: list[str] = []
        metrics = _validate_predictor_values(connection, set(), errors)
    finally:
        connection.close()
    assert metrics["cadd_values"] == 1
    assert metrics["invalid_cadd_values"] == 1
    assert metrics["invalid_revel_values"] == 1
    assert metrics["invalid_loftee_values"] == 1
    assert any("without an enabled adapter" in error for error in errors)
