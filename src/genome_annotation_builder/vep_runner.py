"""Pinned offline Ensembl VEP execution and output verification."""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pysam

from .config import BuildConfig
from .errors import InputValidationError, ResourceError
from .predictor_resources import (
    resolve_predictor_paths,
    verify_predictor_resources,
)
from .util import atomic_output_directory, atomic_write_json, sha256_file
from .vep_parser import VepParser
from .vep_resources import verify_vep_cache


def verify_vep_runtime(config: BuildConfig) -> tuple[str, list[str]]:
    if config.vep_mode == "docker":
        command = [
            "docker",
            "run",
            "--rm",
            "--pull",
            "never",
            "--network",
            "none",
            config.vep_container_image,
            "vep",
            "--help",
        ]
    else:
        executable = shutil.which(config.vep_executable)
        if executable is None:
            raise ResourceError(f"VEP executable not found: {config.vep_executable}")
        command = [executable, "--help"]
    completed = subprocess.run(command, check=False, capture_output=True, text=True, timeout=120)
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        raise ResourceError(f"VEP version command failed: {output[-2000:]}")
    match = re.search(r"ensembl-vep\s*:\s*([^\s]+)", output)
    if not match:
        raise ResourceError("could not parse Ensembl VEP version")
    version = match.group(1)
    if version != config.vep_expected_version:
        raise ResourceError(
            f"VEP version mismatch: expected {config.vep_expected_version}, got {version}"
        )
    return version, command


def _vep_arguments(
    config: BuildConfig,
    input_path: str,
    output_path: str,
    fasta_path: str,
    cache_path: str,
    predictor_paths: dict[str, str] | None = None,
) -> list[str]:
    arguments = [
        "vep",
        "--input_file",
        input_path,
        "--output_file",
        output_path,
        "--format",
        "vcf",
        "--vcf",
        "--compress_output",
        "bgzip",
        "--force_overwrite",
        "--no_stats",
        "--offline",
        "--cache",
        "--dir_cache",
        cache_path,
        "--cache_version",
        config.vep_cache_version,
        "--species",
        "homo_sapiens",
        "--assembly",
        "GRCh38",
        "--fasta",
        fasta_path,
        "--check_ref",
        "--everything",
        "--mane",
        "--canonical",
        "--appris",
        "--tsl",
        "--ccds",
        "--hgvs",
        "--hgvsg",
        "--transcript_version",
        "--gene_version",
        "--symbol",
        "--numbers",
        "--domains",
        "--protein",
        "--biotype",
        "--allele_number",
        "--variant_class",
        "--sift",
        "b",
        "--polyphen",
        "b",
        "--af_gnomadg",
        "--af_gnomade",
        "--fork",
        str(config.vep_forks),
        "--buffer_size",
        str(config.vep_buffer_size),
    ]
    if not config.predictors.enabled:
        return arguments
    if predictor_paths is None:
        raise ResourceError("validated predictor paths are required for VEP plugin execution")
    arguments.extend(["--dir_plugins", config.vep_plugin_dir])

    def required(name: str) -> str:
        value = predictor_paths.get(name)
        if value is None:
            raise ResourceError(f"missing validated VEP predictor path: {name}")
        return value

    if config.predictors.loftee.enabled:
        arguments.extend(
            [
                "--plugin",
                "LoF,"
                f"loftee_path:{required('loftee_path')},"
                f"human_ancestor_fa:{required('loftee_human_ancestor')},"
                f"gerp_bigwig:{required('loftee_gerp_bigwig')},"
                f"conservation_file:{required('loftee_conservation_sqlite')}",
            ]
        )
    if config.predictors.cadd.enabled:
        arguments.extend(
            [
                "--plugin",
                f"CADD,snv={required('cadd_snv')},indels={required('cadd_indel')}",
            ]
        )
    if config.predictors.revel.enabled:
        arguments.extend(["--plugin", f"REVEL,file={required('revel_scores')}"])
    if config.predictors.alphamissense.enabled:
        arguments.extend(
            [
                "--plugin",
                f"AlphaMissense,file={required('alphamissense_scores')},transcript_match=1",
            ]
        )
    if config.predictors.spliceai.enabled:
        arguments.extend(
            [
                "--plugin",
                "SpliceAI,"
                f"snv={required('spliceai_snv')},"
                f"indel={required('spliceai_indel_sentinel')},split_output=1",
            ]
        )
    return arguments


def _predictor_path_strings(config: BuildConfig, *, docker: bool) -> dict[str, str]:
    local = resolve_predictor_paths(config, full_digest=False)
    if not docker:
        return {name: str(path.resolve()) for name, path in local.items()}
    root = config.predictors.root.resolve()
    return {
        name: "/predictors/" + path.resolve().relative_to(root).as_posix()
        for name, path in local.items()
    }


def _docker_perl5lib(config: BuildConfig, paths: dict[str, str]) -> str:
    base = "/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/modules:/plugins"
    if not config.predictors.loftee.enabled:
        return base
    return f"{paths['loftee_path']}:{base}"


def _command(config: BuildConfig, input_vcf: Path, staging: Path) -> list[str]:
    reference = config.reference_fasta.resolve()
    cache = config.vep_cache_dir.resolve()
    if config.vep_mode == "docker":
        predictor_paths = (
            _predictor_path_strings(config, docker=True) if config.predictors.enabled else None
        )
        internal_input = f"/input/{input_vcf.name}"
        internal_output = "/output/annotated.vcf.gz"
        internal_reference = f"/reference/{reference.name}"
        arguments = _vep_arguments(
            config,
            internal_input,
            internal_output,
            internal_reference,
            "/cache",
            predictor_paths,
        )
        docker_arguments = [
            "docker",
            "run",
            "--rm",
            "--pull",
            "never",
            "--network",
            "none",
            "--read-only",
            "--cap-drop",
            "ALL",
            "--security-opt",
            "no-new-privileges",
            "--tmpfs",
            "/tmp:rw,nosuid,noexec,size=4g",
            "--user",
            f"{os.getuid()}:{os.getgid()}",
            "--volume",
            f"{input_vcf.parent.resolve()}:/input:ro",
            "--volume",
            f"{reference.parent}:/reference:ro",
            "--volume",
            f"{cache}:/cache:ro",
        ]
        if config.predictors.enabled:
            assert predictor_paths is not None
            docker_arguments.extend(
                [
                    "--volume",
                    f"{config.predictors.root.resolve()}:/predictors:ro",
                    "--env",
                    f"PERL5LIB={_docker_perl5lib(config, predictor_paths)}",
                ]
            )
        docker_arguments.extend(
            [
                "--volume",
                f"{staging.resolve()}:/output:rw",
                "--workdir",
                "/output",
                config.vep_container_image,
                *arguments,
            ]
        )
        return docker_arguments
    executable = shutil.which(config.vep_executable)
    if executable is None:
        raise ResourceError(f"VEP executable not found: {config.vep_executable}")
    arguments = _vep_arguments(
        config,
        str(input_vcf.resolve()),
        str(staging / "annotated.vcf.gz"),
        str(reference),
        str(cache),
        (_predictor_path_strings(config, docker=False) if config.predictors.enabled else None),
    )
    arguments[0] = executable
    return arguments


def verify_predictor_plugin_runtime(config: BuildConfig) -> dict[str, Any]:
    """Confirm that every enabled plugin and LOFTEE dependency can be loaded."""

    if not config.predictors.enabled:
        return {"status": "not_enabled", "modules": []}
    modules = [
        module
        for plugin, module in (
            ("loftee", "LoF"),
            ("cadd", "CADD"),
            ("revel", "REVEL"),
            ("alphamissense", "AlphaMissense"),
            ("spliceai", "SpliceAI"),
        )
        if getattr(config.predictors, plugin).enabled
    ]
    ancestor_test_code = "1"
    ancestor_report: dict[str, Any] | None = None
    runtime_environment = os.environ.copy()
    if config.predictors.loftee.enabled:
        modules.append("Bio::DB::HTS::Faidx")
        ancestor_test_code = (
            "my $path = $ENV{'GAB_LOFTEE_ANCESTOR'}; "
            "die 'missing LOFTEE ancestor path' unless defined($path); "
            "my $faidx = Bio::DB::HTS::Faidx->new($path); "
            "my ($sequence) = $faidx->get_sequence('1:1-1'); "
            "die 'LOFTEE ancestor fetch failed' "
            "unless defined($sequence) && length($sequence) == 1;"
        )
        ancestor_report = {
            "status": "verified",
            "backend": "Bio::DB::HTS::Faidx",
            "probe": "1:1-1",
        }
    if config.vep_mode == "docker":
        paths = _predictor_path_strings(config, docker=True)
        docker_ancestor_environment = (
            ["--env", f"GAB_LOFTEE_ANCESTOR={paths['loftee_human_ancestor']}"]
            if config.predictors.loftee.enabled
            else []
        )
        command = [
            "docker",
            "run",
            "--rm",
            "--pull",
            "never",
            "--network",
            "none",
            "--read-only",
            "--cap-drop",
            "ALL",
            "--security-opt",
            "no-new-privileges",
            "--user",
            f"{os.getuid()}:{os.getgid()}",
            "--volume",
            f"{config.predictors.root.resolve()}:/predictors:ro",
            "--env",
            f"PERL5LIB={_docker_perl5lib(config, paths)}",
            *docker_ancestor_environment,
            "--entrypoint",
            "perl",
            config.vep_container_image,
            *(f"-M{module}" for module in modules),
            "-e",
            ancestor_test_code,
        ]
    else:
        perl = shutil.which("perl")
        if perl is None:
            raise ResourceError("Perl is required to load executable-mode VEP plugins")
        paths = _predictor_path_strings(config, docker=False)
        if config.predictors.loftee.enabled:
            runtime_environment["GAB_LOFTEE_ANCESTOR"] = paths["loftee_human_ancestor"]
        include_paths = [config.vep_plugin_dir]
        if config.predictors.loftee.enabled:
            include_paths.insert(0, paths["loftee_path"])
        command = [
            perl,
            *(value for path in include_paths for value in ("-I", path)),
            *(f"-M{module}" for module in modules),
            "-e",
            ancestor_test_code,
        ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=False,
        timeout=120,
        env=runtime_environment,
    )
    if completed.returncode != 0:
        stderr = (
            completed.stderr.decode("utf-8", errors="replace")
            if isinstance(completed.stderr, bytes)
            else str(completed.stderr or "")
        )
        stdout = (
            completed.stdout.decode("utf-8", errors="replace")
            if isinstance(completed.stdout, bytes)
            else str(completed.stdout or "")
        )
        raise ResourceError(
            "VEP predictor plugin runtime check failed: " + (stderr or stdout)[-4000:]
        )
    report = {
        "status": "verified",
        "modules": modules,
        "runtime": config.vep_container_image if config.vep_mode == "docker" else "perl",
    }
    if ancestor_report is not None:
        report["ancestral_lookup"] = ancestor_report
    return report


def _execution_environment(config: BuildConfig, staging: Path) -> dict[str, str]:
    environment = {**os.environ, "TMPDIR": str(staging)}
    if (
        config.vep_mode == "executable"
        and config.predictors.enabled
        and config.predictors.loftee.enabled
    ):
        paths = _predictor_path_strings(config, docker=False)
        existing = environment.get("PERL5LIB", "")
        components = [paths["loftee_path"], config.vep_plugin_dir]
        if existing:
            components.append(existing)
        environment["PERL5LIB"] = ":".join(components)
    return environment


def _vcf_identity(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    count = 0
    source = pysam.VariantFile(str(path))
    try:
        if list(source.header.samples) != ["SAMPLE"]:
            raise InputValidationError("VEP stage input/output must contain exactly SAMPLE")
        for record in source:
            sample = record.samples["SAMPLE"]
            gt = sample.get("GT")
            fields = (
                record.contig,
                str(record.pos),
                record.ref,
                ",".join(record.alts or ()),
                "/".join("." if allele is None else str(allele) for allele in (gt or ())),
                "1" if sample.phased else "0",
                str(record.info.get("GAB_SOURCE_ID", "")),
                str(record.info.get("GAB_SOURCE_DIGEST", "")),
                str(record.info.get("GAB_SOURCE_ORDINAL", "")),
                str(record.info.get("GAB_ALT_INDEX", "")),
            )
            digest.update("\0".join(fields).encode("utf-8"))
            digest.update(b"\n")
            count += 1
    finally:
        source.close()
    return {"record_count": count, "record_identity_sha256": digest.hexdigest()}


def _tail_text(path: Path, maximum_bytes: int = 8000) -> str:
    with path.open("rb") as handle:
        size = path.stat().st_size
        if size > maximum_bytes:
            handle.seek(-maximum_bytes, os.SEEK_END)
        return handle.read().decode("utf-8", errors="replace")


def _predictor_stderr_failures(config: BuildConfig, paths: list[Path]) -> list[str]:
    """Return fail-closed plugin warnings that VEP may otherwise treat as success."""

    if not config.predictors.enabled:
        return []
    failures: list[str] = []
    for path in paths:
        if not path.is_file():
            continue
        try:
            with path.open("r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    lowered = line.lower()
                    loftee_failure = config.predictors.loftee.enabled and (
                        "use of uninitialized value $faidx" in lowered
                        or "samtools faidx" in lowered
                        or ("warning" in lowered and "lof.pm" in lowered)
                        or ("warning" in lowered and "ancestr" in lowered)
                    )
                    if loftee_failure:
                        failures.append(line.rstrip()[:1000])
                        if len(failures) >= 20:
                            return failures
        except OSError as exc:
            raise InputValidationError(
                f"cannot inspect VEP predictor warning log {path}: {exc}"
            ) from exc
    return failures


def run_vep_annotation(
    config: BuildConfig,
    input_vcf: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    input_vcf = Path(input_vcf)
    destination = Path(destination)
    cache_manifest = verify_vep_cache(config)
    predictor_manifest = verify_predictor_resources(config)
    predictor_runtime = verify_predictor_plugin_runtime(config)
    version, _version_command = verify_vep_runtime(config)
    input_identity = _vcf_identity(input_vcf)
    started = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    with atomic_output_directory(destination) as staging:
        command = _command(config, input_vcf, staging)
        stdout_log = staging / "vep.stdout.log"
        stderr_log = staging / "vep.stderr.log"
        with stdout_log.open("wb") as stdout, stderr_log.open("wb") as stderr:
            completed = subprocess.run(
                command,
                check=False,
                stdout=stdout,
                stderr=stderr,
                env=_execution_environment(config, staging),
            )
        if completed.returncode != 0:
            detail = _tail_text(stderr_log) or _tail_text(stdout_log)
            raise InputValidationError(
                f"Ensembl VEP failed with exit {completed.returncode}: {detail}"
            )
        predictor_failures = _predictor_stderr_failures(
            config,
            [stderr_log, staging / "annotated.vcf.gz_warnings.txt"],
        )
        if predictor_failures:
            raise InputValidationError(
                "Ensembl VEP emitted a fail-closed predictor warning: "
                + " | ".join(predictor_failures)
            )
        annotated = staging / "annotated.vcf.gz"
        if not annotated.is_file():
            raise InputValidationError("Ensembl VEP did not create annotated.vcf.gz")
        index = subprocess.run(
            [config.bcftools, "index", "--tbi", "--force", str(annotated)],
            check=False,
            capture_output=True,
            text=True,
        )
        if index.returncode != 0:
            raise InputValidationError(f"failed to index VEP output: {index.stderr}")
        output_identity = _vcf_identity(annotated)
        if output_identity != input_identity:
            raise InputValidationError(
                "Ensembl VEP changed variant/genotype identity or record order: "
                f"input={input_identity}, output={output_identity}"
            )
        handle = pysam.VariantFile(str(annotated))
        try:
            parser = VepParser.from_header(handle.header)
            if not parser.enabled:
                raise InputValidationError("VEP output lacks a valid CSQ header")
            parser.require_pipeline_fields()
            parser.require_predictor_fields(config.predictors.enabled_plugins)
            if list(handle.header.samples) != ["SAMPLE"]:
                raise InputValidationError("VEP changed the projected sample columns")
        finally:
            handle.close()
        report = {
            **config.checkpoint_identity,
            "status": "PASS",
            "started_at": started,
            "completed_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
            "vep_version": version,
            "cache_version": config.vep_cache_version,
            "cache_tree_sha256": cache_manifest["tree_sha256"],
            "cache_metadata": cache_manifest["cache_metadata"],
            "container_image": (
                config.vep_container_image if config.vep_mode == "docker" else None
            ),
            "predictors": predictor_manifest,
            "predictor_runtime": predictor_runtime,
            **output_identity,
            "input_sha256": sha256_file(input_vcf),
            "outputs": {
                "annotated.vcf.gz": sha256_file(annotated),
                "annotated.vcf.gz.tbi": sha256_file(Path(f"{annotated}.tbi")),
                "vep.stdout.log": sha256_file(stdout_log),
                "vep.stderr.log": sha256_file(stderr_log),
            },
            "command_semantics": _vep_arguments(
                config,
                "INPUT.vcf.gz",
                "OUTPUT.vcf.gz",
                "GRCh38.fa",
                "VEP_CACHE",
                (
                    {
                        name: f"PREDICTORS/{name}"
                        for name in resolve_predictor_paths(config, full_digest=False)
                    }
                    if config.predictors.enabled
                    else None
                ),
            ),
            "stdout_tail": _tail_text(stdout_log, 4000),
            "stderr_tail": _tail_text(stderr_log, 4000),
        }
        atomic_write_json(staging / "vep_summary.json", report)
    return report
