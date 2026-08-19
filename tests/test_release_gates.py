from __future__ import annotations

import importlib
import io
import stat
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
try:
    privacy_scan = importlib.import_module("privacy_scan")
    distribution_gate = importlib.import_module("verify_distributions")
    reproducibility_gate = importlib.import_module("verify_reproducible_distributions")
finally:
    sys.path.remove(str(SCRIPT_DIR))
scan_payload = privacy_scan.scan_payload
REQUIRED_SUFFIXES = distribution_gate.REQUIRED_SUFFIXES
_validate = distribution_gate._validate
compare_distributions = reproducibility_gate.compare_distributions
ReproducibilityError = reproducibility_gate.ReproducibilityError

PUBLIC_CONFIG_EXAMPLES = (
    "config/config.example.yaml",
    "config/release-c-adapters.example.json",
    "config/release-c-matrix.example.json",
    "config/release-c-resource-metadata.example.json",
    "config/release-c-resources.example.json",
    "config/release-c.example.yaml",
)


def _wheel(path: Path, *, extra_name: str = "genome_annotation_builder/module.py", extra=b"ok"):
    with zipfile.ZipFile(path, "w") as archive:
        for suffix in REQUIRED_SUFFIXES:
            archive.writestr(f"genome_annotation_builder/{suffix}", b"pinned public fixture")
        archive.writestr(extra_name, extra)


def _reproducible_wheel(path: Path, payload: bytes) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        info = zipfile.ZipInfo("genome_annotation_builder/module.py", (2020, 1, 1, 0, 0, 0))
        info.external_attr = (stat.S_IFREG | 0o644) << 16
        archive.writestr(info, payload)


def _source_distribution(path: Path, payload: bytes) -> None:
    with tarfile.open(path, "w:gz") as archive:
        info = tarfile.TarInfo("genome_annotation_builder-1.0/src/module.py")
        info.size = len(payload)
        info.mode = 0o644
        archive.addfile(info, io.BytesIO(payload))


def test_privacy_scan_rejects_binary_payload() -> None:
    errors = scan_payload(relative="src/unexpected.bin", payload=b"a\0b", size=3, denied=[])
    assert errors == ["binary payload is not allowed in public source: src/unexpected.bin"]


def test_gitignore_allows_only_explicit_public_config_examples(tmp_path: Path) -> None:
    repository_root = Path(__file__).resolve().parents[1]
    (tmp_path / ".gitignore").write_bytes((repository_root / ".gitignore").read_bytes())
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)

    for relative in PUBLIC_CONFIG_EXAMPLES:
        completed = subprocess.run(
            ["git", "check-ignore", "--no-index", "-q", relative],
            cwd=tmp_path,
            check=False,
        )
        assert completed.returncode == 1, relative

    for relative in (
        "config/config.yaml",
        "config/production.json",
        "config/patient.json",
        "private/data.json",
        "output/sample.genome",
        "work/run.log",
        "resources/reference.fa",
    ):
        completed = subprocess.run(
            ["git", "check-ignore", "--no-index", "-q", relative],
            cwd=tmp_path,
            check=False,
        )
        assert completed.returncode == 0, relative


def test_distribution_gate_accepts_minimal_public_wheel(tmp_path: Path) -> None:
    wheel = tmp_path / "package.whl"
    _wheel(wheel)
    _validate(wheel, [])


def test_distribution_gate_requires_every_release_c_profile(tmp_path: Path) -> None:
    wheel = tmp_path / "package.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        for suffix in REQUIRED_SUFFIXES[:-1]:
            archive.writestr(
                f"genome_annotation_builder/{suffix}",
                b"pinned public fixture",
            )
    with pytest.raises(SystemExit, match="missing required members"):
        _validate(wheel, [])


def test_distribution_gate_rejects_local_deny_value(tmp_path: Path) -> None:
    wheel = tmp_path / "package.whl"
    _wheel(wheel, extra=b"deployment-only-subject")
    with pytest.raises(SystemExit, match="local deny-list entry"):
        _validate(wheel, [b"deployment-only-subject"])


def test_distribution_gate_rejects_wheel_symlink(tmp_path: Path) -> None:
    wheel = tmp_path / "package.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        link = zipfile.ZipInfo("genome_annotation_builder/link")
        link.create_system = 3
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(link, "module.py")
    with pytest.raises(SystemExit, match="symbolic link"):
        _validate(wheel, [])


def test_distribution_gate_allows_only_explicit_synthetic_vcf(tmp_path: Path) -> None:
    wheel = tmp_path / "package.whl"
    _wheel(
        wheel,
        extra_name="tests/fixtures/synthetic-positive.vcf",
        extra=b"##fileformat=VCFv4.3\n",
    )
    _validate(wheel, [])

    non_synthetic = tmp_path / "non-synthetic.whl"
    _wheel(
        non_synthetic,
        extra_name="tests/fixtures/subject.vcf",
        extra=b"##fileformat=VCFv4.3\n",
    )
    with pytest.raises(SystemExit, match="forbidden members"):
        _validate(non_synthetic, [])


def test_distribution_reproducibility_accepts_identical_payloads(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    for directory in (first, second):
        _reproducible_wheel(directory / "package.whl", b"same wheel")
        _source_distribution(directory / "package.tar.gz", b"same source")

    report = compare_distributions(first, second)
    assert report["status"] == "PASS"
    assert report["wheel"]["archive_byte_identical"] is True
    assert report["source_distribution"]["payload_identical"] is True


def test_distribution_reproducibility_rejects_payload_drift(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    _reproducible_wheel(first / "package.whl", b"same wheel")
    _reproducible_wheel(second / "package.whl", b"same wheel")
    _source_distribution(first / "package.tar.gz", b"first source")
    _source_distribution(second / "package.tar.gz", b"changed source")

    with pytest.raises(ReproducibilityError, match="source-distribution payload differs"):
        compare_distributions(first, second)
