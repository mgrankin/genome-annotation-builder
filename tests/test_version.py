from __future__ import annotations

import tomllib
from pathlib import Path

import yaml

from genome_annotation_builder import PIPELINE_VERSION, __version__

ROOT = Path(__file__).resolve().parents[1]


def test_release_version_is_consistent_across_build_boundaries() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    lock = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))
    environment = yaml.safe_load((ROOT / "workflow/envs/gab.yaml").read_text(encoding="utf-8"))
    locked = {
        package["name"]: package["version"]
        for package in lock["package"]
        if package["name"] == "genome-annotation-builder"
    }

    assert __version__ == PIPELINE_VERSION
    assert project["project"]["version"] == PIPELINE_VERSION
    assert locked == {"genome-annotation-builder": PIPELINE_VERSION}
    assert environment["name"] == f"genome-annotation-builder-{PIPELINE_VERSION}"
    assert "setuptools>=77.0.3" in project["build-system"]["requires"]
