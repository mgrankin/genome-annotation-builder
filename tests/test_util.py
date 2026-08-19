from __future__ import annotations

from pathlib import Path

import pytest

from genome_annotation_builder.util import source_tree_digest


def test_source_tree_digest_ignores_only_runtime_bytecode(tmp_path: Path) -> None:
    package = tmp_path / "package"
    package.mkdir()
    source = package / "module.py"
    source.write_text("VALUE = 1\n", encoding="utf-8")
    original = source_tree_digest(package)

    cache = package / "__pycache__"
    cache.mkdir()
    (cache / "module.cpython-311.pyc").write_bytes(b"runtime-only")
    (package / "orphan.pyc").write_bytes(b"runtime-only")
    assert source_tree_digest(package) == original

    source.write_text("VALUE = 2\n", encoding="utf-8")
    assert source_tree_digest(package) != original


def test_source_tree_digest_includes_non_python_vendor_patches(tmp_path: Path) -> None:
    package = tmp_path / "package"
    patch_directory = package / "vendor_patches"
    patch_directory.mkdir(parents=True)
    patch = patch_directory / "fix.patch"
    patch.write_text("first patch bytes\n", encoding="utf-8")
    original = source_tree_digest(package)
    patch.write_text("different patch bytes\n", encoding="utf-8")
    assert source_tree_digest(package) != original


def test_source_tree_digest_rejects_symlinked_implementation(tmp_path: Path) -> None:
    target = tmp_path / "real.py"
    target.write_text("VALUE = 1\n", encoding="utf-8")
    package = tmp_path / "package"
    package.mkdir()
    (package / "linked.py").symlink_to(target)
    with pytest.raises(ValueError, match="symbolic link"):
        source_tree_digest(package)
