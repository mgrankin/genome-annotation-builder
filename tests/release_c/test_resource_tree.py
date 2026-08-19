from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from genome_annotation_builder.release_c.portable_model import PortableEvidenceError
from genome_annotation_builder.release_c.resource_tree import (
    create_resource_tree_manifest,
    validate_resource_tree_manifest,
)


def populate_tree(root: Path) -> None:
    (root / "empty").mkdir(parents=True)
    (root / "nested" / "deeper").mkdir(parents=True)
    (root / "catalog.json").write_text('{"loci": []}\n', encoding="utf-8")
    (root / "nested" / "deeper" / "index.bin").write_bytes(b"synthetic-index\n")


def test_manifest_is_path_free_deterministic_and_exact(tmp_path: Path) -> None:
    first_root = tmp_path / "first" / "resources"
    second_root = tmp_path / "second" / "resources"
    populate_tree(first_root)
    populate_tree(second_root)
    first_manifest = tmp_path / "first-tree.json"
    second_manifest = tmp_path / "second-tree.json"

    first = create_resource_tree_manifest(first_root, first_manifest)
    second = create_resource_tree_manifest(second_root, second_manifest)

    assert first_manifest.read_bytes() == second_manifest.read_bytes()
    assert first["manifest_sha256"] == second["manifest_sha256"]
    assert first["tree_identity_sha256"] == second["tree_identity_sha256"]
    assert first["files"] == 2
    assert first["directories"] == 3
    assert first["logical_bytes"] == 29
    payload = json.loads(first_manifest.read_text(encoding="utf-8"))
    assert payload["directories"] == ["empty", "nested", "nested/deeper"]
    assert [row["path"] for row in payload["files"]] == [
        "catalog.json",
        "nested/deeper/index.bin",
    ]
    retained = first_manifest.read_text(encoding="utf-8")
    assert str(first_root) not in retained
    assert str(second_root) not in retained

    validated = validate_resource_tree_manifest(
        first_root,
        first_manifest,
        expected_manifest_sha256=first["manifest_sha256"],
        expected_manifest_bytes=first["manifest_bytes"],
    )
    assert validated.portable_summary() == first
    assert set(validated.signatures) == {
        first_root,
        first_root / "empty",
        first_root / "nested",
        first_root / "nested" / "deeper",
        first_root / "catalog.json",
        first_root / "nested" / "deeper" / "index.bin",
        first_manifest,
    }


@pytest.mark.parametrize("mutation", ["modify", "add_file", "add_empty_directory"])
def test_validation_rejects_any_tree_inventory_or_content_change(
    tmp_path: Path, mutation: str
) -> None:
    root = tmp_path / "resources"
    populate_tree(root)
    manifest = tmp_path / "tree.json"
    created = create_resource_tree_manifest(root, manifest)

    if mutation == "modify":
        (root / "catalog.json").write_text('{"loci": ["changed"]}\n', encoding="utf-8")
    elif mutation == "add_file":
        (root / "undeclared.bin").write_bytes(b"undeclared\n")
    else:
        (root / "undeclared-empty").mkdir()

    with pytest.raises(PortableEvidenceError):
        validate_resource_tree_manifest(
            root,
            manifest,
            expected_manifest_sha256=created["manifest_sha256"],
            expected_manifest_bytes=created["manifest_bytes"],
        )


def test_creation_rejects_symbolic_links_and_special_files(tmp_path: Path) -> None:
    symlink_root = tmp_path / "symlink-resources"
    symlink_root.mkdir()
    (symlink_root / "real.txt").write_text("content\n", encoding="utf-8")
    (symlink_root / "alias.txt").symlink_to(symlink_root / "real.txt")
    with pytest.raises(PortableEvidenceError, match="symbolic link"):
        create_resource_tree_manifest(symlink_root, tmp_path / "symlink-tree.json")

    fifo_root = tmp_path / "fifo-resources"
    fifo_root.mkdir()
    (fifo_root / "real.txt").write_text("content\n", encoding="utf-8")
    os.mkfifo(fifo_root / "pipe")
    with pytest.raises(PortableEvidenceError, match="special file"):
        create_resource_tree_manifest(fifo_root, tmp_path / "fifo-tree.json")


def test_creation_rejects_private_relative_member_and_existing_output(tmp_path: Path) -> None:
    private_root = tmp_path / "private-resources"
    private_root.mkdir()
    (private_root / "privateid.txt").write_text("content\n", encoding="utf-8")
    with pytest.raises(PortableEvidenceError, match="forbidden private binding"):
        create_resource_tree_manifest(
            private_root,
            tmp_path / "private-tree.json",
            forbidden_bindings=("privateid",),
        )

    root = tmp_path / "resources"
    populate_tree(root)
    manifest = tmp_path / "tree.json"
    create_resource_tree_manifest(root, manifest)
    original = manifest.read_bytes()
    with pytest.raises(PortableEvidenceError, match="refusing to replace"):
        create_resource_tree_manifest(root, manifest)
    assert manifest.read_bytes() == original


def test_validation_rejects_manifest_tampering(tmp_path: Path) -> None:
    root = tmp_path / "resources"
    populate_tree(root)
    manifest = tmp_path / "tree.json"
    created = create_resource_tree_manifest(root, manifest)
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["tree_identity_sha256"] = "0" * 64
    manifest.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")

    with pytest.raises(PortableEvidenceError, match="manifest digest changed"):
        validate_resource_tree_manifest(
            root,
            manifest,
            expected_manifest_sha256=created["manifest_sha256"],
            expected_manifest_bytes=created["manifest_bytes"],
        )
    with pytest.raises(PortableEvidenceError, match="tree identity changed"):
        validate_resource_tree_manifest(root, manifest)


def test_cli_emits_only_portable_summary_and_sanitized_failure(tmp_path: Path) -> None:
    private_token = "privateid"
    root = tmp_path / private_token / "resources"
    populate_tree(root)
    manifest = tmp_path / "tree.json"
    module = "genome_annotation_builder.release_c.resource_tree"
    created = subprocess.run(
        [
            sys.executable,
            "-m",
            module,
            "create",
            "--root",
            str(root),
            "--output",
            str(manifest),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert created.returncode == 0
    assert json.loads(created.stdout)["status"] == "PASS"
    assert str(root) not in created.stdout
    assert private_token not in created.stdout

    (root / "extra.txt").write_text("changed\n", encoding="utf-8")
    failed = subprocess.run(
        [
            sys.executable,
            "-m",
            module,
            "validate",
            "--root",
            str(root),
            "--manifest",
            str(manifest),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert failed.returncode == 2
    assert failed.stdout == ""
    assert "Traceback" not in failed.stderr
    assert str(root) not in failed.stderr
    assert private_token not in failed.stderr

    blocked_parent = tmp_path / f"{private_token}-blocked-parent"
    blocked_parent.write_text("not a directory\n", encoding="utf-8")
    unexpected = subprocess.run(
        [
            sys.executable,
            "-m",
            module,
            "create",
            "--root",
            str(root),
            "--output",
            str(blocked_parent / "tree.json"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert unexpected.returncode == 2
    assert unexpected.stdout == ""
    assert "Traceback" not in unexpected.stderr
    assert str(blocked_parent) not in unexpected.stderr
    assert private_token not in unexpected.stderr
