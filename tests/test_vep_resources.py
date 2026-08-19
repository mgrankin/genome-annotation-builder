from __future__ import annotations

import io
import tarfile
from pathlib import Path

import pytest

from genome_annotation_builder.errors import ResourceError
from genome_annotation_builder.vep_resources import _bsd_sum, _validate_tar_members


def test_bsd_sum_matches_known_vector() -> None:
    # The value is independently produced by the POSIX/GNU `sum` utility.
    assert _bsd_sum(io.BytesIO(b"genome-annotation-builder\n")) == 50075


def _write_tar(path: Path, name: str, *, symlink: bool = False) -> None:
    payload = b"cache\n"
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo(name)
        if symlink:
            member.type = tarfile.SYMTYPE
            member.linkname = "/etc/passwd"
            archive.addfile(member)
        else:
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))


def test_cache_archive_member_validation(tmp_path: Path) -> None:
    valid = tmp_path / "valid.tar.gz"
    _write_tar(valid, "homo_sapiens/116_GRCh38/info.txt")
    assert _validate_tar_members(valid) == {"archive_files": 1, "archive_directories": 0}


@pytest.mark.parametrize(
    ("name", "symlink"),
    [
        ("homo_sapiens/../escape", False),
        ("wrong_root/cache.txt", False),
        ("homo_sapiens/116_GRCh38/link", True),
    ],
)
def test_cache_archive_rejects_unsafe_members(tmp_path: Path, name: str, symlink: bool) -> None:
    archive = tmp_path / "unsafe.tar.gz"
    _write_tar(archive, name, symlink=symlink)
    with pytest.raises(ResourceError):
        _validate_tar_members(archive)
