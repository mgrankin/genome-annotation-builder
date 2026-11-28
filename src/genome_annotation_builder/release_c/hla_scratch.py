"""Retire the pinned HLA*LA caller's two Picard mate files in fresh matrix attempts."""

import logging
import os
import stat
from contextlib import ExitStack
from pathlib import Path
from typing import Any

Snapshot = tuple[Path, int, dict[str, os.stat_result]]
MATES = ("R_1.fastq", "R_2.fastq")


def new_workspace(path: Path) -> os.stat_result | None:
    """Claim a new matrix workspace using the same creation mode as the bindings builder."""
    if os.path.realpath(path) != os.fspath(path):
        return None
    try:
        path.mkdir(parents=True)
        return path.lstat()
    except OSError:
        return None


def _directory(path: str | Path, handles: ExitStack, *, parent: int | None = None) -> int:
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent)
    handles.callback(os.close, fd)
    return fd


def snapshot_mates(
    workspace: Path, identity: os.stat_result, artifact: dict[str, Any], handles: ExitStack
) -> Snapshot | None:
    """Open produced files immediately after a new native execution. Never follow a link."""
    try:
        if os.path.realpath(workspace) != os.fspath(workspace) or artifact.get("components") != {
            "bestguess_g_group": str(workspace / "SUBJECT/hla/R1_bestguess_G.txt")
        }:
            return None
        workspace_fd = _directory(workspace, handles)
        if not os.path.samestat(os.fstat(workspace_fd), identity):
            return None
        subject = _directory("SUBJECT", handles, parent=workspace_fd)
        source = os.stat("extraction.bam", dir_fd=subject, follow_symlinks=False)
        if (
            not stat.S_ISREG(source.st_mode)
            or source.st_size == 0
            or source.st_dev != identity.st_dev
            or source.st_dev != os.fstat(subject).st_dev
        ):
            return None
        mates = {}
        for name in MATES:
            fd = os.open(name, os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW, dir_fd=subject)
            handles.callback(os.close, fd)
            st = os.fstat(fd)
            if not stat.S_ISREG(st.st_mode) or st.st_nlink != 1 or st.st_dev != source.st_dev:
                return None
            mates[name] = st
    except OSError:
        return None
    return workspace / "SUBJECT", subject, {"extraction.bam": source, **mates}


def remove_validated_mates(snapshot: Snapshot | None) -> None:
    if snapshot is None:
        return
    path, subject, files = snapshot
    try:
        if os.path.realpath(path) != os.fspath(path) or not os.path.samestat(
            path.stat(), os.fstat(subject)
        ):
            return
        for name, initial in files.items():
            observed = os.stat(name, dir_fd=subject, follow_symlinks=False)
            if not os.path.samestat(initial, observed) or (
                initial.st_size,
                initial.st_ctime_ns,
                initial.st_nlink,
            ) != (observed.st_size, observed.st_ctime_ns, observed.st_nlink):
                return
        for name in MATES:
            os.unlink(name, dir_fd=subject)
    except OSError:
        logging.getLogger(__name__).warning("Could not retire verified HLA*LA mate scratch")
