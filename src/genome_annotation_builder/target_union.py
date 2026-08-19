"""Deterministic, disk-backed union of checksum-bound target catalogues."""

from __future__ import annotations

import csv
import sqlite3
from collections.abc import Sequence
from pathlib import Path

from .errors import ResourceError
from .target_callability import CATALOG_COLUMNS


def merge_target_catalogs(
    sources: Sequence[tuple[str, Path]], destination: str | Path
) -> dict[str, int]:
    """Merge catalogues in caller-defined order while enforcing global target IDs."""

    destination = Path(destination)
    database_path = destination.with_suffix(destination.suffix + ".union.sqlite")
    database = sqlite3.connect(str(database_path))
    database.execute("PRAGMA journal_mode=OFF")
    database.execute("PRAGMA synchronous=OFF")
    database.execute("CREATE TABLE targets (target_id TEXT PRIMARY KEY)")
    counts: dict[str, int] = {}
    try:
        with destination.open("w", encoding="utf-8", newline="") as output:
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            writer.writerow(CATALOG_COLUMNS)
            for label, path in sources:
                if label in counts:
                    raise ResourceError(f"duplicate target-catalogue source label: {label}")
                count = 0
                with path.open("r", encoding="utf-8", newline="") as source:
                    reader = csv.DictReader(source, delimiter="\t")
                    if tuple(reader.fieldnames or ()) != CATALOG_COLUMNS:
                        raise ResourceError(
                            f"{label} target header must be exactly: " + "\t".join(CATALOG_COLUMNS)
                        )
                    for row in reader:
                        if None in row or any(value is None for value in row.values()):
                            raise ResourceError(f"{label} target row has wrong column count")
                        target_id = str(row["target_id"])
                        try:
                            database.execute("INSERT INTO targets VALUES (?)", (target_id,))
                        except sqlite3.IntegrityError as exc:
                            raise ResourceError(
                                f"duplicate target_id while combining catalogues: {target_id}"
                            ) from exc
                        writer.writerow([row[key] for key in CATALOG_COLUMNS])
                        count += 1
                counts[label] = count
        database.commit()
    finally:
        database.close()
        if database_path.exists():
            database_path.unlink()
    return counts
