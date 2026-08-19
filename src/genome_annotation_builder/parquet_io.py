"""Bounded-memory deterministic Parquet dataset writers."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq


class PartitionedParquetWriter:
    """Write one stable Parquet file per partition with bounded row groups."""

    def __init__(
        self,
        root: str | Path,
        schema: pa.Schema,
        *,
        partition_field: str | None = "chrom",
        batch_size: int = 50_000,
    ) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=False)
        self.schema = schema
        self.partition_field = partition_field
        self.batch_size = batch_size
        self._buffers: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self._writers: dict[str, pq.ParquetWriter] = {}
        self.rows = 0
        self.rows_by_partition: dict[str, int] = defaultdict(int)
        self._closed = False

    def _partition(self, row: dict[str, Any]) -> str:
        if self.partition_field is None:
            return "_unpartitioned"
        value = str(row[self.partition_field])
        if not value or any(
            char not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-"
            for char in value
        ):
            raise ValueError(f"unsafe Parquet partition value: {value!r}")
        return value

    def write(self, row: dict[str, Any]) -> None:
        if self._closed:
            raise RuntimeError("writer is closed")
        partition = self._partition(row)
        buffer = self._buffers[partition]
        buffer.append(row)
        self.rows += 1
        self.rows_by_partition[partition] += 1
        if len(buffer) >= self.batch_size:
            self._flush(partition)

    def _writer(self, partition: str) -> pq.ParquetWriter:
        existing = self._writers.get(partition)
        if existing is not None:
            return existing
        if self.partition_field is None:
            directory = self.root
        else:
            directory = self.root / f"{self.partition_field}={partition}"
            directory.mkdir(parents=True, exist_ok=True)
        writer = pq.ParquetWriter(
            directory / "part-00000.parquet",
            self.schema,
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
            write_statistics=True,
            version="2.6",
        )
        self._writers[partition] = writer
        return writer

    def _flush(self, partition: str) -> None:
        rows = self._buffers[partition]
        if not rows:
            return
        table = pa.Table.from_pylist(rows, schema=self.schema)
        self._writer(partition).write_table(table, row_group_size=self.batch_size)
        rows.clear()

    def close(self) -> None:
        if self._closed:
            return
        for partition in sorted(self._buffers):
            self._flush(partition)
        for writer in self._writers.values():
            writer.close()
        if not self._writers:
            pq.write_table(
                pa.Table.from_pylist([], schema=self.schema),
                self.root / "part-00000.parquet",
                compression="zstd",
                version="2.6",
            )
        pq.write_metadata(self.schema, self.root / "_common_metadata")
        self._closed = True

    def __enter__(self) -> PartitionedParquetWriter:
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if exc_type is None:
            self.close()
        else:
            for writer in self._writers.values():
                writer.close()
