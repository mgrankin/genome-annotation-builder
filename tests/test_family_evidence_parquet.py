from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from genome_annotation_builder.family_evidence import _iter_core_groups


def _dictionary_string() -> pa.DataType:
    return pa.dictionary(pa.int32(), pa.string())


def test_core_iterator_projects_nested_leaves_across_row_groups(tmp_path: Path) -> None:
    """PyArrow cannot project the whole dictionary-backed struct across row groups."""

    schema = pa.schema(
        [
            pa.field("variant_id", pa.string(), nullable=False),
            pa.field("chrom", _dictionary_string(), nullable=False),
            pa.field("pos", pa.int32(), nullable=False),
            pa.field("ref", pa.string(), nullable=False),
            pa.field("alt", pa.string(), nullable=False),
            pa.field(
                "genotype",
                pa.struct(
                    [
                        pa.field("gt", pa.list_(pa.int8()), nullable=False),
                        pa.field("phased", pa.bool_(), nullable=False),
                        pa.field("zygosity", _dictionary_string(), nullable=False),
                        pa.field("allele_balance", pa.float32()),
                    ]
                ),
                nullable=False,
            ),
        ]
    )
    rows = [
        {
            "variant_id": "chr1:10:A:G",
            "chrom": "chr1",
            "pos": 10,
            "ref": "A",
            "alt": "G",
            "genotype": {
                "gt": [0, 1],
                "phased": False,
                "zygosity": "heterozygous",
                "allele_balance": 0.5,
            },
        },
        {
            "variant_id": "chr1:20:C:T",
            "chrom": "chr1",
            "pos": 20,
            "ref": "C",
            "alt": "T",
            "genotype": {
                "gt": [1, 0],
                "phased": True,
                "zygosity": "heterozygous",
                "allele_balance": 0.45,
            },
        },
    ]
    partition = tmp_path / "variants.parquet" / "chrom=chr1"
    partition.mkdir(parents=True)
    parquet_path = partition / "part-00000.parquet"
    with pq.ParquetWriter(parquet_path, schema, compression="zstd") as writer:
        for row in rows:
            writer.write_table(pa.Table.from_pylist([row], schema=schema))

    groups = list(_iter_core_groups(tmp_path, "chr1"))

    assert [pos for pos, _variants in groups] == [10, 20]
    assert groups[0][1]["chr1:10:A:G"].gt == (0, 1)
    assert groups[0][1]["chr1:10:A:G"].phased is False
    assert groups[1][1]["chr1:20:C:T"].gt == (1, 0)
    assert groups[1][1]["chr1:20:C:T"].phased is True
