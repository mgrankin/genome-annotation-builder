from __future__ import annotations

import csv
import json
from datetime import date
from pathlib import Path

import pytest

from genome_annotation_builder.clingen_freeze import (
    CLINGEN_DOSAGE_HEADER,
    CLINGEN_GENE_DISEASE_HEADER,
    freeze_clingen_snapshot,
)
from genome_annotation_builder.clingen_resources import (
    sync_clingen_resources,
    verify_clingen_resources,
    verify_vendor_clingen_snapshot,
)
from genome_annotation_builder.config import load_config
from genome_annotation_builder.errors import ResourceError
from genome_annotation_builder.util import canonical_json_bytes


def _gene_disease_csv(path: Path, *, classification: str = "Definitive") -> None:
    separator = ["+" * len(value) for value in CLINGEN_GENE_DISEASE_HEADER]
    rows = [
        ["CLINGEN GENE DISEASE VALIDITY CURATIONS", *([""] * 9)],
        ["FILE CREATED: 2026-01-01", *([""] * 9)],
        ["WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity", *([""] * 9)],
        separator,
        list(CLINGEN_GENE_DISEASE_HEADER),
        separator,
        [
            "GENE1",
            "HGNC:1",
            "synthetic condition",
            "MONDO:0000001",
            "AD",
            "SOP4",
            classification,
            "https://search.clinicalgenome.org/kb/gene-validity/CGGCIEX:assertion_1",
            "2017-02-10T00:00:00",
            "Synthetic Gene Curation Expert Panel",
        ],
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, lineterminator="\n").writerows(rows)


def _dosage_tsv(path: Path, *, hi_score: str = "30") -> None:
    values = {field: "" for field in CLINGEN_DOSAGE_HEADER}
    values.update(
        {
            "Gene Symbol": "GENE1",
            "Gene ID": "1",
            "cytoBand": "1p36.33",
            "Genomic Location": "chr1:100-200",
            "Haploinsufficiency Score": hi_score,
            "Haploinsufficiency Description": (
                "Gene associated with autosomal recessive phenotype"
                if hi_score == "30"
                else "No evidence available"
            ),
            "Haploinsufficiency PMID1": "123; 456\u200b",
            "Triplosensitivity Score": "Not yet evaluated",
            "Triplosensitivity Description": "Not yet evaluated",
            "Date Last Evaluated": "2025-12-31",
            "Haploinsufficiency Disease ID": "MONDO:0000001",
        }
    )
    path.write_text(
        "\n".join(
            [
                "#ClinGen Gene Curation Results",
                "#01 Jan,2026",
                "#Genomic Locations are reported on GRCh38 (hg38): GCF_000001405.36",
                "#https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen",
                "#to create link: https://example.invalid/Gene Symbol",
                "#" + "\t".join(CLINGEN_DOSAGE_HEADER),
                "\t".join(values[field] for field in CLINGEN_DOSAGE_HEADER),
                "",
            ]
        ),
        encoding="utf-8",
    )


def test_clingen_freeze_is_typed_lossless_and_semantically_bounded(tmp_path: Path) -> None:
    gene = tmp_path / "gene.csv"
    dosage = tmp_path / "dosage.tsv"
    _gene_disease_csv(gene)
    _dosage_tsv(dosage)
    document = freeze_clingen_snapshot(
        gene,
        dosage,
        gene_disease_snapshot=date(2026, 1, 1),
        dosage_snapshot=date(2026, 1, 1),
        accessed_at=date(2026, 1, 2),
        gene_disease_url="https://example.invalid/gene.csv",
        dosage_url="https://example.invalid/dosage.tsv",
    )
    assert document["license"] == "CC0-1.0"
    assert "not evidence that the variant is pathogenic" in document["semantic_boundary"]
    gene_row = document["gene_disease_validity"][0]
    assert gene_row["assertion_id"] == "CGGCIEX:assertion_1"
    assert gene_row["classification_date"] == "2017-02-10T00:00:00"
    dosage_row = document["dosage_sensitivity"][0]
    assert dosage_row["haploinsufficiency_score"] == 30
    assert dosage_row["triplosensitivity_score"] is None
    assert dosage_row["haploinsufficiency_pmids"] == [123, 456]
    assert dosage_row["haploinsufficiency_pmid_source_values"] == ["123; 456\u200b"]
    assert canonical_json_bytes(document) == canonical_json_bytes(document)


def test_clingen_freeze_rejects_unknown_classification_and_score(tmp_path: Path) -> None:
    gene = tmp_path / "gene.csv"
    dosage = tmp_path / "dosage.tsv"
    _gene_disease_csv(gene, classification="Probably true")
    _dosage_tsv(dosage)
    with pytest.raises(ResourceError, match="unknown ClinGen classification"):
        freeze_clingen_snapshot(
            gene,
            dosage,
            gene_disease_snapshot=date(2026, 1, 1),
            dosage_snapshot=date(2026, 1, 1),
            accessed_at=date(2026, 1, 2),
            gene_disease_url="https://example.invalid/gene.csv",
            dosage_url="https://example.invalid/dosage.tsv",
        )

    _gene_disease_csv(gene)
    _dosage_tsv(dosage, hi_score="4")
    with pytest.raises(ResourceError, match="score/description mismatch"):
        freeze_clingen_snapshot(
            gene,
            dosage,
            gene_disease_snapshot=date(2026, 1, 1),
            dosage_snapshot=date(2026, 1, 1),
            accessed_at=date(2026, 1, 2),
            gene_disease_url="https://example.invalid/gene.csv",
            dosage_url="https://example.invalid/dosage.tsv",
        )


def _config(tmp_path: Path) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(
        f"""
project:
  output_dir: {tmp_path / "output"}
  work_dir: {tmp_path / "work"}
  bundle_name: synthetic.genome
input: {{vcf: synthetic.vcf.gz, sample: SAMPLE, genome_build: GRCh38}}
reference: {{fasta: synthetic.fa, contig_style: chr}}
resources:
  cache_dir: {tmp_path / "resources"}
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
annotation:
  vep: {{enabled: true}}
  clingen:
    enabled: true
    snapshot: "2026-08-15"
    root: {tmp_path / "resources" / "clingen" / "2026-08-15"}
    gene_disease_parquet: {tmp_path / "resources" / "clingen" / "2026-08-15" / "normalized" / "gene_disease_validity.parquet"}
    dosage_parquet: {tmp_path / "resources" / "clingen" / "2026-08-15" / "normalized" / "dosage_sensitivity.parquet"}
    manifest: {tmp_path / "resources" / "clingen" / "2026-08-15" / "resource_manifest.json"}
    license: CC0-1.0
    redistribution_allowed: true
""",
        encoding="utf-8",
    )
    return path


def test_vendored_clingen_snapshot_and_materialized_indexes_verify(tmp_path: Path) -> None:
    vendor = verify_vendor_clingen_snapshot()
    assert vendor["gene_disease_rows"] == 3663
    assert vendor["dosage_rows"] == 1526
    config = load_config(_config(tmp_path))
    first = sync_clingen_resources(config)
    second = sync_clingen_resources(config)
    assert first["status"] == "verified"
    assert second == first
    assert first["outputs"]["gene_disease_validity"]["rows"] == 3663
    assert first["outputs"]["dosage_sensitivity"]["rows"] == 1526

    replay_root = tmp_path / "independent-replay"
    replay_root.mkdir()
    replay_config = load_config(_config(replay_root))
    replay = sync_clingen_resources(replay_config)
    assert replay["manifest_sha256"] == first["manifest_sha256"]
    assert replay_config.clingen.manifest.read_bytes() == config.clingen.manifest.read_bytes()
    assert replay_config.clingen.gene_disease_parquet.read_bytes() == (
        config.clingen.gene_disease_parquet.read_bytes()
    )
    assert replay_config.clingen.dosage_parquet.read_bytes() == (
        config.clingen.dosage_parquet.read_bytes()
    )

    manifest = json.loads(config.clingen.manifest.read_text(encoding="utf-8"))
    manifest["semantic_boundary"] = "a gene match is a diagnosis"
    config.clingen.manifest.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ResourceError, match="manifest boundary mismatch"):
        verify_clingen_resources(config)
