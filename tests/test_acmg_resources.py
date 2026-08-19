from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from genome_annotation_builder.acmg_resources import (
    sync_acmg_sf,
    verify_acmg_sf,
    verify_vendor_acmg_snapshot,
)
from genome_annotation_builder.constants import (
    ACMG_SF_CANONICAL_SHA256,
    ACMG_SF_GENE_COUNT,
)
from genome_annotation_builder.errors import ResourceError


def test_vendored_acmg_snapshot_is_complete_and_pinned() -> None:
    snapshot = verify_vendor_acmg_snapshot()
    assert snapshot["counts"]["genes"] == ACMG_SF_GENE_COUNT
    assert snapshot["counts"] == {
        "genes": 84,
        "policy_gene_disease_pairs": 100,
        "clingen_gene_disease_rows": 177,
        "reporting_guidance_rows": 84,
        "frequency_specification_rows": 104,
    }
    additions = {
        row["gene_symbol"]
        for row in snapshot["policy_gene_disease_pairs"]
        if row["introduced_version"] == "3.3"
    }
    assert additions == {"ABCD1", "CYP27A1", "PLN"}
    assert len(ACMG_SF_CANONICAL_SHA256) == 64


def test_acmg_snapshot_installs_offline_and_detects_tampering(tmp_path: Path) -> None:
    config = SimpleNamespace(resource_cache=tmp_path / "resources", allow_download=False)
    manifest = sync_acmg_sf(config)  # type: ignore[arg-type]
    assert manifest["counts"]["genes"] == 84
    assert manifest["source_archive_present"] is False
    assert verify_acmg_sf(config) == manifest  # type: ignore[arg-type]

    canonical = config.resource_cache / "acmg_sf" / "3.3" / "acmg_sf_v3.3_20260815.json"
    document = json.loads(canonical.read_text(encoding="utf-8"))
    document["policy_gene_disease_pairs"][0]["gene_symbol"] = "TAMPERED"
    canonical.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(ResourceError, match="SHA-256"):
        verify_acmg_sf(config)  # type: ignore[arg-type]
