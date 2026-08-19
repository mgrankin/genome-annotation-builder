"""Project ClinGen gene-level context without creating variant-level claims."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .clingen_resources import verify_clingen_resources
from .config import BuildConfig
from .errors import ResourceError
from .util import sha256_file


def _project_table(source: Path, symbols: set[str], destination: Path) -> pa.Table:
    try:
        table = pq.read_table(source)
    except (OSError, pa.ArrowException) as exc:
        raise ResourceError(f"cannot read ClinGen source table {source}: {exc}") from exc
    if not symbols:
        return table.slice(0, 0)
    mask = pc.is_in(table["gene_symbol"], value_set=pa.array(sorted(symbols), type=pa.string()))
    selected = table.filter(mask)
    if selected.num_rows:
        pq.write_table(
            selected,
            destination,
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
            write_statistics=True,
            version="2.6",
        )
    return selected


def project_clingen_gene_context(
    config: BuildConfig,
    gene_symbols: set[str],
    destination: str | Path,
) -> dict[str, Any]:
    """Materialize only curations for exact HGNC symbols present in one bundle.

    These rows remain gene-level context.  They are deliberately not copied to
    ``clinical_evidence.parquet`` and are never interpreted as assertions about
    a particular allele or subject.
    """

    if not config.clingen.enabled:
        return {"status": "not_enabled"}
    resource = verify_clingen_resources(config)
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    gene_path = destination / "clingen_gene_disease.parquet"
    dosage_path = destination / "clingen_dosage_sensitivity.parquet"
    gene = _project_table(config.clingen.gene_disease_parquet, gene_symbols, gene_path)
    dosage = _project_table(config.clingen.dosage_parquet, gene_symbols, dosage_path)
    matched = set(gene["gene_symbol"].to_pylist()) | set(dosage["gene_symbol"].to_pylist())
    outputs: dict[str, dict[str, Any]] = {}
    if gene.num_rows:
        outputs[gene_path.name] = {"rows": gene.num_rows, "sha256": sha256_file(gene_path)}
    if dosage.num_rows:
        outputs[dosage_path.name] = {
            "rows": dosage.num_rows,
            "sha256": sha256_file(dosage_path),
        }
    return {
        "status": "PASS",
        "snapshot": resource["snapshot"],
        "manifest_sha256": resource["manifest_sha256"],
        "canonical_sha256": resource["canonical_sha256"],
        "semantic_boundary": resource["semantic_boundary"],
        "input_gene_symbols": len(gene_symbols),
        "matched_gene_symbols": len(matched),
        "unmatched_gene_symbols": len(gene_symbols - matched),
        "gene_disease_rows": gene.num_rows,
        "dosage_rows": dosage.num_rows,
        "outputs": outputs,
    }
