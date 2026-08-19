#!/usr/bin/env python3
"""
Reference reader for `.genome/1.x` bundles.

A deliberately minimal implementation that demonstrates the spec is
implementable. Production readers will be more sophisticated; this one
exists to:

  1. Validate that a bundle is conformant
  2. Expose the bundle's data through a simple Python API
  3. Serve as the canonical reference behaviour for ambiguous cases

Usage:
  # Validate a bundle
  python reader.py /path/to/bundle/

  # Use as a library
  from reader import GenomeBundle
  bundle = GenomeBundle('/path/to/bundle/')
  bundle.validate()
  variants = bundle.query("SELECT * FROM variants WHERE pathogenicity.is_actionable = true")
  for row in variants:
      print(row)
"""

import sys
import json
import hashlib
import re
from pathlib import Path
from typing import Iterator, Optional


SUPPORTED_SCHEMA_MAJOR = 1
AGENT_EXPLORATION_REQUIRED_PATHS = (
    "field_glossary.json",
    "inventory.json",
    "raw_artifacts.json",
    "transformations.json",
    "omissions.json",
    "relationships.json",
    "readmygenome.md",
    "readmygenome/SKILL.md",
    "readmygenome/scripts/query.py",
    "readmygenome/reference/query_catalog.md",
    "readmygenome/reference/synonyms.json",
    "readmygenome/evals/canonical_questions.json",
)


class GenomeBundle:
    """Minimal reader for a `.genome/1.x` bundle."""

    def __init__(self, path):
        self.path = Path(path)
        if not self.path.is_dir():
            raise ValueError(f"{path} is not a directory")

        manifest_file = self.path / "manifest.json"
        if not manifest_file.exists():
            raise ValueError(f"manifest.json not found in {path}")
        self.manifest = json.loads(manifest_file.read_text())

        schema_file = self.path / "schema.json"
        if not schema_file.exists():
            raise ValueError(f"schema.json not found in {path}")
        self.schema = json.loads(schema_file.read_text())

        self._con = None  # lazy DuckDB connection

    # ---------- validation ----------

    def validate(self) -> dict:
        """Validate the bundle. Returns a structured report.

        Raises ValueError if the bundle is non-conformant in a way that
        would make further reading unsafe (missing required files, schema
        version mismatch). Returns warnings for non-fatal issues.
        """
        errors, warnings = [], []

        sv = self.manifest.get("schema_version")
        match = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)", str(sv or ""))
        if not match or int(match.group(1)) != SUPPORTED_SCHEMA_MAJOR:
            errors.append(
                f"schema_version {sv!r} is not a supported 1.x semantic version"
            )

        schema_version = self.schema.get("version")
        if schema_version != sv:
            errors.append(
                f"schema.json version {schema_version!r} does not match manifest {sv!r}"
            )

        for required in ("manifest.json", "schema.json"):
            if not (self.path / required).exists():
                errors.append(f"missing required file: {required}")

        if not (self.path / "variants.parquet").is_dir():
            errors.append("missing required directory: variants.parquet/")

        profiles = self.manifest.get("profiles") or {}
        if profiles.get("agent_exploration"):
            for required in AGENT_EXPLORATION_REQUIRED_PATHS:
                if not (self.path / required).exists():
                    errors.append(
                        f"Agent Exploration Profile missing required path: {required}"
                    )

        # Hash verification
        for rel, meta in self.manifest.get("files", {}).items():
            file_path = self.path / rel
            if not file_path.exists():
                errors.append(f"file declared in manifest is missing: {rel}")
                continue
            expected = meta.get("sha256")
            if not expected:
                continue
            actual = self._sha256(file_path)
            if actual != expected:
                errors.append(f"hash mismatch for {rel}: expected {expected[:12]}…, got {actual[:12]}…")

        # Snapshot age warning
        gen_at = self.manifest.get("generated_at")
        if gen_at:
            try:
                from datetime import datetime, timezone
                gen = datetime.fromisoformat(gen_at.replace("Z", "+00:00"))
                age_days = (datetime.now(timezone.utc) - gen).days
                if age_days > 180:
                    warnings.append(
                        f"bundle is {age_days} days old; interpretive content may be stale"
                    )
            except Exception:
                pass

        if errors:
            raise ValueError(f"bundle is non-conformant:\n  - " + "\n  - ".join(errors))

        return {"valid": True, "warnings": warnings}

    # ---------- queries ----------

    def query(self, sql: str, params=None) -> list:
        """Run a SQL query against the bundle. Returns list of dicts."""
        con = self._connection()
        if params:
            cur = con.execute(sql, params)
        else:
            cur = con.execute(sql)
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]

    def actionable_variants(self) -> list:
        """Convenience: return variants meeting the bundle's actionable rule."""
        return self.query("""
            SELECT variant_id, gene.symbol AS gene, consequence.hgvsp,
                   pathogenicity.clinvar_significance,
                   pathogenicity.clinvar_review_stars,
                   predictions.revel
            FROM variants
            WHERE pathogenicity.is_actionable = true
              AND clinical_grade = true
            ORDER BY predictions.revel DESC NULLS LAST
        """)

    def pgx_phenotypes(self) -> list:
        """Convenience: return all gene-level pharmacogenomic phenotype calls."""
        if not (self.path / "pharmacogenomics.parquet").exists():
            return []
        return self.query("""
            SELECT gene_symbol, diplotype, phenotype, cpic_level, affected_drugs
            FROM pharmacogenomics
            ORDER BY gene_symbol
        """)

    def variants_in_gene(self, gene_symbol: str) -> list:
        """Convenience: return all variants in a gene."""
        return self.query("""
            SELECT variant_id, consequence.hgvsp,
                   pathogenicity.clinvar_significance,
                   pathogenicity.is_actionable,
                   predictions.revel
            FROM variants
            WHERE gene.symbol = ?
            ORDER BY predictions.revel DESC NULLS LAST
        """, [gene_symbol])

    def rule(self, name: str) -> Optional[dict]:
        """Return the manifest's rule definition for a pre-computed field."""
        return self.manifest.get("consensus_definitions", {}).get(name)

    @property
    def snapshot_date(self) -> Optional[str]:
        return self.manifest.get("generated_at")

    @property
    def annotation_sources(self) -> dict:
        return self.manifest.get("annotation_sources", {})

    # ---------- internals ----------

    def _connection(self):
        if self._con is None:
            try:
                import duckdb
            except ImportError:
                raise RuntimeError("duckdb required: pip install duckdb")
            self._con = duckdb.connect()

            variants_path = self.path / "variants.parquet"
            if variants_path.exists():
                self._con.execute(f"""
                    CREATE VIEW variants AS
                    SELECT * FROM read_parquet('{variants_path}/**/*.parquet',
                                               hive_partitioning=true)
                """)

            pgx_path = self.path / "pharmacogenomics.parquet"
            if pgx_path.exists():
                self._con.execute(
                    f"CREATE VIEW pharmacogenomics AS SELECT * FROM read_parquet('{pgx_path}')"
                )

            gene_index_path = self.path / "gene_index.parquet"
            if gene_index_path.exists():
                self._con.execute(
                    f"CREATE VIEW gene_index AS SELECT * FROM read_parquet('{gene_index_path}')"
                )

            prs_path = self.path / "prs.parquet"
            if prs_path.exists():
                self._con.execute(
                    f"CREATE VIEW prs AS SELECT * FROM read_parquet('{prs_path}')"
                )

        return self._con

    @staticmethod
    def _sha256(path: Path) -> str:
        h = hashlib.sha256()
        if path.is_dir():
            for p in sorted(path.rglob("*")):
                if p.is_file():
                    h.update(p.relative_to(path).as_posix().encode())
                    h.update(b"\0")
                    with p.open("rb") as f:
                        while chunk := f.read(1024 * 1024):
                            h.update(chunk)
        else:
            with path.open("rb") as f:
                while chunk := f.read(1024 * 1024):
                    h.update(chunk)
        return h.hexdigest()


# ---------- CLI ----------

def main():
    if len(sys.argv) != 2:
        print("usage: python reader.py <bundle_path>", file=sys.stderr)
        sys.exit(2)

    bundle = GenomeBundle(sys.argv[1])
    try:
        report = bundle.validate()
        print(f"✓ bundle is conformant to .genome/{bundle.manifest['schema_version']}")
        print(f"  generated: {bundle.snapshot_date}")
        print(f"  pipeline:  {bundle.manifest.get('pipeline_version')}")
        for w in report["warnings"]:
            print(f"  ⚠ {w}")

        actionable = bundle.actionable_variants()
        pgx = bundle.pgx_phenotypes()
        print(f"\n  {len(actionable)} actionable variant(s)")
        print(f"  {len(pgx)} pharmacogenomic call(s)")

    except ValueError as e:
        print(f"✗ {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
