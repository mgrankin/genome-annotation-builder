"""Static, generic assets for the .genome Agent Exploration Profile."""

from __future__ import annotations

from typing import Any

README = """# Read my genome bundle

This directory is a `.genome/1.1` snapshot. Begin with `manifest.json`,
`inventory.json`, and `omissions.json`. Validate integrity before querying:

```bash
python readmygenome/scripts/query.py --bundle . overview
```

Use only bounded named queries. Never infer a reference genotype or a negative
test from an empty variant lookup. Use the `callability` operation and inspect
its explicit conclusion and evidence rows. Annotation domains listed as omitted
make no positive or negative claim.

This is research data, not a diagnosis. Clinically consequential results need
appropriate confirmation and professional interpretation.
"""


SKILL = """---
name: readmygenome
description: Safely explore one local .genome/1.1 bundle with bounded, evidence-aware queries.
---

# Read My Genome

Use this guide only for the bundle containing it.

## Required workflow

1. Run `python readmygenome/scripts/query.py --bundle . overview`. Stop if
   integrity validation fails.
2. Read `manifest.json`, `inventory.json`, `omissions.json`, and the applicable
   glossary entries before interpreting a field.
3. Query with the bundled named operations. Keep `--limit` small and iterate;
   never load an entire Parquet table into model context.
4. Search all relevant available tables. Continue with synonyms, rsIDs, gene
   symbols, and related identifiers after an empty retrieval.
5. Before any absence statement, run `callability` for the exact locus. Only
   `reference_observed` supports a threshold-qualified reference statement.
   Apply the callability query to the whole allele. A producer-declared PAR limit
   can yield unknown even when source evidence records a passing reference GT.
   `overlapping_variant_observed`, `callable_without_explicit_reference`,
   `unknown_uncallable`, `unknown_no_call`, and `unknown` do not.
6. Keep observed genotype, pipeline-derived values, external annotations,
   reviewed interpretations, and general knowledge visibly separate.
7. Cite table, stable row identifier, annotation source/snapshot, evidence
   snapshot, and review date when available.
8. Surface conflicts and uncertainty. Do not derive a PGx phenotype from
   variant rows. Use `pharmacogenomics.parquet` only for published calls and
   inspect the PGx gene/locus/guidance audits for ambiguity, missing evidence,
   complex-gene exclusion, or another explicit withholding reason.
9. Treat family evidence as a subject-only confidence/review layer. It never
   changes GT, and `not_individually_persisted` LOCO evidence is not proof of
   consistency. Do not infer or reconstruct a relative's genotype.
10. Treat GWAS matches as population-association retrieval only. A carried
    effect allele is not a personal diagnosis, causal claim, penetrance
    estimate, or clinical recommendation; inspect the exact source row and
    snapshot before summarizing it.
11. Treat ClinGen gene-disease and dosage rows as gene-level context only. A
    gene match does not classify the subject's variant and is not a diagnosis;
    an absent row is unavailable context rather than a negative result.
12. Explain results directly while retaining clinical-confirmation and
   medical-action boundaries.

## Named operations

`overview`, `search`, `variant`, `gene`, `clinical`, `pgx`, `prs`, `gwas`,
`callability`, and `evidence` return JSON by default. They are read-only and
bounded. See `reference/query_catalog.md` for examples.

Operational files are non-authoritative. The manifest-listed tables, their
hashes, evidence classes, and declared rules remain authoritative.
"""


QUERY_CATALOG = """# Query catalogue

All examples automatically validate bundle integrity and return bounded JSON.

| Goal | Command | Evidence boundary |
|---|---|---|
| Orient to snapshot | `query.py overview` | Reports available tables and omissions |
| Search an rsID | `query.py variant rs123` | Empty is a retrieval outcome |
| Search a normalized allele | `query.py variant chr1:100:A:G` | Exact GRCh38 allele |
| Search a gene | `query.py gene GENE1` | Variants plus optional ClinGen gene-level context; no variant classification is inferred |
| Broad lookup | `query.py search TERM` | Searches bounded core fields |
| Clinical condition | `query.py clinical TERM` | Searches available variant evidence, policy rows, and gene-level context without inferring a diagnosis |
| Medication/PGx | `query.py pgx DRUG` | Returns published calls or explicit gene/guidance audit limitations |
| Polygenic trait | `query.py prs TRAIT` | Requires PRS table and matching ancestry context |
| GWAS association | `query.py gwas TRAIT_OR_RSID` | Population association only; exact allele and source audit required |
| Locus absence | `query.py callability chr1:100:A:G` | Never infer from empty variants alone |
| Source/family evidence | `query.py evidence chr1:100:A:G` | Maps a normalized row to immutable source and optional subject-only family evidence |

Use `--limit 1..500`; the default is 50. Use `--format table` only for a small
human-readable result. Named operations do not accept arbitrary SQL.
"""


SYNONYMS: dict[str, Any] = {
    "schema_version": "1.0.0",
    "authority": "operational_guidance",
    "conditions": {
        "myocardial infarction": ["heart attack", "MI"],
        "type 2 diabetes mellitus": ["type 2 diabetes", "T2D", "T2DM"],
    },
    "drugs": {
        "acetylsalicylic acid": ["aspirin"],
        "clopidogrel": ["Plavix"],
    },
    "traits": {
        "low-density lipoprotein cholesterol": ["LDL cholesterol", "LDL-C"],
        "body mass index": ["BMI"],
    },
    "genes": {},
    "warning": "Synonyms expand retrieval only and never override bundle evidence.",
}


def canonical_questions() -> dict[str, Any]:
    topics = [
        (
            "orientation",
            "What data and annotation snapshots are available?",
            ["manifest", "inventory"],
        ),
        (
            "important_findings",
            "What clinical-grade findings are present?",
            ["clinical_findings", "variants"],
        ),
        (
            "gene_lookup",
            "What variants and gene-level context are available for a requested gene?",
            ["variants", "clingen_gene_disease", "clingen_dosage_sensitivity"],
        ),
        ("variant_lookup", "Is a requested normalized variant present?", ["variants"]),
        ("rsid_lookup", "Is a requested rsID present?", ["variants"]),
        (
            "condition_lookup",
            "What evidence concerns a requested condition?",
            [
                "clinical_findings",
                "clinical_evidence",
                "clingen_gene_disease",
                "clingen_dosage_sensitivity",
            ],
        ),
        (
            "phenotype_lookup",
            "What evidence concerns a requested phenotype?",
            ["clinical_findings", "interpretations"],
        ),
        (
            "medication_review",
            "Are there gene-level PGx calls for a requested medication?",
            ["pharmacogenomics", "pgx_gene_call_audit", "pgx_guidance_audit"],
        ),
        (
            "pgx_gene",
            "What is the reported PGx diplotype for a requested gene?",
            ["pharmacogenomics", "pgx_gene_call_audit"],
        ),
        (
            "carrier_status",
            "Is a carrier-status finding reported?",
            ["clinical_findings", "clinical_evidence"],
        ),
        ("prs_lookup", "Is a polygenic score available for a requested trait?", ["prs"]),
        (
            "gwas_lookup",
            "Are trait associations recorded for a requested variant?",
            ["gwas_associations", "gwas_catalog_evidence", "variants"],
        ),
        (
            "rare_variants",
            "Which variants satisfy the declared rarity rule?",
            ["variants", "manifest"],
        ),
        (
            "damaging_variants",
            "Which variants satisfy the declared prediction consensus?",
            ["variants", "manifest"],
        ),
        (
            "clinvar_conflicts",
            "Which reported ClinVar interpretations have conflicts?",
            ["variants", "clinical_evidence"],
        ),
        (
            "exact_reference",
            "Is an exact locus explicitly observed as reference?",
            ["variants", "callability", "reference_evidence"],
        ),
        (
            "callable_only",
            "Is a requested locus callable without explicit reference evidence?",
            ["callable_regions", "reference_evidence"],
        ),
        (
            "no_call",
            "Does a requested locus have incomplete genotype evidence?",
            ["reference_evidence"],
        ),
        (
            "unknown_locus",
            "What can be said when no locus evidence is available?",
            ["variants", "callable_regions", "reference_evidence"],
        ),
        (
            "source_provenance",
            "What source record produced a normalized variant?",
            ["variant_provenance", "variants"],
        ),
        (
            "family_evidence",
            "What subject-only family evidence is available for a normalized variant?",
            ["variants", "family_evidence"],
        ),
        (
            "normalization",
            "Was a variant changed by normalization?",
            ["variant_provenance", "transformations"],
        ),
        (
            "non_core",
            "Were symbolic or non-core-contig calls retained separately?",
            ["non_core_variants", "omissions"],
        ),
        (
            "ancestry_limits",
            "What ancestry limitations apply to population or score results?",
            ["manifest", "omissions", "prs"],
        ),
        ("snapshot_age", "How old are interpretation resources?", ["manifest"]),
        (
            "new_literature",
            "Can current literature be distinguished from this snapshot?",
            ["manifest", "clinical_evidence"],
        ),
        (
            "guideline_update",
            "Does a newer guideline change a bundle claim?",
            ["manifest", "transformations"],
        ),
        (
            "broad_query",
            "Summarize available evidence without loading complete tables.",
            ["inventory", "manifest"],
        ),
        (
            "empty_retrieval",
            "What should happen after an empty first search?",
            ["synonyms", "inventory"],
        ),
        (
            "clinical_summary",
            "Produce a clinician-facing evidence summary with confirmation boundaries.",
            ["clinical_findings", "clinical_evidence", "manifest"],
        ),
        ("integrity_failure", "What happens when a declared file hash is wrong?", ["manifest"]),
    ]
    cases = []
    for index, (topic, question, tables) in enumerate(topics, start=1):
        cases.append(
            {
                "case_id": f"AE-{index:03d}",
                "topic": topic,
                "question": question,
                "expected_tables_searched": tables,
                "required_findings": ["Report only rows and availability supported by the bundle."],
                "claims_that_must_not_be_made": [
                    "Do not turn an empty retrieval or omitted domain into a negative genetic claim."
                ],
                "expected_callability_behaviour": (
                    "Check exact reference/no-call evidence before any locus-specific absence statement."
                ),
                "expected_citation_fields": [
                    "table",
                    "row_identifier",
                    "generated_at",
                    "annotation_source_or_omission",
                ],
            }
        )
    return {"profile_version": "1.0.0", "cases": cases}
