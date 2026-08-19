# `.genome/1.1` — Bundle specification

**Version:** 1.1.0
**Status:** Public release
**Date:** August 2026
**Editor:** David Ball, genome.computer
**License:** Apache 2.0

---

## 1. Introduction

### 1.1 Purpose

This document specifies the `.genome/1.1` bundle format: an opinionated, versioned container for a single individual's genomic data structured to be read correctly by general-purpose software agents (including large language models) without requiring out-of-band knowledge or specialised tooling.

The bundle is designed around three principles, defined normatively in §3:

1. **Self-describing fields.** Every value's meaning is recoverable from the data alone.
2. **Self-contained correctness.** Every interpretive claim is computable from the bundle without external lookup.
3. **Explicit reference frames.** Every value that depends on a reference (allele, threshold, guideline list) carries that reference adjacent.
4. **Explicit evidence and absence boundaries.** The bundle distinguishes observed data, pipeline-derived values, external annotations, and reviewed interpretations, and states when a missing row has no negative meaning.

A `.genome` bundle is additive to, not a replacement for, a Variant Call Format (VCF) file. The bundle is intended to ship alongside the originating VCF, which remains the sovereignty record.

### 1.2 Scope

This specification covers:

- The on-disk layout of a `.genome/1.1` bundle (§4)
- The required and optional files within a bundle (§5)
- The schema of each Parquet table (§6)
- The structure and required fields of `manifest.json` (§7)
- The variant identifier and normalisation requirements (§8)
- The pre-computed interpretive fields and their declared rules (§9)
- Evidence, callability, provenance, omission, and relationship records (§6–§7)
- The optional Agent Exploration Profile for bundles intended for direct use by general-purpose agents (§12)
- Versioning and compatibility (§10)
- Validation requirements for compliant readers (§11)

This specification does NOT cover:

- Variant calling, imputation, or annotation pipelines
- Clinical interpretation or medical decision support
- Bundle generation tooling (the `genome-convert` reference implementation is separate)
- Network protocols or transport mechanisms
- Privacy, security, or authentication

### 1.3 Conformance terms

The key words **MUST**, **MUST NOT**, **REQUIRED**, **SHALL**, **SHALL NOT**, **SHOULD**, **SHOULD NOT**, **RECOMMENDED**, **MAY**, and **OPTIONAL** in this document are to be interpreted as described in [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

A bundle is **conformant** to `.genome/1.1` if and only if it satisfies every **MUST** and **REQUIRED** in this document. A reader is **conformant** if it implements every **MUST** in §11. Conformance to the Agent Exploration Profile is separate and requires the additional rules in §12.

### 1.4 Relationship to other specifications

- **VCF (Variant Call Format).** The bundle's `variants.parquet` is conceptually equivalent to a VCF's variant records, restructured for typed columnar access. Bundle generation typically takes an annotated VCF as input.
- **GRCh38.** v1.1 specifies GRCh38 (hg38) as the only supported reference assembly. GRCh37/hg19 is out of scope.
- **HGNC, Ensembl, RefSeq.** Gene symbols use HGNC. Transcript identifiers use Ensembl as canonical, RefSeq supported in `consequence.all_annotations`.
- **ClinVar, ACMG, CPIC, gnomAD.** Versions of these external sources MUST be recorded in `manifest.json.annotation_sources` (§7).
- **CPIC, PharmCAT.** Pharmacogenomic gene-level calls in `pharmacogenomics.parquet` follow the CPIC data model. PharmCAT is the reference haplotype caller; equivalent callers are acceptable.
- **Apache Parquet.** All tabular data is stored as Apache Parquet. The format version MUST be Parquet 2.6 or later. Compression MUST be zstd.

---

## 2. Terminology

- **Bundle.** A directory or tar archive conforming to §4 and containing the files specified in §5.
- **Variant.** A single row in `variants.parquet`, identified by a normalised `variant_id` (§8).
- **Pre-computed field.** A boolean or categorical value computed at bundle generation time against a rule declared in `manifest.json.consensus_definitions` (§9).
- **Evidence class.** One of `observed`, `pipeline_derived`, `external_annotation`, `reviewed_interpretation`, or `bundle_metadata`; identifies where a value or claim originated.
- **Callability.** Evidence about whether an assay and calling pipeline could reliably observe a genotype at a locus or interval. Callability is distinct from the presence or absence of a variant annotation row.
- **No matching record.** A retrieval outcome. It means only that the searched tables returned no matching row; it is not a negative genotype result unless site-level callability and an explicit reference or negative call support that conclusion.
- **Snapshot.** The set of interpretive values in a bundle at the moment it was generated. Snapshots are immutable; updated interpretations require a new bundle.
- **Recomputability boundary.** The conceptual division between durable variant-level data (which does not change) and dated interpretive content (which is regenerated from the same variants under updated rules). See §10.3.
- **Conformant reader.** Software that reads a `.genome/1.x` bundle and satisfies every **MUST** in §11.

---

## 3. Design invariants

These four invariants are normative. A bundle that satisfies the file-layout and schema requirements but violates an invariant is non-conformant.

### 3.1 Invariant I-1: Self-describing fields

Every value in the bundle MUST be interpretable from the schema and the value alone, without reference to any external document, header binding, or implicit field order. Values that are members of a fixed enumeration MUST be stored in dictionary-encoded columns whose value set is declared in `schema.json`.

### 3.2 Invariant I-2: Self-contained correctness

Every pre-computed interpretive field (§9) MUST be derivable from a rule fully specified in `manifest.json.consensus_definitions`, including the exact threshold, the source list and version, and the cryptographic hash of any external authoritative list (e.g. ACMG SF gene list).

### 3.3 Invariant I-3: Explicit reference frames

Every value whose meaning depends on a reference MUST carry that reference in the same row. Specifically:

- An effect size MUST co-occur with `effect_type` and `effect_allele` (§6.1.10).
- A pharmacogenomic phenotype MUST be reported at gene level with explicit diplotype (§6.2); variant-level rows MUST NOT carry phenotype claims.
- A rarity classification MUST be relative to a population identified in the manifest (§9.1).

### 3.4 Invariant I-4: Explicit evidence and absence boundaries

Every interpretation record MUST state its evidence class, review status, evidence snapshot, and whether confirmation is required. A bundle MUST NOT encode the absence of a row in `variants.parquet`, `clinical_findings.parquet`, `clinical_evidence.parquet`, or `pharmacogenomics.parquet` as proof of a reference genotype, negative test, or absence of disease.

A negative site-level genetic claim requires either an explicit genotype call or callability evidence that identifies the locus, assay scope, and confidence semantics. When that evidence is unavailable, readers MUST report “no matching record found” and the unavailable evidence, not “negative” or “not present.”

---

## 4. Bundle layout

A `.genome/1.1` bundle is a directory whose root contains the files specified in §5. The directory MAY be packaged as a tar archive for transport.

```
<bundle_root>/
  manifest.json              [REQUIRED]
  schema.json                [REQUIRED]
  variants.parquet/          [REQUIRED — directory, hive-partitioned]
    chrom=chr1/
      part-0000.parquet
    chrom=chr2/
      part-0000.parquet
    ...
  pharmacogenomics.parquet   [OPTIONAL]
  gene_index.parquet         [OPTIONAL]
  prs.parquet                [OPTIONAL]
  clinical_findings.parquet  [OPTIONAL]
  clinical_evidence.parquet  [OPTIONAL]
  interpretations.parquet    [OPTIONAL]
  gwas_associations.parquet  [OPTIONAL]
  callability.parquet        [OPTIONAL]
  callable_regions.parquet   [OPTIONAL]
  field_glossary.json        [OPTIONAL]
  inventory.json             [OPTIONAL]
  raw_artifacts.json         [OPTIONAL]
  transformations.json       [OPTIONAL]
  omissions.json             [OPTIONAL]
  relationships.json         [OPTIONAL]
  readmygenome.md             [OPTIONAL]
  readmygenome/              [OPTIONAL]
  findings.json              [OPTIONAL]
  report.pdf                 [OPTIONAL]
  reads.cram                 [OPTIONAL]
  reads.cram.crai            [OPTIONAL — REQUIRED if reads.cram present]
```

The bundle root MAY contain other files; conformant readers MUST ignore unknown files. Conformant readers MUST NOT depend on the absence of unknown files.

The bundle root MUST NOT contain symbolic links, device files, or named pipes.

---

## 5. Files

### 5.1 `manifest.json` (REQUIRED)

JSON document conforming to the manifest schema in §7. Encoding MUST be UTF-8. Indentation is OPTIONAL but RECOMMENDED.

### 5.2 `schema.json` (REQUIRED)

JSON Schema document defining the structure of the Parquet tables and the manifest. The shipped `schema.json` MUST be byte-identical to the version published at `https://github.com/genome-computer/genome-spec` for the declared `schema_version`.

### 5.3 `variants.parquet/` (REQUIRED)

Directory containing variant data, hive-partitioned by `chrom`. Each partition contains one or more Parquet files. See §6.1 for schema.

The partition column `chrom` MUST take values from the set `{chr1, chr2, ..., chr22, chrX, chrY, chrM}`. Variants on alternate, unplaced, or decoy contigs MUST NOT appear in the bundle and MUST be reported in the conversion log.

### 5.4 `pharmacogenomics.parquet` (OPTIONAL)

Single Parquet file containing gene-level pharmacogenomic calls. See §6.2 for schema. If absent, the bundle makes no pharmacogenomic claims; readers MUST NOT infer phenotypes from `variants.parquet` rows alone.

### 5.5 `gene_index.parquet` (OPTIONAL but RECOMMENDED)

Single Parquet file with one row per gene appearing in `variants.parquet`. See §6.3 for schema. Used by readers to plan gene-scoped queries efficiently.

### 5.6 `prs.parquet` (OPTIONAL)

Single Parquet file containing polygenic scores. See §6.4 for schema. Rows in `prs.parquet` are exploratory by default; their `clinical_grade` is determined per-row by the rule in `manifest.json.consensus_definitions.clinical_grade`.

### 5.7 `findings.json` (OPTIONAL)

JSON document containing a resolved, partner-facing summary of clinically significant findings. The schema of this file is API-surface, not standard-bound, and is OPTIONAL. Future revisions MAY normalise this file's structure.

### 5.8 `report.pdf` (OPTIONAL)

A human-readable PDF report. Out of scope for this specification.

### 5.9 `reads.cram` and `reads.cram.crai` (OPTIONAL)

Aligned sequencing reads in CRAM format with index. Provided as a sovereignty layer for downstream re-analysis. If `reads.cram` is present, `reads.cram.crai` MUST also be present.

### 5.10 Evidence, interpretation, and callability tables (OPTIONAL)

`clinical_findings.parquet`, `clinical_evidence.parquet`, `interpretations.parquet`, `gwas_associations.parquet`, `callability.parquet`, and `callable_regions.parquet` are additive v1.1 tables defined in §6.5–§6.10. If present, every file MUST follow its declared schema and MUST be listed in `manifest.json.files`.

Their absence means the bundle makes no claim in that evidence domain. It MUST NOT be interpreted as a negative clinical, genetic, pharmacogenomic, or trait result.

### 5.11 Self-description and provenance documents (OPTIONAL)

`field_glossary.json`, `inventory.json`, `raw_artifacts.json`, `transformations.json`, `omissions.json`, and `relationships.json` are defined in §7.5–§7.10. They are OPTIONAL for core v1.1 conformance and REQUIRED by the Agent Exploration Profile (§12).

### 5.12 Agent guidance (OPTIONAL)

`readmygenome.md` is a short human- and agent-readable entry point. `readmygenome/` MAY contain an operating guide, read-only query tools, a query catalogue, synonyms, and evaluations. Core readers MUST treat these as non-authoritative operational aids: genomic claims remain authoritative only when grounded in the manifest-listed data and evidence files.

---

## 6. Table schemas

This section specifies the column schema of each Parquet table. The authoritative schema, including JSON Schema validation rules, is in `schema.json` (§5.2). This section is a human-readable rendering for implementers; in case of discrepancy, `schema.json` governs.

All Parquet files MUST be Parquet 2.6 or later, with zstd compression. All Parquet files SHOULD include row-group statistics on `pos` (where applicable) and SHOULD include bloom filters on `gene.symbol`, `rsid`, and `pathogenicity.clinvar_id`.

### 6.1 `variants.parquet`

One row per (sample, normalised variant). v1.0 is single-sample.

#### 6.1.1 Identification

| Column | Type | Required | Description |
|---|---|---|---|
| `variant_id` | string | yes | Normalised identifier, format `chr{N}:{pos}:{ref}:{alt}`. See §8. |
| `chrom` | dictionary<string> | yes | Chromosome (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`). Partition column. |
| `pos` | int32 | yes | 1-based position on GRCh38. |
| `ref` | string | yes | Reference allele, IUPAC nucleotide string. |
| `alt` | string | yes | Alternate allele. Always biallelic post-normalisation. |
| `variant_type` | dictionary<string> | yes | One of `SNV`, `INS`, `DEL`, `MNV`, `SV`. |
| `rsid` | string (nullable) | no | dbSNP rsID when assigned. |

#### 6.1.2 `quality` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `call_confidence` | dictionary<string> | yes | `high`, `medium`, `low` — derived from DP and GQ thresholds in manifest. Reflects technical quality only. |
| `origin` | dictionary<string> | yes | `observed` or `imputed`. |
| `imputation_r2` | float32 (nullable) | no | Imputation accuracy (0–1). Populated iff `origin = imputed`. |
| `imputation_panel` | dictionary<string> (nullable) | no | Reference panel name. |
| `dp` | int16 (nullable) | no | Read depth. |
| `gq` | int8 (nullable) | no | Genotype quality. |

**Conformance:** `call_confidence`, `origin`, and `imputation_r2` MUST NOT be conflated into a single tier value. Readers MUST NOT treat `imputed` as equivalent to low-quality; high-R² imputed calls are often higher confidence than marginal directly-genotyped calls.

#### 6.1.3 `genotype` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `gt` | list<int8> | yes | Allele indices, e.g. `[0, 1]` for het, `[1, 1]` for hom alt. |
| `phased` | bool | yes | Whether the genotype is phased. |
| `zygosity` | dictionary<string> | yes | `hom_ref`, `het`, `hom_alt`. |
| `allele_balance` | float32 (nullable) | no | Allele balance for hets (sanity check). |

#### 6.1.4 `gene` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `symbol` | dictionary<string> (nullable) | no | HGNC gene symbol. Null for intergenic. |
| `ensembl_id` | string (nullable) | no | Ensembl gene ID. |
| `hgnc_id` | string (nullable) | no | HGNC ID. |
| `biotype` | dictionary<string> (nullable) | no | e.g. `protein_coding`, `lncRNA`. |
| `is_canonical` | bool (nullable) | no | True iff the row's `consequence` uses the canonical transcript. |

#### 6.1.5 `consequence` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `most_severe` | dictionary<string> (nullable) | no | VEP most-severe-consequence call. |
| `impact` | dictionary<string> (nullable) | no | `HIGH`, `MODERATE`, `LOW`, `MODIFIER`. |
| `hgvsc` | string (nullable) | no | Coding-sequence notation. |
| `hgvsp` | string (nullable) | no | Protein-level notation. |
| `all_annotations` | list<struct> (nullable) | no | Per-transcript annotations. |

**Conformance:** VEP CSQ pipe-delimited strings from upstream pipelines MUST be parsed at bundle generation time and decomposed into typed columns. Pipe-delimited strings MUST NOT appear in the bundle.

#### 6.1.6 `population_frequency` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `gnomad_af` | float32 (nullable) | no | Overall gnomAD allele frequency. |
| `gnomad_af_popmax` | float32 (nullable) | no | Highest frequency across major populations. |
| `gnomad_af_nfe`, `_eas`, `_afr`, `_sas`, `_amr` | float32 (nullable) | no | Per-population frequencies. |
| `is_rare` | bool | yes | Pre-computed per the rule in manifest. See §9.1. |

#### 6.1.7 `pathogenicity` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `clinvar_significance` | dictionary<string> (nullable) | no | One of: `Pathogenic`, `Likely_pathogenic`, `Uncertain_significance`, `Likely_benign`, `Benign`, `Conflicting_interpretations`. |
| `clinvar_has_conflicts` | bool | yes | True iff ClinVar submissions disagree. |
| `clinvar_conflict_summary` | dictionary<string> (nullable) | no | Compact description of conflict, e.g. `P_vs_LB`. |
| `clinvar_review_stars` | int8 (nullable) | no | ClinVar review confidence (0–4). |
| `clinvar_submitters_count` | int16 (nullable) | no | Number of ClinVar submissions. |
| `clinvar_id` | string (nullable) | no | ClinVar variant ID for cross-reference. |
| `is_actionable` | bool | yes | Pre-computed per the rule in manifest. See §9.2. |

**Conformance:** Significance and conflict status MUST be separate fields. A `Pathogenic` value with `clinvar_has_conflicts = true` MUST NOT be reported by a conformant reader as equivalent to a `Pathogenic` value with `clinvar_has_conflicts = false`.

#### 6.1.8 `predictions` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `cadd_phred` | float32 (nullable) | no | CADD PHRED-scaled score. |
| `revel` | float32 (nullable) | no | REVEL ensemble score (0–1). |
| `alphamissense` | float32 (nullable) | no | AlphaMissense score (0–1). |
| `spliceai_max` | float32 (nullable) | no | Max SpliceAI score across four consequences. |
| `polyphen` | float32 (nullable) | no | PolyPhen-2 score. |
| `sift` | float32 (nullable) | no | SIFT score. |
| `consensus_damaging` | bool | yes | Pre-computed per the rule in manifest. See §9.3. |

#### 6.1.9 `pharmacogenomics` (struct, variant-level)

| Field | Type | Required | Description |
|---|---|---|---|
| `is_pgx` | bool | yes | Whether the variant is part of a PGx-relevant gene. |
| `contributes_to_gene` | dictionary<string> (nullable) | no | Gene symbol. |

**Conformance:** Variant-level PGx fields MUST NOT include star alleles, phenotypes, drug recommendations, or CPIC levels. Those live exclusively in `pharmacogenomics.parquet` (§6.2). Readers MUST NOT derive phenotype claims from variant-level rows alone.

#### 6.1.10 `trait_associations` (struct)

| Field | Type | Required | Description |
|---|---|---|---|
| `is_gwas_hit` | bool | yes | Whether the variant has any GWAS association. |
| `traits` | list<string> (nullable) | no | Trait names. |
| `effect_size` | float32 (nullable) | no | Effect magnitude. |
| `effect_type` | dictionary<string> (nullable) | no | `beta`, `odds_ratio`, or `hazard_ratio`. |
| `effect_allele` | string (nullable) | no | Allele the effect is measured against. MUST equal `ref` or `alt`. |
| `ci_lower`, `ci_upper` | float32 (nullable) | no | 95% CI bounds. |
| `p_value` | float32 (nullable) | no | Association p-value. |
| `sample_size` | int32 (nullable) | no | GWAS sample size. |
| `study_pmids` | list<string> (nullable) | no | PubMed IDs of source studies. |

**Conformance (I-3):** If `effect_size` is non-null, `effect_type` AND `effect_allele` MUST also be non-null. A row containing `effect_size` without both `effect_type` and `effect_allele` is non-conformant. Bundle generators MUST reject input where these conditions cannot be satisfied; conformant readers MUST treat such rows as malformed and refuse to report the effect.

#### 6.1.11 Top-level row fields

| Column | Type | Required | Description |
|---|---|---|---|
| `clinical_grade` | bool | yes | Pre-computed per the rule in manifest. See §9.4. |

### 6.2 `pharmacogenomics.parquet`

One row per (gene, diplotype). Source of truth for pharmacogenomic phenotypes.

| Column | Type | Required | Description |
|---|---|---|---|
| `gene_symbol` | dictionary<string> | yes | HGNC gene symbol. |
| `diplotype` | string | yes | Two-haplotype call, e.g. `*1/*4`. |
| `phenotype` | dictionary<string> | yes | e.g. `intermediate metaboliser`, `poor metaboliser`. |
| `activity_score` | float32 (nullable) | no | CPIC activity score where applicable. |
| `copy_number` | int8 (nullable) | no | Diploid copy number for genes where CNV affects function. |
| `cpic_level` | dictionary<string> | yes | `A`, `B`, `C`, or `D`. |
| `affected_drugs` | list<string> | yes | Drug names with CPIC guidance for this gene. |
| `guideline_url` | string (nullable) | no | CPIC guideline URL. |
| `caller` | dictionary<string> | yes | Name of the haplotype caller used (e.g. `pharmcat-2.10`). |
| `caller_version` | string | yes | Caller version. |

### 6.3 `gene_index.parquet`

One row per gene appearing in `variants.parquet`.

| Column | Type | Required | Description |
|---|---|---|---|
| `gene_symbol` | dictionary<string> | yes | HGNC gene symbol. |
| `chrom` | dictionary<string> | yes | Chromosome. |
| `start_pos` | int32 | yes | Lowest variant position in this gene. |
| `end_pos` | int32 | yes | Highest variant position in this gene. |
| `variant_count` | int32 | yes | Total variants in this gene. |
| `actionable_count` | int32 | yes | Variants meeting `is_actionable`. |
| `pgx_count` | int32 | yes | Variants where `pharmacogenomics.is_pgx = true`. |

### 6.4 `prs.parquet`

One row per scored trait.

| Column | Type | Required | Description |
|---|---|---|---|
| `trait` | string | yes | Trait name. |
| `score_value` | float32 | yes | Raw polygenic score. |
| `percentile` | float32 | yes | Position in reference population (0–100). |
| `reference_population` | dictionary<string> | yes | Population the percentile is calibrated against. |
| `training_source` | string | yes | Source identifier (e.g. PGS Catalog ID `PGS000001`). |
| `training_date` | date | yes | Date of source summary statistics. |
| `contributing_variants` | list<string> (nullable) | no | List of `variant_id`s that contributed. |
| `model_type` | dictionary<string> | yes | e.g. `linear`, `lassosum`, `prscs`. |

### 6.5 `clinical_findings.parquet`

One row per pipeline-selected finding. This table is a retrieval index, not proof of expert review.

| Column | Type | Required | Description |
|---|---|---|---|
| `finding_id` | string | yes | Stable identifier within the bundle. |
| `variant_id` | string (nullable) | no | Related variant identifier, when the finding is variant-level. |
| `gene_symbol` | string (nullable) | no | HGNC symbol. |
| `condition` | string | yes | Normalised condition label. |
| `claim_type` | dictionary<string> | yes | One of the claim types in §6.9. |
| `classification` | string (nullable) | no | Source classification; null means the source did not provide one. |
| `clinical_grade` | bool | yes | Computed under the manifest rule. |
| `evidence_ids` | list<string> | yes | Related `clinical_evidence.evidence_id` values. Empty means no linked evidence rows were retained. |

### 6.6 `clinical_evidence.parquet`

One row per retained evidence assertion supporting, qualifying, or conflicting with a finding or interpretation.

| Column | Type | Required | Description |
|---|---|---|---|
| `evidence_id` | string | yes | Stable identifier within the bundle. |
| `finding_id` | string (nullable) | no | Related finding. |
| `variant_id` | string (nullable) | no | Related variant. |
| `source` | string | yes | Source name, e.g. `ClinVar`, `CPIC`, or `ACMG_SF`. |
| `source_record_id` | string (nullable) | no | Identifier in the source system. |
| `source_version` | string | yes | Version or snapshot date. |
| `assertion` | string | yes | Normalised assertion represented by this row. |
| `review_status` | string (nullable) | no | Source-specific review status; null means unavailable. |
| `retrieved_at` | timestamp | yes | Evidence snapshot time. |

### 6.7 `callability.parquet`

One row per queried site or panel target for which the generator retained locus-specific analytical evidence.

| Column | Type | Required | Description |
|---|---|---|---|
| `callability_id` | string | yes | Stable identifier within the bundle. |
| `variant_id` | string (nullable) | no | Exact variant identifier when the target is allele-specific. |
| `gene_symbol` | string (nullable) | no | HGNC symbol when assigned. |
| `chrom` | dictionary<string> | yes | GRCh38 chromosome. |
| `pos` | int32 | yes | One-based GRCh38 position. |
| `reference_observed` | bool (nullable) | no | True only when an explicit reference call was retained; null means not established. |
| `callable` | bool (nullable) | no | Whether the pipeline judged the locus callable; null means not assessed. |
| `call_confidence` | dictionary<string> (nullable) | no | `high`, `medium`, `low`, or null when not assessed. |
| `depth` | int32 (nullable) | no | Retained read depth. |
| `genotype_quality` | float32 (nullable) | no | Retained genotype quality. |
| `assay_scope` | string | yes | Scope that produced the callability evidence. |
| `evidence_scope` | dictionary<string> | yes | `site`, `target_panel`, or `derived_interval`. |

A row with `callable = true` but `reference_observed = null` does not establish a reference genotype. A missing row has no negative meaning.

### 6.8 `callable_regions.parquet`

One row per contiguous interval declared callable by the originating assay and pipeline.

| Column | Type | Required | Description |
|---|---|---|---|
| `region_id` | string | yes | Stable identifier within the bundle. |
| `chrom` | dictionary<string> | yes | GRCh38 chromosome. |
| `start_pos` | int32 | yes | One-based inclusive start. |
| `end_pos` | int32 | yes | One-based inclusive end. |
| `callable` | bool | yes | Interval-level pipeline assessment. |
| `minimum_depth` | int32 (nullable) | no | Applied minimum-depth threshold. |
| `assay_scope` | string | yes | Assay or processing scope. |
| `derivation` | string | yes | Named rule or pipeline stage that produced the interval. |

Interval callability can support locus reasoning only when the queried locus lies within the interval and the assay scope is applicable. It does not itself prove a homozygous-reference genotype.

### 6.9 `interpretations.parquet`

One row per structured interpretation. Interpretation rows make review status and model matching explicit rather than requiring a reader to infer clinical meaning from booleans.

| Column | Type | Required | Description |
|---|---|---|---|
| `interpretation_id` | string | yes | Stable identifier within the bundle. |
| `finding_id` | string (nullable) | no | Related finding. |
| `variant_id` | string (nullable) | no | Related variant. |
| `gene_symbol` | string (nullable) | no | HGNC symbol. |
| `interpretation_status` | dictionary<string> | yes | `algorithmic`, `draft`, `expert_reviewed`, or `guideline_based`. |
| `claim_type` | dictionary<string> | yes | `variant_classification`, `carrier_status`, `disease_risk`, `pharmacogenomic`, `trait_association`, or `computational_prediction`. |
| `condition` | string (nullable) | no | Normalised condition or trait. |
| `inheritance_model` | string (nullable) | no | Applied inheritance model. |
| `required_genotype` | string (nullable) | no | Genotype required by the claim model. |
| `observed_genotype` | string (nullable) | no | Bundle-observed genotype used by the interpretation. |
| `genotype_matches_model` | bool (nullable) | no | Explicit model-match result; null means not evaluated. |
| `penetrance` | string (nullable) | no | Source-described penetrance, including uncertainty. |
| `effect_size` | float64 (nullable) | no | Effect estimate; when populated, its type and effect allele MUST be linked through evidence. |
| `analytical_confidence` | dictionary<string> (nullable) | no | `high`, `medium`, `low`, or null when not assessed. |
| `classification_review_status` | string (nullable) | no | Review status of the classification source. |
| `gene_disease_validity` | string (nullable) | no | Source-described gene–disease validity. |
| `actionability_basis` | string (nullable) | no | Guideline or rule supporting actionability. |
| `confirmation_required` | bool | yes | Whether confirmatory testing is required before clinical use. |
| `reviewed_by` | string (nullable) | no | Reviewer identity or role; null for unreviewed records. |
| `reviewed_at` | timestamp (nullable) | no | Review time; null for unreviewed records. |
| `evidence_snapshot` | string | yes | Version or date of the evidence set. |
| `evidence_ids` | list<string> | yes | Related evidence identifiers. |

`expert_reviewed` MUST NOT be used unless `reviewed_by` and `reviewed_at` are populated. `guideline_based` MUST name the guideline in `actionability_basis` or linked evidence. `algorithmic` MUST NOT be presented as expert review.

### 6.10 `gwas_associations.parquet`

One row per retained study association. Rows represent external population evidence, not a personalised diagnosis or treatment recommendation.

| Column | Type | Required | Description |
|---|---|---|---|
| `association_id` | string | yes | Stable identifier within the bundle. |
| `variant_id` | string | yes | Related variant. |
| `trait` | string | yes | Normalised study trait. |
| `effect_allele` | string | yes | Allele to which the effect applies. |
| `effect_size` | float64 | yes | Study effect estimate. |
| `effect_type` | dictionary<string> | yes | `beta`, `odds_ratio`, or `hazard_ratio`. |
| `p_value` | float64 (nullable) | no | Study p-value. |
| `sample_size` | int64 (nullable) | no | Study sample size. |
| `study_pmids` | list<string> | yes | Supporting publications. |
| `source_version` | string | yes | Catalogue snapshot. |

---

## 7. Manifest and self-description documents

`manifest.json` is the reproducibility anchor. It MUST be present in every bundle and MUST contain at minimum:

```json
{
  "schema_version": "1.1.0",
  "pipeline_version": "<string identifying the bundle generator>",
  "genome_build": "GRCh38",
  "generated_at": "<ISO 8601 datetime>",
  "annotation_sources": {
    "vep": "<version>",
    "clinvar": "<version or date>",
    "gnomad": "<version>",
    "cadd": "<version>",
    "revel": "<version>",
    "alphamissense": "<version>",
    "spliceai": "<version>",
    "acmg_sf_list": "<version>",
    "cpic_guidelines": "<version or date>"
  },
  "consensus_definitions": {
    "is_rare": { ... },
    "is_actionable": { ... },
    "consensus_damaging": { ... },
    "clinical_grade": { ... }
  },
  "normalization": "<exact normalisation command, e.g. 'bcftools norm -f GRCh38.fa --atomize'>",
  "files": {
    "<relative_path>": {
      "sha256": "<64 hex chars>",
      "rows": <integer or null>,
      "bytes": <integer or null>
    }
  }
}
```

### 7.1 Required keys

`schema_version`, `pipeline_version`, `genome_build`, `generated_at`, `annotation_sources`, `consensus_definitions`, `normalization`, and `files` are **REQUIRED**. Bundles missing any of these are non-conformant.

### 7.2 `consensus_definitions`

Each rule MUST be a JSON object containing at minimum:

- `rule`: a human-readable string describing the rule
- The thresholds, source identifiers, and (where applicable) hashes used by the rule

See §9 for the required content of each rule object.

### 7.3 `files`

The `files` object MUST contain a SHA256 hash for every file in the bundle (excluding `manifest.json` itself). For `variants.parquet/` (a directory), the hash is computed by hashing every contained file in sorted relative-path order, separated by null bytes (the algorithm is implemented in the reference reader).

Conformant readers MUST verify file hashes before treating a bundle as authoritative.

### 7.4 Optional keys

Bundles MAY include `source_bundle_sha256` for an immutable migration lineage and `profiles.agent_exploration` for the profile version declared in §12. Bundles MAY include additional manifest keys for provider-specific metadata. Conformant readers MUST ignore unknown keys.

### 7.5 `field_glossary.json`

Machine-readable documentation for every field present in every retained Parquet table. The document MUST map each `<table>.<field_path>` to an object containing:

- `type` and whether the field is nullable
- `null_meaning`
- `source` and `derivation`
- `allowed_values`, or null when unconstrained
- `evidence_class`
- `observed_or_inferred`
- applicable `genome_build`
- `version_or_snapshot`
- `confidence_semantics`
- `absence_meaning`
- `related_evidence`

If the Agent Exploration Profile is declared, every queryable Parquet field MUST have exactly one glossary entry. “See README” or equivalent out-of-band documentation is not sufficient.

### 7.6 `inventory.json`

Complete logical inventory of the bundle. Each entry MUST include relative path, media type, logical role, required/optional status, and row count where applicable. It MUST also state whether the entry is authoritative data, derived data, metadata, operational guidance, or evaluation material. File integrity remains authoritative in `manifest.json.files`; `inventory.json` MUST NOT attempt a self-referential hash of itself or `manifest.json`.

### 7.7 `raw_artifacts.json`

Records the disposition of every source or sovereignty artifact used to create the bundle. Allowed dispositions are `embedded`, `separate`, `referenced`, `deleted`, and `never_received`. Each record MUST state artifact type, disposition, reason, integrity identifier when retained, and whether its absence limits re-analysis, callability, or confirmation.

This document describes disposition; it MUST NOT imply that a deleted or unreceived artifact was analysed beyond the recorded pipeline inputs.

### 7.8 `transformations.json`

Ordered pipeline-stage records. Each stage MUST include a stable stage identifier, software name and version, parameters or configuration digest, input identifiers, output identifiers, start/completion time when available, and status. External reference snapshots used by a stage MUST be identified.

### 7.9 `omissions.json`

Records source fields, record classes, contigs, annotations, or artifacts deliberately not retained. Each omission MUST include scope, reason, pipeline stage, consequence for interpretation, and whether recovery requires a retained artifact, re-annotation, or a new source upload.

An omission record is not a genetic finding and MUST NOT be interpreted as absence of the omitted biological feature.

### 7.10 `relationships.json`

Declares supported joins among tables and metadata. Each relationship MUST identify source table and fields, target table and fields, cardinality, whether the join is complete or partial, and the meaning of an unmatched row. Recommended joins include variant-to-finding, finding-to-evidence, interpretation-to-evidence, PGx-to-guideline evidence, PRS-to-source, and callability-to-variant or interval.

---

## 8. Variant identifier and normalisation

### 8.1 Identifier format

`variant_id` MUST be the string `{chrom}:{pos}:{ref}:{alt}` where:

- `chrom` ∈ `{chr1, ..., chr22, chrX, chrY, chrM}` (UCSC style)
- `pos` is 1-based
- `ref` and `alt` are uppercase IUPAC nucleotide strings (`A`, `C`, `G`, `T`)

### 8.2 Normalisation

Every variant MUST be left-normalised and atomised. The reference normalisation command is:

```
bcftools norm -f GRCh38.fa --atomize
```

The `manifest.json.normalization` field MUST record the exact command used. Equivalent normalisation by other tools is acceptable provided the resulting `(pos, ref, alt)` triples are identical to the bcftools output for the same input.

### 8.3 Multiallelic sites

Multiallelic sites MUST be split into biallelic rows. After splitting, each variant has exactly one `alt` allele.

### 8.4 Excluded variants

Variants on alternate, unplaced, or decoy contigs MUST NOT appear in the bundle. The conversion log SHOULD record the count of excluded variants by exclusion reason.

---

## 9. Pre-computed interpretive fields

Every pre-computed field MUST be defined by a rule in `manifest.json.consensus_definitions`. The rule MUST be sufficient to reproduce the field's value from the bundle's other contents.

### 9.1 `is_rare`

Required rule object structure:

```json
{
  "rule": "<human-readable description>",
  "threshold": <number, e.g. 0.01>,
  "population": "<population key, e.g. 'gnomad_af_popmax'>",
  "null_treatment": "<one of: 'rare', 'not_rare', 'null'>"
}
```

Default rule (RECOMMENDED but not REQUIRED):

```json
{
  "rule": "gnomad_af_popmax < 0.01; null AF treated as rare",
  "threshold": 0.01,
  "population": "gnomad_af_popmax",
  "null_treatment": "rare"
}
```

### 9.2 `is_actionable`

Required rule object structure:

```json
{
  "rule": "<human-readable description>",
  "gene_list": "<identifier, e.g. 'ACMG_SF_v3.2'>",
  "gene_list_sha256": "<64 hex chars>",
  "significance_set": ["<allowed clinvar_significance values>"],
  "min_review_stars": <integer, 0-4>
}
```

Default rule (RECOMMENDED):

```json
{
  "rule": "gene ∈ ACMG SF v3.2 AND clinvar_significance ∈ {Pathogenic, Likely_pathogenic} AND clinvar_review_stars >= 2",
  "gene_list": "ACMG_SF_v3.2",
  "gene_list_sha256": "<hash of the canonical 78-gene list>",
  "significance_set": ["Pathogenic", "Likely_pathogenic"],
  "min_review_stars": 2
}
```

The `gene_list_sha256` MUST be the SHA256 of a canonical text representation of the gene list (one symbol per line, sorted alphabetically, trailing newline). Bundle generators MUST publish the canonical form of any non-standard gene list used.

### 9.3 `consensus_damaging`

Required rule object structure:

```json
{
  "rule": "<human-readable description>",
  "applies_to": ["<consequence values this rule is computed for>"],
  "tools": [
    {"name": "<tool>", "field": "<field path>", "threshold": <number>}
  ],
  "aggregation": "<one of: 'at_least_n', 'all', 'majority'>",
  "n": <integer, required iff aggregation = 'at_least_n'>
}
```

Default rule (RECOMMENDED):

```json
{
  "rule": "≥2 of REVEL/AlphaMissense/CADD exceed thresholds, applied to missense variants only",
  "applies_to": ["missense_variant"],
  "tools": [
    {"name": "REVEL", "field": "predictions.revel", "threshold": 0.75},
    {"name": "AlphaMissense", "field": "predictions.alphamissense", "threshold": 0.564},
    {"name": "CADD", "field": "predictions.cadd_phred", "threshold": 20.0}
  ],
  "aggregation": "at_least_n",
  "n": 2
}
```

### 9.4 `clinical_grade`

Required rule object structure:

```json
{
  "rule": "<human-readable description>",
  "qualifying_sources": [
    {"source": "<identifier>", "minimum_evidence": "<criterion>"}
  ]
}
```

Default rule (RECOMMENDED):

```json
{
  "rule": "true iff backed by ACMG SF, CPIC, or ClinVar at >=2 review stars",
  "qualifying_sources": [
    {"source": "ACMG_SF_v3.2", "minimum_evidence": "gene_membership"},
    {"source": "CPIC", "minimum_evidence": "level_A_or_B"},
    {"source": "ClinVar", "minimum_evidence": "review_stars>=2"}
  ]
}
```

PRS-derived rows in `prs.parquet` and GWAS-only signals in `trait_associations` are NOT clinical-grade by default.

### 9.5 Conformance

A bundle whose `manifest.json.consensus_definitions` does not match the actual computed values in the Parquet tables is non-conformant. Bundle generators SHOULD verify rule conformance before emitting `manifest.json`.

---

## 10. Versioning and compatibility

### 10.1 Semantic versioning

`schema_version` follows semantic versioning:

- **Patch** (1.0.x) — clarifications, typo fixes, additional non-normative examples. Readers and bundles remain compatible.
- **Minor** (1.x.0) — additive changes (new optional columns, new optional files). Readers MUST ignore unknown columns and files; bundles MUST NOT remove or change existing fields.
- **Major** (x.0.0) — breaking changes (removed fields, semantic changes, new requirements). Readers MUST check `schema_version` and reject bundles whose major version they do not support.

### 10.2 Forward compatibility

A v1.0 reader MUST accept v1.x bundles where x > 0 by ignoring unknown columns and files. A v1.0 reader MAY reject v2.x or later bundles and SHOULD do so explicitly rather than misparse.

An Agent Exploration Profile declaration does not change the core major version. Readers that do not implement the profile MUST continue to read the core v1 tables and ignore profile assets.

### 10.3 Recomputability boundary

The bundle distinguishes durable variant-level data from dated interpretive content:

- **Inside the boundary** (durable): `variant_id`, `chrom`, `pos`, `ref`, `alt`, `genotype`, raw values within `quality`, raw values within `population_frequency`, raw scores within `predictions`, raw `clinvar_*` fields, and retained assay-level callability observations.
- **Outside the boundary** (dated snapshot): all pre-computed booleans (`is_rare`, `is_actionable`, `consensus_damaging`, `clinical_grade`), `pharmacogenomics.parquet` calls, `prs.parquet` scores, external evidence joins, and structured interpretations.

Re-interpretation is the operation of regenerating outside-boundary content from inside-boundary content under updated rules. Conformant bundle generators SHOULD support re-interpretation as a primary operation. Conformant readers MAY surface the snapshot date and offer the user a choice between snapshot-consistent and latest-consistent interpretation modes.

---

## 11. Reader conformance

A conformant `.genome/1.x` reader MUST:

1. **Validate before reading.** Verify the bundle directory exists, required files are present, `manifest.json` parses, and the declared `schema_version` is supported.
2. **Verify file hashes.** Compute SHA256 of every file declared in `manifest.json.files` and confirm matches before treating bundle contents as authoritative.
3. **Reject non-conformant bundles.** Refuse to operate on bundles missing required fields, with hash mismatches, or with violated invariants. MUST NOT silently coerce.
4. **Cite rule provenance.** When surfacing the value of any pre-computed field to a user, include the rule definition from `manifest.json.consensus_definitions` and the source version it was computed against.
5. **Surface snapshot date.** When surfacing any interpretive value, make the bundle's `generated_at` available to the user.
6. **Distinguish bundle claims from external knowledge.** When a reader supplements bundle content with external information (e.g. general knowledge about a gene), the reader MUST make the distinction visible to the user.
7. **Respect the recomputability boundary.** If the reader supports re-interpretation, it MUST NOT mutate the original bundle; re-interpretation produces a new bundle with new identifiers.
8. **Treat missing rows as retrieval outcomes.** A reader MUST NOT turn an empty query result into a negative genetic claim without applicable explicit genotype or callability evidence.
9. **Preserve evidence classes.** A reader MUST distinguish observed data, pipeline-derived values, external annotations, and reviewed interpretations in user-facing answers.
10. **Bound queries.** Agent-oriented readers MUST use projected, filtered, bounded queries and MUST NOT load complete Parquet tables into language-model context.

A conformant reader SHOULD:

- Surface ClinVar conflict status (`clinvar_has_conflicts`) prominently when reporting `clinvar_significance`.
- Default clinical-context queries to `clinical_grade = true` filtering.
- Warn the user when the bundle is materially old (RECOMMENDED threshold: 180 days) and the question depends on rules that change frequently.
- Refuse to derive pharmacogenomic phenotype claims from variant-level data.

---

## 12. Agent Exploration Profile

The Agent Exploration Profile is an optional conformance profile for a `.genome/1.1` bundle intended to be downloaded and explored directly with Codex, Claude Code, or another general-purpose software agent. A bundle declares the profile by setting `manifest.json.profiles.agent_exploration` to a semantic version.

### 12.1 Required profile assets

A conformant Agent Exploration Profile bundle MUST contain:

- `field_glossary.json`, `inventory.json`, `raw_artifacts.json`, `transformations.json`, `omissions.json`, and `relationships.json`
- `readmygenome.md`
- `readmygenome/SKILL.md`
- a read-only query CLI with JSON output by default
- a tested query catalogue
- a synonym resource for conditions, drugs, traits, and genes
- an expected-answer evaluation suite

If structured interpretations are present, they MUST use `interpretations.parquet`. If the bundle can make any locus-specific absence statement, it MUST include applicable `callability.parquet`, `callable_regions.parquet`, or explicit reference genotype evidence.

### 12.2 Operating guide requirements

The bundled operating guide MUST direct an agent to:

1. Validate integrity and orient to manifest, inventory, schema, and evidence scope before substantive querying.
2. Never load whole Parquet tables into model context.
3. Use bounded, iterative, read-only DuckDB queries.
4. Search all relevant tables rather than treating a single empty lookup as final.
5. Check callability before interpreting absence.
6. Separate observed data, pipeline-derived values, external annotations, reviewed interpretations, and general knowledge.
7. Continue with synonyms and related identifiers when an initial retrieval is empty.
8. Cite table, row identifier, source version, evidence snapshot, and interpretation or review date where applicable.
9. Explain findings directly while preserving confirmation and medical-action boundaries.

### 12.3 Query CLI requirements

The CLI MUST provide automatic bundle-root discovery and named operations for overview, free-text search, variant, gene, clinical condition, pharmacogenomics drug, polygenic trait, callability, and evidence lookup. It MUST support partitioned Parquet, emit JSON by default, report query counts and truncation, provide an optional human-readable table format, work without network access for basic exploration, and reject mutating SQL.

### 12.4 Query catalogue and evaluations

The catalogue MUST cover genome orientation, important findings, gene and variant lookup, disease and phenotype lookup, medication and PGx review, carrier status, polygenic traits, rare damaging variants, negative findings, callability, conflicting interpretations, ancestry limitations, comparison with new literature or guidelines, and a clinician-facing summary.

The evaluation suite MUST contain at least 30 canonical questions. Each case MUST identify expected tables searched, required findings, claims that must not be made, expected callability behaviour, and expected citation fields. It MUST include broad-query and empty-retrieval regression cases.

### 12.5 Profile conformance and authority

Operational guidance, synonym resources, and evaluation fixtures are non-authoritative. They MUST NOT override bundle data, evidence snapshots, the binding schema, or manifest rules. A profile-conformant agent that cannot execute the bundled tools MAY still explore the bundle, but MUST preserve the same bounded-query, evidence, absence, and citation rules.

---

## 13. Examples

### 13.1 Minimal conformant manifest

```json
{
  "schema_version": "1.1.0",
  "pipeline_version": "genome-convert-1.1.0",
  "genome_build": "GRCh38",
  "generated_at": "2026-04-20T14:32:00Z",
  "annotation_sources": {
    "vep": "111",
    "clinvar": "2026-04",
    "gnomad": "4.1",
    "cadd": "1.7",
    "revel": "1.3",
    "alphamissense": "1.0",
    "spliceai": "1.3",
    "acmg_sf_list": "v3.2",
    "cpic_guidelines": "2026-03"
  },
  "consensus_definitions": {
    "is_rare": {
      "rule": "gnomad_af_popmax < 0.01; null AF treated as rare",
      "threshold": 0.01,
      "population": "gnomad_af_popmax",
      "null_treatment": "rare"
    },
    "is_actionable": {
      "rule": "gene ∈ ACMG SF v3.2 AND clinvar_significance ∈ {Pathogenic, Likely_pathogenic} AND clinvar_review_stars >= 2",
      "gene_list": "ACMG_SF_v3.2",
      "gene_list_sha256": "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
      "significance_set": ["Pathogenic", "Likely_pathogenic"],
      "min_review_stars": 2
    },
    "consensus_damaging": {
      "rule": "≥2 of REVEL/AlphaMissense/CADD exceed thresholds, missense only",
      "applies_to": ["missense_variant"],
      "tools": [
        {"name": "REVEL", "field": "predictions.revel", "threshold": 0.75},
        {"name": "AlphaMissense", "field": "predictions.alphamissense", "threshold": 0.564},
        {"name": "CADD", "field": "predictions.cadd_phred", "threshold": 20.0}
      ],
      "aggregation": "at_least_n",
      "n": 2
    },
    "clinical_grade": {
      "rule": "true iff backed by ACMG SF, CPIC, or ClinVar at >=2 review stars",
      "qualifying_sources": [
        {"source": "ACMG_SF_v3.2", "minimum_evidence": "gene_membership"},
        {"source": "CPIC", "minimum_evidence": "level_A_or_B"},
        {"source": "ClinVar", "minimum_evidence": "review_stars>=2"}
      ]
    }
  },
  "normalization": "bcftools norm -f GRCh38.fa --atomize",
  "files": {
    "schema.json": {"sha256": "...", "bytes": 12345},
    "variants.parquet": {"sha256": "...", "rows": 4523891, "bytes": 480000000},
    "pharmacogenomics.parquet": {"sha256": "...", "rows": 47},
    "gene_index.parquet": {"sha256": "...", "rows": 19847},
    "prs.parquet": {"sha256": "...", "rows": 23}
  }
}
```

### 13.2 Canonical query: actionable findings

```sql
SELECT variant_id,
       gene.symbol AS gene,
       consequence.hgvsp,
       pathogenicity.clinvar_significance,
       pathogenicity.clinvar_review_stars,
       predictions.revel
FROM read_parquet('variants.parquet/**/*.parquet', hive_partitioning=true)
WHERE pathogenicity.is_actionable = true
  AND clinical_grade = true
ORDER BY predictions.revel DESC NULLS LAST;
```

This query returns in <500ms on a typical WGS bundle. It cannot return false positives from CSQ misparsing, ClinVar substring matching, or ACMG list approximation, because none of those failure modes is reachable through the schema.

---

## 14. Open issues for future versions

The following are explicitly out of scope for v1.1 and will be addressed in subsequent versions:

- **Multi-sample bundles.** v1.1 or v1.2.
- **Polygenic-score versioning across multiple training sources for the same trait.** v1.1.
- **Phasing serialisation across `genotype.gt`.** v1.1.
- **Primary-findings rule (rare-disease diagnostic context).** v1.1.
- **Structural-variant representation beyond the existing `SV` discriminator.** v1.2 or later.
- **GRCh37/hg19 support.** Not currently planned.

---

## 15. Acknowledgements

This specification incorporates feedback from external technical reviewers whose critiques shaped §3 (invariant framing), §10.3 (recomputability boundary), §6.1.10 (effect-allele requirement), and §11 (reader conformance). Specific contributions are recorded in `RATIONALE.md` in this repository.

---

## 16. Change log

### 1.1.0 (August 2026)

Additive evidence, interpretation, callability, provenance, omission, and relationship records; explicit missing-row semantics; and the optional Agent Exploration Profile. Existing v1.0 fields are unchanged.

### 1.0.0 (April 2026)

Initial public release. See `CHANGELOG.md` for the development history including v0.9 → v1.0 changes.

---

## 17. References

- RFC 2119 — Key words for use in RFCs to Indicate Requirement Levels
- VCF Specification v4.4 — https://samtools.github.io/hts-specs/
- Apache Parquet Specification — https://parquet.apache.org/docs/
- ClinVar — https://www.ncbi.nlm.nih.gov/clinvar/
- ACMG Secondary Findings v3.2 — Miller DT et al., *Genetics in Medicine* (2023)
- CPIC — https://cpicpgx.org/
- gnomAD v4.1 — https://gnomad.broadinstitute.org/
- PharmCAT — https://pharmcat.org/

---

*This specification is published at https://github.com/genome-computer/genome-spec under the Apache 2.0 license. Correspondence: research@genome.computer.*
