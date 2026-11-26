# Genome Annotation Builder

`genome-annotation-builder` is a reproducible, evidence-preserving pipeline for
turning a GRCh38 reference-confidence VCF/gVCF into a `.genome/1.1` bundle.
It treats an absent VCF row, a no-call, an interval judged callable, and an
explicit homozygous-reference observation as four different facts.

The project is designed for public use with arbitrary samples. It contains no
real genomes, sample identifiers, local filesystem paths, or derived personal
metrics. A default-deny `.gitignore` prevents common genomic formats and local
runtime directories from entering Git.

## Current scope

The first implementation milestone builds and validates the lossless
small-variant core:

- streams a single-sample GRCh38 reference-confidence VCF;
- separates called biological alleles from symbolic reference blocks and
  no-call evidence;
- decomposes called multiallelic genotypes without converting missing calls to
  reference;
- normalizes biological alleles against the exact configured GRCh38 FASTA;
- writes typed, chromosome-partitioned Parquet plus a compact reference-
  confidence ledger;
- optionally runs Ensembl VEP 116.0 entirely offline against the checksum-
  pinned GRCh38 cache, verifies that no allele or genotype changed, and parses
  every transcript annotation into typed structures;
- materializes a gene query index from all annotated transcript memberships;
  selects a deterministic MANE/canonical representative within each gene,
  chooses the primary gene by the representative consequence, and retains the
  independently most-severe consequence across every transcript;
- optionally materializes the immutable ClinGen 2026-08-15 gene-disease
  validity and GRCh38 dosage-sensitivity snapshots, then projects only exact
  HGNC-symbol matches into two gene-context extensions; these rows never enter
  variant-level clinical evidence and never imply pathogenicity or diagnosis;
- optionally installs a dated, checksum-pinned NCBI ClinVar GRCh38 snapshot,
  normalizes it against the configured FASTA, performs an exact-allele streaming
  join, and preserves review status and conflicting assertions separately;
- vendors a deterministic, logo-free ClinGen/ACMG SF v3.3 policy snapshot with
  all 84 genes, disease-level reportability, inheritance, reporting guidance,
  MANE transcripts, and frequency specifications;
- optionally evaluates disease-specific ACMG secondary-finding candidates with
  exact-allele ClinVar evidence, inheritance/genotype and gene-specific gates;
  every result remains algorithmic, `clinical_grade=false`, and requires
  confirmation;
- emits normative `clinical_evidence.parquet`, `clinical_findings.parquet`, and
  `interpretations.parquet` retrieval tables while retaining lossless ClinVar
  fields in `extensions/clinvar_assertions.parquet`;
- optionally consumes a checksum-pinned catalogue of actually queried loci and
  emits normative `callability.parquet` plus a lossless target-to-source proof
  extension; carried variants and no-calls always supersede reference blocks,
  and interval callability alone never becomes a reference genotype;
- optionally ingests one or more checksum-pinned PGS Catalog format-2.0
  GRCh38-harmonized additive models, orients every biallelic SNV against the
  exact FASTA, and scores only explicit passing variant/reference genotypes;
  per-locus and per-model audit extensions preserve unsupported, other-ALT,
  no-call, low-confidence, interval-only, and unknown states;
- emits normative `prs.parquet` only when the declared missing-weight gate
  passes and a digest-bound empirical calibration is applicable to its named
  reference population; an uncalibrated raw score remains audit-only;
- scores locus audits in bounded batches and aggregates from disk; models with
  more than 10,000 contributing variant IDs retain the exact count and full
  locus-level mapping while the optional embedded ID list is null;
- optionally prepares the exact PharmCAT 3.4.0 positions and a frozen CPIC
  pair-view snapshot, automatically adds all official alleles to the target
  callability union, and materializes a position-complete matcher VCF;
- runs the immutable PharmCAT container with no network, a read-only root,
  dropped capabilities, and none of its missing-to-reference or research-mode
  options; structurally complex genes are excluded from VCF-only matching;
- emits normative `pharmacogenomics.parquet` only for complete, unambiguous,
  non-complex calls with an exact non-provisional CPIC A-D mapping, while
  retaining every locus, withheld gene call, and structured guidance decision
  in typed audit extensions;
- optionally projects a digest-pinned normalized family joint VCF and PED into
  `extensions/family_evidence.parquet`; every row must preserve the accepted
  subject `GT` and phase exactly, while relative identifiers and genotypes are
  never persisted;
- optionally joins a child-independent LOCO exception ledger, retaining exact
  de-novo/conflict/uninformative rows but never interpreting absence from that
  exception-only file as site-level consistency;
- optionally freezes and normalizes one dated NHGRI-EBI GWAS Catalog release,
  then projects only exact GRCh38 coordinate + rsID-alias + effect-allele
  matches; typed associations and a lossless source-row audit remain separate;
- optionally installs checksum-pinned offline LOFTEE, CADD, REVEL,
  AlphaMissense, and SpliceAI resources, verifies every upstream and derived
  object, passes only validated paths to VEP, and records consequence-aware
  score completeness without treating a missing score as benign;
- optionally executes the separately staged, receipt-bound Release C for
  read-backed structural variants, copy-number variants, repeat expansions,
  HLA, complex pharmacogenomics, mitochondrial calls, and bounded local-read
  evidence;
- validates every Release C input immediately before its caller, executes an
  immutable profile without network access, retains failed attempts, and
  accepts output only through exact artifact, schema, reference, and row-count
  reconciliation;
- runs technical datasets sequentially, preserves their evidence separately,
  and computes domain-specific replicate concordance without majority voting,
  silent row removal, or interpreting a missing row as a negative result;
- assembles an atomic `.genome/1.1` directory with provenance and omissions;
- performs stronger row and cross-artifact validation than the upstream
  reference reader;
- exposes bounded, read-only queries suitable for the Agent Exploration
  Profile.

The base VEP adapter includes cache-provided dbSNP and gnomAD frequencies. The
ClinVar/ACMG, predictor, GWAS, and read-backed layers create versioned evidence
chains, not diagnoses or expert-reviewed reports. Callers or model classes
outside the explicit adapters require a separately checksum-pinned profile. A
missing source remains missing; it is never interpreted as benign, reference,
or an absent structural event.

## Read-backed Release C

Release C is intentionally a separate extension stage because it consumes an
indexed GRCh38 alignment rather than inferring structural or complex-locus
evidence from a small-variant VCF. The public implementation ships generic
caller profiles, typed adapters, sequential matrix orchestration, subject-stage
assembly, atomic bundle extension, validation, and bounded queries. It contains
no real genome, subject identifier, input path, or production result.

After bundle validation, `gab release-c project-release build` creates
identity-bound per-bundle and project-wide JSON/TSV/HTML annotation QC,
machine-readable bundle validation reports, manifests, an audit log, and
`SHA256SUMS`. QC covers reconciliation, variant/genotype classes, Ti/Tv, VEP,
ClinVar, gnomAD, applicable predictors, difficult regions, PGx/PGS,
callability, family evidence, technical-replicate concordance, and read-backed
categories. It exports aggregate counters only—never loci, alleles, genotypes,
private paths, or technical dataset identifiers—and uses a bounded-memory
Parquet scan suitable for whole-genome bundles. For cohorts of three or more,
fixed cross-bundle median rules stop publication on an unexplained QC outlier;
an investigated outlier can proceed only with a portable explanation supplied
from an ignored local file.

Start with the public examples in `config/release-c.example.yaml`,
`config/release-c-matrix.example.json`, and
`config/release-c-resource-metadata.example.json`; see
[`docs/release-c.md`](docs/release-c.md) for the fail-closed lifecycle and
resource-freeze procedure.

## Reproducibility boundary

The core format is pinned to genome.computer `genome-spec` commit
`7354ae288784234f854562b7dceb87f067fad550`. Normative files are fetched from
that immutable revision and accepted only if their SHA-256 values match the
compiled-in allow-list. Every build records input, reference, resource,
software, command, and output digests.

## Quick start

Requirements: Python 3.11+, `bcftools`, and a GRCh38 FASTA with `.fai`.
Docker is additionally required for enabled container adapters such as VEP or
PharmCAT. The optional Snakemake path uses Conda/Miniforge and the checked-in
`workflow/envs/gab.yaml` environment.

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[dev,workflow]'
cp config/config.example.yaml config/config.yaml
# Edit only config/config.yaml; it is ignored by Git.
gab resources sync --config config/config.yaml
# If annotation.vep.enabled is true, this verifies and installs the public
# Ensembl cache archive (approximately 27.6 GB compressed):
gab resources vep-sync --config config/config.yaml
# If any annotation.predictors adapter is enabled, review its source terms,
# download/prepare every exact release, and then independently re-verify it:
gab resources predictors-sync --config config/config.yaml
gab resources predictors-verify --config config/config.yaml
# Run fixed public positive controls through the same immutable VEP/plugins;
# preflight refuses enabled predictors until this semantic gate verifies:
gab resources predictors-doctor --config config/config.yaml
gab resources predictors-doctor-verify --config config/config.yaml
# If annotation.clinvar.enabled is true, this installs and exact-reference
# normalizes the pinned dated NCBI snapshot:
gab resources clinvar-sync --config config/config.yaml
# Install the vendored ACMG SF v3.3 / ClinGen 2026-08-15 snapshot. This works
# offline; with resources.allow_download=true, the exact source workbook is
# additionally archived and cross-checked.
gab resources acmg-sf-sync --config config/config.yaml
# If annotation.clingen.enabled is true, deterministically materialize and
# verify the vendored CC0 gene-validity and dosage-sensitivity indexes:
gab resources clingen-sync --config config/config.yaml
gab resources clingen-verify --config config/config.yaml
# If annotation.gwas.enabled is true, download/CRC-check the exact dated source
# archive, normalize it deterministically, and verify its Parquet + manifest:
gab resources gwas-sync --config config/config.yaml
# If annotation.pgs.enabled is true, verify every locally acquired scoring and
# calibration file against the configured digest and exact GRCh38 FASTA:
gab resources pgs-verify --config config/config.yaml
# If annotation.pgx.enabled is true, pull/verify only the exact configured OCI
# digest, install/verify the exact PharmCAT positions, and verify the already
# frozen, locally configured CPIC snapshot:
gab resources pharmcat-sync --config config/config.yaml
gab resources pgx-sync --config config/config.yaml
gab preflight --config config/config.yaml
gab build --config config/config.yaml
gab validate /path/to/output/sample.genome
gab query /path/to/output/sample.genome overview
```

The same restart-safe sequence can be driven by the lint-clean workflow. The
top-level rule intentionally has no declared output: it invokes `gab build` on
every workflow run, so an existing immutable bundle is revalidated instead of
being trusted merely because `manifest.json` exists.

```bash
conda config --set channel_priority strict
GAB_CONFIG=config/config.yaml snakemake \
  --snakefile workflow/Snakefile \
  --cores 1 \
  --use-conda
```

A build reuses an existing bundle only when its recorded build key matches the
current input, implementation, configuration and resources. To rebuild a legacy
bundle or run with a different configuration, choose a new output path or set
`output.overwrite` explicitly. The build validates the staged bundle before
promoting it with an atomic rename.

The example configuration pins the official VEP 116.0 Linux/amd64 container
by immutable OCI digest. On another architecture, use the matching immutable
release digest and retain the exact runtime identity in the build manifest.

## Multiple samples

The portable format and builder are deliberately single-sample: create one
immutable configuration and one bundle per subject. A cohort operator may run
those configurations sequentially to bound disk, RAM, and container pressure;
no cross-sample annotation concurrency is required for quality. Shared evidence
such as a normalized joint callset is consumed through an immutable,
checksum-bound adapter, while each portable bundle retains only the subject's
own genotype and privacy-minimized derived statuses. Restarting the same config
revalidates an existing final bundle or reuses digest-verified work stages.

## Offline prediction resources

All prediction adapters are disabled by default. Enabling an adapter requires
VEP 116.0 from the pinned Docker image, its built-in `/plugins` directory, the
exact licence/terms text compiled for the pinned release, and
`redistribution_allowed: false`. Executable-mode VEP remains available for the
base adapter, but predictor plugins are rejected there because this release
cannot bind arbitrary host plugin code to an immutable runtime identity. This
project distributes
the adapter code and immutable identities only; it does not distribute the
large databases or grant rights to use or redistribute them. Review each
upstream resource independently before downloading or publishing derived data.

`predictors-sync` uses restartable downloads, exact byte counts and SHA-256
checks (plus publisher MD5 where available), structurally validates indexed
inputs, and atomically writes a deterministic relative-path private manifest.
Local installation time is deliberately excluded from its content identity, so
the same exact upstream and derived bytes reproduce the same manifest SHA-256.
It also creates and verifies deterministic derivatives where the VEP plugin
requires them:

- LOFTEE is pinned to one source commit. Its upstream `LoF.pm` invokes an
  external `samtools faidx` that is absent from the immutable VEP image, so the
  builder applies the shipped `biodbhts-faidx-v1` patch with GNU patch 2.7.6.
  The source archive, patch, exact preimage, exact patched `LoF.pm`, complete
  source tree, ancestor FASTA, conservation bigWig, and read-only conservation
  SQLite database are all independently checksum-bound. The patch uses the
  image's existing `Bio::DB::HTS::Faidx`; no host `samtools` is mounted;
- CADD v1.7 uses separate GRCh38 SNV and small-indel tabix resources;
- REVEL v1.3 is filtered to rows with a GRCh38 coordinate, stably sorted,
  BGZF-compressed, transcript-preserving, and tabix-indexed;
- AlphaMissense uses transcript-matched GRCh38 scores and a locally generated
  tabix index. The exact hosted file carries a stricter CC BY-NC-SA 4.0 notice
  than the upstream repository's CC BY 4.0 statement, so the stricter boundary
  is recorded;
- the pinned Ensembl SpliceAI MANE v1.4 release contains masked SNV scores only.
  The adapter supplies a labelled zero-record indel file solely to satisfy the
  VEP plugin interface and writes `SpliceAI indel predictions` to
  `omissions.json`. A null indel score is unavailable evidence, never a zero-
  effect prediction.

The manifest is accepted only when its exact plugin and artifact sets, paths,
digests, versions, assembly, coverage, structural checks, and terms agree with
the configured release boundary. Preflight additionally loads every requested
Perl plugin inside the network-disabled pinned runtime and, for LOFTEE, reads a
real base from the indexed ancestor FASTA through `Bio::DB::HTS::Faidx`.
An uninitialized ancestral lookup, attempted `samtools faidx`, or related
LOFTEE warning makes VEP fail even if VEP itself exits zero. Final bundle validation
cross-checks enabled plugin arguments, checksum provenance, coverage claims,
typed score domains, and omissions. Predictor source databases and manifests
remain outside the portable `.genome` bundle; only annotations and immutable
digest references are retained.

Resource checksums and module loading are necessary but not sufficient. The
mandatory predictor doctor sends three fixed public GRCh38 controls through the
same immutable VEP 116.0 runtime and asserts exact values: CADD PHRED
`21.0/25.4/23.4`, LOFTEE `HC`, REVEL `0.675`, AlphaMissense `0.5749`, and
SpliceAI maximum delta score `0.03`. It also requires byte-for-byte allele,
genotype, ploidy, and order identity across VEP. These positive controls prove
that the configured lookup/extraction path functions for known rows; they do
not establish genome-wide performance or clinical validity. Preflight and the
portable validator bind the control definition, result rows, resource
manifest, and report SHA-256.

## ClinGen gene-level context boundary

The optional ClinGen adapter vendors a canonical, logo-free CC0 snapshot made
from the official gene-disease validity download dated 2026-08-15 and the
official archived GRCh38 dosage-sensitivity table dated 2026-08-14. The source
byte counts, SHA-256 values, URLs, row counts, allowed validity
classifications, dosage-score semantics, and attribution are compiled into the
builder. Independent resource materializations are byte-identical.

Projection is by exact, case-sensitive HGNC symbol from typed VEP transcript
annotations. Matching rows are written only to
`extensions/clingen_gene_disease.parquet` and
`extensions/clingen_dosage_sensitivity.parquet`. They are not copied to
`clinical_evidence.parquet`, do not set variant pathogenicity or actionability,
and do not create a subject diagnosis. A missing match means unavailable
context, not evidence that a gene lacks a disease relationship or dosage
mechanism.

The validator reconciles projected counts and hashes with the immutable
resource manifest, requires every projected symbol to exist exactly in
`gene_index.parquet`, verifies every classification and special dosage score,
and rejects any ClinGen gene-context row represented as variant-level clinical
evidence. `gab query ... gene` and `gab query ... clinical` expose the rows with
the same semantic warning.

## Target catalogue

`callability.target_catalog` is optional. When enabled, `path` must name a
UTF-8 tab-separated file and `expected_sha256` must match its exact bytes. The
header is deliberately strict:

```text
target_id chrom pos ref alt gene_symbol target_source target_source_version evidence_scope
```

The actual delimiter is a tab. `target_id`, source, source version, chromosome,
position, and evidence scope are required. `ref` and `alt` must either both be
blank for a locus-only target or both contain an already normalized concrete
GRCh38 allele. `gene_symbol` may be blank. Allowed evidence scopes are `site`,
`target_panel`, and `derived_interval`. Target IDs must be globally unique.

Every allele is independently passed through `bcftools norm` against the exact
configured FASTA; a catalogue entry that would move or change is rejected, not
silently rewritten. The materializer is bounded-memory and uses a work-local
disk-backed join. For each target it applies this precedence:

```text
carried core or symbolic/non-core variant
  > overlapping no-call
  > complete reference evidence covering the full normalized REF span
  > covering callable/uncallable interval
  > unknown
```

`callability.parquet` contains the normative result. The namespaced
`extensions/target_callability_evidence.parquet` row retains the catalogue file
and row digests, source record identifiers, observed ploidy, thresholds, and a
human-readable reason. A `callable=true, reference_observed=null` result is not
a negative genotype claim.

## Polygenic scoring

`annotation.pgs.models` accepts locally acquired PGS Catalog scoring files; no
scoring data is bundled with this repository. Acquire files from the official
[PGS Catalog downloads](https://www.pgscatalog.org/downloads/) page and review
the model-specific licence plus
[EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/). Each model must declare its PGS
identifier, exact SHA-256, snapshot/version, training date, trait, licence,
redistribution decision, model type, and maximum missing absolute-weight
fraction. The current adapter deliberately supports only additive format-2.0
rows with official harmonized GRCh38 coordinates and concrete biallelic SNV
alleles. Haplotype, diplotype, dosage-specific, inferred-other-allele, failed
harmonization, indel, and structurally complex rows are retained as unsupported
weight rather than guessed. A non-palindromic reverse-complement that matches
the FASTA is retained with `strand_flipped=true`; A/T and C/G rows remain
unsupported unless a future adapter supplies independent strand
disambiguation.

PGS target rows are automatically merged with the optional user target
catalogue and passed through the same reference-confidence materializer. An
effect dosage is calculated only from an exact passing carried genotype or a
complete passing explicit-reference observation. In particular:

- no-call, unknown, interval-only, low-confidence, unsupported, and different-
  ALT observations receive no dosage and no contribution;
- an effect allele equal to REF is handled with the observed ploidy, including
  a valid haploid call, rather than assuming diploidy;
- weights are not re-scaled and missing genotypes are not imputed;
- every model reports expected/found reference/alternate/missing/ambiguous
  loci and expected/used/missing absolute weight.

Normative `.genome/1.1` requires a percentile. The builder therefore publishes
a row in `prs.parquet` only when an explicit calibration TSV is SHA-bound to
the same scoring file, identifies its population/source/date/missingness scope,
and supplies a monotone empirical score-to-percentile CDF. Otherwise the raw
partial score and the reason for non-publication remain in
`extensions/pgs_model_audit.parquet`; all locus decisions remain in
`extensions/pgs_locus_audit.parquet`. PGS results are exploratory population
estimates, not diagnoses.

Calibration files are UTF-8 TSV (optionally gzip-compressed) with this exact
shape; values must be monotone and contain at least two points:

```text
#pgs_id=PGS000001
#scoring_file_sha256=<64 lowercase hex characters>
#reference_population=<documented population/cohort>
#calibration_source=<source or derivation identifier>
#calibration_source_version=<immutable version>
#calibration_date=YYYY-MM-DD
#max_missing_weight_fraction=0.0
score	percentile
<raw score>	<0..100>
```

## GWAS Catalog associations

`annotation.gwas` pins one dated NHGRI-EBI GWAS Catalog bulk release by URL,
archive byte count, SHA-256, member name/size, and ZIP CRC. `gwas-sync` streams
the strict tabular member into a deterministically sorted Parquet resource and
records every normalization count. The full source archive and normalized
resource remain outside sample bundles.

Projection is intentionally conservative. A row is joined only when all three
identities agree: exact GRCh38 coordinate, an exact rsID alias from the carried
VCF allele, and a single concrete effect allele equal to REF or ALT. Multi-locus
rows, non-core chromosomes, unknown/complex risk alleles, invalid P values, and
missing traits are counted as omissions instead of guessed. Effect-allele
dosage is computed from the accepted subject GT with its observed ploidy.

The legacy bulk column is named `OR or BETA` but does not itself identify the
effect scale. A normative `gwas_associations.parquet` row is therefore emitted
only when the accompanying text explicitly identifies beta direction, odds
ratio, or hazard ratio. Exact matches with ambiguous or missing scale/value are
retained only in `extensions/gwas_catalog_evidence.parquet`, including the raw
effect and confidence-interval text plus the source-row number and digest. Each
`association_id` commits to the frozen snapshot, normalized source digest,
unique source-row number, and matched variant. Consequently, semantically
duplicate upstream rows remain two auditable source observations without an ID
collision or silent deduplication. Variant-level `trait_associations` is a
compact retrieval summary backed by those tables.

Use `gab query BUNDLE gwas TERM` for bounded searches by trait, rsID, variant,
study accession, or PMID. These are literature-curated population associations,
not personal diagnoses, causal findings, penetrance estimates, or clinical
recommendations. The configured terms and `redistribution_allowed` decision are
recorded in the manifest; review the source terms independently before sharing
derived bundles.

## Pharmacogenomics

`annotation.pgx` binds four independent identities: the PharmCAT release and
commit, its official positions VCF, the exact container image digest, and a
locally frozen CPIC `pair_view` snapshot. The helper
`scripts/freeze_cpic_pairs.py` creates a canonical JSON snapshot; review the
applicable source terms, record the exact CPIC **data release**, date, SHA-256,
and redistribution decision in the private configuration. The freezer resolves
that release's official Git tag and commit, hashes the OpenAPI document actually
served, requires an exact `Content-Range`, and accepts the row set only after two
identical complete retrievals. This deliberately keeps the CPIC data release,
the serving PostgREST version, and the snapshot digest as separate identities.
Neither licensed source data nor sample-derived reports belongs in Git.

Freeze CPIC before `pgx-sync`, because synchronization downloads only the exact
small official PharmCAT positions asset and then verifies both configured data
boundaries. It never chooses a live CPIC snapshot implicitly:

```bash
python scripts/freeze_cpic_pairs.py \
  --source-version EXACT_CPIC_DATA_RELEASE \
  --snapshot-date YYYY-MM-DD \
  --output /private/resources/cpic/cpic_pairs.json
# Copy the printed SHA-256 and reviewed source terms into config/config.yaml.
gab resources pgx-sync --config config/config.yaml
```

The builder first merges all PharmCAT alleles into the same target catalogue
used for reference-confidence materialization. For each multiallelic source
record it reconstructs the original genotype only when all carried alleles are
exact targets with passing DP/GQ/FILTER evidence. A reference genotype is
created only from complete, passing, explicit homozygous-reference evidence.
The following remain no-calls:

- missing or incomplete genotypes;
- interval callability without a locus-specific reference observation;
- a different overlapping ALT;
- low-confidence variant or reference evidence;
- unknown evidence;
- positions in `CYP2D6`, `HLA-A`, `HLA-B`, or `MT-RNR1` that require a
  validated CNV/SV/HLA-aware outside call.

The materialized VCF uses a non-identifying internal sample label and is passed
to the exact PharmCAT image with networking disabled, a read-only container
root, no added privileges, and no missing/absent/unspecified-to-reference or
research-mode flags. Raw reports stay in the private restart-safe work area.
The portable bundle contains:

- `extensions/pgx_locus_audit.parquet` — evidence and decision for every
  official position;
- `extensions/pgx_gene_call_audit.parquet` — every PharmCAT gene call and its
  publication/withholding reason;
- `extensions/pgx_guidance_audit.parquet` — structured drug guidance linked to
  the canonical report;
- optional `pharmacogenomics.parquet` — only one complete unambiguous
  non-complex diplotype/phenotype with a non-provisional exact CPIC A-D mapping.

The strongest CPIC level is used for the compact core row; mappings at other
levels remain visible in the audit. These results are research annotations, not
prescribing instructions, and require professional review and appropriate
clinical confirmation.

## Family evidence

`annotation.family_evidence` is optional and generic. When enabled, it requires
an indexed, coordinate-sorted, normalized family joint VCF, a standard six-
column PED, explicit source version/date, and exact SHA-256 values for both
artifacts. Caller-specific FORMAT and difficult-region INFO tags are mapped in
private configuration; no deployment-specific tag names or identifiers are
compiled into the adapter.

The selected subject and every configured parent must occur literally in both
the PED and joint VCF. Complete carried family alleles are reconciled to the
accepted core through the immutable normalization-provenance ledger. This
supports a complex source allele that atomizes into multiple concrete
biallelic rows while requiring exact source identity, target cardinality,
projected genotype, ploidy, phased allele order, and phase state. Incomplete
subject genotypes such as `1/.` remain explicit no-calls: they are counted,
including their partial ALT evidence, but never become a portable variant or a
reference claim. Family evidence is monotonic; unexplained disagreement aborts
the build instead of refining or replacing the accepted genotype.

The resulting `extensions/family_evidence.parquet` contains only subject-local
or minimally derived values:

- the subject genotype, phase state, optional phase set, and source confidence;
- observed paternal, maternal, both-parent, neither-parent, or unresolved
  origin and conservative autosomal Mendelian status;
- independent-caller and technical-replicate concordance categories;
- difficult-region flags and configured de-novo, mosaicism, or dropout review
  indicators;
- immutable source digests and an evidence digest that excludes relative
  genotypes.

Parent/child identifiers and their genotypes never enter the extension. A
source tag containing any PED identifier is rejected. The manifest records
zero relative identifiers, zero relative genotypes, and zero genotype
refinements, and the validator independently joins every family row back to
`variants.parquet`.

An optional `loco` block accepts a checksum-pinned, uniquely keyed exception
TSV plus its child-independent summary JSON. Source rows may arrive in any
order; the adapter records ordering inversions and deterministically sorts them
before joining. Exact exceptions and their source child genotype are projected
through the same normalization provenance. The current contract retains source
`de_novo_candidate`, `mendelian_conflict`, and
`uninformative_parent_not_callable` states. If that source genotype differs
from the accepted core genotype, the effective state is
`source_genotype_discordant`: the source status remains available for audit,
but cannot promote a de-novo/conflict claim or alter the accepted genotype.
Because the input is not a complete per-site ledger, an absent row is recorded
as `not_individually_persisted`, not “consistent.” Aggregate consistency
metrics remain provenance only. This prevents an exception-file lookup from
becoming a false negative answer.

## Safety and interpretation

This software is research infrastructure, not a medical device. It does not
diagnose disease and does not replace clinical confirmation or professional
interpretation. Generated bundles contain sensitive genetic information and
must be protected accordingly. See [SECURITY.md](SECURITY.md) and
[docs/evidence-model.md](docs/evidence-model.md).

## Development

```bash
python -m pip install -e '.[dev]'
pytest
ruff check .
python scripts/privacy_scan.py --tracked
```

The public implementation contract is documented in
[docs/implementation-spec.md](docs/implementation-spec.md). The upstream
format remains authoritative where the two overlap.

## License

Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
