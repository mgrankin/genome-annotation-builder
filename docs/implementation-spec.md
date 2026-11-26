# Public implementation specification

This document defines the generic pipeline contract. It deliberately contains
no deployment-specific sample names, paths, pedigree, or results.

## Input contract

The primary input is one coordinate-sorted, indexed, single-sample GRCh38
reference-confidence VCF/gVCF. A multisample file is accepted only when the
sample is selected explicitly. The configured FASTA, its index, VCF contig
dictionary, and REF alleles must agree.

The input is immutable. Its SHA-256 and structural summary are recorded before
any transformation. A build stops on truncated input, invalid genotype allele
indices, unsupported ploidy, reference mismatch, duplicate normalized allele,
or incomplete output reconciliation.

## Evidence states

The ingestion layer records locus-level genotype evidence and the concrete alleles observed in those genotypes:

- `called_variant`: a concrete A/C/G/T ALT allele carried by the selected
  genotype;
- `explicit_reference_site`: a complete `0` genotype at a concrete variant
  site;
- `reference_block`: a homozygous/haploid-reference symbolic gVCF block with
  an inclusive `END`;
- `no_call`: a record or interval with an incomplete genotype;
- `symbolic_nonreference`: a symbolic non-reference call that cannot be stored
  as a small variant;
- `unsupported`: malformed, non-small-variant, or unsupported-ploidy evidence,
  which is fatal unless an explicit omission policy exists.

An absent record has no evidence state. A callable interval is not itself proof
of a homozygous-reference genotype. A reference answer at a requested locus is
permitted only when a complete explicit reference observation covers the locus
and no overlapping non-reference/no-call evidence supersedes it.

## Allele decomposition and normalization

Concrete ALT alleles present in a complete sample genotype enter
`variants.parquet`. An incomplete genotype retains its locus-level `no_call`
and false callability. A concrete ALT explicitly present in its GT also goes
through the normal normalization and VEP stages. The resulting annotations
are published in `extensions/observed_partial_variants.parquet`. For example,
`./2` projects only ALT 2, with GT `[null, 1]`. Unknown slots remain null
and their position and phase separator are preserved. Neither `0/.` nor
`./.` contributes an observed ALT, including when native dosage metadata
is positive.

The extension carries the normalized variant ID, coordinates, REF, ALT,
nullable GT slots, phase and typed gene and transcript consequences. Its
original source GT, source ID, digest, ordinal and selected ALT index link
to the no-call evidence. An allele that decomposes into multiple atoms
produces one row per atom; the source-allele accounting counts it once.
These records do not supply a complete dosage to clinical or trait models.

For each selected ALT in a complete genotype, the original genotype is projected
to a biallelic genotype while preserving allele order and phase: the selected
allele maps to `1`; reference and other ALT alleles map to `0`. Thus a `1/2`
call produces two traceable heterozygous rows, while an uncarried ALT produces
no variant row.

Each projected row receives a digest over the complete source record and
selected sample fields. `bcftools norm` left-aligns and minimizes it against the
configured FASTA. The digest survives normalization and is used for exact
source/output reconciliation. Symbolic alleles never enter `variants.parquet`.

Core columns use the exact physical widths and dictionary encodings declared by
`.genome/1.1`. If an observed DP exceeds signed `int16` or GQ exceeds signed
`int8`, the nullable core convenience field is left null instead of wrapping or
clipping; the exact integer remains in `extensions/variant_provenance.parquet`
and the build summary counts the projection.

## Bundle layers

The core `.genome/1.1` bundle contains:

- byte-identical upstream `schema.json`;
- chromosome-partitioned `variants.parquet/`;
- `manifest.json` with immutable resource and transformation provenance;
- field glossary, inventory, raw artifact, transformation, omission, and
  relationship metadata;
- a bounded read-only query interface and Agent Exploration Profile files.

The lossless extension contains:

- `extensions/reference_confidence_regions.parquet/` with reference blocks,
  explicit reference sites, no-calls, filters, DP/GQ, ploidy, and source
  digests;
- `extensions/variant_provenance.parquet/` mapping every output allele to its
  source record and projected genotype;
- `extensions/clinvar_assertions.parquet/` retaining the complete typed
  exact-allele ClinVar assertion fields whenever ClinVar is enabled;
- optional `extensions/family_evidence.parquet/` retaining only subject-local
  genotype/phase and privacy-minimized pedigree-derived evidence;
- `extensions/reconciliation.json` proving that every accepted carried allele
  produced exactly one output row and every excluded record has a reason.

Optional annotation tables are emitted only when actually built. Their absence
is listed in `omissions.json`; empty or missing sources never satisfy positive
clinical, rarity, prediction, pharmacogenomic, or trait rules.

`callable_regions.parquet` is produced with a bounded-memory coordinate sweep.
Duplicate and overlapping source records become deterministic, disjoint
segments. If evidence layers disagree, non-callable wins conservatively; the
`minimum_depth` column records the configured threshold, not an observed depth.
Observed DP/MIN_DP remains in the source-linked reference-confidence ledger.

### Target-level callability

When a frozen target catalogue is configured, every target is checksum-bound,
required to use a unique identifier, checked against the exact GRCh38 FASTA,
and independently proven to be left-aligned/minimal by `bcftools norm`. The
catalogue may contain allele-specific or locus-only rows and identifies its
source, source version, and evidence scope. A non-normalized or duplicate row
is fatal.

Materialization joins the catalogue to carried core alleles, symbolic/non-core
calls, the complete reference/no-call ledger, and callable intervals. It uses
the strict precedence `carried variant > no-call > complete reference evidence
> covering interval > unknown`. Reference evidence must cover the complete
normalized target REF span. Multiple overlapping reference records are
combined conservatively: callability requires all selected records to pass,
and retained DP/GQ/confidence use the least favourable complete evidence.

The normative result is `callability.parquet`. The one-to-one namespaced
`extensions/target_callability_evidence.parquet` row preserves catalogue file
and row digests, all selected source identifiers, applicable observed ploidy,
thresholds, precedence, resolution, and reason. Validators prove the core and
extension rows are identical on shared fields and that every source identifier
resolves. An interval-only row may say callable, but its
`reference_observed` remains null.

## PGS Catalog scoring boundary

Configured PGS resources are local, licence-explicit, and SHA-256 pinned. The
adapter requires PGS Catalog scoring-file format 2.0 and official harmonized
GRCh38 coordinates. It independently queries the configured FASTA to decide
which concrete allele is REF; neither the effect-allele label nor a harmonized
coordinate is treated as a genotype. The current implementation supports
additive, concrete biallelic SNV rows. Unsupported model rows remain in the
weight denominator and audit instead of being silently removed.

Allele orientation is explicit. A direct REF match records
`strand_flipped=false`; if neither source allele matches but a non-palindromic
reverse complement does, both alleles are complemented and
`strand_flipped=true` is retained with the source alleles. Palindromic A/T and
C/G rows are strand-ambiguous without independent frequency/source evidence
and remain unsupported rather than being guessed.

Every supported scoring row becomes a namespaced target and is merged with the
optional user target catalogue. The ordinary target-callability precedence is
therefore reused before scoring. Exact carried alleles and complete passing
reference observations may yield a dosage. A different overlapping ALT,
explicit no-call, low-confidence reference, interval-only callability, missing
evidence, unsupported ploidy, or structurally represented call yields no
dosage. Effect-allele dosage is ploidy-aware; missing genotypes are not imputed
and scores are not re-scaled.

`extensions/pgs_locus_audit.parquet` records the scoring-row digest, normalized
target, target-catalogue digest, callability proof, effect orientation, dosage,
contribution, and reason. `extensions/pgs_model_audit.parquet` records expected,
supported, alternate-called, reference-called, missing, ambiguous,
low-confidence, and unsupported locus counts; expected/used/missing absolute
weight; partial raw score; licence; and publication decision. Locus rows are
written in bounded batches and model accounting is aggregated from the Parquet
audit with a disk-backed engine. The model audit stores the exact number of
distinct contributing variants. Its convenience ID list is complete only up
to 10,000 IDs; above that bound it is null and the locus audit remains the
authoritative unbounded mapping. The nullable normative
`contributing_variants` field follows the same rule.

The normative `prs.parquet` table is stricter than raw-score computation. A row
is emitted only when model missing weight does not exceed its declared limit
and a checksum-pinned monotone empirical calibration is bound to the identical
scoring-file digest, applicable missingness limit, reference population,
source, version, and date. Missing or inapplicable calibration leaves the raw
score audit-only because `.genome/1.1` requires a non-null percentile.

## GWAS Catalog association boundary

The optional GWAS adapter accepts exactly one configured dated NHGRI-EBI bulk
release. Synchronization validates the archive URL, byte count, SHA-256, sole
member name and size, ZIP CRC, and strict header before producing a sorted
typed Parquet resource plus a manifest. The raw and normalized resources are
referenced by digest and are not embedded wholesale in a sample bundle.

A carried allele matches a source row only when the normalized GRCh38
coordinate, an rsID alias, and the concrete effect allele all agree. Effect
allele dosage is recomputed from the unchanged subject genotype and ploidy.
Unsupported multi-locus/non-core rows and ambiguous risk alleles are counted,
not coerced. The stream join is coordinate sorted and bounded-memory.

`gwas_associations.parquet` contains only rows with an explicit effect type and
numeric value. `extensions/gwas_catalog_evidence.parquet` retains every exact
match, including ambiguous legacy `OR or BETA` rows, raw effect/CI text, source
row number/digest, and the reason a normative row was withheld. The stable
association identifier commits to snapshot, normalized source digest, unique
source-row number, and matched variant. Upstream rows with identical semantic
content therefore remain independently addressable instead of colliding or
being silently deduplicated. The compact variant summary is
validator-reconciled to these two tables. All GWAS output is
population-association evidence and is never marked diagnostic or clinical.

## PharmCAT and CPIC boundary

The small-variant PGx adapter is bound to an exact PharmCAT release commit,
official positions-file SHA-256, immutable OCI image digest, and a canonical
locally frozen CPIC `pair_view` JSON snapshot. The CPIC snapshot records the
official data-release tag and commit separately from the actually served
PostgREST/OpenAPI identity. It also records the date, exact fields, exact
`Content-Range`, canonical row-set SHA-256, and two-pass stability gate; licence
and redistribution decisions remain private configuration. Resource verification
parses every row and accepts normative
recommendation mappings only when `usedforrecommendation=Yes`,
`provisional=false`, and the CPIC level is exactly A, B, C, or D. Composite,
provisional, retired, and non-recommendation rows remain provenance but do not
open the publication gate.

Every PharmCAT source ALT is independently normalized with the configured
`bcftools` and exact FASTA before it enters the target catalogue. The original
PharmCAT locus/allele and its normalized target variant ID are both retained;
this permits the matcher VCF to preserve PharmCAT's expected representation
while evidence lookup uses the same normalized GRCh38 identity as the rest of
the bundle. A changed normalization is explicit in the locus audit.

The target-evidence join reconstructs each original multiallelic genotype. A
matcher-eligible variant call requires exact passing carried alleles and
consistent genotype slots, ploidy, phase, DP, and GQ. A matcher-eligible
reference call requires complete passing explicit-reference evidence for every
ALT. Other-ALT, no-call, low-confidence, interval-only, unknown, or internally
inconsistent evidence becomes an explicit missing genotype, never `0/0`.

`CYP2D6`, `HLA-A`, `HLA-B`, and `MT-RNR1` are structurally complex and are
excluded from VCF-only matching. If such a gene occurs in the frozen positions
boundary, a gene-audit placeholder records that a separately validated
CNV/SV/HLA-aware outside call is required; absence from the matcher report is
not a negative result.

PharmCAT runs from the exact local image with `--pull never`, network disabled,
a read-only root, dropped capabilities, no new privileges, host UID/GID, bounded
processes/memory, and read-only input. Missing/absent/unspecified-to-reference
and research-mode options are prohibited. Successful process exit is
insufficient: the runner requires exactly one calls-only report and JSON report,
checks the embedded version, canonicalizes nondeterministic timestamp data, and
records raw/canonical digests.

The normative `pharmacogenomics.parquet` projection requires one unambiguous
complete non-complex source diplotype, a phenotype, and at least one qualifying
CPIC mapping. The compact row uses the strongest exact CPIC level and only drugs
at that level. All mappings, withheld calls, and structured guidance remain in
`extensions/pgx_gene_call_audit.parquet` and
`extensions/pgx_guidance_audit.parquet`; every matcher-locus decision remains in
`extensions/pgx_locus_audit.parquet`.

## Pedigree-aware family evidence boundary

The optional family adapter accepts a SHA-256-bound normalized multisample
joint VCF and six-column PED. The configured subject and each known parent must
exist literally in both. Caller-specific confidence, reason, phase-set,
independent-caller, technical-replicate, de-novo, and difficult-region tags are
mapped explicitly by configuration and verified against the VCF header.

The adapter streams the family VCF and chromosome-partitioned core ledger in
coordinate order. Every complete carried subject allele must reconcile through
the immutable normalization-provenance ledger. An unchanged source allele maps
to its one accepted core row; a complex source allele may map to multiple
concrete biallelic rows only when provenance proves exact source identity,
target cardinality, and projected GT/phase for every row. A family-only allele,
core-only allele, duplicate, undocumented representation mismatch, or
GT/phase disagreement is fatal. Incomplete subject genotypes are explicit
no-calls and are omitted from carried-allele materialization while total and
partial-ALT counts remain in provenance. They never imply either a variant or
reference. Carried symbolic/non-ACGT records remain outside this small-variant
extension and are counted explicitly.

Diploid phase is compared strictly from the `/` or `|` separator. A haploid GT
has no separator and therefore cannot encode a source phase relationship;
htslib's exposed haploid phase bit is retained only as the normalized/core
representation and is never interpreted as evidence of biological phasing.

Autosomal Mendelian compatibility is derived only when the subject and both
parents have complete diploid biallelic genotypes. Other ploidies,
sex/mitochondrial chromosomes, a missing parent, and partial calls remain
`not_assessed` or unresolved. Parental origin describes observed ALT presence,
not biological proof beyond the recorded pedigree. Caller-provided conflict or
de-novo classes may lower confidence or trigger review, but cannot modify GT.

`extensions/family_evidence.parquet` stores the subject GT/phase, optional phase
set, family confidence/reasons, inheritance and observed parental-origin state,
Mendelian status, caller/replicate concordance, difficult-region context, and
configured de-novo/mosaicism/dropout indicators. It stores neither relative
identifiers nor relative genotypes. PED identifiers found inside a configured
source version or retained family tag are fatal. Row identifiers and evidence
digests commit only the subject-local/minimally derived payload and immutable
source digests.

An optional LOCO input consists of an immutable child-independent summary and a
uniquely keyed exception ledger. Summary counters, fraction, and
child-independent-model flag are cross-checked; the entire TSV is read, its
status counts must reproduce the summary, source-order inversions are counted,
and rows are deterministically sorted before the join. Exact exceptions are
projected through the same normalization provenance and may annotate a core
row. The extension retains both source status and source-child-GT concordance.
The producer may omit zero-valued summary counters. The adapter reconciles the
carried, informative and exception totals and checks the embedded review digest
when supplied. An uninformative child projection (`child_gt=.`) retains the
`uninformative_child_genotype` status and has null genotype concordance.
If a complete source child genotype disagrees with the accepted core genotype, the
effective status is `source_genotype_discordant`; the original status remains
auditable but cannot be promoted into a de-novo or Mendelian-conflict claim.
The accepted core genotype is never changed. The absence of an exception is
deliberately represented as `not_individually_persisted`, because an
exception-only file is not a complete site-level consistency ledger.

The bundle validator rechecks the extension's exact physical schema, digest
boundary, identifier uniqueness, normalization-projection accounting, LOCO
source/effective-status state machine, privacy counters, and one-to-one
coordinate/GT/phase equality with `variants.parquet`. The bounded `evidence`
query returns this subject-only row alongside ordinary variant provenance.

## Offline VEP annotation boundary

When enabled, functional annotation uses Ensembl VEP 116.0 with indexed human
GRCh38 cache 116. The container reference is an immutable OCI digest; the
archive is checked against Ensembl's published size/BSD sum and the pipeline's
SHA-256 allow-list, and the extracted tree receives a separate SHA-256 digest.
The container runs without a network and with read-only input, reference, and
cache mounts.

Before and after VEP, an order-sensitive digest covers CHROM, POS, REF, ALT,
GT, phase, and source-provenance identifiers. Any change is fatal. CSQ is not
stored as an opaque interpretation surface: every transcript is percent-
decoded into a typed list. Transcript selection is gene-aware and
deterministic: MANE Select, MANE Plus Clinical, canonical, APPRIS principal,
and consequence severity choose one representative within each gene; the
primary gene is then chosen by the representative consequence before
transcript rank. `most_severe` is calculated independently across all
transcripts. This prevents a nearby gene's upstream/downstream MANE annotation
from overriding a coding consequence in another gene.

Cache-provided gnomAD genome frequencies take precedence over exome
frequencies for core population fields. Population max considers every
available cache population group, including groups without dedicated
`.genome/1.1` columns. Missing frequency or predictor values make no benign or
rarity claim.

The optional prediction boundary consists of exactly five explicit offline
adapters: LOFTEE at one source commit, CADD v1.7 GRCh38 SNV and small-indel
tables, REVEL v1.3, the exact hosted AlphaMissense GRCh38 release, and Ensembl's
SpliceAI 1.3.1 masked MANE v1.4 GRCh38 SNV release. Every enabled resource has a
compiled-in URL, byte count, SHA-256 and, where published, MD5. Unfrozen
identities fail closed. All downloaded files, generated indexes, extracted
source trees, and deterministic derivatives live outside both the public Git
tree and the portable sample bundle.

Resource synchronization is restartable but promotion is atomic. REVEL rows
without a GRCh38 position are counted and omitted; retained `1-22,X,Y` rows are
stably sorted without dropping transcript identifiers, BGZF-compressed, and
tabix-indexed. LOFTEE archives reject links and path traversal; its ancillary
bigWig, indexed ancestor FASTA, and conservation SQLite database receive
independent structural checks. The exact pinned upstream `LoF.pm` is additionally
transformed by the shipped `biodbhts-faidx-v1` unified diff because the pinned
VEP image has `Bio::DB::HTS::Faidx` but no external `samtools`. GNU patch 2.7.6
must apply with zero fuzz to the fixed preimage; the patch SHA-256, original and
patched `LoF.pm` SHA-256 values, transformation semantics, and resulting source
tree digest are mandatory manifest fields. The unmodified extracted tree is not
rewritten or accepted as the runtime tree. AlphaMissense receives a verified local tabix
index. CADD and SpliceAI publisher indexes are opened and queried rather than
trusted by filename.

The official Ensembl SpliceAI release in this boundary contains SNVs only. A
deterministic, labelled, zero-record indel VCF is generated solely because the
VEP plugin requires both arguments. Its existence is not coverage: resource
metadata fixes `indel_scores_available=false`, and every resulting bundle must
declare SpliceAI indel prediction as an omission. A future indel resource
requires a new official, licensed, checksum-pinned adapter and cannot be
substituted through configuration.

VEP receives predictor resources through a read-only mount and sanitized
`PREDICTORS/...` command semantics. Predictor adapters require the pinned Docker
runtime and its `/plugins` tree; arbitrary executable-mode host plugin code is
rejected because it lacks an immutable code identity. Preflight loads every
enabled Perl module inside the same pinned, network-disabled runtime. LOFTEE's
runtime check also reads `1:1-1` from the indexed ancestor resource through
`Bio::DB::HTS::Faidx`; it does not merely import the module. VEP output must declare all
required CSQ fields. CADD is finite and non-negative; REVEL, AlphaMissense,
SpliceAI, SIFT, PolyPhen, and frequency values are finite and bounded to their
defined domains. Completeness denominators are consequence-aware: LOFTEE is
assessed only for predicted loss-of-function consequences, CADD only for its
supported SNV/indel classes, REVEL and AlphaMissense only for missense SNVs,
and the pinned SpliceAI data only for MANE-overlapping SNVs. Missing applicable
scores are counted as missing evidence and cast no damaging/benign vote.
Known LOFTEE ancestral-lookup warnings are fatal even when the VEP process exits
zero, preventing a silently disabled filter from producing a PASS checkpoint.

The private predictor manifest records relative paths only. Portable
`raw_artifacts.json` repeats one common manifest digest plus every upstream and
derived artifact digest, terms, conservative non-redistribution decision, and
coverage. `transformations.json` records the matching plugin arguments and
coverage; `omissions.json` records every disabled adapter and the mandatory
SpliceAI indel gap. The final validator cross-checks these three documents and
rejects predictor values for an adapter that was not enabled.

A checksum-valid database and importable Perl module can still be wired
incorrectly. Therefore every enabled predictor configuration must pass the
fixed public semantic-control set through the same immutable VEP runtime before
preflight succeeds. The controls assert CADD PHRED `21.0`, `25.4`, and `23.4`;
LOFTEE `HC`; REVEL `0.675`; AlphaMissense `0.5749`; and SpliceAI maximum delta
score `0.03`, while requiring unchanged alleles, genotypes, ploidy, and record
order. The control definition, resource manifest, result rows, and report
SHA-256 are cross-linked through raw artifacts, transformations, and manifest
provenance. This is a positive wiring test, not a genome-wide performance or
clinical-validity claim.

## ClinGen gene-context boundary

The ClinGen adapter is independent of the ACMG Secondary Findings policy
snapshot. It freezes two official sources into one deterministic CC0 canonical
object: 3,663 gene-disease validity assertions dated 2026-08-15 and 1,526
GRCh38 gene dosage-sensitivity rows from the 2026-08-14 archive. Source bytes,
SHA-256, URLs, headers, creation dates, enums, row counts, coordinates, special
dosage scores (`30`, `40`, and nullable `Not yet evaluated`), and lossless PMID
source cells are validated before the canonical object is accepted.

`clingen-sync` deterministically materializes typed Parquet resource indexes;
two independent materializations must be byte-identical. Bundle projection
uses exact, case-sensitive HGNC symbols taken from VEP transcript annotations.
Only non-empty projections are emitted as
`extensions/clingen_gene_disease.parquet` and
`extensions/clingen_dosage_sensitivity.parquet`. An absent projection receives
an explicit omission and is never interpreted as negative evidence.

These are gene-level knowledge rows, not variant assertions. They never enter
`clinical_evidence.parquet`, never set pathogenicity, actionability, or
`clinical_grade`, and never create a diagnosis. Validation requires exact joins
to `gene_index.parquet`, reconciles row counts and output hashes to provenance,
checks classification and dosage semantics, and rejects any ClinGen row leaked
into variant-level clinical evidence.

## ClinVar annotation boundary

When enabled, the ClinVar adapter accepts only the compiled-in dated NCBI
GRCh38 snapshot. The downloaded archive, publisher MD5 sidecar, tabix index,
normalized VCF, and normalized index are checksum-verified. Normalization uses
the exact configured GRCh38 FASTA, retains primary contigs only, and fails on a
REF mismatch.

Matching is coordinate-sorted and exact on CHROM/POS/REF/ALT; it does not join
by position alone. Every assertion remains available in the lossless ClinVar
extension with raw significance, review status, conditions, submitters,
accessions, and snapshot date. A corresponding row in normative
`clinical_evidence.parquet` provides the standard source/version/assertion
retrieval surface. Their evidence identifiers must match one-to-one. Aggregate
significance, review stars, and conflict state are separate fields. These
annotations are external assertions, not diagnoses or expert review performed
by this pipeline.

## ACMG Secondary Findings boundary

The actionability knowledge boundary is ACMG SF v3.3 (PMID 40568962, DOI
`10.1016/j.gim.2025.101454`) plus a ClinGen curation workbook frozen on
2026-08-15. The workbook is accepted only at its compiled-in byte count and
SHA-256. A deterministic extractor removes presentation assets and produces a
canonical CC0 JSON snapshot containing:

- all 84 policy genes and 100 policy gene-disease rows;
- 177 ClinGen disease-level reportability rows;
- one reporting-guidance row for every policy gene;
- 104 MANE/BA1/BS1/DAFT frequency-specification rows;
- an explicit `under_construction` state for variant-specific guidance.

The sorted 84-gene set has an independent digest, and the policy, ClinGen
reporting, and guidance gene sets must match exactly. This prevents a currently
incomplete convenience CSV (81 unique genes at the frozen date) from silently
becoming the policy boundary. ClinGen curated content is CC0 with attribution.

Membership in an ACMG gene is not itself an actionable finding. When the
adapter is enabled, the versioned evaluator requires exact disease identifiers,
matching VEP and ClinVar genes, P/LP significance, sufficient review, no
conflict, high analytical confidence, compatible inheritance/genotype, and all
implemented variant-class/gene-specific rules. Unresolved phase, mechanism,
transcript, or second-allele questions receive explicit review statuses instead
of being promoted.

Each retained policy decision is written to a namespaced audit extension and
materialized as a normative finding, policy evidence row, and structured
interpretation. The interpretation status is `algorithmic`, `clinical_grade`
is false, `confirmation_required` is true, and reviewer fields are null. The
validator checks every variant, finding, evidence, interpretation, and policy
link. Only decisions that pass every implemented gate set variant-level
`is_actionable=true`; this remains a retrieval candidate rather than a clinical
result.

## Workflow execution boundary

The checked-in Snakemake workflow uses a versioned Conda environment and
checksum/resource marker files outside the final bundle. Snakemake never owns
or pre-creates the `.genome` directory: only the builder's atomic promotion may
create that immutable artifact. The top-level outputless rule deliberately runs
on every invocation. A first run builds and validates the bundle; a repeated run
replays verified checkpoints or validates the existing final artifact. Source,
workflow, lockfile, configuration, VCF, reference, and resource markers are all
explicit inputs.

## Validation gates

A final bundle is promoted only after all of the following pass:

1. normative file revision and checksum checks;
2. manifest schema and file digest checks;
3. exact Parquet physical widths/dictionary encodings and row validation
   against the relevant `.genome/1.1` definitions;
4. A/C/G/T allele and GRCh38 REF checks;
5. genotype/ploidy/zygosity/allele-balance invariants;
6. one-to-one source/output reconciliation;
7. reference-confidence interval ordering and state consistency, including
   unique non-overlapping callable regions;
8. positive, other-ALT/symbolic overlap, explicit-reference, callable-only,
   no-call, low-confidence reference, haploid, and unknown query truth-table
   tests;
9. Agent Exploration required-path and bounded-query checks;
10. clinical finding/evidence/interpretation and lossless-extension referential
    integrity;
11. target catalogue normalization, digest, precedence, one-to-one projection,
    and source-link validation;
12. PGS model/locus accounting, target-proof joins, ploidy-aware dosage,
    missing-weight, calibration, and normative-publication validation;
13. complete PharmCAT position/allele accounting, normalized target mapping,
    per-locus evidence recomputation, forbidden-option/runtime gates, complex-
    gene withholding, CPIC publication semantics, and core/audit linkage;
14. manifest implementation-tree, CPIC release/OpenAPI/row-set, resource,
    container, and licence provenance boundaries;
15. family source/PED/LOCO accounting, subject-only schema, one-to-one
    coordinate/GT/phase reconciliation, and zero-relative-data invariants;
16. GWAS source identity, exact coordinate/rsID/effect-allele joins, genotype
    dosage, typed/lossless reconciliation, compact-summary consistency, and
    population-only interpretation boundaries;
17. predictor artifact-set/digest/licence boundaries, VEP plugin command
    semantics, applicable score domains, coverage claims, and explicit disabled
    or unsupported-class omissions;
18. per-bundle and project-wide annotation QC reconciliation, including
    variant/genotype distributions, Ti/Tv, VEP/ClinVar/gnomAD/predictor
    completeness, difficult-region, PGx/PGS, callability, family,
    technical-replicate, and read-backed summaries, with semantic replay of
    JSON/TSV/HTML identities;
19. privacy scan of the public Git index.

## Release sequence

- Release A: lossless small-variant and reference-confidence core, pinned
  offline VEP functional/population annotation, exact-allele ClinVar evidence,
  conformance, provenance, callability, validation, and query profile.
- Release B: additional frozen predictor and gene-disease resources, reviewed
  clinical rules/findings, the conservative small-variant PharmCAT/CPIC layer,
  PGS/GWAS, family evidence, and transparent rule definitions.
- Release C: read-backed SV/CNV, repeat and complex-locus calling, mitochondrial
  analysis, and other assay-specific extensions.

The releases are cumulative. Later interpretation never overwrites the raw
genotype or hides an earlier evidence class.
