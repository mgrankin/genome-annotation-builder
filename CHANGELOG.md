# Changelog

## Unreleased

- Reconcile CNVpytor's four-decimal TSV fractions against full-precision HDF5
  values with an exactly inclusive half-unit decimal boundary. This avoids a
  binary floating-point false mismatch at `0.00005` while values outside the
  caller's documented rendering precision continue to fail closed.
- Preserve standard IUPAC repeat motifs emitted from a digest-pinned
  ExpansionHunter variant catalog. Ambiguous motifs remain typed rows, are
  forced to `REVIEW` with `CATALOG:IUPAC_AMBIGUOUS_MOTIF`, and are counted
  separately; non-IUPAC symbols still fail closed.
- Correct the pinned mitochondrial source-contig policy to match the exact
  shipped GRCh38 reference (`MT`). The workflow continues to emit canonical
  portable mitochondrial coordinates as `chrM`; ship the byte-identical policy
  example and verify its digest against the immutable caller profile.
- Decode GATK `AS_FilterStatus` without losing allele boundaries when GATK uses
  pipes between alleles and commas between multiple filters for one allele.
  Preserve standards-compliant Number=A input and fail closed on every
  genuinely irreconcilable cardinality.
- Run both SpecHLA stages through the pinned image's `env-execute` bootstrap.
  The OCI runner deliberately replaces image entrypoints, while the SpecHLA
  wrappers require the Conda environment initialized by that bootstrap;
  command provenance now records the bridge explicitly and remains path-free.
  Bind the repaired SpecHLA resource tree containing all twelve Bowtie2 mapping
  index members required by the two baseline mapping FASTA files.
- Correct the HLA*LA resource boundary so `--customGraphDir` mounts the exact
  parent directory containing the named `PRG_MHC_GRCh38_withIMGT` graph.
  HLA*LA appends the graph name itself; mounting the graph contents directly
  makes the caller reject an otherwise valid resource tree. Re-freeze the
  exact parent-tree manifest and caller-profile identity.
- Correct the CNVpytor 1.3.2 combined HDF5 audit to require the exact signals
  consumed by `call_2d`: corrected/masked parity RD plus unflagged binned parity
  and unique RD. CNVpytor does not emit a GC-corrected unique-RD histogram.
  Re-freeze the caller-profile identity so the mounted audit implementation is
  still checked byte-for-byte before execution. Re-audit the 10-kb dataset
  after the 100-kb call mutates the shared HDF5, so both published audits bind
  the same final `.pytor` SHA-256.

## 0.2.0 - 2026-08-17

- Add generic deterministic multi-bundle release reporting with full enriched
  bundle replay, project `SHA256SUMS`, portable manifests, JSON/TSV/HTML QC,
  a deterministic path-free JSONL validation log, per-bundle machine-readable
  validation reports, and explicit technical-replicate concordance metrics.
- Add bounded-memory per-bundle and project-wide annotation QC for record
  reconciliation, variant/genotype classes, Ti/Tv, indel lengths, VEP,
  ClinVar, gnomAD, applicable predictors, difficult regions, PGx/PGS,
  callability, family evidence, technical replicates, and read-backed domains;
  bind and semantically replay its JSON, long-form TSV, and HTML reports, and
  reject unexplained cross-bundle QC outliers under fixed public median rules.
- Add a generic, receipt-bound Release C for read-backed SV, CNV, repeat,
  HLA, complex-PGx, mitochondrial, and bounded local-read evidence.
- Ship nine immutable caller profiles with exact OCI, executable, GRCh38,
  resource-file, and resource-tree identities, plus fail-closed catalog and
  resource-manifest validators. Matrix execution rejects legacy or hand-written
  freezes without a canonical manifest identity and explicit privacy boundary.
- Add sequential matrix execution, immediate pre-call input rehashing,
  restart-safe caller lifecycles, scope-reconciled typed adapters, and retained
  failed attempts without caller voting or missing-to-negative conversion.
- Validate every Manta BND pair before applying the portable core-contig scope,
  omit both endpoints together when either is non-core, and publish exact
  aggregate raw/typed/omitted reconciliation without exposing loci.
- Revalidate and reuse a completed immutable caller lifecycle when only its
  strict adapter failed, instead of rerunning the expensive caller on retry.
- Represent DELLY's documented `RDCN=-1` missing/invalid sentinel as a null
  copy number with an explicit aggregate count and row limitation; reject every
  other negative copy-number value.
- Bind every caller receipt, including zero-row runs, to an immutable strict-
  adapter summary in the final extension manifest; publish exact aggregate
  raw/typed/scope-omission arithmetic in bundle omissions and annotation QC
  without exposing loci, alleles, subjects, or host paths.
- Add technical-replicate concordance and per-subject extension staging while
  preserving every portable typed source row and reporting single-source
  domains explicitly.
- Add atomic `.genome/1.1` extension, seven Arrow schemas, complete provenance,
  validation, bounded `read-backed` queries, synthetic controls, and public
  examples containing no personal data or environment-specific paths.

## 0.1.18 - 2026-08-16

- Treat the phase bit exposed by htslib for a separatorless haploid genotype as
  representational during family normalization reconciliation. Haploid allele,
  genotype, ploidy, source identity, and target-cardinality checks remain
  fail-closed; diploid `/` versus `|` semantics remain strict.
- Add end-to-end regressions proving that a complex haploid source allele can
  project to multiple normalized rows without a false phase conflict, while an
  actual diploid phase change is still rejected.

## 0.1.17 - 2026-08-16

- Include the immutable predictor `control_definition_sha256` directly in the
  portable semantic-control transformation row, matching the manifest and raw
  artifact boundaries checked by the fail-closed bundle validator.
- Make GWAS `association_id` commit to both the normalized source-record digest
  and its unique source row number. This preserves every row from a frozen
  snapshot even when upstream contains semantically duplicate rows, without
  ID collisions or silent deduplication.
- Add builder-level predictor provenance and repeated-GWAS-source-row
  regressions, including deterministic replay and normative/lossless ID
  reconciliation.

## 0.1.16 - 2026-08-16

- Treat incomplete subject genotypes in the family joint VCF as explicit
  no-calls. Count both all incomplete records and the subset carrying partial
  ALT evidence, but never convert `1/.` or `./1` into a portable variant or a
  reference claim.
- Project family and leave-one-child-out annotations through the immutable
  variant-provenance ledger when one complex source allele normalizes to
  multiple concrete biallelic rows. Validate source identity, projected
  genotype and phase, target cardinality, and digest accounting end to end.
- Accept uniquely keyed LOCO exception ledgers in arbitrary source order,
  record the number of ordering inversions, and sort the immutable rows
  deterministically before the streaming join.
- Preserve the original LOCO status and child genotype separately. When the
  LOCO child genotype differs from the accepted core genotype, publish the
  effective state `source_genotype_discordant`; never promote that row to a
  de-novo or Mendelian-conflict claim and never modify the accepted genotype.
- Add end-to-end complex-allele, partial-genotype, out-of-order LOCO, and
  genotype-discordance reconciliation regressions.

## 0.1.15 - 2026-08-16

- Read only the required `genotype.gt` and `genotype.phased` Parquet leaves
  while materializing family evidence. This avoids PyArrow's unsupported
  conversion of a dictionary-backed nested struct spanning multiple row
  groups, without changing any variant or genotype semantics.
- Add a multi-row-group regression matching the production failure boundary.

## 0.1.14 - 2026-08-16

- Rebuild the normalized VCF header in the public canonical core-contig order
  `chr1` through `chr22`, `chrX`, `chrY`, `chrM` before the final external sort.
  Preserve non-core reference contigs after the core set in their original
  relative order.
- Add a full-stream fail-closed normalization gate for header order, coordinate
  monotonicity, concrete biallelic representation, and duplicate query alleles.
  This prevents a valid but lexically ordered reference dictionary such as
  `1,10,...,19,2` from violating downstream streaming ClinVar/GWAS merge joins.
- Add a regression using a lexically ordered synthetic GRCh38 header and a
  streaming ClinVar lookup, proving the normalized output is naturally ordered
  without changing the carried genotype or normalized allele set.

## 0.1.13 - 2026-08-16

- Preserve the selected-ALT projected genotype explicitly before normalization.
- Convert bcftools atomization's temporary `CONCRETE,*` representation into
  one concrete biallelic row while restoring and validating GT/phase against
  source allele provenance; retain the raw atomized BCF for audit.
- Add a phased complex-allele regression proving one-to-many normalization,
  source/output reconciliation, and absence of spanning-star alleles from the
  portable variant ledger.

## 0.1.12 - 2026-08-15

- Add a deterministic CC0 ClinGen gene-level adapter covering all 3,663 dated
  gene-disease validity assertions and 1,526 archived GRCh38 dosage-sensitivity
  rows, with exact source SHA-256, typed schemas, lossless special-score/PMID
  semantics, offline sync/verify commands, and byte-identical replay tests.
- Project only exact HGNC-symbol matches into dedicated gene-context
  extensions; keep them out of variant-level clinical evidence and emit zero
  pathogenicity, actionability, clinical-grade, or diagnosis claims.
- Add provenance, inventory, relationship, omission, bounded-query, and strict
  validator gates for ClinGen, including rejection of orphan symbols, semantic
  boundary drift, count/hash mismatches, or gene-context leakage into clinical
  evidence.
- Add a mandatory five-adapter predictor doctor using fixed public GRCh38
  controls with exact CADD, LOFTEE, REVEL, AlphaMissense, and SpliceAI expected
  values; bind its definition/results/report to preflight and portable
  provenance while preserving VEP allele/genotype identity.
- Make the predictor resource manifest byte-reproducible by removing local
  installation time from schema v3; identical upstream/derived bytes now retain
  one stable manifest identity and semantic-control run boundary.
- Pin and package a zero-fuzz LOFTEE `biodbhts-faidx-v1` source patch with
  exact upstream, patch, preimage, patched-file, and source-tree SHA-256
  identities. Use the VEP image's `Bio::DB::HTS::Faidx` instead of an absent
  external `samtools` for ancestral-allele lookup.
- Probe the real indexed ancestor FASTA inside the immutable, network-disabled
  VEP runtime and fail closed on LOFTEE ancestral-lookup warnings even when VEP
  exits zero.
- Track every packaged implementation file recursively in the Snakemake DAG so
  a vendor patch or normative asset change invalidates the build checkpoint.
- Add a fail-closed wheel/sdist release gate and source-manifest exclusions so
  private/runtime directories, genomic payloads, Python bytecode, and local
  cache artifacts cannot enter a public distribution archive.
- Require a PEP 639-capable Setuptools build backend so the SPDX licence
  expression and licence-file metadata are reproducible in isolated builds.

## 0.1.11

- Add a strict dated NHGRI-EBI GWAS Catalog resource boundary with pinned URL,
  archive size/SHA-256, sole ZIP member identity/CRC, exact TSV header,
  deterministic normalization, typed Parquet schema, and restart-safe sync and
  verification commands.
- Match carried alleles only by exact GRCh38 coordinate, rsID alias, and effect
  allele equal to REF or ALT; calculate effect-allele dosage from the unchanged
  subject genotype and preserve every exact source-row match in a lossless audit.
- Emit normative `gwas_associations.parquet` only when effect type and value are
  explicit. Keep ambiguous legacy `OR or BETA` rows audit-only, preserve raw
  effect/CI text, and orient explicitly labelled beta decreases without guessing.
- Populate the compact variant `trait_associations` structure from a deterministic
  typed association, while retaining all traits/PMIDs and reconciling it against
  both typed and lossless tables during validation.
- Add bounded GWAS queries by trait, rsID, normalized variant, study accession,
  and PMID, including typed associations reached through their lossless source
  row; state the population-only, non-diagnostic interpretation boundary.
- Record source terms, redistribution decision, raw-artifact disposition,
  transformations, normalization exclusions, untyped exact matches, and exact
  bundle/resource counts in manifests and omissions.
- Add Snakemake GWAS resource synchronization and end-to-end resource,
  projection, bundle, validation, query, dosage, and false-`or` regression gates.

## 0.1.10

- Add an optional SHA-256-bound family joint VCF/PED adapter with configurable
  FORMAT/INFO tag mapping and an optional child-independent LOCO exception
  ledger.
- Reconcile every portable family-evidence row one-to-one with the accepted
  subject variant, genotype, ploidy, and phase; any attempted family-driven
  genotype refinement is fatal.
- Emit `extensions/family_evidence.parquet` with subject-only inheritance,
  parental-origin, phase-set, Mendelian, caller, replicate, difficult-region,
  de-novo, mosaicism/dropout, and exact LOCO exception statuses.
- Keep all relative identifiers and genotypes out of the extension and reject
  configured family tag values that contain a PED identifier.
- Treat absence from an exception-only LOCO ledger as
  `not_individually_persisted`, never as per-site Mendelian consistency.
- Validate physical schema, immutable family source boundaries, privacy
  invariants, one-to-one variant linkage, and exact subject GT/phase equality;
  expose the extension through the bounded evidence query.
- Add end-to-end trio, LOCO, parental-origin, de-novo, difficult-region,
  technical-replicate, privacy-leak, and genotype-mismatch regression gates.

## 0.1.9

- Add an immutable PharmCAT 3.4.0 and frozen CPIC pair-snapshot boundary with
  exact positions, release, container, source, licence, date, and SHA-256
  verification.
- Bind a CPIC freeze to the official data-release tag and commit, the served
  OpenAPI/PostgREST identity, exact pagination, and two identical complete
  retrieval passes; keep those identities separate from the snapshot SHA-256.
- Materialize every official PharmCAT locus from target-level evidence and
  preserve a complete typed audit. Only exact passing carried genotypes or
  complete explicit-reference records reach the matcher; no-call, other-ALT,
  low-confidence, interval-only, and unknown evidence never becomes reference.
- Run the digest-pinned PharmCAT container offline with a read-only root,
  dropped capabilities, host UID/GID, bounded resources, exact output gates,
  and none of PharmCAT's missing-to-reference or research-mode options.
- Withhold structurally complex CYP2D6, HLA-A, HLA-B, and MT-RNR1 VCF-only
  matches pending separately validated CNV/SV/HLA-aware outside calls.
- Project only complete, unambiguous, non-complex gene calls with an exact
  non-provisional CPIC A-D recommendation mapping into normative
  `pharmacogenomics.parquet`; retain all gene decisions and structured guidance
  in audit extensions.
- Recompute every PGx locus from bundled callability/variant evidence during
  validation and cross-check normative calls, guidance, resource provenance,
  counts, and digests.
- Add a lint-clean Conda-backed Snakemake driver whose outputless top-level rule
  leaves atomic bundle ownership to the builder and revalidates immutable
  completed bundles on every workflow invocation.

## 0.1.8

- Add strict, licence-explicit, checksum-pinned PGS Catalog format-2.0 resource
  verification and automatic union of normalized PGS loci with the target
  callability catalogue.
- Add ploidy-aware additive scoring from exact carried or explicit-reference
  genotype evidence; no-call, other-ALT, low-confidence, interval-only,
  unknown, and unsupported loci never receive an invented dosage.
- Preserve source and oriented alleles, explicitly audit safe non-palindromic
  strand flips, and reject strand-ambiguous palindromic SNPs.
- Add per-locus and per-model PGS audit extensions with complete locus and
  absolute-weight accounting, source/calibration digests, and reasons.
- Aggregate genome-wide PGS audits from disk instead of retaining all loci in
  memory. Embed contributing variant IDs only up to a documented bound; retain
  the authoritative unbounded mapping in the locus audit.
- Publish normative `prs.parquet` only for scores passing declared missing-
  weight limits and an applicable digest-bound empirical population
  calibration; retain valid uncalibrated scores audit-only.
- Add end-to-end calibrated/uncalibrated, REF-effect, ALT-effect, haploid,
  no-call, interval-only, other-ALT, low-confidence, and unknown regressions.

## 0.1.7

- Add an optional checksum-pinned, normalized target-catalogue input and
  materialize normative `callability.parquet` rows for arbitrary ClinVar,
  ACMG, PGx, PGS, GWAS, or user-defined loci.
- Preserve a row-level target-to-source proof chain in
  `extensions/target_callability_evidence.parquet`, including target/resource
  digests, carried variants, symbolic/non-core calls, no-calls, complete
  reference records, covering intervals, ploidy, DP/GQ thresholds, and reasons.
- Enforce conservative `variant > no-call > complete reference > interval >
  unknown` precedence, exact physical schemas, referential integrity, bounded
  disk-backed normalization validation, and positive/reference/no-call/unknown
  query truth tables.
- Search clinical findings by gene and identifiers, and adopt SPDX package
  metadata without setuptools license deprecation warnings.

## 0.1.6

- Match the normative `.genome/1.1` physical widths and dictionary encodings
  across core variants, callable regions, and the gene index; preserve source
  DP/GQ outside normative widths in the provenance extension.
- Make all read-only query results directly JSON serializable and normalize
  timezone-aware values to deterministic UTC strings.
- Add explicit physical-schema regression tests.

## 0.1.5

- Emit the normative `.genome/1.1` `clinical_evidence.parquet`,
  `clinical_findings.parquet`, and `interpretations.parquet` schemas for ACMG
  policy candidates while keeping every result algorithmic, non-clinical-grade,
  confirmation-required, and linked to versioned evidence.
- Preserve lossless exact-allele ClinVar records separately in the namespaced
  `extensions/clinvar_assertions.parquet` table and validate the one-to-one
  projection into normative clinical evidence.
- Make the gene index physically conformant (`dictionary<string>` and `int32`)
  and add the required `pgx_count`, with exact variant/gene accounting.
- Validate all finding, evidence, interpretation, policy, variant, and extension
  links and expose the complete chain through bounded read-only queries.
- Declare the timezone runtime needed by DuckDB when retrieving normative UTC
  clinical-evidence timestamps.

## 0.1.4

- Make gene-index actionability counts gene-specific by joining the exact
  variant/policy-gene decisions, preventing an actionable allele in one gene
  from being attributed to every overlapping transcript gene.

## 0.1.3

- Integrate the pinned ACMG SF v3.3/ClinGen snapshot into build preflight,
  restart keys, resource provenance, and bundle manifests.
- Add a disease-aware secondary-findings evaluator that requires exact-allele
  ClinVar P/LP evidence, sufficient review, no conflict, high analytical
  confidence, compatible inheritance/genotype, and implemented gene-specific
  rules; gene membership alone never qualifies.
- Preserve every automated decision in a typed extension, require clinical
  confirmation, and keep unresolved compound heterozygotes and mechanism- or
  transcript-specific cases in explicit review states.

## 0.1.2

- Search variant-level ClinVar aggregates, typed raw assertions, and optional
  reviewed findings together in the bounded clinical query.
- Distinguish external ClinVar annotation from a versioned actionability or
  pipeline-reviewed clinical-grade interpretation policy.
- Add a deterministic ACMG SF v3.3 / ClinGen resource extractor, offline
  vendored snapshot, source-workbook cross-check, CC0 attribution, and strict
  84-gene/multi-table completeness gates.

## 0.1.1

- Build `callable_regions.parquet` with a bounded-memory sweep-line so duplicate
  and overlapping source records produce deterministic, disjoint intervals.
- Interpret `minimum_depth` as the applied callability threshold and use
  conservative non-callable precedence where source intervals overlap.

## 0.1.0 — unreleased

- Add checksum-pinned `.genome/1.1` core resources.
- Add strict reference-confidence VCF classification and source accounting.
- Add genotype-aware multiallelic projection, GRCh38 normalization, atomization,
  and chromosome-partitioned Parquet output.
- Add compact reference/no-call, non-core, and variant-provenance extensions.
- Add atomic bundle assembly, explicit omissions, Agent Exploration Profile,
  bounded read-only queries, and strict validation.
- Add deterministic synthetic tests, privacy scanning, and CI.
- Add checksum-pinned, network-isolated Ensembl VEP 116.0 execution and
  restart-safe cache installation.
- Add typed per-transcript CSQ parsing, deterministic transcript selection,
  independent most-severe consequence, gnomAD frequency parsing, and a derived
  gene index.
- Add pre/post annotation allele/genotype identity gates and VEP-specific
  synthetic tests.
- Add an immutable dated NCBI ClinVar resource installer, exact-reference
  normalization, exact-allele streaming join, typed assertion evidence, and
  conflict-preserving aggregation.
