# Release C: checksum-bound read-backed extensions

Release C adds evidence that cannot be reconstructed from a small-variant VCF:
structural variants, copy-number changes, repeat expansions, HLA typing,
mitochondrial heteroplasmy, and bounded local read review. It does not modify a
small-variant genotype. Each caller observation remains a separately
provenanced row, including observations that disagree.

This layer is designed for arbitrary GRCh38 short-read WGS samples. Subject
identifiers, local paths, alignments, resource databases, raw caller output,
and runtime logs remain private. Only path-free typed evidence, exact hashes,
and audit receipts are portable.

## Security and evidence boundary

Every caller run has the following chain:

1. A private manifest binds one or more technical datasets to one logical
   subject, one exact GRCh38 FASTA, and an explicit per-caller alignment role.
2. The input gate immediately re-hashes the selected BAM/CRAM and index,
   validates reference/alignment contigs, sort order, read groups, sample
   cardinality, indexed reads, and free space.
3. A public caller profile is compiled with private mount bindings into a
   restart-safe run plan. Static tools and resource trees must match declared
   digests before the plan is published.
4. The OCI runner uses a locally present digest-pinned image with no network,
   a read-only root, a non-root user, dropped capabilities, bounded resources,
   one writable workspace, and no shell interpolation.
5. A caller receipt binds the PASS execution report to the exact input-gate
   identity, caller profile, resources, and raw artifact digest.
6. A strict adapter converts raw caller output into a typed Parquet domain.
7. The extension stage verifies every row against its input gate, receipt, and
   immutable strict-adapter summary, strips private bindings, and atomically
   publishes a portable extension set. Zero-row adapter runs remain explicit,
   so an empty result cannot be confused with a caller that never ran.

An absent variant-only caller row is never reference evidence. Caller
concordance does not use majority voting and never deletes or rewrites portable
typed source rows. The portable coordinate scope is GRCh38 `chr1`–`chr22`,
`chrX`, `chrY`, and `chrM`. Raw Manta records on alternate, decoy, or unplaced
contigs remain bound by the raw-artifact digest and are reported as explicit
scope omissions. Every BND pair is checked for reciprocal IDs, coordinates,
confidence intervals, and orientation before that scope is applied; if either
endpoint is outside the portable scope, both records are omitted together.
The adapter fails unless typed plus explicitly omitted records exactly reconcile
to the raw source count. The aggregate parsing summaries and their identities
are embedded in `release_c_extensions.json`. Scope omissions are also copied to
the final `omissions.json` with machine-readable raw, typed, non-BND, BND-record,
and reciprocal-BND-pair counts, but never with loci, alleles, sample identifiers,
or host paths.

DELLY's pinned source initializes `FORMAT/RDCN` to `-1` when a read-depth copy
number cannot be estimated (for example, when both flanking control counts are
unavailable), and its own downstream code treats that value as missing/invalid.
The adapter therefore maps exactly `-1` to a null copy number, records an
aggregate sentinel count, and preserves an explicit row limitation. Any other
negative `RDCN` remains invalid and fails closed.

## Shipped profiles

Run this after installation to prove that all profiles survived packaging and
still match the release boundary:

```bash
gab release-c profile-catalog
```

The catalog includes digest-pinned profiles for:

- Manta germline SV discovery;
- DELLY short-read SV discovery;
- CNVpytor combined read-depth and B-allele-frequency calling. The in-container
  HDF5 audit follows CNVpytor 1.3.2 `call_2d` signal semantics exactly: the
  corrected/masked parity-RD histogram and the unflagged binned parity/unique-RD
  histograms must all be present; a nonexistent GC-corrected unique-RD signal is
  not assumed. Because both bin sizes share one mutable HDF5, the portable
  10-kb audit is repeated after the 100-kb call and both published audits must
  bind the final `.pytor` bytes. Four-decimal TSV fields reconcile to the
  full-precision HDF5 value at an inclusive half-rendered-unit boundary;
  decimal comparison keeps that boundary exact and larger differences fail
  closed;
- ExpansionHunter with separate female and male ploidy models. Repeat motifs
  may use standard IUPAC DNA ambiguity symbols when the digest-pinned variant
  catalog does so; those rows are retained losslessly, forced to `REVIEW`, and
  carry `CATALOG:IUPAC_AMBIGUOUS_MOTIF` rather than being silently omitted;
- HLA*LA and SpecHLA short-read HLA typing. HLA*LA's exact resource-tree mount
  is the parent directory containing the named `PRG_MHC_GRCh38_withIMGT`
  subdirectory, because `--customGraphDir` is a parent path and HLA*LA appends
  `--graph` internally. SpecHLA's extraction and calling stages explicitly
  pass through the pinned image's `env-execute` bootstrap before invoking the
  public wrappers, because the runner replaces the image entrypoint. Its exact
  resource tree must also include the complete six-file Bowtie2 index for each
  of the two mapping FASTA prefixes;
- GATK Mutect2 mitochondrial normal/shifted calling. The pinned GRCh38 input
  reference names its mitochondrial contig `MT`; the workflow remaps and emits
  canonical portable mitochondrial coordinates as `chrM`. Copy
  `config/mt-call-policy.example.json` unchanged when using the shipped caller
  profile, because the profile verifies this policy byte-for-byte.
- bounded local pileup capture in a minimal pysam OCI runtime.

Cyrius and Aldy support is output-only and licence-gated. Their adapters can
validate locally acquired caller output, but the project does not download,
redistribute, or silently accept those caller licences.

## Private manifest

Copy `config/release-c.example.yaml` to an ignored local path and replace every
placeholder with the actual files. The configuration may contain identifiers
and paths; it must not be committed.

For a normal run, hash every bound file while creating the manifest:

```bash
gab release-c manifest create \
  --config local/release-c.yaml \
  --output local/release-c-manifest.json
```

For very large alignments, `--reuse-supplied-alignment-sha256` may defer only
the alignment re-hash. The resulting manifest is explicitly marked partial,
and the input gate must perform the full re-hash immediately before every
first caller execution. Indexes, reference files, and small VCF inputs are
never deferred.

Validate the private schema without opening files:

```bash
gab release-c manifest validate-structure \
  --manifest local/release-c-manifest.json
```

Every logical subject must have exactly one technical dataset labelled
`primary`; any additional independently sequenced dataset is labelled
`technical_replicate`. This label selects provenance semantics only. It does
not permit copying a call from the primary source into a replicate or dropping
a discordant replicate row.

## One caller lifecycle

First compile an exact private mount specification from the input manifest,
the selected public profile, and a private resource-role map copied from
`config/release-c-resources.example.json`:

```bash
gab release-c binding-builder \
  --manifest local/release-c-manifest.json \
  --dataset-id TECHNICAL_DATASET_01 \
  --profile-id manta-germline-1.6.0 \
  --resource-bindings local/release-c-resources.json \
  --resource-manifest local/resources.json \
  --workspace work/caller-output/dataset-01/manta \
  --output work/bindings/dataset-01/manta.json \
  --user 1000:1000 \
  --minimum-free-bytes 50000000000
```

The builder derives alignment, index, BAF, and reference roles from the
private input manifest. Static binaries, catalogues, wrapper sources, and
resource trees must be declared by role in the resource map. Extra, missing,
or profile-digest-mismatched roles fail before a run plan is published.
The resource manifest used by `matrix-runner` must embed the exact output of
`gab release-c profile-catalog`; a freeze created for another code/resource
release is rejected before any caller starts. Shipped profiles pin the exact
GRCh38 FASTA/index. HLA*LA, SpecHLA, CNVpytor, and mitochondrial directory
mounts additionally pin the SHA-256 of their complete tree manifests.
Matrix execution accepts only the canonical output of `resource-manifest
create`, including its verified `manifest_identity_sha256` and explicit
`private_paths_persisted=false` boundary. A legacy or hand-written freeze is
rejected even when it happens to contain the current profile catalog.

After independently validating the deployment's tools, resource trees, and
positive/negative controls, copy
`config/release-c-resource-metadata.example.json` to an ignored local path and
record their path-free versions, identities, receipts, and every remaining
limitation. Bind that metadata deterministically to the installed profiles:

```bash
gab release-c resource-manifest create \
  --metadata local/release-c-resource-metadata.json \
  --output local/tool-resource-freeze.json

gab release-c resource-manifest validate \
  --manifest local/tool-resource-freeze.json
```

This command does not turn an assertion into a scientific control: the
metadata must point by digest to controls that were actually run. It prevents
profile drift, path leakage, silent warning removal, and accidental reuse of a
freeze from another release. `production_dataset_execution_claimed` must stay
false because production completion is proven later by lifecycle receipts.

The restart-safe lifecycle command performs the input gate, profile
compilation, isolated OCI execution, receipt finalization, and validation as
one no-overwrite checkpoint chain. Re-running the same command validates and
reuses completed checkpoints; a retained failure requires a new work
directory.

```bash
gab release-c caller-lifecycle \
  --manifest local/release-c-manifest.json \
  --dataset-id TECHNICAL_DATASET_01 \
  --profile-id manta-germline-1.6.0 \
  --bindings local/manta-bindings.json \
  --resource-manifest local/resources.json \
  --work-directory work/lifecycles/dataset-01/manta
```

The private work directory contains the executable plan, logs, and raw caller
artifacts. Its `lifecycle-summary.json` is path-free and identifier-free, but
it should still pass the project privacy gate before publication.

After a lifecycle reaches `PASS`, materialize its receipt-bound raw artifact
with the adapter selected by the exact shipped profile:

```bash
gab release-c profile-adapter materialize \
  --profile-id manta-germline-1.6.0 \
  --lifecycle-directory work/lifecycles/dataset-01/manta \
  --bindings work/bindings/dataset-01/manta.json \
  --technical-source-role primary \
  --output-directory work/typed/dataset-01/manta

gab release-c profile-adapter validate \
  --profile-id manta-germline-1.6.0 \
  --lifecycle-directory work/lifecycles/dataset-01/manta \
  --bindings work/bindings/dataset-01/manta.json \
  --technical-source-role primary \
  --output-directory work/typed/dataset-01/manta
```

The adapter re-hashes every raw artifact component against the caller receipt,
checks every typed row against the exact caller run, reference, alignment,
resource manifest, and full raw-artifact identity, then publishes a
no-overwrite directory. HLA adapters additionally require a private options
file copied from `config/release-c-adapters.example.json`; its database and
nomenclature identities become evidence fields, while its local path does not.

The following outline uses generic private paths. Exact mount roles are listed
inside the selected profile and must all be supplied in `bindings.json`.

```bash
gab release-c input-gate \
  --manifest local/release-c-manifest.json \
  --dataset-id TECHNICAL_DATASET_01 \
  --caller manta \
  --output work/gates/manta.json \
  --minimum-free-bytes 50000000000

gab release-c caller-profile prepare \
  --profile-id manta-germline-1.6.0 \
  --bindings local/manta-bindings.json \
  --input-gate work/gates/manta.json \
  --resource-manifest local/resources.json \
  --output-directory work/plans/manta

gab release-c oci-runner run \
  --boundary work/plans/manta/prepared-boundary.json \
  --execution-spec work/plans/manta/execution-spec.json \
  --report work/reports/manta-execution.json \
  --log-dir work/logs/manta

gab release-c caller-receipt finalize \
  --boundary work/plans/manta/prepared-boundary.json \
  --input-gate work/gates/manta.json \
  --resource-manifest local/resources.json \
  --execution-report work/reports/manta-execution.json \
  --artifact-spec work/plans/manta/artifact-spec.json \
  --output work/receipts/manta.json
```

Use the caller-specific adapter to materialize typed evidence, then stage all
domains for one logical subject:

```bash
gab release-c extension-stage stage \
  --private-manifest local/release-c-manifest.json \
  --dataset-id TECHNICAL_DATASET_01 \
  --domain structural_variants=work/typed/manta.parquet \
  --input-gate work/gates/manta.json \
  --caller-receipt work/receipts/manta.json \
  --adapter-summary work/typed/adapter-summary.json \
  --output output/release-c-extensions

gab release-c extension-stage validate \
  --stage output/release-c-extensions
```

All final paths use no-overwrite publication. A failed or interrupted run
retains its private workspace for diagnosis but cannot masquerade as a final
portable extension.

## Enriching a .genome/1.1 bundle

After the subject-bound extension stage passes, create a new enriched bundle.
The base bundle is never modified. `copy` is the portable default;
`hardlink` is available for an immutable base and same-filesystem output when
storage is constrained. The chosen materialization mode and exact base
manifest digest are recorded in transformation provenance.

```bash
gab release-c bundle-extension build \
  --base-bundle output/base/sample.genome \
  --release-c-stage output/release-c-extensions \
  --output output/enriched/sample.genome \
  --reference-fasta local/GRCh38.fa \
  --materialization copy \
  --minimum-free-bytes 100000000000

gab release-c bundle-extension validate \
  --bundle output/enriched/sample.genome \
  --reference-fasta local/GRCh38.fa
```

Publication occurs only after the core `.genome/1.1` validator, every typed
Release-C schema, caller-receipt linkage, field-glossary coverage, omission
accounting, manifest/inventory digests, and a bounded `read-backed` query all
pass. The builder replaces the old generic read-backed omission with explicit
per-domain omissions, so an unconfigured caller can never look like a normal
or negative result. Explicit portable-scope omissions reconcile to retained raw
artifacts, and caller/adapter receipts remain present even for zero typed rows.
Caller and technical-replicate rows remain distinct.

Inside an enriched bundle, search all read-backed domains with:

```bash
python readmygenome/scripts/query.py read-backed chr7:117559593 --bundle .
python readmygenome/scripts/query.py read-backed HLA-B --bundle .
python readmygenome/scripts/query.py read-backed CYP2D6 --bundle .
```

The query is output-bounded and returns the source caller, technical-source
role, concordance state, confidence tier, limitations, and confirmation flag.
It does not infer a diagnosis or a negative result from absence.

## Bounded local read-back

Local pileup review is intentionally restricted to loci that need extra
evidence. The target selector reads a validated `.genome/1.1` bundle and
prioritizes clinical findings, de-novo and mosaicism candidates, suspected
dropout, pedigree/LOCO conflicts, technical-replicate discordance, and
independent-caller discordance. It keeps every clinical finding, applies a
deterministic hash order within each priority tier, and emits the selected
alleles in canonical genomic order. The default limit is 20,000 targets and
the hard safety limit is 100,000.

```bash
gab release-c readback-targets create \
  --bundle output/sample.genome \
  --output local/readback-targets.tsv \
  --receipt local/readback-targets.receipt.json \
  --maximum-targets 20000

gab release-c readback-targets validate \
  --targets local/readback-targets.tsv \
  --receipt local/readback-targets.receipt.json \
  --bundle output/sample.genome
```

Builder 0.2.8 publishes a typed, zero-row clinical findings table when a
completed ACMG evaluation produces no decisions. For an older native bundle
that omitted the table, the selector checks the registered inventory,
transformation and policy provenance before setting
`source_tables.clinical_findings_state` to `legacy_evaluated_empty` and
`clinical_findings_sha256` to null. An evaluated zero leaves family review
signals eligible for selection. Each variant linked by a finding remains
mandatory, including findings from nonactionable decisions.

The selector rejects a source whose family-evidence cardinality no longer
matches `variants.parquet`, refuses to truncate clinical findings, and uses
no-overwrite publication. The receipt binds the exact bundle manifest,
clinical table, family-evidence tree, selection policy, and target bytes while
persisting neither a subject identifier nor a local path. The target TSV and
receipt still reveal genotype-specific coordinates and therefore remain
private.

Bind the target TSV to the `readback-bounded-pysam-0.24.0-v2` profile, run one
caller lifecycle for each independent technical alignment, and feed the raw
capture JSON to `readback-adapter`. A target absent from the capture, a
no-call, or insufficient fragment evidence is never converted to a reference
genotype. A canonical header-only target TSV still runs the caller against the
indexed alignment and reference. The caller produces a raw capture with zero
records. The adapter records zero typed rows, zero targets with no eligible
reads and the raw artifact digest. The optional evidence table is omitted.
The v1 catalog and its original capture remain available for historical runs.

When two independent alignments use the same target manifest, annotate their
finite local observations without averaging or voting:

```bash
gab release-c readback-concordance \
  --input-parquet work/typed/source-01/readback.parquet \
  --input-parquet work/typed/source-02/readback.parquet \
  --output-parquet work/concordance/readback.parquet \
  --reference-sha256 REPLACE_WITH_EXACT_GRCH38_SHA256
```

The comparison uses Wilson intervals over eligible REF/ALT fragments. Missing,
unsupported, low-depth, or policy-mismatched observations remain explicitly
not comparable; non-overlapping intervals are marked discordant and both rows
are retained for review.

### Mitochondrial replicate fractions

For each exact mitochondrial allele, `gab release-c mt-concordance` compares the
two technical observations using Wilson intervals on the target ALT depth divided
by the sum of the native record's allele depths. A multiallelic record contributes
its reference depth once and the depth of every ALT, including filtered sibling
alleles. The caller's FORMAT/AF and the original counts remain in the output.
A singleton stays source-specific.

The subject stage verifies the complete typed adapter output. The standalone CLI
also recognizes a complete typed partition from its adjacent `adapter-summary.json`
and checks its file digest and row inventory. For other inputs the supplied
allele depths must account for FORMAT/DP to resolve the denominator. A missing
sibling with unexplained depth gives `unresolved_allelic_depth`; low depth after
resolution gives `insufficient_allelic_depth`. Unassigned reads can make DP
larger than the sum of AD in a complete native partition.

## Sequential caller matrix

The default matrix runs seven fixed profiles plus the independently confirmed
XX/XY ExpansionHunter profile for every technical dataset. Copy
`config/release-c-matrix.example.json` to an ignored local path to select a
smaller validated subset. Profiles are always executed one at a time.

```bash
gab release-c matrix-runner \
  --manifest local/release-c-manifest.json \
  --resource-bindings local/release-c-resources.json \
  --resource-manifest local/tool-resource-freeze.json \
  --adapter-options local/release-c-adapters.json \
  --matrix-config local/release-c-matrix.json \
  --work-root work/release-c-matrix \
  --runtime docker \
  --user 1000:1000 \
  --minimum-free-bytes 100000000000 \
  --indexed-probe-bases 67108864
```

Before every binding, caller, and adapter boundary the runner checks the free
space floor. Each job performs binding recomputation, full alignment re-hash,
64-MiB indexed random-access probe, network-disabled caller lifecycle, caller
receipt, and strict typed adapter. A one-base probe is suitable only for a
plan-compilation smoke test, not a production matrix.
Power-loss replay validates and reuses complete checkpoints. A caller failure
is retained; a fresh attempt is created only after review and explicit
`--retry-failed`, without deleting the prior attempt. If the caller lifecycle
completed and only strict adapter materialization failed, the retry revalidates
and reuses that immutable caller receipt in the same attempt; it does not rerun
the expensive caller.

For a fresh HLA*LA run, the matrix can reclaim the two temporary paired FASTQ
files produced by Picard. It tracks them in the new attempt workspace and
removes them after validating the typed HLA result. The extracted BAM, unpaired
reads, raw HLA table and typed result stay available for inspection. Failed
attempts, existing workspaces and validation-only replays retain their files.

After completion, validate every checkpoint without creating new work:

```bash
gab release-c matrix-runner \
  --manifest local/release-c-manifest.json \
  --resource-bindings local/release-c-resources.json \
  --resource-manifest local/tool-resource-freeze.json \
  --adapter-options local/release-c-adapters.json \
  --matrix-config local/release-c-matrix.json \
  --work-root work/release-c-matrix \
  --runtime docker \
  --user 1000:1000 \
  --minimum-free-bytes 100000000000 \
  --indexed-probe-bases 67108864 \
  --validate-only
```

The plan and runtime ledgers are private. The final summary contains aggregate
counts and public profile IDs only: no dataset/subject identifiers, local
paths, or private-manifest digest. It explicitly records zero removed and zero
collapsed portable typed caller/replicate rows, while adapter parsing summaries
separately account for any raw records outside the declared portable scope.

After the full matrix is complete, assemble one logical subject by explicitly
listing all and only that subject's technical dataset IDs. The assembler
revalidates the matrix identity, typed files, gates, receipts, and common
reference; it refuses datasets belonging to different logical subjects.

```bash
gab release-c subject-stage build \
  --manifest local/release-c-manifest.json \
  --matrix-work-root work/release-c-matrix \
  --dataset-id TECHNICAL_DATASET_01 \
  --dataset-id TECHNICAL_DATASET_02 \
  --work-directory work/subject-stage-attempt-01 \
  --output-directory output/subject-release-c \
  --minimum-free-bytes 100000000000

gab release-c subject-stage validate \
  --stage output/subject-release-c \
  --reference-sha256 REPLACE_WITH_EXACT_GRCH38_SHA256
```

Manta/DELLY, CNV configurations, HLA callers, and technical replicates are
compared with domain-specific rules. Every portable typed source row is
retained. A single source remains explicitly not assessed for concordance;
disagreement becomes review evidence and never triggers majority voting or an
intersection drop.
The private work directory records dataset IDs and therefore stays ignored;
the published extension stage contains neither them nor the logical-subject
identifier.

After every enriched bundle has passed validation, build deterministic
project-level manifests, checksums, and JSON/TSV/HTML QC reports. Bundle IDs
are supplied by the local operator; the command never derives or embeds host
paths.

```bash
gab release-c project-release build \
  --bundle subject-01=output/subject-01.genome \
  --bundle subject-02=output/subject-02.genome \
  --output-root output \
  --reference-fasta local/GRCh38.fa \
  --generated-at 2026-01-01T00:00:00Z \
  --minimum-free-bytes 100000000000

gab release-c project-release validate \
  --bundle subject-01=output/subject-01.genome \
  --bundle subject-02=output/subject-02.genome \
  --output-root output \
  --reference-fasta local/GRCh38.fa
```

The build publishes `SHA256SUMS`, `manifests/`, `reports/`, and a deterministic
path-free `logs/project_release.jsonl` without overwriting differing bytes. A
replay validates and reuses byte-identical metadata. The full validator
rechecks every enriched bundle, its pinned `.genome/1.1` schema, query smoke,
Release-C receipts, project manifest identities, the semantic identity of the
JSONL audit log, and every top-level checksum.

Each bundle also receives an identity-bound
`reports/bundle_validation/<ID>.json` containing its declared file hashes,
complete core and Release-C validator result, annotation-source provenance,
query regression, warnings, and explicit empty error set. Project QC aggregates
typed concordant, discordant, not-comparable, source-specific, and technical-
replicate row counts without persisting technical dataset identifiers.

The same report now performs a bounded-memory scan of each validated bundle and
publishes per-bundle plus project-wide annotation QC. The JSON and long-form TSV
include source/output reconciliation, SNV/INS/DEL/MNV/SV and genotype counts,
Ti/Tv, fixed indel-length bins, VEP consequence and gene/transcript/HGVS
coverage, ClinVar match/conflict/review distributions, gnomAD and applicable-
class predictor completeness, difficult-region stratification, PGx and PGS
accounting, target callability states, subject-only family states, technical-
replicate concordance, read-backed categories, and strict-adapter raw/typed/
scope-omission reconciliation. The HTML is a compact
overview of the same identity-bound counters. No locus, allele, genotype, or
technical dataset identifier is copied into project QC. Validation recomputes
the aggregate arithmetic and, on a full replay, recomputes every bundle QC
document from its Parquet data; editing and re-signing only one report therefore
does not make it valid.

For projects with at least three bundles, fixed public rules compare selected
QC fractions and counts with the cohort median. An unexplained outlier stops
publication. After investigation, the operator may supply an ignored local
`--qc-explanations` JSON document containing only the flagged portable bundle
ID, metric, uppercase explanation code, and a short non-sensitive explanation.
Stale explanations and explanations for non-outliers are rejected; the local
file path is never persisted.

## Technical replicates

Run each technical dataset independently through input gate, caller, receipt,
and adapter. Concordance modules compare only independently bound alignments.
A replay of the same alignment is not an independent replicate.

For multiple technical datasets belonging to one logical subject:

- keep every portable typed source row and source-alignment digest;
- label exact, interval-compatible, discordant, not-comparable, and
  source-specific states explicitly;
- do not average allele fractions or copy-number estimates;
- do not interpret a missing row as a negative result;
- stage the combined typed rows only after source-specific validation passes.

## Operational rules

- Keep private configs, bindings, manifests, alignments, resources, caller
  output, and logs outside Git.
- Use the exact reference and contig convention declared in the manifest.
- Select the ExpansionHunter profile that matches the independently confirmed
  sex/ploidy model.
- Validate every resource tree with `resource-tree` before profile compilation.
- Do not enable licence-gated callers without explicit licence acceptance.
- Re-run the input gate after any input, index, reference, or resource change.
- Never publish an extension-stage directory unless its independent validator
  returns `PASS`.
