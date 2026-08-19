# Evidence model

The builder uses a monotonic evidence model: transformations may add evidence,
normalization, or interpretation, but cannot erase the source observation.

| State | Genotype claim | Callability claim | Allowed absence answer |
|---|---|---|---|
| carried concrete ALT | non-reference | record-local | no |
| complete `0` at site | reference at exact site | record-local | yes, exact locus only |
| complete `0` symbolic block | reference throughout represented block | threshold-dependent | yes, when thresholds pass |
| incomplete genotype | unknown | no | no |
| interval-only callability | none | callable | no |
| no covering record | none | unknown | no |

Overlaps are resolved conservatively: a carried variant wins over reference;
no-call wins over interval-derived reference for an absence query; conflicting
source records fail validation unless their provenance documents an intentional
merge.

For configured target catalogues, the full order is `carried core or
symbolic/non-core variant > overlapping no-call > complete reference evidence
> covering interval > unknown`. The target proof extension records the exact
catalogue row digest, every selected evidence identifier, observed ploidy,
DP/GQ thresholds, and the resolution reason. A different carried ALT or a
symbolic event overlapping the target is not an explicit reference result.

Polygenic scoring is a consumer of this same proof, not a shortcut around it.
Only `non_reference_observed` for the exact target allele or a passing
`reference_observed` state can produce a dosage. The following states always
produce `dosage=null` and `contribution=null`: explicit no-call, another ALT,
low-confidence reference, callable interval without a genotype, uncallable
interval, unknown locus, unsupported allele/model row, or ambiguous ploidy.
The model audit keeps their absolute weight in the missing denominator.

A valid raw PGS is still not a normative percentile. Publication requires an
independently digest-bound calibration CDF for the same scoring file and a
named reference population. Absence of calibration means “valid uncalibrated
raw score retained in audit,” not average risk. PGS output is never
`clinical_grade` by default.

GWAS projection is an exact external-evidence join, not phenotype inference.
The carried allele must match the frozen resource on GRCh38 coordinate, rsID
alias, and effect allele; effect dosage is then calculated from the unchanged
subject GT. The legacy bulk effect column does not establish whether a value is
a beta or ratio. Rows lacking an explicit textual effect type remain in the
lossless audit but do not enter the normative typed association table. An exact
match means only that a population study reported an association; it does not
establish personal phenotype, causality, penetrance, or clinical actionability.
Semantically duplicate rows in the frozen upstream release remain distinct
source observations: their stable identifiers include the unique source-row
number as well as the normalized row digest and matched variant. They are not
silently collapsed into one assertion.

Pharmacogenomic matching is another consumer of the target proof. The official
PharmCAT representation and normalized evidence target are linked explicitly.
Only an exact passing carried genotype or complete passing explicit-reference
evidence is eligible for the matcher. A no-call, other ALT, low-confidence call,
interval-only state, unknown state, or structurally complex locus is written as
missing in the matcher input and remains fully described in the locus audit.

A PharmCAT gene result is not automatically a normative PGx result. Publication
requires exactly one complete non-complex diplotype, a phenotype, and an exact
non-provisional CPIC A-D recommendation mapping. Ambiguity, missing positions,
missing phenotype, absent normative mapping, or need for an outside complex-
locus caller leaves the gene audit-only. This means “not published under the
recorded gate,” never “normal metabolizer” or “standard dose.”

Computational predictor scores are independent, applicability-bounded evidence.
LOFTEE applies to predicted loss-of-function transcript consequences and its
HC/LC result is accepted only when the pinned source patch and real indexed
ancestral-allele lookup pass in the immutable VEP runtime; REVEL and
AlphaMissense apply to missense SNVs; CADD applies only to the pre-scored SNV
and small-indel classes in its frozen release; the pinned Ensembl SpliceAI
resource applies only to MANE-overlapping SNVs. A score may be absent because a
variant is outside that scope, because no applicable transcript matched, or
because the frozen resource has no row. These states are unavailable evidence,
not benign votes. In particular, the operational zero-record SpliceAI indel
input is an interface sentinel and never evidence that an indel has zero splice
effect.

The `consensus_damaging` convenience Boolean uses only populated REVEL,
AlphaMissense, and CADD scores and requires the declared number of threshold
votes. Null scores cast no vote and do not reduce the threshold. LOFTEE flags
remain transcript evidence rather than being silently converted into the same
numeric vote. Portable predictor claims are valid only when VEP command
semantics, the common resource-manifest digest, every artifact digest, coverage
metadata, and explicit omissions all agree.

Clinical evidence classes are kept separate from computational prediction,
population frequency/association, family evidence, and raw read support. A negative or
missing value in one class cannot be used as evidence from another class.

Family evidence is monotonic and privacy-minimized. It may describe observed
parental origin, Mendelian compatibility, phase, independent-caller or
technical-replicate concordance, difficult context, and review candidates, but
it cannot rewrite the accepted subject allele, genotype, ploidy, or phase. Each
family row is joined one-to-one to `variants.parquet`, and validation requires
exact subject `GT`/phase equality. Relative identifiers and relative genotypes
remain in the separately retained digest-bound family source, never in the
portable subject bundle.

The join follows immutable normalization provenance rather than assuming one
family source row per normalized output row. A complex source allele may
therefore annotate multiple atomized rows only when every projection is proved
and reconciled. An incomplete source genotype remains a no-call even when one
ALT allele is visible; it cannot create either a carried-variant or reference
claim.

A leave-one-out review file can have exception-only semantics. In that case an
exact retained exception is evidence for its recorded review category, while a
missing row means only “not individually persisted.” It does not prove
Mendelian consistency. Cohort-level consistency fractions are provenance and
cannot be applied to an individual absent locus. Source ordering is not
semantic and is canonicalized deterministically after uniqueness and summary
accounting checks. When a retained LOCO child genotype disagrees with the
accepted core genotype, the effective evidence is
`source_genotype_discordant`; the original source status remains audit data and
the accepted genotype is unchanged.

Clinical retrieval uses four deliberately separate layers:

- `clinical_evidence.parquet` records versioned external assertions and policy
  evidence in the normative `.genome/1.1` shape;
- `extensions/clinvar_assertions.parquet` preserves richer source fields without
  redefining the normative table;
- `clinical_findings.parquet` indexes pipeline-selected candidates;
- `interpretations.parquet` states the applied model, review status, evidence
  snapshot, and confirmation requirement explicitly.

ClinGen gene-level context is a separate knowledge layer:

- `extensions/clingen_gene_disease.parquet` contains expert-curated
  gene-disease validity assertions;
- `extensions/clingen_dosage_sensitivity.parquet` contains GRCh38 gene dosage
  sensitivity context.

Neither table classifies a carried allele or subject. They do not enter
`clinical_evidence.parquet`; exact gene membership is not pathogenicity,
actionability, or diagnosis evidence. No matching row means unavailable
context, not a negative finding.

An algorithmic interpretation is never expert review. In automated builds the
reviewer fields remain null, `clinical_grade` remains false, and orthogonal or
professional confirmation remains required even when `is_actionable` marks a
policy candidate.

The normative `.genome/1.1` schema requires `is_rare` to be Boolean. This
pipeline uses its conservative `null_treatment=rare` option: an allele without
a usable population frequency is retained by rarity filters, but that Boolean
is not evidence that rarity was observed. Consumers must inspect the nullable
frequency fields before making a population-frequency claim.
