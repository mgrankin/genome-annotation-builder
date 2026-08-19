# Contributing

Contributions must preserve three invariants:

1. missing evidence is never promoted to a reference or benign conclusion;
2. normalized output alleles and genotypes remain traceable to immutable source
   records;
3. no real genomic or identifying data enters Git history, logs, fixtures, or
   bug reports.

Run the unit/integration suite, lint checks, Git-index privacy scan, and
wheel/sdist payload verification before submitting a change. New
transformations need a deterministic synthetic test covering positive,
reference, no-call, and uncallable cases.

Public knowledge snapshots are immutable release inputs. Updating the ACMG SF
snapshot requires retaining the new source workbook outside Git, updating its
compiled-in byte/SHA identity, regenerating the logo-free canonical JSON with
`scripts/extract_acmg_sf.py`, documenting license/attribution and policy IDs,
and reviewing the full semantic diff. Never update a floating resource in
place under an existing pipeline version.
