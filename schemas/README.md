# Schemas

Normative genome.computer files are not edited or silently vendored. Run
`gab resources sync` to download the immutable revision listed in
`genome_annotation_builder.constants`, verify every SHA-256, and place it in the
configured private resource cache. The exact `schema.json` bytes are copied
into each staged bundle and re-hashed before promotion.

Project-owned extension schemas are implemented as exact Arrow schemas in the
package and use their own namespace; they never replace the upstream core
schema. In particular, `extensions/target_callability_evidence.parquet`
provides the lossless proof behind each normative `callability.parquet` row.
Polygenic scoring adds `extensions/pgs_locus_audit.parquet` and
`extensions/pgs_model_audit.parquet`; only a score passing their completeness
and calibration gates is projected into normative `prs.parquet`. Model rows
carry an exact contributing-variant count plus a completeness flag for the
bounded convenience list; the locus audit is always the full relationship.
