"""Immutable format and resource boundary used by the pipeline."""

from __future__ import annotations

PIPELINE_VERSION = "0.2.1"
SCHEMA_VERSION = "1.1.0"
GENOME_BUILD = "GRCh38"

GENOME_SPEC_REPOSITORY = "https://github.com/genome-computer/genome-spec"
GENOME_SPEC_COMMIT = "7354ae288784234f854562b7dceb87f067fad550"
GENOME_SPEC_RAW_ROOT = (
    f"https://raw.githubusercontent.com/genome-computer/genome-spec/{GENOME_SPEC_COMMIT}"
)

# Keys are paths at the pinned upstream revision. Values are exact SHA-256.
PINNED_UPSTREAM_FILES: dict[str, str] = {
    "SPEC.md": "e2aa151190bbb37a8304254dbf88dbb4d29eab1113655113c771ccffe7d1d0b1",
    "schema.json": "6d07eaa4bacca4bb2c474eacd8e1a658470475fc4a46438432fb4509d0f11615",
    "profiles/agent-exploration.json": (
        "128ab7bc42994dfeae8df7d38f14cfacfc7adb6e0ccb9cb65dfde84ea52ed6f9"
    ),
    "scripts/genome-convert": ("cab078eee144f24c69236ed03b4291f75c0fc2d384ff529e3eb5e6e5afb42d59"),
    "reference/reader.py": ("cd3262a02a66d47a31c8ecce3c9f2fc1101ac51772bd9e27167cd188bb19e4e5"),
    "LICENSE": "f1da28603af0a2e6c12ca2478ecd6bc36e4f4725e8a0414de453e93c0823a5c9",
}

BIOLOGICAL_BASES = frozenset("ACGT")
SYMBOLIC_REFERENCE_ALLELES = frozenset({"<*>", "<NON_REF>"})
STANDARD_CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"])

# PharmCAT/CPIC boundary for the first PGx implementation. PharmCAT release
# assets and the multi-architecture image are immutable at these identities.
PHARMCAT_VERSION = "3.4.0"
PHARMCAT_RELEASE_COMMIT = "41f21b71886b961b24a00344ef2dad2387aaa0c9"
PHARMCAT_RELEASE_DATE = "2026-07-14"
PHARMCAT_LICENSE = "MPL-2.0"
PHARMCAT_POSITIONS_URL = (
    "https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat_positions_3.4.0.vcf.bgz"
)
PHARMCAT_POSITIONS_BYTES = 11_027
PHARMCAT_POSITIONS_SHA256 = "01a28fe31675c9e8a8b18bcab4afb82a62abdfea6606856aeca0e0378e7d6b58"
PHARMCAT_POSITIONS_RECORDS = 1_226
PHARMCAT_POSITIONS_ALLELES = 1_288
PHARMCAT_DOCKER_IMAGE = (
    "pgkb/pharmcat@sha256:bc498d55c33094002bd7e14221a6b3fcd7243a587675369bc1a5ecf6ac0e5319"
)
PHARMCAT_DOCKER_AMD64_MANIFEST = (
    "sha256:5f742b225204dcfb63d1b8ba4ddd4b8e41d266ca692f2b027582ba31992b8a71"
)
STRUCTURALLY_COMPLEX_PGX_GENES = frozenset({"CYP2D6", "HLA-A", "HLA-B", "MT-RNR1"})

# Ensembl release 116 is the functional-annotation boundary for pipeline 0.1.
# The archive's SHA-256 is added to this allow-list after downloading the
# release artifact and independently checking Ensembl's published CHECKSUMS
# entry.  The older BSD checksum remains recorded because it is Ensembl's
# authoritative checksum for this archive.
VEP_VERSION = "116.0"
VEP_CACHE_VERSION = "116"
VEP_CACHE_ASSEMBLY = "GRCh38"
VEP_CACHE_ARCHIVE = "homo_sapiens_vep_116_GRCh38.tar.gz"
VEP_CACHE_URL = (
    f"https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/{VEP_CACHE_ARCHIVE}"
)
VEP_CACHE_CHECKSUMS_URL = (
    "https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/CHECKSUMS"
)
VEP_CACHE_CHECKSUMS_SHA256 = "9d4df34db16f24d96b0936ef61405759494ccbd21b5daa5c7e12f8583521556d"
VEP_CACHE_ARCHIVE_BYTES = 27_644_657_162
VEP_CACHE_ARCHIVE_BSD_SUM = 56_036
VEP_CACHE_ARCHIVE_1K_BLOCKS = 26_996_736
VEP_CACHE_ARCHIVE_SHA256 = "014c7dd9bb5ad06665866d62eb80f31ca761197bcb9b59280300676f996e600d"
VEP_DOCKER_AMD64 = (
    "ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11"
)

# Optional VEP predictor resources.  These are intentionally not vendored: the
# files are large and several have independent licence terms.  ``resources
# predictors-sync`` downloads exact release assets into the private resource
# cache, verifies both the upstream checksum (where one is published) and our
# stronger SHA-256 pin, then records only relative paths in a local manifest.

# LOFTEE GRCh38 branch, frozen at one exact source commit.  The source is MIT;
# the ancillary Broad-hosted GRCh38 data do not carry an explicit redistribution
# grant in the upstream README, so the builder conservatively treats them as
# non-redistributable inputs.
LOFTEE_COMMIT = "a46b502a68c812c8ae0c5a5721c0603fe81cae8d"
LOFTEE_SOURCE_ARCHIVE = f"loftee-{LOFTEE_COMMIT}.tar.gz"
LOFTEE_SOURCE_URL = f"https://codeload.github.com/konradjk/loftee/tar.gz/{LOFTEE_COMMIT}"
LOFTEE_SOURCE_BYTES = 1_605_763
LOFTEE_SOURCE_SHA256 = "fee0d112d5d5f6814bb69cafe8d3f92f3a29af0e626c5338bc9737fb88879c7d"
# The pinned upstream LoF.pm shells out to ``samtools faidx`` for ancestral
# allele lookup.  The immutable VEP 116 image intentionally does not ship that
# executable, and upstream LOFTEE only emits a Perl warning when the lookup
# fails.  Apply one reviewable, byte-pinned patch which uses the
# Bio::DB::HTS::Faidx module already present in that same image.  Both sides of
# the transformation are fixed so an upstream drift or fuzzy patch can never be
# accepted silently.
LOFTEE_PATCH_ID = "biodbhts-faidx-v1"
LOFTEE_PATCH_FILE = "loftee_biodbhts_faidx_v1.patch"
LOFTEE_PATCH_SHA256 = "29970270b19159b25fc11a8ab5f7ed62ce084d60d1e7e02efd1ed111b896e74a"
LOFTEE_LOF_PM_ORIGINAL_SHA256 = "649b15cbb2670284f0dfdbe6e3aa1601cdb1fec629e5ba686b412cddd05efe35"
LOFTEE_LOF_PM_PATCHED_SHA256 = "ffbb60bf93f7458f9b064fd7202ab76054cd5c8f5665062f46c880909d2a82c1"
LOFTEE_PATCHED_SOURCE_DIR = f"source-{LOFTEE_PATCH_ID}"
LOFTEE_LICENSE = "MIT (source); ancillary GRCh38 resource redistribution terms not stated"
LOFTEE_DATA_ROOT_URL = "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38"
LOFTEE_GERP_FILE = "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
LOFTEE_GERP_BYTES = 12_617_579_354
LOFTEE_GERP_SHA256 = "8801e57ce8effbef9248b122caee16da338bfa8319dd41392c9ce4e58ea6cfe8"
LOFTEE_ANCESTOR_FILE = "human_ancestor.fa.gz"
LOFTEE_ANCESTOR_BYTES = 885_231_320
LOFTEE_ANCESTOR_SHA256 = "624b0e0f8e1c4e7c8639edf17df03f213621d1f56f6fa1c8025fad56082b157b"
LOFTEE_ANCESTOR_FAI_BYTES = 736
LOFTEE_ANCESTOR_FAI_SHA256 = "703e9a2011886f90679e99d0333ca91afc81df726db91c76d57b9e2d6f182ec5"
LOFTEE_ANCESTOR_GZI_BYTES = 764_520
LOFTEE_ANCESTOR_GZI_SHA256 = "a31e4ee63e519a0da4a8e2750dabbca177129b7f312b2a993c17a99981e6a3bf"
LOFTEE_SQL_GZ_FILE = "loftee.sql.gz"
LOFTEE_SQL_GZ_BYTES = 29_960_948
LOFTEE_SQL_GZ_SHA256 = "527182b8a4f81d0f68d06fa311c0d22730423f2d419859c782f477ab67656441"

# CADD v1.7, GRCh38.  Upstream publishes MD5SUMs; SHA-256 is additionally
# frozen here after downloading the exact objects.  Commercial use requires a
# separate CADD licence and is never inferred by this software.
CADD_VERSION = "1.7"
CADD_LICENSE = "CADD v1.7 terms; non-commercial use unless separately licensed"
CADD_ROOT_URL = f"https://kircherlab.bihealth.org/download/CADD/v{CADD_VERSION}/GRCh38"
CADD_SNV_FILE = "whole_genome_SNVs.tsv.gz"
CADD_SNV_BYTES = 87_473_403_655
CADD_SNV_MD5 = "88577a55f1cd519d44e0f415ba248eb9"
CADD_SNV_SHA256 = "11582f442d3db5bc529a3786e9ddfe27860642c2447d64e468387ad9c050645a"
CADD_SNV_TBI_BYTES = 2_761_840
CADD_SNV_TBI_MD5 = "347df8fac17ea374c4598f4f44c7ce8b"
CADD_SNV_TBI_SHA256 = "bf21382da1faa72dfee9dd0905bc294b9e074f70b31daaed370f19b0db4ce8f4"
CADD_INDEL_FILE = "gnomad.genomes.r4.0.indel.tsv.gz"
CADD_INDEL_BYTES = 1_257_151_321
CADD_INDEL_MD5 = "4b9c685c96d396af4d001c2f7dd9d8f9"
CADD_INDEL_SHA256 = "19d5a784a99e5b0602089e6c527dbf9bba35d774fcc04b0d4669b0e00f1e8c8a"
CADD_INDEL_TBI_BYTES = 1_899_705
CADD_INDEL_TBI_MD5 = "85f3d2daa9202c5915c0ce0f1c749a66"
CADD_INDEL_TBI_SHA256 = "dd2a1fca35816a82bc18a68e13f4fc9fce5536a2e7ed7a11e36863c589019f7e"

# REVEL v1.3 all-chromosome archive from the immutable Zenodo record.  The
# record does not declare a machine-readable redistribution licence; retain a
# conservative non-commercial/non-redistributable boundary in output bundles.
REVEL_VERSION = "1.3"
REVEL_ARCHIVE = "revel-v1.3_all_chromosomes.zip"
REVEL_URL = "https://zenodo.org/api/records/7072866/files/revel-v1.3_all_chromosomes.zip/content"
REVEL_ARCHIVE_BYTES = 667_102_707
REVEL_ARCHIVE_MD5 = "3ea2bc33e6b5455fc7e9899da863b5fe"
REVEL_ARCHIVE_SHA256 = "815801d832e4ab76922dddfd901d7f24a755df57078a7138aa14f615da86a8a2"
REVEL_MEMBER = "revel_with_transcript_ids"
REVEL_MEMBER_BYTES = 6_496_023_640
REVEL_LICENSE = "REVEL upstream terms; academic/non-commercial use, redistribution not granted"

# AlphaMissense genomic predictions.  The official repository states CC BY
# 4.0, while this exact hosted file embeds a stricter CC BY-NC-SA 4.0 notice.
# The pipeline records that conflict and applies the stricter file notice.
ALPHAMISSENSE_VERSION = "2023"
ALPHAMISSENSE_FILE = "AlphaMissense_hg38.tsv.gz"
ALPHAMISSENSE_URL = f"https://storage.googleapis.com/dm_alphamissense/{ALPHAMISSENSE_FILE}"
ALPHAMISSENSE_BYTES = 642_961_469
ALPHAMISSENSE_MD5 = "9fd167735f16a1b87da6eb3e4c25fcb5"
ALPHAMISSENSE_SHA256 = "0516cfd71c0767ac8f9c469252d429000e94e02c008b6e3a46d4b4646fcd3475"
ALPHAMISSENSE_LICENSE = (
    "CC BY-NC-SA 4.0 notice embedded in exact file; upstream repository states CC BY 4.0"
)

# Ensembl-generated masked SpliceAI scores for every possible SNV overlapping
# MANE v1.4 transcripts.  Ensembl explicitly does not supply indel scores.  The
# adapter therefore creates a labelled empty indel sentinel for the VEP plugin
# and reports indel coverage as unavailable rather than sourcing an unpinned
# third-party mirror.
SPLICEAI_VERSION = "1.3.1"
SPLICEAI_SCOPE = "Ensembl MANE v1.4 masked SNV scores, GRCh38"
SPLICEAI_SNV_FILE = "spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz"
SPLICEAI_ROOT_URL = "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins"
SPLICEAI_SNV_URL = f"{SPLICEAI_ROOT_URL}/{SPLICEAI_SNV_FILE}"
SPLICEAI_SNV_BYTES = 28_643_031_420
SPLICEAI_SNV_SHA256 = "56a3d88ed14fdd8a7d8884ffdd1c745a64102258d3ccc2584916efd1099b062e"
SPLICEAI_SNV_TBI_BYTES = 1_266_506
SPLICEAI_SNV_TBI_SHA256 = "a64c62434c554f5c98a10f84dab8e5544f45b145cd28acf10411ec350e49c20d"
SPLICEAI_LICENSE = "Ensembl terms plus upstream SpliceAI terms; redistribution review required"

# Immutable dated ClinVar snapshot. NCBI's rolling `clinvar.vcf.gz` alias is
# intentionally not used by the pipeline.
CLINVAR_SNAPSHOT = "20260808"
CLINVAR_ARCHIVE = f"clinvar_{CLINVAR_SNAPSHOT}.vcf.gz"
CLINVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
CLINVAR_URL = f"{CLINVAR_BASE_URL}/{CLINVAR_ARCHIVE}"
CLINVAR_TBI_URL = f"{CLINVAR_URL}.tbi"
CLINVAR_MD5_URL = f"{CLINVAR_URL}.md5"
CLINVAR_ARCHIVE_BYTES = 193_118_277
CLINVAR_ARCHIVE_MD5 = "9ffe83758f93f2623f39e83ef1659e05"
CLINVAR_ARCHIVE_SHA256 = "a99917628ff65443740fa010bbe8fcc59f4a2bdce8fcfccf439c22bba8928975"
CLINVAR_TBI_BYTES = 610_187
CLINVAR_TBI_SHA256 = "bd1e1b19a57d5511ff11820da518690107b227e3cfa7b3a3f1cdab0995367e05"
CLINVAR_MD5_FILE_SHA256 = "fc4bb09dd1adea7cff8de0370d6e3629a9f59d00c2bd4b52b7f9860fb7e80a2d"

# ACMG Secondary Findings v3.3, accompanied by the live ClinGen curation
# workbook frozen on 2026-08-15. ClinGen curated content is CC0; the workbook
# is parsed into a deterministic, logo-free JSON snapshot before vendoring.
ACMG_SF_VERSION = "3.3"
ACMG_SF_SNAPSHOT = "20260815"
ACMG_SF_POLICY_PMID = "40568962"
ACMG_SF_POLICY_DOI = "10.1016/j.gim.2025.101454"
ACMG_SF_WORKBOOK_URL = (
    "https://docs.google.com/spreadsheets/d/"
    "1ecSUe0bJ-2gAMXBsrf5U1wMEDOmiIjtR_A5388uj7HY/export?format=xlsx"
)
ACMG_SF_WORKBOOK_NAME = "clingen_secondary_findings_resource_20260815.xlsx"
ACMG_SF_WORKBOOK_BYTES = 1_119_389
ACMG_SF_WORKBOOK_SHA256 = "ac29d564a27dd08a7fca0d0cb111d833a5275c248322e89e71e54b4ed296067c"
ACMG_SF_GENE_COUNT = 84
ACMG_SF_POLICY_ROW_COUNT = 100
ACMG_SF_CLINGEN_ROW_COUNT = 177
ACMG_SF_GUIDANCE_ROW_COUNT = 84
ACMG_SF_FREQUENCY_ROW_COUNT = 104
ACMG_SF_GENE_SET_SHA256 = "95b5e5692ac5741125cd1f121f2727541796e08cdee87fa5cffde4fffeac7da0"
ACMG_SF_CANONICAL_SHA256 = "e0292fddbcefa006c51d1a81a389ee60b2e1eb3db3b199e82e86f5765d02930a"

# ClinGen expert-curated gene-disease validity and GRCh38 gene dosage
# sensitivity.  The live gene-validity export and the short-retention dosage
# archive are transformed into one deterministic, vendored CC0 snapshot.  This
# makes an old pipeline release reproducible after either upstream URL changes
# or expires.  The source byte identities remain recorded for audit.
CLINGEN_GENE_LEVEL_SNAPSHOT = "2026-08-15"
CLINGEN_GENE_LEVEL_FILE = "clingen_gene_level_20260815.json"
CLINGEN_GENE_LEVEL_CANONICAL_BYTES = 3_571_392
CLINGEN_GENE_LEVEL_CANONICAL_SHA256 = (
    "3e3d8c95715bc8729999545df93bb5e8d0904de221c6d49cacafef61830356cf"
)
CLINGEN_GENE_DISEASE_SNAPSHOT = "2026-08-15"
CLINGEN_GENE_DISEASE_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
CLINGEN_GENE_DISEASE_BYTES = 1_117_506
CLINGEN_GENE_DISEASE_SHA256 = "474dc140863af4e19d2b6c300e76c26bcb2fcb6299906c08725011f82a16a4e4"
CLINGEN_GENE_DISEASE_ROWS = 3_663
CLINGEN_DOSAGE_SNAPSHOT = "2026-08-14"
CLINGEN_DOSAGE_URL = (
    "https://ftp.clinicalgenome.org/archive/20260814/ClinGen_gene_curation_list_GRCh38.tsv"
)
CLINGEN_DOSAGE_BYTES = 247_559
CLINGEN_DOSAGE_SHA256 = "bc615d879a81cb16130e52659510b4a6734e354c3db59c91c6bab5fc66619a92"
CLINGEN_DOSAGE_ROWS = 1_526
CLINGEN_LICENSE = "CC0-1.0"
CLINGEN_ATTRIBUTION = (
    "Curated content obtained from the Clinical Genome Resource (ClinGen); accessed 2026-08-15."
)

# Dated NHGRI-EBI GWAS Catalog literature-curated association release.  The
# rolling ``latest`` alias is used only to discover releases; production uses
# the immutable dated FTP path below.  The archive contains one TSV member.
GWAS_CATALOG_SNAPSHOT = "2026-08-03"
GWAS_CATALOG_ARCHIVE = "gwas-catalog-associations_ontology-annotated-full.zip"
GWAS_CATALOG_MEMBER = "gwas-catalog-download-associations-alt-full.tsv"
GWAS_CATALOG_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/08/03/" + GWAS_CATALOG_ARCHIVE
)
GWAS_CATALOG_ARCHIVE_BYTES = 72_788_405
GWAS_CATALOG_ARCHIVE_SHA256 = "3ed9a1335926181d36334fe66fbf6754faf5120c780bc385241c5ad3fafcfe47"
GWAS_CATALOG_MEMBER_BYTES = 738_628_769
GWAS_CATALOG_LICENSE = "EMBL-EBI Terms of Use"
GWAS_CATALOG_ATTRIBUTION = "NHGRI-EBI GWAS Catalog, release 2026-08-03"
