from pathlib import Path

import pytest

from genome_annotation_builder.config import load_config
from genome_annotation_builder.errors import ConfigurationError

BASE = """
project:
  output_dir: output
  work_dir: work
  bundle_name: synthetic.genome
input:
  vcf: synthetic.vcf.gz
  sample: SAMPLE_ID
  genome_build: GRCh38
reference:
  fasta: synthetic.fa
  expected_sha256: null
  contig_style: bare
resources:
  cache_dir: resources
  allow_download: false
  genome_spec_commit: 7354ae288784234f854562b7dceb87f067fad550
"""


def test_load_config_resolves_relative_paths(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(BASE, encoding="utf-8")
    config = load_config(path)
    assert config.input_vcf == tmp_path / "synthetic.vcf.gz"
    assert config.bundle_path == tmp_path / "output" / "synthetic.genome"
    assert config.quality.high_min_dp == 20
    assert config.include_sample_identifier is False


def test_unknown_key_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(BASE + "unexpected: true\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="unknown keys"):
        load_config(path)


def test_non_grch38_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(BASE.replace("GRCh38", "GRCh37"), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="only GRCh38"):
        load_config(path)


def test_enabled_vep_container_must_be_digest_pinned(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  vep:
    enabled: true
    mode: docker
    container_image: ensemblorg/ensembl-vep:release_116.0
    cache_version: "116"
    expected_version: "116.0"
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="pinned"):
        load_config(path)


def test_acmg_requires_vep_and_clinvar(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  acmg_sf:
    enabled: true
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="requires both"):
        load_config(path)


def test_acmg_review_threshold_is_bounded(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  acmg_sf:
    min_clinvar_review_stars: 5
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="integer from 0 to 4"):
        load_config(path)


def test_predictors_require_vep_and_preserve_private_relative_paths(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    predictor = """
annotation:
  predictors:
    root: private-predictors
    cadd:
      enabled: true
      license: "CADD v1.7 terms; non-commercial use unless separately licensed"
      redistribution_allowed: false
"""
    path.write_text(BASE + predictor, encoding="utf-8")
    with pytest.raises(ConfigurationError, match=r"requires annotation\.vep"):
        load_config(path)

    predictor = predictor.replace(
        "annotation:\n",
        "annotation:\n  vep:\n    enabled: true\n    mode: docker\n    plugin_dir: /plugins\n",
    )
    path.write_text(BASE + predictor, encoding="utf-8")
    config = load_config(path)
    assert config.predictors.enabled_plugins == ("cadd",)
    assert config.predictors.root == tmp_path / "private-predictors"
    assert config.predictors.manifest == (
        tmp_path / "private-predictors" / "predictor_resource_manifest.json"
    )


def test_predictor_resource_redistribution_cannot_be_asserted(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  vep:
    enabled: true
    mode: docker
  predictors:
    alphamissense:
      enabled: true
      license: "CC BY-NC-SA 4.0 notice embedded in exact file; upstream repository states CC BY 4.0"
      redistribution_allowed: true
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="cannot be asserted"):
        load_config(path)


@pytest.mark.parametrize(
    ("vep", "message"),
    [
        (
            "enabled: true\n    mode: executable\n    plugin_dir: /plugins",
            "pinned Docker VEP runtime",
        ),
        (
            "enabled: true\n    mode: docker\n    plugin_dir: /custom-plugins",
            "plugin directory /plugins",
        ),
    ],
)
def test_predictors_require_plugin_code_from_the_pinned_container(
    tmp_path: Path, vep: str, message: str
) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + f"""
annotation:
  vep:
    {vep}
  predictors:
    cadd:
      enabled: true
      license: "CADD v1.7 terms; non-commercial use unless separately licensed"
      redistribution_allowed: false
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match=message):
        load_config(path)


def test_gwas_catalog_configuration_is_dated_and_digest_pinned(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  gwas:
    enabled: true
    snapshot: "2026-08-03"
    archive_expected_sha256: 3ed9a1335926181d36334fe66fbf6754faf5120c780bc385241c5ad3fafcfe47
    license: EMBL-EBI Terms of Use
    redistribution_allowed: false
""",
        encoding="utf-8",
    )
    config = load_config(path)
    assert config.gwas.enabled is True
    assert config.gwas.snapshot == "2026-08-03"
    assert config.gwas.redistribution_allowed is False

    path.write_text(
        path.read_text(encoding="utf-8").replace("2026-08-03", "2026-08-04", 1),
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="requires GWAS Catalog snapshot"):
        load_config(path)


def test_enabled_target_catalog_requires_digest(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
callability:
  target_catalog:
    enabled: true
    path: synthetic-targets.tsv
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="requires expected_sha256"):
        load_config(path)


def test_target_catalog_configuration_is_typed(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
callability:
  target_catalog:
    enabled: true
    path: synthetic-targets.tsv
    expected_sha256: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
    assay_scope: test_assay
""",
        encoding="utf-8",
    )
    config = load_config(path)
    assert config.target_catalog_enabled is True
    assert config.target_catalog_path == tmp_path / "synthetic-targets.tsv"
    assert config.target_catalog_assay_scope == "test_assay"


def test_family_evidence_configuration_is_digest_bound_and_generic(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  family_evidence:
    enabled: true
    joint_vcf: joint.vcf.gz
    joint_vcf_expected_sha256: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
    pedigree: family.ped
    pedigree_expected_sha256: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
    source_version: family-callset-v1
    snapshot_date: "2026-08-15"
    tags:
      confidence_format: FC
      reasons_format: FR
      difficult_info_flags: [DIFFICULT]
      de_novo_confidence_values: [DE_NOVO_CANDIDATE]
      conflict_confidence_values: [MENDELIAN_CONFLICT]
      parental_mosaicism_reason_tokens: []
      parental_dropout_reason_tokens: []
      child_dropout_reason_tokens: []
""",
        encoding="utf-8",
    )
    config = load_config(path)
    assert config.family_evidence.enabled is True
    assert config.family_evidence.joint_vcf == tmp_path / "joint.vcf.gz"
    assert config.family_evidence.tags.confidence_format == "FC"
    assert config.family_evidence.tags.difficult_info_flags == ("DIFFICULT",)


def test_enabled_family_evidence_requires_source_digests(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  family_evidence:
    enabled: true
    joint_vcf: joint.vcf.gz
    pedigree: family.ped
    source_version: family-callset-v1
    snapshot_date: "2026-08-15"
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="joint_vcf_expected_sha256"):
        load_config(path)


def test_family_loco_configuration_is_all_or_none(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  family_evidence:
    enabled: false
    loco:
      review_tsv: review.tsv
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="all-or-none"):
        load_config(path)


def test_pgs_configuration_is_digest_and_licence_explicit(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  pgs:
    enabled: true
    models:
      - pgs_id: PGS000001
        scoring_file: score.txt.gz
        expected_sha256: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        source_version: "2026-08-15"
        training_date: "2015-02-18"
        trait: Example trait
        model_type: linear
        license: CC0-1.0
        redistribution_allowed: false
        max_missing_weight_fraction: 0.01
""",
        encoding="utf-8",
    )
    config = load_config(path)
    assert config.pgs_enabled is True
    assert len(config.pgs_models) == 1
    assert config.pgs_models[0].pgs_id == "PGS000001"
    assert config.pgs_models[0].max_missing_weight_fraction == pytest.approx(0.01)
    assert config.pgs_models[0].redistribution_allowed is False


def test_pgs_requires_explicit_redistribution_decision(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(
        BASE
        + """
annotation:
  pgs:
    enabled: true
    models:
      - pgs_id: PGS000001
        scoring_file: score.txt.gz
        expected_sha256: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        source_version: "2026-08-15"
        training_date: "2015-02-18"
        trait: Example trait
        license: CC0-1.0
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="redistribution_allowed"):
        load_config(path)


PGX = """
annotation:
  pgx:
    enabled: true
    pharmcat:
      version: "3.4.0"
      release_commit: 41f21b71886b961b24a00344ef2dad2387aaa0c9
      positions_vcf: pharmcat_positions.vcf.bgz
      expected_sha256: 01a28fe31675c9e8a8b18bcab4afb82a62abdfea6606856aeca0e0378e7d6b58
      container_image: pgkb/pharmcat@sha256:bc498d55c33094002bd7e14221a6b3fcd7243a587675369bc1a5ecf6ac0e5319
      license: MPL-2.0
      redistribution_allowed: false
    cpic:
      pairs_json: cpic_pairs.json
      expected_sha256: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      source_version: "1.60.0"
      snapshot_date: "2026-08-15"
      license: reviewed-resource-terms
      redistribution_allowed: false
    docker_executable: docker
    timeout_seconds: 1800
"""


def test_pgx_configuration_is_immutable_and_licence_explicit(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(BASE + PGX, encoding="utf-8")
    config = load_config(path)
    assert config.pgx.enabled is True
    assert config.pgx.version == "3.4.0"
    assert config.pgx.positions_vcf == tmp_path / "pharmcat_positions.vcf.bgz"
    assert config.pgx.cpic_source_version == "1.60.0"
    assert config.pgx.redistribution_allowed is False
    assert config.pgx.cpic_redistribution_allowed is False


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        ('version: "3.4.0"', 'version: "3.3.0"', "requires PharmCAT"),
        (
            "bc498d55c33094002bd7e14221a6b3fcd7243a587675369bc1a5ecf6ac0e5319",
            "a" * 64,
            "container_image",
        ),
        ("license: MPL-2.0", 'license: ""', "license must be explicit"),
    ],
)
def test_pgx_rejects_unpinned_or_implicit_boundaries(
    tmp_path: Path, old: str, new: str, message: str
) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(BASE + PGX.replace(old, new, 1), encoding="utf-8")
    with pytest.raises(ConfigurationError, match=message):
        load_config(path)
