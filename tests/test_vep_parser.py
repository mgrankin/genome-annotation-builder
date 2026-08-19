from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pysam
import pytest

from genome_annotation_builder.arrow_schemas import variant_schema
from genome_annotation_builder.errors import InputValidationError
from genome_annotation_builder.variant_build import _blank_variant_row
from genome_annotation_builder.vep_parser import VepParser, apply_vep_annotation

CSQ_FIELDS = (
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
    "HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
    "Existing_variation|DISTANCE|STRAND|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|"
    "APPRIS|TSL|CCDS|ENSP|HGNC_ID|SIFT|PolyPhen|DOMAINS|LoF|LoF_filter|LoF_flags|"
    "LoF_info|ALLELE_NUM|gnomADg_AF|gnomADg_NFE_AF|gnomADg_AFR_AF|gnomADg_EAS_AF|"
    "gnomADg_SAS_AF|gnomADg_AMR_AF|CADD_PHRED|REVEL|am_pathogenicity|"
    "SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL"
)


def _csq_row(**values: str) -> str:
    return "|".join(values.get(field, "") for field in CSQ_FIELDS.split("|"))


def _record(tmp_path: Path, rows: list[str]) -> tuple[pysam.VariantFile, pysam.VariantRecord]:
    path = tmp_path / "annotated.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=1000>\n"
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP annotations. Format: {CSQ_FIELDS}">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t10\t.\tA\tG\t.\tPASS\tCSQ={','.join(rows)}\n",
        encoding="utf-8",
    )
    handle = pysam.VariantFile(str(path))
    return handle, next(handle)


def test_transcript_selection_and_most_severe_are_independent(tmp_path: Path) -> None:
    mane = _csq_row(
        Allele="G",
        Consequence="missense_variant",
        IMPACT="MODERATE",
        SYMBOL="GENE1",
        Gene="ENSG000001.5",
        Feature_type="Transcript",
        Feature="ENST_MANE.2",
        BIOTYPE="protein_coding",
        EXON="2/7",
        HGVSc="ENST_MANE.2%3Ac.100A%3EG",
        HGVSp="ENSP_MANE.1%3Ap.Lys34Arg",
        Existing_variation="rs123&COSV1",
        STRAND="1",
        CANONICAL="YES",
        MANE_SELECT="NM_000001.3",
        APPRIS="P1",
        TSL="1",
        HGNC_ID="42",
        SIFT="deleterious(0.01)",
        PolyPhen="probably_damaging(0.99)",
        DOMAINS="Pfam:PF00001&PROSITE:PS00001",
        ALLELE_NUM="1",
        gnomADg_AF="0.002",
        gnomADg_NFE_AF="0.003",
        gnomADg_AFR_AF="0.010",
        CADD_PHRED="25.4",
        REVEL="0.81",
        am_pathogenicity="0.70",
        SpliceAI_pred_DS_AG="0.11",
        SpliceAI_pred_DS_DL="0.42",
    )
    severe_noncanonical = _csq_row(
        Allele="G",
        Consequence="stop_gained",
        IMPACT="HIGH",
        SYMBOL="GENE1",
        Gene="ENSG000001.5",
        Feature_type="Transcript",
        Feature="ENST_OTHER.1",
        BIOTYPE="protein_coding",
        HGVSc="ENST_OTHER.1%3Ac.88A%3EG",
        HGVSp="ENSP_OTHER.1%3Ap.Lys30Ter",
        Existing_variation="rs123",
        STRAND="1",
        ALLELE_NUM="1",
    )
    handle, record = _record(tmp_path, [mane, severe_noncanonical])
    try:
        annotation = VepParser.from_header(handle.header).parse(record)
    finally:
        handle.close()

    assert annotation is not None
    assert annotation["rsid"] == "rs123"
    assert annotation["gene"] == {
        "symbol": "GENE1",
        "ensembl_id": "ENSG000001.5",
        "hgnc_id": "HGNC:42",
        "biotype": "protein_coding",
        "is_canonical": True,
    }
    consequence = annotation["consequence"]
    assert consequence["most_severe"] == "stop_gained"
    assert consequence["impact"] == "HIGH"
    assert consequence["hgvsc"] == "ENST_MANE.2:c.100A>G"
    assert consequence["hgvsp"] == "ENSP_MANE.1:p.Lys34Arg"
    assert len(consequence["all_annotations"]) == 2
    assert consequence["all_annotations"][0]["mane_select"] == "NM_000001.3"
    assert consequence["all_annotations"][0]["domains"] == [
        "Pfam:PF00001",
        "PROSITE:PS00001",
    ]
    frequency = annotation["population_frequency"]
    assert frequency["gnomad_af"] == 0.002
    assert frequency["gnomad_af_popmax"] == 0.01
    assert frequency["is_rare"] is False
    predictions = annotation["predictions"]
    assert predictions["cadd_phred"] == 25.4
    assert predictions["spliceai_max"] == 0.42
    assert predictions["sift"] == 0.01
    assert predictions["polyphen"] == 0.99
    assert predictions["consensus_damaging"] is False


def test_annotation_is_arrow_typed_and_missing_scores_cast_no_vote(tmp_path: Path) -> None:
    row = _csq_row(
        Allele="G",
        Consequence="missense_variant",
        IMPACT="MODERATE",
        SYMBOL="GENE2",
        Gene="ENSG000002",
        Feature_type="Transcript",
        Feature="ENST000002",
        BIOTYPE="protein_coding",
        CANONICAL="YES",
        ALLELE_NUM="1",
        CADD_PHRED="22",
    )
    handle, record = _record(tmp_path, [row])
    try:
        annotation = VepParser.from_header(handle.header).parse(record)
    finally:
        handle.close()
    assert annotation is not None
    assert annotation["population_frequency"]["is_rare"] is True
    assert annotation["predictions"]["consensus_damaging"] is False

    variant = _blank_variant_row(
        variant_id="chr1:10:A:G",
        chrom="chr1",
        pos=10,
        ref="A",
        alt="G",
        rsid=None,
        call_confidence="high",
        dp=30,
        gq=50,
        gt=(0, 1),
        phased=False,
        allele_balance=0.5,
    )
    apply_vep_annotation(variant, annotation)
    table = pa.Table.from_pylist([variant], schema=variant_schema())
    assert table.num_rows == 1
    assert table.column("gene")[0].as_py()["symbol"] == "GENE2"


def test_enabled_predictor_fields_are_required_as_one_explicit_gate() -> None:
    parser = VepParser(tuple(CSQ_FIELDS.split("|")))
    parser.require_predictor_fields(("loftee", "cadd", "revel", "alphamissense", "spliceai"))

    without_spliceai = tuple(
        field for field in CSQ_FIELDS.split("|") if not field.startswith("SpliceAI_pred_")
    )
    with pytest.raises(InputValidationError, match="spliceai"):
        VepParser(without_spliceai).require_predictor_fields(("spliceai",))


def test_out_of_range_predictor_score_is_rejected(tmp_path: Path) -> None:
    row = _csq_row(
        Allele="G",
        Consequence="missense_variant",
        IMPACT="MODERATE",
        SYMBOL="GENE",
        Gene="ENSG000001",
        Feature_type="Transcript",
        Feature="ENST000001",
        BIOTYPE="protein_coding",
        ALLELE_NUM="1",
        REVEL="1.01",
    )
    handle, record = _record(tmp_path, [row])
    try:
        with pytest.raises(InputValidationError, match="REVEL score"):
            VepParser.from_header(handle.header).parse(record)
    finally:
        handle.close()
