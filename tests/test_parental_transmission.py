"""The recorded trio can resolve transmission even when both parents carry ALT."""

from dataclasses import replace
from itertools import product

from genome_annotation_builder.family_evidence import (
    FamilyVariant,
    PedigreeMember,
    _pedigree_status,
)

CHILD = PedigreeMember("family", "child", "father", "mother", 0)
SITE = FamilyVariant(
    variant_id="chr1:10:C:T",
    chrom="chr1",
    pos=10,
    ref="C",
    alt="T",
    subject_gt=(0, 1),
    subject_phased=False,
    father_gt=(1, 1),
    mother_gt=(0, 1),
    phase_set=None,
    family_confidence=None,
    family_reasons=(),
    caller_support=None,
    replicate_support=None,
    de_novo_score=None,
    difficult_region_flags=(),
)


def test_complete_trio_transmission_table() -> None:
    paternal = "inherited_paternal"
    maternal = "inherited_maternal"
    both = "inherited_biparental"
    ambiguous = "inherited_parent_ambiguous"
    conflict = "mendelian_conflict"
    no_alt = "not_assessed"
    # Father and mother ALT counts index the rows. The columns are 0/0, 0/1 and 1/1 children.
    expected = {
        (0, 0): (no_alt, conflict, conflict),
        (0, 1): (no_alt, maternal, conflict),
        (0, 2): (conflict, maternal, conflict),
        (1, 0): (no_alt, paternal, conflict),
        (1, 1): (no_alt, ambiguous, both),
        (1, 2): (conflict, maternal, both),
        (2, 0): (conflict, paternal, conflict),
        (2, 1): (conflict, paternal, both),
        (2, 2): (conflict, conflict, both),
    }
    carriers = {
        (False, False): "neither_parent_observed",
        (True, False): "paternal_observed",
        (False, True): "maternal_observed",
        (True, True): "both_parents_observed",
    }
    genotypes = list(product((0, 1), repeat=2))
    for father, mother, child in product(genotypes, repeat=3):
        state = expected[sum(father), sum(mother)][sum(child)]
        for phased in (False, True):
            record = replace(
                SITE, father_gt=father, mother_gt=mother, subject_gt=child, subject_phased=phased
            )
            assert _pedigree_status(record, CHILD) == (
                "conflict" if state == conflict else "consistent",
                carriers[1 in father, 1 in mother],
                state,
            )


def test_transmission_requires_complete_autosomal_diploid_trio() -> None:
    unavailable = ("not_assessed", "unknown", "not_assessed")
    unresolved = ("parent_or_subject_no_call", "unknown", "unresolved_parent_no_call")
    for chrom in ("chrX", "chrY", "chrM"):
        assert _pedigree_status(replace(SITE, chrom=chrom), CHILD) == unavailable
    assert _pedigree_status(replace(SITE, subject_gt=(1,)), CHILD) == unavailable
    assert _pedigree_status(SITE, replace(CHILD, father=None)) == unavailable
    assert _pedigree_status(SITE, replace(CHILD, father=None, mother=None)) == (
        "not_applicable_founder",
        "not_applicable",
        "founder_observation",
    )
    for parental in (None, (1,), (1, None), (None, 1), (1, 2)):
        assert _pedigree_status(replace(SITE, father_gt=parental), CHILD) == unresolved
        assert _pedigree_status(replace(SITE, mother_gt=parental), CHILD) == unresolved
