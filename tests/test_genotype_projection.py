import pytest

from arrayevidence.genotype_projection import (
    project_multiallelic_genotype,
)
from arrayevidence.multiallelic import MultiallelicMatch
from arrayevidence.vcf import VCFRecord


def record(gt, alts=("C", "G")):
    return VCFRecord(
        chrom="1",
        pos=100,
        variant_id="synthetic",
        ref="A",
        alts=alts,
        gt=gt,
    )


def match(index=2):
    return MultiallelicMatch(
        status="PARTIAL_MATCH",
        reason="SYNTHETIC",
        matched_alt_indices=(index,),
    )


@pytest.mark.parametrize(
    "gt,status,dosage,other",
    [
        ("0/0", "ASSAY_REF_ONLY", 0, False),
        ("0/2", "ASSAY_ALT_PRESENT", 1, False),
        ("2/2", "ASSAY_ALT_PRESENT", 2, False),
        ("0/1", "OTHER_ALT_PRESENT", 0, True),
        ("1/2", "OTHER_ALT_PRESENT", 1, True),
        ("1/1", "OTHER_ALT_PRESENT", 0, True),
        ("2|0", "ASSAY_ALT_PRESENT", 1, False),
        ("2", "ASSAY_ALT_PRESENT", 1, False),
    ],
)
def test_projection(gt, status, dosage, other):
    result = project_multiallelic_genotype(
        match(), record(gt)
    )

    assert result.status == status
    assert result.assay_alt_dosage == dosage
    assert result.other_alt_present is other


@pytest.mark.parametrize("gt", ["./.", "0/.", "."])
def test_missing_genotypes(gt):
    result = project_multiallelic_genotype(
        match(), record(gt)
    )

    assert result.status == "NO_CALL"
    assert result.assay_alt_dosage is None


@pytest.mark.parametrize("gt", ["0/3", "-1/2", "A/2", "0/2|1"])
def test_invalid_genotypes(gt):
    result = project_multiallelic_genotype(
        match(), record(gt)
    )

    assert result.status == "UNRESOLVED"
    assert result.assay_alt_dosage is None


def test_unvalidated_match_cannot_project():
    invalid_match = MultiallelicMatch(
        status="MISMATCH",
        reason="SYNTHETIC",
        matched_alt_indices=(),
    )

    result = project_multiallelic_genotype(
        invalid_match, record("0/2")
    )

    assert result.status == "UNRESOLVED"


def test_alt_index_one():
    result = project_multiallelic_genotype(
        match(index=1), record("0/1")
    )

    assert result.status == "ASSAY_ALT_PRESENT"
    assert result.assay_alt_dosage == 1


def test_triploid_genotype():
    result = project_multiallelic_genotype(
        match(), record("0/2/2")
    )

    assert result.status == "ASSAY_ALT_PRESENT"
    assert result.assay_alt_dosage == 2
    assert result.ploidy == 3


def test_invalid_target_index():
    result = project_multiallelic_genotype(
        match(index=3), record("0/2")
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "INVALID_MATCHED_ALT_INDEX"


def test_no_implicit_diploid_assumption():
    result = project_multiallelic_genotype(
        match(), record("2")
    )

    assert result.ploidy == 1
    assert result.assay_alt_dosage == 1
