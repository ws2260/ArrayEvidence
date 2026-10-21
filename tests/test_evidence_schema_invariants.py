"""Cross-field invariants for versioned technical evidence rows."""

import pytest

from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AlleleCandidateRow,
    AssayEvidenceRow,
)


def make_row(**overrides):
    values = dict(
        schema_version=SCHEMA_VERSION,
        sample_id="synthetic_sample",
        assay_id="synthetic_assay",
        manifest_chrom="1",
        manifest_pos=5,
        coordinate_status="VCF_ABSENT",
        coordinate_evidence="ABSENT",
        coordinate_candidate_count=0,
        coordinate_genotype_evidence="NOT_EVALUATED",
        allele_evidence="MATCH",
        allele_reason="ONE_EQUIVALENT_VCF_ALLELE",
        indel_evidence="UNIQUE",
        projection_status="ASSAY_ALT_PRESENT",
        projection_reason="MATCHED_ALT_PROJECTED",
        assay_alt_dosage=1,
        other_alt_present=False,
        ploidy=2,
    )
    values.update(overrides)
    return AssayEvidenceRow(**values)


@pytest.mark.parametrize(
    "overrides",
    [
        {"assay_alt_dosage": None},
        {"assay_alt_dosage": 0},
        {"other_alt_present": True},
        {"other_alt_present": None},
        {"ploidy": None},
    ],
)
def test_assay_alt_present_requires_evaluable_matching_genotype(overrides):
    with pytest.raises(ValueError):
        make_row(**overrides)


@pytest.mark.parametrize(
    "overrides",
    [
        {"assay_alt_dosage": 1},
        {"assay_alt_dosage": None},
        {"other_alt_present": True},
        {"ploidy": None},
    ],
)
def test_assay_ref_only_requires_zero_dosage(overrides):
    with pytest.raises(ValueError):
        make_row(
            **{
                "projection_status": "ASSAY_REF_ONLY",
                "assay_alt_dosage": 0,
                **overrides,
            }
        )


@pytest.mark.parametrize(
    "overrides",
    [
        {"other_alt_present": False},
        {"other_alt_present": None},
        {"assay_alt_dosage": None},
        {"ploidy": None},
    ],
)
def test_other_alt_present_requires_evaluable_genotype(overrides):
    with pytest.raises(ValueError):
        make_row(
            **{
                "projection_status": "OTHER_ALT_PRESENT",
                "other_alt_present": True,
                **overrides,
            }
        )


@pytest.mark.parametrize(
    "status",
    ["NO_CALL", "NOT_EVALUATED", "UNRESOLVED"],
)
def test_non_evaluable_projection_rejects_dosage(status):
    with pytest.raises(ValueError):
        make_row(projection_status=status)


def test_dosage_cannot_exceed_ploidy():
    with pytest.raises(ValueError):
        make_row(assay_alt_dosage=3, ploidy=2)


def test_valid_other_alt_with_zero_matched_dosage():
    row = make_row(
        projection_status="OTHER_ALT_PRESENT",
        assay_alt_dosage=0,
        other_alt_present=True,
    )
    assert row.assay_alt_dosage == 0
    assert row.other_alt_present is True


def test_valid_no_call_has_unknown_projection_values():
    row = make_row(
        projection_status="NO_CALL",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )
    assert row.assay_alt_dosage is None


def test_unresolved_allele_cannot_have_projected_dosage():
    with pytest.raises(ValueError):
        make_row(allele_evidence="UNRESOLVED")


def test_absent_indel_cannot_have_projected_dosage():
    with pytest.raises(ValueError):
        make_row(indel_evidence="ABSENT")


def test_candidate_source_coordinate_can_differ_from_manifest():
    candidate = AlleleCandidateRow(
        schema_version=SCHEMA_VERSION,
        sample_id="synthetic_sample",
        assay_id="synthetic_assay",
        candidate_index=1,
        vcf_chrom="1",
        vcf_pos=6,
        vcf_id="synthetic_variant",
        vcf_ref="A",
        vcf_alt="AA",
        vcf_alt_index=1,
        vcf_gt="0/1",
        candidate_status="UNIQUE",
    )

    assert candidate.vcf_pos == 6
    assert candidate.vcf_alt_index == 1
