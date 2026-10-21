"""Synthetic contract tests for versioned technical evidence."""

import pytest

from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AlleleCandidateRow,
    AssayEvidenceRow,
)


def make_assay_row(**overrides):
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


def make_candidate_row(**overrides):
    values = dict(
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
    values.update(overrides)
    return AlleleCandidateRow(**values)


def test_schema_version_is_explicit():
    assert SCHEMA_VERSION == "1.0.0"


def test_assay_row_preserves_shifted_indel_evidence():
    row = make_assay_row()

    assert row.coordinate_status == "VCF_ABSENT"
    assert row.coordinate_evidence == "ABSENT"
    assert row.indel_evidence == "UNIQUE"
    assert row.allele_evidence == "MATCH"
    assert row.assay_alt_dosage == 1


def test_missing_evidence_does_not_imply_hom_ref():
    row = make_assay_row(
        allele_evidence="UNRESOLVED",
        indel_evidence="ABSENT",
        projection_status="NOT_EVALUATED",
        projection_reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )

    assert row.assay_alt_dosage is None
    assert row.projection_status == "NOT_EVALUATED"


def test_candidate_preserves_shifted_vcf_source():
    candidate = make_candidate_row()

    assert candidate.vcf_pos == 6
    assert candidate.vcf_alt == "AA"
    assert candidate.vcf_alt_index == 1


@pytest.mark.parametrize(
    "invalid_index",
    [0, -1],
)
def test_candidate_rejects_invalid_alt_index(invalid_index):
    with pytest.raises(ValueError):
        make_candidate_row(vcf_alt_index=invalid_index)


def test_assay_row_rejects_negative_dosage():
    with pytest.raises(ValueError):
        make_assay_row(assay_alt_dosage=-1)


def test_assay_row_requires_identifiers():
    with pytest.raises(ValueError):
        make_assay_row(sample_id="")


def test_candidate_requires_identifiers():
    with pytest.raises(ValueError):
        make_candidate_row(assay_id="")


def test_assay_row_rejects_unsupported_schema_version():
    with pytest.raises(ValueError):
        make_assay_row(schema_version="999.0.0")
