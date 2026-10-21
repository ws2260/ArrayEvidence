"""Synthetic parent-child evidence consistency contract."""

import pytest

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AlleleCandidateRow,
    AssayEvidenceRow,
)


def assay(**overrides):
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
        allele_reason="REFERENCE_NORMALIZED_INDEL_MATCH",
        indel_evidence="UNIQUE",
        projection_status="ASSAY_ALT_PRESENT",
        projection_reason="MATCHED_ALT_PROJECTED",
        assay_alt_dosage=1,
        other_alt_present=False,
        ploidy=2,
    )
    values.update(overrides)
    return AssayEvidenceRow(**values)


def candidate(index=1, **overrides):
    values = dict(
        schema_version=SCHEMA_VERSION,
        sample_id="synthetic_sample",
        assay_id="synthetic_assay",
        candidate_index=index,
        vcf_chrom="1",
        vcf_pos=6,
        vcf_id=f"synthetic_variant_{index}",
        vcf_ref="A",
        vcf_alt="AA",
        vcf_alt_index=1,
        vcf_gt="0/1",
        candidate_status="UNIQUE",
    )
    values.update(overrides)
    return AlleleCandidateRow(**values)


def assert_invalid(tmp_path, parent, children):
    with pytest.raises(ValueError):
        export_evidence_tsv((tmp_path / "bundle"), [parent], children)

    assert not (tmp_path / "bundle" / "assay_evidence.tsv").exists()
    assert not (tmp_path / "bundle" / "allele_candidates.tsv").exists()


def unresolved_indel(status):
    return assay(
        allele_evidence="UNRESOLVED",
        allele_reason="NO_UNIQUE_REFERENCE_MATCH",
        indel_evidence=status,
        projection_status="NOT_EVALUATED",
        projection_reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )


def ambiguous_indel():
    return assay(
        allele_evidence="UNRESOLVED",
        allele_reason="MULTIPLE_EQUIVALENT_VCF_ALLELES",
        indel_evidence="MULTIPLE",
        projection_status="UNRESOLVED",
        projection_reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )


def snv_parent(status):
    return assay(
        coordinate_status="VCF_CALLED",
        coordinate_evidence="PRESENT",
        coordinate_candidate_count=1,
        coordinate_genotype_evidence="CALLED",
        allele_evidence=status,
        allele_reason="SYNTHETIC_SNV_EVIDENCE",
        indel_evidence="NOT_EVALUATED",
        projection_status="NOT_EVALUATED",
        projection_reason="NOT_EVALUATED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )


def test_unique_indel_requires_one_candidate(tmp_path):
    assert_invalid(tmp_path, assay(), [])


def test_unique_indel_rejects_wrong_candidate_status(tmp_path):
    assert_invalid(
        tmp_path,
        assay(),
        [candidate(candidate_status="MULTIPLE")],
    )


def test_multiple_indel_requires_at_least_two_candidates(tmp_path):
    assert_invalid(
        tmp_path,
        ambiguous_indel(),
        [candidate(candidate_status="MULTIPLE")],
    )


def test_multiple_indel_rejects_mixed_candidate_status(tmp_path):
    assert_invalid(
        tmp_path,
        ambiguous_indel(),
        [
            candidate(1, candidate_status="MULTIPLE"),
            candidate(2, candidate_status="UNIQUE"),
        ],
    )


def test_multiple_indel_rejects_resolved_projection(tmp_path):
    # This state is rejected at the schema boundary, before export.
    with pytest.raises(
        ValueError,
        match="Evaluable projection requires established allele evidence",
    ):
        assay(
            allele_evidence="UNRESOLVED",
            indel_evidence="MULTIPLE",
            projection_status="ASSAY_ALT_PRESENT",
            assay_alt_dosage=1,
            other_alt_present=False,
            ploidy=2,
        )

    assert not (tmp_path / "bundle" / "assay_evidence.tsv").exists()
    assert not (tmp_path / "bundle" / "allele_candidates.tsv").exists()


@pytest.mark.parametrize("status", ["ABSENT", "UNRESOLVED"])
def test_nonmatched_indel_rejects_candidates(tmp_path, status):
    assert_invalid(
        tmp_path,
        unresolved_indel(status),
        [candidate()],
    )


def test_snv_match_requires_one_candidate(tmp_path):
    assert_invalid(tmp_path, snv_parent("MATCH"), [])


def test_snv_partial_match_requires_one_candidate(tmp_path):
    assert_invalid(
        tmp_path,
        snv_parent("PARTIAL_MATCH"),
        [],
    )


def test_snv_match_rejects_wrong_candidate_status(tmp_path):
    assert_invalid(
        tmp_path,
        snv_parent("MATCH"),
        [candidate(candidate_status="PARTIAL_MATCH")],
    )


@pytest.mark.parametrize(
    "status",
    ["MISMATCH", "UNRESOLVED", "NOT_EVALUATED"],
)
def test_nonmatched_snv_rejects_candidates(tmp_path, status):
    assert_invalid(
        tmp_path,
        snv_parent(status),
        [candidate(candidate_status="MATCH")],
    )


def test_valid_unique_indel(tmp_path):
    export_evidence_tsv(
        (tmp_path / "bundle"),
        [assay()],
        [candidate()],
    )
    assert (tmp_path / "bundle" / "allele_candidates.tsv").is_file()


def test_valid_multiple_indel(tmp_path):
    export_evidence_tsv(
        (tmp_path / "bundle"),
        [ambiguous_indel()],
        [
            candidate(1, candidate_status="MULTIPLE"),
            candidate(2, candidate_status="MULTIPLE"),
        ],
    )
    assert (tmp_path / "bundle" / "allele_candidates.tsv").is_file()


@pytest.mark.parametrize("status", ["MATCH", "PARTIAL_MATCH"])
def test_valid_snv_match(tmp_path, status):
    export_evidence_tsv(
        (tmp_path / "bundle"),
        [snv_parent(status)],
        [candidate(candidate_status=status)],
    )
    assert (tmp_path / "bundle" / "allele_candidates.tsv").is_file()


def test_valid_absent_indel(tmp_path):
    export_evidence_tsv(
        (tmp_path / "bundle"),
        [unresolved_indel("ABSENT")],
        [],
    )
    assert (tmp_path / "bundle" / "assay_evidence.tsv").is_file()
