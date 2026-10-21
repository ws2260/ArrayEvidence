"""Relational integrity requirements for technical evidence export."""

import pytest

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AlleleCandidateRow,
    AssayEvidenceRow,
)


def assay(sample="sample_A", name="assay_A", **overrides):
    values = dict(
        schema_version=SCHEMA_VERSION,
        sample_id=sample,
        assay_id=name,
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


def candidate(
    sample="sample_A",
    name="assay_A",
    index=1,
    **overrides,
):
    values = dict(
        schema_version=SCHEMA_VERSION,
        sample_id=sample,
        assay_id=name,
        candidate_index=index,
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


def assert_rejected_without_outputs(tmp_path, assays, candidates):
    with pytest.raises(ValueError):
        export_evidence_tsv((tmp_path / "bundle"), assays, candidates)

    assert not (tmp_path / "bundle" / "assay_evidence.tsv").exists()
    assert not (tmp_path / "bundle" / "allele_candidates.tsv").exists()


def test_duplicate_assay_key_rejected(tmp_path):
    assert_rejected_without_outputs(
        tmp_path,
        [assay(), assay()],
        [],
    )


def test_orphan_candidate_rejected(tmp_path):
    assert_rejected_without_outputs(
        tmp_path,
        [assay()],
        [candidate(name="unknown_assay")],
    )


def test_candidate_from_other_sample_rejected(tmp_path):
    assert_rejected_without_outputs(
        tmp_path,
        [assay()],
        [candidate(sample="sample_B")],
    )


def test_duplicate_candidate_index_rejected(tmp_path):
    assert_rejected_without_outputs(
        tmp_path,
        [assay()],
        [
            candidate(index=1, vcf_id="first"),
            candidate(index=1, vcf_id="second"),
        ],
    )


@pytest.mark.parametrize("indices", [
    [2],
    [1, 3],
    [2, 3],
])
def test_nonconsecutive_candidate_indices_rejected(
    tmp_path, indices
):
    assert_rejected_without_outputs(
        tmp_path,
        [assay()],
        [
            candidate(index=index, vcf_id=f"variant_{index}")
            for index in indices
        ],
    )


def test_two_equivalent_indel_candidates_are_valid(tmp_path):
    ambiguous = assay(
        allele_evidence="UNRESOLVED",
        allele_reason="MULTIPLE_EQUIVALENT_VCF_ALLELES",
        indel_evidence="MULTIPLE",
        projection_status="UNRESOLVED",
        projection_reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )

    export_evidence_tsv(
        (tmp_path / "bundle"),
        [ambiguous],
        [
            candidate(index=2, vcf_id="second", candidate_status="MULTIPLE"),
            candidate(index=1, vcf_id="first", candidate_status="MULTIPLE"),
        ],
    )

    lines = (
        tmp_path / "bundle" / "allele_candidates.tsv"
    ).read_text(encoding="utf-8").splitlines()

    assert len(lines) == 3
    assert "first" in lines[1]
    assert "second" in lines[2]


def test_multiple_samples_can_reuse_same_assay_id(tmp_path):
    export_evidence_tsv(
        (tmp_path / "bundle"),
        [
            assay(sample="sample_A"),
            assay(sample="sample_B"),
        ],
        [
            candidate(sample="sample_A"),
            candidate(sample="sample_B"),
        ],
    )

    lines = (
        tmp_path / "bundle" / "assay_evidence.tsv"
    ).read_text(encoding="utf-8").splitlines()

    assert len(lines) == 3
