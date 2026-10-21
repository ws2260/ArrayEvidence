"""Synthetic contract tests for reconciliation-to-evidence conversion."""

import pytest

from arrayevidence.evidence_builder import build_evidence_rows
from arrayevidence.reconcile import AlleleSource, ReconciliationResult
from arrayevidence.vcf import VCFRecord


def record(
    *,
    pos=6,
    variant_id="synthetic_variant",
    ref="A",
    alts=("AA",),
    gt="0/1",
):
    return VCFRecord(
        chrom="1",
        pos=pos,
        variant_id=variant_id,
        ref=ref,
        alts=tuple(alts),
        gt=gt,
    )


def result(**overrides):
    values = dict(
        assay_id="synthetic_assay",
        coordinate=("1", 5),
        status="VCF_ABSENT",
        candidate_record_count=0,
        coordinate_evidence="ABSENT",
        genotype_evidence="NOT_EVALUATED",
        allele_evidence="MATCH",
        allele_reason="ONE_EQUIVALENT_VCF_ALLELE",
        projection_status="ASSAY_ALT_PRESENT",
        projection_reason="MATCHED_ALT_PROJECTED",
        assay_alt_dosage=1,
        other_alt_present=False,
        ploidy=2,
        matched_alt_indices=(1,),
        indel_evidence="UNIQUE",
        allele_sources=(AlleleSource(record(), 1),),
    )
    values.update(overrides)
    return ReconciliationResult(**values)


def test_shifted_indel_preserves_coordinate_and_source_separately():
    assays, candidates = build_evidence_rows(
        "synthetic_sample",
        [result()],
    )

    assert len(assays) == 1
    assert len(candidates) == 1

    assay = assays[0]
    candidate = candidates[0]

    assert assay.sample_id == "synthetic_sample"
    assert assay.manifest_chrom == "1"
    assert assay.manifest_pos == 5
    assert assay.coordinate_status == "VCF_ABSENT"
    assert assay.coordinate_evidence == "ABSENT"
    assert assay.coordinate_candidate_count == 0
    assert assay.indel_evidence == "UNIQUE"
    assert assay.projection_status == "ASSAY_ALT_PRESENT"
    assert assay.assay_alt_dosage == 1

    assert candidate.vcf_pos == 6
    assert candidate.vcf_ref == "A"
    assert candidate.vcf_alt == "AA"
    assert candidate.vcf_alt_index == 1
    assert candidate.vcf_gt == "0/1"
    assert candidate.candidate_status == "UNIQUE"


def test_ambiguous_indels_preserve_all_candidates_without_projection():
    first = record(pos=4, variant_id="first", ref="T", alts=("TTT",))
    second = record(pos=5, variant_id="second", ref="G", alts=("GTT",))

    unresolved = result(
        allele_evidence="UNRESOLVED",
        allele_reason="MULTIPLE_EQUIVALENT_VCF_ALLELES",
        indel_evidence="MULTIPLE",
        projection_status="UNRESOLVED",
        projection_reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
        matched_alt_indices=(1, 1),
        allele_sources=(
            AlleleSource(first, 1),
            AlleleSource(second, 1),
        ),
    )

    assays, candidates = build_evidence_rows(
        "synthetic_sample",
        [unresolved],
    )

    assert len(assays) == 1
    assert assays[0].assay_alt_dosage is None
    assert assays[0].projection_status == "UNRESOLVED"
    assert [c.vcf_id for c in candidates] == ["first", "second"]
    assert [c.candidate_index for c in candidates] == [1, 2]
    assert all(c.candidate_status == "MULTIPLE" for c in candidates)


def test_absent_indel_does_not_create_candidate_or_hom_ref():
    absent = result(
        allele_evidence="UNRESOLVED",
        allele_reason="NO_EQUIVALENT_VCF_ALLELE",
        indel_evidence="ABSENT",
        projection_status="NOT_EVALUATED",
        projection_reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
        matched_alt_indices=(),
        allele_sources=(),
    )

    assays, candidates = build_evidence_rows(
        "synthetic_sample",
        [absent],
    )

    assert len(assays) == 1
    assert candidates == []
    assert assays[0].assay_alt_dosage is None
    assert assays[0].coordinate_genotype_evidence == "NOT_EVALUATED"


def test_multiallelic_snv_preserves_second_alt_and_phased_gt():
    source = record(
        pos=5,
        ref="G",
        alts=("C", "A"),
        gt="0|2",
    )

    snv = result(
        status="VCF_CALLED",
        candidate_record_count=1,
        coordinate_evidence="PRESENT",
        genotype_evidence="CALLED",
        allele_evidence="PARTIAL_MATCH",
        allele_reason="ASSAY_ALT_PRESENT_IN_MULTIALLELIC_VCF",
        indel_evidence="NOT_EVALUATED",
        matched_alt_indices=(2,),
        allele_sources=(AlleleSource(source, 2),),
    )

    assays, candidates = build_evidence_rows(
        "synthetic_sample",
        [snv],
    )

    assert assays[0].coordinate_candidate_count == 1
    assert candidates[0].vcf_alt == "A"
    assert candidates[0].vcf_alt_index == 2
    assert candidates[0].vcf_gt == "0|2"
    assert candidates[0].candidate_status == "PARTIAL_MATCH"


def test_unmapped_assay_preserves_null_coordinates():
    unmapped = result(
        coordinate=None,
        status="UNMAPPED_ASSAY",
        coordinate_evidence="UNMAPPED",
        genotype_evidence="NOT_EVALUATED",
        allele_evidence="NOT_EVALUATED",
        allele_reason="UNMAPPED_ASSAY",
        indel_evidence="NOT_EVALUATED",
        projection_status="NOT_EVALUATED",
        projection_reason="NOT_EVALUATED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
        matched_alt_indices=(),
        allele_sources=(),
    )

    assays, candidates = build_evidence_rows(
        "synthetic_sample",
        [unmapped],
    )

    assert assays[0].manifest_chrom is None
    assert assays[0].manifest_pos is None
    assert candidates == []


def test_multiple_assays_preserve_input_order():
    first = result(assay_id="assay_B")
    second = result(
        assay_id="assay_A",
        allele_sources=(AlleleSource(record(variant_id="other"), 1),),
    )

    assays, candidates = build_evidence_rows(
        "synthetic_sample",
        [first, second],
    )

    assert [row.assay_id for row in assays] == ["assay_B", "assay_A"]
    assert [row.assay_id for row in candidates] == ["assay_B", "assay_A"]


def test_empty_sample_id_is_rejected():
    with pytest.raises(ValueError):
        build_evidence_rows("", [result()])


def test_invalid_source_alt_index_is_rejected():
    invalid = result(
        allele_sources=(AlleleSource(record(alts=("AA",)), 2),),
    )

    with pytest.raises(ValueError):
        build_evidence_rows("synthetic_sample", [invalid])


def test_empty_results_produce_empty_tables():
    assays, candidates = build_evidence_rows(
        "synthetic_sample",
        [],
    )

    assert assays == []
    assert candidates == []
