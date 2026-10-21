"""Synthetic regression tests inspired by complex multi-assay loci.

No participant data or proprietary manifest content is included.
"""

from collections import Counter

from arrayevidence.assay_alleles import ValidatedAssayAlleles
from arrayevidence.validated_match import match_validated_alleles
from arrayevidence.vcf import VCFRecord


def make_record(ref, alt, gt="0/0"):
    return VCFRecord(
        chrom="1",
        pos=100,
        variant_id="synthetic",
        ref=ref,
        alts=(alt,),
        gt=gt,
    )


def make_validated_alleles(ref, alt):
    return ValidatedAssayAlleles(
        status="VALIDATED",
        reason="REFERENCE_VALIDATED",
        orientation="FORWARD",
        alleles=tuple(sorted((ref, alt))),
        reference_base=ref,
    )


def test_same_coordinate_and_reference_do_not_imply_match():
    assay = make_validated_alleles("T", "G")
    record = make_record("T", "C")

    result = match_validated_alleles(assay, record)

    assert result.status == "MISMATCH"
    assert result.reason == "ALLELE_SET_DIFFERENT"


def test_homozygous_reference_does_not_override_allele_mismatch():
    assay = make_validated_alleles("T", "G")
    record = make_record("T", "C", gt="0/0")

    assert record.call_status == "HOM_REF"

    result = match_validated_alleles(assay, record)

    assert result.status == "MISMATCH"


def test_repeated_assays_are_not_independent_coordinates():
    assays = [make_validated_alleles("T", "G") for _ in range(10)]
    record = make_record("T", "C")

    results = [
        match_validated_alleles(assay, record)
        for assay in assays
    ]

    counts = Counter(result.status for result in results)

    assert counts["MISMATCH"] == 10

    # Ten assay-level mismatches can arise at one genomic position.
    assay_coordinates = [("1", 100)] * len(assays)

    assert len(set(assay_coordinates)) == 1


def test_distinct_alt_alleles_are_not_interchangeable():
    assay_g = make_validated_alleles("T", "G")
    assay_c = make_validated_alleles("T", "C")
    record = make_record("T", "C")

    assert match_validated_alleles(
        assay_g, record
    ).status == "MISMATCH"

    assert match_validated_alleles(
        assay_c, record
    ).status == "MATCH"


def test_reference_disagreement_is_not_an_allele_match():
    assay = make_validated_alleles("T", "G")
    record = make_record("C", "G")

    result = match_validated_alleles(assay, record)

    assert result.status == "MISMATCH"
    assert result.reason == "VCF_REF_DISAGREES_WITH_REFERENCE"
