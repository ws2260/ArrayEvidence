"""Synthetic end-to-end reconciliation tests for multiallelic SNVs."""

import pytest

from arrayevidence.manifest import Assay
from arrayevidence.reconcile import reconcile
from arrayevidence.vcf import VCFRecord


class SyntheticReference:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, reference, start, end):
        if reference != "1":
            raise KeyError(reference)
        return self.sequence[start:end]

    def get_reference_length(self, reference):
        if reference != "1":
            raise KeyError(reference)
        return len(self.sequence)


def make_assay():
    return Assay(
        ilmn_id="synthetic_assay",
        name="synthetic",
        genome_build="38",
        chrom="1",
        map_info="5",
        snp="[A/G]",
        source_seq="CC[A/G]AT",
    )


def make_record(gt, alts=("C", "G")):
    return VCFRecord(
        chrom="1",
        pos=5,
        variant_id="synthetic",
        ref="G",
        alts=alts,
        gt=gt,
    )


@pytest.mark.parametrize(
    "gt,status,dosage,other",
    [
        ("0/0", "ASSAY_REF_ONLY", 0, False),
        ("0/2", "ASSAY_ALT_PRESENT", 1, False),
        ("1/2", "OTHER_ALT_PRESENT", 1, True),
        ("0/1", "OTHER_ALT_PRESENT", 0, True),
        ("./.", "NO_CALL", None, None),
    ],
)
def test_integrated_projection(gt, status, dosage, other):
    reference = SyntheticReference("AACCGATGCT")

    # Reference is G, so assay ALT is A.
    # VCF ALT2 is A; ALT1 is C.
    record = make_record(gt, alts=("C", "A"))

    result = reconcile(
        [make_assay()],
        [record],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.allele_evidence == "PARTIAL_MATCH"
    assert result.matched_alt_indices == (2,)
    assert result.projection_status == status
    assert result.assay_alt_dosage == dosage
    assert result.other_alt_present is other


def test_multiallelic_reference_mismatch():
    reference = SyntheticReference("AACCGATGCT")

    record = VCFRecord(
        chrom="1",
        pos=5,
        variant_id="synthetic",
        ref="T",
        alts=("C", "A"),
        gt="0/2",
    )

    result = reconcile(
        [make_assay()],
        [record],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.allele_evidence == "MISMATCH"
    assert result.projection_status == "NOT_EVALUATED"


def test_multiple_records_not_silently_resolved():
    reference = SyntheticReference("AACCGATGCT")

    result = reconcile(
        [make_assay()],
        [
            make_record("0/2", alts=("C", "A")),
            make_record("0/1", alts=("A",)),
        ],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.candidate_record_count == 2
    assert result.allele_evidence == "UNRESOLVED"
    assert result.allele_reason == "MULTIPLE_VCF_CANDIDATES"
    assert result.projection_status == "NOT_EVALUATED"


def test_biallelic_behavior_preserved():
    reference = SyntheticReference("AACCGATGCT")

    result = reconcile(
        [make_assay()],
        [make_record("0/1", alts=("A",))],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.allele_evidence == "MATCH"
    assert result.projection_status == "NOT_EVALUATED"
    assert result.matched_alt_indices == ()


def test_without_reference_does_not_claim_partial_match():
    result = reconcile(
        [make_assay()],
        [make_record("0/2", alts=("C", "A"))],
    )[0]

    assert result.allele_evidence == "NOT_EVALUATED"
    assert result.projection_status == "NOT_EVALUATED"
