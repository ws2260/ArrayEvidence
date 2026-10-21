"""Synthetic tests for allele-aware indel genotype projection."""

import pytest

from arrayevidence.indel_events import SequenceEvent
from arrayevidence.indel_index import IndelCandidateIndex
from arrayevidence.indel_projection import project_indel_candidate
from arrayevidence.vcf import VCFRecord


class SyntheticReference:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, chrom, start, end):
        if chrom != "1":
            raise KeyError(chrom)
        if start < 0 or end > len(self.sequence):
            raise ValueError("Out of bounds")
        return self.sequence[start:end]

    def get_reference_length(self, chrom):
        if chrom != "1":
            raise KeyError(chrom)
        return len(self.sequence)


def make_record(gt, alts=("CT", "CA")):
    return VCFRecord(
        chrom="1",
        pos=2,
        variant_id="synthetic",
        ref="C",
        alts=alts,
        gt=gt,
    )


def find_insertion(records, inserted="A"):
    index = IndelCandidateIndex(
        records,
        SyntheticReference("ACGT"),
    )
    return index.find(
        SequenceEvent("1", 2, "", inserted, "INSERTION")
    )


@pytest.mark.parametrize(
    "gt,status,dosage,other,ploidy",
    [
        ("0/2", "ASSAY_ALT_PRESENT", 1, False, 2),
        ("1/2", "OTHER_ALT_PRESENT", 1, True, 2),
        ("0/1", "OTHER_ALT_PRESENT", 0, True, 2),
        ("0/0", "ASSAY_REF_ONLY", 0, False, 2),
        ("2/2", "ASSAY_ALT_PRESENT", 2, False, 2),
        ("./.", "NO_CALL", None, None, None),
        ("2", "ASSAY_ALT_PRESENT", 1, False, 1),
    ],
)
def test_unique_indel_projection(gt, status, dosage, other, ploidy):
    result = find_insertion([make_record(gt)])

    assert result.status == "UNIQUE"
    assert result.alleles[0].alt_index == 2

    projection = project_indel_candidate(result)

    assert projection.status == status
    assert projection.assay_alt_dosage == dosage
    assert projection.other_alt_present is other
    assert projection.ploidy == ploidy


def test_ambiguous_records_cannot_be_projected():
    result = find_insertion([
        make_record("0/2"),
        make_record("0/1", alts=("CA",)),
    ])

    assert result.status == "MULTIPLE"

    projection = project_indel_candidate(result)

    assert projection.status == "UNRESOLVED"
    assert projection.assay_alt_dosage is None
    assert projection.other_alt_present is None


def test_duplicate_matching_alts_cannot_be_projected():
    result = find_insertion([
        make_record("1/2", alts=("CA", "CA")),
    ])

    assert result.status == "MULTIPLE"

    projection = project_indel_candidate(result)

    assert projection.status == "UNRESOLVED"
    assert projection.assay_alt_dosage is None


def test_absent_indel_cannot_be_projected():
    result = find_insertion([make_record("0/1")], inserted="G")

    assert result.status == "ABSENT"

    projection = project_indel_candidate(result)

    assert projection.status == "UNRESOLVED"
    assert projection.assay_alt_dosage is None
