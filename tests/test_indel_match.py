import pytest

from arrayevidence.indel_events import SequenceEvent
from arrayevidence.indel_match import match_indel
from arrayevidence.manifest_indel import ManifestIndelResult
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


def manifest_event(start, deleted="", inserted=""):
    status = "INSERTION" if inserted else "DELETION"

    return ManifestIndelResult(
        status="VALIDATED",
        reason="SYNTHETIC_REFERENCE_VALIDATION",
        orientation="FORWARD",
        event=SequenceEvent(
            "1", start, deleted, inserted, status
        ),
    )


def vcf(pos, ref, alt):
    return VCFRecord(
        chrom="1",
        pos=pos,
        variant_id="synthetic",
        ref=ref,
        alts=(alt,),
        gt="0/1",
    )


def test_exact_insertion_match():
    reference = SyntheticReference("ACGT")

    result = match_indel(
        manifest_event(2, inserted="T"),
        vcf(2, "C", "CT"),
        reference,
    )

    assert result.status == "MATCH"


def test_exact_deletion_match():
    reference = SyntheticReference("ACGT")

    result = match_indel(
        manifest_event(2, deleted="G"),
        vcf(2, "CG", "C"),
        reference,
    )

    assert result.status == "MATCH"


def test_repeat_equivalent_insertion():
    reference = SyntheticReference("CAAAAG")

    result = match_indel(
        manifest_event(4, inserted="A"),
        vcf(2, "A", "AA"),
        reference,
    )

    assert result.status == "MATCH"


def test_repeat_equivalent_deletion():
    reference = SyntheticReference("CAAAAG")

    result = match_indel(
        manifest_event(4, deleted="A"),
        vcf(2, "AA", "A"),
        reference,
    )

    assert result.status == "MATCH"


def test_different_insertion():
    reference = SyntheticReference("ACGT")

    result = match_indel(
        manifest_event(2, inserted="A"),
        vcf(2, "C", "CT"),
        reference,
    )

    assert result.status == "MISMATCH"


def test_vcf_reference_mismatch():
    reference = SyntheticReference("ACGT")

    result = match_indel(
        manifest_event(2, inserted="T"),
        vcf(2, "G", "GT"),
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "VCF_REF_DISAGREES_WITH_REFERENCE"


def test_unsupported_complex_vcf():
    reference = SyntheticReference("ACGT")

    result = match_indel(
        manifest_event(2, inserted="T"),
        vcf(2, "CG", "TA"),
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "VCF_NOT_SIMPLE_INDEL"


def test_unvalidated_manifest():
    reference = SyntheticReference("ACGT")

    manifest = ManifestIndelResult(
        "UNRESOLVED",
        "NO_FLANK_ALIGNMENT",
        None,
        None,
    )

    result = match_indel(
        manifest,
        vcf(2, "C", "CT"),
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "MANIFEST_EVENT_NOT_VALIDATED"


def test_multiallelic_vcf_not_supported():
    reference = SyntheticReference("ACGT")

    record = VCFRecord(
        chrom="1",
        pos=2,
        variant_id="synthetic",
        ref="C",
        alts=("CT", "CA"),
        gt="0/1",
    )

    result = match_indel(
        manifest_event(2, inserted="T"),
        record,
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "MULTIALLELIC_VCF_NOT_SUPPORTED"


def test_different_chromosome():
    reference = SyntheticReference("ACGT")

    record = VCFRecord(
        chrom="2",
        pos=2,
        variant_id="synthetic",
        ref="C",
        alts=("CT",),
        gt="0/1",
    )

    result = match_indel(
        manifest_event(2, inserted="T"),
        record,
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "CHROMOSOME_MISMATCH"
