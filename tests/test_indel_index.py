from arrayevidence.indel_events import SequenceEvent
from arrayevidence.indel_index import IndelCandidateIndex
from arrayevidence.vcf import VCFRecord


class SyntheticReference:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, reference, start, end):
        if reference != "1":
            raise KeyError(reference)
        if start < 0 or end > len(self.sequence):
            raise ValueError("Out of bounds")
        return self.sequence[start:end]

    def get_reference_length(self, reference):
        if reference != "1":
            raise KeyError(reference)
        return len(self.sequence)


def record(pos, ref, alts, gt="0/1"):
    return VCFRecord(
        chrom="1",
        pos=pos,
        variant_id="synthetic",
        ref=ref,
        alts=tuple(alts),
        gt=gt,
    )


def test_unique_insertion():
    reference = SyntheticReference("ACGT")
    candidate = record(2, "C", ["CT"])
    index = IndelCandidateIndex([candidate], reference)

    result = index.find(
        SequenceEvent("1", 2, "", "T", "INSERTION")
    )

    assert result.status == "UNIQUE"
    assert result.records == (candidate,)


def test_equivalent_repeat_representations_preserved():
    reference = SyntheticReference("CAAAAG")

    first = record(2, "AA", ["A"], "0/0")
    second = record(3, "AA", ["A"], "0/0")

    index = IndelCandidateIndex(
        [first, second],
        reference,
    )

    result = index.find(
        SequenceEvent("1", 4, "A", "", "DELETION")
    )

    assert result.status == "MULTIPLE"
    assert result.records == (first, second)


def test_no_equivalent_allele():
    reference = SyntheticReference("ACGT")

    index = IndelCandidateIndex(
        [record(2, "C", ["CT"])],
        reference,
    )

    result = index.find(
        SequenceEvent("1", 2, "", "A", "INSERTION")
    )

    assert result.status == "ABSENT"


def test_invalid_vcf_reference_excluded():
    reference = SyntheticReference("ACGT")

    index = IndelCandidateIndex(
        [record(2, "G", ["GT"])],
        reference,
    )

    result = index.find(
        SequenceEvent("1", 2, "", "T", "INSERTION")
    )

    assert result.status == "ABSENT"
    assert index.skipped_alleles == 1


def test_multiallelic_record_can_be_indexed():
    reference = SyntheticReference("ACGT")

    candidate = record(
        2, "C", ["CT", "CA"]
    )

    index = IndelCandidateIndex(
        [candidate],
        reference,
    )

    result = index.find(
        SequenceEvent("1", 2, "", "A", "INSERTION")
    )

    assert result.status == "UNIQUE"
    assert result.records == (candidate,)


def test_unsupported_manifest_event():
    reference = SyntheticReference("ACGT")
    index = IndelCandidateIndex([], reference)

    result = index.find(
        SequenceEvent("1", 2, "C", "T", "SUBSTITUTION")
    )

    assert result.status == "UNRESOLVED"


def test_duplicate_records_not_discarded():
    reference = SyntheticReference("ACGT")

    first = record(2, "C", ["CT"], "0/0")
    second = record(2, "C", ["CT"], "0/1")

    index = IndelCandidateIndex(
        [first, second],
        reference,
    )

    result = index.find(
        SequenceEvent("1", 2, "", "T", "INSERTION")
    )

    assert result.status == "MULTIPLE"
    assert len(result.records) == 2
    assert result.records[0].gt == "0/0"
    assert result.records[1].gt == "0/1"


def test_multiallelic_candidate_preserves_alt_index():
    """The matched ALT must be identified independently of the record."""
    reference = SyntheticReference("ACGT")
    candidate = record(
        2, "C", ["CT", "CA"], gt="0/2"
    )

    index = IndelCandidateIndex([candidate], reference)

    result = index.find(
        SequenceEvent("1", 2, "", "A", "INSERTION")
    )

    assert result.status == "UNIQUE"
    assert result.records == (candidate,)

    # VCF ALT indices are 1-based; 0 denotes REF.
    assert len(result.alleles) == 1
    assert result.alleles[0].record == candidate
    assert result.alleles[0].alt_index == 2


def test_equivalent_alts_in_one_record_remain_distinct():
    """Two matching ALTs are two allele candidates, not one."""
    reference = SyntheticReference("ACGT")
    candidate = record(
        2, "C", ["CT", "CT"], gt="1/2"
    )

    index = IndelCandidateIndex([candidate], reference)

    result = index.find(
        SequenceEvent("1", 2, "", "T", "INSERTION")
    )

    assert result.status == "MULTIPLE"

    assert result.records == (candidate, candidate)
    assert tuple(
        allele.alt_index for allele in result.alleles
    ) == (1, 2)
