from arrayevidence.manifest import Assay
from arrayevidence.reference import (
    reverse_complement,
    validate_orientation,
)


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


def assay(source_seq, pos=5):
    return Assay(
        ilmn_id="synthetic",
        name="synthetic",
        genome_build="38",
        chrom="1",
        map_info=str(pos),
        snp="[A/G]",
        source_seq=source_seq,
    )


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACG") == "CGTT"


def test_forward_orientation():
    ref = SyntheticReference("AACCGATGCT")
    result = validate_orientation(
        assay("CC[A/G]AT", pos=5),
        ref,
    )
    assert result.status == "FORWARD"


def test_reverse_orientation():
    ref = SyntheticReference("AACCGATGCT")
    result = validate_orientation(
        assay("AT[A/G]GG", pos=5),
        ref,
    )
    assert result.status == "REVERSE"


def test_no_match():
    ref = SyntheticReference("AACCGATGCT")
    result = validate_orientation(
        assay("AA[A/G]AA", pos=5),
        ref,
    )
    assert result.status == "NO_MATCH"


def test_unsupported_source_sequence():
    ref = SyntheticReference("AACCGATGCT")
    result = validate_orientation(
        assay("CCAAT"),
        ref,
    )
    assert result.status == "UNRESOLVED"
