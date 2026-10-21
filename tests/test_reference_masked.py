from arrayevidence.manifest import Assay
from arrayevidence.reference import validate_orientation_masked


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


def assay(source_seq):
    return Assay(
        ilmn_id="synthetic",
        name="synthetic",
        genome_build="38",
        chrom="1",
        map_info="5",
        snp="[A/G]",
        source_seq=source_seq,
    )


def test_masked_forward():
    ref = SyntheticReference("AACCGATGCT")

    result = validate_orientation_masked(
        assay("CN[A/G]AT"),
        ref,
        min_informative_bases=3,
    )

    assert result.status == "FORWARD"
    assert result.reason == "MASKED_FORWARD_MATCH"


def test_insufficient_informative_bases():
    ref = SyntheticReference("AACCGATGCT")

    result = validate_orientation_masked(
        assay("NN[A/G]NN"),
        ref,
        min_informative_bases=3,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "INSUFFICIENT_INFORMATIVE_BASES"


def test_masked_reverse():
    ref = SyntheticReference("AACCGATGCT")

    result = validate_orientation_masked(
        assay("AN[A/G]GG"),
        ref,
        min_informative_bases=3,
    )

    assert result.status == "REVERSE"


def test_masked_mismatch():
    ref = SyntheticReference("AACCGATGCT")

    result = validate_orientation_masked(
        assay("AA[A/G]AA"),
        ref,
        min_informative_bases=3,
    )

    assert result.status == "NO_MATCH"
