from arrayevidence.assay_alleles import validate_assay_alleles
from arrayevidence.manifest import Assay


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


def assay(source_seq, snp="[A/G]"):
    return Assay(
        ilmn_id="synthetic",
        name="synthetic",
        genome_build="38",
        chrom="1",
        map_info="5",
        snp=snp,
        source_seq=source_seq,
    )


def test_forward_validated_alleles():
    reference = SyntheticReference("AACCGATGCT")

    result = validate_assay_alleles(
        assay("CC[A/G]AT"),
        reference,
        min_informative_bases=3,
    )

    assert result.status == "VALIDATED"
    assert result.orientation == "FORWARD"
    assert result.alleles == ("A", "G")
    assert result.reference_base == "G"


def test_reverse_validated_alleles():
    reference = SyntheticReference("AACCGATGCT")

    result = validate_assay_alleles(
        assay("AT[C/T]GG", snp="[C/T]"),
        reference,
        min_informative_bases=3,
    )

    assert result.status == "VALIDATED"
    assert result.orientation == "REVERSE"
    assert result.alleles == ("A", "G")
    assert result.reference_base == "G"


def test_unresolved_when_reference_allele_absent():
    reference = SyntheticReference("AACCGATGCT")

    result = validate_assay_alleles(
        assay("CC[A/C]AT", snp="[A/C]"),
        reference,
        min_informative_bases=3,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "REFERENCE_BASE_NOT_IN_ASSAY_ALLELES"


def test_unresolved_when_orientation_not_established():
    reference = SyntheticReference("AACCGATGCT")

    result = validate_assay_alleles(
        assay("NN[A/G]NN"),
        reference,
        min_informative_bases=3,
    )

    assert result.status == "UNRESOLVED"
