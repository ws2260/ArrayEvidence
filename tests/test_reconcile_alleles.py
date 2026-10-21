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


def assay(source_seq="CC[A/G]AT"):
    return Assay(
        ilmn_id="synthetic",
        name="synthetic",
        genome_build="38",
        chrom="1",
        map_info="5",
        snp="[A/G]",
        source_seq=source_seq,
    )


def record(ref="G", alt="A", gt="0/1"):
    return VCFRecord(
        chrom="1",
        pos=5,
        variant_id="synthetic",
        ref=ref,
        alts=(alt,),
        gt=gt,
    )


def test_matching_alleles():
    reference = SyntheticReference("AACCGATGCT")

    result = reconcile(
        [assay()],
        [record()],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.coordinate_evidence == "PRESENT"
    assert result.genotype_evidence == "CALLED"
    assert result.allele_evidence == "MATCH"


def test_mismatching_alleles():
    reference = SyntheticReference("AACCGATGCT")

    result = reconcile(
        [assay()],
        [record(alt="C")],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.allele_evidence == "MISMATCH"
    assert result.allele_reason == "ALLELE_SET_DIFFERENT"


def test_no_reference_does_not_claim_match():
    result = reconcile(
        [assay()],
        [record()],
    )[0]

    assert result.allele_evidence == "NOT_EVALUATED"


def test_multiple_candidates_remain_unresolved():
    reference = SyntheticReference("AACCGATGCT")

    result = reconcile(
        [assay()],
        [record(), record(alt="C")],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.coordinate_evidence == "MULTIPLE"
    assert result.allele_evidence == "UNRESOLVED"
    assert result.allele_reason == "MULTIPLE_VCF_CANDIDATES"


def test_no_call_is_separate_from_allele_match():
    reference = SyntheticReference("AACCGATGCT")

    result = reconcile(
        [assay()],
        [record(gt="./.")],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.genotype_evidence == "NO_CALL"
    assert result.allele_evidence == "MATCH"


def test_unvalidated_assay_remains_unresolved():
    reference = SyntheticReference("AACCGATGCT")

    result = reconcile(
        [assay(source_seq="NN[A/G]NN")],
        [record()],
        reference=reference,
        min_informative_bases=3,
    )[0]

    assert result.allele_evidence == "UNRESOLVED"
    assert result.allele_reason == "INSUFFICIENT_INFORMATIVE_BASES"
