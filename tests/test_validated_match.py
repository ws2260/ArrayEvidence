from arrayevidence.assay_alleles import ValidatedAssayAlleles
from arrayevidence.validated_match import match_validated_alleles
from arrayevidence.vcf import VCFRecord


def validated(alleles=("A", "G"), reference_base="G"):
    return ValidatedAssayAlleles(
        status="VALIDATED",
        reason="REFERENCE_AND_FLANKS_CONCORDANT",
        orientation="FORWARD",
        alleles=alleles,
        reference_base=reference_base,
    )


def record(ref="G", alts=("A",)):
    return VCFRecord(
        chrom="1",
        pos=100,
        variant_id="synthetic",
        ref=ref,
        alts=alts,
        gt="0/1",
    )


def test_validated_match():
    result = match_validated_alleles(validated(), record())
    assert result.status == "MATCH"
    assert result.reason == "REFERENCE_VALIDATED_ALLELE_SET"


def test_vcf_reference_disagreement():
    result = match_validated_alleles(
        validated(),
        record(ref="A", alts=("G",)),
    )
    assert result.status == "MISMATCH"
    assert result.reason == "VCF_REF_DISAGREES_WITH_REFERENCE"


def test_different_alternate_allele():
    result = match_validated_alleles(
        validated(),
        record(ref="G", alts=("C",)),
    )
    assert result.status == "MISMATCH"
    assert result.reason == "ALLELE_SET_DIFFERENT"


def test_unvalidated_assay():
    assay = ValidatedAssayAlleles(
        status="UNRESOLVED",
        reason="INSUFFICIENT_INFORMATIVE_BASES",
        orientation="UNRESOLVED",
        alleles=(),
        reference_base=None,
    )

    result = match_validated_alleles(assay, record())
    assert result.status == "UNRESOLVED"


def test_multiallelic_record():
    result = match_validated_alleles(
        validated(),
        record(alts=("A", "C")),
    )
    assert result.status == "UNRESOLVED"
    assert result.reason == "UNSUPPORTED_VCF_ALLELES"


def test_invalid_vcf_alleles():
    result = match_validated_alleles(
        validated(),
        record(ref="G", alts=("G",)),
    )
    assert result.status == "UNRESOLVED"
    assert result.reason == "INVALID_VCF_ALLELES"
