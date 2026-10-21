from arrayevidence.alleles import match_alleles
from arrayevidence.manifest import Assay
from arrayevidence.vcf import VCFRecord


def assay(snp, chrom="1", pos=100):
    return Assay(
        ilmn_id="synthetic_probe",
        name="synthetic_variant",
        genome_build="38",
        chrom=chrom,
        map_info=str(pos),
        snp=snp,
        source_seq="",
    )


def record(ref, alt, chrom="1", pos=100):
    return VCFRecord(
        chrom=chrom,
        pos=pos,
        variant_id="synthetic_variant",
        ref=ref,
        alts=(alt,),
        gt="0/1",
    )


def test_direct_allele_match():
    result = match_alleles(
        assay("[A/G]"),
        record("A", "G"),
    )
    assert result.status == "MATCH"
    assert result.reason == "DIRECT_ALLELE_SET"


def test_reverse_strand_requires_validation():
    result = match_alleles(
        assay("[A/G]"),
        record("T", "C"),
    )
    assert result.status == "UNRESOLVED"
    assert result.reason == "POSSIBLE_REVERSE_STRAND"


def test_different_alleles():
    result = match_alleles(
        assay("[A/G]"),
        record("A", "C"),
    )
    assert result.status == "MISMATCH"


def test_different_coordinate():
    result = match_alleles(
        assay("[A/G]"),
        record("A", "G", pos=200),
    )
    assert result.reason == "DIFFERENT_COORDINATE"


def test_multiallelic_vcf_unresolved():
    vcf = VCFRecord(
        chrom="1",
        pos=100,
        variant_id="synthetic",
        ref="A",
        alts=("G", "C"),
        gt="0/2",
    )

    result = match_alleles(assay("[A/G]"), vcf)

    assert result.status == "UNRESOLVED"
    assert result.reason == "UNSUPPORTED_VCF_ALLELES"
