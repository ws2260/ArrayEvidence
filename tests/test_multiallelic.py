from arrayevidence.assay_alleles import ValidatedAssayAlleles
from arrayevidence.multiallelic import match_multiallelic_snv
from arrayevidence.vcf import VCFRecord


def assay(ref="A", alt="G", status="VALIDATED"):
    return ValidatedAssayAlleles(
        status=status,
        reason="SYNTHETIC",
        orientation="FORWARD",
        alleles=tuple(sorted((ref, alt))),
        reference_base=ref,
    )


def record(ref="A", alts=("C", "G"), gt="1/2"):
    return VCFRecord(
        chrom="1",
        pos=100,
        variant_id="synthetic",
        ref=ref,
        alts=alts,
        gt=gt,
    )


def test_assay_alt_present_at_second_alt_index():
    result = match_multiallelic_snv(assay(), record())

    assert result.status == "PARTIAL_MATCH"
    assert result.matched_alt_indices == (2,)


def test_assay_alt_present_at_first_alt_index():
    result = match_multiallelic_snv(
        assay(),
        record(alts=("G", "C")),
    )

    assert result.status == "PARTIAL_MATCH"
    assert result.matched_alt_indices == (1,)


def test_assay_alt_absent():
    result = match_multiallelic_snv(
        assay(alt="T"),
        record(),
    )

    assert result.status == "MISMATCH"
    assert result.reason == "ASSAY_ALT_NOT_IN_VCF"


def test_reference_mismatch():
    result = match_multiallelic_snv(
        assay(),
        record(ref="T", alts=("C", "G")),
    )

    assert result.status == "MISMATCH"
    assert result.reason == "VCF_REF_DISAGREES_WITH_REFERENCE"


def test_indel_containing_record_is_unresolved():
    result = match_multiallelic_snv(
        assay(),
        record(alts=("C", "AT")),
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "UNSUPPORTED_VCF_ALLELES"


def test_duplicate_alt_is_unresolved():
    result = match_multiallelic_snv(
        assay(),
        record(alts=("G", "G")),
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "DUPLICATE_VCF_ALLELES"


def test_unvalidated_assay_is_unresolved():
    result = match_multiallelic_snv(
        assay(status="UNRESOLVED"),
        record(),
    )

    assert result.status == "UNRESOLVED"


def test_biallelic_record_not_handled_here():
    result = match_multiallelic_snv(
        assay(),
        record(alts=("G",)),
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "NOT_MULTIALLELIC"


def test_partial_match_does_not_depend_on_genotype():
    for gt in ("0/0", "0/1", "1/2", "./."):
        result = match_multiallelic_snv(
            assay(),
            record(gt=gt),
        )

        assert result.status == "PARTIAL_MATCH"
        assert result.matched_alt_indices == (2,)
