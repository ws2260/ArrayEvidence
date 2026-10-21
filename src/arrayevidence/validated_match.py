"""Reference-validated matching between assay alleles and VCF records."""

from dataclasses import dataclass

from .assay_alleles import ValidatedAssayAlleles
from .vcf import VCFRecord


@dataclass(frozen=True)
class ValidatedAlleleMatch:
    status: str
    reason: str


def match_validated_alleles(
    assay_alleles: ValidatedAssayAlleles,
    record: VCFRecord,
) -> ValidatedAlleleMatch:
    """Compare validated reference-forward SNV alleles with a VCF record.

    This function assumes that assay and VCF coordinates have already
    been checked and that both use the same reference assembly.
    """

    if assay_alleles.status != "VALIDATED":
        return ValidatedAlleleMatch(
            "UNRESOLVED",
            "ASSAY_ALLELES_NOT_VALIDATED",
        )

    if len(record.ref) != 1 or len(record.alts) != 1:
        return ValidatedAlleleMatch(
            "UNRESOLVED",
            "UNSUPPORTED_VCF_ALLELES",
        )

    ref = record.ref.upper()
    alt = record.alts[0].upper()

    if (
        ref not in "ACGT"
        or alt not in "ACGT"
        or ref == alt
    ):
        return ValidatedAlleleMatch(
            "UNRESOLVED",
            "INVALID_VCF_ALLELES",
        )

    if ref != assay_alleles.reference_base:
        return ValidatedAlleleMatch(
            "MISMATCH",
            "VCF_REF_DISAGREES_WITH_REFERENCE",
        )

    if {ref, alt} == set(assay_alleles.alleles):
        return ValidatedAlleleMatch(
            "MATCH",
            "REFERENCE_VALIDATED_ALLELE_SET",
        )

    return ValidatedAlleleMatch(
        "MISMATCH",
        "ALLELE_SET_DIFFERENT",
    )
