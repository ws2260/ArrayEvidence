"""Conservative reference-validated matching for multiallelic SNVs."""

from dataclasses import dataclass

from .assay_alleles import ValidatedAssayAlleles
from .vcf import VCFRecord


@dataclass(frozen=True)
class MultiallelicMatch:
    status: str
    reason: str
    matched_alt_indices: tuple[int, ...]


def match_multiallelic_snv(
    assay: ValidatedAssayAlleles,
    record: VCFRecord,
) -> MultiallelicMatch:
    """Compare validated biallelic assay alleles to a multiallelic SNV.

    ALT indices are 1-based, following VCF genotype conventions.

    A partial allele match does not establish that the genotype can be
    assigned to the assay. Genotype interpretation is a separate step.

    Assay and VCF coordinates and reference assembly must already agree.
    """

    if assay.status != "VALIDATED":
        return MultiallelicMatch(
            "UNRESOLVED",
            "ASSAY_ALLELES_NOT_VALIDATED",
            (),
        )

    if len(record.alts) < 2:
        return MultiallelicMatch(
            "UNRESOLVED",
            "NOT_MULTIALLELIC",
            (),
        )

    ref = record.ref.upper()
    alts = tuple(alt.upper() for alt in record.alts)

    if (
        len(ref) != 1
        or ref not in "ACGT"
        or any(len(alt) != 1 or alt not in "ACGT" for alt in alts)
    ):
        return MultiallelicMatch(
            "UNRESOLVED",
            "UNSUPPORTED_VCF_ALLELES",
            (),
        )

    if len(set((ref, *alts))) != len(alts) + 1:
        return MultiallelicMatch(
            "UNRESOLVED",
            "DUPLICATE_VCF_ALLELES",
            (),
        )

    if ref != assay.reference_base:
        return MultiallelicMatch(
            "MISMATCH",
            "VCF_REF_DISAGREES_WITH_REFERENCE",
            (),
        )

    assay_alts = set(assay.alleles) - {ref}

    if len(assay_alts) != 1:
        return MultiallelicMatch(
            "UNRESOLVED",
            "INVALID_ASSAY_ALLELE_SET",
            (),
        )

    matched_indices = tuple(
        index
        for index, alt in enumerate(alts, start=1)
        if alt in assay_alts
    )

    if not matched_indices:
        return MultiallelicMatch(
            "MISMATCH",
            "ASSAY_ALT_NOT_IN_VCF",
            (),
        )

    return MultiallelicMatch(
        "PARTIAL_MATCH",
        "ASSAY_ALT_PRESENT_IN_MULTIALLELIC_VCF",
        matched_indices,
    )
