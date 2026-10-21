"""Reference-validated allele representation for biallelic SNV assays."""

import re
from dataclasses import dataclass

from .manifest import Assay
from .reference import (
    ReferenceReader,
    reverse_complement,
    validate_orientation_masked,
)


SNP_PATTERN = re.compile(r"^\[([ACGT])/([ACGT])\]$", re.I)
SOURCE_PATTERN = re.compile(r"\[([ACGT])/([ACGT])\]", re.I)


@dataclass(frozen=True)
class ValidatedAssayAlleles:
    status: str
    reason: str
    orientation: str
    alleles: tuple[str, ...]
    reference_base: str | None


def validate_assay_alleles(
    assay: Assay,
    reference: ReferenceReader,
    flank_length: int = 20,
    min_informative_bases: int = 20,
) -> ValidatedAssayAlleles:
    """Resolve an assay's biallelic SNV alleles to reference orientation."""

    unresolved = lambda reason, orientation="UNRESOLVED": ValidatedAssayAlleles(
        status="UNRESOLVED",
        reason=reason,
        orientation=orientation,
        alleles=(),
        reference_base=None,
    )

    if assay.coordinate is None:
        return unresolved("UNMAPPED_ASSAY")

    snp_match = SNP_PATTERN.fullmatch(assay.snp.upper())
    source_matches = SOURCE_PATTERN.findall(assay.source_seq.upper())

    if snp_match is None or len(source_matches) != 1:
        return unresolved("UNSUPPORTED_ASSAY_ALLELES")

    orientation = validate_orientation_masked(
        assay,
        reference,
        flank_length=flank_length,
        min_informative_bases=min_informative_bases,
    )

    if orientation.status not in {"FORWARD", "REVERSE"}:
        return unresolved(
            orientation.reason,
            orientation.status,
        )

    source_alleles = set(source_matches[0])

    if orientation.status == "REVERSE":
        source_alleles = {
            reverse_complement(allele)
            for allele in source_alleles
        }

    chrom, pos = assay.coordinate

    try:
        reference_base = reference.fetch(
            chrom, pos - 1, pos
        ).upper()
    except (KeyError, ValueError):
        return unresolved("REFERENCE_UNAVAILABLE", orientation.status)

    if reference_base not in source_alleles:
        return unresolved(
            "REFERENCE_BASE_NOT_IN_ASSAY_ALLELES",
            orientation.status,
        )

    return ValidatedAssayAlleles(
        status="VALIDATED",
        reason="REFERENCE_AND_FLANKS_CONCORDANT",
        orientation=orientation.status,
        alleles=tuple(sorted(source_alleles)),
        reference_base=reference_base,
    )
