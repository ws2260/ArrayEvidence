"""Conservative genotype projection for matched indel alleles.

Indel equivalence is established upstream by reference-normalized
candidate discovery. This module projects an observed VCF genotype
only when exactly one matching ALT allele has been identified.

It does not infer array intensity, genotype accuracy, or copy number.
"""

from .genotype_projection import (
    GenotypeProjection,
    project_multiallelic_genotype,
)
from .indel_index import IndelCandidateResult
from .multiallelic import MultiallelicMatch


def project_indel_candidate(
    candidate: IndelCandidateResult,
) -> GenotypeProjection:
    """Project a uniquely matched indel ALT onto its VCF genotype."""

    if candidate.status != "UNIQUE" or len(candidate.alleles) != 1:
        return GenotypeProjection(
            status="UNRESOLVED",
            reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
            assay_alt_dosage=None,
            other_alt_present=None,
            ploidy=None,
        )

    allele = candidate.alleles[0]

    match = MultiallelicMatch(
        status="PARTIAL_MATCH",
        reason="REFERENCE_NORMALIZED_INDEL_MATCH",
        matched_alt_indices=(allele.alt_index,),
    )

    return project_multiallelic_genotype(
        match,
        allele.record,
    )
