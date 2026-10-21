"""Conservative genotype projection for multiallelic SNV records.

Projection describes observed VCF allele indices relative to a
reference-validated assay ALT. It does not infer array intensities,
copy number, or an independently validated assay genotype.
"""

from dataclasses import dataclass

from .multiallelic import MultiallelicMatch
from .vcf import VCFRecord


@dataclass(frozen=True)
class GenotypeProjection:
    status: str
    reason: str
    assay_alt_dosage: int | None
    other_alt_present: bool | None
    ploidy: int | None


def project_multiallelic_genotype(
    match: MultiallelicMatch,
    record: VCFRecord,
) -> GenotypeProjection:
    if match.status != "PARTIAL_MATCH":
        return GenotypeProjection(
            "UNRESOLVED",
            "ALLELE_MATCH_NOT_ESTABLISHED",
            None,
            None,
            None,
        )

    if (
        len(match.matched_alt_indices) != 1
        or not 1 <= match.matched_alt_indices[0] <= len(record.alts)
    ):
        return GenotypeProjection(
            "UNRESOLVED",
            "INVALID_MATCHED_ALT_INDEX",
            None,
            None,
            None,
        )

    gt = record.gt

    if not gt or gt == ".":
        return GenotypeProjection(
            "NO_CALL", "MISSING_GENOTYPE", None, None, None
        )

    if "/" in gt and "|" in gt:
        return GenotypeProjection(
            "UNRESOLVED", "MIXED_GENOTYPE_SEPARATORS",
            None, None, None
        )

    alleles = gt.replace("|", "/").split("/")

    if any(allele == "." for allele in alleles):
        return GenotypeProjection(
            "NO_CALL", "MISSING_GENOTYPE_ALLELE",
            None, None, None
        )

    try:
        indices = tuple(int(allele) for allele in alleles)
    except ValueError:
        return GenotypeProjection(
            "UNRESOLVED", "INVALID_GENOTYPE_INDEX",
            None, None, None
        )

    if any(index < 0 or index > len(record.alts) for index in indices):
        return GenotypeProjection(
            "UNRESOLVED", "GENOTYPE_INDEX_OUT_OF_RANGE",
            None, None, None
        )

    target = match.matched_alt_indices[0]
    dosage = sum(index == target for index in indices)
    other_alt = any(index not in (0, target) for index in indices)

    if other_alt:
        return GenotypeProjection(
            "OTHER_ALT_PRESENT",
            "GENOTYPE_CONTAINS_UNASSAYED_ALT",
            dosage,
            True,
            len(indices),
        )

    if dosage:
        return GenotypeProjection(
            "ASSAY_ALT_PRESENT",
            "TARGET_ALT_OBSERVED",
            dosage,
            False,
            len(indices),
        )

    return GenotypeProjection(
        "ASSAY_REF_ONLY",
        "ALL_CALLED_ALLELES_REFERENCE",
        0,
        False,
        len(indices),
    )
