"""Conservative allele matching for biallelic SNV assays.

Coordinate matching alone does not establish allele equivalence.
"""

import re
from dataclasses import dataclass

from .manifest import Assay
from .vcf import VCFRecord


SNP_PATTERN = re.compile(r"^\[([ACGT])/([ACGT])\]$")

COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class AlleleMatch:
    status: str
    reason: str


def match_alleles(assay: Assay, record: VCFRecord) -> AlleleMatch:
    """Compare a manifest SNP definition with a VCF record.

    This is a conservative preliminary matcher. It does not infer
    strand orientation or validate alleles against the reference.
    """

    if assay.coordinate is None:
        return AlleleMatch("UNRESOLVED", "UNMAPPED_ASSAY")

    if assay.coordinate != record.coordinate:
        return AlleleMatch("MISMATCH", "DIFFERENT_COORDINATE")

    match = SNP_PATTERN.fullmatch(assay.snp.upper())

    if match is None:
        return AlleleMatch("UNRESOLVED", "UNSUPPORTED_ASSAY_ALLELES")

    if len(record.ref) != 1 or len(record.alts) != 1:
        return AlleleMatch("UNRESOLVED", "UNSUPPORTED_VCF_ALLELES")

    if record.ref not in "ACGT" or record.alts[0] not in "ACGT":
        return AlleleMatch("UNRESOLVED", "UNSUPPORTED_VCF_ALLELES")

    assay_alleles = set(match.groups())
    vcf_alleles = {record.ref, record.alts[0]}

    if assay_alleles == vcf_alleles:
        return AlleleMatch("MATCH", "DIRECT_ALLELE_SET")

    complemented = {
        allele.translate(COMPLEMENT)
        for allele in assay_alleles
    }

    if complemented == vcf_alleles:
        return AlleleMatch("UNRESOLVED", "POSSIBLE_REVERSE_STRAND")

    return AlleleMatch("MISMATCH", "ALLELE_SET_DIFFERENT")
