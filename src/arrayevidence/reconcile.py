"""Reference-aware assay-to-VCF reconciliation.

Coordinate evidence, allele equivalence, and genotype projection
are deliberately separate evidence dimensions.

The caller must ensure that manifest, VCF, and FASTA use the same
reference assembly and chromosome naming convention.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

from .assay_alleles import validate_assay_alleles
from .indel_index import IndelCandidateIndex
from .indel_normalize import left_normalize_indel
from .indel_projection import project_indel_candidate
from .manifest_indel import reconstruct_manifest_indel
from .genotype_projection import project_multiallelic_genotype
from .manifest import Assay
from .multiallelic import match_multiallelic_snv
from .reference import ReferenceReader
from .validated_match import match_validated_alleles
from .vcf import VCFRecord


@dataclass(frozen=True)
class ReconciliationResult:
    assay_id: str
    coordinate: tuple[str, int] | None
    status: str
    candidate_record_count: int
    coordinate_evidence: str
    genotype_evidence: str
    allele_evidence: str = "NOT_EVALUATED"
    allele_reason: str = "REFERENCE_NOT_PROVIDED"
    projection_status: str = "NOT_EVALUATED"
    projection_reason: str = "NOT_EVALUATED"
    assay_alt_dosage: int | None = None
    other_alt_present: bool | None = None
    ploidy: int | None = None
    matched_alt_indices: tuple[int, ...] = ()
    indel_evidence: str = "NOT_EVALUATED"


def reconcile(
    assays: Iterable[Assay],
    records: Iterable[VCFRecord],
    reference: ReferenceReader | None = None,
    flank_length: int = 20,
    min_informative_bases: int = 20,
    indel_flank_length: int = 15,
    indel_min_flank_length: int = 10,
) -> list[ReconciliationResult]:

    records = tuple(records)
    assays = tuple(assays)
    indel_index = (
        IndelCandidateIndex(records, reference)
        if reference is not None
        and any(a.snp.upper() in {'[I/D]', '[D/I]'} for a in assays)
        else None
    )

    by_coordinate = defaultdict(list)

    for record in records:
        by_coordinate[record.coordinate].append(record)

    results = []

    for assay in assays:
        coordinate = assay.coordinate

        allele_evidence = "NOT_EVALUATED"
        allele_reason = "REFERENCE_NOT_PROVIDED"

        projection_status = "NOT_EVALUATED"
        projection_reason = "NOT_EVALUATED"

        assay_alt_dosage = None
        other_alt_present = None
        ploidy = None
        matched_alt_indices = ()
        indel_evidence = "NOT_EVALUATED"

        if coordinate is None:
            candidates = []
            status = "UNMAPPED_ASSAY"
            coordinate_evidence = "UNMAPPED"
            genotype_evidence = "NOT_EVALUATED"
            allele_reason = "UNMAPPED_ASSAY"

        else:
            candidates = by_coordinate.get(coordinate, [])

            if not candidates:
                status = "VCF_ABSENT"
                coordinate_evidence = "ABSENT"
                genotype_evidence = "NOT_EVALUATED"
                allele_reason = "VCF_ABSENT"

            else:
                coordinate_evidence = (
                    "MULTIPLE"
                    if len(candidates) > 1
                    else "PRESENT"
                )

                genotype_evidence = (
                    "CALLED"
                    if any(
                        record.call_status != "NO_CALL"
                        for record in candidates
                    )
                    else "NO_CALL"
                )

                status = (
                    "VCF_CALLED"
                    if genotype_evidence == "CALLED"
                    else "VCF_NO_CALL"
                )

                if len(candidates) > 1:
                    allele_evidence = "UNRESOLVED"
                    allele_reason = "MULTIPLE_VCF_CANDIDATES"

                elif (
                    reference is not None
                    and assay.snp.upper() not in {"[I/D]", "[D/I]"}
                ):
                    record = candidates[0]

                    validated = validate_assay_alleles(
                        assay,
                        reference,
                        flank_length=flank_length,
                        min_informative_bases=min_informative_bases,
                    )

                    if len(record.alts) > 1:
                        match = match_multiallelic_snv(
                            validated,
                            record,
                        )

                        allele_evidence = match.status
                        matched_alt_indices = match.matched_alt_indices

                        allele_reason = (
                            validated.reason
                            if validated.status != "VALIDATED"
                            else match.reason
                        )

                        if match.status == "PARTIAL_MATCH":
                            projection = project_multiallelic_genotype(
                                match,
                                record,
                            )

                            projection_status = projection.status
                            projection_reason = projection.reason
                            assay_alt_dosage = (
                                projection.assay_alt_dosage
                            )
                            other_alt_present = (
                                projection.other_alt_present
                            )
                            ploidy = projection.ploidy

                    else:
                        match = match_validated_alleles(
                            validated,
                            record,
                        )

                        allele_evidence = match.status

                        allele_reason = (
                            validated.reason
                            if validated.status != "VALIDATED"
                            else match.reason
                        )

        if (
            indel_index is not None
            and assay.snp.upper() in {"[I/D]", "[D/I]"}
        ):
            reconstructed = reconstruct_manifest_indel(
                assay,
                reference,
                flank_length=indel_flank_length,
                min_flank_length=indel_min_flank_length,
            )

            if reconstructed.status != "VALIDATED" or reconstructed.event is None:
                indel_evidence = "UNRESOLVED"
                allele_evidence = "UNRESOLVED"
                allele_reason = reconstructed.reason
                projection_status = "NOT_EVALUATED"
                projection_reason = "MANIFEST_INDEL_NOT_VALIDATED"
                assay_alt_dosage = None
                other_alt_present = None
                ploidy = None
                matched_alt_indices = ()
            else:
                normalized = left_normalize_indel(
                    reconstructed.event, reference
                )

                if normalized.status != "NORMALIZED" or normalized.event is None:
                    indel_evidence = "UNRESOLVED"
                    allele_evidence = "UNRESOLVED"
                    allele_reason = normalized.reason
                    projection_status = "NOT_EVALUATED"
                    projection_reason = "INDEL_NORMALIZATION_FAILED"
                    assay_alt_dosage = None
                    other_alt_present = None
                    ploidy = None
                    matched_alt_indices = ()
                else:
                    candidate = indel_index.find(normalized.event)
                    indel_evidence = candidate.status
                    allele_reason = candidate.reason
                    matched_alt_indices = tuple(
                        a.alt_index for a in candidate.alleles
                    )

                    if candidate.status == "UNIQUE":
                        allele_evidence = "MATCH"
                        projection = project_indel_candidate(candidate)
                        projection_status = projection.status
                        projection_reason = projection.reason
                        assay_alt_dosage = projection.assay_alt_dosage
                        other_alt_present = projection.other_alt_present
                        ploidy = projection.ploidy
                    else:
                        allele_evidence = "UNRESOLVED"
                        projection_status = (
                            "UNRESOLVED"
                            if candidate.status == "MULTIPLE"
                            else "NOT_EVALUATED"
                        )
                        projection_reason = "UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED"
                        assay_alt_dosage = None
                        other_alt_present = None
                        ploidy = None

        results.append(
            ReconciliationResult(
                assay_id=assay.ilmn_id,
                coordinate=coordinate,
                status=status,
                candidate_record_count=len(candidates),
                coordinate_evidence=coordinate_evidence,
                genotype_evidence=genotype_evidence,
                allele_evidence=allele_evidence,
                allele_reason=allele_reason,
                projection_status=projection_status,
                projection_reason=projection_reason,
                assay_alt_dosage=assay_alt_dosage,
                other_alt_present=other_alt_present,
                ploidy=ploidy,
                matched_alt_indices=matched_alt_indices,
                indel_evidence=indel_evidence,
            )
        )

    return results
