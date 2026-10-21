"""Pure conversion from reconciliation results to versioned evidence rows.

This module does not perform allele matching, genotype calling,
clinical annotation, file I/O, or participant identifier generation.

Input order is preserved. Every candidate retains its original
VCF record representation and 1-based ALT index.
"""

from collections.abc import Iterable

from .evidence_schema import (
    SCHEMA_VERSION,
    AlleleCandidateRow,
    AssayEvidenceRow,
)
from .reconcile import ReconciliationResult


def build_evidence_rows(
    sample_id: str,
    results: Iterable[ReconciliationResult],
) -> tuple[list[AssayEvidenceRow], list[AlleleCandidateRow]]:
    """Convert reconciled assays into two linked technical evidence tables.

    One AssayEvidenceRow is created per input result. Each established
    source allele becomes an AlleleCandidateRow, including every
    equivalent allele in an ambiguous indel match.

    Missing source alleles never imply homozygous reference genotypes.
    """

    if not isinstance(sample_id, str) or not sample_id.strip():
        raise ValueError("sample_id must be a non-empty string")

    assay_rows: list[AssayEvidenceRow] = []
    candidate_rows: list[AlleleCandidateRow] = []

    for result in results:
        coordinate = result.coordinate

        if coordinate is None:
            manifest_chrom = None
            manifest_pos = None
        else:
            manifest_chrom, manifest_pos = coordinate

        assay_rows.append(
            AssayEvidenceRow(
                schema_version=SCHEMA_VERSION,
                sample_id=sample_id,
                assay_id=result.assay_id,
                manifest_chrom=manifest_chrom,
                manifest_pos=manifest_pos,
                coordinate_status=result.status,
                coordinate_evidence=result.coordinate_evidence,
                coordinate_candidate_count=result.candidate_record_count,
                coordinate_genotype_evidence=result.genotype_evidence,
                allele_evidence=result.allele_evidence,
                allele_reason=result.allele_reason,
                indel_evidence=result.indel_evidence,
                projection_status=result.projection_status,
                projection_reason=result.projection_reason,
                assay_alt_dosage=result.assay_alt_dosage,
                other_alt_present=result.other_alt_present,
                ploidy=result.ploidy,
            )
        )

        candidate_status = (
            result.indel_evidence
            if result.indel_evidence in {"UNIQUE", "MULTIPLE"}
            else result.allele_evidence
        )

        for candidate_index, source in enumerate(
            result.allele_sources,
            start=1,
        ):
            record = source.record
            alt_index = source.alt_index

            if (
                type(alt_index) is not int
                or alt_index < 1
                or alt_index > len(record.alts)
            ):
                raise ValueError(
                    f"Invalid ALT index {alt_index!r} for "
                    f"assay {result.assay_id!r}"
                )

            candidate_rows.append(
                AlleleCandidateRow(
                    schema_version=SCHEMA_VERSION,
                    sample_id=sample_id,
                    assay_id=result.assay_id,
                    candidate_index=candidate_index,
                    vcf_chrom=record.chrom,
                    vcf_pos=record.pos,
                    vcf_id=record.variant_id,
                    vcf_ref=record.ref,
                    vcf_alt=record.alts[alt_index - 1],
                    vcf_alt_index=alt_index,
                    vcf_gt=record.gt,
                    candidate_status=candidate_status,
                )
            )

    return assay_rows, candidate_rows
