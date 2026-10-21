"""Reference-validated equivalence of simple manifest and VCF indels."""

from dataclasses import dataclass

from .indel_events import minimal_sequence_event
from .indel_normalize import left_normalize_indel
from .manifest_indel import ManifestIndelResult
from .reference import ReferenceReader
from .vcf import VCFRecord


@dataclass(frozen=True)
class IndelMatch:
    status: str
    reason: str


def match_indel(
    manifest: ManifestIndelResult,
    record: VCFRecord,
    reference: ReferenceReader,
) -> IndelMatch:
    """Compare independently validated, normalized simple indel events.

    Assumes the caller has verified reference assembly compatibility.
    """

    if manifest.status != "VALIDATED" or manifest.event is None:
        return IndelMatch(
            "UNRESOLVED",
            "MANIFEST_EVENT_NOT_VALIDATED",
        )

    if len(record.alts) != 1:
        return IndelMatch(
            "UNRESOLVED",
            "MULTIALLELIC_VCF_NOT_SUPPORTED",
        )

    if record.chrom != manifest.event.chrom:
        return IndelMatch(
            "UNRESOLVED",
            "CHROMOSOME_MISMATCH",
        )

    ref = record.ref.upper()
    alt = record.alts[0].upper()

    if (
        not ref
        or not alt
        or not set(ref + alt) <= set("ACGT")
    ):
        return IndelMatch(
            "UNRESOLVED",
            "UNSUPPORTED_VCF_ALLELES",
        )

    # Validate the entire VCF REF allele against the reference.
    try:
        observed = reference.fetch(
            record.chrom,
            record.pos - 1,
            record.pos - 1 + len(ref),
        ).upper()
    except (KeyError, ValueError):
        return IndelMatch(
            "UNRESOLVED",
            "VCF_REFERENCE_UNAVAILABLE",
        )

    if observed != ref:
        return IndelMatch(
            "UNRESOLVED",
            "VCF_REF_DISAGREES_WITH_REFERENCE",
        )

    vcf_event = minimal_sequence_event(
        record.chrom,
        record.pos,
        ref,
        alt,
    )

    if vcf_event.status not in {"INSERTION", "DELETION"}:
        return IndelMatch(
            "UNRESOLVED",
            "VCF_NOT_SIMPLE_INDEL",
        )

    manifest_normalized = left_normalize_indel(
        manifest.event,
        reference,
    )

    vcf_normalized = left_normalize_indel(
        vcf_event,
        reference,
    )

    if (
        manifest_normalized.status != "NORMALIZED"
        or vcf_normalized.status != "NORMALIZED"
    ):
        return IndelMatch(
            "UNRESOLVED",
            "NORMALIZATION_FAILED",
        )

    if manifest_normalized.event == vcf_normalized.event:
        return IndelMatch(
            "MATCH",
            "NORMALIZED_SEQUENCE_EVENT_EQUIVALENT",
        )

    return IndelMatch(
        "MISMATCH",
        "NORMALIZED_SEQUENCE_EVENT_DIFFERENT",
    )
