"""Reference-aware left normalization of simple insertion/deletion events.

Coordinates are zero-based interbase positions. This module does not
interpret Illumina assay alleles or complex substitutions.
"""

from dataclasses import dataclass

from .indel_events import SequenceEvent
from .reference import ReferenceReader


@dataclass(frozen=True)
class NormalizedIndel:
    status: str
    reason: str
    event: SequenceEvent | None


def left_normalize_indel(
    event: SequenceEvent,
    reference: ReferenceReader,
) -> NormalizedIndel:
    """Left-normalize a simple insertion or deletion.

    The reference must use the same assembly and contig naming
    convention as the event.

    Returns UNRESOLVED if reference validation fails.
    """

    if event.status not in {"INSERTION", "DELETION"}:
        return NormalizedIndel(
            "UNRESOLVED",
            "NOT_SIMPLE_INDEL",
            None,
        )

    chrom = event.chrom
    start = event.start

    if start < 0:
        return NormalizedIndel(
            "UNRESOLVED",
            "INVALID_COORDINATE",
            None,
        )

    try:
        reference_length = reference.get_reference_length(chrom)
    except (KeyError, ValueError):
        return NormalizedIndel(
            "UNRESOLVED",
            "REFERENCE_CONTIG_UNAVAILABLE",
            None,
        )

    if start > reference_length:
        return NormalizedIndel(
            "UNRESOLVED",
            "COORDINATE_OUT_OF_RANGE",
            None,
        )

    if event.status == "INSERTION":
        sequence = event.inserted.upper()

        if (
            not sequence
            or event.deleted
            or not set(sequence) <= set("ACGT")
        ):
            return NormalizedIndel(
                "UNRESOLVED",
                "INVALID_INSERTION",
                None,
            )

    else:
        sequence = event.deleted.upper()

        if (
            not sequence
            or event.inserted
            or not set(sequence) <= set("ACGT")
        ):
            return NormalizedIndel(
                "UNRESOLVED",
                "INVALID_DELETION",
                None,
            )

        if start + len(sequence) > reference_length:
            return NormalizedIndel(
                "UNRESOLVED",
                "DELETION_OUT_OF_RANGE",
                None,
            )

        observed = reference.fetch(
            chrom,
            start,
            start + len(sequence),
        ).upper()

        if observed != sequence:
            return NormalizedIndel(
                "UNRESOLVED",
                "DELETION_REFERENCE_MISMATCH",
                None,
            )

    # Shift left while the preceding reference base equals
    # the last base of the inserted/deleted sequence.
    #
    # Rotate the event sequence when shifting.
    while start > 0:
        previous = reference.fetch(
            chrom,
            start - 1,
            start,
        ).upper()

        if previous != sequence[-1]:
            break

        sequence = sequence[-1] + sequence[:-1]
        start -= 1

    if event.status == "INSERTION":
        normalized = SequenceEvent(
            chrom=chrom,
            start=start,
            deleted="",
            inserted=sequence,
            status="INSERTION",
        )
    else:
        normalized = SequenceEvent(
            chrom=chrom,
            start=start,
            deleted=sequence,
            inserted="",
            status="DELETION",
        )

    return NormalizedIndel(
        "NORMALIZED",
        "REFERENCE_VALIDATED",
        normalized,
    )
