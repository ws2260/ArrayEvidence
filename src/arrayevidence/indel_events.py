"""Minimal sequence-change representation for VCF alleles.

This module trims common REF/ALT context. It does not perform
repeat-aware left alignment or establish assay equivalence.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class SequenceEvent:
    chrom: str
    start: int
    deleted: str
    inserted: str
    status: str

    @property
    def length_change(self) -> int:
        return len(self.inserted) - len(self.deleted)


def minimal_sequence_event(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> SequenceEvent:
    """Return a minimally trimmed, zero-based interbase event.

    VCF POS is one-based. The returned start is zero-based,
    immediately before the first changed reference base.

    This function requires explicit A/C/G/T alleles.
    """

    ref = ref.upper()
    alt = alt.upper()

    if (
        pos < 1
        or not ref
        or not alt
        or not set(ref + alt) <= set("ACGT")
    ):
        return SequenceEvent(
            chrom, pos - 1, ref, alt, "UNRESOLVED"
        )

    start = pos - 1

    # Remove shared prefix.
    while ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        start += 1

    # Remove shared suffix.
    while ref and alt and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    if not ref and not alt:
        status = "NO_CHANGE"
    elif not ref:
        status = "INSERTION"
    elif not alt:
        status = "DELETION"
    elif len(ref) == len(alt):
        status = "SUBSTITUTION"
    else:
        status = "COMPLEX_INDEL"

    return SequenceEvent(
        chrom=chrom,
        start=start,
        deleted=ref,
        inserted=alt,
        status=status,
    )
