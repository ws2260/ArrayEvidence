"""Reference-based orientation validation for SNV assay flanks."""

import re
from dataclasses import dataclass
from typing import Protocol

from .manifest import Assay


SOURCE_PATTERN = re.compile(r"\[([ACGT])/([ACGT])\]", re.I)

COMPLEMENT = str.maketrans(
    "ACGTacgt",
    "TGCAtgca",
)


class ReferenceReader(Protocol):
    def fetch(self, reference: str, start: int, end: int) -> str:
        ...

    def get_reference_length(self, reference: str) -> int:
        ...


@dataclass(frozen=True)
class OrientationResult:
    status: str
    reason: str
    forward_matches: int
    reverse_matches: int
    compared_bases: int


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def validate_orientation(
    assay: Assay,
    reference: ReferenceReader,
    flank_length: int = 20,
) -> OrientationResult:
    """Compare SourceSeq flanks against both reference orientations.

    SourceSeq contains a bracketed variant. The reference position
    is 1-based, whereas reference.fetch uses 0-based half-open
    coordinates.

    Only exact, complete flank matches establish orientation.
    """

    if flank_length < 1:
        raise ValueError("flank_length must be positive")

    if assay.coordinate is None:
        return OrientationResult(
            "UNRESOLVED", "UNMAPPED_ASSAY", 0, 0, 0
        )

    matches = list(SOURCE_PATTERN.finditer(assay.source_seq))

    if len(matches) != 1:
        return OrientationResult(
            "UNRESOLVED", "UNSUPPORTED_SOURCE_SEQUENCE", 0, 0, 0
        )

    match = matches[0]

    left = assay.source_seq[:match.start()].upper()
    right = assay.source_seq[match.end():].upper()

    if not left or not right:
        return OrientationResult(
            "UNRESOLVED", "INSUFFICIENT_FLANKS", 0, 0, 0
        )

    if not re.fullmatch(r"[ACGT]+", left + right):
        return OrientationResult(
            "UNRESOLVED", "NON_ACGT_FLANKS", 0, 0, 0
        )

    left = left[-flank_length:]
    right = right[:flank_length]

    chrom, pos = assay.coordinate

    try:
        chrom_length = reference.get_reference_length(chrom)

        if pos - 1 < len(left) or pos + len(right) > chrom_length:
            return OrientationResult(
                "UNRESOLVED", "REFERENCE_BOUNDARY", 0, 0, 0
            )

        forward_left = reference.fetch(
            chrom, pos - 1 - len(left), pos - 1
        ).upper()

        forward_right = reference.fetch(
            chrom, pos, pos + len(right)
        ).upper()

    except (KeyError, ValueError):
        return OrientationResult(
            "UNRESOLVED", "REFERENCE_UNAVAILABLE", 0, 0, 0
        )

    forward = forward_left + forward_right

    reverse = reverse_complement(right) + reverse_complement(left)

    observed = left + right

    forward_ok = observed == forward
    reverse_ok = reverse == forward

    forward_matches = sum(
        a == b for a, b in zip(observed, forward)
    )

    reverse_matches = sum(
        a == b for a, b in zip(reverse, forward)
    )

    compared = len(forward)

    if forward_ok and reverse_ok:
        status, reason = "AMBIGUOUS", "BOTH_ORIENTATIONS_MATCH"
    elif forward_ok:
        status, reason = "FORWARD", "EXACT_FORWARD_FLANK_MATCH"
    elif reverse_ok:
        status, reason = "REVERSE", "EXACT_REVERSE_FLANK_MATCH"
    else:
        status, reason = "NO_MATCH", "NO_EXACT_FLANK_MATCH"

    return OrientationResult(
        status,
        reason,
        forward_matches,
        reverse_matches,
        compared,
    )


def validate_orientation_masked(
    assay: Assay,
    reference: ReferenceReader,
    flank_length: int = 20,
    min_informative_bases: int = 20,
) -> OrientationResult:
    """Validate orientation while treating N as missing flank evidence.

    This is a separate, masked-evidence tier. N never counts as
    a matching base, and all informative positions must agree.
    """

    if flank_length < 1 or min_informative_bases < 1:
        raise ValueError("Flank and informative-base limits must be positive")

    if assay.coordinate is None:
        return OrientationResult(
            "UNRESOLVED", "UNMAPPED_ASSAY", 0, 0, 0
        )

    matches = list(SOURCE_PATTERN.finditer(assay.source_seq))

    if len(matches) != 1:
        return OrientationResult(
            "UNRESOLVED", "UNSUPPORTED_SOURCE_SEQUENCE", 0, 0, 0
        )

    match = matches[0]

    left = assay.source_seq[:match.start()].upper()
    right = assay.source_seq[match.end():].upper()

    if not left or not right:
        return OrientationResult(
            "UNRESOLVED", "INSUFFICIENT_FLANKS", 0, 0, 0
        )

    if not re.fullmatch(r"[ACGTN]+", left + right):
        return OrientationResult(
            "UNRESOLVED", "UNSUPPORTED_FLANK_CHARACTERS", 0, 0, 0
        )

    left = left[-flank_length:]
    right = right[:flank_length]

    chrom, pos = assay.coordinate

    try:
        chrom_length = reference.get_reference_length(chrom)

        if pos - 1 < len(left) or pos + len(right) > chrom_length:
            return OrientationResult(
                "UNRESOLVED", "REFERENCE_BOUNDARY", 0, 0, 0
            )

        ref_left = reference.fetch(
            chrom, pos - 1 - len(left), pos - 1
        ).upper()

        ref_right = reference.fetch(
            chrom, pos, pos + len(right)
        ).upper()

    except (KeyError, ValueError):
        return OrientationResult(
            "UNRESOLVED", "REFERENCE_UNAVAILABLE", 0, 0, 0
        )

    forward_query = left + right

    # Reverse orientation swaps the flanks.
    reverse_query = (
        reverse_complement(right)
        + reverse_complement(left)
    )

    reference_sequence = ref_left + ref_right

    def evaluate(query):
        informative = 0
        matching = 0

        for observed, expected in zip(query, reference_sequence):
            if observed == "N":
                continue

            if observed not in "ACGT":
                continue

            informative += 1

            if observed == expected:
                matching += 1

        return informative, matching

    forward_info, forward_matches = evaluate(forward_query)
    reverse_info, reverse_matches = evaluate(reverse_query)

    forward_ok = (
        forward_info >= min_informative_bases
        and forward_matches == forward_info
    )

    reverse_ok = (
        reverse_info >= min_informative_bases
        and reverse_matches == reverse_info
    )

    compared = max(forward_info, reverse_info)

    if forward_ok and reverse_ok:
        status, reason = "AMBIGUOUS", "MASKED_BOTH_ORIENTATIONS"
    elif forward_ok:
        status, reason = "FORWARD", "MASKED_FORWARD_MATCH"
    elif reverse_ok:
        status, reason = "REVERSE", "MASKED_REVERSE_MATCH"
    elif compared < min_informative_bases:
        status, reason = "UNRESOLVED", "INSUFFICIENT_INFORMATIVE_BASES"
    else:
        status, reason = "NO_MATCH", "MASKED_FLANK_MISMATCH"

    return OrientationResult(
        status,
        reason,
        forward_matches,
        reverse_matches,
        compared,
    )
