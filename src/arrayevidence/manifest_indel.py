"""Conservative reconstruction of indel events from assay SourceSeq.

Requires a unique exact alignment of both flanks to the reference.
This module does not establish VCF allele equivalence.
"""

from dataclasses import dataclass
import re

from .indel_events import SequenceEvent
from .manifest import Assay
from .reference import ReferenceReader, reverse_complement


BRACKET = re.compile(r"\[([^\]]+)\]")


@dataclass(frozen=True)
class ManifestIndelResult:
    status: str
    reason: str
    orientation: str | None
    event: SequenceEvent | None


def reconstruct_manifest_indel(
    assay: Assay,
    reference: ReferenceReader,
    flank_length: int = 15,
    min_flank_length: int = 10,
    search_radius: int = 100,
    max_reference_gap: int = 30,
) -> ManifestIndelResult:
    """Reconstruct a simple reference-oriented indel event.

    All coordinates are zero-based interbase coordinates.
    """

    if assay.coordinate is None:
        return ManifestIndelResult(
            "UNRESOLVED", "UNMAPPED_ASSAY", None, None
        )

    if assay.snp.upper() not in {"[I/D]", "[D/I]"}:
        return ManifestIndelResult(
            "UNRESOLVED", "NOT_INDEL_ASSAY", None, None
        )

    matches = list(BRACKET.finditer(assay.source_seq))

    if len(matches) != 1:
        return ManifestIndelResult(
            "UNRESOLVED", "INVALID_SOURCE_BRACKET", None, None
        )

    match = matches[0]
    tokens = match.group(1).upper().split("/")

    if len(tokens) != 2:
        return ManifestIndelResult(
            "UNRESOLVED", "INVALID_ALLELE_TOKENS", None, None
        )

    if not all(
        token == "-" or (
            token and set(token) <= set("ACGT")
        )
        for token in tokens
    ):
        return ManifestIndelResult(
            "UNRESOLVED", "UNSUPPORTED_ALLELE_TOKENS", None, None
        )

    alleles = tuple(
        "" if token == "-" else token
        for token in tokens
    )

    if alleles[0] == alleles[1]:
        return ManifestIndelResult(
            "UNRESOLVED", "IDENTICAL_ALLELES", None, None
        )

    left = assay.source_seq[:match.start()].upper()
    right = assay.source_seq[match.end():].upper()

    left = left[-flank_length:]
    right = right[:flank_length]

    if (
        len(left) < min_flank_length
        or len(right) < min_flank_length
        or not set(left + right) <= set("ACGT")
    ):
        return ManifestIndelResult(
            "UNRESOLVED", "UNSUPPORTED_FLANKS", None, None
        )

    chrom, pos = assay.coordinate

    try:
        reference_length = reference.get_reference_length(chrom)

        window_start = max(0, pos - search_radius - 1)
        window_end = min(
            reference_length,
            pos + search_radius,
        )

        window = reference.fetch(
            chrom,
            window_start,
            window_end,
        ).upper()

    except (KeyError, ValueError):
        return ManifestIndelResult(
            "UNRESOLVED", "REFERENCE_UNAVAILABLE", None, None
        )

    alignments = set()

    for orientation, flank_left, flank_right in (
        ("FORWARD", left, right),
        (
            "REVERSE",
            reverse_complement(right),
            reverse_complement(left),
        ),
    ):
        for gap in range(max_reference_gap + 1):
            limit = (
                len(window)
                - len(flank_left)
                - gap
                - len(flank_right)
                + 1
            )

            for index in range(max(0, limit)):
                middle_start = index + len(flank_left)
                right_start = middle_start + gap

                if (
                    window[index:middle_start] == flank_left
                    and window[
                        right_start:right_start + len(flank_right)
                    ] == flank_right
                ):
                    alignments.add((
                        orientation,
                        window_start + middle_start,
                        window[middle_start:right_start],
                    ))

    if not alignments:
        return ManifestIndelResult(
            "UNRESOLVED", "NO_FLANK_ALIGNMENT", None, None
        )

    if len(alignments) != 1:
        return ManifestIndelResult(
            "UNRESOLVED", "AMBIGUOUS_FLANK_ALIGNMENT", None, None
        )

    orientation, boundary, reference_allele = next(
        iter(alignments)
    )

    if orientation == "REVERSE":
        alleles = tuple(
            reverse_complement(allele)
            for allele in alleles
        )

    if reference_allele not in alleles:
        return ManifestIndelResult(
            "UNRESOLVED",
            "REFERENCE_ALLELE_NOT_IN_MANIFEST",
            orientation,
            None,
        )

    alternatives = [
        allele
        for allele in alleles
        if allele != reference_allele
    ]

    if len(alternatives) != 1:
        return ManifestIndelResult(
            "UNRESOLVED",
            "NO_UNIQUE_ALTERNATE",
            orientation,
            None,
        )

    alternate = alternatives[0]

    if reference_allele == "" and alternate:
        event_type = "INSERTION"
    elif reference_allele and alternate == "":
        event_type = "DELETION"
    else:
        return ManifestIndelResult(
            "UNRESOLVED",
            "COMPLEX_ALLELE_REPRESENTATION",
            orientation,
            None,
        )

    event = SequenceEvent(
        chrom=chrom,
        start=boundary,
        deleted=reference_allele,
        inserted=alternate,
        status=event_type,
    )

    return ManifestIndelResult(
        "VALIDATED",
        "UNIQUE_REFERENCE_ALIGNMENT",
        orientation,
        event,
    )
