import pytest

from arrayevidence.indel_events import SequenceEvent
from arrayevidence.indel_normalize import left_normalize_indel


class SyntheticReference:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, reference, start, end):
        if reference != "1":
            raise KeyError(reference)
        return self.sequence[start:end]

    def get_reference_length(self, reference):
        if reference != "1":
            raise KeyError(reference)
        return len(self.sequence)


def event(start, deleted="", inserted=""):
    status = "INSERTION" if inserted else "DELETION"
    return SequenceEvent(
        "1", start, deleted, inserted, status
    )


def test_homopolymer_insertion_equivalence():
    reference = SyntheticReference("CAAAAG")

    first = left_normalize_indel(
        event(2, inserted="A"),
        reference,
    )
    second = left_normalize_indel(
        event(5, inserted="A"),
        reference,
    )

    assert first.status == "NORMALIZED"
    assert first.event == second.event
    assert first.event.start == 1


def test_homopolymer_deletion_equivalence():
    reference = SyntheticReference("CAAAAG")

    first = left_normalize_indel(
        event(1, deleted="A"),
        reference,
    )
    second = left_normalize_indel(
        event(4, deleted="A"),
        reference,
    )

    assert first.event == second.event
    assert first.event.start == 1


def test_tandem_repeat_insertion_rotation():
    reference = SyntheticReference("CATATG")

    first = left_normalize_indel(
        event(1, inserted="AT"),
        reference,
    )
    second = left_normalize_indel(
        event(3, inserted="AT"),
        reference,
    )

    assert first.event == second.event
    assert first.event.start == 1


def test_tandem_repeat_deletion_rotation():
    reference = SyntheticReference("CATATG")

    first = left_normalize_indel(
        event(1, deleted="AT"),
        reference,
    )
    second = left_normalize_indel(
        event(3, deleted="AT"),
        reference,
    )

    assert first.event == second.event
    assert first.event.start == 1


def test_nonrepeat_insertion_not_shifted():
    reference = SyntheticReference("ACGT")

    result = left_normalize_indel(
        event(2, inserted="T"),
        reference,
    )

    assert result.event.start == 2
    assert result.event.inserted == "T"


def test_deletion_reference_mismatch():
    reference = SyntheticReference("ACGT")

    result = left_normalize_indel(
        event(1, deleted="G"),
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "DELETION_REFERENCE_MISMATCH"


def test_deletion_out_of_range():
    reference = SyntheticReference("ACGT")

    result = left_normalize_indel(
        event(3, deleted="GT"),
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "DELETION_OUT_OF_RANGE"


def test_invalid_coordinate():
    reference = SyntheticReference("ACGT")

    result = left_normalize_indel(
        event(-1, inserted="A"),
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "INVALID_COORDINATE"


def test_unsupported_complex_event():
    reference = SyntheticReference("ACGT")

    complex_event = SequenceEvent(
        "1", 1, "C", "AT", "COMPLEX_INDEL"
    )

    result = left_normalize_indel(
        complex_event,
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "NOT_SIMPLE_INDEL"


def test_insertion_at_reference_start():
    reference = SyntheticReference("ACGT")

    result = left_normalize_indel(
        event(0, inserted="A"),
        reference,
    )

    assert result.status == "NORMALIZED"
    assert result.event.start == 0


def test_deletion_at_reference_start():
    reference = SyntheticReference("ACGT")

    result = left_normalize_indel(
        event(0, deleted="A"),
        reference,
    )

    assert result.status == "NORMALIZED"
    assert result.event.start == 0


def test_unknown_reference_contig():
    reference = SyntheticReference("ACGT")

    result = left_normalize_indel(
        SequenceEvent("2", 1, "", "A", "INSERTION"),
        reference,
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "REFERENCE_CONTIG_UNAVAILABLE"
