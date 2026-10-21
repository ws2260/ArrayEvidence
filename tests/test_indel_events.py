import pytest

from arrayevidence.indel_events import minimal_sequence_event


@pytest.mark.parametrize(
    "pos,ref,alt,start,deleted,inserted,status",
    [
        (100, "A", "AG", 100, "", "G", "INSERTION"),
        (100, "AC", "A", 100, "C", "", "DELETION"),
        (100, "AC", "ACG", 101, "", "G", "INSERTION"),
        (100, "ACG", "AC", 101, "G", "", "DELETION"),
        (100, "ACG", "ATG", 100, "C", "T", "SUBSTITUTION"),
        (100, "AC", "GT", 99, "AC", "GT", "SUBSTITUTION"),
        (100, "AC", "ATG", 100, "C", "TG", "COMPLEX_INDEL"),
        (100, "A", "A", 100, "", "", "NO_CHANGE"),
    ],
)
def test_minimal_events(
    pos, ref, alt, start, deleted, inserted, status
):
    result = minimal_sequence_event(
        "1", pos, ref, alt
    )

    assert result.start == start
    assert result.deleted == deleted
    assert result.inserted == inserted
    assert result.status == status


def test_length_change():
    insertion = minimal_sequence_event(
        "1", 100, "A", "ATTT"
    )

    deletion = minimal_sequence_event(
        "1", 100, "ATTT", "A"
    )

    assert insertion.length_change == 3
    assert deletion.length_change == -3


@pytest.mark.parametrize(
    "pos,ref,alt",
    [
        (0, "A", "AT"),
        (100, "A", "<DEL>"),
        (100, "A", "*"),
        (100, "N", "A"),
        (100, "", "A"),
    ],
)
def test_unsupported_alleles(pos, ref, alt):
    result = minimal_sequence_event(
        "1", pos, ref, alt
    )

    assert result.status == "UNRESOLVED"


def test_shared_suffix_trimming():
    result = minimal_sequence_event(
        "1", 100, "ATC", "AGC"
    )

    assert result.start == 100
    assert result.deleted == "T"
    assert result.inserted == "G"


def test_different_anchoring_can_produce_same_event():
    first = minimal_sequence_event(
        "1", 100, "A", "AG"
    )

    second = minimal_sequence_event(
        "1", 99, "CA", "CAG"
    )

    assert first == second


def test_repeat_equivalence_not_assumed():
    first = minimal_sequence_event(
        "1", 100, "A", "AA"
    )

    second = minimal_sequence_event(
        "1", 101, "A", "AA"
    )

    assert first != second
