
"""Synthetic tests for the pysam FASTA reference adapter."""

from pathlib import Path

import pysam
import pytest

from arrayevidence.fasta_adapter import PysamReferenceReader


@pytest.fixture
def indexed_fasta(tmp_path: Path) -> Path:
    path = tmp_path / "synthetic.fa"

    path.write_text(
        ">1" + chr(10)
        + "ACGTACGTACGT" + chr(10)
        + ">X" + chr(10)
        + "TTGGAACC" + chr(10),
        encoding="ascii",
    )

    pysam.faidx(str(path))

    return path


def test_reference_fetch_uses_zero_based_half_open_coordinates(
    indexed_fasta,
):
    with PysamReferenceReader(indexed_fasta) as reference:
        assert reference.fetch("1", 0, 4) == "ACGT"
        assert reference.fetch("1", 4, 8) == "ACGT"
        assert reference.fetch("X", 2, 6) == "GGAA"


def test_reference_lengths(indexed_fasta):
    with PysamReferenceReader(indexed_fasta) as reference:
        assert reference.get_reference_length("1") == 12
        assert reference.get_reference_length("X") == 8


def test_chromosome_names_are_not_automatically_modified(
    indexed_fasta,
):
    with PysamReferenceReader(indexed_fasta) as reference:
        with pytest.raises(KeyError):
            reference.fetch("chr1", 0, 4)

        with pytest.raises(KeyError):
            reference.get_reference_length("chr1")


def test_invalid_coordinates_are_rejected(indexed_fasta):
    with PysamReferenceReader(indexed_fasta) as reference:
        with pytest.raises(ValueError):
            reference.fetch("1", -1, 4)

        with pytest.raises(ValueError):
            reference.fetch("1", 5, 4)

        with pytest.raises(ValueError):
            reference.fetch("1", 0, 13)


def test_missing_fasta_is_rejected(tmp_path):
    path = tmp_path / "missing.fa"

    with pytest.raises(FileNotFoundError):
        PysamReferenceReader(path)


def test_missing_fasta_index_is_rejected(tmp_path):
    path = tmp_path / "unindexed.fa"
    path.write_text(
        ">1" + chr(10) + "ACGT" + chr(10),
        encoding="ascii",
    )

    with pytest.raises(FileNotFoundError):
        PysamReferenceReader(path)


def test_context_manager_closes_reference(indexed_fasta):
    reference = PysamReferenceReader(indexed_fasta)

    with reference:
        assert reference.fetch("1", 0, 2) == "AC"

    with pytest.raises(ValueError, match="closed"):
        reference.fetch("1", 0, 2)


def test_explicit_close_is_idempotent(indexed_fasta):
    reference = PysamReferenceReader(indexed_fasta)

    reference.close()
    reference.close()

    with pytest.raises(ValueError, match="closed"):
        reference.get_reference_length("1")
