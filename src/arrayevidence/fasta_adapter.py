"""Read-only indexed FASTA adapter for reference-based evidence."""

from pathlib import Path

import pysam


class PysamReferenceReader:
    """Access an indexed FASTA using zero-based, half-open coordinates.

    Chromosome names must match the FASTA index exactly. An existing
    .fai index is required; this class never creates or updates one.
    """

    def __init__(self, path: str | Path) -> None:
        fasta_path = Path(path)

        if not fasta_path.is_file():
            raise FileNotFoundError(
                "Reference FASTA file does not exist"
            )

        index_path = Path(str(fasta_path) + ".fai")

        if not index_path.is_file():
            raise FileNotFoundError(
                "Reference FASTA index (.fai) does not exist"
            )

        self._fasta = pysam.FastaFile(str(fasta_path))
        self._closed = False
        self._references = frozenset(self._fasta.references)

    def _require_open(self) -> None:
        if self._closed:
            raise ValueError("Reference FASTA reader is closed")

    def _require_chrom(self, chrom: str) -> None:
        if chrom not in self._references:
            raise KeyError(
                "Chromosome is absent from the reference FASTA"
            )

    def get_reference_length(self, chrom: str) -> int:
        """Return the length of a reference sequence."""

        self._require_open()
        self._require_chrom(chrom)

        return self._fasta.get_reference_length(chrom)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch [start, end) using zero-based coordinates."""

        self._require_open()
        self._require_chrom(chrom)

        if not isinstance(start, int) or not isinstance(end, int):
            raise TypeError(
                "Reference coordinates must be integers"
            )

        length = self._fasta.get_reference_length(chrom)

        if start < 0 or end < start or end > length:
            raise ValueError(
                "Reference coordinates are outside valid bounds"
            )

        return self._fasta.fetch(chrom, start, end)

    def close(self) -> None:
        """Close the reference handle; safe to call repeatedly."""

        if not self._closed:
            self._fasta.close()
            self._closed = True

    def __enter__(self) -> "PysamReferenceReader":
        self._require_open()
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()
