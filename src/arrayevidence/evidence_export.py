"""Deterministic TSV serialization of versioned technical evidence.

This module does not perform genotype calling, allele matching,
clinical annotation, de-identification, or access authorization.

Callers are responsible for ensuring that the supplied sample
identifiers and evidence are authorized for the intended output.
"""

import csv
import ctypes
import errno
import io
import os
import shutil
import tempfile
from dataclasses import fields
from pathlib import Path
from typing import Iterable

from .evidence_schema import AlleleCandidateRow, AssayEvidenceRow


NULL_TOKEN = "NULL"

ASSAY_FILENAME = "assay_evidence.tsv"
CANDIDATE_FILENAME = "allele_candidates.tsv"


def _serialize_value(value: object) -> str:
    """Serialize one scalar without conflating missingness and zero."""

    if value is None:
        return NULL_TOKEN

    if isinstance(value, str):
        if value == NULL_TOKEN:
            raise ValueError(
                "Literal 'NULL' is reserved for missing values"
            )

        if any(character in value for character in "\t\r\n"):
            raise ValueError(
                "TSV values must not contain tabs or line breaks"
            )

        return value

    if type(value) in {int, bool}:
        return str(value)

    raise TypeError(
        f"Unsupported evidence value type: {type(value).__name__}"
    )


def _prepare_table(rows, row_type):
    """Validate and serialize all rows before any file is written."""

    names = tuple(field.name for field in fields(row_type))
    prepared = []

    for row in rows:
        if not isinstance(row, row_type):
            raise TypeError(
                f"Expected {row_type.__name__}, "
                f"received {type(row).__name__}"
            )

        prepared.append(
            tuple(
                _serialize_value(getattr(row, name))
                for name in names
            )
        )

    return names, prepared


def _validate_table_integrity(
    assay_rows: list[AssayEvidenceRow],
    candidate_rows: list[AlleleCandidateRow],
) -> None:
    """Validate assay keys and candidate relationships before writing.

    Assay key:
        (sample_id, assay_id)

    Candidate key:
        (sample_id, assay_id, candidate_index)

    Candidate indices must be consecutive and 1-based within each
    assay. Assays without candidate rows are permitted.
    """

    assay_keys: set[tuple[str, str]] = set()

    for row in assay_rows:
        key = (row.sample_id, row.assay_id)

        if key in assay_keys:
            raise ValueError(
                f"Duplicate assay evidence key: {key!r}"
            )

        assay_keys.add(key)

    candidate_indices: dict[
        tuple[str, str], set[int]
    ] = {}

    for row in candidate_rows:
        key = (row.sample_id, row.assay_id)

        if key not in assay_keys:
            raise ValueError(
                f"Candidate references missing assay evidence: {key!r}"
            )

        indices = candidate_indices.setdefault(key, set())

        if row.candidate_index in indices:
            raise ValueError(
                f"Duplicate candidate index {row.candidate_index} "
                f"for assay key {key!r}"
            )

        indices.add(row.candidate_index)

    for key, indices in candidate_indices.items():
        expected = set(range(1, len(indices) + 1))

        if indices != expected:
            raise ValueError(
                f"Nonconsecutive candidate indices for "
                f"assay key {key!r}: "
                f"observed={sorted(indices)!r}, "
                f"expected={sorted(expected)!r}"
            )


def _validate_evidence_consistency(
    assay_rows: list[AssayEvidenceRow],
    candidate_rows: list[AlleleCandidateRow],
) -> None:
    """Check that candidate provenance agrees with parent evidence.

    This validates the current reconciliation evidence contract.
    It does not independently establish allele equivalence or
    genotype correctness.
    """

    candidates_by_assay: dict[
        tuple[str, str], list[AlleleCandidateRow]
    ] = {}

    for candidate in candidate_rows:
        key = (candidate.sample_id, candidate.assay_id)
        candidates_by_assay.setdefault(key, []).append(candidate)

    for assay in assay_rows:
        key = (assay.sample_id, assay.assay_id)
        candidates = candidates_by_assay.get(key, [])

        indel_status = assay.indel_evidence
        allele_status = assay.allele_evidence

        if indel_status == "UNIQUE":
            expected_count = 1
            expected_status = "UNIQUE"

        elif indel_status == "MULTIPLE":
            if len(candidates) < 2:
                raise ValueError(
                    f"MULTIPLE indel evidence requires at least "
                    f"two candidates for {key!r}"
                )

            if assay.projection_status not in {
                "UNRESOLVED",
                "NOT_EVALUATED",
            }:
                raise ValueError(
                    f"MULTIPLE indel evidence cannot have "
                    f"resolved projection for {key!r}"
                )

            if allele_status != "UNRESOLVED":
                raise ValueError(
                    f"MULTIPLE indel evidence requires "
                    f"unresolved allele evidence for {key!r}"
                )

            expected_count = len(candidates)
            expected_status = "MULTIPLE"

        elif indel_status in {"ABSENT", "UNRESOLVED"}:
            expected_count = 0
            expected_status = None

        elif indel_status == "NOT_EVALUATED":
            if allele_status in {"MATCH", "PARTIAL_MATCH"}:
                expected_count = 1
                expected_status = allele_status
            else:
                expected_count = 0
                expected_status = None

        else:
            raise ValueError(
                f"Unsupported indel evidence status "
                f"{indel_status!r} for {key!r}"
            )

        if len(candidates) != expected_count:
            raise ValueError(
                f"Evidence candidate count mismatch for {key!r}: "
                f"indel_evidence={indel_status!r}, "
                f"allele_evidence={allele_status!r}, "
                f"expected={expected_count}, "
                f"observed={len(candidates)}"
            )

        if expected_status is not None:
            for candidate in candidates:
                if candidate.candidate_status != expected_status:
                    raise ValueError(
                        f"Candidate status mismatch for {key!r}: "
                        f"expected={expected_status!r}, "
                        f"observed={candidate.candidate_status!r}"
                    )


def _render_tsv(header, rows):
    """Render TSV entirely in memory using Unix line endings."""

    buffer = io.StringIO(newline="")
    writer = csv.writer(
        buffer,
        delimiter="\t",
        lineterminator="\n",
        quoting=csv.QUOTE_MINIMAL,
    )

    writer.writerow(header)
    writer.writerows(rows)

    return buffer.getvalue().encode("utf-8")


def _publish_directory_no_replace(
    staging: Path,
    destination: Path,
) -> None:
    """Atomically publish a new directory without overwriting.

    Linux renameat2 with RENAME_NOREPLACE ensures that another
    process cannot create the destination between checking and
    publication.

    This implementation intentionally requires Linux/WSL.
    """

    if os.name != "posix":
        raise OSError(
            "Atomic no-replace bundle publication requires Linux"
        )

    libc = ctypes.CDLL(None, use_errno=True)

    renameat2 = getattr(libc, "renameat2", None)

    if renameat2 is None:
        raise OSError(
            "renameat2 is unavailable; refusing unsafe publication"
        )

    renameat2.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    renameat2.restype = ctypes.c_int

    AT_FDCWD = -100
    RENAME_NOREPLACE = 1

    result = renameat2(
        AT_FDCWD,
        os.fsencode(staging),
        AT_FDCWD,
        os.fsencode(destination),
        RENAME_NOREPLACE,
    )

    if result != 0:
        error = ctypes.get_errno()

        if error == errno.EEXIST:
            raise FileExistsError(
                error,
                "Evidence bundle destination already exists",
                str(destination),
            )

        raise OSError(
            error,
            os.strerror(error),
            str(destination),
        )


def export_evidence_tsv(
    output_dir: str | Path,
    assay_rows: Iterable[AssayEvidenceRow],
    candidate_rows: Iterable[AlleleCandidateRow],
) -> tuple[Path, Path]:
    """Write deterministic assay and allele-candidate TSV tables.

    Sorting keys:
      Assays: sample_id, assay_id, then complete serialized row.
      Candidates: sample_id, assay_id, candidate_index,
                  then complete serialized row.

    The complete serialized row breaks ties deterministically.

    All input rows are validated and rendered before staging.
    A completed directory is published without overwriting an
    existing destination. Requires Linux/WSL renameat2 support.

    Publication is atomic on a supported local filesystem, but
    power-loss durability is not guaranteed.
    """

    # Materialize and type-check inputs before relational checks.
    assay_rows = list(assay_rows)
    candidate_rows = list(candidate_rows)

    assay_header, assays = _prepare_table(
        assay_rows,
        AssayEvidenceRow,
    )

    candidate_header, candidates = _prepare_table(
        candidate_rows,
        AlleleCandidateRow,
    )

    # All consistency checks precede output directory creation.
    _validate_table_integrity(assay_rows, candidate_rows)
    _validate_evidence_consistency(assay_rows, candidate_rows)

    assay_columns = {
        name: index
        for index, name in enumerate(assay_header)
    }
    candidate_columns = {
        name: index
        for index, name in enumerate(candidate_header)
    }

    assays.sort(
        key=lambda row: (
            row[assay_columns["sample_id"]],
            row[assay_columns["assay_id"]],
            row,
        )
    )

    candidates.sort(
        key=lambda row: (
            row[candidate_columns["sample_id"]],
            row[candidate_columns["assay_id"]],
            int(row[candidate_columns["candidate_index"]]),
            row,
        )
    )

    assay_data = _render_tsv(assay_header, assays)
    candidate_data = _render_tsv(
        candidate_header,
        candidates,
    )

    destination = Path(output_dir)

    if destination.exists() or destination.is_symlink():
        raise FileExistsError(
            f"Evidence bundle destination already exists: {destination}"
        )

    parent = destination.parent

    if not parent.is_dir():
        raise FileNotFoundError(
            f"Evidence bundle parent directory does not exist: {parent}"
        )

    # Stage beside the destination to keep the rename on the
    # same filesystem. Never write directly into the final bundle.
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{destination.name}.staging-",
            dir=parent,
        )
    )

    try:
        (staging / ASSAY_FILENAME).write_bytes(assay_data)
        (staging / CANDIDATE_FILENAME).write_bytes(candidate_data)

        # Publish both completed files as one visible directory.
        _publish_directory_no_replace(staging, destination)

    finally:
        # After successful publication staging no longer exists.
        # After failure, remove only the staging directory.
        if staging.exists():
            shutil.rmtree(staging)

    return (
        destination / ASSAY_FILENAME,
        destination / CANDIDATE_FILENAME,
    )
