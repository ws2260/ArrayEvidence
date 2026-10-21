"""Safe publication requirements for complete evidence bundles."""

from pathlib import Path

import pytest

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AssayEvidenceRow,
)


def empty_evidence_export(destination):
    return export_evidence_tsv(destination, [], [])


def test_new_bundle_directory_is_published(tmp_path):
    destination = tmp_path / "evidence_bundle"

    paths = empty_evidence_export(destination)

    assert destination.is_dir()
    assert paths == (
        destination / "assay_evidence.tsv",
        destination / "allele_candidates.tsv",
    )
    assert all(path.is_file() for path in paths)


def test_existing_empty_directory_is_not_overwritten(tmp_path):
    destination = tmp_path / "existing"
    destination.mkdir()

    with pytest.raises(FileExistsError):
        empty_evidence_export(destination)

    assert list(destination.iterdir()) == []


def test_existing_bundle_is_not_overwritten(tmp_path):
    destination = tmp_path / "existing"
    destination.mkdir()

    original = destination / "assay_evidence.tsv"
    original.write_bytes(b"ORIGINAL_CONTENT\n")

    with pytest.raises(FileExistsError):
        empty_evidence_export(destination)

    assert original.read_bytes() == b"ORIGINAL_CONTENT\n"
    assert not (destination / "allele_candidates.tsv").exists()


def test_failed_validation_creates_no_bundle(tmp_path):
    destination = tmp_path / "invalid"

    with pytest.raises((TypeError, ValueError)):
        export_evidence_tsv(
            destination,
            ["not_an_assay_row"],
            [],
        )

    assert not destination.exists()
    assert list(tmp_path.iterdir()) == []


def test_write_failure_leaves_no_published_bundle(
    tmp_path, monkeypatch
):
    destination = tmp_path / "evidence_bundle"

    original_write_bytes = Path.write_bytes
    calls = 0

    def fail_second_write(self, data):
        nonlocal calls
        calls += 1

        if calls == 2:
            raise OSError("Synthetic second-file write failure")

        return original_write_bytes(self, data)

    monkeypatch.setattr(
        Path,
        "write_bytes",
        fail_second_write,
    )

    with pytest.raises(
        OSError,
        match="Synthetic second-file write failure",
    ):
        empty_evidence_export(destination)

    assert not destination.exists()
    assert list(tmp_path.iterdir()) == []


def test_bundle_contains_only_two_expected_files(tmp_path):
    destination = tmp_path / "bundle"

    empty_evidence_export(destination)

    assert sorted(path.name for path in destination.iterdir()) == [
        "allele_candidates.tsv",
        "assay_evidence.tsv",
    ]
