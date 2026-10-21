"""Failure handling and publication-boundary tests."""

import errno
from pathlib import Path

import pytest

import arrayevidence.evidence_export as exporter


def test_destination_created_during_publication_is_preserved(
    tmp_path, monkeypatch
):
    destination = tmp_path / "bundle"
    original_publish = exporter._publish_directory_no_replace

    def competing_publication(staging, target):
        assert target == destination
        assert staging.is_dir()
        assert sorted(p.name for p in staging.iterdir()) == [
            "allele_candidates.tsv",
            "assay_evidence.tsv",
        ]

        # Simulate another process creating the destination
        # after the initial existence check.
        target.mkdir()
        (target / "other_process.txt").write_text(
            "DO NOT MODIFY",
            encoding="utf-8",
        )

        return original_publish(staging, target)

    monkeypatch.setattr(
        exporter,
        "_publish_directory_no_replace",
        competing_publication,
    )

    with pytest.raises(FileExistsError):
        exporter.export_evidence_tsv(destination, [], [])

    assert destination.is_dir()
    assert (
        destination / "other_process.txt"
    ).read_text(encoding="utf-8") == "DO NOT MODIFY"

    assert sorted(p.name for p in destination.iterdir()) == [
        "other_process.txt",
    ]

    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "bundle",
    ]


def test_unexpected_publication_failure_cleans_staging(
    tmp_path, monkeypatch
):
    destination = tmp_path / "bundle"

    def failed_publish(staging, target):
        assert staging.is_dir()
        assert not target.exists()
        raise OSError(
            errno.EIO,
            "Synthetic publication I/O failure",
        )

    monkeypatch.setattr(
        exporter,
        "_publish_directory_no_replace",
        failed_publish,
    )

    with pytest.raises(OSError, match="Synthetic publication"):
        exporter.export_evidence_tsv(destination, [], [])

    assert not destination.exists()
    assert list(tmp_path.iterdir()) == []


def test_missing_parent_is_not_created(tmp_path):
    destination = tmp_path / "missing_parent" / "bundle"

    with pytest.raises(FileNotFoundError):
        exporter.export_evidence_tsv(destination, [], [])

    assert not destination.exists()
    assert not destination.parent.exists()


def test_existing_file_destination_is_preserved(tmp_path):
    destination = tmp_path / "bundle"
    destination.write_bytes(b"EXISTING_FILE\n")

    with pytest.raises(FileExistsError):
        exporter.export_evidence_tsv(destination, [], [])

    assert destination.read_bytes() == b"EXISTING_FILE\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_existing_symlink_destination_is_preserved(tmp_path):
    target = tmp_path / "target"
    target.mkdir()

    destination = tmp_path / "bundle"
    destination.symlink_to(target, target_is_directory=True)

    with pytest.raises(FileExistsError):
        exporter.export_evidence_tsv(destination, [], [])

    assert destination.is_symlink()
    assert destination.resolve() == target.resolve()
    assert list(target.iterdir()) == []


def test_successful_publication_leaves_no_staging_directory(
    tmp_path,
):
    destination = tmp_path / "bundle"

    exporter.export_evidence_tsv(destination, [], [])

    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "bundle",
    ]
    assert sorted(p.name for p in destination.iterdir()) == [
        "allele_candidates.tsv",
        "assay_evidence.tsv",
    ]
