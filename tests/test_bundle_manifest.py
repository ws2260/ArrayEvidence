"""Synthetic contract tests for deterministic bundle metadata."""

import hashlib
import json
from pathlib import Path

import pytest

from arrayevidence.evidence_export import export_evidence_tsv


EXPECTED_FILES = {
    "assay_evidence.tsv",
    "allele_candidates.tsv",
}


def read_manifest(bundle):
    return json.loads(
        (bundle / "bundle_manifest.json").read_text(
            encoding="utf-8"
        )
    )


def test_export_publishes_three_complete_files(tmp_path):
    bundle = tmp_path / "bundle"

    export_evidence_tsv(bundle, [], [])

    assert {p.name for p in bundle.iterdir()} == (
        EXPECTED_FILES | {"bundle_manifest.json"}
    )


def test_manifest_has_exact_minimal_schema(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    manifest = read_manifest(bundle)

    assert set(manifest) == {
        "bundle_format_version",
        "evidence_schema_version",
        "files",
    }
    assert manifest["bundle_format_version"] == "1.0.0"
    assert manifest["evidence_schema_version"] == "1.0.0"
    assert set(manifest["files"]) == EXPECTED_FILES


def test_manifest_hashes_exact_published_bytes(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    manifest = read_manifest(bundle)

    for filename in EXPECTED_FILES:
        data = (bundle / filename).read_bytes()
        entry = manifest["files"][filename]

        assert set(entry) == {
            "sha256",
            "size_bytes",
            "row_count",
        }

        assert entry["sha256"] == hashlib.sha256(
            data
        ).hexdigest()

        assert entry["size_bytes"] == len(data)

        # Empty exports contain a header but zero data rows.
        assert entry["row_count"] == 0


def test_manifest_is_deterministic_across_destinations(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"

    export_evidence_tsv(first, [], [])
    export_evidence_tsv(second, [], [])

    assert (
        first / "bundle_manifest.json"
    ).read_bytes() == (
        second / "bundle_manifest.json"
    ).read_bytes()


def test_manifest_uses_canonical_json_and_newline(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    data = (bundle / "bundle_manifest.json").read_bytes()
    manifest = json.loads(data)

    expected = (
        json.dumps(
            manifest,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        + "\n"
    ).encode("utf-8")

    assert data == expected


def test_manifest_does_not_embed_identifiers_or_paths(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    data = (bundle / "bundle_manifest.json").read_text(
        encoding="utf-8"
    )

    assert str(tmp_path) not in data
    assert "sample_id" not in data
    assert "assay_id" not in data
    assert "created_at" not in data
    assert "reference_assembly" not in data


def test_manifest_is_included_before_publication(
    tmp_path, monkeypatch
):
    import arrayevidence.evidence_export as exporter

    original_publish = exporter._publish_directory_no_replace

    def inspect_staging(staging, destination):
        assert {p.name for p in staging.iterdir()} == {
            "assay_evidence.tsv",
            "allele_candidates.tsv",
            "bundle_manifest.json",
        }

        return original_publish(staging, destination)

    monkeypatch.setattr(
        exporter,
        "_publish_directory_no_replace",
        inspect_staging,
    )

    export_evidence_tsv(tmp_path / "bundle", [], [])


def test_existing_destination_still_cannot_be_overwritten(
    tmp_path,
):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    sentinel = bundle / "sentinel.txt"
    sentinel.write_text("KEEP", encoding="utf-8")

    with pytest.raises(FileExistsError):
        export_evidence_tsv(bundle, [], [])

    assert sentinel.read_text(encoding="utf-8") == "KEEP"
    assert sorted(p.name for p in bundle.iterdir()) == [
        "sentinel.txt"
    ]
