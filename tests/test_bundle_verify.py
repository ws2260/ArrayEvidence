"""Synthetic contract tests for independent bundle verification."""

import hashlib
import json

import pytest

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.bundle_verify import verify_bundle


def make_bundle(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])
    return bundle


def rewrite_manifest(bundle, transform):
    path = bundle / "bundle_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    transform(manifest)
    path.write_bytes(
        (
            json.dumps(
                manifest,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
            )
            + "\n"
        ).encode("utf-8")
    )


def test_valid_bundle_passes(tmp_path):
    result = verify_bundle(make_bundle(tmp_path))
    assert result.valid is True
    assert result.errors == ()


def test_verification_does_not_modify_files(tmp_path):
    bundle = make_bundle(tmp_path)
    before = {
        p.name: p.read_bytes()
        for p in bundle.iterdir()
    }

    verify_bundle(bundle)

    after = {
        p.name: p.read_bytes()
        for p in bundle.iterdir()
    }
    assert before == after


def test_missing_bundle_fails(tmp_path):
    result = verify_bundle(tmp_path / "missing")
    assert result.valid is False
    assert result.errors


def test_missing_data_file_fails(tmp_path):
    bundle = make_bundle(tmp_path)
    (bundle / "assay_evidence.tsv").unlink()

    result = verify_bundle(bundle)
    assert result.valid is False


def test_extra_file_fails(tmp_path):
    bundle = make_bundle(tmp_path)
    (bundle / "unexpected.txt").write_text("x")

    result = verify_bundle(bundle)
    assert result.valid is False


def test_modified_tsv_fails(tmp_path):
    bundle = make_bundle(tmp_path)
    with (bundle / "assay_evidence.tsv").open("ab") as handle:
        handle.write(b"tampered\n")

    result = verify_bundle(bundle)
    assert result.valid is False


def test_wrong_hash_fails(tmp_path):
    bundle = make_bundle(tmp_path)

    def change(manifest):
        manifest["files"]["assay_evidence.tsv"]["sha256"] = "0" * 64

    rewrite_manifest(bundle, change)
    assert verify_bundle(bundle).valid is False


def test_wrong_byte_size_fails(tmp_path):
    bundle = make_bundle(tmp_path)

    def change(manifest):
        manifest["files"]["assay_evidence.tsv"]["size_bytes"] += 1

    rewrite_manifest(bundle, change)
    assert verify_bundle(bundle).valid is False


def test_wrong_row_count_fails(tmp_path):
    bundle = make_bundle(tmp_path)

    def change(manifest):
        manifest["files"]["assay_evidence.tsv"]["row_count"] += 1

    rewrite_manifest(bundle, change)
    assert verify_bundle(bundle).valid is False


@pytest.mark.parametrize("field", [
    "bundle_format_version",
    "evidence_schema_version",
])
def test_unsupported_version_fails(tmp_path, field):
    bundle = make_bundle(tmp_path)

    def change(manifest):
        manifest[field] = "999.0.0"

    rewrite_manifest(bundle, change)
    assert verify_bundle(bundle).valid is False


def test_unexpected_manifest_field_fails(tmp_path):
    bundle = make_bundle(tmp_path)

    def change(manifest):
        manifest["participant_id"] = "not_allowed"

    rewrite_manifest(bundle, change)
    assert verify_bundle(bundle).valid is False


def test_noncanonical_json_fails(tmp_path):
    bundle = make_bundle(tmp_path)
    path = bundle / "bundle_manifest.json"
    manifest = json.loads(path.read_bytes())
    path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    assert verify_bundle(bundle).valid is False


def test_malformed_json_fails(tmp_path):
    bundle = make_bundle(tmp_path)
    (bundle / "bundle_manifest.json").write_text(
        "{invalid",
        encoding="utf-8",
    )

    assert verify_bundle(bundle).valid is False


def test_symlinked_data_file_fails(tmp_path):
    bundle = make_bundle(tmp_path)
    path = bundle / "assay_evidence.tsv"
    target = tmp_path / "outside.tsv"
    target.write_bytes(path.read_bytes())

    path.unlink()
    path.symlink_to(target)

    assert verify_bundle(bundle).valid is False


def test_symlinked_bundle_directory_fails(tmp_path):
    bundle = make_bundle(tmp_path)
    alias = tmp_path / "alias"
    alias.symlink_to(bundle, target_is_directory=True)

    assert verify_bundle(alias).valid is False


def test_invalid_metadata_type_fails(tmp_path):
    bundle = make_bundle(tmp_path)

    def change(manifest):
        manifest["files"]["assay_evidence.tsv"]["size_bytes"] = True

    rewrite_manifest(bundle, change)
    assert verify_bundle(bundle).valid is False


def test_manifest_filename_set_is_exact(tmp_path):
    bundle = make_bundle(tmp_path)

    def change(manifest):
        manifest["files"]["extra.tsv"] = {
            "sha256": hashlib.sha256(b"").hexdigest(),
            "size_bytes": 0,
            "row_count": 0,
        }

    rewrite_manifest(bundle, change)
    assert verify_bundle(bundle).valid is False
