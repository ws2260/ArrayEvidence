"""Adversarial synthetic tests for independent bundle verification."""

import hashlib
import json

import pytest

from arrayevidence.bundle_verify import verify_bundle
from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AssayEvidenceRow,
)


def make_bundle(tmp_path, assay_rows=()):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, assay_rows, [])
    return bundle


def canonical_bytes(manifest):
    return (
        json.dumps(
            manifest,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        + "\n"
    ).encode("utf-8")


def change_manifest(bundle, modify):
    path = bundle / "bundle_manifest.json"
    manifest = json.loads(path.read_bytes())
    modify(manifest)
    path.write_bytes(canonical_bytes(manifest))


def replace_tsv_and_rehash(bundle, filename, data, row_count):
    (bundle / filename).write_bytes(data)

    def modify(manifest):
        entry = manifest["files"][filename]
        entry["sha256"] = hashlib.sha256(data).hexdigest()
        entry["size_bytes"] = len(data)
        entry["row_count"] = row_count

    change_manifest(bundle, modify)


def test_duplicate_json_keys_rejected(tmp_path):
    bundle = make_bundle(tmp_path)
    path = bundle / "bundle_manifest.json"

    original = path.read_text(encoding="utf-8")
    assert original.startswith('{"bundle_format_version":')

    modified = original.replace(
        '{"bundle_format_version":',
        '{"bundle_format_version":"1.0.0",'
        '"bundle_format_version":',
        1,
    )

    path.write_text(modified, encoding="utf-8")

    assert verify_bundle(bundle).valid is False


@pytest.mark.parametrize("field", [
    "bundle_format_version",
    "evidence_schema_version",
])
def test_boolean_version_rejected(tmp_path, field):
    bundle = make_bundle(tmp_path)

    change_manifest(
        bundle,
        lambda manifest: manifest.__setitem__(field, True),
    )

    assert verify_bundle(bundle).valid is False


@pytest.mark.parametrize("field", [
    "size_bytes",
    "row_count",
])
def test_negative_metadata_rejected(tmp_path, field):
    bundle = make_bundle(tmp_path)

    def modify(manifest):
        manifest["files"]["assay_evidence.tsv"][field] = -1

    change_manifest(bundle, modify)

    assert verify_bundle(bundle).valid is False


def test_embedded_nul_rejected_even_with_matching_hash(tmp_path):
    bundle = make_bundle(tmp_path)
    path = bundle / "assay_evidence.tsv"

    data = path.read_bytes()
    assert b"\t" in data

    corrupted = data.replace(b"\t", b"\x00", 1)

    replace_tsv_and_rehash(
        bundle,
        "assay_evidence.tsv",
        corrupted,
        row_count=0,
    )

    assert verify_bundle(bundle).valid is False


def test_blank_tsv_record_rejected_even_with_matching_hash(tmp_path):
    bundle = make_bundle(tmp_path)
    path = bundle / "assay_evidence.tsv"

    corrupted = path.read_bytes() + b"\n"

    replace_tsv_and_rehash(
        bundle,
        "assay_evidence.tsv",
        corrupted,
        row_count=1,
    )

    assert verify_bundle(bundle).valid is False


def test_quoted_tsv_field_rejected_even_with_matching_hash(tmp_path):
    bundle = make_bundle(tmp_path)
    path = bundle / "assay_evidence.tsv"

    data = path.read_bytes()
    corrupted = data.replace(
        b"schema_version",
        b'"schema_version"',
        1,
    )

    replace_tsv_and_rehash(
        bundle,
        "assay_evidence.tsv",
        corrupted,
        row_count=0,
    )

    assert verify_bundle(bundle).valid is False


def test_nonempty_synthetic_bundle_verifies(tmp_path):
    row = AssayEvidenceRow(
        schema_version=SCHEMA_VERSION,
        sample_id="synthetic_sample",
        assay_id="synthetic_assay",
        manifest_chrom="1",
        manifest_pos=5,
        coordinate_status="VCF_ABSENT",
        coordinate_evidence="ABSENT",
        coordinate_candidate_count=0,
        coordinate_genotype_evidence="NOT_EVALUATED",
        allele_evidence="NOT_EVALUATED",
        allele_reason="REFERENCE_NOT_PROVIDED",
        indel_evidence="NOT_EVALUATED",
        projection_status="NOT_EVALUATED",
        projection_reason="NOT_EVALUATED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )

    bundle = make_bundle(tmp_path, [row])

    result = verify_bundle(bundle)
    assert result.valid is True
    assert result.errors == ()

    manifest = json.loads(
        (bundle / "bundle_manifest.json").read_bytes()
    )

    assert manifest["files"]["assay_evidence.tsv"]["row_count"] == 1
    assert manifest["files"]["allele_candidates.tsv"]["row_count"] == 0


def test_corrupted_utf8_rejected_with_matching_hash(tmp_path):
    bundle = make_bundle(tmp_path)
    path = bundle / "assay_evidence.tsv"

    corrupted = path.read_bytes() + b"\xff\n"

    replace_tsv_and_rehash(
        bundle,
        "assay_evidence.tsv",
        corrupted,
        row_count=1,
    )

    assert verify_bundle(bundle).valid is False
