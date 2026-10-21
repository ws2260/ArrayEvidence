"""Targeted synthetic tests for strict bundle parsing."""

import hashlib
import json

from arrayevidence.bundle_verify import verify_bundle
from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AssayEvidenceRow,
)


def make_nonempty_bundle(tmp_path):
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

    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [row], [])
    return bundle


def update_integrity_metadata(bundle, filename):
    path = bundle / filename
    data = path.read_bytes()

    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_bytes())

    entry = manifest["files"][filename]
    entry["sha256"] = hashlib.sha256(data).hexdigest()
    entry["size_bytes"] = len(data)

    manifest_path.write_bytes(
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


def test_duplicate_json_key_is_rejected_by_parser(tmp_path):
    """A strict parser must reject duplicate keys, not silently merge."""
    import arrayevidence.bundle_verify as verifier

    assert hasattr(verifier, "_parse_manifest_strict")

    duplicate = (
        b'{"bundle_format_version":"1.0.0",'
        b'"bundle_format_version":"1.0.0"}'
    )

    try:
        verifier._parse_manifest_strict(duplicate)
    except ValueError:
        pass
    else:
        raise AssertionError("Duplicate JSON keys were accepted")


def test_double_quote_inside_data_field_is_rejected(tmp_path):
    bundle = make_nonempty_bundle(tmp_path)
    path = bundle / "assay_evidence.tsv"

    data = path.read_bytes()
    assert b"synthetic_sample" in data

    path.write_bytes(
        data.replace(
            b"synthetic_sample",
            b'synthetic_"sample',
            1,
        )
    )

    update_integrity_metadata(bundle, "assay_evidence.tsv")

    result = verify_bundle(bundle)
    assert result.valid is False


def test_nul_inside_data_field_is_rejected(tmp_path):
    bundle = make_nonempty_bundle(tmp_path)
    path = bundle / "assay_evidence.tsv"

    original = path.read_bytes()
    assert b"synthetic_sample" in original

    # Construct a real NUL byte without ambiguous escape sequences.
    corrupted = original.replace(
        b"synthetic_sample",
        b"synthetic_" + bytes([0]) + b"sample",
        1,
    )

    assert corrupted != original
    assert bytes([0]) in corrupted

    path.write_bytes(corrupted)

    # Recompute checksums and byte sizes so that the verifier
    # must detect the invalid TSV content independently.
    update_integrity_metadata(
        bundle,
        "assay_evidence.tsv",
    )

    # Confirm the on-disk fixture actually contains a NUL byte.
    actual = path.read_bytes()
    assert bytes([0]) in actual

    result = verify_bundle(bundle)

    assert result.valid is False
    assert any(
        "TSV structure" in error
        for error in result.errors
    )
