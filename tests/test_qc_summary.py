import json
from pathlib import Path

import pytest

from arrayevidence.qc import summarize_bundle, render_qc_json


ASSAY_HEADER = (
    "schema_version\tsample_id\tassay_id\tmanifest_chrom\t"
    "manifest_pos\tcoordinate_status\tcoordinate_evidence\t"
    "coordinate_candidate_count\tcoordinate_genotype_evidence\t"
    "allele_evidence\tallele_reason\tindel_evidence\t"
    "projection_status\tprojection_reason\tassay_alt_dosage\t"
    "other_alt_present\tploidy"
)

CANDIDATE_HEADER = (
    "schema_version\tsample_id\tassay_id\tcandidate_index\t"
    "vcf_chrom\tvcf_pos\tvcf_id\tvcf_ref\tvcf_alt\t"
    "vcf_alt_index\tvcf_gt\tcandidate_status"
)


def test_qc_summary_requires_verified_bundle(tmp_path):
    with pytest.raises(ValueError):
        summarize_bundle(tmp_path / "nonexistent")


def test_qc_json_is_deterministic():
    summary = {
        "assay_rows": 2,
        "candidate_rows": 1,
        "coordinate_status": {
            "VCF_CALLED": 1,
            "VCF_ABSENT": 1,
        },
    }

    first = render_qc_json(summary)
    second = render_qc_json(summary)

    assert first == second
    assert first.endswith("\n")
    assert json.loads(first) == summary


def test_qc_json_uses_canonical_key_order():
    summary = {
        "z": 1,
        "a": {"y": 2, "b": 3},
    }

    rendered = render_qc_json(summary)

    assert rendered == (
        '{"a":{"b":3,"y":2},"z":1}\n'
    )


def test_qc_json_rejects_nonfinite_values():
    with pytest.raises(ValueError):
        render_qc_json({"rate": float("nan")})


def test_qc_summary_does_not_accept_unverified_tsv_directory(tmp_path):
    bundle = tmp_path / "unverified"
    bundle.mkdir()

    (bundle / "assay_evidence.tsv").write_text(
        ASSAY_HEADER + "\n",
        encoding="utf-8",
    )

    (bundle / "allele_candidates.tsv").write_text(
        CANDIDATE_HEADER + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError):
        summarize_bundle(bundle)


def test_qc_summary_rejects_modified_verified_bundle(
    tmp_path,
    monkeypatch,
):
    from arrayevidence.evidence_export import export_evidence_tsv

    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    assay_file = bundle / "assay_evidence.tsv"

    with assay_file.open("a", encoding="utf-8") as handle:
        handle.write("\n")

    with pytest.raises(ValueError):
        summarize_bundle(bundle)


def test_empty_verified_bundle_has_zero_counts(tmp_path):
    from arrayevidence.evidence_export import export_evidence_tsv

    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    summary = summarize_bundle(bundle)

    assert summary["assay_rows"] == 0
    assert summary["candidate_rows"] == 0
    assert summary["coordinate_status"] == {}
    assert summary["allele_evidence"] == {}
    assert summary["projection_status"] == {}


def test_qc_summary_excludes_sensitive_metadata(tmp_path):
    from arrayevidence.evidence_export import export_evidence_tsv

    bundle = tmp_path / "private_sample_identifier"
    export_evidence_tsv(bundle, [], [])

    summary = summarize_bundle(bundle)
    rendered = render_qc_json(summary)

    assert "private_sample_identifier" not in rendered
    assert "sample_id" not in rendered
    assert "assay_id" not in rendered
    assert "vcf_gt" not in rendered
    assert "vcf_pos" not in rendered
    assert str(tmp_path) not in rendered
