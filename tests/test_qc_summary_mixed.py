"""Synthetic mixed-category QC tests."""

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.qc import summarize_bundle

from test_evidence_export import assay_row, candidate_row


ASSAY_FIELDS = (
    "coordinate_status",
    "coordinate_genotype_evidence",
    "allele_evidence",
    "allele_reason",
    "indel_evidence",
    "projection_status",
    "projection_reason",
)


def test_mixed_qc_counts(tmp_path):
    first = assay_row(
        sample_id="SYNTHETIC",
        assay_id="ASSAY_001",
    )

    second = assay_row(
        sample_id="SYNTHETIC",
        assay_id="ASSAY_002",
    )

    first_candidate = candidate_row(
        sample_id="SYNTHETIC",
        assay_id="ASSAY_001",
    )

    second_candidate = candidate_row(
        sample_id="SYNTHETIC",
        assay_id="ASSAY_002",
    )

    bundle = tmp_path / "bundle"

    export_evidence_tsv(
        bundle,
        [first, second],
        [first_candidate, second_candidate],
    )

    summary = summarize_bundle(bundle)

    assert summary["assay_rows"] == 2
    assert summary["candidate_rows"] == 2

    for field in ASSAY_FIELDS:
        assert sum(summary[field].values()) == 2

    assert sum(summary["candidate_status"].values()) == 2


def test_qc_category_totals_match_denominators(tmp_path):
    row = assay_row(
        sample_id="SYNTHETIC",
        assay_id="ASSAY_001",
    )

    candidate = candidate_row(
        sample_id="SYNTHETIC",
        assay_id="ASSAY_001",
    )

    bundle = tmp_path / "bundle"

    export_evidence_tsv(
        bundle,
        [row],
        [candidate],
    )

    summary = summarize_bundle(bundle)

    for field in ASSAY_FIELDS:
        assert (
            sum(summary[field].values())
            == summary["assay_rows"]
        )

    assert (
        sum(summary["candidate_status"].values())
        == summary["candidate_rows"]
    )
