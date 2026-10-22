"""Synthetic mixed-category QC validation."""

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import AssayEvidenceRow
from arrayevidence.qc import summarize_bundle, render_qc_json

from test_evidence_export import assay_row, candidate_row


def test_distinct_evidence_categories_and_denominators(tmp_path):
    matched = assay_row(
        sample_id="SYNTHETIC",
        assay_id="PRIVATE_SYNTHETIC_ASSAY_73921",
    )

    candidate = candidate_row(
        sample_id="SYNTHETIC",
        assay_id="PRIVATE_SYNTHETIC_ASSAY_73921",
    )

    absent = AssayEvidenceRow(
        schema_version="1.0.0",
        sample_id="SYNTHETIC",
        assay_id="PRIVATE_SYNTHETIC_ASSAY_84620",
        manifest_chrom="1",
        manifest_pos=200,
        coordinate_status="VCF_ABSENT",
        coordinate_evidence="VCF_ABSENT",
        coordinate_candidate_count=0,
        coordinate_genotype_evidence="NOT_EVALUATED",
        allele_evidence="NOT_EVALUATED",
        allele_reason="VCF_ABSENT",
        indel_evidence="NOT_EVALUATED",
        projection_status="NOT_EVALUATED",
        projection_reason="NOT_EVALUATED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )

    bundle = tmp_path / "bundle"

    export_evidence_tsv(
        bundle,
        [matched, absent],
        [candidate],
    )

    before = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    summary = summarize_bundle(bundle)

    assert summary["assay_rows"] == 2
    assert summary["candidate_rows"] == 1

    assert summary["allele_evidence"] == {
        "MATCH": 1,
        "NOT_EVALUATED": 1,
    }

    assert summary["indel_evidence"] == {
        "NOT_EVALUATED": 1,
        "UNIQUE": 1,
    }

    assert summary["coordinate_status"]["VCF_ABSENT"] >= 1
    assert summary["allele_reason"]["VCF_ABSENT"] == 1

    for field in (
        "coordinate_status",
        "coordinate_genotype_evidence",
        "allele_evidence",
        "allele_reason",
        "indel_evidence",
        "projection_status",
        "projection_reason",
    ):
        assert sum(summary[field].values()) == 2

    assert sum(summary["candidate_status"].values()) == 1

    rendered = render_qc_json(summary)

    assert "SYNTHETIC" not in rendered
    assert "PRIVATE_SYNTHETIC_ASSAY_73921" not in rendered
    assert "PRIVATE_SYNTHETIC_ASSAY_84620" not in rendered
    assert "sample_id" not in rendered
    assert "assay_id" not in rendered
    assert "vcf_gt" not in rendered

    after = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    assert after == before
