"""Synthetic nonempty QC tests, independent of participant data."""

from pathlib import Path

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import AssayEvidenceRow
from arrayevidence.qc import summarize_bundle


def test_nonempty_qc_counts_synthetic_assay(tmp_path):
    row = AssayEvidenceRow(
        schema_version="1.0.0",
        sample_id="SYNTHETIC",
        assay_id="SYNTHETIC_ASSAY",
        manifest_chrom="1",
        manifest_pos=100,
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
    export_evidence_tsv(bundle, [row], [])

    before = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    summary = summarize_bundle(bundle)

    assert summary["assay_rows"] == 1
    assert summary["candidate_rows"] == 0
    assert summary["coordinate_status"] == {"VCF_ABSENT": 1}
    assert summary["allele_evidence"] == {"NOT_EVALUATED": 1}
    assert summary["projection_status"] == {"NOT_EVALUATED": 1}
    assert summary["candidate_status"] == {}

    after = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    assert after == before
