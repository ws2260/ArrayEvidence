"""Optional local integration test using authorized AIDA inputs."""

from collections import Counter
from pathlib import Path

import pytest

from arrayevidence.manifest import read_assays
from arrayevidence.vcf import read_vcf
from arrayevidence.reconcile import reconcile


@pytest.mark.integration
def test_aida_coordinate_baseline():
    root = Path.home() / "research-data/AIDA"

    manifest = (
        root / "resources/illumina_gsa_v3/csv"
        / "GSA-24v3-0_A2.csv"
    )
    vcf_dir = root / "processed/vcf/dragen_array_v1.4"

    if not manifest.exists() or not vcf_dir.exists():
        pytest.skip("AIDA inputs unavailable")

    vcfs = sorted(vcf_dir.glob("*.vcf.gz"))
    if len(vcfs) != 48:
        pytest.skip("Expected AIDA VCF cohort unavailable")

    results = reconcile(
        read_assays(manifest),
        read_vcf(vcfs[0]),
    )

    observed = Counter(r.status for r in results)

    assert observed == {
        "UNMAPPED_ASSAY": 210,
        "VCF_ABSENT": 11150,
        "VCF_CALLED": 638148,
        "VCF_NO_CALL": 4519,
    }
