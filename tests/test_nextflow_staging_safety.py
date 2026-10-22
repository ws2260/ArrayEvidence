"""Safety checks for Nextflow validated staging."""

from pathlib import Path


WORKFLOW = (
    Path(__file__).resolve().parents[1]
    / "workflows"
    / "from_vcf.nf"
)


def test_staging_does_not_enable_overwrite():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "process STAGE_EVIDENCE {" in source
    assert "overwrite: true" not in source


def test_publication_remains_separate():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "process PUBLISH_EVIDENCE {" not in source
    assert "publish_evidence.py" not in source


def test_staging_preflight_guard_exists():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "STAGING_PREFLIGHT" in source


def test_staging_preflight_runs_before_build():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "STAGING_PREFLIGHT" in source
    assert source.index("STAGING_PREFLIGHT") < source.index(
        "bundles = BUILD_EVIDENCE("
    )



def test_workflow_invokes_reservation_helper():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "reserve_staging.py" in source


def test_reservation_uses_nextflow_session_identity():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "workflow.sessionId" in source
    assert "workflow.resume" in source


def test_reservation_precedes_build_submission():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "reserve_staging.py" in source
    assert source.index("reserve_staging.py") < source.index(
        "bundles = BUILD_EVIDENCE("
    )



def test_staging_uses_atomic_finalization():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "finalize_staging.py" in source


def test_staging_does_not_publish_directory_directly():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "publishDir params.staging_dir" not in source


def test_nextflow_declares_configurable_python_interpreter():
    from pathlib import Path
    workflow = (
        Path(__file__).resolve().parents[1]
        / "workflows" / "from_vcf.nf"
    ).read_text(encoding="utf-8")
    assert "params.python_bin = 'python3'" in workflow


def test_nextflow_uses_configured_python_for_both_helpers():
    from pathlib import Path
    workflow = (
        Path(__file__).resolve().parents[1]
        / "workflows" / "from_vcf.nf"
    ).read_text(encoding="utf-8")
    assert "params.python_bin.toString()" in workflow
    assert '${params.python_bin} "${helper}"' in workflow
