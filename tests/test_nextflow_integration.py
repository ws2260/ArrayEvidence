"""Opt-in synthetic Nextflow integration tests."""

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest

from test_cli_build_end_to_end import synthetic_inputs


def run_nextflow(repo, manifest, vcf, reference, executable, stage, resume=None):
    command = [
        "conda", "run", "-n", "arrayevidence-nextflow",
        "nextflow", "run", str(repo / "workflows/from_vcf.nf"),
        "--manifest", str(manifest),
        "--vcf", str(vcf),
        "--reference", str(reference),
        "--sample_id", "PSEUDONYM_001",
        "--arrayevidence_bin", str(executable),
        "--staging_dir", str(stage),
    ]
    if resume:
        command.extend(["-resume", resume])

    return subprocess.run(
        command,
        cwd=repo,
        capture_output=True,
        text=True,
        timeout=180,
    )


def file_hashes(root):
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*")) if p.is_file()
    }


@pytest.mark.skipif(
    os.environ.get("ARRAYEVIDENCE_RUN_NEXTFLOW") != "1",
    reason="Set ARRAYEVIDENCE_RUN_NEXTFLOW=1 to run Nextflow integration",
)
def test_nextflow_atomic_staging_and_completed_resume(tmp_path):
    assert shutil.which("conda"), "Conda executable unavailable"

    repo = Path(__file__).resolve().parents[1]
    executable = Path(sys.executable).parent / "arrayevidence"
    assert executable.is_file(), "ArrayEvidence CLI unavailable"

    manifest, vcf, reference = synthetic_inputs(tmp_path)
    stage = tmp_path / "staging"

    first = run_nextflow(
        repo, manifest, vcf, reference, executable, stage
    )
    assert first.returncode == 0, (first.stdout + first.stderr)[-5000:]

    marker = stage / ".arrayevidence-session"
    validated = stage / "validated"
    bundle = validated / "bundle"
    qc = validated / "qc.json"

    assert marker.is_file()
    assert bundle.is_dir()
    assert qc.is_file()
    assert not list(stage.glob(".arrayevidence-stage-*"))

    session = marker.read_text().strip()
    assert session

    verified = subprocess.run(
        [str(executable), "verify", str(bundle)],
        capture_output=True, text=True,
    )
    assert verified.returncode == 0, verified.stderr

    independent_qc = subprocess.run(
        [str(executable), "qc", str(bundle)],
        capture_output=True, text=True, check=True,
    )
    assert json.loads(qc.read_text()) == json.loads(independent_qc.stdout)

    before = file_hashes(stage)

    resumed = run_nextflow(
        repo, manifest, vcf, reference, executable, stage, resume=session
    )
    log = resumed.stdout + resumed.stderr

    assert resumed.returncode != 0
    assert "STAGING_PREFLIGHT" in log
    assert "BUILD_EVIDENCE (" not in log
    assert file_hashes(stage) == before



@pytest.mark.skipif(
    os.environ.get("ARRAYEVIDENCE_RUN_NEXTFLOW") != "1",
    reason="Set ARRAYEVIDENCE_RUN_NEXTFLOW=1 to run Nextflow integration",
)
def test_nextflow_interrupted_stage_can_resume(tmp_path):
    assert shutil.which("conda"), "Conda executable unavailable"

    repo = Path(__file__).resolve().parents[1]
    executable = Path(sys.executable).parent / "arrayevidence"
    assert executable.is_file()

    manifest, vcf, reference = synthetic_inputs(tmp_path)
    stage = tmp_path / "staging"

    source = (repo / "workflows/from_vcf.nf").read_text()
    needle = '''    test -s qc.json

    ${params.python_bin} "${helper}"'''
    assert source.count(needle) == 1

    # Use a unique temporary workflow within workflows/ so that
    # projectDir-relative helper resolution remains unchanged.
    import uuid
    temporary = (
        repo / "workflows"
        / ("_resume_integration_" + uuid.uuid4().hex + ".nf")
    )

    try:
        interrupted_source = source.replace(
            needle,
            '''    test -s qc.json
    echo SYNTHETIC_CONTROLLED_FAILURE >&2
    false

    ${params.python_bin} "${helper}"''',
            1,
        )
        temporary.write_text(interrupted_source, encoding="utf-8")

        # The existing runner accepts a workflow path through its
        # repository argument only indirectly, so use a local runner.
        def run(resume=None):
            command = [
                "conda", "run", "-n", "arrayevidence-nextflow",
                "nextflow", "run", str(temporary),
                "--manifest", str(manifest),
                "--vcf", str(vcf),
                "--reference", str(reference),
                "--sample_id", "PSEUDONYM_001",
                "--arrayevidence_bin", str(executable),
                "--staging_dir", str(stage),
            ]
            if resume:
                command.extend(["-resume", resume])
            return subprocess.run(
                command,
                cwd=repo,
                capture_output=True,
                text=True,
                timeout=180,
            )

        failed = run()
        failed_log = failed.stdout + failed.stderr

        assert failed.returncode != 0
        assert "SYNTHETIC_CONTROLLED_FAILURE" in failed_log
        assert (stage / ".arrayevidence-session").is_file()
        assert not (stage / "validated").exists()

        session = (stage / ".arrayevidence-session").read_text().strip()
        assert session

        temporary.write_text(source, encoding="utf-8")

        resumed = run(resume=session)
        assert resumed.returncode == 0, (
            resumed.stdout + resumed.stderr
        )[-5000:]

        bundle = stage / "validated" / "bundle"
        qc = stage / "validated" / "qc.json"

        assert bundle.is_dir()
        assert qc.is_file()
        assert not list(stage.glob(".arrayevidence-stage-*"))

        verification = subprocess.run(
            [str(executable), "verify", str(bundle)],
            capture_output=True, text=True,
        )
        assert verification.returncode == 0, verification.stderr

        independent_qc = subprocess.run(
            [str(executable), "qc", str(bundle)],
            capture_output=True, text=True, check=True,
        )
        assert json.loads(qc.read_text()) == json.loads(
            independent_qc.stdout
        )

    finally:
        temporary.unlink(missing_ok=True)
