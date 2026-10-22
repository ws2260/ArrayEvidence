"""Synthetic tests for atomic staging finalization."""

import json
import subprocess
import sys
from pathlib import Path

HELPER = Path(__file__).resolve().parents[1] / "scripts/finalize_staging.py"


def finalize(bundle, qc, stage, executable, session="session-A"):
    return subprocess.run(
        [
            sys.executable, str(HELPER),
            str(bundle), str(qc), str(stage), str(executable),
            session,
        ],
        capture_output=True, text=True,
    )


def test_finalizer_rejects_missing_bundle(tmp_path):
    stage = tmp_path / "stage"
    stage.mkdir()
    qc = tmp_path / "qc.json"
    qc.write_text(json.dumps({"synthetic": True}))

    result = finalize(
        tmp_path / "missing_bundle", qc, stage, sys.executable
    )
    assert result.returncode != 0
    assert not (stage / "validated").exists()


def test_finalizer_rejects_existing_destination(tmp_path):
    stage = tmp_path / "stage"
    stage.mkdir()
    validated = stage / "validated"
    validated.mkdir()
    sentinel = validated / "DO_NOT_OVERWRITE"
    sentinel.write_text("PRESERVE")

    result = finalize(
        tmp_path / "missing_bundle",
        tmp_path / "missing_qc.json",
        stage,
        sys.executable,
    )

    assert result.returncode != 0
    assert sentinel.read_text() == "PRESERVE"


def test_finalizer_requires_reserved_staging_directory(tmp_path):
    result = finalize(
        tmp_path / "missing_bundle",
        tmp_path / "missing_qc.json",
        tmp_path / "missing_stage",
        sys.executable,
    )
    assert result.returncode != 0
    assert not (tmp_path / "missing_stage").exists()


def test_finalizer_script_exists():
    assert HELPER.is_file()



def test_finalizer_requires_session_argument(tmp_path):
    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / ".arrayevidence-session").write_text("session-A\\n")

    result = subprocess.run(
        [
            sys.executable, str(HELPER),
            str(tmp_path / "missing_bundle"),
            str(tmp_path / "missing_qc.json"),
            str(stage),
            sys.executable,
        ],
        capture_output=True, text=True,
    )

    assert result.returncode != 0
    assert "session" in result.stderr.lower()
    assert not (stage / "validated").exists()



def test_finalizer_rejects_wrong_session(tmp_path):
    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / ".arrayevidence-session").write_text("session-A\\n")

    result = subprocess.run(
        [
            sys.executable, str(HELPER),
            str(tmp_path / "missing_bundle"),
            str(tmp_path / "missing_qc.json"),
            str(stage),
            sys.executable,
            "session-B",
        ],
        capture_output=True, text=True,
    )

    assert result.returncode != 0
    assert "another session" in result.stderr.lower()
    assert not (stage / "validated").exists()



def test_real_synthetic_bundle_finalization(tmp_path):
    import hashlib

    root = Path(__file__).resolve().parents[1]
    from test_cli_build_end_to_end import synthetic_inputs

    manifest, vcf, reference = synthetic_inputs(tmp_path)
    executable = Path(sys.executable).parent / "arrayevidence"
    assert executable.is_file(), "ArrayEvidence CLI is not installed"

    bundle = tmp_path / "bundle"
    qc = tmp_path / "qc.json"
    stage = tmp_path / "stage"

    build = subprocess.run(
        [
            str(executable), "build",
            "--manifest", str(manifest),
            "--vcf", str(vcf),
            "--reference", str(reference),
            "--sample-id", "PSEUDONYM_001",
            "--output", str(bundle),
        ],
        capture_output=True, text=True,
    )
    assert build.returncode == 0, build.stderr

    report = subprocess.run(
        [str(executable), "qc", str(bundle)],
        capture_output=True, text=True, check=True,
    )
    qc.write_text(report.stdout, encoding="utf-8")

    reserve = subprocess.run(
        [
            sys.executable, str(root / "scripts/reserve_staging.py"),
            str(stage), "session-A", "false",
        ],
        capture_output=True, text=True,
    )
    assert reserve.returncode == 0, reserve.stderr

    result = finalize(bundle, qc, stage, executable)
    assert result.returncode == 0, result.stderr

    validated = stage / "validated"
    assert (validated / "bundle").is_dir()
    assert json.loads((validated / "qc.json").read_text()) == json.loads(
        qc.read_text()
    )

    verified = subprocess.run(
        [str(executable), "verify", str(validated / "bundle")],
        capture_output=True, text=True,
    )
    assert verified.returncode == 0, verified.stderr

    def hashes():
        return {
            str(f.relative_to(stage)): hashlib.sha256(
                f.read_bytes()
            ).hexdigest()
            for f in stage.rglob("*") if f.is_file()
        }

    before = hashes()
    duplicate = finalize(bundle, qc, stage, executable)

    assert duplicate.returncode != 0
    assert hashes() == before
    assert not list(stage.glob(".arrayevidence-stage-*"))
