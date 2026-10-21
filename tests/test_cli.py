"""Synthetic contract tests for the ArrayEvidence CLI."""

import subprocess
import sys
from importlib.metadata import entry_points
from pathlib import Path

from arrayevidence.evidence_export import export_evidence_tsv


def run_cli(*arguments):
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "arrayevidence",
            *map(str, arguments),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def test_verify_valid_bundle(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    result = run_cli("verify", bundle)

    assert result.returncode == 0
    assert "PASS: Bundle verification succeeded." in result.stdout
    assert result.stderr == ""


def test_verify_invalid_bundle(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    (bundle / "assay_evidence.tsv").write_bytes(
        b"tampered\n"
    )

    result = run_cli("verify", bundle)

    assert result.returncode == 1
    assert "FAIL: Bundle verification failed." in result.stderr
    assert "SHA-256 mismatch" in result.stderr


def test_verify_missing_bundle(tmp_path):
    missing = tmp_path / "missing"

    result = run_cli("verify", missing)

    assert result.returncode == 1
    assert "FAIL: Bundle verification failed." in result.stderr


def test_verify_does_not_modify_bundle(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    before = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    result = run_cli("verify", bundle)

    after = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    assert result.returncode == 0
    assert before == after


def test_help_exits_successfully():
    result = run_cli("--help")

    assert result.returncode == 0
    assert "verify" in result.stdout
    assert result.stderr == ""


def test_verify_help_exits_successfully():
    result = run_cli("verify", "--help")

    assert result.returncode == 0
    assert "bundle" in result.stdout.lower()
    assert result.stderr == ""


def test_missing_command_returns_usage_error():
    result = run_cli()

    assert result.returncode == 2
    assert "usage:" in result.stderr.lower()


def test_missing_bundle_argument_returns_usage_error():
    result = run_cli("verify")

    assert result.returncode == 2
    assert "usage:" in result.stderr.lower()


def test_unknown_command_returns_usage_error():
    result = run_cli("unknown")

    assert result.returncode == 2
    assert "invalid choice" in result.stderr.lower()


def test_console_script_is_declared():
    scripts = entry_points(group="console_scripts")
    matches = [
        script
        for script in scripts
        if script.name == "arrayevidence"
    ]

    assert len(matches) == 1
    assert matches[0].value == "arrayevidence.cli:main"


def test_cli_does_not_print_evidence_rows(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    result = run_cli("verify", bundle)

    combined = result.stdout + result.stderr

    assert result.returncode == 0
    assert "sample_id" not in combined
    assert "assay_id" not in combined
    assert "vcf_gt" not in combined


def test_verify_rejects_extra_positional_arguments(tmp_path):
    bundle = tmp_path / "bundle"
    export_evidence_tsv(bundle, [], [])

    result = run_cli("verify", bundle, "unexpected")

    assert result.returncode == 2
    assert "usage:" in result.stderr.lower()
