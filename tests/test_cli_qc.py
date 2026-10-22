"""Synthetic CLI tests for read-only aggregate QC."""

import json
from pathlib import Path

from arrayevidence.cli import main
from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.qc import summarize_bundle


def test_qc_valid_bundle_outputs_deterministic_json(
    tmp_path, capsys
):
    bundle = tmp_path / "private_synthetic_identifier"
    export_evidence_tsv(bundle, [], [])

    before = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    first_code = main(["qc", str(bundle)])
    first = capsys.readouterr()

    second_code = main(["qc", str(bundle)])
    second = capsys.readouterr()

    assert first_code == second_code == 0
    assert first.out == second.out
    assert first.err == second.err == ""

    report = json.loads(first.out)

    assert report == summarize_bundle(bundle)
    assert report["assay_rows"] == 0
    assert report["candidate_rows"] == 0
    assert first.out.endswith("\n")

    assert "private_synthetic_identifier" not in first.out
    assert str(tmp_path) not in first.out

    after = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    assert before == after


def test_qc_missing_bundle_reports_sanitized_error(
    tmp_path, capsys
):
    private_path = tmp_path / "PRIVATE_PATH_98123"

    code = main(["qc", str(private_path)])
    output = capsys.readouterr()

    assert code == 1
    assert output.out == ""
    assert output.err == (
        "FAIL [QC_ERROR]: Evidence QC failed.\n"
    )
    assert "PRIVATE_PATH_98123" not in output.err
    assert str(tmp_path) not in output.err


def test_qc_corrupted_bundle_reports_sanitized_error(
    tmp_path, capsys
):
    bundle = tmp_path / "private_bundle"
    export_evidence_tsv(bundle, [], [])

    with (bundle / "assay_evidence.tsv").open(
        "a", encoding="utf-8"
    ) as handle:
        handle.write("\n")

    code = main(["qc", str(bundle)])
    output = capsys.readouterr()

    assert code == 1
    assert output.out == ""
    assert output.err == (
        "FAIL [QC_ERROR]: Evidence QC failed.\n"
    )


def test_qc_unexpected_exception_is_sanitized(
    monkeypatch, tmp_path, capsys
):
    import arrayevidence.qc as qc

    secret = "PRIVATE_INTERNAL_VALUE_72981"

    def broken_summary(_bundle):
        raise RuntimeError(secret)

    monkeypatch.setattr(qc, "summarize_bundle", broken_summary)

    code = main(["qc", str(tmp_path)])
    output = capsys.readouterr()

    assert code == 1
    assert output.out == ""
    assert output.err == (
        "FAIL [QC_ERROR]: Evidence QC failed.\n"
    )
    assert secret not in output.err


def test_qc_requires_bundle_argument():
    import pytest

    with pytest.raises(SystemExit) as exc:
        main(["qc"])

    assert exc.value.code == 2


def test_qc_help_is_available():
    import pytest

    with pytest.raises(SystemExit) as exc:
        main(["qc", "--help"])

    assert exc.value.code == 0
