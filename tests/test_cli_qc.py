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


def test_qc_output_file_matches_stdout(tmp_path, capsys):
    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    before = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    assert main(["qc", str(bundle)]) == 0
    expected = capsys.readouterr().out

    output = tmp_path / "qc.json"
    assert main([
        "qc", str(bundle), "--output", str(output)
    ]) == 0

    captured = capsys.readouterr()

    assert captured.out == ""
    assert captured.err == ""
    assert output.read_bytes() == expected.encode("utf-8")

    after = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }
    assert before == after


def test_qc_output_refuses_existing_file(tmp_path, capsys):
    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    output = tmp_path / "existing.json"
    output.write_text("DO_NOT_OVERWRITE")

    code = main([
        "qc", str(bundle), "--output", str(output)
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.\n"
    )
    assert output.read_text() == "DO_NOT_OVERWRITE"


def test_qc_output_invalid_bundle_creates_no_file(
    tmp_path, capsys
):
    output = tmp_path / "qc.json"

    code = main([
        "qc", str(tmp_path / "missing_bundle"),
        "--output", str(output),
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_ERROR]: Evidence QC failed.\n"
    )
    assert not output.exists()


def test_qc_output_rejects_directory(tmp_path, capsys):
    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    output = tmp_path / "existing_directory"
    output.mkdir()

    code = main([
        "qc", str(bundle), "--output", str(output)
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.\n"
    )


def test_qc_output_help_documents_option(capsys):
    import pytest

    with pytest.raises(SystemExit) as exc:
        main(["qc", "--help"])

    captured = capsys.readouterr()

    assert exc.value.code == 0
    assert "--output" in captured.out


def test_qc_output_link_failure_cleans_temporary_file(
    tmp_path, capsys, monkeypatch
):
    import os

    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    destination = tmp_path / "qc.json"

    def fail_link(_source, _destination):
        raise OSError("PRIVATE_INTERNAL_ERROR")

    monkeypatch.setattr(os, "link", fail_link)

    code = main([
        "qc", str(bundle),
        "--output", str(destination),
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.\n"
    )
    assert "PRIVATE_INTERNAL_ERROR" not in captured.err
    assert not destination.exists()
    assert list(tmp_path.glob(".arrayevidence-qc-*.tmp")) == []


def test_qc_output_existing_file_leaves_no_temporary_file(
    tmp_path, capsys
):
    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    destination = tmp_path / "qc.json"
    destination.write_bytes(b"EXISTING_CONTENT")

    code = main([
        "qc", str(bundle),
        "--output", str(destination),
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.\n"
    )
    assert destination.read_bytes() == b"EXISTING_CONTENT"
    assert list(tmp_path.glob(".arrayevidence-qc-*.tmp")) == []


def test_qc_output_cleanup_failure_is_sanitized(
    tmp_path, capsys, monkeypatch
):
    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    destination = tmp_path / "qc.json"

    original_unlink = Path.unlink

    def fail_temporary_cleanup(self, *args, **kwargs):
        if self.name.startswith(".arrayevidence-qc-"):
            raise PermissionError("PRIVATE_CLEANUP_ERROR")
        return original_unlink(self, *args, **kwargs)

    monkeypatch.setattr(Path, "unlink", fail_temporary_cleanup)

    code = main([
        "qc", str(bundle),
        "--output", str(destination),
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_OUTPUT_ERROR]: Evidence QC output cleanup failed.\n"
    )
    assert "PRIVATE_CLEANUP_ERROR" not in captured.err
    assert destination.is_file()
    assert json.loads(destination.read_text()) == summarize_bundle(bundle)


def test_qc_output_rejects_destination_inside_bundle(
    tmp_path, capsys
):
    from arrayevidence.bundle_verify import verify_bundle

    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    destination = bundle / "qc.json"

    before = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    assert verify_bundle(bundle).valid

    code = main([
        "qc", str(bundle),
        "--output", str(destination),
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.\n"
    )

    assert not destination.exists()
    assert verify_bundle(bundle).valid

    after = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }
    assert before == after


def test_qc_output_rejects_symlink_into_bundle(
    tmp_path, capsys
):
    from arrayevidence.bundle_verify import verify_bundle

    bundle = tmp_path / "synthetic_bundle"
    export_evidence_tsv(bundle, [], [])

    alias = tmp_path / "bundle_alias"
    alias.symlink_to(bundle, target_is_directory=True)

    destination = alias / "qc.json"

    before = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }

    code = main([
        "qc", str(bundle),
        "--output", str(destination),
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert captured.out == ""
    assert captured.err == (
        "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.\n"
    )
    assert not destination.exists()
    assert verify_bundle(bundle).valid

    after = {
        path.name: path.read_bytes()
        for path in bundle.iterdir()
    }
    assert before == after
