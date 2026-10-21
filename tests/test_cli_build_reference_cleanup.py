"""Synthetic CLI reference-cleanup error classification tests."""

from unittest.mock import patch

from arrayevidence.cli import main


class FailingCleanupReference:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        raise RuntimeError("SYNTHETIC_PRIVATE_CLEANUP_MARKER_987")


def build_args(tmp_path):
    return [
        "build",
        "--manifest", str(tmp_path / "synthetic.csv"),
        "--vcf", str(tmp_path / "synthetic.vcf"),
        "--reference", str(tmp_path / "synthetic.fa"),
        "--sample-id", "PSEUDONYM_001",
        "--output", str(tmp_path / "bundle"),
    ]


def test_reference_cleanup_failure_reports_reference_error(
    tmp_path,
    capsys,
):
    with (
        patch(
            "arrayevidence.cli.read_assays",
            return_value=iter([object()]),
        ),
        patch(
            "arrayevidence.cli.read_vcf",
            return_value=iter([object()]),
        ),
        patch(
            "arrayevidence.cli.PysamReferenceReader",
            return_value=FailingCleanupReference(),
        ),
        patch(
            "arrayevidence.cli.reconcile",
            return_value=[object()],
        ) as reconcile,
        patch(
            "arrayevidence.cli.build_evidence_rows",
        ) as build_rows,
        patch(
            "arrayevidence.cli.export_evidence_tsv",
        ) as export,
    ):
        code = main(build_args(tmp_path))

    output = capsys.readouterr()
    combined = output.out + output.err

    assert code == 1
    reconcile.assert_called_once()
    build_rows.assert_not_called()
    export.assert_not_called()

    assert "REFERENCE_ERROR" in output.err
    assert "RECONCILIATION_ERROR" not in output.err
    assert "SYNTHETIC_PRIVATE_CLEANUP_MARKER_987" not in combined
    assert "PSEUDONYM_001" not in combined
    assert "Traceback" not in combined
    assert not (tmp_path / "bundle").exists()
