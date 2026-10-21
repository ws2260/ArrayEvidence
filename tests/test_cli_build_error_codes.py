"""Synthetic stage-specific CLI error-reporting tests."""

from unittest.mock import patch

import pytest

from arrayevidence.cli import main


class FakeReference:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return None


def build_args(tmp_path):
    return [
        "build",
        "--manifest", str(tmp_path / "synthetic.csv"),
        "--vcf", str(tmp_path / "synthetic.vcf"),
        "--reference", str(tmp_path / "synthetic.fa"),
        "--sample-id", "PSEUDONYM_001",
        "--output", str(tmp_path / "bundle"),
    ]


STAGES = [
    ("read_assays", "MANIFEST_READ_ERROR"),
    ("read_vcf", "VCF_READ_ERROR"),
    ("PysamReferenceReader", "REFERENCE_ERROR"),
    ("reconcile", "RECONCILIATION_ERROR"),
    ("build_evidence_rows", "EVIDENCE_ERROR"),
    ("export_evidence_tsv", "EXPORT_ERROR"),
]


@pytest.mark.parametrize("failure_stage,error_code", STAGES)
def test_build_reports_sanitized_stage_error(
    tmp_path,
    capsys,
    failure_stage,
    error_code,
):
    secret = "SYNTHETIC_PRIVATE_MARKER_987"

    defaults = {
        "read_assays": iter([object()]),
        "read_vcf": iter([object()]),
        "PysamReferenceReader": FakeReference(),
        "reconcile": [object()],
        "build_evidence_rows": ([object()], [object()]),
        "export_evidence_tsv": (object(), object()),
    }

    patches = {}

    for name, result in defaults.items():
        patches[name] = patch(
            f"arrayevidence.cli.{name}",
            return_value=result,
        )

    active = []

    try:
        for name, context in patches.items():
            mocked = context.__enter__()
            active.append(context)

            if name == failure_stage:
                mocked.side_effect = RuntimeError(secret)

        code = main(build_args(tmp_path))

    finally:
        for context in reversed(active):
            context.__exit__(None, None, None)

    output = capsys.readouterr()
    combined = output.out + output.err

    assert code == 1
    assert error_code in output.err
    assert secret not in combined
    assert "PSEUDONYM_001" not in combined
    assert "Traceback" not in combined
    assert not (tmp_path / "bundle").exists()


def test_missing_manifest_reports_manifest_code(tmp_path, capsys):
    code = main(build_args(tmp_path))
    output = capsys.readouterr()

    assert code == 1
    assert "MANIFEST_READ_ERROR" in output.err
    assert "PSEUDONYM_001" not in output.err
    assert not (tmp_path / "bundle").exists()
