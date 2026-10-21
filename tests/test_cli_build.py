
"""Synthetic CLI contract tests for evidence-bundle construction."""

from pathlib import Path
from unittest.mock import patch

import pytest

from arrayevidence.cli import main


class FakeReference:
    def __init__(self, path):
        self.path = path
        self.closed = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.closed = True


def build_args(tmp_path, sample_id="PSEUDONYM_001"):
    return [
        "build",
        "--manifest", str(tmp_path / "synthetic.csv"),
        "--vcf", str(tmp_path / "synthetic.vcf"),
        "--reference", str(tmp_path / "synthetic.fa"),
        "--sample-id", sample_id,
        "--output", str(tmp_path / "bundle"),
    ]


def test_build_orchestrates_existing_pipeline(tmp_path, capsys):
    args = build_args(tmp_path)
    reference = FakeReference(tmp_path / "synthetic.fa")

    assays = [object()]
    records = [object()]
    results = [object()]
    assay_rows = [object()]
    candidate_rows = [object()]

    with (
        patch(
            "arrayevidence.cli.read_assays",
            return_value=iter(assays),
            create=True,
        ) as read_assays,
        patch(
            "arrayevidence.cli.read_vcf",
            return_value=iter(records),
            create=True,
        ) as read_vcf,
        patch(
            "arrayevidence.cli.PysamReferenceReader",
            return_value=reference,
            create=True,
        ) as reference_reader,
        patch(
            "arrayevidence.cli.reconcile",
            return_value=results,
            create=True,
        ) as reconcile,
        patch(
            "arrayevidence.cli.build_evidence_rows",
            return_value=(assay_rows, candidate_rows),
            create=True,
        ) as build_rows,
        patch(
            "arrayevidence.cli.export_evidence_tsv",
            return_value=(
                tmp_path / "bundle" / "assay_evidence.tsv",
                tmp_path / "bundle" / "allele_candidates.tsv",
            ),
            create=True,
        ) as export,
    ):
        code = main(args)

    assert code == 0
    assert reference.closed

    read_assays.assert_called_once_with(
        str(tmp_path / "synthetic.csv")
    )
    read_vcf.assert_called_once_with(
        str(tmp_path / "synthetic.vcf")
    )
    reference_reader.assert_called_once_with(
        str(tmp_path / "synthetic.fa")
    )

    reconcile.assert_called_once_with(
        assays,
        records,
        reference=reference,
    )

    build_rows.assert_called_once_with(
        "PSEUDONYM_001",
        results,
    )

    export.assert_called_once_with(
        str(tmp_path / "bundle"),
        assay_rows,
        candidate_rows,
    )

    output = capsys.readouterr()
    assert "PSEUDONYM_001" not in output.out
    assert "PSEUDONYM_001" not in output.err


def test_build_accepts_operator_supplied_pseudonym(tmp_path):
    args = build_args(
        tmp_path,
        sample_id="ANONYMOUS_CODE_42",
    )

    with (
        patch(
            "arrayevidence.cli.read_assays",
            return_value=iter([]),
            create=True,
        ),
        patch(
            "arrayevidence.cli.read_vcf",
            return_value=iter([]),
            create=True,
        ),
        patch(
            "arrayevidence.cli.PysamReferenceReader",
            return_value=FakeReference("synthetic.fa"),
            create=True,
        ),
        patch(
            "arrayevidence.cli.reconcile",
            return_value=[],
            create=True,
        ),
        patch(
            "arrayevidence.cli.build_evidence_rows",
            return_value=([], []),
            create=True,
        ) as build_rows,
        patch(
            "arrayevidence.cli.export_evidence_tsv",
            return_value=(Path("a"), Path("b")),
            create=True,
        ),
    ):
        assert main(args) == 0

    build_rows.assert_called_once_with(
        "ANONYMOUS_CODE_42",
        [],
    )


def test_build_refuses_existing_output_before_reading(
    tmp_path,
    capsys,
):
    (tmp_path / "bundle").mkdir()

    with patch(
        "arrayevidence.cli.read_assays",
        create=True,
    ) as read_assays:
        code = main(build_args(tmp_path))

    assert code == 1
    read_assays.assert_not_called()

    output = capsys.readouterr()
    assert "PSEUDONYM_001" not in output.err


def test_build_sanitizes_input_errors(tmp_path, capsys):
    private_marker = "PRIVATE_SAMPLE_IDENTIFIER_123"

    with patch(
        "arrayevidence.cli.read_assays",
        side_effect=ValueError(private_marker),
        create=True,
    ):
        code = main(build_args(tmp_path))

    assert code == 1

    output = capsys.readouterr()

    assert private_marker not in output.out
    assert private_marker not in output.err


@pytest.mark.parametrize(
    "omitted",
    [
        "--manifest",
        "--vcf",
        "--reference",
        "--sample-id",
        "--output",
    ],
)
def test_build_requires_all_arguments(tmp_path, omitted):
    args = build_args(tmp_path)
    index = args.index(omitted)
    del args[index:index + 2]

    with pytest.raises(SystemExit) as error:
        main(args)

    assert error.value.code == 2


def test_build_rejects_blank_sample_id(tmp_path, capsys):
    code = main(
        build_args(tmp_path, sample_id="   ")
    )

    assert code == 2

    output = capsys.readouterr()
    assert "PSEUDONYM_001" not in output.err
