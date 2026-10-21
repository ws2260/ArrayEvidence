"""Synthetic operational-failure tests for build CLI."""

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from arrayevidence.cli import main


def arguments(tmp_path):
    return [
        "build",
        "--manifest", str(tmp_path / "manifest.csv"),
        "--vcf", str(tmp_path / "sample.vcf"),
        "--reference", str(tmp_path / "reference.fa"),
        "--sample-id", "PSEUDONYM_001",
        "--output", str(tmp_path / "bundle"),
    ]


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


def test_unexpected_exception_is_sanitized(tmp_path, capsys):
    secret = "SYNTHETIC_PRIVATE_MARKER_987"

    with patch(
        "arrayevidence.cli.read_assays",
        side_effect=RuntimeError(secret),
    ):
        code = main(arguments(tmp_path))

    output = capsys.readouterr()

    assert code == 1
    assert secret not in output.out + output.err
    assert "Traceback" not in output.out + output.err


def test_missing_manifest_fails_without_bundle(tmp_path):
    args = arguments(tmp_path)

    result = run_cli(*args)

    assert result.returncode == 1
    assert not (tmp_path / "bundle").exists()
    assert "Traceback" not in result.stderr
    assert "PSEUDONYM_001" not in result.stderr


def test_malformed_vcf_fails_without_bundle(tmp_path):
    manifest = tmp_path / "manifest.csv"
    vcf = tmp_path / "sample.vcf"

    manifest.write_text(
        "[Assay]\n"
        "IlmnID,Name,GenomeBuild,Chr,MapInfo,SNP,SourceSeq\n"
        "synthetic,rs1,38,1,30,[A/G],AAAA[A/G]CCCC\n",
        encoding="utf-8",
    )

    vcf.write_text(
        "not a valid VCF\n",
        encoding="utf-8",
    )

    result = run_cli(*arguments(tmp_path))

    assert result.returncode == 1
    assert not (tmp_path / "bundle").exists()
    assert "Traceback" not in result.stderr


def test_missing_fasta_index_fails_without_bundle(tmp_path):
    manifest = tmp_path / "manifest.csv"
    vcf = tmp_path / "sample.vcf"
    reference = tmp_path / "reference.fa"

    manifest.write_text(
        "[Assay]\n"
        "IlmnID,Name,GenomeBuild,Chr,MapInfo,SNP,SourceSeq\n"
        "synthetic,rs1,38,1,30,[A/G],AAAA[A/G]CCCC\n",
        encoding="utf-8",
    )

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        "\tFORMAT\tSYNTHETIC_SAMPLE\n"
        "1\t30\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )

    reference.write_text(
        ">1\n" + "A" * 100 + "\n",
        encoding="ascii",
    )

    assert not Path(str(reference) + ".fai").exists()

    result = run_cli(*arguments(tmp_path))

    assert result.returncode == 1
    assert not (tmp_path / "bundle").exists()
    assert "Traceback" not in result.stderr
