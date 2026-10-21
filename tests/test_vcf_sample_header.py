"""Synthetic VCF sample-header validation contract."""

import gzip
from pathlib import Path

import pytest

from arrayevidence.vcf import read_vcf


HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
    "INFO\tFORMAT\tSAMPLE_A\n"
)

RECORD = (
    "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
)


def write_vcf(path: Path, content: str):
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(content)
    else:
        path.write_text(content, encoding="utf-8")
    return path


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_single_sample_vcf_is_accepted(tmp_path, suffix):
    path = write_vcf(
        tmp_path / ("synthetic" + suffix),
        HEADER + RECORD,
    )

    records = list(read_vcf(path))

    assert len(records) == 1
    assert records[0].gt == "0/1"


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_multiple_samples_are_rejected(tmp_path, suffix):
    content = (
        HEADER.replace(
            "SAMPLE_A\n",
            "SAMPLE_A\tSAMPLE_B\n",
        )
        + RECORD.rstrip("\n")
        + "\t0/0\n"
    )

    path = write_vcf(
        tmp_path / ("multisample" + suffix),
        content,
    )

    with pytest.raises(ValueError, match="single.sample"):
        list(read_vcf(path))


def test_site_only_vcf_is_rejected(tmp_path):
    content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
        "FILTER\tINFO\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\n"
    )

    path = write_vcf(tmp_path / "site_only.vcf", content)

    with pytest.raises(ValueError, match="header|single.sample"):
        list(read_vcf(path))


def test_missing_chrom_header_is_rejected(tmp_path):
    path = write_vcf(
        tmp_path / "missing_header.vcf",
        "##fileformat=VCFv4.2\n" + RECORD,
    )

    with pytest.raises(ValueError, match="header"):
        list(read_vcf(path))


def test_malformed_chrom_header_is_rejected(tmp_path):
    content = HEADER.replace(
        "#CHROM\tPOS\tID",
        "#CHROM\tWRONG\tID",
    )

    path = write_vcf(
        tmp_path / "malformed_header.vcf",
        content + RECORD,
    )

    with pytest.raises(ValueError, match="header"):
        list(read_vcf(path))


def test_empty_vcf_without_header_is_rejected(tmp_path):
    path = write_vcf(
        tmp_path / "empty.vcf",
        "",
    )

    with pytest.raises(ValueError, match="header"):
        list(read_vcf(path))


def test_header_only_single_sample_vcf_is_accepted(tmp_path):
    path = write_vcf(
        tmp_path / "header_only.vcf",
        HEADER,
    )

    assert list(read_vcf(path)) == []


def test_record_with_extra_sample_column_is_rejected(tmp_path):
    path = write_vcf(
        tmp_path / "extra_record_column.vcf",
        HEADER + RECORD.rstrip("\n") + "\t0/0\n",
    )

    with pytest.raises(ValueError, match="sample|column"):
        list(read_vcf(path))


def test_record_with_missing_sample_column_is_rejected(tmp_path):
    record = RECORD.rsplit("\t", 1)[0] + "\n"
    path = write_vcf(
        tmp_path / "missing_record_column.vcf",
        HEADER + record,
    )

    with pytest.raises(ValueError, match="sample|column"):
        list(read_vcf(path))


def test_duplicate_chrom_header_is_rejected(tmp_path):
    path = write_vcf(
        tmp_path / "duplicate_header.vcf",
        HEADER + HEADER.splitlines(keepends=True)[1] + RECORD,
    )

    with pytest.raises(ValueError, match="header"):
        list(read_vcf(path))


def test_empty_sample_name_is_rejected(tmp_path):
    header = HEADER.replace(
        "SAMPLE_A" + chr(10),
        chr(10),
    )

    path = write_vcf(
        tmp_path / "empty_sample.vcf",
        header + RECORD,
    )

    with pytest.raises(ValueError, match="header"):
        list(read_vcf(path))


def test_metadata_after_chrom_header_is_rejected(tmp_path):
    path = write_vcf(
        tmp_path / "late_metadata.vcf",
        HEADER + "##source=synthetic" + chr(10) + RECORD,
    )

    with pytest.raises(ValueError, match="header"):
        list(read_vcf(path))


def test_missing_gt_format_field_is_rejected(tmp_path):
    record = RECORD.replace(
        chr(9) + "GT" + chr(9),
        chr(9) + "DP" + chr(9),
    )

    path = write_vcf(
        tmp_path / "missing_gt.vcf",
        HEADER + record,
    )

    with pytest.raises(ValueError, match="GT field"):
        list(read_vcf(path))


def test_missing_genotype_value_is_rejected(tmp_path):
    record = RECORD.rsplit(chr(9), 1)[0] + chr(9) + chr(10)

    path = write_vcf(
        tmp_path / "missing_gt_value.vcf",
        HEADER + record,
    )

    with pytest.raises(ValueError, match="genotype value"):
        list(read_vcf(path))
