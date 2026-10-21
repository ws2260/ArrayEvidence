"""Real CLI build/verify tests using synthetic files only."""

import csv
import subprocess
import sys
from pathlib import Path

import pysam


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


def read_table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def synthetic_inputs(tmp_path):
    reference = tmp_path / "synthetic.fa"
    manifest = tmp_path / "synthetic_manifest.csv"
    vcf = tmp_path / "synthetic.vcf"

    sequence = (
        "ACGT" * 15
        + "G"
        + "TGCA" * 15
    )

    position = 61
    left = sequence[position - 21:position - 1]
    right = sequence[position:position + 20]

    assert len(left) == 20
    assert len(right) == 20
    assert sequence[position - 1] == "G"

    reference.write_text(
        ">1\n" + sequence + "\n",
        encoding="ascii",
    )
    pysam.faidx(str(reference))

    manifest.write_text(
        "[Heading]\n"
        "Descriptor,Value\n"
        "[Assay]\n"
        "IlmnID,Name,GenomeBuild,Chr,MapInfo,SNP,SourceSeq\n"
        f"synthetic_multi,rs_synthetic,38,1,{position},"
        f"[G/A],{left}[G/A]{right}\n"
        f"synthetic_absent,rs_absent,38,1,80,"
        f"[A/C],{left}[A/C]{right}\n"
        "[Controls]\n"
        "Control,Value\n",
        encoding="utf-8",
    )

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        "\tFORMAT\tORIGINAL_SYNTHETIC_VCF_SAMPLE\n"
        f"1\t{position}\tsynthetic_variant\tG\tC,A\t."
        "\tPASS\t.\tGT\t0|2\n",
        encoding="utf-8",
    )

    return manifest, vcf, reference


def build_args(manifest, vcf, reference, output):
    return (
        "build",
        "--manifest", manifest,
        "--vcf", vcf,
        "--reference", reference,
        "--sample-id", "PSEUDONYM_001",
        "--output", output,
    )


def test_build_verify_determinism_and_no_overwrite(tmp_path):
    manifest, vcf, reference = synthetic_inputs(tmp_path)

    first = tmp_path / "bundle_first"
    second = tmp_path / "bundle_second"

    args_first = build_args(manifest, vcf, reference, first)
    args_second = build_args(manifest, vcf, reference, second)

    result = run_cli(*args_first)

    assert result.returncode == 0, (
        result.stdout,
        result.stderr,
    )

    assert "PSEUDONYM_001" not in result.stdout
    assert "ORIGINAL_SYNTHETIC_VCF_SAMPLE" not in (
        result.stdout + result.stderr
    )

    expected_files = {
        "assay_evidence.tsv",
        "allele_candidates.tsv",
        "bundle_manifest.json",
    }

    assert {p.name for p in first.iterdir()} == expected_files

    assays = read_table(first / "assay_evidence.tsv")
    candidates = read_table(first / "allele_candidates.tsv")

    assert len(assays) == 2

    by_id = {row["assay_id"]: row for row in assays}

    multi = by_id["synthetic_multi"]
    absent = by_id["synthetic_absent"]

    assert multi["allele_evidence"] == "PARTIAL_MATCH"
    assert multi["projection_status"] == "ASSAY_ALT_PRESENT"
    assert multi["assay_alt_dosage"] == "1"

    assert absent["coordinate_status"] == "VCF_ABSENT"
    assert absent["coordinate_genotype_evidence"] == "NOT_EVALUATED"
    assert absent["assay_alt_dosage"] == "NULL"

    assert len(candidates) == 1
    assert candidates[0]["assay_id"] == "synthetic_multi"
    assert candidates[0]["vcf_alt_index"] == "2"
    assert candidates[0]["vcf_alt"] == "A"
    assert candidates[0]["vcf_gt"] == "0|2"

    for row in assays + candidates:
        assert row["sample_id"] == "PSEUDONYM_001"

    all_bytes = b"".join(
        (first / name).read_bytes()
        for name in sorted(expected_files)
    )

    assert b"ORIGINAL_SYNTHETIC_VCF_SAMPLE" not in all_bytes

    verified = run_cli("verify", first)
    assert verified.returncode == 0, (
        verified.stdout,
        verified.stderr,
    )

    second_result = run_cli(*args_second)
    assert second_result.returncode == 0, (
        second_result.stdout,
        second_result.stderr,
    )

    for name in expected_files:
        assert (first / name).read_bytes() == (
            second / name
        ).read_bytes()

    before = {
        name: (first / name).read_bytes()
        for name in expected_files
    }

    repeated = run_cli(*args_first)

    assert repeated.returncode == 1
    assert "PSEUDONYM_001" not in (
        repeated.stdout + repeated.stderr
    )

    after = {
        name: (first / name).read_bytes()
        for name in expected_files
    }

    assert before == after
