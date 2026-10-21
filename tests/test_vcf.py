from pathlib import Path

from arrayevidence.vcf import read_vcf


def test_vcf_genotype_states(tmp_path: Path):
    path = tmp_path / "synthetic.vcf"

    path.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\n"
        "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
        "1\t300\trs3\tG\tA\t.\tPASS\t.\tGT\t./.\n"
        "X\t400\trs4\tA\tC\t.\tPASS\t.\tGT\t.\n"
        "1\t500\trs5\tA\tC,G\t.\tPASS\t.\tGT\t0/2\n",
        encoding="utf-8",
    )

    records = list(read_vcf(path))

    assert len(records) == 5
    assert records[0].call_status == "HOM_REF"
    assert records[1].call_status == "CALLED_NONREF"
    assert records[2].call_status == "NO_CALL"
    assert records[3].call_status == "NO_CALL"
    assert records[4].alts == ("C", "G")
    assert records[4].call_status == "CALLED_NONREF"
