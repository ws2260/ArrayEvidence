import argparse
import random
from pathlib import Path

import pysam


def generate(output, n):
    output.mkdir(parents=True, exist_ok=False)
    rng = random.Random(20261023)

    # Nonrepetitive synthetic reference; assays spaced 60 bp apart.
    length = 60 * (n + 2)
    sequence = "".join(rng.choices("ACGT", k=length))

    reference = output / "synthetic.fa"
    reference.write_text(
        ">1\n" + sequence + "\n",
        encoding="ascii",
    )
    pysam.faidx(str(reference))

    manifest = output / "synthetic_manifest.csv"
    vcf = output / "synthetic.vcf"

    with manifest.open("w", encoding="utf-8") as m, \
         vcf.open("w", encoding="utf-8") as v:

        m.write(
            "[Heading]\nDescriptor,Value\n"
            "[Assay]\n"
            "IlmnID,Name,GenomeBuild,Chr,MapInfo,SNP,SourceSeq\n"
        )

        v.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            "\tFORMAT\tSYNTHETIC_SAMPLE\n"
        )

        for i in range(n):
            pos = 61 + 60 * i
            ref = sequence[pos - 1]
            alt = next(base for base in "ACGT" if base != ref)

            left = sequence[pos - 21:pos - 1]
            right = sequence[pos:pos + 20]

            m.write(
                f"synthetic_{i},rs_synthetic_{i},38,1,{pos},"
                f"[{ref}/{alt}],{left}[{ref}/{alt}]{right}\n"
            )

            if i % 5 != 0:
                gt = ("0/0", "0/1", "1/1")[i % 3]
                v.write(
                    f"1\t{pos}\tsynthetic_{i}\t{ref}\t{alt}\t."
                    f"\tPASS\t.\tGT\t{gt}\n"
                )

        m.write("[Controls]\nControl,Value\n")

    print(f"Generated {n} synthetic assays")
    print(f"VCF records: {n - (n + 4) // 5}")
    print(f"Output: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--assays", type=int, default=1000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.assays < 1:
        parser.error("--assays must be positive")

    generate(args.output, args.assays)
