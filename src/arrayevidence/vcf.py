"""Minimal, deterministic VCF parsing for assay reconciliation."""

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass(frozen=True)
class VCFRecord:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alts: tuple[str, ...]
    gt: str

    @property
    def coordinate(self) -> tuple[str, int]:
        return self.chrom, self.pos

    @property
    def call_status(self) -> str:
        if not self.gt or self.gt == ".":
            return "NO_CALL"

        alleles = self.gt.replace("|", "/").split("/")

        if "." in alleles:
            return "NO_CALL"

        if all(allele == "0" for allele in alleles):
            return "HOM_REF"

        return "CALLED_NONREF"


def read_vcf(path: str | Path) -> Iterator[VCFRecord]:
    """Read genotype records from a single-sample VCF."""
    path = Path(path)

    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                raise ValueError("Expected a single-sample VCF record")

            chrom, pos, variant_id, ref, alt = fields[:5]
            format_keys = fields[8].split(":")
            sample_values = fields[9].split(":")

            if "GT" not in format_keys:
                raise ValueError("VCF record has no GT field")

            gt_index = format_keys.index("GT")

            if gt_index >= len(sample_values):
                raise ValueError("VCF record has no genotype value")

            yield VCFRecord(
                chrom=chrom,
                pos=int(pos),
                variant_id=variant_id,
                ref=ref,
                alts=tuple(alt.split(",")) if alt != "." else (),
                gt=sample_values[gt_index],
            )
