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
    """Read genotype records from a strictly single-sample VCF.

    Require exactly one sample in the #CHROM header. Reject malformed
    headers and record widths rather than silently selecting a sample.
    """

    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open

    required_columns = (
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    )

    header_seen = False

    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                if header_seen:
                    raise ValueError(
                        "Malformed VCF header: metadata after #CHROM"
                    )
                continue

            if line.startswith("#"):
                if header_seen:
                    raise ValueError(
                        "Malformed VCF header: duplicate #CHROM"
                    )

                columns = line.rstrip('\r\n').split('\t')

                if tuple(columns[:9]) != required_columns:
                    raise ValueError(
                        "Malformed VCF header: invalid columns"
                    )

                if len(columns) != 10:
                    raise ValueError(
                        "Expected a single-sample VCF header"
                    )

                if not columns[9].strip():
                    raise ValueError(
                        "Malformed VCF header: empty sample name"
                    )

                header_seen = True
                continue

            if not header_seen:
                raise ValueError(
                    "Missing or malformed VCF header"
                )

            fields = line.rstrip('\r\n').split('\t')

            if len(fields) != 10:
                raise ValueError(
                    "VCF record sample column count mismatch"
                )

            chrom, pos, variant_id, ref, alt = fields[:5]
            format_keys = fields[8].split(":")
            sample_values = fields[9].split(":")

            if "GT" not in format_keys:
                raise ValueError("VCF record has no GT field")

            gt_index = format_keys.index("GT")

            if (
                gt_index >= len(sample_values)
                or not sample_values[gt_index]
            ):
                raise ValueError(
                    "VCF record has no genotype value"
                )

            yield VCFRecord(
                chrom=chrom,
                pos=int(pos),
                variant_id=variant_id,
                ref=ref,
                alts=tuple(alt.split(",")) if alt != "." else (),
                gt=sample_values[gt_index],
            )

    if not header_seen:
        raise ValueError("Missing VCF header")
