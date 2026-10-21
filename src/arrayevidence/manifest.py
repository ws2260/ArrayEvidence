"""Section-aware parser for Illumina CSV assay manifests."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


REQUIRED_COLUMNS = {
    "IlmnID",
    "Name",
    "GenomeBuild",
    "Chr",
    "MapInfo",
    "SNP",
    "SourceSeq",
}


@dataclass(frozen=True)
class Assay:
    ilmn_id: str
    name: str
    genome_build: str
    chrom: str
    map_info: str
    snp: str
    source_seq: str
    ilmn_strand: str = ""
    source_strand: str = ""
    ref_strand: str = ""

    @property
    def coordinate(self) -> tuple[str, int] | None:
        """Return a mapped coordinate, or None for an unmapped assay."""
        if self.chrom in {"", "0", "NA"}:
            return None

        try:
            position = int(self.map_info)
        except ValueError:
            return None

        if position <= 0:
            return None

        return self.chrom, position


def read_assays(path: str | Path) -> Iterator[Assay]:
    """Yield assay definitions from the [Assay] section only."""
    path = Path(path)

    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        in_assay_section = False
        columns = None
        found_assay_section = False

        for row in reader:
            if not row:
                continue

            first = row[0].strip()

            if first.startswith("[") and first.endswith("]"):
                if in_assay_section:
                    break

                in_assay_section = first == "[Assay]"
                if in_assay_section:
                    found_assay_section = True
                continue

            if not in_assay_section:
                continue

            if columns is None:
                columns = {name: i for i, name in enumerate(row)}
                missing = REQUIRED_COLUMNS - columns.keys()

                if missing:
                    raise ValueError(
                        f"Missing manifest columns: {sorted(missing)}"
                    )
                continue

            def field(name: str) -> str:
                index = columns[name]
                if index >= len(row):
                    raise ValueError(
                        f"Malformed assay row: missing {name}"
                    )
                return row[index].strip()

            yield Assay(
                ilmn_id=field("IlmnID"),
                name=field("Name"),
                genome_build=field("GenomeBuild"),
                chrom=field("Chr"),
                map_info=field("MapInfo"),
                snp=field("SNP"),
                source_seq=field("SourceSeq"),
                ilmn_strand=field("IlmnStrand") if "IlmnStrand" in columns else "",
                source_strand=field("SourceStrand") if "SourceStrand" in columns else "",
                ref_strand=field("RefStrand") if "RefStrand" in columns else "",
            )

        if not found_assay_section:
            raise ValueError("Manifest has no [Assay] section")
