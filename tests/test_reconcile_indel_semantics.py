"""Semantic regression tests for indel reconciliation.

Coordinate evidence and normalized allele evidence are separate.
Indel assays must not be processed through the SNV validation path.
"""

from unittest.mock import patch

from arrayevidence.manifest import Assay
from arrayevidence.reconcile import reconcile
from arrayevidence.vcf import VCFRecord


class SyntheticReference:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, chrom, start, end):
        if chrom != "1":
            raise KeyError(chrom)
        if start < 0 or end > len(self.sequence):
            raise ValueError("Out of bounds")
        return self.sequence[start:end]

    def get_reference_length(self, chrom):
        if chrom != "1":
            raise KeyError(chrom)
        return len(self.sequence)


def assay(source_seq="ACGT[-/TT]GCTA", pos=5):
    return Assay(
        ilmn_id="semantic_indel",
        name="semantic_indel",
        genome_build="38",
        chrom="1",
        map_info=str(pos),
        snp="[I/D]",
        source_seq=source_seq,
    )


def record(pos, ref, alts, gt="0/1"):
    return VCFRecord(
        chrom="1",
        pos=pos,
        variant_id="semantic",
        ref=ref,
        alts=tuple(alts),
        gt=gt,
    )


def run(assays, records, sequence="ACGTGCTA"):
    return reconcile(
        assays,
        records,
        reference=SyntheticReference(sequence),
        flank_length=4,
        min_informative_bases=4,
        indel_flank_length=4,
        indel_min_flank_length=4,
    )


def test_shifted_indel_preserves_coordinate_evidence():
    result = run(
        [assay("CGTC[-/A]AAGC", pos=5)],
        [record(6, "A", ["AA"])],
        sequence="CGTCAAGCTA",
    )[0]

    assert result.status == "VCF_ABSENT"
    assert result.coordinate_evidence == "ABSENT"
    assert result.candidate_record_count == 0
    assert result.genotype_evidence == "NOT_EVALUATED"

    assert result.indel_evidence == "UNIQUE"
    assert result.allele_evidence == "MATCH"
    assert result.projection_status == "ASSAY_ALT_PRESENT"
    assert result.assay_alt_dosage == 1


def test_indel_does_not_invoke_snv_validation():
    # Assay and VCF record deliberately share a coordinate,
    # triggering the old coordinate-based SNV validation branch.
    with patch(
        "arrayevidence.reconcile.validate_assay_alleles",
        side_effect=AssertionError(
            "SNV validation must not run for an indel assay"
        ),
    ):
        result = run(
            [assay(pos=4)],
            [record(4, "T", ["TTT"])],
        )[0]

    assert result.indel_evidence == "UNIQUE"
    assert result.allele_evidence == "MATCH"
    assert result.projection_status == "ASSAY_ALT_PRESENT"


def test_ambiguous_indel_does_not_inherit_genotype_projection():
    result = run(
        [assay(pos=4)],
        [
            record(4, "T", ["TTT"], "0/1"),
            record(4, "T", ["TTT"], "0/0"),
        ],
    )[0]

    assert result.coordinate_evidence == "MULTIPLE"
    assert result.genotype_evidence == "CALLED"

    assert result.indel_evidence == "MULTIPLE"
    assert result.allele_evidence == "UNRESOLVED"
    assert result.projection_status == "UNRESOLVED"
    assert result.assay_alt_dosage is None
    assert result.other_alt_present is None
    assert result.ploidy is None
