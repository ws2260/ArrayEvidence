"""Synthetic integration tests for reference-normalized indel reconciliation."""

import pytest

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
        ilmn_id="synthetic_indel",
        name="synthetic",
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
        variant_id="synthetic",
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


@pytest.mark.parametrize(
    "gt,status,dosage",
    [
        ("0/1", "ASSAY_ALT_PRESENT", 1),
        ("1/1", "ASSAY_ALT_PRESENT", 2),
        ("0/0", "ASSAY_REF_ONLY", 0),
        ("./.", "NO_CALL", None),
    ],
)
def test_exact_insertion_match(gt, status, dosage):
    result = run(
        [assay()],
        [record(4, "T", ["TTT"], gt)],
    )[0]

    assert result.indel_evidence == "UNIQUE"
    assert result.allele_evidence == "MATCH"
    assert result.matched_alt_indices == (1,)
    assert result.projection_status == status
    assert result.assay_alt_dosage == dosage


def test_shifted_equivalent_insertion():
    # Reference contains an A homopolymer.
    # Manifest event and VCF representation use different positions.
    result = run(
        [assay("CGTC[-/A]AAGC", pos=5)],
        [record(6, "A", ["AA"], "0/1")],
        sequence="CGTCAAGCTA",
    )[0]

    assert result.indel_evidence == "UNIQUE"
    assert result.allele_evidence == "MATCH"
    assert result.coordinate_evidence == "ABSENT"
    assert result.projection_status == "ASSAY_ALT_PRESENT"


def test_multiallelic_indel_uses_correct_alt():
    result = run(
        [assay()],
        [record(4, "T", ["TA", "TTT"], "0/2")],
    )[0]

    assert result.indel_evidence == "UNIQUE"
    assert result.matched_alt_indices == (2,)
    assert result.projection_status == "ASSAY_ALT_PRESENT"
    assert result.assay_alt_dosage == 1


def test_ambiguous_equivalent_records_not_projected():
    result = run(
        [assay()],
        [
            record(4, "T", ["TTT"], "0/1"),
            record(4, "T", ["TTT"], "0/0"),
        ],
    )[0]

    assert result.indel_evidence == "MULTIPLE"
    assert result.allele_evidence == "UNRESOLVED"
    assert result.projection_status == "UNRESOLVED"
    assert result.assay_alt_dosage is None


def test_absent_equivalent_allele_not_hom_ref():
    result = run(
        [assay()],
        [record(4, "T", ["TA"], "0/0")],
    )[0]

    assert result.indel_evidence == "ABSENT"
    assert result.allele_evidence == "UNRESOLVED"
    assert result.projection_status == "NOT_EVALUATED"
    assert result.assay_alt_dosage is None


def test_unresolved_manifest_not_projected():
    result = run(
        [assay("CCCC[-/TT]GGGG")],
        [record(4, "T", ["TTT"], "0/1")],
    )[0]

    assert result.indel_evidence == "UNRESOLVED"
    assert result.allele_evidence == "UNRESOLVED"
    assert result.projection_status == "NOT_EVALUATED"
    assert result.assay_alt_dosage is None
