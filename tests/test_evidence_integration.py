"""End-to-end synthetic reconciliation and evidence-row integration."""

from arrayevidence.evidence_builder import build_evidence_rows
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
            raise ValueError("Reference interval out of range")
        return self.sequence[start:end]

    def get_reference_length(self, chrom):
        if chrom != "1":
            raise KeyError(chrom)
        return len(self.sequence)


def assay(name, pos, snp, source_seq):
    return Assay(
        ilmn_id=name,
        name=name,
        genome_build="38",
        chrom="1",
        map_info=str(pos),
        snp=snp,
        source_seq=source_seq,
    )


def record(pos, ref, alts, gt="0/1", variant_id="."):
    return VCFRecord(
        chrom="1",
        pos=pos,
        variant_id=variant_id,
        ref=ref,
        alts=tuple(alts),
        gt=gt,
    )


def pipeline(assays, records, sequence, sample_id="synthetic_sample"):
    reconciled = reconcile(
        assays,
        records,
        reference=SyntheticReference(sequence),
        flank_length=4,
        min_informative_bases=4,
        indel_flank_length=4,
        indel_min_flank_length=4,
    )
    return build_evidence_rows(sample_id, reconciled)


def test_shifted_indel_survives_full_pipeline():
    a = assay(
        "shifted_assay",
        5,
        "[I/D]",
        "CGTC[-/A]AAGC",
    )
    r = record(
        6,
        "A",
        ("AA",),
        variant_id="shifted_vcf",
    )

    assays, candidates = pipeline(
        [a], [r], "CGTCAAGCTA"
    )

    assert len(assays) == 1
    assert len(candidates) == 1

    evidence = assays[0]
    source = candidates[0]

    assert evidence.coordinate_status == "VCF_ABSENT"
    assert evidence.coordinate_evidence == "ABSENT"
    assert evidence.coordinate_candidate_count == 0
    assert evidence.indel_evidence == "UNIQUE"
    assert evidence.allele_evidence == "MATCH"
    assert evidence.projection_status == "ASSAY_ALT_PRESENT"
    assert evidence.assay_alt_dosage == 1

    assert source.sample_id == "synthetic_sample"
    assert source.assay_id == "shifted_assay"
    assert source.vcf_id == "shifted_vcf"
    assert source.vcf_pos == 6
    assert source.vcf_ref == "A"
    assert source.vcf_alt == "AA"
    assert source.vcf_alt_index == 1
    assert source.vcf_gt == "0/1"
    assert source.candidate_status == "UNIQUE"


def test_ambiguous_indel_keeps_both_sources():
    a = assay(
        "ambiguous_assay",
        5,
        "[I/D]",
        "ACGT[-/TT]GCTA",
    )
    r1 = record(
        4, "T", ("TTT",),
        variant_id="first",
    )
    r2 = record(
        4, "T", ("TTT",),
        variant_id="second",
    )

    assays, candidates = pipeline(
        [a], [r1, r2], "ACGTGCTA"
    )

    assert len(assays) == 1
    assert assays[0].indel_evidence == "MULTIPLE"
    assert assays[0].projection_status == "UNRESOLVED"
    assert assays[0].assay_alt_dosage is None

    assert len(candidates) == 2
    assert [c.vcf_id for c in candidates] == [
        "first", "second"
    ]
    assert [c.candidate_index for c in candidates] == [1, 2]
    assert all(
        c.candidate_status == "MULTIPLE"
        for c in candidates
    )


def test_multiallelic_snv_preserves_alt_two():
    a = assay(
        "multiallelic_assay",
        5,
        "[G/A]",
        "ACGT[G/A]CTA",
    )
    r = record(
        5,
        "G",
        ("C", "A"),
        gt="0|2",
        variant_id="multiallelic_vcf",
    )

    assays, candidates = pipeline(
        [a], [r], "ACGTGCTA"
    )

    assert len(assays) == 1
    assert assays[0].allele_evidence == "PARTIAL_MATCH"
    assert assays[0].projection_status == "ASSAY_ALT_PRESENT"
    assert assays[0].assay_alt_dosage == 1

    assert len(candidates) == 1
    assert candidates[0].vcf_alt_index == 2
    assert candidates[0].vcf_alt == "A"
    assert candidates[0].vcf_gt == "0|2"
    assert candidates[0].candidate_status == "PARTIAL_MATCH"


def test_absent_indel_does_not_become_hom_ref():
    a = assay(
        "absent_assay",
        5,
        "[I/D]",
        "ACGT[-/TT]GCTA",
    )

    assays, candidates = pipeline(
        [a], [], "ACGTGCTA"
    )

    assert len(assays) == 1
    assert assays[0].indel_evidence == "ABSENT"
    assert assays[0].assay_alt_dosage is None
    assert assays[0].projection_status == "NOT_EVALUATED"
    assert candidates == []
