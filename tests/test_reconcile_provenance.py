"""Synthetic provenance contract for reconciled assay alleles."""

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


def assay(*, snp, source_seq, pos=5):
    return Assay(
        ilmn_id="synthetic_assay",
        name="synthetic_assay",
        genome_build="38",
        chrom="1",
        map_info=str(pos),
        snp=snp,
        source_seq=source_seq,
    )


def record(*, pos, ref, alts, gt="0/1", variant_id="."):
    return VCFRecord(
        chrom="1",
        pos=pos,
        variant_id=variant_id,
        ref=ref,
        alts=tuple(alts),
        gt=gt,
    )


def run(assays, records, reference):
    return reconcile(
        assays,
        records,
        reference=reference,
        flank_length=4,
        min_informative_bases=4,
        indel_flank_length=4,
        indel_min_flank_length=4,
    )


def test_biallelic_snv_match_preserves_source_alt_one():
    reference = SyntheticReference("ACGTGCTA")
    a = assay(snp="[G/A]", source_seq="ACGT[G/A]CTA")
    r = record(pos=5, ref="G", alts=("A",), variant_id="snv1")

    result = run([a], [r], reference)[0]

    assert result.allele_evidence == "MATCH"
    assert len(result.allele_sources) == 1
    assert result.allele_sources[0].record == r
    assert result.allele_sources[0].alt_index == 1


def test_multiallelic_snv_preserves_matched_alt_index():
    reference = SyntheticReference("ACGTGCTA")
    a = assay(snp="[G/A]", source_seq="ACGT[G/A]CTA")
    r = record(
        pos=5,
        ref="G",
        alts=("C", "A"),
        gt="0/2",
        variant_id="multi1",
    )

    result = run([a], [r], reference)[0]

    assert result.allele_evidence == "PARTIAL_MATCH"
    assert len(result.allele_sources) == 1
    assert result.allele_sources[0].record == r
    assert result.allele_sources[0].alt_index == 2


def test_shifted_indel_preserves_actual_vcf_record():
    reference = SyntheticReference("CGTCAAGCTA")
    a = assay(
        snp="[I/D]",
        source_seq="CGTC[-/A]AAGC",
        pos=5,
    )
    r = record(
        pos=6,
        ref="A",
        alts=("AA",),
        gt="0/1",
        variant_id="shifted1",
    )

    result = run([a], [r], reference)[0]

    assert result.coordinate_evidence == "ABSENT"
    assert result.indel_evidence == "UNIQUE"
    assert result.allele_evidence == "MATCH"
    assert len(result.allele_sources) == 1
    assert result.allele_sources[0].record == r
    assert result.allele_sources[0].alt_index == 1


def test_ambiguous_indels_preserve_all_sources_without_projection():
    reference = SyntheticReference("ACGTGCTA")
    a = assay(
        snp="[I/D]",
        source_seq="ACGT[-/TT]GCTA",
        pos=5,
    )
    r1 = record(pos=4, ref="T", alts=("TTT",), variant_id="indel1")
    r2 = record(pos=4, ref="T", alts=("TTT",), variant_id="indel2")

    result = run([a], [r1, r2], reference)[0]

    assert result.indel_evidence == "MULTIPLE"
    assert result.allele_evidence == "UNRESOLVED"
    assert result.projection_status == "UNRESOLVED"
    assert [
        (source.record.variant_id, source.alt_index)
        for source in result.allele_sources
    ] == [("indel1", 1), ("indel2", 1)]


def test_absent_indel_has_no_allele_sources():
    reference = SyntheticReference("ACGTGCTA")
    a = assay(
        snp="[I/D]",
        source_seq="ACGT[-/TT]GCTA",
        pos=5,
    )

    result = run([a], [], reference)[0]

    assert result.indel_evidence == "ABSENT"
    assert result.allele_sources == ()


def test_snv_mismatch_has_no_confirmed_allele_sources():
    reference = SyntheticReference("ACGTGCTA")
    a = assay(snp="[G/A]", source_seq="ACGT[G/A]CTA")
    r = record(pos=5, ref="G", alts=("C",), variant_id="mismatch1")

    result = run([a], [r], reference)[0]

    assert result.allele_evidence == "MISMATCH"
    assert result.allele_sources == ()
