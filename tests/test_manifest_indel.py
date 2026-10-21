from arrayevidence.manifest import Assay
from arrayevidence.manifest_indel import reconstruct_manifest_indel


class SyntheticReference:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, reference, start, end):
        if reference != "1":
            raise KeyError(reference)
        return self.sequence[start:end]

    def get_reference_length(self, reference):
        if reference != "1":
            raise KeyError(reference)
        return len(self.sequence)


def make_assay(source_seq, pos, snp="[I/D]"):
    return Assay(
        ilmn_id="synthetic_indel",
        name="synthetic",
        genome_build="38",
        chrom="1",
        map_info=str(pos),
        snp=snp,
        source_seq=source_seq,
    )


def reconstruct(assay, sequence, gap=10):
    return reconstruct_manifest_indel(
        assay,
        SyntheticReference(sequence),
        flank_length=4,
        min_flank_length=4,
        search_radius=20,
        max_reference_gap=gap,
    )


def test_forward_insertion():
    result = reconstruct(
        make_assay("ACGT[-/TT]GCTA", 5),
        "ACGTGCTA",
    )

    assert result.status == "VALIDATED"
    assert result.orientation == "FORWARD"
    assert result.event.status == "INSERTION"
    assert result.event.start == 4
    assert result.event.inserted == "TT"


def test_forward_deletion():
    result = reconstruct(
        make_assay("ACGT[-/TT]GCTA", 5),
        "ACGTTTGCTA",
    )

    assert result.status == "VALIDATED"
    assert result.event.status == "DELETION"
    assert result.event.start == 4
    assert result.event.deleted == "TT"


def test_reverse_insertion():
    result = reconstruct(
        make_assay("TAGC[-/AA]ACGT", 5),
        "ACGTGCTA",
    )

    assert result.status == "VALIDATED"
    assert result.orientation == "REVERSE"
    assert result.event.status == "INSERTION"
    assert result.event.inserted == "TT"


def test_invalid_manifest_tokens():
    result = reconstruct(
        make_assay("ACGT[I/D]GCTA", 5),
        "ACGTGCTA",
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "UNSUPPORTED_ALLELE_TOKENS"


def test_no_flank_alignment():
    result = reconstruct(
        make_assay("CCCC[-/TT]GGGG", 5),
        "ACGTGCTA",
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "NO_FLANK_ALIGNMENT"


def test_not_indel_assay():
    result = reconstruct(
        make_assay("ACGT[-/TT]GCTA", 5, snp="[A/G]"),
        "ACGTGCTA",
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "NOT_INDEL_ASSAY"


def test_unmapped_assay():
    assay = make_assay("ACGT[-/TT]GCTA", 5)
    assay = Assay(
        ilmn_id=assay.ilmn_id,
        name=assay.name,
        genome_build=assay.genome_build,
        chrom="0",
        map_info="0",
        snp=assay.snp,
        source_seq=assay.source_seq,
    )

    result = reconstruct(
        assay,
        "ACGTGCTA",
    )

    assert result.status == "UNRESOLVED"
    assert result.reason == "UNMAPPED_ASSAY"
