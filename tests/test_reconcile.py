from arrayevidence.manifest import Assay
from arrayevidence.vcf import VCFRecord
from arrayevidence.reconcile import reconcile


def assay(name, chrom, pos):
    return Assay(
        ilmn_id=name,
        name=name,
        genome_build="38",
        chrom=chrom,
        map_info=str(pos),
        snp="[A/G]",
        source_seq="AAA[A/G]CCC",
    )


def record(chrom, pos, gt):
    return VCFRecord(
        chrom=chrom,
        pos=pos,
        variant_id=f"synthetic_{pos}",
        ref="A",
        alts=("G",),
        gt=gt,
    )


def test_reconciliation_states():
    assays = [
        assay("a1", "1", 100),
        assay("a2", "1", 100),
        assay("a3", "1", 200),
        assay("a4", "1", 300),
        assay("a5", "0", 0),
    ]

    records = [
        record("1", 100, "0/1"),
        record("1", 200, "./."),
    ]

    results = reconcile(assays, records)

    assert [r.status for r in results] == [
        "VCF_CALLED",
        "VCF_CALLED",
        "VCF_NO_CALL",
        "VCF_ABSENT",
        "UNMAPPED_ASSAY",
    ]

    assert results[0].candidate_record_count == 1
    assert results[1].candidate_record_count == 1
    assert results[3].candidate_record_count == 0


def test_multiple_records_at_same_position():
    assays = [assay("a1", "1", 100)]

    records = [
        record("1", 100, "./."),
        record("1", 100, "0/1"),
    ]

    results = reconcile(assays, records)

    assert results[0].status == "VCF_CALLED"
    assert results[0].candidate_record_count == 2


def test_separate_coordinate_and_genotype_evidence():
    assays = [
        assay("mapped_called", "1", 100),
        assay("mapped_no_call", "1", 200),
        assay("mapped_absent", "1", 300),
        assay("unmapped", "0", 0),
        assay("multiple_candidates", "1", 400),
    ]

    records = [
        record("1", 100, "0/0"),
        record("1", 200, "./."),
        record("1", 400, "./."),
        record("1", 400, "0/1"),
    ]

    results = {
        r.assay_id: r
        for r in reconcile(assays, records)
    }

    assert results["mapped_called"].coordinate_evidence == "PRESENT"
    assert results["mapped_called"].genotype_evidence == "CALLED"

    assert results["mapped_no_call"].coordinate_evidence == "PRESENT"
    assert results["mapped_no_call"].genotype_evidence == "NO_CALL"

    assert results["mapped_absent"].coordinate_evidence == "ABSENT"
    assert results["mapped_absent"].genotype_evidence == "NOT_EVALUATED"

    assert results["unmapped"].coordinate_evidence == "UNMAPPED"
    assert results["unmapped"].genotype_evidence == "NOT_EVALUATED"

    assert results["multiple_candidates"].coordinate_evidence == "MULTIPLE"
    assert results["multiple_candidates"].genotype_evidence == "CALLED"
    assert results["multiple_candidates"].candidate_record_count == 2
