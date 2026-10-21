"""Synthetic deterministic TSV evidence export contract."""

import csv
from dataclasses import fields

import pytest

from arrayevidence.evidence_export import export_evidence_tsv
from arrayevidence.evidence_schema import (
    SCHEMA_VERSION,
    AlleleCandidateRow,
    AssayEvidenceRow,
)


def assay_row(sample_id="sample_B", assay_id="assay_B", **overrides):
    values = dict(
        schema_version=SCHEMA_VERSION,
        sample_id=sample_id,
        assay_id=assay_id,
        manifest_chrom="1",
        manifest_pos=5,
        coordinate_status="VCF_ABSENT",
        coordinate_evidence="ABSENT",
        coordinate_candidate_count=0,
        coordinate_genotype_evidence="NOT_EVALUATED",
        allele_evidence="MATCH",
        allele_reason="REFERENCE_NORMALIZED_INDEL_MATCH",
        indel_evidence="UNIQUE",
        projection_status="ASSAY_ALT_PRESENT",
        projection_reason="MATCHED_ALT_PROJECTED",
        assay_alt_dosage=1,
        other_alt_present=False,
        ploidy=2,
    )
    values.update(overrides)
    return AssayEvidenceRow(**values)


def candidate_row(
    sample_id="sample_B",
    assay_id="assay_B",
    candidate_index=1,
    **overrides,
):
    values = dict(
        schema_version=SCHEMA_VERSION,
        sample_id=sample_id,
        assay_id=assay_id,
        candidate_index=candidate_index,
        vcf_chrom="1",
        vcf_pos=6,
        vcf_id="synthetic_variant",
        vcf_ref="A",
        vcf_alt="AA",
        vcf_alt_index=1,
        vcf_gt="0/1",
        candidate_status="UNIQUE",
    )
    values.update(overrides)
    return AlleleCandidateRow(**values)


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_export_creates_two_versioned_tables(tmp_path):
    export_evidence_tsv(
        (tmp_path / "bundle"),
        [assay_row()],
        [candidate_row()],
    )

    assay_path = tmp_path / "bundle" / "assay_evidence.tsv"
    candidate_path = tmp_path / "bundle" / "allele_candidates.tsv"

    assert assay_path.is_file()
    assert candidate_path.is_file()

    assays = read_tsv(assay_path)
    candidates = read_tsv(candidate_path)

    assert len(assays) == 1
    assert len(candidates) == 1

    assert assays[0]["schema_version"] == SCHEMA_VERSION
    assert candidates[0]["schema_version"] == SCHEMA_VERSION
    assert assays[0]["sample_id"] == "sample_B"
    assert candidates[0]["vcf_pos"] == "6"
    assert candidates[0]["vcf_gt"] == "0/1"


def test_headers_follow_schema_field_order(tmp_path):
    export_evidence_tsv((tmp_path / "bundle"), [], [])

    assay_header = (
        tmp_path / "bundle" / "assay_evidence.tsv"
    ).read_text(encoding="utf-8").splitlines()[0]

    candidate_header = (
        tmp_path / "bundle" / "allele_candidates.tsv"
    ).read_text(encoding="utf-8").splitlines()[0]

    assert assay_header.split("\t") == [
        field.name for field in fields(AssayEvidenceRow)
    ]
    assert candidate_header.split("\t") == [
        field.name for field in fields(AlleleCandidateRow)
    ]


def test_null_zero_and_false_remain_distinct(tmp_path):
    absent = assay_row(
        allele_evidence="UNRESOLVED",
        allele_reason="NO_EQUIVALENT_VCF_ALLELE",
        indel_evidence="ABSENT",
        projection_status="NOT_EVALUATED",
        projection_reason="NOT_EVALUATED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )

    ref_only = assay_row(
        sample_id="sample_A",
        assay_id="assay_A",
        projection_status="ASSAY_REF_ONLY",
        assay_alt_dosage=0,
        other_alt_present=False,
    )

    export_evidence_tsv(
        (tmp_path / "bundle"),
        [absent, ref_only],
        [candidate_row(sample_id="sample_A", assay_id="assay_A")],
    )

    rows = read_tsv(tmp_path / "bundle" / "assay_evidence.tsv")
    by_sample = {row["sample_id"]: row for row in rows}

    assert by_sample["sample_B"]["assay_alt_dosage"] == "NULL"
    assert by_sample["sample_B"]["other_alt_present"] == "NULL"
    assert by_sample["sample_A"]["assay_alt_dosage"] == "0"
    assert by_sample["sample_A"]["other_alt_present"] == "False"


def test_unmapped_coordinates_use_null(tmp_path):
    unmapped = assay_row(
        manifest_chrom=None,
        manifest_pos=None,
        coordinate_status="UNMAPPED_ASSAY",
        coordinate_evidence="UNMAPPED",
        allele_evidence="NOT_EVALUATED",
        allele_reason="UNMAPPED_ASSAY",
        indel_evidence="NOT_EVALUATED",
        projection_status="NOT_EVALUATED",
        projection_reason="NOT_EVALUATED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )

    export_evidence_tsv((tmp_path / "bundle"), [unmapped], [])

    row = read_tsv(tmp_path / "bundle" / "assay_evidence.tsv")[0]

    assert row["manifest_chrom"] == "NULL"
    assert row["manifest_pos"] == "NULL"


def test_export_is_deterministic_across_input_order(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"

    assays = [
        assay_row("sample_B", "assay_B"),
        assay_row("sample_A", "assay_A"),
    ]
    candidates = [
        candidate_row("sample_B", "assay_B"),
        candidate_row("sample_A", "assay_A"),
    ]

    export_evidence_tsv(first, assays, candidates)
    export_evidence_tsv(
        second,
        list(reversed(assays)),
        list(reversed(candidates)),
    )

    for filename in ("assay_evidence.tsv", "allele_candidates.tsv"):
        assert (first / filename).read_bytes() == (
            second / filename
        ).read_bytes()


def test_candidate_index_controls_order_within_assay(tmp_path):
    candidates = [
        candidate_row(candidate_index=2, vcf_id="second"),
        candidate_row(candidate_index=1, vcf_id="first"),
    ]

    ambiguous = assay_row(
        allele_evidence="UNRESOLVED",
        allele_reason="MULTIPLE_EQUIVALENT_VCF_ALLELES",
        indel_evidence="MULTIPLE",
        projection_status="UNRESOLVED",
        projection_reason="UNIQUE_INDEL_ALLELE_NOT_ESTABLISHED",
        assay_alt_dosage=None,
        other_alt_present=None,
        ploidy=None,
    )
    candidates = [
        candidate_row(
            candidate_index=row.candidate_index,
            vcf_id=row.vcf_id,
            candidate_status="MULTIPLE",
        )
        for row in candidates
    ]

    export_evidence_tsv(
        (tmp_path / "bundle"),
        [ambiguous],
        candidates,
    )

    rows = read_tsv(tmp_path / "bundle" / "allele_candidates.tsv")

    assert [row["vcf_id"] for row in rows] == [
        "first", "second"
    ]


@pytest.mark.parametrize("bad_value", [
    "embedded\ttab",
    "embedded\nnewline",
    "embedded\rcarriage_return",
])
def test_unsafe_tsv_values_rejected_before_files_created(
    tmp_path, bad_value
):
    invalid = assay_row(assay_id=bad_value)

    with pytest.raises(ValueError):
        export_evidence_tsv((tmp_path / "bundle"), [invalid], [])

    assert not (tmp_path / "bundle" / "assay_evidence.tsv").exists()
    assert not (tmp_path / "bundle" / "allele_candidates.tsv").exists()


def test_unsafe_candidate_value_rejected_before_files_created(tmp_path):
    invalid = candidate_row(vcf_id="bad\tidentifier")

    with pytest.raises(ValueError):
        export_evidence_tsv(
            (tmp_path / "bundle"),
            [assay_row()],
            [invalid],
        )

    assert not (tmp_path / "bundle" / "assay_evidence.tsv").exists()
    assert not (tmp_path / "bundle" / "allele_candidates.tsv").exists()


def test_utf8_and_unix_line_endings(tmp_path):
    export_evidence_tsv(
        (tmp_path / "bundle"),
        [assay_row(assay_id="synthetic_α")],
        [candidate_row(assay_id="synthetic_α")],
    )

    data = (tmp_path / "bundle" / "assay_evidence.tsv").read_bytes()

    assert "synthetic_α".encode("utf-8") in data
    assert b"\r\n" not in data
    assert data.endswith(b"\n")
