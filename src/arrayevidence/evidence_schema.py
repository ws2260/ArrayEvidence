"""Versioned, immutable technical evidence row models.

These models describe evidence; they do not perform allele matching,
genotype calling, clinical annotation, or data export.

No missing genotype or allele evidence is inferred.
"""

from dataclasses import dataclass


SCHEMA_VERSION = "1.0.0"


def _require_identifier(value: str, field: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a non-empty string")


def _require_version(value: str) -> None:
    if value != SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported schema version: {value!r}; "
            f"expected {SCHEMA_VERSION!r}"
        )


def _require_positive_integer(value: int, field: str) -> None:
    if type(value) is not int or value < 1:
        raise ValueError(f"{field} must be a positive integer")


@dataclass(frozen=True)
class AssayEvidenceRow:
    """One technical evidence row for one sample and one assay."""

    schema_version: str
    sample_id: str
    assay_id: str

    manifest_chrom: str | None
    manifest_pos: int | None

    coordinate_status: str
    coordinate_evidence: str
    coordinate_candidate_count: int
    coordinate_genotype_evidence: str

    allele_evidence: str
    allele_reason: str
    indel_evidence: str

    projection_status: str
    projection_reason: str
    assay_alt_dosage: int | None
    other_alt_present: bool | None
    ploidy: int | None

    def __post_init__(self) -> None:
        _require_version(self.schema_version)
        _require_identifier(self.sample_id, "sample_id")
        _require_identifier(self.assay_id, "assay_id")

        if self.manifest_chrom is not None:
            _require_identifier(self.manifest_chrom, "manifest_chrom")

        if self.manifest_pos is not None:
            _require_positive_integer(
                self.manifest_pos, "manifest_pos"
            )

        if (self.manifest_chrom is None) != (
            self.manifest_pos is None
        ):
            raise ValueError(
                "manifest_chrom and manifest_pos must both be "
                "present or both be None"
            )

        if (
            type(self.coordinate_candidate_count) is not int
            or self.coordinate_candidate_count < 0
        ):
            raise ValueError(
                "coordinate_candidate_count must be "
                "a non-negative integer"
            )

        if self.assay_alt_dosage is not None:
            if (
                type(self.assay_alt_dosage) is not int
                or self.assay_alt_dosage < 0
            ):
                raise ValueError(
                    "assay_alt_dosage must be a "
                    "non-negative integer or None"
                )

        if self.ploidy is not None:
            _require_positive_integer(self.ploidy, "ploidy")

        if (
            self.assay_alt_dosage is not None
            and self.ploidy is not None
            and self.assay_alt_dosage > self.ploidy
        ):
            raise ValueError(
                "assay_alt_dosage cannot exceed ploidy"
            )

        if (
            self.other_alt_present is not None
            and type(self.other_alt_present) is not bool
        ):
            raise ValueError(
                "other_alt_present must be bool or None"
            )

        evaluable_statuses = {
            "ASSAY_ALT_PRESENT",
            "ASSAY_REF_ONLY",
            "OTHER_ALT_PRESENT",
        }
        non_evaluable_statuses = {
            "NO_CALL",
            "NOT_EVALUATED",
            "UNRESOLVED",
        }

        if self.projection_status in evaluable_statuses:
            if (
                self.assay_alt_dosage is None
                or self.other_alt_present is None
                or self.ploidy is None
            ):
                raise ValueError(
                    "Evaluable projection requires dosage, "
                    "other_alt_present, and ploidy"
                )

            if self.projection_status == "ASSAY_ALT_PRESENT":
                if (
                    self.assay_alt_dosage < 1
                    or self.other_alt_present is not False
                ):
                    raise ValueError(
                        "ASSAY_ALT_PRESENT requires positive dosage "
                        "and no other ALT"
                    )

            elif self.projection_status == "ASSAY_REF_ONLY":
                if (
                    self.assay_alt_dosage != 0
                    or self.other_alt_present is not False
                ):
                    raise ValueError(
                        "ASSAY_REF_ONLY requires zero dosage "
                        "and no other ALT"
                    )

            elif self.projection_status == "OTHER_ALT_PRESENT":
                if self.other_alt_present is not True:
                    raise ValueError(
                        "OTHER_ALT_PRESENT requires another ALT"
                    )

            if self.allele_evidence not in {
                "MATCH",
                "PARTIAL_MATCH",
            }:
                raise ValueError(
                    "Evaluable projection requires established "
                    "allele evidence"
                )

            if self.indel_evidence in {
                "ABSENT",
                "UNRESOLVED",
                "MULTIPLE",
            }:
                raise ValueError(
                    "Unresolved indel evidence cannot support "
                    "genotype projection"
                )

        elif self.projection_status in non_evaluable_statuses:
            if (
                self.assay_alt_dosage is not None
                or self.other_alt_present is not None
                or self.ploidy is not None
            ):
                raise ValueError(
                    "Non-evaluable projection must not contain "
                    "dosage, other ALT, or ploidy values"
                )

        else:
            raise ValueError(
                f"Unsupported projection status: "
                f"{self.projection_status!r}"
            )


@dataclass(frozen=True)
class AlleleCandidateRow:
    """One source VCF ALT allele candidate for one sample and assay."""

    schema_version: str
    sample_id: str
    assay_id: str

    candidate_index: int

    vcf_chrom: str
    vcf_pos: int
    vcf_id: str
    vcf_ref: str
    vcf_alt: str
    vcf_alt_index: int
    vcf_gt: str

    candidate_status: str

    def __post_init__(self) -> None:
        _require_version(self.schema_version)
        _require_identifier(self.sample_id, "sample_id")
        _require_identifier(self.assay_id, "assay_id")

        _require_positive_integer(
            self.candidate_index, "candidate_index"
        )
        _require_identifier(self.vcf_chrom, "vcf_chrom")
        _require_positive_integer(self.vcf_pos, "vcf_pos")
        _require_identifier(self.vcf_ref, "vcf_ref")
        _require_identifier(self.vcf_alt, "vcf_alt")

        _require_positive_integer(
            self.vcf_alt_index, "vcf_alt_index"
        )

        if not isinstance(self.vcf_id, str):
            raise ValueError("vcf_id must be a string")

        if not isinstance(self.vcf_gt, str):
            raise ValueError("vcf_gt must be a string")

        _require_identifier(
            self.candidate_status, "candidate_status"
        )
