"""Reference-normalized index of simple VCF indel alleles.

The index preserves every VCF record. Equivalent records are not
silently deduplicated, and genotype evidence is not merged.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

from .indel_events import SequenceEvent, minimal_sequence_event
from .indel_normalize import left_normalize_indel
from .reference import ReferenceReader
from .vcf import VCFRecord


@dataclass(frozen=True)
class IndelCandidateResult:
    status: str
    reason: str
    records: tuple[VCFRecord, ...]


class IndelCandidateIndex:
    def __init__(
        self,
        records: Iterable[VCFRecord],
        reference: ReferenceReader,
    ):
        self.reference = reference
        self._by_event = defaultdict(list)
        self.indexed_alleles = 0
        self.skipped_alleles = 0

        for record in records:
            for alt in record.alts:
                event = self._validated_event(record, alt)

                if event is None:
                    self.skipped_alleles += 1
                    continue

                self._by_event[event].append(record)
                self.indexed_alleles += 1

    def _validated_event(
        self,
        record: VCFRecord,
        alt: str,
    ) -> SequenceEvent | None:
        ref = record.ref.upper()
        alt = alt.upper()

        if (
            record.pos < 1
            or not ref
            or not alt
            or not set(ref + alt) <= set("ACGT")
        ):
            return None

        try:
            observed = self.reference.fetch(
                record.chrom,
                record.pos - 1,
                record.pos - 1 + len(ref),
            ).upper()
        except (KeyError, ValueError):
            return None

        if observed != ref:
            return None

        event = minimal_sequence_event(
            record.chrom,
            record.pos,
            ref,
            alt,
        )

        if event.status not in {"INSERTION", "DELETION"}:
            return None

        normalized = left_normalize_indel(
            event,
            self.reference,
        )

        if normalized.status != "NORMALIZED":
            return None

        return normalized.event

    def find(
        self,
        event: SequenceEvent,
    ) -> IndelCandidateResult:
        """Find VCF records equivalent to a validated manifest event.

        The caller must independently establish that the manifest event
        is reference-validated and uses the same reference assembly.
        """

        if event.status not in {"INSERTION", "DELETION"}:
            return IndelCandidateResult(
                "UNRESOLVED",
                "UNSUPPORTED_MANIFEST_EVENT",
                (),
            )

        normalized = left_normalize_indel(
            event,
            self.reference,
        )

        if normalized.status != "NORMALIZED":
            return IndelCandidateResult(
                "UNRESOLVED",
                "MANIFEST_NORMALIZATION_FAILED",
                (),
            )

        candidates = tuple(
            self._by_event.get(normalized.event, ())
        )

        if not candidates:
            return IndelCandidateResult(
                "ABSENT",
                "NO_EQUIVALENT_VCF_ALLELE",
                (),
            )

        if len(candidates) == 1:
            return IndelCandidateResult(
                "UNIQUE",
                "ONE_EQUIVALENT_VCF_RECORD",
                candidates,
            )

        return IndelCandidateResult(
            "MULTIPLE",
            "MULTIPLE_EQUIVALENT_VCF_RECORDS",
            candidates,
        )
