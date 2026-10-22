"""Read-only aggregate QC for verified ArrayEvidence bundles.

This module intentionally excludes sample identifiers, assay identifiers,
genomic coordinates, genotypes, and input/output paths from its reports.
"""

from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path

from .bundle_verify import verify_bundle


ASSAY_CATEGORIES = (
    "coordinate_status",
    "coordinate_genotype_evidence",
    "allele_evidence",
    "allele_reason",
    "indel_evidence",
    "projection_status",
    "projection_reason",
)

CANDIDATE_CATEGORIES = (
    "candidate_status",
)


def render_qc_json(summary: dict) -> str:
    """Render a deterministic, strict JSON document."""
    return json.dumps(
        summary,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _count_table(
    path: Path,
    categories: tuple[str, ...],
) -> tuple[int, dict[str, dict[str, int]]]:
    counts = {
        field: Counter()
        for field in categories
    }

    total = 0

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not set(categories).issubset(reader.fieldnames or []):
            raise ValueError("Unexpected evidence table schema")

        for row in reader:
            total += 1

            for field in categories:
                counts[field][row[field]] += 1

    return total, {
        field: dict(sorted(counter.items()))
        for field, counter in counts.items()
    }


def summarize_bundle(bundle_dir: str | Path) -> dict:
    """Summarize one verified bundle without modifying it."""
    bundle = Path(bundle_dir)

    verification = verify_bundle(bundle)

    if not verification.valid:
        raise ValueError("Bundle verification failed")

    assay_total, assay_counts = _count_table(
        bundle / "assay_evidence.tsv",
        ASSAY_CATEGORIES,
    )

    candidate_total, candidate_counts = _count_table(
        bundle / "allele_candidates.tsv",
        CANDIDATE_CATEGORIES,
    )

    return {
        "assay_rows": assay_total,
        "candidate_rows": candidate_total,
        **assay_counts,
        **candidate_counts,
    }
