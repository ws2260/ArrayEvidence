"""Command-line interface for ArrayEvidence."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import sys

from .bundle_verify import verify_bundle
from .manifest import read_assays
from .vcf import read_vcf
from .fasta_adapter import PysamReferenceReader
from .reconcile import reconcile
from .evidence_builder import build_evidence_rows
from .evidence_export import export_evidence_tsv


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="arrayevidence",
        description=(
            "Deterministic technical evidence for "
            "genotyping arrays."
        ),
    )

    commands = parser.add_subparsers(
        dest="command",
        required=True,
    )

    verify_parser = commands.add_parser(
        "verify",
        help="Verify the integrity of an evidence bundle.",
        description=(
            "Check bundle structure, checksums, and "
            "TSV integrity without modifying the bundle."
        ),
    )

    verify_parser.add_argument(
        "bundle",
        help="Path to the evidence bundle directory.",
    )

    build_parser = commands.add_parser(
        "build",
        help="Build an immutable technical evidence bundle.",
        description=(
            "Reconcile an assay manifest against a single-sample "
            "VCF and indexed reference FASTA."
        ),
    )

    build_parser.add_argument(
        "--manifest",
        required=True,
        help="Path to the assay manifest.",
    )
    build_parser.add_argument(
        "--vcf",
        required=True,
        help="Path to the single-sample VCF.",
    )
    build_parser.add_argument(
        "--reference",
        required=True,
        help="Path to the indexed reference FASTA.",
    )
    build_parser.add_argument(
        "--sample-id",
        required=True,
        help="Operator-supplied pseudonymous sample identifier.",
    )
    build_parser.add_argument(
        "--output",
        required=True,
        help="New evidence bundle directory; must not exist.",
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the CLI and return a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "build":
        from pathlib import Path

        if not args.sample_id.strip():
            print(
                "ERROR: --sample-id must be non-empty.",
                file=sys.stderr,
            )
            return 2

        destination = Path(args.output)

        if destination.exists() or destination.is_symlink():
            print(
                "FAIL: Evidence bundle destination already exists.",
                file=sys.stderr,
            )
            return 1

        try:
            assays = list(read_assays(args.manifest))
            records = list(read_vcf(args.vcf))

            with PysamReferenceReader(args.reference) as reference:
                results = reconcile(
                    assays,
                    records,
                    reference=reference,
                )

            assay_rows, candidate_rows = build_evidence_rows(
                args.sample_id,
                results,
            )

            export_evidence_tsv(
                args.output,
                assay_rows,
                candidate_rows,
            )

        except Exception:
            print(
                "FAIL: Evidence bundle construction failed.",
                file=sys.stderr,
            )
            return 1

        print("PASS: Evidence bundle constructed.")
        return 0

    if args.command == "verify":
        result = verify_bundle(args.bundle)

        if result.valid:
            print("PASS: Bundle verification succeeded.")
            return 0

        print(
            "FAIL: Bundle verification failed.",
            file=sys.stderr,
        )

        for error in result.errors:
            print(f"  - {error}", file=sys.stderr)

        return 1

    parser.error("Unsupported command.")
    return 2
