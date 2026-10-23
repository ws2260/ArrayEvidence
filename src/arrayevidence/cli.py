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

    qc_parser = commands.add_parser(
        "qc",
        help="Summarize a verified evidence bundle.",
        description=(
            "Produce deterministic aggregate QC JSON "
            "without modifying the evidence bundle."
        ),
    )

    qc_parser.add_argument(
        "bundle",
        help="Path to the evidence bundle directory.",
    )
    qc_parser.add_argument(
        "--output",
        help="Write QC JSON to a new file instead of stdout.",
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

        stage = "BUILD_ERROR"

        try:
            stage = "MANIFEST_READ_ERROR"
            assays = list(read_assays(args.manifest))

            stage = "VCF_READ_ERROR"
            records = list(read_vcf(args.vcf))

            stage = "REFERENCE_ERROR"
            with PysamReferenceReader(args.reference) as reference:
                stage = "RECONCILIATION_ERROR"
                results = reconcile(
                    assays,
                    records,
                    reference=reference,
                )
                stage = "REFERENCE_ERROR"

            stage = "EVIDENCE_ERROR"
            assay_rows, candidate_rows = build_evidence_rows(
                args.sample_id,
                results,
            )

            stage = "EXPORT_ERROR"
            export_evidence_tsv(
                args.output,
                assay_rows,
                candidate_rows,
            )

        except Exception:
            print(
                f"FAIL [{stage}]: Evidence bundle construction failed.",
                file=sys.stderr,
            )
            return 1

        print("PASS: Evidence bundle constructed.")
        return 0

    if args.command == "qc":
        from . import qc

        try:
            summary = qc.summarize_bundle(args.bundle)
            rendered = qc.render_qc_json(summary)
        except Exception:
            print(
                "FAIL [QC_ERROR]: Evidence QC failed.",
                file=sys.stderr,
            )
            return 1

        if args.output is None:
            sys.stdout.write(rendered)
            return 0

        from pathlib import Path
        import os
        import tempfile

        destination = Path(args.output)

        try:
            bundle_path = Path(args.bundle).resolve(strict=True)
            output_path = destination.resolve(strict=False)

            if output_path == bundle_path or bundle_path in output_path.parents:
                raise ValueError("QC output must be outside the evidence bundle")
        except (OSError, RuntimeError, ValueError):
            print(
                "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.",
                file=sys.stderr,
            )
            return 1

        temporary = None
        publication_failed = False
        cleanup_failed = False

        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="",
                dir=destination.parent,
                prefix=".arrayevidence-qc-",
                suffix=".tmp",
                delete=False,
            ) as handle:
                temporary = Path(handle.name)
                handle.write(rendered)

            os.link(temporary, destination)
        except OSError:
            publication_failed = True
        finally:
            if temporary is not None:
                try:
                    temporary.unlink(missing_ok=True)
                except OSError:
                    cleanup_failed = True

        if publication_failed:
            print(
                "FAIL [QC_OUTPUT_ERROR]: Evidence QC output failed.",
                file=sys.stderr,
            )
            return 1

        if cleanup_failed:
            print(
                "FAIL [QC_OUTPUT_ERROR]: Evidence QC output cleanup failed.",
                file=sys.stderr,
            )
            return 1

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
