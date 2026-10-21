"""Command-line interface for ArrayEvidence."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import sys

from .bundle_verify import verify_bundle


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

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the CLI and return a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

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
