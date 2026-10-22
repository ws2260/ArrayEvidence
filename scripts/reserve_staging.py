"""Reserve a staging directory for a Nextflow session."""

import argparse
import os
from pathlib import Path

MARKER = ".arrayevidence-session"


def reserve(destination, session, resume):
    path = Path(destination).absolute()
    parent = path.parent

    if not session or "/" in session or "\\x00" in session:
        raise ValueError("Invalid session identifier")

    if not parent.is_dir() or parent.is_symlink():
        raise ValueError("Invalid staging parent")

    if resume:
        if path.is_symlink() or not path.is_dir():
            raise ValueError("No valid staging directory to resume")

        marker = path / MARKER
        if marker.is_symlink() or not marker.is_file():
            raise ValueError("Staging ownership marker missing")

        if marker.read_text(encoding="utf-8").strip() != session:
            raise ValueError("Staging belongs to another session")

        # Fail closed: publishDir overwrite:false silently skips existing
        # validated outputs, even when resumed tasks recompute successfully.
        validated = path / "validated"
        if validated.exists() or validated.is_symlink():
            raise ValueError(
                "Validated staging output already exists; "
                "resume cannot safely reuse this destination"
            )

        return

    # Atomic reservation: fails if the destination already exists.
    path.mkdir(mode=0o700)

    marker = path / MARKER
    try:
        with marker.open("x", encoding="utf-8") as handle:
            handle.write(session + "\n")
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        # Only remove our own empty reservation, never user data.
        marker.unlink(missing_ok=True)
        try:
            path.rmdir()
        except OSError:
            pass
        raise


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("destination")
    ap.add_argument("session")
    ap.add_argument("resume", choices=("true", "false"))
    args = ap.parse_args()

    reserve(args.destination, args.session, args.resume == "true")
    print("PASS: Staging reservation accepted")
