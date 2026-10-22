"""Atomically finalize a reserved staging directory."""

import argparse
from pathlib import Path

from publish_evidence import publish


def finalize(bundle, qc, staging_dir, executable, session):
    stage = Path(staging_dir).absolute()

    if stage.is_symlink() or not stage.is_dir():
        raise ValueError("Reserved staging directory is missing or invalid")

    marker = stage / ".arrayevidence-session"
    if marker.is_symlink() or not marker.is_file():
        raise ValueError("Staging ownership marker is missing")

    owner = marker.read_text(encoding="utf-8").strip()
    if not owner:
        raise ValueError("Staging ownership marker is empty")

    if owner != session:
        raise ValueError("Staging belongs to another session")

    destination = stage / "validated"

    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Validated staging output already exists")

    publish(bundle, qc, destination, executable)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle")
    ap.add_argument("qc")
    ap.add_argument("staging_dir")
    ap.add_argument("executable")
    ap.add_argument("session")
    args = ap.parse_args()

    finalize(
        args.bundle,
        args.qc,
        args.staging_dir,
        args.executable,
        args.session,
    )
    print("PASS: Staging finalized atomically")
