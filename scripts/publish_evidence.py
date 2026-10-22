"""Atomic, no-overwrite publication on supported Linux filesystems."""
import argparse
import ctypes
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

libc = ctypes.CDLL(None, use_errno=True)
renameat2 = libc.renameat2
renameat2.argtypes = [
    ctypes.c_int, ctypes.c_char_p,
    ctypes.c_int, ctypes.c_char_p, ctypes.c_uint,
]
renameat2.restype = ctypes.c_int

def publish(bundle, qc, destination, executable):
    dest = Path(destination).absolute()
    if not dest.parent.is_dir():
        raise ValueError("Destination parent does not exist")
    if dest.exists() or dest.is_symlink():
        raise FileExistsError("Destination already exists")

    stage = Path(tempfile.mkdtemp(
        prefix=".arrayevidence-stage-", dir=dest.parent
    ))
    try:
        shutil.copytree(bundle, stage / "bundle")
        shutil.copy2(qc, stage / "qc.json")
        json.loads((stage / "qc.json").read_text())

        subprocess.run(
            [executable, "verify", str(stage / "bundle")],
            check=True, capture_output=True, text=True,
        )

        rc = renameat2(
            -100, os.fsencode(stage),
            -100, os.fsencode(dest), 1
        )
        if rc != 0:
            err = ctypes.get_errno()
            raise OSError(err, os.strerror(err), str(dest))
    finally:
        if stage.exists():
            shutil.rmtree(stage)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle")
    ap.add_argument("qc")
    ap.add_argument("destination")
    ap.add_argument("executable")
    a = ap.parse_args()
    publish(a.bundle, a.qc, a.destination, a.executable)
