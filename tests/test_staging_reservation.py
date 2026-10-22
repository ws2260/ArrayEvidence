"""Session-owned staging reservation tests (synthetic only)."""

import subprocess
import sys
from pathlib import Path

HELPER = Path(__file__).resolve().parents[1] / "scripts" / "reserve_staging.py"


def reserve(path, session, resume=False):
    return subprocess.run(
        [
            sys.executable, str(HELPER),
            str(path), session,
            "true" if resume else "false",
        ],
        capture_output=True, text=True,
    )


def test_fresh_reservation_and_matching_resume(tmp_path):
    stage = tmp_path / "stage"
    assert reserve(stage, "session-A").returncode == 0
    assert stage.is_dir()
    assert reserve(stage, "session-A", resume=True).returncode == 0


def test_fresh_run_cannot_reuse_existing_stage(tmp_path):
    stage = tmp_path / "stage"
    assert reserve(stage, "session-A").returncode == 0
    assert reserve(stage, "session-A").returncode != 0


def test_different_session_cannot_resume(tmp_path):
    stage = tmp_path / "stage"
    assert reserve(stage, "session-A").returncode == 0
    assert reserve(stage, "session-B", resume=True).returncode != 0


def test_unowned_directory_rejected(tmp_path):
    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / "sentinel").write_text("PRESERVE")
    assert reserve(stage, "session-A", resume=True).returncode != 0
    assert (stage / "sentinel").read_text() == "PRESERVE"


def test_symlink_rejected(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    link = tmp_path / "link"
    link.symlink_to(target, target_is_directory=True)
    assert reserve(link, "session-A").returncode != 0
    assert reserve(link, "session-A", resume=True).returncode != 0


def test_missing_parent_rejected(tmp_path):
    stage = tmp_path / "missing" / "stage"
    assert reserve(stage, "session-A").returncode != 0
    assert not stage.exists()


def test_marker_contains_exact_session_id(tmp_path):
    stage = tmp_path / "stage"
    assert reserve(stage, "session-A").returncode == 0

    marker = stage / ".arrayevidence-session"
    assert marker.read_text(encoding="utf-8") == "session-A\n"


def test_matching_resume_rejects_existing_contents(tmp_path):
    stage = tmp_path / "stage"
    assert reserve(stage, "session-A").returncode == 0

    evidence = stage / "validated" / "qc.json"
    evidence.parent.mkdir()
    evidence.write_text('{"synthetic": true}\\n')

    result = reserve(stage, "session-A", resume=True)

    assert result.returncode != 0
    assert "Validated staging output already exists" in result.stderr
    assert evidence.read_text() == '{"synthetic": true}\\n'



def test_concurrent_fresh_reservations(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    stage = tmp_path / "stage"

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(
            lambda i: reserve(stage, f"session-{i}"),
            range(8),
        ))

    assert sum(r.returncode == 0 for r in results) == 1

    marker = stage / ".arrayevidence-session"
    owner = marker.read_text(encoding="utf-8").strip()
    assert owner in {f"session-{i}" for i in range(8)}

    assert reserve(stage, owner, resume=True).returncode == 0


def test_corrupted_marker_rejected(tmp_path):
    stage = tmp_path / "stage"
    assert reserve(stage, "session-A").returncode == 0

    marker = stage / ".arrayevidence-session"
    marker.write_text("session-B\\n")

    assert reserve(stage, "session-A", resume=True).returncode != 0
    assert marker.read_text() == "session-B\\n"



def test_resume_rejects_existing_validated_output(tmp_path):
    stage = tmp_path / "stage"
    assert reserve(stage, "session-A").returncode == 0

    validated = stage / "validated"
    validated.mkdir()
    sentinel = validated / "DO_NOT_OVERWRITE.txt"
    sentinel.write_text("PRESERVE\\n")

    result = reserve(stage, "session-A", resume=True)

    assert result.returncode != 0, (
        "Resume must not silently reuse an existing validated output"
    )
    assert sentinel.read_text() == "PRESERVE\\n"
