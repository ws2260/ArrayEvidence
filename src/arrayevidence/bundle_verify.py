"""Read-only integrity verification for ArrayEvidence TSV bundles."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from dataclasses import fields
import hashlib
import io
import json
from pathlib import Path
import re

from .evidence_schema import (
    SCHEMA_VERSION,
    AlleleCandidateRow,
    AssayEvidenceRow,
)
from .evidence_export import (
    ASSAY_FILENAME,
    BUNDLE_FORMAT_VERSION,
    BUNDLE_MANIFEST_FILENAME,
    CANDIDATE_FILENAME,
)


EXPECTED_FILES = frozenset(
    {
        ASSAY_FILENAME,
        CANDIDATE_FILENAME,
        BUNDLE_MANIFEST_FILENAME,
    }
)

EXPECTED_DATA_FILES = frozenset(
    {
        ASSAY_FILENAME,
        CANDIDATE_FILENAME,
    }
)

SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")

MAX_MANIFEST_BYTES = 65536
MAX_ERRORS = 20


@dataclass(frozen=True)
class BundleVerificationResult:
    valid: bool
    errors: tuple[str, ...]


def _parse_manifest_strict(data: bytes) -> object:
    """Parse UTF-8 JSON, rejecting duplicate keys at every level."""

    def reject_duplicate_keys(pairs):
        result = {}

        for key, value in pairs:
            if key in result:
                raise ValueError(
                    "Duplicate JSON object key is not permitted."
                )
            result[key] = value

        return result

    return json.loads(
        data.decode("utf-8", errors="strict"),
        object_pairs_hook=reject_duplicate_keys,
        parse_constant=lambda value: (_ for _ in ()).throw(
            ValueError("Nonstandard JSON numeric constant.")
        ),
    )


def _canonical_json(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        + "\n"
    ).encode("utf-8")


def _is_nonnegative_integer(value: object) -> bool:
    return type(value) is int and value >= 0


def _check_manifest_structure(
    manifest: object,
    errors: list[str],
) -> bool:
    if not isinstance(manifest, dict):
        errors.append("Manifest must be a JSON object.")
        return False

    expected_keys = {
        "bundle_format_version",
        "evidence_schema_version",
        "files",
    }

    if set(manifest) != expected_keys:
        errors.append("Unexpected or missing manifest fields.")
        return False

    if manifest["bundle_format_version"] != BUNDLE_FORMAT_VERSION:
        errors.append("Unsupported bundle format version.")

    if manifest["evidence_schema_version"] != SCHEMA_VERSION:
        errors.append("Unsupported evidence schema version.")

    entries = manifest["files"]

    if not isinstance(entries, dict):
        errors.append("Manifest files field must be an object.")
        return False

    if set(entries) != EXPECTED_DATA_FILES:
        errors.append("Unexpected or missing manifest file entries.")
        return False

    for filename in sorted(EXPECTED_DATA_FILES):
        entry = entries[filename]

        if not isinstance(entry, dict):
            errors.append(f"Invalid metadata object: {filename}")
            continue

        if set(entry) != {"sha256", "size_bytes", "row_count"}:
            errors.append(f"Invalid metadata fields: {filename}")
            continue

        digest = entry["sha256"]

        if (
            not isinstance(digest, str)
            or SHA256_PATTERN.fullmatch(digest) is None
        ):
            errors.append(f"Invalid SHA-256 metadata: {filename}")

        if not _is_nonnegative_integer(entry["size_bytes"]):
            errors.append(f"Invalid byte-size metadata: {filename}")

        if not _is_nonnegative_integer(entry["row_count"]):
            errors.append(f"Invalid row-count metadata: {filename}")

    return not errors


def _count_tsv_rows(
    data: bytes,
    expected_columns: tuple[str, ...],
) -> int:
    text = data.decode("utf-8", errors="strict")

    if "\x00" in text:
        raise ValueError("TSV contains a NUL character.")

    if '"' in text:
        raise ValueError("TSV contains unsupported quoting.")

    if not text.endswith("\n"):
        raise ValueError("TSV must end with a newline.")

    if "\r" in text:
        raise ValueError("TSV contains carriage returns.")

    reader = csv.reader(
        io.StringIO(text, newline=""),
        delimiter="\t",
        quoting=csv.QUOTE_NONE,
    )

    rows = iter(reader)

    try:
        header = next(rows)
    except StopIteration as exc:
        raise ValueError("TSV header is missing.") from exc

    if tuple(header) != expected_columns:
        raise ValueError("TSV header does not match schema.")

    count = 0

    for row in rows:
        if len(row) != len(expected_columns):
            raise ValueError("TSV row has incorrect column count.")

        count += 1

    return count


def verify_bundle(
    bundle_dir: str | Path,
) -> BundleVerificationResult:
    """Verify bundle structure and integrity without modifying files.

    This is not a cryptographic authentication or authorization check.
    """
    errors: list[str] = []
    bundle = Path(bundle_dir)

    if bundle.is_symlink():
        return BundleVerificationResult(
            False,
            ("Bundle directory must not be a symbolic link.",),
        )

    if not bundle.is_dir():
        return BundleVerificationResult(
            False,
            ("Bundle directory does not exist or is not a directory.",),
        )

    try:
        members = list(bundle.iterdir())
    except OSError:
        return BundleVerificationResult(
            False,
            ("Cannot enumerate bundle directory.",),
        )

    names = {member.name for member in members}

    if names != EXPECTED_FILES:
        errors.append("Bundle has missing or unexpected files.")

    for member in members:
        if member.is_symlink():
            errors.append(f"Symbolic link not permitted: {member.name}")
        elif not member.is_file():
            errors.append(f"Non-file bundle member: {member.name}")

    if errors:
        return BundleVerificationResult(False, tuple(errors[:MAX_ERRORS]))

    manifest_path = bundle / BUNDLE_MANIFEST_FILENAME

    try:
        if manifest_path.stat().st_size > MAX_MANIFEST_BYTES:
            raise ValueError("Manifest exceeds size limit.")

        manifest_bytes = manifest_path.read_bytes()

        if len(manifest_bytes) > MAX_MANIFEST_BYTES:
            raise ValueError("Manifest exceeds size limit.")

        manifest = _parse_manifest_strict(manifest_bytes)

        if manifest_bytes != _canonical_json(manifest):
            raise ValueError("Manifest is not canonical JSON.")

    except (OSError, UnicodeError, ValueError, TypeError):
        return BundleVerificationResult(
            False,
            ("Manifest is missing, malformed, or noncanonical.",),
        )

    if not _check_manifest_structure(manifest, errors):
        return BundleVerificationResult(
            False,
            tuple(errors[:MAX_ERRORS]),
        )

    schemas = {
        ASSAY_FILENAME: tuple(
            field.name for field in fields(AssayEvidenceRow)
        ),
        CANDIDATE_FILENAME: tuple(
            field.name for field in fields(AlleleCandidateRow)
        ),
    }

    for filename in sorted(EXPECTED_DATA_FILES):
        entry = manifest["files"][filename]

        try:
            data = (bundle / filename).read_bytes()

            digest = hashlib.sha256(data).hexdigest()

            if digest != entry["sha256"]:
                errors.append(f"SHA-256 mismatch: {filename}")

            if len(data) != entry["size_bytes"]:
                errors.append(f"Byte-size mismatch: {filename}")

            count = _count_tsv_rows(data, schemas[filename])

            if count != entry["row_count"]:
                errors.append(f"Row-count mismatch: {filename}")

        except (OSError, UnicodeError, ValueError, csv.Error):
            errors.append(f"Cannot validate TSV structure: {filename}")

        if len(errors) >= MAX_ERRORS:
            break

    return BundleVerificationResult(
        valid=not errors,
        errors=tuple(errors[:MAX_ERRORS]),
    )
