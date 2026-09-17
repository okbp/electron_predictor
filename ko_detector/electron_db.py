"""Electron donor / acceptor categories per genome from the curated database TSV.

The database (``data/electron/electron_donor_acceptor_DATABASE_with_genome.tsv``) has one row per
genome x role x compound, with ``donor_category`` filled on ``electron_donor`` rows and
``acceptor_category`` on ``electron_acceptor`` rows. It is saved by Excel, so it may be cp932
rather than UTF-8 and may contain NUL characters.
"""

from __future__ import annotations

import csv
import io
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)

ENCODINGS = ("utf-8-sig", "cp932")
REQUIRED_COLUMNS = ("genome_id", "role", "donor_category", "acceptor_category", "consensus")
# database role value -> (role, category column)
ROLE_COLUMNS = {
    "electron_donor": ("donor", "donor_category"),
    "electron_acceptor": ("acceptor", "acceptor_category"),
}
GENOME_ID_RE = re.compile(r"GC[AF]_\d+\.\d+")


@dataclass(frozen=True)
class ElectronRecord:
    genome_id: str
    role: str          # "donor" / "acceptor"
    category: str      # donor_category or acceptor_category
    compound: str      # compound_name
    consensus: str     # "used" / "not_used" / "CONFLICTING"
    confidence: str    # confidence_grade (A-D)


def _decode(raw: bytes, path: Path) -> str:
    for encoding in ENCODINGS:
        try:
            text = raw.decode(encoding)
        except UnicodeDecodeError:
            continue
        if encoding != ENCODINGS[0]:
            logger.info("%s is not UTF-8; read as %s", path, encoding)
        if "\x00" in text:
            logger.info("%s: removed %d NUL character(s)", path, text.count("\x00"))
            text = text.replace("\x00", "")
        return text
    raise ValueError(f"{path}: cannot decode as {' or '.join(ENCODINGS)}")


def read_electron_db(path: Path) -> Dict[str, List[ElectronRecord]]:
    """Records per genome. Every genome listed in the database is a key, even without a category."""
    reader = csv.DictReader(io.StringIO(_decode(path.read_bytes(), path), newline=""), delimiter="\t")
    missing = [column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or ())]
    if missing:
        raise ValueError(f"{path}: missing column(s): {', '.join(missing)}")

    records: Dict[str, List[ElectronRecord]] = {}
    invalid_ids = set()
    unknown_roles = 0
    for row in reader:
        genome_id = (row["genome_id"] or "").strip()
        if not GENOME_ID_RE.fullmatch(genome_id):
            invalid_ids.add(genome_id)
            continue
        genome_records = records.setdefault(genome_id, [])
        role_columns = ROLE_COLUMNS.get((row["role"] or "").strip())
        if role_columns is None:
            unknown_roles += 1
            continue
        role, column = role_columns
        category = (row[column] or "").strip()
        if not category:
            continue
        genome_records.append(
            ElectronRecord(
                genome_id=genome_id,
                role=role,
                category=category,
                compound=(row.get("compound_name") or "").strip(),
                consensus=(row["consensus"] or "").strip(),
                confidence=(row.get("confidence_grade") or "").strip(),
            )
        )
    if invalid_ids:
        logger.warning("%s: skipped rows with invalid genome_id %s", path.name, sorted(invalid_ids)[:5])
    if unknown_roles:
        logger.warning("%s: skipped %d row(s) whose role is not electron_donor / electron_acceptor", path.name, unknown_roles)
    return records
