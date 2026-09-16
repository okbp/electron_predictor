"""Parse the KO reference sheets and read/write the KO configuration file."""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .models import KoEntry

logger = logging.getLogger(__name__)

CONFIG_COLUMNS = (
    "ko",
    "role",
    "functional_group",
    "substrate",
    "gene_enzyme",
    "priority",
    "detectable_via_kegg",
    "recommended_tool",
    "caveats",
    "ko_text",
    "source",
)
REQUIRED_CONFIG_COLUMNS = ("ko", "functional_group", "gene_enzyme")

# Accepted header names in the reference sheet, per config field.
REFERENCE_COLUMNS: Dict[str, Sequence[str]] = {
    "functional_group": ("Functional group",),
    "substrate": ("Electron donor", "Electron acceptor", "Substrate"),
    "gene_enzyme": ("Gene / enzyme", "Gene/enzyme", "Gene"),
    "ko": ("KO",),
    "priority": ("Priority",),
    "detectable_via_kegg": ("Detectable via KEGG",),
    "caveats": ("Degeneracy and caveats", "Caveats"),
    "recommended_tool": ("Recommended tool",),
}
REQUIRED_REFERENCE_FIELDS = ("functional_group", "gene_enzyme", "ko")
# The substrate column's header tells which sheet this is.
ROLE_BY_SUBSTRATE_HEADER = {"Electron donor": "donor", "Electron acceptor": "acceptor"}
ROLES = ("donor", "acceptor")

KO_ID_RE = re.compile(r"K\d{5}")
# A single KO ("K17218") or a range ("K17222-K17227", en dash and "K17222-17227" accepted).
KO_TOKEN_RE = re.compile(r"\bK(\d{5})\b(?:\s*[-‐-―~]\s*K?(\d{5})\b)?")
# What may appear between KOs in a KO list ("K00368 / K15864", "K03385, K15876 / K00362, K00363").
KO_LIST_FILLER_RE = re.compile(r"[\s,;/]+|\b(?:and|or)\b")
MAX_KO_RANGE = 100


def parse_ko_field(text: str) -> List[str]:
    """Extract KO identifiers from a KO cell, expanding ranges."""
    kos: List[str] = []
    for match in KO_TOKEN_RE.finditer(text):
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else start
        if end < start:
            raise ValueError(f"invalid KO range {match.group(0)!r}")
        if end - start >= MAX_KO_RANGE:
            raise ValueError(f"KO range too large {match.group(0)!r}")
        for number in range(start, end + 1):
            ko = f"K{number:05d}"
            if ko not in kos:
                kos.append(ko)
    return kos


def is_ko_list(text: str) -> bool:
    """True when the cell holds only KOs, ranges and separators.

    A description that merely mentions a KO — e.g. "Closely related to psrABC (K08352); poor KO
    resolution" on the phsABC row — is not a KO list for that row.
    """
    return bool(text.strip()) and not KO_LIST_FILLER_RE.sub("", KO_TOKEN_RE.sub("", text))


def _resolve_reference_columns(fieldnames: Optional[Sequence[str]], path: Path) -> Tuple[Dict[str, str], str]:
    present = {name.strip(): name for name in fieldnames or ()}
    columns: Dict[str, str] = {}
    role = ""
    for key, candidates in REFERENCE_COLUMNS.items():
        for candidate in candidates:
            if candidate in present:
                columns[key] = present[candidate]
                if key == "substrate":
                    role = ROLE_BY_SUBSTRATE_HEADER.get(candidate, "")
                break
    missing = [key for key in REQUIRED_REFERENCE_FIELDS if key not in columns]
    if missing:
        raise ValueError(f"{path}: missing reference column(s) for {', '.join(missing)}")
    return columns, role


def read_reference(path: Path) -> List[KoEntry]:
    """Read a reference sheet; each KO of each data row becomes one entry."""
    entries: List[KoEntry] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns, role = _resolve_reference_columns(reader.fieldnames, path)
        for row in reader:
            values = {key: "" for key in REFERENCE_COLUMNS}
            values.update({key: (row.get(column) or "").strip() for key, column in columns.items()})
            # Blank separator rows and footnotes only fill the first column.
            if not values["functional_group"] or not (values["gene_enzyme"] or values["ko"]):
                continue
            source = f"{path.name}:{reader.line_num}"
            try:
                kos = parse_ko_field(values["ko"])
            except ValueError as exc:
                raise ValueError(f"{source}: {exc}") from None
            if kos and not is_ko_list(values["ko"]):
                logger.warning(
                    "%s: KO field %r is a description, not a KO list (mentions %s); kept with empty ko for %s",
                    source, values["ko"], ", ".join(kos), values["gene_enzyme"],
                )
                kos = []
            elif not kos:
                logger.warning(
                    "%s: no KO in %r for %s (kept with empty ko)", source, values["ko"], values["gene_enzyme"]
                )
            for ko in kos or [""]:
                entries.append(
                    KoEntry(
                        ko=ko,
                        functional_group=values["functional_group"],
                        substrate=values["substrate"],
                        gene_enzyme=values["gene_enzyme"],
                        priority=values["priority"],
                        detectable_via_kegg=values["detectable_via_kegg"],
                        recommended_tool=values["recommended_tool"],
                        caveats=values["caveats"],
                        ko_text=values["ko"],
                        source=source,
                        role=role,
                    )
                )
    return entries


def read_references(paths: Iterable[Path]) -> List[KoEntry]:
    entries: List[KoEntry] = []
    for path in paths:
        entries.extend(read_reference(path))
    return entries


def write_config(entries: Sequence[KoEntry], path: Path, comments: Sequence[str] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.write("# KO configuration for ko_detector\n")
        for comment in comments:
            handle.write(f"# {comment}\n")
        handle.write("# Rows with an empty 'ko' are not searched; enter a KO to enable them.\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(CONFIG_COLUMNS)
        for entry in entries:
            writer.writerow([getattr(entry, column) for column in CONFIG_COLUMNS])


def read_config(path: Path) -> List[KoEntry]:
    entries: List[KoEntry] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        missing = [column for column in REQUIRED_CONFIG_COLUMNS if column not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path}: missing config column(s): {', '.join(missing)}")
        for row in reader:
            values = {column: (row.get(column) or "").strip() for column in CONFIG_COLUMNS}
            values["ko"] = values["ko"].upper()
            if values["ko"] and not KO_ID_RE.fullmatch(values["ko"]):
                raise ValueError(f"{path}: invalid KO {values['ko']!r} ({values['gene_enzyme']})")
            if values["role"] and values["role"] not in ROLES:
                raise ValueError(f"{path}: invalid role {values['role']!r} ({values['gene_enzyme']})")
            if not any(values.values()):
                continue
            entries.append(KoEntry(**values))
    return entries


def read_configs(paths: Iterable[Path]) -> List[KoEntry]:
    """Entries of several configuration files in order; an entry repeated verbatim is kept once."""
    entries: List[KoEntry] = []
    for path in paths:
        entries.extend(read_config(path))
    return list(dict.fromkeys(entries))


def searchable_kos(entries: Iterable[KoEntry]) -> List[str]:
    """Unique non-empty KOs in configuration order."""
    return list(dict.fromkeys(entry.ko for entry in entries if entry.ko))


def entries_by_ko(entries: Iterable[KoEntry]) -> Dict[str, List[KoEntry]]:
    grouped: Dict[str, List[KoEntry]] = {}
    for entry in entries:
        if entry.ko:
            grouped.setdefault(entry.ko, []).append(entry)
    return grouped


def join_unique(values: Iterable[str], separator: str = "; ") -> str:
    return separator.join(dict.fromkeys(value for value in values if value))
