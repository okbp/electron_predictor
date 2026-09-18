"""Phenotypes and electron donor / acceptor categories per genome from ``phenotype_data.tsv``.

The file has one row per organism, joined to genomes by ``genome_accession``. Each trait column
(``oxygen``, ``optimal_temperature``, ``optimal_pH``, ``habitat``) has companion columns
(``<trait>_status``, ``_source_tier``, ``_source``, ``_evidence``, ``_url``, ``_other_values``).
``electron_donors`` / ``electron_acceptors`` list tokens such as ``S2O3(2-) thiosulfate`` separated
by ``; `` and are turned into the donor / acceptor categories shown in the report through three
reviewed tables in ``config/``:

- ``electron_compound_categories.tsv``: (role, token) -> category;
- ``phenotype_corrections.tsv``: tokens listed on the wrong side (moved to the other role);
- ``habitat_classes.tsv``: habitat value (HTML tags removed) -> environment class.
"""

from __future__ import annotations

import collections
import csv
import dataclasses
import io
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

ENCODINGS = ("utf-8-sig", "cp932")
TRAITS = ("oxygen", "optimal_temperature", "optimal_pH", "habitat")
ROLES = ("donor", "acceptor")
OXYGEN_LEVELS = ("obligate aerobic", "aerobic", "microaerophilic", "facultative anaerobic", "anaerobic", "obligate anaerobic")
OXYGEN_ALIASES = {"facultative aerobic": "facultative anaerobic"}   # user decision
REQUIRED_COLUMNS = ("organism_id", "organism", "genome_accession", "electron_donors", "electron_acceptors") + TRAITS
COMPOUND_TABLE = "electron_compound_categories.tsv"
CORRECTIONS_TABLE = "phenotype_corrections.tsv"
HABITAT_TABLE = "habitat_classes.tsv"
GENOME_ID_RE = re.compile(r"GC[AF]_\d+\.\d+")
NUMBER_RE = re.compile(r"^\s*(-?\d+(?:\.\d+)?)\s*(?:-\s*(-?\d+(?:\.\d+)?))?\s*$")
BRACKETED_RE = re.compile(r"([^\[\];]+?)\s*\[([^\]]*)\]")
TAG_RE = re.compile(r"<[^>]+>")
EVIDENCE_MAX = 400


@dataclass(frozen=True)
class TraitValue:
    value: str                       # as written (oxygen lower-cased, habitat without HTML tags)
    klass: str = ""                  # oxygen level or habitat class; empty for numbers
    number: Optional[float] = None   # temperature / pH; midpoint of a range
    status: str = ""
    tier: str = ""
    source: str = ""
    evidence: str = ""
    url: str = ""
    other_values: str = ""

    @property
    def tier_rank(self) -> int:
        match = re.match(r"\d", self.tier)
        return int(match.group()) if match else 9


@dataclass
class GenomePhenotype:
    genome_id: str
    organisms: List[str]
    traits: Dict[str, TraitValue]


@dataclass(frozen=True)
class ElectronRecord:
    genome_id: str
    role: str        # "donor" / "acceptor"
    category: str
    compound: str    # token as written, e.g. "S2O3(2-) thiosulfate"
    consensus: str   # "used" / "CONFLICTING" / "not_used" / "audit_error"
    source: str      # e.g. "literature pipeline (33 papers)"


@dataclass(frozen=True)
class PhenotypeTables:
    token_categories: Dict[Tuple[str, str], str]
    corrections: Dict[Tuple[str, str, str], str]
    habitat_classes: Dict[str, str]


def _decode(raw: bytes, path: Path) -> str:
    for encoding in ENCODINGS:
        try:
            text = raw.decode(encoding)
        except UnicodeDecodeError:
            continue
        if encoding != ENCODINGS[0]:
            logger.info("%s is not UTF-8; read as %s", path, encoding)
        return text.replace("\x00", "")
    raise ValueError(f"{path}: cannot decode as {' or '.join(ENCODINGS)}")


def _read_table(path: Path, required: Tuple[str, ...]) -> List[Dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        missing = [column for column in required if column not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path}: missing column(s): {', '.join(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def read_tables(config_dir: Path) -> PhenotypeTables:
    token_categories: Dict[Tuple[str, str], str] = {}
    for row in _read_table(config_dir / COMPOUND_TABLE, ("role", "token", "category")):
        if row["role"] not in ROLES:
            raise ValueError(f"{COMPOUND_TABLE}: invalid role {row['role']!r}")
        token_categories[(row["role"], row["token"])] = row["category"] or "Other"
    corrections: Dict[Tuple[str, str, str], str] = {}
    for row in _read_table(config_dir / CORRECTIONS_TABLE, ("genome_accession", "token", "from_role", "to_role")):
        if row["from_role"] not in ROLES or row["to_role"] not in ROLES:
            raise ValueError(f"{CORRECTIONS_TABLE}: invalid role in {row}")
        corrections[(row["genome_accession"], row["token"], row["from_role"])] = row["to_role"]
    habitat_classes = {
        row["habitat"]: row["class"] for row in _read_table(config_dir / HABITAT_TABLE, ("habitat", "class")) if row["class"]
    }
    return PhenotypeTables(token_categories, corrections, habitat_classes)


def is_shifted(row: Dict[str, str]) -> bool:
    """A row whose cells slid by a few columns (a source tier where a status belongs)."""
    return bool(re.match(r"^\d \(", row.get("oxygen_status") or "")) or (row.get("optimal_pH_kind") or "") not in ("", "optimum")


def parse_number(text: str) -> Optional[Tuple[float, float]]:
    """``25`` -> (25, 25); ``55-60`` -> (55, 60); ``-2-30`` -> (-2, 30)."""
    match = NUMBER_RE.match(text or "")
    if not match or not (text or "").strip():
        return None
    low = float(match.group(1))
    return low, float(match.group(2)) if match.group(2) else low


def clean_habitat(text: str) -> str:
    return re.sub(r"\s+", " ", TAG_RE.sub("", text or "")).strip()


def split_tokens(text: str) -> List[str]:
    return [token.strip() for token in (text or "").split(";") if token.strip()]


def bracketed_tokens(text: str) -> Dict[str, str]:
    """``S0 elemental sulfur [1 paper]; Fe(II) [1 paper; putative only]`` -> token -> bracket text."""
    return {token.strip(): note.strip() for token, note in BRACKETED_RE.findall(text or "")}


def token_sources(tokens: List[str], text: str) -> Dict[str, str]:
    """Split ``"tok1: source; tok2: source (a; b)"`` by the known tokens (sources may contain ';')."""
    found = []
    for token in tokens:
        match = re.search(r"(?:^|;\s*)" + re.escape(token) + r":", text or "")
        if match:
            found.append((match.end(), match.start(), token))
    found.sort()
    sources: Dict[str, str] = {}
    for i, (end, _, token) in enumerate(found):
        stop = found[i + 1][1] if i + 1 < len(found) else len(text)
        sources[token] = text[end:stop].strip().rstrip(";").strip()
    return sources


def is_audit_error(source: str) -> bool:
    """Only a failed audit and no other support (a literature or audited-OK source keeps the token)."""
    return "audited ERROR" in source and "literature" not in source and "audited OK" not in source


def _trait(row: Dict[str, str], trait: str, tables: PhenotypeTables, unknown_habitats: collections.Counter) -> Optional[TraitValue]:
    raw = (row.get(trait) or "").strip()
    if not raw:
        return None
    meta = {
        "status": (row.get(f"{trait}_status") or "").strip(),
        "tier": (row.get(f"{trait}_source_tier") or "").strip(),
        "source": (row.get(f"{trait}_source") or "").strip(),
        "evidence": (row.get(f"{trait}_evidence") or "").strip()[:EVIDENCE_MAX],
        "url": (row.get(f"{trait}_url") or "").strip(),
        "other_values": (row.get(f"{trait}_other_values") or "").strip(),
    }
    if trait == "oxygen":
        level = OXYGEN_ALIASES.get(raw.lower(), raw.lower())
        if level not in OXYGEN_LEVELS:
            logger.warning("%s: unknown oxygen value %r", row.get("organism_id"), raw)
            return None
        return TraitValue(value=raw, klass=level, **meta)
    if trait == "habitat":
        value = clean_habitat(raw)
        klass = tables.habitat_classes.get(value)
        if klass is None:
            unknown_habitats[value] += 1
            klass = "other"
        return TraitValue(value=value, klass=klass, **meta)
    pair = parse_number(raw)
    if pair is None:
        logger.warning("%s: %s %r is not a number or range", row.get("organism_id"), trait, raw)
        return None
    return TraitValue(value=raw, number=round((pair[0] + pair[1]) / 2, 6), **meta)   # 6.8, not 6.800000000000001


def _electron(row: Dict[str, str], genome_id: str, tables: PhenotypeTables, unknown_tokens: collections.Counter) -> List[ElectronRecord]:
    records: List[ElectronRecord] = []
    conflicting = set(split_tokens(row.get("electron_donor_acceptor_conflicting_reports", "")))

    def add(role: str, token: str, consensus: str, source: str) -> None:
        role = tables.corrections.get((genome_id, token, role), role)
        category = tables.token_categories.get((role, token))
        if category is None:
            unknown_tokens[(role, token)] += 1
            category = "Other"
        records.append(ElectronRecord(genome_id, role, category, token, consensus, source))

    for role in ROLES:
        listed = split_tokens(row.get(f"electron_{role}s", ""))
        sources = token_sources(listed, row.get(f"electron_{role}s_source", ""))
        not_used = bracketed_tokens(row.get(f"reported_not_used_as_{role}", ""))
        for token in listed:
            source = sources.get(token, "")
            if is_audit_error(source):
                consensus = "audit_error"
            elif token in conflicting or token in not_used:
                consensus = "CONFLICTING"
            else:
                consensus = "used"
            add(role, token, consensus, source)
        for token, note in not_used.items():
            if token not in listed:
                add(role, token, "not_used", f"reported not used ({note})")
    return records


def _merge_trait(kept: Tuple[TraitValue, str], new: Tuple[TraitValue, str]) -> Tuple[TraitValue, str]:
    """Two organism rows share a genome: keep the better-sourced (value, organism), note the other value."""
    if kept[0].value == new[0].value:
        return kept
    (best, organism), (other, other_organism) = (kept, new) if kept[0].tier_rank <= new[0].tier_rank else (new, kept)
    note = f"{other.value} [{other_organism}]"
    return dataclasses.replace(best, other_values="; ".join(v for v in (best.other_values, note) if v)), organism


def read_phenotype(path: Path, tables: PhenotypeTables) -> Tuple[Dict[str, GenomePhenotype], Dict[str, List[ElectronRecord]]]:
    """Phenotypes and electron records per genome accession (both dicts have the same keys)."""
    reader = csv.DictReader(io.StringIO(_decode(path.read_bytes(), path), newline=""), delimiter="\t")
    missing = [column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or ())]
    if missing:
        raise ValueError(f"{path}: missing column(s): {', '.join(missing)}")

    phenotypes: Dict[str, GenomePhenotype] = {}
    electron: Dict[str, List[ElectronRecord]] = {}
    shifted: List[str] = []
    without_accession = 0
    invalid: List[str] = []
    unknown_habitats: collections.Counter = collections.Counter()
    unknown_tokens: collections.Counter = collections.Counter()
    owners: Dict[Tuple[str, str], str] = {}     # (genome, trait) -> organism the kept value comes from
    for row in reader:
        if is_shifted(row):
            shifted.append(f"{row.get('organism_id')} {row.get('organism')}")
            continue
        genome_id = (row.get("genome_accession") or "").strip()
        if not genome_id:
            without_accession += 1
            continue
        if not GENOME_ID_RE.fullmatch(genome_id):
            invalid.append(genome_id)
            continue
        organism = (row.get("organism") or "").strip()
        traits = {}
        for trait in TRAITS:
            value = _trait(row, trait, tables, unknown_habitats)
            if value is not None:
                traits[trait] = value
        existing = phenotypes.get(genome_id)
        if existing is None:
            phenotypes[genome_id] = GenomePhenotype(genome_id, [organism], traits)
        else:
            logger.info("%s: several organism rows (%s, %s); taking their union", genome_id, "; ".join(existing.organisms), organism)
            existing.organisms.append(organism)
            for trait, value in traits.items():
                if trait in existing.traits:
                    kept = (existing.traits[trait], owners[(genome_id, trait)])
                    existing.traits[trait], owners[(genome_id, trait)] = _merge_trait(kept, (value, organism))
                else:
                    existing.traits[trait] = value
        for trait in phenotypes[genome_id].traits:
            owners.setdefault((genome_id, trait), organism)
        records = electron.setdefault(genome_id, [])
        seen = {(r.role, r.compound) for r in records}
        for record in _electron(row, genome_id, tables, unknown_tokens):
            if (record.role, record.compound) not in seen:     # a union; the first mention wins
                seen.add((record.role, record.compound))
                records.append(record)

    if shifted:
        logger.warning("%s: skipped column-shifted row(s): %s", path.name, "; ".join(shifted))
    if invalid:
        logger.warning("%s: skipped rows with invalid genome_accession %s", path.name, invalid[:5])
    if unknown_habitats:
        logger.warning("%d habitat value(s) not in config/%s (shown as other): %s", len(unknown_habitats), HABITAT_TABLE,
                       "; ".join(list(unknown_habitats)[:5]))
    if unknown_tokens:
        logger.warning("%d electron token(s) not in config/%s (shown as Other): %s", len(unknown_tokens), COMPOUND_TABLE,
                       "; ".join(f"{role}:{token}" for role, token in list(unknown_tokens)[:5]))
    logger.info("%s: %d genome(s) linked, %d row(s) without genome_accession", path.name, len(phenotypes), without_accession)
    return phenotypes, electron
