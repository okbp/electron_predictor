"""Data structures shared across ko_detector modules."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

STATUS_OK = "ok"
STATUS_MULTIPLE_FILES = "multiple_files"
STATUS_NO_FILE = "no_file"
STATUS_ERROR = "error"
# Genomes with these statuses count towards the summary denominator.
VALID_STATUSES = (STATUS_OK, STATUS_MULTIPLE_FILES)


@dataclass(frozen=True)
class KoEntry:
    """One KO taken from one row of the reference sheet.

    ``ko`` is empty when the row has no KEGG KO (e.g. "Poor KEGG coverage");
    such entries are kept so they can be filled in by hand but are not searched.
    """

    ko: str
    functional_group: str
    substrate: str
    gene_enzyme: str
    priority: str = ""
    detectable_via_kegg: str = ""
    recommended_tool: str = ""
    caveats: str = ""
    ko_text: str = ""
    source: str = ""


@dataclass(frozen=True)
class Hit:
    """One gene assigned to a target KO by KofamScan."""

    genome_id: str
    ko: str
    gene: str
    threshold: str
    score: str
    evalue: str
    significant: bool


@dataclass
class GenomeResult:
    genome_id: str
    status: str
    files: List[Path] = field(default_factory=list)
    n_rows: int = 0
    hits: List[Hit] = field(default_factory=list)
    message: str = ""

    @property
    def is_valid(self) -> bool:
        return self.status in VALID_STATUSES

    def hits_by_ko(self) -> Dict[str, List[Hit]]:
        grouped: Dict[str, List[Hit]] = {}
        for hit in self.hits:
            grouped.setdefault(hit.ko, []).append(hit)
        return grouped
