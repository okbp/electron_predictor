"""Aggregate per-genome results into per-KO counts."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import List, Sequence

from .models import GenomeResult, KoEntry
from .reference import entries_by_ko, join_unique, searchable_kos


@dataclass(frozen=True)
class KoSummary:
    ko: str
    entries: Sequence[KoEntry]
    n_genomes: int
    n_genomes_significant: int
    total_genomes: int

    @property
    def fraction(self) -> float:
        return self.n_genomes / self.total_genomes if self.total_genomes else 0.0

    @property
    def gene_enzyme(self) -> str:
        return join_unique(entry.gene_enzyme for entry in self.entries)

    @property
    def functional_group(self) -> str:
        return join_unique(entry.functional_group for entry in self.entries)

    @property
    def substrate(self) -> str:
        return join_unique(entry.substrate for entry in self.entries)


def summarize(entries: Sequence[KoEntry], results: Sequence[GenomeResult]) -> List[KoSummary]:
    """Count genomes carrying each KO; KOs found in no genome are included with 0."""
    valid = [result for result in results if result.is_valid]
    present: Counter = Counter()
    present_significant: Counter = Counter()
    for result in valid:
        for ko, hits in result.hits_by_ko().items():
            present[ko] += 1
            if any(hit.significant for hit in hits):
                present_significant[ko] += 1

    grouped = entries_by_ko(entries)
    return [
        KoSummary(
            ko=ko,
            entries=grouped[ko],
            n_genomes=present[ko],
            n_genomes_significant=present_significant[ko],
            total_genomes=len(valid),
        )
        for ko in searchable_kos(entries)
    ]
