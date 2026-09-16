"""Resolve genome lineages (domain … species) from the NCBI assembly summary and new_taxdump.

The files are the ones downloaded by ``scripts/setup_taxonomy.py``. Each dump is read line by
line and only the taxids that are needed are kept, so the whole taxonomy is never held in memory.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Set, Tuple

logger = logging.getLogger(__name__)

TAXONOMY_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
ASSEMBLY_SUMMARY_FILE = "assembly_summary_refseq.txt"
TAXDUMP_DIR = "new_taxdump"
TAXDUMP_FILES = ("merged.dmp", "taxidlineage.dmp", "nodes.dmp", "rankedlineage.dmp")
# Dumps released before 2025 call the top rank "superkingdom".
RANK_ALIASES = {"superkingdom": "domain"}


@dataclass(frozen=True)
class GenomeTaxonomy:
    genome_id: str
    taxid: str
    organism_name: str
    # (name, taxid) for each of TAXONOMY_RANKS; ("", "") when the lineage has no such rank.
    lineage: Tuple[Tuple[str, str], ...]

    def rank(self, rank: str) -> Tuple[str, str]:
        return self.lineage[TAXONOMY_RANKS.index(rank)]


def taxonomy_available(taxonomy_dir: Path) -> bool:
    return (taxonomy_dir / ASSEMBLY_SUMMARY_FILE).is_file() and all(
        (taxonomy_dir / TAXDUMP_DIR / name).is_file() for name in TAXDUMP_FILES
    )


def _dmp_fields(line: str) -> List[str]:
    return line.rstrip("\n").rstrip("|").rstrip("\t").split("\t|\t")


def _read_dmp(path: Path, ids: Set[str] = None) -> Iterator[List[str]]:
    """Fields of each line, optionally only for lines whose first column is in ``ids``."""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if ids is not None and line[: line.find("\t")] not in ids:
                continue
            yield _dmp_fields(line)


def read_assembly_summary(path: Path, genome_ids: Iterable[str]) -> Dict[str, Tuple[str, str]]:
    """Map accession -> (taxid, organism_name); another version of the same accession is used as a fallback."""
    wanted = set(genome_ids)
    by_base = {genome_id.split(".")[0]: genome_id for genome_id in wanted}
    exact: Dict[str, Tuple[str, str]] = {}
    other_version: Dict[str, Tuple[str, str]] = {}
    columns = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                if line.startswith("#assembly_accession"):
                    header = line[1:].rstrip("\n").split("\t")
                    columns = (header.index("assembly_accession"), header.index("taxid"), header.index("organism_name"))
                continue
            if columns is None:
                raise ValueError(f"{path}: header line '#assembly_accession ...' not found")
            fields = line.rstrip("\n").split("\t")
            accession, taxid, organism = (fields[i] for i in columns)
            if accession in wanted:
                exact[accession] = (taxid, organism)
            else:
                genome_id = by_base.get(accession.split(".")[0])
                if genome_id is not None:
                    other_version.setdefault(genome_id, (taxid, organism))
    other_version.update(exact)
    return other_version


def load_genome_taxonomy(taxonomy_dir: Path, genome_ids: Iterable[str]) -> Dict[str, GenomeTaxonomy]:
    """Lineage of each genome found in the assembly summary (genomes that are not found are left out)."""
    genome_ids = list(genome_ids)
    dump = taxonomy_dir / TAXDUMP_DIR
    assemblies = read_assembly_summary(taxonomy_dir / ASSEMBLY_SUMMARY_FILE, genome_ids)

    merged = {fields[0]: fields[1] for fields in _read_dmp(dump / "merged.dmp")}
    current = {genome_id: merged.get(taxid, taxid) for genome_id, (taxid, _) in assemblies.items()}
    needed = set(current.values())
    # taxidlineage gives the ancestors, nodes their ranks, rankedlineage their names.
    ancestors = {fields[0]: fields[1].split() for fields in _read_dmp(dump / "taxidlineage.dmp", needed)}
    all_ids = set(needed)
    for chain in ancestors.values():
        all_ids.update(chain)
    ranks = {fields[0]: RANK_ALIASES.get(fields[2], fields[2]) for fields in _read_dmp(dump / "nodes.dmp", all_ids)}
    names = {fields[0]: fields[1] for fields in _read_dmp(dump / "rankedlineage.dmp", all_ids)}

    result: Dict[str, GenomeTaxonomy] = {}
    for genome_id in genome_ids:
        if genome_id not in assemblies:
            continue
        taxid = current[genome_id]
        by_rank: Dict[str, Tuple[str, str]] = {}
        for tid in ancestors.get(taxid, []) + [taxid]:
            rank = ranks.get(tid)
            if rank in TAXONOMY_RANKS:
                by_rank[rank] = (names.get(tid, ""), tid)
        result[genome_id] = GenomeTaxonomy(
            genome_id=genome_id,
            taxid=taxid,
            organism_name=assemblies[genome_id][1],
            lineage=tuple(by_rank.get(rank, ("", "")) for rank in TAXONOMY_RANKS),
        )

    missing = len(set(genome_ids)) - len(result)
    if missing:
        logger.warning("%d genome(s) not found in %s", missing, ASSEMBLY_SUMMARY_FILE)
    unresolved = sum(1 for taxid in needed if taxid not in ranks)
    if unresolved:
        logger.warning("%d taxid(s) not found in %s/nodes.dmp", unresolved, TAXDUMP_DIR)
    return result
