"""Discover GCF genome directories and read KofamScan ``*.kolist_gene.tsv`` files."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import AbstractSet, List, Tuple

from .models import (
    STATUS_ERROR,
    STATUS_MULTIPLE_FILES,
    STATUS_NO_FILE,
    STATUS_OK,
    GenomeResult,
    Hit,
)

logger = logging.getLogger(__name__)

GENOME_DIR_PREFIX = "GCF_"
DEFAULT_PATTERN = "*.kolist_gene.tsv"
# Column order used when a file has no header line.
KOFAM_COLUMNS = ("KO", "gene name", "thrshld", "score", "E-value", "significant")


def find_genome_dirs(input_dir: Path, prefix: str = GENOME_DIR_PREFIX) -> List[Path]:
    """Genome directories directly under ``input_dir`` (not recursive)."""
    return sorted(path for path in input_dir.iterdir() if path.is_dir() and path.name.startswith(prefix))


def read_kolist(
    path: Path,
    target_kos: AbstractSet[str],
    genome_id: str,
    significant_only: bool = False,
) -> Tuple[int, List[Hit]]:
    """Return the number of data rows and the hits for ``target_kos``.

    Every row is treated as a gene the genome carries; ``significant_only``
    keeps only rows marked ``*`` (score above the KofamScan threshold).
    """
    index = {name: i for i, name in enumerate(KOFAM_COLUMNS)}
    n_rows = 0
    hits: List[Hit] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0].strip() == "KO":
                index = {name.strip(): i for i, name in enumerate(fields)}
                if "gene name" not in index:
                    raise ValueError(f"{path}: header has no 'gene name' column")
                continue
            n_rows += 1

            def get(name: str) -> str:
                i = index.get(name)
                return fields[i].strip() if i is not None and i < len(fields) else ""

            ko = get("KO")
            if ko not in target_kos:
                continue
            significant = get("significant") == "*"
            if significant_only and not significant:
                continue
            hits.append(
                Hit(
                    genome_id=genome_id,
                    ko=ko,
                    gene=get("gene name"),
                    threshold=get("thrshld"),
                    score=get("score"),
                    evalue=get("E-value"),
                    significant=significant,
                )
            )
    return n_rows, hits


def scan_genome(
    genome_dir: Path,
    target_kos: AbstractSet[str],
    pattern: str = DEFAULT_PATTERN,
    significant_only: bool = False,
) -> GenomeResult:
    genome_id = genome_dir.name
    files = sorted(path for path in genome_dir.glob(pattern) if path.is_file())
    if not files:
        return GenomeResult(genome_id, STATUS_NO_FILE, message=f"no file matching {pattern}")

    status = STATUS_OK if len(files) == 1 else STATUS_MULTIPLE_FILES
    result = GenomeResult(genome_id, status, files=files)
    seen = set()
    try:
        for path in files:
            n_rows, hits = read_kolist(path, target_kos, genome_id, significant_only)
            result.n_rows += n_rows
            for hit in hits:
                key = (hit.ko, hit.gene)
                if key not in seen:
                    seen.add(key)
                    result.hits.append(hit)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        return GenomeResult(genome_id, STATUS_ERROR, files=files, message=str(exc))
    if status == STATUS_MULTIPLE_FILES:
        result.message = f"{len(files)} files merged"
    return result


def scan_directory(
    input_dir: Path,
    target_kos: AbstractSet[str],
    pattern: str = DEFAULT_PATTERN,
    significant_only: bool = False,
) -> List[GenomeResult]:
    genome_dirs = find_genome_dirs(input_dir)
    logger.info("found %d %s* directories in %s", len(genome_dirs), GENOME_DIR_PREFIX, input_dir)
    results: List[GenomeResult] = []
    for i, genome_dir in enumerate(genome_dirs, 1):
        result = scan_genome(genome_dir, target_kos, pattern, significant_only)
        if result.status != STATUS_OK:
            logger.warning("%s: %s (%s)", result.genome_id, result.status, result.message)
        results.append(result)
        if i % 500 == 0:
            logger.info("processed %d / %d genomes", i, len(genome_dirs))
    return results
