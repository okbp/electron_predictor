"""Write scan results as TSV/JSON and read them back for report rendering."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence

from .models import GenomeResult, Hit, KoEntry
from .reference import entries_by_ko, join_unique, searchable_kos, write_config
from .summary import KoSummary

STATUS_FILE = "genome_status.tsv"
HITS_FILE = "genome_ko_hits.tsv"
MATRIX_FILE = "genome_ko_matrix.tsv"
SUMMARY_FILE = "ko_summary.tsv"
CONFIG_COPY_FILE = "ko_config_used.tsv"
RUN_INFO_FILE = "run_info.json"

STATUS_COLUMNS = ("genome_id", "status", "n_rows", "n_target_kos", "n_hits", "files", "message")
HIT_COLUMNS = (
    "genome_id",
    "ko",
    "gene",
    "threshold",
    "score",
    "evalue",
    "significant",
    "gene_enzyme",
    "functional_group",
)
SUMMARY_COLUMNS = (
    "ko",
    "gene_enzyme",
    "functional_group",
    "substrate",
    "n_genomes",
    "n_genomes_significant",
    "fraction",
    "total_genomes",
)
FILE_SEPARATOR = ";"


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _read_tsv(path: Path) -> Iterator[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_status(path: Path, results: Sequence[GenomeResult]) -> None:
    rows = (
        (
            result.genome_id,
            result.status,
            result.n_rows,
            len(result.hits_by_ko()),
            len(result.hits),
            FILE_SEPARATOR.join(str(file) for file in result.files),
            result.message,
        )
        for result in results
    )
    _write_tsv(path, STATUS_COLUMNS, rows)


def write_hits(path: Path, entries: Sequence[KoEntry], results: Sequence[GenomeResult]) -> None:
    grouped = entries_by_ko(entries)
    rows = (
        (
            hit.genome_id,
            hit.ko,
            hit.gene,
            hit.threshold,
            hit.score,
            hit.evalue,
            "*" if hit.significant else "",
            join_unique(entry.gene_enzyme for entry in grouped.get(hit.ko, ())),
            join_unique(entry.functional_group for entry in grouped.get(hit.ko, ())),
        )
        for result in results
        if result.is_valid
        for hit in result.hits
    )
    _write_tsv(path, HIT_COLUMNS, rows)


def write_matrix(path: Path, entries: Sequence[KoEntry], results: Sequence[GenomeResult]) -> None:
    """Genome x KO table; values are the number of genes assigned to the KO (0 = absent)."""
    kos = searchable_kos(entries)
    rows = []
    for result in results:
        if not result.is_valid:
            continue
        by_ko = result.hits_by_ko()
        rows.append([result.genome_id] + [len(by_ko.get(ko, ())) for ko in kos])
    _write_tsv(path, ["genome_id"] + kos, rows)


def write_summary(path: Path, summaries: Sequence[KoSummary]) -> None:
    rows = (
        (
            summary.ko,
            summary.gene_enzyme,
            summary.functional_group,
            summary.substrate,
            summary.n_genomes,
            summary.n_genomes_significant,
            f"{summary.fraction:.4f}",
            summary.total_genomes,
        )
        for summary in summaries
    )
    _write_tsv(path, SUMMARY_COLUMNS, rows)


def write_scan_outputs(
    output_dir: Path,
    entries: Sequence[KoEntry],
    results: Sequence[GenomeResult],
    summaries: Sequence[KoSummary],
    run_info: Dict[str, object],
) -> Dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {name: output_dir / name for name in (STATUS_FILE, HITS_FILE, MATRIX_FILE, SUMMARY_FILE)}
    write_status(paths[STATUS_FILE], results)
    write_hits(paths[HITS_FILE], entries, results)
    write_matrix(paths[MATRIX_FILE], entries, results)
    write_summary(paths[SUMMARY_FILE], summaries)

    paths[CONFIG_COPY_FILE] = output_dir / CONFIG_COPY_FILE
    write_config(entries, paths[CONFIG_COPY_FILE], comments=["copy of the configuration used for this scan"])
    paths[RUN_INFO_FILE] = output_dir / RUN_INFO_FILE
    paths[RUN_INFO_FILE].write_text(json.dumps(run_info, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return paths


def read_run_info(output_dir: Path) -> Dict[str, object]:
    path = output_dir / RUN_INFO_FILE
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_results(output_dir: Path) -> List[GenomeResult]:
    """Rebuild genome results from ``genome_status.tsv`` and ``genome_ko_hits.tsv``."""
    results: Dict[str, GenomeResult] = {}
    for row in _read_tsv(output_dir / STATUS_FILE):
        results[row["genome_id"]] = GenomeResult(
            genome_id=row["genome_id"],
            status=row["status"],
            files=[Path(file) for file in row["files"].split(FILE_SEPARATOR) if file],
            n_rows=int(row["n_rows"] or 0),
            message=row["message"],
        )
    for row in _read_tsv(output_dir / HITS_FILE):
        result = results.get(row["genome_id"])
        if result is None:
            raise ValueError(f"{HITS_FILE}: genome {row['genome_id']!r} is not in {STATUS_FILE}")
        result.hits.append(
            Hit(
                genome_id=row["genome_id"],
                ko=row["ko"],
                gene=row["gene"],
                threshold=row["threshold"],
                score=row["score"],
                evalue=row["evalue"],
                significant=row["significant"] == "*",
            )
        )
    return list(results.values())
