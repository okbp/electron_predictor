"""Command line interface: ``python -m ko_detector {build-config,scan,render-html}``."""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO

from . import __version__
from .html_report import REPORT_FILES, write_reports
from .kofam import DEFAULT_PATTERN, GENOME_DIR_PREFIX, scan_directory
from .models import GenomeResult, KoEntry
from .reference import read_configs, read_references, searchable_kos, write_config
from .results_io import (
    CONFIG_COPY_FILE,
    GENOME_TAXONOMY_FILE,
    read_genome_taxonomy,
    read_results,
    read_run_info,
    write_scan_outputs,
)
from .summary import KoSummary, summarize
from .taxonomy import GenomeTaxonomy, load_genome_taxonomy, taxonomy_available

logger = logging.getLogger("ko_detector")

# scan writes to <DEFAULT_RESULTS_DIR>/<input directory name> when -o is omitted.
DEFAULT_RESULTS_DIR = Path("data") / "results"
# Taxonomy downloaded by scripts/setup_taxonomy.py; used for the tree when present.
DEFAULT_TAXONOMY_DIR = Path("data") / "taxonomy"


def _load_entries(configs: Optional[List[Path]], references: Optional[List[Path]]) -> List[KoEntry]:
    if configs:
        entries = read_configs(configs)
        logger.info("loaded %d config rows from %s", len(entries), ", ".join(map(str, configs)))
    else:
        entries = read_references(references or [])
        logger.info("parsed %d KO entries from reference sheet(s)", len(entries))
    return entries


def _print_summary(results: Sequence[GenomeResult], summaries: Sequence[KoSummary], stream: TextIO) -> None:
    n_valid = sum(result.is_valid for result in results)
    print(
        f"Genomes: {n_valid} analysed, {len(results) - n_valid} excluded "
        f"(of {len(results)} {GENOME_DIR_PREFIX}* directories)",
        file=stream,
    )
    print("ko\tn_genomes\tn_significant\tgene_enzyme", file=stream)
    for summary in summaries:
        print(
            f"{summary.ko}\t{summary.n_genomes}/{summary.total_genomes}\t"
            f"{summary.n_genomes_significant}\t{summary.gene_enzyme}",
            file=stream,
        )


def _report_titles(args: argparse.Namespace) -> Dict[str, str]:
    titles = {"ja": args.title, "en": args.title_en}
    return {lang: title for lang, title in titles.items() if title}


def _taxonomy_dir(args: argparse.Namespace) -> Optional[Path]:
    if args.no_taxonomy:
        return None
    if args.taxonomy_dir is not None:
        if not taxonomy_available(args.taxonomy_dir):
            raise ValueError(f"taxonomy data not found in {args.taxonomy_dir} (run scripts/setup_taxonomy.py)")
        return args.taxonomy_dir
    if taxonomy_available(DEFAULT_TAXONOMY_DIR):
        return DEFAULT_TAXONOMY_DIR
    logger.info("no taxonomy data in %s; reports are written without the tree", DEFAULT_TAXONOMY_DIR)
    return None


def _load_taxonomy(args: argparse.Namespace, results: Sequence[GenomeResult]) -> Optional[Dict[str, GenomeTaxonomy]]:
    directory = _taxonomy_dir(args)
    if directory is None:
        return None
    logger.info("reading taxonomy from %s", directory)
    return load_genome_taxonomy(directory, [result.genome_id for result in results if result.is_valid])


def cmd_build_config(args: argparse.Namespace) -> int:
    entries = read_references(args.reference)
    write_config(entries, args.output, comments=[f"generated from: {', '.join(map(str, args.reference))}"])
    n_kos = len(searchable_kos(entries))
    n_empty = sum(1 for entry in entries if not entry.ko)
    print(f"Wrote {args.output}: {len(entries)} rows, {n_kos} unique KOs, {n_empty} rows without KO")
    return 0


def cmd_scan(args: argparse.Namespace) -> int:
    if not args.input_dir.is_dir():
        logger.error("input directory not found: %s", args.input_dir)
        return 1
    if args.output_dir is None:
        args.output_dir = DEFAULT_RESULTS_DIR / args.input_dir.resolve().name
    entries = _load_entries(args.config, args.reference)
    target_kos = set(searchable_kos(entries))
    if not target_kos:
        logger.error("no KO to search in the configuration")
        return 1

    results = scan_directory(args.input_dir, target_kos, args.pattern, args.significant_only)
    if not results:
        logger.error("no %s* directories found in %s", GENOME_DIR_PREFIX, args.input_dir)
        return 1

    summaries = summarize(entries, results)
    taxonomy = _load_taxonomy(args, results)
    run_info = {
        "tool": f"ko_detector {__version__}",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "input_dir": str(args.input_dir),
        "pattern": args.pattern,
        "significant_only": args.significant_only,
        "config": [str(path) for path in args.config] if args.config else None,
        "reference": [str(path) for path in args.reference] if args.reference else None,
        "taxonomy": taxonomy is not None,
        "n_genome_dirs": len(results),
        "n_genomes_analysed": sum(result.is_valid for result in results),
    }
    paths = write_scan_outputs(args.output_dir, entries, results, summaries, run_info, taxonomy)
    if not args.no_html:
        reports = write_reports(args.output_dir, entries, results, run_info, _report_titles(args), taxonomy)
        paths.update({path.name: path for path in reports.values()})

    _print_summary(results, summaries, sys.stdout)
    for path in paths.values():
        logger.info("wrote %s", path)
    return 0


def cmd_render_html(args: argparse.Namespace) -> int:
    entries = read_configs(args.config or [args.results_dir / CONFIG_COPY_FILE])
    results = read_results(args.results_dir)
    taxonomy_file = args.results_dir / GENOME_TAXONOMY_FILE
    if not args.no_taxonomy and args.taxonomy_dir is None and taxonomy_file.is_file():
        taxonomy = read_genome_taxonomy(taxonomy_file)
    else:
        taxonomy = _load_taxonomy(args, results)
    output_dir = args.output_dir or args.results_dir
    reports = write_reports(
        output_dir, entries, results, read_run_info(args.results_dir), _report_titles(args), taxonomy
    )
    for path in reports.values():
        print(f"Wrote {path}")
    return 0


def _add_title_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--title", help=f"title of the Japanese report ({REPORT_FILES['ja']})")
    parser.add_argument("--title-en", help=f"title of the English report ({REPORT_FILES['en']})")


def _add_taxonomy_arguments(parser: argparse.ArgumentParser, default_help: str) -> None:
    parser.add_argument("--taxonomy-dir", type=Path,
                        help=f"NCBI taxonomy directory for the tree (default: {default_help})")
    parser.add_argument("--no-taxonomy", action="store_true", help="do not add the taxonomy tree")


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-v", "--verbose", action="store_true", help="show progress messages")

    parser = argparse.ArgumentParser(
        prog="ko_detector",
        description="Detect reference KOs in KofamScan results of GCF genomes.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser(
        "build-config", parents=[common], help="convert reference TSV(s) into a KO configuration file"
    )
    p.add_argument("-r", "--reference", type=Path, action="append", required=True,
                   help="reference TSV (repeatable)")
    p.add_argument("-o", "--output", type=Path, required=True, help="configuration file to write")
    p.set_defaults(func=cmd_build_config)

    p = subparsers.add_parser("scan", parents=[common], help="detect KOs in every GCF_* directory")
    p.add_argument("-i", "--input-dir", type=Path, required=True,
                   help=f"directory containing {GENOME_DIR_PREFIX}* genome directories")
    p.add_argument("-o", "--output-dir", type=Path,
                   help=f"directory for result files (default: {DEFAULT_RESULTS_DIR}/<input-dir name>)")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("-c", "--config", type=Path, action="append",
                        help="KO configuration file (from build-config); repeat for donor and acceptor")
    source.add_argument("-r", "--reference", type=Path, action="append",
                        help="reference TSV read directly instead of a configuration file (repeatable)")
    p.add_argument("--pattern", default=DEFAULT_PATTERN, help=f"KofamScan file glob (default: {DEFAULT_PATTERN})")
    p.add_argument("--significant-only", action="store_true",
                   help="count only rows marked '*' (default: every row counts as present)")
    p.add_argument("--no-html", action="store_true", help=f"do not write {' / '.join(REPORT_FILES.values())}")
    _add_taxonomy_arguments(p, f"{DEFAULT_TAXONOMY_DIR} when present")
    _add_title_arguments(p)
    p.set_defaults(func=cmd_scan)

    p = subparsers.add_parser("render-html", parents=[common], help="(re)build the HTML reports from scan results")
    p.add_argument("-i", "--results-dir", type=Path, required=True, help="output directory of a previous scan")
    p.add_argument("-c", "--config", type=Path, action="append",
                   help=f"KO configuration, repeatable (default: <results-dir>/{CONFIG_COPY_FILE})")
    p.add_argument("-o", "--output-dir", type=Path,
                   help=f"directory for {' / '.join(REPORT_FILES.values())} (default: <results-dir>)")
    _add_taxonomy_arguments(p, f"<results-dir>/{GENOME_TAXONOMY_FILE}, else {DEFAULT_TAXONOMY_DIR} when present")
    _add_title_arguments(p)
    p.set_defaults(func=cmd_render_html)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )
    try:
        return args.func(args)
    except (OSError, ValueError) as exc:
        logger.error("%s", exc)
        return 1
