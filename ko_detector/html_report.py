"""Render self-contained HTML reports (genome x KO presence table) in each UI language."""

from __future__ import annotations

import html
import json
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .models import GenomeResult, KoEntry
from .taxonomy import TAXONOMY_RANKS, GenomeTaxonomy

TEMPLATE_DIR = Path(__file__).parent / "templates"
TEMPLATE_PATH = TEMPLATE_DIR / "report.html"
MESSAGES_DIR = TEMPLATE_DIR / "messages"
LANG_PLACEHOLDER = "__REPORT_LANG__"
TITLE_PLACEHOLDER = "__REPORT_TITLE__"
DATA_PLACEHOLDER = "__REPORT_DATA__"
# UI language -> report file name
REPORT_FILES = {"ja": "report.html", "en": "report_en.html"}


def load_messages(lang: str) -> Dict[str, str]:
    return json.loads((MESSAGES_DIR / f"{lang}.json").read_text(encoding="utf-8"))


def _attach_taxonomy(genomes: List[Dict[str, object]], taxonomy: Dict[str, GenomeTaxonomy]) -> Dict[str, object]:
    """Build the taxonomy tree of the listed genomes and point each genome at its species node.

    Nodes are ``[name, rank index, taxid, parent index]``. Every genome gets a node at every rank;
    a missing rank becomes a placeholder node (empty name and taxid) shared under the same parent.
    """
    nodes: List[List[object]] = []
    index: Dict[Tuple[int, str], int] = {}
    for genome in genomes:
        info = taxonomy.get(str(genome["id"]))
        parent = -1
        for rank_index in range(len(TAXONOMY_RANKS)):
            name, taxid = info.lineage[rank_index] if info else ("", "")
            key = (parent, taxid)
            if key not in index:
                index[key] = len(nodes)
                nodes.append([name if taxid else "", rank_index, taxid, parent])
            parent = index[key]
        genome["tax"] = parent
        if info:
            genome["org"] = info.organism_name
            genome["taxid"] = info.taxid
    return {"ranks": list(TAXONOMY_RANKS), "nodes": nodes}


def build_payload(
    entries: Sequence[KoEntry],
    results: Sequence[GenomeResult],
    lang: str = "ja",
    title: Optional[str] = None,
    run_info: Optional[Dict[str, object]] = None,
    taxonomy: Optional[Dict[str, GenomeTaxonomy]] = None,
) -> Dict[str, object]:
    run_info = run_info or {}
    messages = load_messages(lang)
    ko_rows = Counter(entry.ko for entry in entries if entry.ko)
    columns = [
        {
            "ko": entry.ko,
            "role": entry.role,
            "group": entry.functional_group,
            "substrate": entry.substrate,
            "gene": entry.gene_enzyme,
            "priority": entry.priority,
            "tool": entry.recommended_tool,
            "caveats": entry.caveats,
            "source": entry.source,
            "shared": ko_rows[entry.ko] > 1,
        }
        for entry in entries
        if entry.ko
    ]
    genomes: List[Dict[str, object]] = []
    excluded: List[Dict[str, str]] = []
    n_without_hits = 0
    for result in results:
        if not result.is_valid:
            excluded.append({"id": result.genome_id, "status": result.status, "message": result.message})
            continue
        hits = {
            ko: [[hit.gene, hit.score, hit.threshold, hit.evalue, hit.significant] for hit in ko_hits]
            for ko, ko_hits in result.hits_by_ko().items()
        }
        # Genomes carrying none of the target KOs are counted but not listed in the report.
        if not hits:
            n_without_hits += 1
            continue
        genomes.append({"id": result.genome_id, "hits": hits})
    undetectable = [
        {
            "group": entry.functional_group,
            "substrate": entry.substrate,
            "gene": entry.gene_enzyme,
            "koText": entry.ko_text,
            "tool": entry.recommended_tool,
        }
        for entry in entries
        if not entry.ko
    ]
    return {
        "lang": lang,
        "title": title or messages["title"],
        "i18n": messages,
        "generated": run_info.get("generated_at", ""),
        "source": run_info.get("input_dir", ""),
        "significantOnly": bool(run_info.get("significant_only", False)),
        "columns": columns,
        "genomes": genomes,
        "genomesWithoutHits": n_without_hits,
        "excluded": excluded,
        "undetectable": undetectable,
        "taxonomy": _attach_taxonomy(genomes, taxonomy) if taxonomy is not None else None,
    }


def render_html(payload: Dict[str, object]) -> str:
    template = TEMPLATE_PATH.read_text(encoding="utf-8")
    data = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    # "<" only occurs inside JSON strings, so escaping it keeps "</script>" out of the page.
    data = data.replace("<", "\\u003c")
    return (
        template.replace(LANG_PLACEHOLDER, html.escape(str(payload["lang"])))
        .replace(TITLE_PLACEHOLDER, html.escape(str(payload["title"])))
        .replace(DATA_PLACEHOLDER, data)
    )


def write_html(path: Path, payload: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_html(payload), encoding="utf-8")


def write_reports(
    output_dir: Path,
    entries: Sequence[KoEntry],
    results: Sequence[GenomeResult],
    run_info: Optional[Dict[str, object]] = None,
    titles: Optional[Dict[str, str]] = None,
    taxonomy: Optional[Dict[str, GenomeTaxonomy]] = None,
) -> Dict[str, Path]:
    """Write one report per language; ``titles`` overrides the default title per language."""
    titles = titles or {}
    paths: Dict[str, Path] = {}
    for lang, filename in REPORT_FILES.items():
        paths[lang] = output_dir / filename
        write_html(paths[lang], build_payload(entries, results, lang, titles.get(lang), run_info, taxonomy))
    return paths
