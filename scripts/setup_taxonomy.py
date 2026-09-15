#!/usr/bin/env python3
"""Download NCBI new_taxdump and the RefSeq assembly summary into data/taxonomy.

Standard library only. Data that already exists is kept unless --force is given.

    python3 scripts/setup_taxonomy.py                       # set up everything missing
    python3 scripts/setup_taxonomy.py --force               # download everything again
    python3 scripts/setup_taxonomy.py --only taxdump        # one target only
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import tarfile
import time
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Sequence

logger = logging.getLogger("setup_taxonomy")

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = ROOT / "data" / "taxonomy"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
TAXDUMP_DIR = "new_taxdump"
ASSEMBLY_SUMMARY_FILE = "assembly_summary_refseq.txt"
INFO_FILE = "download_info.json"
# new_taxdump counts as set up when these files exist.
TAXDUMP_REQUIRED = ("nodes.dmp", "names.dmp", "rankedlineage.dmp")
TARGETS = ("taxdump", "assembly_summary")
CHUNK_SIZE = 1024 * 1024
RETRIES = 3
TIMEOUT = 60


def fetch_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=TIMEOUT) as response:
        return response.read().decode("utf-8")


def parse_md5(text: str) -> str:
    """Return the checksum from an NCBI ``.md5`` file ("<md5>  <file name>")."""
    token = text.split()[0].lower() if text.split() else ""
    if len(token) != 32 or any(c not in "0123456789abcdef" for c in token):
        raise ValueError(f"unexpected md5 file content: {text!r}")
    return token


def _download_once(url: str, dest: Path, part: Path, label: str) -> Dict[str, object]:
    digest = hashlib.md5()
    size = 0
    with urllib.request.urlopen(url, timeout=TIMEOUT) as response, part.open("wb") as out:
        total = int(response.headers.get("Content-Length") or 0)
        last_modified = response.headers.get("Last-Modified", "")
        next_report = 0.1
        while True:
            chunk = response.read(CHUNK_SIZE)
            if not chunk:
                break
            out.write(chunk)
            digest.update(chunk)
            size += len(chunk)
            if total and size / total >= next_report:
                logger.info("  %s: %3.0f%% (%d / %d MB)", label, 100 * size / total, size >> 20, total >> 20)
                next_report += 0.1
    if total and size != total:
        raise OSError(f"incomplete download of {url}: {size} of {total} bytes")
    part.replace(dest)
    return {"url": url, "last_modified": last_modified, "bytes": size, "md5": digest.hexdigest()}


def download(url: str, dest: Path, retries: int = RETRIES, label: Optional[str] = None) -> Dict[str, object]:
    """Download ``url`` to ``dest`` through a ``.part`` file; return url, Last-Modified, size and md5."""
    part = dest.with_name(dest.name + ".part")
    for attempt in range(1, retries + 1):
        try:
            return _download_once(url, dest, part, label or dest.name)
        except urllib.error.HTTPError:
            part.unlink(missing_ok=True)
            raise
        except OSError as exc:
            part.unlink(missing_ok=True)
            if attempt == retries:
                raise
            wait = 5 * attempt
            logger.warning("download failed (%s); retrying in %d s", exc, wait)
            time.sleep(wait)
    raise AssertionError("unreachable")


def check_member(member: tarfile.TarInfo) -> None:
    path = Path(member.name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe path in archive: {member.name!r}")
    if not (member.isfile() or member.isdir()):
        raise ValueError(f"unsupported entry type in archive: {member.name!r}")


def _replace_dir(new: Path, dest: Path) -> None:
    old = dest.with_name(dest.name + ".old")
    if old.exists():
        shutil.rmtree(old)
    if dest.exists():
        dest.rename(old)
    new.rename(dest)
    if old.exists():
        shutil.rmtree(old)


def extract_taxdump(archive_path: Path, dest: Path) -> None:
    """Extract into a temporary directory first so a failure never leaves a half-written ``dest``."""
    tmp = dest.with_name(dest.name + ".tmp")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    try:
        with tarfile.open(archive_path, "r:gz") as archive:
            for member in archive:
                check_member(member)
                archive.extract(member, tmp)
        missing = [name for name in TAXDUMP_REQUIRED if not (tmp / name).is_file()]
        if missing:
            raise ValueError(f"{archive_path.name} does not contain {', '.join(missing)}")
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    _replace_dir(tmp, dest)


def setup_taxdump(output_dir: Path, url: str, force: bool, keep_archive: bool) -> Optional[Dict[str, object]]:
    dest = output_dir / TAXDUMP_DIR
    if not force and all((dest / name).is_file() for name in TAXDUMP_REQUIRED):
        logger.info("skip taxdump: %s already exists (use --force to download again)", dest)
        return None

    archive = output_dir / url.rsplit("/", 1)[-1]
    expected = parse_md5(fetch_text(url + ".md5"))
    logger.info("downloading %s", url)
    info = download(url, archive)
    try:
        if info["md5"] != expected:
            raise ValueError(f"md5 mismatch for {archive.name}: expected {expected}, got {info['md5']}")
        logger.info("extracting %s into %s", archive.name, dest)
        extract_taxdump(archive, dest)
    finally:
        if not keep_archive:
            archive.unlink(missing_ok=True)
    return info


def setup_assembly_summary(output_dir: Path, url: str, force: bool) -> Optional[Dict[str, object]]:
    dest = output_dir / ASSEMBLY_SUMMARY_FILE
    if not force and dest.is_file():
        logger.info("skip assembly summary: %s already exists (use --force to download again)", dest)
        return None

    tmp = dest.with_name(dest.name + ".download")
    logger.info("downloading %s", url)
    info = download(url, tmp, label=dest.name)
    with tmp.open(encoding="utf-8", errors="replace") as handle:
        first_line = handle.readline()
    if not first_line.startswith("#"):
        tmp.unlink()
        raise ValueError(f"{url} does not look like an assembly summary (first line: {first_line[:80]!r})")
    tmp.replace(dest)
    return info


def update_info(output_dir: Path, updates: Dict[str, Dict[str, object]]) -> None:
    path = output_dir / INFO_FILE
    info = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    downloaded_at = datetime.now().astimezone().isoformat(timespec="seconds")
    for target, values in updates.items():
        info[target] = dict(values, downloaded_at=downloaded_at)
    path.write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Download NCBI new_taxdump and assembly_summary_refseq.txt.")
    parser.add_argument("-o", "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                        help=f"destination directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--only", choices=TARGETS, action="append",
                        help="set up only this target (repeatable; default: all)")
    parser.add_argument("--force", action="store_true", help="download again even if the data already exists")
    parser.add_argument("--keep-archive", action="store_true", help="keep new_taxdump.tar.gz after extraction")
    parser.add_argument("--taxdump-url", default=TAXDUMP_URL, help="new_taxdump.tar.gz URL (md5 read from <url>.md5)")
    parser.add_argument("--assembly-summary-url", default=ASSEMBLY_SUMMARY_URL, help="assembly summary URL")
    parser.add_argument("-q", "--quiet", action="store_true", help="show warnings and errors only")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO, format="%(message)s")
    targets = args.only or list(TARGETS)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    updates: Dict[str, Dict[str, object]] = {}
    try:
        if "taxdump" in targets:
            info = setup_taxdump(args.output_dir, args.taxdump_url, args.force, args.keep_archive)
            if info:
                updates["taxdump"] = info
        if "assembly_summary" in targets:
            info = setup_assembly_summary(args.output_dir, args.assembly_summary_url, args.force)
            if info:
                updates["assembly_summary"] = info
    except (OSError, ValueError, tarfile.TarError) as exc:
        logger.error("error: %s", exc)
        return 1
    finally:
        if updates:
            update_info(args.output_dir, updates)
    logger.info("done: %s", args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
