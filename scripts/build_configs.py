#!/usr/bin/env python3
"""Build the KO configuration files for the electron donor and acceptor reference sheets.

Standard library only; run from anywhere. Wraps ``python3 -m ko_detector build-config`` with the
fixed file names used in ``config/``.

    python3 scripts/build_configs.py                   # both sheets
    python3 scripts/build_configs.py --only acceptor   # one sheet
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ko_detector.reference import read_reference, searchable_kos, write_config  # noqa: E402

DEFAULT_CONFIG_DIR = ROOT / "config"
# role -> (reference sheet, configuration file), both inside the config directory
SHEETS = {
    "donor": ("chemolithoautotroph_donor_acceptor_donor_ko_reference.tsv", "donor_ko_config.tsv"),
    "acceptor": ("chemolithoautotroph_donor_acceptor_acceptor_ko_reference.tsv", "acceptor_ko_config.tsv"),
}

logger = logging.getLogger("build_configs")


def build(role: str, config_dir: Path) -> Tuple[Path, int, int, List[str]]:
    reference_name, config_name = SHEETS[role]
    reference = config_dir / reference_name
    output = config_dir / config_name
    entries = read_reference(reference)
    wrong = sorted({entry.role for entry in entries} - {role})
    if wrong:
        raise ValueError(f"{reference}: expected an 'Electron {role}' sheet, found role(s) {', '.join(wrong) or '(none)'}")
    write_config(entries, output, comments=[f"generated from: {reference.name}"])
    without_ko = [f"{entry.functional_group} / {entry.gene_enzyme} ({entry.ko_text})" for entry in entries if not entry.ko]
    return output, len(entries), len(searchable_kos(entries)), without_ko


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build donor / acceptor KO configuration files from the reference sheets.")
    parser.add_argument("--only", choices=tuple(SHEETS), action="append",
                        help="build only this sheet (repeatable; default: all)")
    parser.add_argument("--config-dir", type=Path, default=DEFAULT_CONFIG_DIR,
                        help=f"directory holding the reference sheets and configurations (default: {DEFAULT_CONFIG_DIR})")
    parser.add_argument("-v", "--verbose", action="store_true", help="show parser warnings for rows without KO")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.WARNING if args.verbose else logging.ERROR, format="%(levelname)s: %(message)s")
    try:
        for role in args.only or list(SHEETS):
            output, n_rows, n_kos, without_ko = build(role, args.config_dir)
            print(f"{role}: wrote {output} ({n_rows} rows, {n_kos} unique KOs, {len(without_ko)} rows without KO)")
            for item in without_ko:
                print(f"  no KO: {item}")
    except (OSError, ValueError) as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
