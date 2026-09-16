#!/usr/bin/env bash
# Copy only the KofamScan result files (*kolist_gene.tsv) of every GCF_* directory, keeping the
# GCF_* directory structure.
#
#   scripts/copy_kolist.sh SRC_DIR DEST_DIR [--dry-run]
#
# - Only GCF_* directories directly under SRC_DIR are used (same rule as ko_detector; no recursion).
# - GCF_* directories without a *kolist_gene.tsv file are not created in DEST_DIR.
# - Uses rsync when available (re-running copies only changed files); otherwise find + cp.
#   Set COPY_KOLIST_NO_RSYNC=1 to force find + cp.
set -euo pipefail

usage() { echo "usage: $0 SRC_DIR DEST_DIR [--dry-run]" >&2; exit 2; }
[ $# -ge 2 ] && [ $# -le 3 ] || usage
src=${1%/}
dest=${2%/}
dry_run=${3:-}
[ -z "$dry_run" ] || [ "$dry_run" = "--dry-run" ] || usage
[ -d "$src" ] || { echo "error: source directory not found: $src" >&2; exit 1; }

list_files() {
  find "$src" -mindepth 2 -maxdepth 2 -type f -path "$src/GCF_*/*kolist_gene.tsv"
}

n_files=$(list_files | wc -l | tr -d ' ')
echo "found $n_files *kolist_gene.tsv file(s) in $src/GCF_*/"

if [ -n "$dry_run" ]; then
  list_files | head -n 5 | sed "s|^$src/|  would copy: $dest/|"
  [ "$n_files" -le 5 ] || echo "  ... and $((n_files - 5)) more"
  exit 0
fi

mkdir -p "$dest"
if command -v rsync >/dev/null 2>&1 && [ -z "${COPY_KOLIST_NO_RSYNC:-}" ]; then
  # Anchored filters: GCF_* directories at the top, their *kolist_gene.tsv files, nothing else.
  rsync -a --prune-empty-dirs \
    --include='/GCF_*/' --include='/GCF_*/*kolist_gene.tsv' --exclude='*' \
    "$src/" "$dest/"
else
  list_files | while IFS= read -r file; do
    rel=${file#"$src"/}
    mkdir -p "$dest/${rel%/*}"
    cp -p "$file" "$dest/$rel"
  done
fi

n_copied=$(find "$dest" -mindepth 2 -maxdepth 2 -type f -path "$dest/GCF_*/*kolist_gene.tsv" | wc -l | tr -d ' ')
echo "done: $dest now holds $n_copied *kolist_gene.tsv file(s)"
