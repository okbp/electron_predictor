# electron_predictor — electron donor KO presence

Starting from a reference table of the KEGG Orthology (KO) entries used as electron donors by
chemolithoautotrophic bacteria, this tool checks whether each KO is present in genomes already
annotated with KofamScan. It outputs per-KO genome counts and an HTML report you can browse.

- Python 3.9+, **standard library only** (no `pip install`)
- Run commands from the **repository root** (`python3 -m ko_detector …`)
- The HTML report is a single self-contained file; open it by double-clicking (no server or network needed)

## Quick start

```bash
# 1. Build the KO configuration from the reference table (first time, and whenever the reference changes)
python3 -m ko_detector build-config \
    -r config/chemolithoautotroph_donor_acceptor_KO_reference.tsv \
    -o config/ko_config.tsv

# 2. Scan the genomes
python3 -m ko_detector scan \
    -i /Users/okabeppuyouko/work/GMO/refseq_reference_genomes_only_gene \
    -c config/ko_config.tsv \
    -o data/results/refseq_reference_genomes

# 3. Open the report
open data/results/refseq_reference_genomes/report.html
```

Scanning the 23,434 RefSeq genomes takes about 50 seconds.

## Input

### Reference table (`config/chemolithoautotroph_donor_acceptor_KO_reference.tsv`)

Tab-separated. Required columns: `Functional group`, `Electron donor`, `Gene / enzyme`, `KO`.
Optional columns: `Priority`, `Detectable via KEGG`, `Degeneracy and caveats`, `Recommended tool`.

The `KO` column is interpreted as follows.

| Written as | Interpreted as |
| --- | --- |
| `K17218` | a single KO |
| `K17229, K17230` | several KOs, comma-separated |
| `K17222-K17227` | a range, expanded (6 KOs). Expanded exactly as written; KOs not listed are not added |
| `Poor KEGG coverage`, etc. | no KO. Kept in the configuration with an empty KO and not used for detection |

Blank rows and the trailing notes (Scope, Priority, …) are skipped.
When the same KO appears in several rows (e.g. K10944–K10946 under both amoABC and pmoABC), each row gets its own column.

### KO configuration (`config/ko_config.tsv`)

A one-KO-per-row TSV written by `build-config`; you can edit it by hand.
Filling in a KO on a row whose KO is empty (e.g. hxlAB, marked `See KEGG`) makes it part of the next scan.
Lines starting with `#` are comments.

### Genomes (KofamScan results)

```
<input directory>/
├── GCF_000164905.1/
│   └── GCF_000164905.1.faa.kofam.kolist_gene.tsv
├── GCF_000005845.2/
│   └── ...
```

- Only directories starting with `GCF_` **directly under** the input directory are used (no recursion)
- In each, `*.kolist_gene.tsv` is read (change with `--pattern`)
- Six columns: `KO / gene name / thrshld / score / E-value / significant`
- **Every KO listed in the file counts as present.** Use `--significant-only` to count only rows marked `*`

## Commands

### `build-config` — reference table → KO configuration

| Option | Description |
| --- | --- |
| `-r, --reference` | reference TSV (repeatable) |
| `-o, --output` | configuration file to write |

### `scan` — detect KOs in the genomes

| Option | Description |
| --- | --- |
| `-i, --input-dir` | directory containing the `GCF_*` directories (required) |
| `-o, --output-dir` | output directory; defaults to `data/results/<input directory name>` |
| `-c, --config` | KO configuration (either this or `-r` is required) |
| `-r, --reference` | read reference TSV(s) directly instead of a configuration file |
| `--pattern` | KofamScan file glob (default `*.kolist_gene.tsv`) |
| `--significant-only` | count only rows marked `*` |
| `--no-html` | do not write the HTML reports |
| `--title` / `--title-en` | title of the Japanese / English report |
| `-v, --verbose` | show progress |

When finished, the number of genomes carrying each KO is printed to standard output.

### `render-html` — rebuild only the HTML from scan results

| Option | Description |
| --- | --- |
| `-i, --results-dir` | output directory of a scan (required) |
| `-c, --config` | KO configuration (default `<results-dir>/ko_config_used.tsv`) |
| `-o, --output-dir` | where to write the reports (default `<results-dir>`) |
| `--title` / `--title-en` | report titles |

## Output

| File | Contents |
| --- | --- |
| `report.html` / `report_en.html` | HTML report (Japanese / English), a single file with the data embedded |
| `ko_summary.tsv` | per KO: genomes carrying it, how many of them significant, fraction |
| `genome_ko_matrix.tsv` | genome × KO gene counts (0 = absent); complexes get one column per KO |
| `genome_ko_hits.tsv` | hit details (gene ID, score, threshold, E-value, significant) |
| `genome_status.tsv` | per-genome status (`ok` / `multiple_files` / `no_file` / `error`) |
| `ko_config_used.tsv` | copy of the KO configuration used for the scan |
| `run_info.json` | run date, input and options |

- `no_file` / `error` genomes are left out of the denominator; `multiple_files` genomes are counted using the union of their files
- Genomes carrying none of the target KOs are **not shown in the HTML**, but they are included in the TSVs and in the `ko_summary.tsv` denominator

**To view the report, `report.html` alone is enough.**
To rebuild it with `render-html` you need `genome_status.tsv`, `genome_ko_hits.tsv` and `ko_config_used.tsv`
(the last one is not needed if you pass a configuration with `-c`); `run_info.json` is optional.

## Using the HTML report

- **Overview / Detail**: Overview shows every genome on one screen. Several genomes share each screen pixel, and colour strength shows the share carrying the KO.
  Detail shows one row per genome, with genome IDs once rows are 10px or taller. Switching back to Detail returns to where you left off
- **+ / −**: change the row height (1–16px). Clicking in the overview opens the detail view at that position
- **Tooltips**: hover a cell to see the genome, the KO, its hit genes (score / threshold / E-value) and every KO the genome carries.
  Click to pin the tooltip and open NCBI Datasets / KEGG links (× or Esc to close).
  Hover a heading to see the KO description, how many genomes carry it, and caveats
- **Sorting**: click a KO heading to list genomes carrying that KO first, or "KOs" to sort by number of KOs
- **Filters**: genome ID search, functional groups (which columns are shown), excluding below-threshold hits
- **⇩ TSV**: save the listed genomes × visible columns as 0/1 (with a BOM for Excel)
- "Items without KO", "Excluded genomes" and "About" at the top open on mouse-over

Colours: dark green = present (above threshold), light green = below-threshold hits only, grey = absent. `※` marks a KO shared with another functional group.

## Caveats

- Detection is **KO presence only**. Enzymes sharing KOs — amoABC vs pmoABC, nxrAB vs narGH, the direction of dsrAB
  (sulfur oxidation vs sulfate reduction) — are not distinguished. Check them separately with phylogenetic placement, HydDB, FeGenie, etc.
- Reference items without a KEGG KO (11 items such as cyc2 and dsrEFH) cannot be detected ("Items without KO" in the report)
- `data/` is in `.gitignore`; results are not stored in the repository

## Setup: NCBI taxonomy

Downloads NCBI `new_taxdump` and `assembly_summary_refseq.txt` into `data/taxonomy/` (standard library only).

```bash
python3 scripts/setup_taxonomy.py                 # fetch whatever is missing (~400 MB, a few minutes)
python3 scripts/setup_taxonomy.py --force         # download the latest versions again
python3 scripts/setup_taxonomy.py --only taxdump  # one target only
```

- Existing data is never overwritten unless `--force` is given
- `new_taxdump.tar.gz` is checked against its md5 before extraction and deleted afterwards (keep it with `--keep-archive`)
- Source URL, Last-Modified and md5 are recorded in `data/taxonomy/download_info.json`

## Tests

```bash
python3 -m unittest
```

No network access is needed. Tests that use `data/examples` or the reference table are skipped when those files are absent.

## Layout

```
electron_predictor/
├── README.md
├── CLAUDE.md
├── config/
│   ├── chemolithoautotroph_donor_acceptor_KO_reference.tsv   # survey results (source)
│   └── ko_config.tsv                                          # generated by build-config
├── ko_detector/
│   ├── cli.py            # subcommands
│   ├── reference.py      # reference parsing, configuration read/write
│   ├── kofam.py          # GCF directory discovery, kolist_gene.tsv reading
│   ├── summary.py        # per-KO counts
│   ├── results_io.py     # result TSVs: write and read back
│   ├── html_report.py    # HTML report generation
│   ├── models.py
│   └── templates/
│       ├── report.html   # report template (CSS / JS inline)
│       └── messages/     # UI strings (ja.json / en.json)
├── scripts/
│   └── setup_taxonomy.py
├── tests/
└── data/                 # in .gitignore
    ├── examples/         # sample genome
    ├── results/          # scan output
    └── taxonomy/         # downloaded by setup_taxonomy.py
```
