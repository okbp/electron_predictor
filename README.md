# electron_predictor — electron donor / acceptor KO presence

Starting from reference tables of the KEGG Orthology (KO) entries used as electron donors and electron
acceptors by chemolithoautotrophic bacteria, this tool checks whether each KO is present in genomes already
annotated with KofamScan. Donors and acceptors are analysed in the same scan. It outputs per-KO genome counts and an HTML report you can browse.

- Python 3.9+, **standard library only** (no `pip install`)
- Run commands from the **repository root** (`python3 -m ko_detector …`)
- The HTML report is a single self-contained file; open it by double-clicking (no server or network needed)

## Quick start

```bash
# 0. (once, optional) Download NCBI taxonomy for the taxonomy tree — see "Setup" below
python3 scripts/setup_taxonomy.py

# 1. Build the donor and acceptor KO configurations from the reference tables
#    (first time, and whenever a reference table changes)
python3 scripts/build_configs.py

# 2. Scan the genomes for donor and acceptor KOs
python3 -m ko_detector scan \
    -i /Users/okabeppuyouko/work/GMO/refseq_reference_genomes_only_gene \
    -c config/donor_ko_config.tsv \
    -c config/acceptor_ko_config.tsv \
    -o data/results/refseq_reference_genomes

# 3. Open the report
open data/results/refseq_reference_genomes/report.html
```

Scanning the 23,434 RefSeq genomes takes about 50 seconds, plus about 5 seconds to read the taxonomy.

## Input

### Reference tables

| Sheet | File |
| --- | --- |
| Electron donors | `config/chemolithoautotroph_donor_acceptor_donor_ko_reference.tsv` |
| Electron acceptors | `config/chemolithoautotroph_donor_acceptor_acceptor_ko_reference.tsv` |

Tab-separated (fields quoted by Excel are fine). Required columns: `Functional group`, `Electron donor` or
`Electron acceptor`, `Gene / enzyme`, `KO`. That second column's header sets the role (`donor` / `acceptor`) of every row.
Optional columns: `Priority`, `Detectable via KEGG`, `Degeneracy and caveats`, `Recommended tool`.

The `KO` column is interpreted as follows.

| Written as | Interpreted as |
| --- | --- |
| `K17218` | a single KO |
| `K17229, K17230` | several KOs, comma-separated |
| `K17222-K17227` | a range, expanded (6 KOs). Expanded exactly as written; KOs not listed are not added |
| `K00368 / K15864`, `K03385, K15876 / K00362, K00363` | several KOs separated by `/`, `,` or `;` |
| `Poor KEGG coverage`, etc. | no KO. Kept in the configuration with an empty KO and not used for detection |
| `Closely related to psrABC (K08352); poor KO resolution` | a description, not a KO list: treated as no KO (with a warning), even though it mentions a KO |

Blank rows and the trailing notes (Scope, Priority, …) are skipped.
When the same KO appears in several rows (e.g. K10944–K10946 under both amoABC and pmoABC), each row gets its own column.

### KO configurations (`config/donor_ko_config.tsv`, `config/acceptor_ko_config.tsv`)

One-KO-per-row TSVs written by `scripts/build_configs.py` (or `build-config`), with a `role` column; you can edit them by hand.
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

### `scripts/build_configs.py` — both reference tables → both KO configurations

| Option | Description |
| --- | --- |
| `--only donor` / `--only acceptor` | build one sheet only (repeatable; default both) |
| `--config-dir` | directory holding the reference tables and configurations (default `config/`) |
| `-v, --verbose` | show parser warnings |

It prints the number of KOs per sheet and lists the rows without a KO. It stops with an error if a table's
header does not match its role (e.g. an acceptor table saved under the donor file name).

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
| `-c, --config` | KO configuration; repeat for donor and acceptor (either this or `-r` is required) |
| `-r, --reference` | read reference TSV(s) directly instead of a configuration file |
| `--pattern` | KofamScan file glob (default `*.kolist_gene.tsv`) |
| `--significant-only` | count only rows marked `*` |
| `--no-html` | do not write the HTML reports |
| `--taxonomy-dir` | NCBI taxonomy for the tree (default `data/taxonomy` when present) |
| `--no-taxonomy` | do not resolve lineages; reports have no tree |
| `--electron-db` | electron donor / acceptor database for the category columns (default `data/electron/electron_donor_acceptor_DATABASE_with_genome.tsv` when present) |
| `--no-electron-db` | reports without category columns |
| `--title` / `--title-en` | title of the Japanese / English report |
| `-v, --verbose` | show progress |

When finished, the number of genomes carrying each KO is printed to standard output.

### `render-html` — rebuild only the HTML from scan results

| Option | Description |
| --- | --- |
| `-i, --results-dir` | output directory of a scan (required) |
| `-c, --config` | KO configuration, repeatable (default `<results-dir>/ko_config_used.tsv`) |
| `-o, --output-dir` | where to write the reports (default `<results-dir>`) |
| `--taxonomy-dir` | read lineages from this taxonomy directory (default: `<results-dir>/genome_taxonomy.tsv`, else `data/taxonomy` when present) |
| `--no-taxonomy` | reports without the tree |
| `--electron-db` | read categories from this database (default: `<results-dir>/genome_electron_categories.tsv`, else the default database when present) |
| `--no-electron-db` | reports without category columns |
| `--title` / `--title-en` | report titles |

## Output

| File | Contents |
| --- | --- |
| `report.html` / `report_en.html` | HTML report (Japanese / English), a single file with the data embedded |
| `ko_summary.tsv` | per KO: role (`donor`, `acceptor` or `donor; acceptor`), genomes carrying it, how many of them significant, fraction |
| `genome_ko_matrix.tsv` | genome × KO gene counts (0 = absent); complexes get one column per KO |
| `genome_ko_hits.tsv` | hit details (gene ID, score, threshold, E-value, significant, role) |
| `genome_status.tsv` | per-genome status (`ok` / `multiple_files` / `no_file` / `error`) |
| `genome_electron_categories.tsv` | donor / acceptor categories per genome from the electron database (role, category, compound, consensus, confidence); only when the database was available |
| `genome_taxonomy.tsv` | organism name and lineage (domain … species, with taxids) per genome; only when taxonomy was available |
| `ko_config_used.tsv` | copy of the KO configurations used for the scan (donor and acceptor rows together) |
| `run_info.json` | run date, input and options |

- `no_file` / `error` genomes are left out of the denominator; `multiple_files` genomes are counted using the union of their files
- Genomes carrying none of the target KOs are **not shown in the HTML**, but they are included in the TSVs and in the `ko_summary.tsv` denominator

**To view the report, `report.html` alone is enough.**
To rebuild it with `render-html` you need `genome_status.tsv`, `genome_ko_hits.tsv` and `ko_config_used.tsv`
(the last one is not needed if you pass a configuration with `-c`); `run_info.json` is optional.
For the category columns it uses `genome_electron_categories.tsv` (or the database). For the taxonomy tree it also uses `genome_taxonomy.tsv`; without that file it reads `data/taxonomy` again, and without either the report has no tree.

## Using the HTML report

- **Donor / acceptor categories**: `donor_category` / `acceptor_category` from the electron donor / acceptor database
  (`data/electron/…DATABASE_with_genome.tsv`, matched by `genome_id`) appear as "Donor categories" and "Acceptor categories"
  sections left of the donor and acceptor KO sections, and collapse like them. Colours: donors — Hydrogen sky blue,
  Iron copper, Sulfur yellow, Carbon monoxide navy, Ammonia scarlet, Nitrite teal (Sulfate orange if it appears);
  acceptors — Aerobes sky blue, Nitrate-reducing teal, SRB orange, Sulfur-reducing yellow, FeRB copper,
  Methanogens / Acetogens navy; Other grey. Strong = consensus `used`, faint = `CONFLICTING` only; `not_used` is shown only in
  tooltips. "no DB" marks genomes the database does not list. Category columns are not counted in "KOs"
- **Donor / acceptor sections**: columns are split into "Electron donors" (green) and "Electron acceptors" (orange). Collapse or expand
  each with the band at the top of the headings, the narrow strip of a collapsed section, or the toolbar buttons
  (green = donors, orange = acceptors; filled when expanded, outlined when collapsed).
  KOs in a collapsed section are left out of "KOs" and the TSV export. Results without roles show a single section without the band

- **Taxonomy tree** (when taxonomy is available): rows follow the NCBI classification and a cladogram
  (domain → phylum → class → order → family → genus → species; branch lengths have no meaning) is drawn on the left,
  with phylum names and alternating phylum bands across the matrix. Hover a branch or a phylum name to see the taxon
  at that rank, its genome count and lineage; pin it to open NCBI Taxonomy. Missing ranks appear as "unclassified <rank>"
- **Overview / Detail**: Overview shows every genome on one screen. Several genomes share each screen pixel, and colour strength shows the share carrying the KO.
  Detail shows one row per genome, with genome IDs once rows are 10px or taller. Switching back to Detail returns to where you left off
- **+ / −**: change the row height (1–16px). Clicking in the overview opens the detail view at that position
- **Tooltips**: hover a cell to see the genome, the KO, its hit genes (score / threshold / E-value) and every KO the genome carries.
  Click to pin the tooltip and open NCBI Datasets / KEGG links (× or Esc to close).
  Hover a heading to see the KO description, how many genomes carry it, and caveats
- **Sorting**: click a KO heading to list genomes carrying that KO first, or "KOs" to sort by number of KOs.
  The tree is hidden while sorted this way; the "Taxonomy tree" button returns to taxonomic order
- **Filters**: search by genome ID, organism or taxon name (e.g. `Aquificota`), functional groups (which columns are shown),
  excluding below-threshold hits. The tree is pruned to the remaining genomes
- **⇩ TSV**: save the listed genomes × visible columns as 0/1 (with a BOM for Excel)
- "Items without KO", "Excluded genomes" and "About" at the top open on mouse-over

Colours: dark green = present (above threshold), light green = below-threshold hits only, grey = absent.
Acceptor columns use the same scheme in orange. `※` marks a KO shared with another functional group.

## Caveats

- Detection is **KO presence only**. Enzymes sharing KOs — amoABC vs pmoABC, nxrAB vs narGH, the direction of dsrAB
  (sulfur oxidation vs sulfate reduction) — are not distinguished. Check them separately with phylogenetic placement, HydDB, FeGenie, etc.
- A KO used on both sheets (13 KOs: amo/pmo, mmoX, coxL, narGH/nxrAB, sat, aprAB, dsrAB, cooS) is counted once, with role
  `donor; acceptor`. Its presence does not tell which direction the genome uses
- Reference items without a KEGG KO (11 donor and 11 acceptor items, such as cyc2, dsrEFH, mtrABC, qmoABC) cannot be detected ("Items without KO" in the report)
- In the HTML report, donor and acceptor KOs are separate collapsible sections (a KO on both sheets appears in both)
- `data/` is in `.gitignore`; results are not stored in the repository

## Setup: NCBI taxonomy

Downloads NCBI `new_taxdump` and `assembly_summary_refseq.txt` into `data/taxonomy/` (standard library only).
`scan` and `render-html` use them to draw the taxonomy tree; without them the report is drawn without the tree.

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
│   ├── chemolithoautotroph_donor_acceptor_donor_ko_reference.tsv      # survey results: donors (source)
│   ├── chemolithoautotroph_donor_acceptor_acceptor_ko_reference.tsv   # survey results: acceptors (source)
│   ├── donor_ko_config.tsv                                             # generated by scripts/build_configs.py
│   └── acceptor_ko_config.tsv                                          # generated by scripts/build_configs.py
├── ko_detector/
│   ├── cli.py            # subcommands
│   ├── reference.py      # reference parsing, configuration read/write
│   ├── kofam.py          # GCF directory discovery, kolist_gene.tsv reading
│   ├── summary.py        # per-KO counts
│   ├── taxonomy.py       # genome lineages from NCBI taxonomy
│   ├── results_io.py     # result TSVs: write and read back
│   ├── html_report.py    # HTML report generation
│   ├── models.py
│   └── templates/
│       ├── report.html   # report template (CSS / JS inline)
│       └── messages/     # UI strings (ja.json / en.json)
├── scripts/
│   ├── build_configs.py
│   └── setup_taxonomy.py
├── tests/
└── data/                 # in .gitignore
    ├── examples/         # sample genome
    ├── results/          # scan output
    └── taxonomy/         # downloaded by setup_taxonomy.py
```
