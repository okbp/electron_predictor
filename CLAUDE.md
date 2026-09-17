# CLAUDE.md — electron_predictor (electron donor / acceptor KO presence)

Guidance for Claude Code working in this repository. User-facing instructions are in `README.md`.

## Purpose

Check whether the KOs in the reference tables of electron donors and electron acceptors used by chemolithoautotrophic
bacteria (`config/chemolithoautotroph_donor_acceptor_{donor,acceptor}_ko_reference.tsv`) are present in genomes annotated
with KofamScan, and output per-KO genome counts plus an HTML report.
The main target is 23,434 RefSeq genomes (`/Users/okabeppuyouko/work/GMO/refseq_reference_genomes_only_gene`).

## Settled decisions (set by the user — do not change)

- **Standard library only.** Do not add external packages or assume `pip install`
  (a `pyproject.toml` was added once and removed at the user's request).
- **The package lives at `ko_detector/` in the repository root.** Do not switch to a `src/` layout
  (it was moved there once and reverted, because `python3 -m ko_detector` stopped working from the root).
  Run everything from the root (`python3 -m ko_detector …` / `python3 -m unittest`).
- **Python 3.9 compatible** (the local interpreter is 3.9.6). Use `from __future__ import annotations`;
  no runtime-evaluated `X | Y` types, no `match`, no 3.10+ APIs.
- **Every KO listed in a KofamScan file counts as present**, even when `significant` is empty.
  `--significant-only` is an option.
- **GCF directories are searched only directly under the input directory** (no recursion).
- **KO ranges are expanded exactly as written in the reference.** Do not add KOs that are not listed
  (e.g. do not add soxD, K22622, to `K17222-K17227`).
- **Detection is KO presence only.** No phylogenetic resolution of shared KOs (amo/pmo, nxr/nar, dsr direction).
  Complexes are reported per KO; there is currently no "is the whole complex present" call.
- **Donors and acceptors are analysed in one scan** (`scan -c config/donor_ko_config.tsv -c config/acceptor_ko_config.tsv`).
  Every configuration row has a `role` (`donor` / `acceptor`) taken from the reference header (`Electron donor` /
  `Electron acceptor`). A KO on both sheets is detected once and reported with role `donor; acceptor`.
- **The report splits columns into collapsible role sections** (electron donors / electron acceptors), because the combined
  table is too wide. Both start expanded; the top header band, a collapsed section's stub, and toolbar buttons toggle them.
- Always write **two reports: Japanese `report.html` and English `report_en.html`**.
- Even at ~20,000 genomes the report shows **everything on one page without pagination**, compactly, with details in tooltips.
- Genomes carrying none of the target KOs are **not shown in the HTML** (only their count). They stay in the TSVs and in the denominators.
- **Taxonomy tree like `../autotrophic_genomes/viewer`**: NCBI classification (domain → phylum → class → order → family → genus → species),
  no distance calculation. Drawn on the report canvas, **not with d3 / SVG** (~23,000 nodes + 18,242 leaves would be too many DOM elements).
  Rows follow taxonomic order by default; sorting by a KO or KO count hides the tree and the "Taxonomy tree" button restores it.
  Missing ranks become "unclassified" placeholder nodes so every genome sits at the same depth.

## Modules

| File | Role |
| --- | --- |
| `cli.py` | `build-config` / `scan` / `render-html`. Default `-o` is `data/results/<input directory name>` and default taxonomy is `data/taxonomy` when present (both relative to the working directory) |
| `reference.py` | reference table → `KoEntry`; KO configuration (TSV) read/write |
| `kofam.py` | `GCF_*` discovery, `*.kolist_gene.tsv` reading → `GenomeResult` |
| `summary.py` | per-KO genome counts (`KoSummary`) |
| `electron_db.py` | donor / acceptor categories per genome (`ElectronRecord`) from `data/electron/electron_donor_acceptor_DATABASE_with_genome.tsv` |
| `taxonomy.py` | genome lineages (`GenomeTaxonomy`) from the assembly summary and new_taxdump in `data/taxonomy` |
| `results_io.py` | writes result TSVs / `run_info.json`, and reads them back for `render-html` |
| `html_report.py` | embeds the payload (JSON) and UI strings into the template |
| `models.py` | `KoEntry` / `Hit` / `GenomeResult` and status constants |
| `templates/report.html` | the report (CSS / JS inline) |
| `templates/messages/{ja,en}.json` | UI strings |
| `scripts/build_configs.py` | builds `config/donor_ko_config.tsv` and `config/acceptor_ko_config.tsv` from the two reference tables (fixed names in `SHEETS`) |
| `scripts/setup_taxonomy.py` | downloads NCBI new_taxdump and assembly_summary_refseq.txt (independent of ko_detector) |

## Reference parsing (`reference.py`)

- Columns are looked up by header name (`REFERENCE_COLUMNS`). `Electron donor` / `Electron acceptor` / `Substrate`
  all map to `substrate`; the header that matched sets `role` (`ROLE_BY_SUBSTRATE_HEADER`; `Substrate` gives an empty role).
  `read_config` rejects roles other than `donor` / `acceptor` / empty. `scripts/build_configs.py` refuses a table whose
  role does not match the file it is built into.
- The acceptor table is exported by Excel: fields with commas are quoted. The csv module handles this; do not split on tabs by hand.
- Rows with an empty `Functional group`, or with both `Gene / enzyme` and `KO` empty, are skipped (blank rows, trailing notes).
- `KO_TOKEN_RE` picks up single KOs and ranges in order. Ranges must be fewer than 100 wide; reversed ranges are errors.
  Separators between KOs may be `,` `/` `;` `and` `or` (`K00368 / K15864`, `K03385, K15876 / K00362, K00363`).
- **A KO cell is used only if it is purely a KO list** (`is_ko_list`). A description that mentions a KO —
  phsABC's `Closely related to psrABC (K08352); poor KO resolution` — becomes an empty KO with a warning; otherwise
  K08352 (psrABC's own KO) would be counted as phsABC.
- Rows yielding no KO are **kept in the configuration with an empty KO** and a warning, so a KO can be filled in by hand.
  `scan` ignores them; the HTML lists them under "Items without KO".
- A row with several KOs becomes several entries. The same KO in several rows stays that way (K10944–6 appear in 3 rows).
  Aggregates (`ko_summary.tsv` / `genome_ko_matrix.tsv`) are per KO, with annotations joined by `; `.
- The configuration is written with the csv module (caveats containing `"` get quoted). Lines starting with `#` are comments.

## Scan (`kofam.py` / `summary.py`)

- If a header line exists (first column `KO`), columns are looked up by name; otherwise `KOFAM_COLUMNS` order is assumed.
- `scan` / `render-html` take `-c` repeatedly; `read_configs` concatenates the files and drops rows repeated verbatim.
  The union of KOs is detected in a single pass over the genomes. `ko_config_used.tsv` stores all rows together.
- The same `(KO, gene)` is one hit. Several matching files in one genome are merged (`multiple_files`).
- Status: `ok` / `multiple_files` are counted; `no_file` / `error` are left out of the denominator and the scan continues.
- Summaries always include KOs found in no genome.

## Taxonomy (`taxonomy.py`)

- Inputs (downloaded by `scripts/setup_taxonomy.py`): `assembly_summary_refseq.txt` (accession → taxid, organism name),
  `new_taxdump/merged.dmp` (old → current taxid), `taxidlineage.dmp` (ancestors), `nodes.dmp` (rank of each taxid),
  `rankedlineage.dmp` (name of each taxid, column 2 only). Same approach as the viewer's `build-tree-all.mjs`.
- **Do not take species from the `rankedlineage.dmp` species column**: it is empty when the taxid itself is a species
  (15,791 of 18,242 genomes). Walk the ancestors (plus the taxid itself) and pick the node whose rank is each of the 7 ranks.
- Dump fields are separated by `\t|\t`; `taxidlineage.dmp` lists ancestors space-separated (root first, with a trailing space).
  The top rank is `domain` in current dumps; older dumps say `superkingdom` (`RANK_ALIASES`).
- Each dump is streamed and filtered by taxid (first column checked before splitting), so memory stays small.
  23,434 genomes resolve in about 3–5 seconds. If an accession is missing, another version of the same accession is used.
- `scan` writes `genome_taxonomy.tsv`. `render-html` prefers that file, then `--taxonomy-dir` / `data/taxonomy`;
  `--no-taxonomy` disables the tree. Tests pass `--no-taxonomy` so they never read the real `data/taxonomy`.

## Electron donor / acceptor database (`electron_db.py`)

- One row per genome × role × compound. `role` is `electron_donor` / `electron_acceptor`; the category is in
  `donor_category` / `acceptor_category` respectively. Joined to scan results by `genome_id` (358 of the 586 autotroph genomes).
- **The file is saved by Excel as cp932, not UTF-8** (e.g. μ is `83 CA`) and contains NUL characters. `read_electron_db`
  tries UTF-8 then cp932 and strips NULs; do not assume UTF-8. Rows whose `genome_id` is not a GCF/GCA accession (one row has `0`) are skipped.
- Every listed genome is a key even without categories (so "not in the database" and "no category" stay distinct).
- `scan` writes `genome_electron_categories.tsv` (UTF-8; a listed genome without categories gets an empty row).
  `render-html` prefers that file, then `--electron-db` / the default path; `--no-electron-db` disables it. Tests pass
  `--no-electron-db` so they never read the real database.

## Output / read-back contract (`results_io.py`)

`render-html` rebuilds `GenomeResult`s from `genome_status.tsv` (all genomes) + `genome_ko_hits.tsv` (hits) +
`ko_config_used.tsv` (or `-c`). `run_info.json` is optional; `genome_taxonomy.tsv` is optional (tree).
**When renaming columns or files, change both the writer and the reader** (`CliTest` checks the round trip).
`genome_ko_matrix.tsv` / `ko_summary.tsv` are for users and are not read back.

## HTML report (`html_report.py` / `templates/report.html`)

### Self-contained single file

- The template placeholders `__REPORT_LANG__` / `__REPORT_TITLE__` / `__REPORT_DATA__` are replaced.
  Data goes into `<script type="application/json">` with **every `<` replaced by its JSON Unicode escape
  (backslash + `u003c`)** to guard against `</script>`. `<` only occurs inside JSON strings, so the value is unchanged.
- No external scripts, CSS, fonts or fetch. The report must open from `file://` by double-clicking.
- Payload: `columns` (configuration rows that have a KO, in display order), `genomes` (**only genomes with hits**;
  `hits: {KO: [[gene, score, threshold, evalue, significant], …]}`), `genomesWithoutHits` (count),
  `excluded`, `undetectable`, `i18n` (UI strings), `taxonomy` (null, or `{ranks, nodes}` with nodes
  `[name, rank index, taxid, parent index]`; each genome then also carries `tax` = its species node, `org` and `taxid`).
  Nodes are identified by `(parent, taxid)`; placeholders have empty name and taxid and are shared under the same parent.

### UI strings (i18n)

- UI text lives in `messages/ja.json` and `en.json`; the template refers to it with `t('key', {…})` or
  `data-i18n` / `data-i18n-placeholder` / `data-i18n-aria` / `data-i18n-title`.
- **Keep the keys identical in both languages.** Every key the template uses must exist (checked by `HtmlTest`;
  write the first argument of `t(` as a literal — a variable escapes the check).
- To add a language, add one JSON file and one entry in `REPORT_FILES`.

### Rendering (canvas)

- Do not render tens of thousands of rows × 44 columns as a DOM table (it cannot meet the no-pagination requirement).
  Two canvases are used, `#head` for headings and `#body` for the matrix; scrolling is handled by a transparent
  `#scroller` whose only child `#sizer` provides the size. **Only visible rows are drawn.**
- Presence is expanded at start-up into `levels` (`Uint8Array`, genome × KO; 0 absent / 1 below-threshold only / 2 above threshold).
  Filtering and sorting rebuild `view.order` (the list of genome indices shown) in `recompute()`.
- **Overview** (`state.rowPx == null`): when rows would be thinner than 1px, each screen pixel (bin) aggregates its
  genomes and the share carrying the KO sets the opacity (`binRange()`). If everything fits at ≥1px, rows are drawn
  normally, up to 16px.
- **Detail**: `ROW_STEPS` (1–16px). Genome IDs and KO counts are drawn from `LABEL_MIN_PX` (10px).
- Overview / Detail buttons: switching to Overview stores the centre genome and row height in `lastDetail`; Detail restores it.
  **A stored row height below 10px is replaced by 12px (`CLICK_ZOOM_PX`)** — otherwise shrinking to 1px with "−" and
  switching back reopened an unreadable 1px view.
- Column width `CELL` 10px, gap 2px, 8px between functional groups. Headings: functional group (two staggered rows)
  → gene → KO (vertical) → share bar. Geometry constants are at the top of the script.
- Colours are CSS custom properties (light / dark), cached by `css()` and used on the canvas.
  Check dark mode after changing colours. The page reloads when the colour scheme changes.

### Taxonomy tree

- Built once at start-up: children sorted by name (placeholders last), pre-order DFS gives `taxPos` (row order) and `tPre`.
- `layoutTree(order)` runs in `recompute()` whenever the tree is on: every clade is a contiguous row range
  `[tFirst, tLast]`; a node's centre is the midpoint of its first and last visible children (like `d3.cluster`).
  Filtering therefore prunes the tree automatically.
- Layout on the left: phylum names (`TREE_LABEL_W`), tree (`TREE_W`, one step per rank, domain stub like the viewer),
  then genome IDs / KO counts, then the matrix. Alternating phylum bands span the whole width; a line separates domains.
  Rank names are written vertically in the header.
- **Link opacity depends on on-screen density, per rank**: each rank is one path with alpha = average node spacing / 3px
  (capped at 1, at least 0.2 when drawn), and ranks below 0.08 are skipped. In the overview this drops species and genus
  and fades family; at 12px rows every rank is fully drawn. Genome leaf links are drawn only in the detail view.
  A fixed fade does not work: a single canvas path is composited once, so ~15,000 species links in 700px still render as a
  solid grey block at any constant alpha.
- Hit-testing in the tree area: rank = `ceil((x − tree origin) / TREE_STEP)` (the segment between rank ticks d−1 and d leads
  to the node at rank d, as in the viewer); the phylum label column counts as phylum. The hovered / pinned path is drawn bold.
- Measured with 18,242 genomes: load ~155 ms; tree layout + overview draw ~13 ms; detail / overview switch ~20 ms.

### Category sections

- Payload `electron.categories` (`{donor: [...], acceptor: [...]}`) and per genome `el`: `[role, category index, compound,
  consensus, confidence]`; no `el` = not in the database.
- The template turns categories into columns (`c.cat`) of `donor_category` / `acceptor_category` sections inserted left of the
  donor / acceptor KO sections, right after the payload is parsed (before `HAS_ROLES`, `KOS` etc. are derived). Each has a
  final "no DB" column (level 2 for genomes without `el`).
- Levels: `used` = 2, only `CONFLICTING` = 1 (drawn at 40% opacity), `not_used` = 0 (tooltip only). `lv()` does not hide
  level 1 for categories under "exclude below-threshold hits".
- **Colours and column order are the user's choice**, matched by keyword in `CATEGORY_RULES`; the list order is the column order.
  Donor: Sulfur (sulfur, 黄色) → Ammonia (緋色) → Nitrite (nitr, 青緑色) → Hydrogen (水色) → Carbon monoxide (carbon, 紺色)
  → Iron (銅色) → [sulfate, オレンジ, not in the data yet] → Other (grey) → no DB.
  Acceptor: Aerobes (aerob, 水色) → Nitrate-reducing (青緑色) → SRB (sulfate, オレンジ) → Sulfur-reducing (黄色)
  → Methanogens / Acetogens (紺色) → FeRB (iron, 銅色) → Other → no DB.
  Keywords are checked top-down, so none may be a substring of another category name (currently "sulfur" ∉ "sulfate",
  "iron" ∉ "denitrifiers"); check this when adding a keyword. The donor sheet has "Nitrite oxidizing bacteria" but no
  Nitrate/Sulfate category; nitrite was given the nitrate colour.
- **Category column headings show the category name** (user request), not the short element label: written vertically in
  SANS, truncated with "…" to fit the group + gene + KO rows (`CAT_LABEL_H`). Category sections therefore draw no group label,
  bracket or gene label, and the whole height below the band counts as the column in header hit-testing. Only the "no DB"
  column keeps its short label. The rule's short label (third field of `CATEGORY_RULES`) is no longer displayed.
- Category columns are excluded from "KOs" (`view.visK`), `KO_COUNT` and KO chips; they have their own tooltip rows.

### Role sections

- `HAS_ROLES` (any column has a `role`) adds a `HEAD_ROLE` band at the top of the header; results without roles
  (older `ko_config_used.tsv`) get no band and one section.
- `recompute()` lays out one section per role in `ROLES` order (config order): columns with `GROUP_GAP` between functional
  groups, `ROLE_GAP` between sections (separator line in the middle). A collapsed section becomes a `STUB_W` strip and
  its columns leave `view.cols`, so "KOs", sorting counts and the TSV export only cover expanded sections.
- **Acceptor columns are orange, donor columns green** (user request): `--hit-acc` / `--weak-acc` / `--bar-acc`
  (light and dark) are used for cells, the overview aggregation, header share bars and tooltip chips (`.kc i.sa` / `.wa`)
  whenever `column.role === 'acceptor'`. The legend's orange swatches and "green = donors / orange = acceptors" note are
  `.role-only` and hidden without roles.
- The toolbar section toggles (`#role-toggles`) are coloured by role too (user request): filled with `--btn-don` / `--btn-acc`
  when expanded, outlined in that colour when collapsed. `--btn-acc` (light #b9560c) is darker than the cell orange so white
  text stays readable; dark mode uses the bright cell colours with dark text.
- Hit-testing: header y < `HEAD_ROLE`, or anywhere over a stub (header or body), is the section toggle (`sectionAtX`).
  All other header rows are offset by `HEAD_ROLE`; keep that offset when touching header geometry.

### Tooltips

- Row and column are **computed from coordinates** (`bodyTarget()` / `headTarget()`); no per-cell elements or listeners.
- Hover shows a following tooltip; clicking a cell pins it (`.pinned`, NCBI Datasets / KEGG links); ×, Esc or clicking the same cell unpins.
  A click in the overview does not pin — it switches to the detail view at that position.
- Unpinned tooltips close on scroll.
- Sorting (heading click) resets `scrollTop` to 0, otherwise the top of the new order is off screen.
- `.kc` (KO chip) is `inline-flex`, so whitespace between its children is not rendered. Space them with margins
  (they once rendered as "sqrK17218").

## `scripts/setup_taxonomy.py`

- The default destination is `data/taxonomy` relative to the script location (independent of the working directory).
- Existing data is left alone unless `--force` (`new_taxdump/` counts as present when `nodes.dmp` / `names.dmp` / `rankedlineage.dmp` exist).
- Downloads go through `.part` files with 3 retries (not for HTTP errors). The taxdump is verified against its `.md5`,
  extracted into a temporary directory, checked, then swapped in. Tar members with absolute paths, `..`, or types other
  than regular files / directories are rejected.
- `download_info.json` records URL / Last-Modified / size / md5 / download time.

## Verifying changes

```bash
python3 -m unittest                                   # 40 tests, no network
python3 scripts/build_configs.py                      # donor 38 KOs / acceptor 59 KOs, 11 rows without KO each
python3 -m ko_detector scan -i data/examples -c config/donor_ko_config.tsv -c config/acceptor_ko_config.tsv -o <scratch dir>
python3 -m ko_detector render-html -i data/results/refseq_reference_genomes -o <scratch dir>
```

**After touching the template, render the 23,434-genome results and look at them** (unit tests only check UI string keys and the payload).

- JS syntax: extract the `<script>` from the generated HTML and run `node --check`.
- Appearance: headless Chrome `--screenshot`. For interactions, inject a `window.onerror` hook and an action script
  (button `click()`, dispatched `MouseEvent`s) into a copy of the generated HTML, store results as attributes on
  `documentElement`, and read them with `--dump-dom`.
- Headless Chrome sometimes does not exit: use a separate `--user-data-dir` per run and a `pkill` timer.
  A `mousemove` sent right after scrolling is undone by the scroll event closing the tooltip, so wait a moment first.
- Write verification files to a scratch directory; do not overwrite `data/results/`.

## Pitfalls

- `data/` is in `.gitignore`; results under `data/results/` are not kept in the repository
  (`data/results/` once disappeared for reasons outside the session; a scan takes about 50 seconds to regenerate).
- The default `scan` output directory depends on the input directory name, so the same input goes to different places with and without `-o`.
- KofamScan files are `*.kolist_gene.tsv`; the `_kolist_gene.tsv` spelling does not match the default pattern.
- Counts include over-estimates from shared KOs (e.g. most of the ~5,200 nxrAB genomes are likely narGH).
