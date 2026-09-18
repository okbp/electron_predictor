import contextlib
import io
import json
import re
import tempfile
import unittest
from pathlib import Path

from ko_detector.cli import main
from ko_detector.results_io import read_genome_taxonomy, write_genome_taxonomy
from ko_detector.models import STATUS_OK, GenomeResult
from ko_detector.taxonomy import load_genome_taxonomy, taxonomy_available

ASSEMBLY_SUMMARY = (
    "##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n"
    "#assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\n"
    "GCF_000000001.1\tP1\tB1\t\tna\t1001\t100\tAlpha beta DSM 1\n"   # strain-level taxid
    "GCF_000000002.2\tP2\tB2\t\tna\t999\t200\tGamma delta\n"         # merged taxid, other version
    "GCF_000000003.1\tP3\tB3\t\tna\t300\t300\tEpsilon zeta\n"        # lineage without family
)
# taxid, parent, rank, name
NODES = [
    ("1", "1", "no rank", "root"),
    ("2", "1", "superkingdom", "Bacteria"),
    ("10", "2", "phylum", "Phylum A"),
    ("20", "10", "class", "Class A"),
    ("30", "20", "order", "Order A"),
    ("40", "30", "family", "Family A"),
    ("50", "40", "genus", "Alpha"),
    ("60", "40", "genus", "Gamma"),
    ("70", "30", "genus", "Epsilon"),
    ("100", "50", "species", "Alpha beta"),
    ("200", "60", "species", "Gamma delta"),
    ("300", "70", "species", "Epsilon zeta"),
    ("1001", "100", "strain", "Alpha beta DSM 1"),
]


def make_taxonomy_dir(root: Path) -> Path:
    directory = root / "taxonomy"
    dump = directory / "new_taxdump"
    dump.mkdir(parents=True)
    (directory / "assembly_summary_refseq.txt").write_text(ASSEMBLY_SUMMARY, encoding="utf-8")
    (dump / "merged.dmp").write_text("999\t|\t200\t|\n", encoding="utf-8")
    parents = {taxid: parent for taxid, parent, _, _ in NODES}

    def lineage(taxid):
        chain = []
        while parents[taxid] != taxid:
            taxid = parents[taxid]
            chain.append(taxid)
        return " ".join(reversed(chain))

    (dump / "nodes.dmp").write_text(
        "".join(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\n" for taxid, parent, rank, _ in NODES), encoding="utf-8")
    (dump / "rankedlineage.dmp").write_text(
        "".join(f"{taxid}\t|\t{name}\t|\t\t|\t\t|\n" for taxid, _, _, name in NODES), encoding="utf-8")
    (dump / "taxidlineage.dmp").write_text(
        "".join(f"{taxid}\t|\t{lineage(taxid)} \t|\n" for taxid, _, _, _ in NODES), encoding="utf-8")
    return directory


class LoadTaxonomyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.directory = make_taxonomy_dir(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def test_resolves_lineages(self):
        self.assertTrue(taxonomy_available(self.directory))
        with self.assertLogs("ko_detector.taxonomy", "WARNING"):
            taxonomy = load_genome_taxonomy(
                self.directory, ["GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1", "GCF_999999999.1"])
        self.assertEqual(sorted(taxonomy), ["GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1"])

        strain = taxonomy["GCF_000000001.1"]
        self.assertEqual((strain.taxid, strain.organism_name), ("1001", "Alpha beta DSM 1"))
        self.assertEqual(strain.rank("domain"), ("Bacteria", "2"))          # superkingdom is treated as domain
        self.assertEqual(strain.rank("species"), ("Alpha beta", "100"))     # species above the strain taxid

        merged = taxonomy["GCF_000000002.1"]                                # summary lists version .2
        self.assertEqual(merged.taxid, "200")
        self.assertEqual(merged.rank("genus"), ("Gamma", "60"))

        no_family = taxonomy["GCF_000000003.1"]
        self.assertEqual(no_family.rank("family"), ("", ""))
        self.assertEqual(no_family.rank("genus"), ("Epsilon", "70"))

    def test_taxonomy_tsv_round_trip(self):
        taxonomy = load_genome_taxonomy(self.directory, ["GCF_000000001.1", "GCF_000000003.1"])
        results = [GenomeResult("GCF_000000001.1", STATUS_OK), GenomeResult("GCF_000000003.1", STATUS_OK)]
        path = self.root / "genome_taxonomy.tsv"
        write_genome_taxonomy(path, results, taxonomy)
        self.assertEqual(read_genome_taxonomy(path), taxonomy)

    def test_missing_directory(self):
        self.assertFalse(taxonomy_available(self.root / "nowhere"))


class ScanWithTaxonomyTest(unittest.TestCase):
    def test_scan_and_render_draw_the_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            directory = make_taxonomy_dir(root)
            reference = root / "ref.tsv"
            reference.write_text("Functional group\tElectron donor\tGene / enzyme\tKO\nSOB\tHS-\tsqr\tK17218\n",
                                 encoding="utf-8")
            for genome in ("GCF_000000001.1", "GCF_000000003.1"):
                (root / "in" / genome).mkdir(parents=True)
                (root / "in" / genome / f"{genome}.faa.kofam.kolist_gene.tsv").write_text(
                    "KO\tgene name\tthrshld\tscore\tE-value\tsignificant\nK17218\tWP_1\t1\t2\t1e-5\t*\n",
                    encoding="utf-8")
            out = root / "out"
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["scan", "-i", str(root / "in"), "-o", str(out), "-r", str(reference),
                             "--taxonomy-dir", str(directory), "--no-phenotype"])
            self.assertEqual(code, 0)
            self.assertEqual(sorted(read_genome_taxonomy(out / "genome_taxonomy.tsv")),
                             ["GCF_000000001.1", "GCF_000000003.1"])

            def tree(path):
                page = path.read_text(encoding="utf-8")
                data = json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>', page, re.S).group(1))
                return data["taxonomy"], {g["id"]: g for g in data["genomes"]}

            scanned, genomes = tree(out / "report.html")
            nodes = scanned["nodes"]
            # domain, phylum, class, order shared; then family A / unclassified family, genus and species for each
            self.assertEqual(len(nodes), 4 + 3 + 3)
            species = nodes[genomes["GCF_000000001.1"]["tax"]]
            self.assertEqual((species[0], species[1], species[2]), ("Alpha beta", 6, "100"))
            self.assertIn(["", 4, "", 3], nodes)                            # placeholder family under order
            self.assertEqual(genomes["GCF_000000001.1"]["org"], "Alpha beta DSM 1")

            # render-html reuses genome_taxonomy.tsv without the taxonomy directory
            rendered = root / "rendered"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["render-html", "-i", str(out), "-o", str(rendered), "--no-phenotype"]), 0)
            self.assertEqual(tree(rendered / "report_en.html")[0], scanned)
