import tempfile
import unittest
from pathlib import Path

from ko_detector.kofam import read_kolist, scan_directory
from ko_detector.models import STATUS_MULTIPLE_FILES, STATUS_NO_FILE, STATUS_OK

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "data" / "examples"

KOLIST = (
    "KO\tgene name\tthrshld\tscore\tE-value\tsignificant\n"
    "K17218\tWP_1\t240.27\t363.1\t3.2e-109\t*\n"
    "K17222\tWP_2\t161.27\t248.4\t8e-75\t*\n"
    "K17222\tWP_3\t161.27\t258.1\t8.8e-78\t*\n"
    "K00370\tWP_4\t835.47\t601.2\t3.8e-181\t\n"
    "K99999\tWP_5\t1.00\t2.0\t1e-3\t*\n"
)
TARGETS = {"K17218", "K17222", "K00370", "K11180"}


class ReadKolistTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "x.kolist_gene.tsv"
        self.path.write_text(KOLIST, encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def test_all_rows_count_as_present(self):
        n_rows, hits = read_kolist(self.path, TARGETS, "G")
        self.assertEqual(n_rows, 5)
        self.assertEqual([(h.ko, h.gene) for h in hits],
                         [("K17218", "WP_1"), ("K17222", "WP_2"), ("K17222", "WP_3"), ("K00370", "WP_4")])
        self.assertFalse(hits[3].significant)
        self.assertEqual(hits[0].evalue, "3.2e-109")

    def test_significant_only(self):
        _, hits = read_kolist(self.path, TARGETS, "G", significant_only=True)
        self.assertNotIn("K00370", {h.ko for h in hits})


class ScanDirectoryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        (root / "GCF_A").mkdir()
        (root / "GCF_A" / "GCF_A.faa.kofam.kolist_gene.tsv").write_text(KOLIST, encoding="utf-8")
        (root / "GCF_B").mkdir()
        (root / "GCF_B" / "other.tsv").write_text(KOLIST, encoding="utf-8")
        (root / "GCF_C").mkdir()
        (root / "GCF_C" / "1.kolist_gene.tsv").write_text(KOLIST, encoding="utf-8")
        (root / "GCF_C" / "2.kolist_gene.tsv").write_text(KOLIST, encoding="utf-8")
        (root / "GCA_D").mkdir()
        (root / "GCA_D" / "d.kolist_gene.tsv").write_text(KOLIST, encoding="utf-8")
        (root / "nested" / "GCF_E").mkdir(parents=True)
        (root / "nested" / "GCF_E" / "e.kolist_gene.tsv").write_text(KOLIST, encoding="utf-8")
        self.root = root

    def tearDown(self):
        self.tmp.cleanup()

    def test_scan(self):
        with self.assertLogs("ko_detector.kofam", "WARNING"):
            results = {r.genome_id: r for r in scan_directory(self.root, TARGETS)}
        self.assertEqual(sorted(results), ["GCF_A", "GCF_B", "GCF_C"])
        self.assertEqual(results["GCF_A"].status, STATUS_OK)
        self.assertEqual(results["GCF_B"].status, STATUS_NO_FILE)
        self.assertFalse(results["GCF_B"].is_valid)
        self.assertEqual(results["GCF_C"].status, STATUS_MULTIPLE_FILES)
        self.assertEqual(len(results["GCF_C"].hits), 4)  # duplicates across files merged
        self.assertEqual(results["GCF_C"].n_rows, 10)


@unittest.skipUnless((EXAMPLES / "GCF_000164905.1").is_dir(), "example genome not available")
class ExampleGenomeTest(unittest.TestCase):
    def test_example(self):
        targets = {"K17218", "K17222", "K03518", "K03519", "K03520", "K11180", "K00370"}
        (result,) = scan_directory(EXAMPLES, targets)
        by_ko = result.hits_by_ko()
        self.assertEqual(result.genome_id, "GCF_000164905.1")
        self.assertEqual(len(by_ko["K17222"]), 2)
        self.assertNotIn("K11180", by_ko)
        self.assertTrue({"K17218", "K03518", "K03519", "K03520", "K00370"} <= set(by_ko))
        self.assertFalse(by_ko["K00370"][0].significant)
