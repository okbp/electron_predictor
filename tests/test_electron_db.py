import contextlib
import io
import json
import re
import tempfile
import unittest
from pathlib import Path

from ko_detector.cli import main
from ko_detector.electron_db import read_electron_db
from ko_detector.html_report import build_payload
from ko_detector.models import STATUS_OK, GenomeResult, Hit, KoEntry
from ko_detector.results_io import read_electron_categories, write_electron_categories

HEADER = "genome_id\torganism_name\trole\tcompound_name\tdonor_category\tacceptor_category\tconsensus\tconfidence_grade\tprimary_quote\n"
ROWS = (
    "GCF_000000001.1\tA\telectron_donor\thydrogen\tHydrogen oxidizing bacteria\t\tused\tA\t25 \u03bcg/mL\n"
    "GCF_000000001.1\tA\telectron_donor\tthiosulfate\tSulfur oxidizing bacteria\t\tCONFLICTING\tC\tq\n"
    "GCF_000000001.1\tA\telectron_acceptor\toxygen\t\tAerobes / aerobic respiration\tused\tB\tq\n"
    "GCF_000000001.1\tA\telectron_acceptor\tnitrate\t\tNitrate-reducing bacteria (denitrifiers / DNRA)\tnot_used\tD\tq\n"
    "GCF_000000002.1\tB\telectron_donor\tunknown\t\t\tused\tB\tq\n"
    "0\tbroken\telectron_donor\thydrogen\tHydrogen oxidizing bacteria\t\tused\tA\tq\n"
)


class ReadElectronDbTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def check(self, records):
        self.assertEqual(sorted(records), ["GCF_000000001.1", "GCF_000000002.1"])
        self.assertEqual(records["GCF_000000002.1"], [])          # listed, but without a category
        first = records["GCF_000000001.1"]
        self.assertEqual([(r.role, r.category, r.consensus, r.confidence) for r in first], [
            ("donor", "Hydrogen oxidizing bacteria", "used", "A"),
            ("donor", "Sulfur oxidizing bacteria", "CONFLICTING", "C"),
            ("acceptor", "Aerobes / aerobic respiration", "used", "B"),
            ("acceptor", "Nitrate-reducing bacteria (denitrifiers / DNRA)", "not_used", "D"),
        ])

    def test_reads_cp932_with_nul_characters(self):
        path = self.dir / "db.tsv"
        path.write_bytes((HEADER + ROWS).encode("cp932").replace(b"\tq\n", b"\tq\x00\n", 1))
        with self.assertLogs("ko_detector.electron_db", "INFO") as logs:
            records = read_electron_db(path)
        self.check(records)
        self.assertIn("read as cp932", "\n".join(logs.output))

    def test_reads_utf8(self):
        path = self.dir / "db.tsv"
        path.write_text(HEADER + ROWS, encoding="utf-8")
        with self.assertLogs("ko_detector.electron_db", "WARNING"):
            self.check(read_electron_db(path))

    def test_missing_columns(self):
        path = self.dir / "db.tsv"
        path.write_text("genome_id\trole\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            read_electron_db(path)

    def test_tsv_round_trip_and_payload(self):
        path = self.dir / "db.tsv"
        path.write_text(HEADER + ROWS, encoding="utf-8")
        with self.assertLogs("ko_detector.electron_db", "WARNING"):
            records = read_electron_db(path)
        results = [GenomeResult(g, STATUS_OK, hits=[Hit(g, "K17218", "WP_1", "1", "2", "1e-5", True)])
                   for g in ("GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1")]
        out = self.dir / "genome_electron_categories.tsv"
        write_electron_categories(out, results, records)
        self.assertEqual(read_electron_categories(out), records)

        entries = [KoEntry(ko="K17218", functional_group="SOB", substrate="HS-", gene_enzyme="sqr", role="donor")]
        payload = build_payload(entries, results, "ja", electron=records)
        self.assertEqual(payload["electron"]["categories"], {
            "donor": ["Hydrogen oxidizing bacteria", "Sulfur oxidizing bacteria"],
            "acceptor": ["Aerobes / aerobic respiration", "Nitrate-reducing bacteria (denitrifiers / DNRA)"],
        })
        genomes = {g["id"]: g for g in payload["genomes"]}
        self.assertEqual(genomes["GCF_000000001.1"]["el"][1], ["donor", 1, "thiosulfate", "CONFLICTING", "C"])
        self.assertEqual(genomes["GCF_000000002.1"]["el"], [])
        self.assertNotIn("el", genomes["GCF_000000003.1"])             # not in the database


class ScanWithElectronDbTest(unittest.TestCase):
    def test_scan_and_render_use_the_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            db = root / "db.tsv"
            db.write_bytes((HEADER + ROWS).encode("cp932"))
            reference = root / "ref.tsv"
            reference.write_text("Functional group\tElectron donor\tGene / enzyme\tKO\nSOB\tHS-\tsqr\tK17218\n", encoding="utf-8")
            for genome in ("GCF_000000001.1", "GCF_000000003.1"):
                (root / "in" / genome).mkdir(parents=True)
                (root / "in" / genome / f"{genome}.kolist_gene.tsv").write_text(
                    "KO\tgene name\tthrshld\tscore\tE-value\tsignificant\nK17218\tWP_1\t1\t2\t1e-5\t*\n", encoding="utf-8")
            out = root / "out"
            with contextlib.redirect_stdout(io.StringIO()), self.assertLogs("ko_detector", "WARNING"):
                self.assertEqual(main(["scan", "-i", str(root / "in"), "-o", str(out), "-r", str(reference),
                                       "--no-taxonomy", "--electron-db", str(db)]), 0)
            self.assertEqual(sorted(read_electron_categories(out / "genome_electron_categories.tsv")), ["GCF_000000001.1"])

            def electron(path):
                page = path.read_text(encoding="utf-8")
                data = json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>', page, re.S).group(1))
                return data["electron"], {g["id"]: g for g in data["genomes"]}

            scanned, genomes = electron(out / "report.html")
            self.assertEqual(scanned["categories"]["donor"], ["Hydrogen oxidizing bacteria", "Sulfur oxidizing bacteria"])
            self.assertIn("el", genomes["GCF_000000001.1"])
            self.assertNotIn("el", genomes["GCF_000000003.1"])

            # render-html reuses genome_electron_categories.tsv (the database itself is gone)
            db.unlink()
            rendered = root / "rendered"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["render-html", "-i", str(out), "-o", str(rendered), "--no-taxonomy"]), 0)
            self.assertEqual(electron(rendered / "report_en.html")[0], scanned)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["render-html", "-i", str(out), "-o", str(rendered), "--no-taxonomy", "--no-electron-db"]), 0)
            self.assertIsNone(electron(rendered / "report.html")[0])
