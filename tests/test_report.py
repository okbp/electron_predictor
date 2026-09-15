import contextlib
import io
import json
import os
import re
import tempfile
import unittest
from pathlib import Path

from ko_detector.cli import main
from ko_detector.html_report import REPORT_FILES, TEMPLATE_PATH, build_payload, load_messages, render_html
from ko_detector.models import STATUS_NO_FILE, STATUS_OK, GenomeResult, Hit, KoEntry
from ko_detector.results_io import read_results
from ko_detector.summary import summarize

ENTRIES = [
    KoEntry(ko="K10944", functional_group="AOB", substrate="NH3", gene_enzyme="amoABC", source="r:2"),
    KoEntry(ko="K10535", functional_group="AOB", substrate="NH3", gene_enzyme="hao", source="r:3"),
    KoEntry(ko="K10944", functional_group="Methanotroph", substrate="CH4", gene_enzyme="pmoABC", source="r:4"),
    KoEntry(ko="", functional_group="FeOB", substrate="Fe(II)", gene_enzyme="cyc2", ko_text="Not in KEGG", source="r:5"),
]


def hit(genome, ko, gene, significant=True):
    return Hit(genome, ko, gene, "10", "20", "1e-5", significant)


RESULTS = [
    GenomeResult("GCF_1", STATUS_OK, n_rows=10, hits=[hit("GCF_1", "K10944", "g1"), hit("GCF_1", "K10535", "g2", False)]),
    GenomeResult("GCF_2", STATUS_OK, n_rows=10, hits=[hit("GCF_2", "K10944", "g3", False)]),
    GenomeResult("GCF_3", STATUS_NO_FILE, message="no file"),
]


def embedded_data(page):
    return json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>', page, re.S).group(1))


class SummaryTest(unittest.TestCase):
    def test_summarize(self):
        summaries = {s.ko: s for s in summarize(ENTRIES, RESULTS)}
        self.assertEqual(list(summaries), ["K10944", "K10535"])
        self.assertEqual(summaries["K10944"].n_genomes, 2)
        self.assertEqual(summaries["K10944"].n_genomes_significant, 1)
        self.assertEqual(summaries["K10944"].total_genomes, 2)
        self.assertEqual(summaries["K10944"].gene_enzyme, "amoABC; pmoABC")
        self.assertEqual(summaries["K10535"].n_genomes_significant, 0)


class HtmlTest(unittest.TestCase):
    def test_messages_have_same_keys_in_every_language(self):
        keys = {lang: set(load_messages(lang)) for lang in REPORT_FILES}
        self.assertEqual(keys["ja"], keys["en"])

    def test_template_uses_only_defined_messages(self):
        template = TEMPLATE_PATH.read_text(encoding="utf-8")
        used = set(re.findall(r"data-i18n(?:-placeholder|-aria|-title)?=\"(\w+)\"", template))
        used |= set(re.findall(r"\bt\('(\w+)'", template))
        used |= set(re.findall(r"'((?:stat|col)[A-Z]\w*)'", template))
        self.assertTrue(used)
        self.assertEqual(used - set(load_messages("ja")), set())

    def test_payload_and_escaping(self):
        entries = ENTRIES + [KoEntry(ko="K00001", functional_group="</script><b>", substrate="", gene_enzyme="x")]
        payload = build_payload(entries, RESULTS, "ja", title="A & B")
        self.assertEqual(len(payload["columns"]), 4)
        self.assertTrue(payload["columns"][0]["shared"])
        self.assertEqual([g["id"] for g in payload["genomes"]], ["GCF_1", "GCF_2"])
        self.assertEqual(payload["excluded"][0]["id"], "GCF_3")
        self.assertEqual(len(payload["undetectable"]), 1)

        page = render_html(payload)
        self.assertIn('<html lang="ja">', page)
        self.assertIn("<title>A &amp; B</title>", page)
        self.assertEqual(page.count("</script>"), 2)
        self.assertEqual(embedded_data(page), payload)

    def test_genomes_without_target_kos_are_not_listed(self):
        results = RESULTS + [GenomeResult("GCF_4", STATUS_OK, n_rows=10)]
        payload = build_payload(ENTRIES, results, "ja")
        self.assertEqual([g["id"] for g in payload["genomes"]], ["GCF_1", "GCF_2"])
        self.assertEqual(payload["genomesWithoutHits"], 1)
        self.assertEqual([e["id"] for e in payload["excluded"]], ["GCF_3"])
        self.assertEqual(summarize(ENTRIES, results)[0].total_genomes, 3)

    def test_default_titles_per_language(self):
        ja = build_payload(ENTRIES, RESULTS, "ja")
        en = build_payload(ENTRIES, RESULTS, "en")
        self.assertEqual(ja["title"], "KO 保有状況レポート")
        self.assertEqual(en["title"], "KO Presence Report")
        self.assertIn('<html lang="en">', render_html(en))


class CliTest(unittest.TestCase):
    def test_build_config_scan_and_render(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "ref.tsv"
            reference.write_text(
                "Functional group\tElectron donor\tGene / enzyme\tKO\n"
                "SOB\tHS-\tsqr\tK17218\n"
                "SOB\tS2O3 2-\tsox\tK17222-K17223\n"
                "FeOB\tFe(II)\tcyc2\tPoor KEGG coverage\n",
                encoding="utf-8",
            )
            genome = root / "in" / "GCF_000000001.1"
            genome.mkdir(parents=True)
            (genome / "GCF_000000001.1.faa.kofam.kolist_gene.tsv").write_text(
                "KO\tgene name\tthrshld\tscore\tE-value\tsignificant\n"
                "K17218\tWP_1\t240.27\t363.1\t3.2e-109\t*\n"
                "K17223\tWP_2\t74.20\t35.8\t1.9e-4\t\n",
                encoding="utf-8",
            )
            (root / "in" / "GCF_000000002.1").mkdir()

            config = root / "ko_config.tsv"
            out = root / "out"
            with contextlib.redirect_stdout(io.StringIO()), self.assertLogs("ko_detector", "WARNING"):
                self.assertEqual(main(["build-config", "-r", str(reference), "-o", str(config)]), 0)
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout), self.assertLogs("ko_detector", "WARNING"):
                self.assertEqual(
                    main(["scan", "-i", str(root / "in"), "-o", str(out), "-c", str(config), "--title-en", "My report"]),
                    0,
                )
            self.assertIn("Genomes: 1 analysed, 1 excluded", stdout.getvalue())

            summary = (out / "ko_summary.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual([line.split("\t")[:5] for line in summary[1:]], [
                ["K17218", "sqr", "SOB", "HS-", "1"],
                ["K17222", "sox", "SOB", "S2O3 2-", "0"],
                ["K17223", "sox", "SOB", "S2O3 2-", "1"],
            ])
            matrix = (out / "genome_ko_matrix.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(matrix, ["genome_id\tK17218\tK17222\tK17223", "GCF_000000001.1\t1\t0\t1"])

            ja = embedded_data((out / "report.html").read_text(encoding="utf-8"))
            en = embedded_data((out / "report_en.html").read_text(encoding="utf-8"))
            self.assertEqual((ja["lang"], ja["title"]), ("ja", "KO 保有状況レポート"))
            self.assertEqual((en["lang"], en["title"]), ("en", "My report"))
            self.assertEqual(ja["genomes"], en["genomes"])

            results = {r.genome_id: r for r in read_results(out)}
            self.assertEqual(len(results["GCF_000000001.1"].hits), 2)
            self.assertFalse(results["GCF_000000002.1"].is_valid)

            cwd = os.getcwd()
            os.chdir(root)
            try:
                with contextlib.redirect_stdout(io.StringIO()), self.assertLogs("ko_detector", "WARNING"):
                    self.assertEqual(main(["scan", "-i", "in", "-c", str(config)]), 0)
            finally:
                os.chdir(cwd)
            self.assertTrue((root / "data" / "results" / "in" / "ko_summary.tsv").is_file())

            rendered = root / "rendered"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["render-html", "-i", str(out), "-o", str(rendered)]), 0)
            self.assertEqual(sorted(p.name for p in rendered.iterdir()), ["report.html", "report_en.html"])
            self.assertIn("GCF_000000001.1", (rendered / "report_en.html").read_text(encoding="utf-8"))
