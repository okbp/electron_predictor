import contextlib
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path

from ko_detector.reference import read_config

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "build_configs.py"
_spec = importlib.util.spec_from_file_location("build_configs", SCRIPT)
build_configs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build_configs)

DONOR = "Functional group\tElectron donor\tGene / enzyme\tKO\nSOB\tHS-\tsqr\tK17218\nFeOB\tFe(II)\tcyc2\tPoor KEGG coverage\n"
ACCEPTOR = "Functional group\tElectron acceptor\tGene / enzyme\tKO\nO2\tO2\t\"coxABC (aa3-type, low affinity)\"\tK02274-K02276\n"


class BuildConfigsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        donor_sheet, _ = build_configs.SHEETS["donor"]
        acceptor_sheet, _ = build_configs.SHEETS["acceptor"]
        (self.dir / donor_sheet).write_text(DONOR, encoding="utf-8")
        (self.dir / acceptor_sheet).write_text(ACCEPTOR, encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, *args):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            code = build_configs.main(["--config-dir", str(self.dir), *args])
        return code, stdout.getvalue()

    def test_builds_both_configs(self):
        code, out = self.run_script()
        self.assertEqual(code, 0)
        donor = read_config(self.dir / "donor_ko_config.tsv")
        acceptor = read_config(self.dir / "acceptor_ko_config.tsv")
        self.assertEqual([(e.ko, e.role) for e in donor], [("K17218", "donor"), ("", "donor")])
        self.assertEqual([(e.ko, e.role) for e in acceptor], [("K02274", "acceptor"), ("K02275", "acceptor"), ("K02276", "acceptor")])
        self.assertIn("donor: wrote", out)
        self.assertIn("no KO: FeOB / cyc2 (Poor KEGG coverage)", out)

    def test_only_one_sheet(self):
        code, _ = self.run_script("--only", "acceptor")
        self.assertEqual(code, 0)
        self.assertTrue((self.dir / "acceptor_ko_config.tsv").is_file())
        self.assertFalse((self.dir / "donor_ko_config.tsv").exists())

    def test_rejects_a_sheet_of_the_wrong_role(self):
        donor_sheet, _ = build_configs.SHEETS["donor"]
        (self.dir / donor_sheet).write_text(ACCEPTOR, encoding="utf-8")
        with self.assertLogs("build_configs", "ERROR"):
            code, _ = self.run_script("--only", "donor")
        self.assertEqual(code, 1)
        self.assertFalse((self.dir / "donor_ko_config.tsv").exists())
