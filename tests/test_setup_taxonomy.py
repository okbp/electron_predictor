import hashlib
import importlib.util
import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "setup_taxonomy.py"
_spec = importlib.util.spec_from_file_location("setup_taxonomy", SCRIPT)
setup_taxonomy = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(setup_taxonomy)

TAXDUMP_FILES = {"nodes.dmp": "1\t|\t1\t|\n", "names.dmp": "1\t|\troot\t|\n", "rankedlineage.dmp": "1\t|\troot\t|\n"}
ASSEMBLY_SUMMARY = "##  See README\n#assembly_accession\tbioproject\nGCF_000005845.2\tPRJNA57779\n"


class SetupTaxonomyTest(unittest.TestCase):
    """Runs the script against local file:// URLs, so no network access is needed."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.remote = root / "remote"
        self.remote.mkdir()
        self.out = root / "taxonomy"
        self.archive = self.remote / "new_taxdump.tar.gz"
        self.assembly = self.remote / "assembly_summary_refseq.txt"
        self.assembly.write_text(ASSEMBLY_SUMMARY, encoding="utf-8")
        self.make_archive(TAXDUMP_FILES)

    def tearDown(self):
        self.tmp.cleanup()

    def make_archive(self, files, md5=None):
        with tarfile.open(self.archive, "w:gz") as archive:
            for name, text in files.items():
                data = text.encode("utf-8")
                info = tarfile.TarInfo(name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
        digest = md5 or hashlib.md5(self.archive.read_bytes()).hexdigest()
        Path(f"{self.archive}.md5").write_text(f"{digest}  new_taxdump.tar.gz\n", encoding="utf-8")

    def run_setup(self, *extra):
        argv = [
            "-o", str(self.out),
            "--taxdump-url", self.archive.as_uri(),
            "--assembly-summary-url", self.assembly.as_uri(),
            *extra,
        ]
        with self.assertLogs("setup_taxonomy", "INFO") as logs:
            code = setup_taxonomy.main(argv)
        return code, "\n".join(logs.output)

    def test_downloads_and_extracts(self):
        code, _ = self.run_setup()
        self.assertEqual(code, 0)
        self.assertEqual((self.out / "new_taxdump" / "names.dmp").read_text(encoding="utf-8"), TAXDUMP_FILES["names.dmp"])
        self.assertEqual((self.out / "assembly_summary_refseq.txt").read_text(encoding="utf-8"), ASSEMBLY_SUMMARY)
        self.assertEqual(sorted(p.name for p in self.out.iterdir()),
                         ["assembly_summary_refseq.txt", "download_info.json", "new_taxdump"])
        info = json.loads((self.out / "download_info.json").read_text(encoding="utf-8"))
        self.assertEqual(info["taxdump"]["md5"], hashlib.md5(self.archive.read_bytes()).hexdigest())
        self.assertEqual(info["assembly_summary"]["url"], self.assembly.as_uri())

    def test_existing_data_is_kept_unless_forced(self):
        self.assertEqual(self.run_setup()[0], 0)
        names = self.out / "new_taxdump" / "names.dmp"
        names.write_text("local edit\n", encoding="utf-8")

        code, logs = self.run_setup()
        self.assertEqual(code, 0)
        self.assertIn("skip taxdump", logs)
        self.assertEqual(names.read_text(encoding="utf-8"), "local edit\n")

        code, _ = self.run_setup("--force", "--only", "taxdump", "--keep-archive")
        self.assertEqual(code, 0)
        self.assertEqual(names.read_text(encoding="utf-8"), TAXDUMP_FILES["names.dmp"])
        self.assertTrue((self.out / "new_taxdump.tar.gz").is_file())
        self.assertFalse((self.out / "new_taxdump.old").exists())

    def test_md5_mismatch_leaves_nothing_behind(self):
        self.make_archive(TAXDUMP_FILES, md5="0" * 32)
        code, logs = self.run_setup("--only", "taxdump")
        self.assertEqual(code, 1)
        self.assertIn("md5 mismatch", logs)
        self.assertEqual(list(self.out.iterdir()), [])

    def test_unsafe_archive_member_is_rejected(self):
        self.make_archive(dict(TAXDUMP_FILES, **{"../evil.dmp": "x"}))
        code, logs = self.run_setup("--only", "taxdump")
        self.assertEqual(code, 1)
        self.assertIn("unsafe path", logs)
        self.assertFalse((self.out.parent / "evil.dmp").exists())
        self.assertEqual(list(self.out.iterdir()), [])

    def test_parse_md5(self):
        self.assertEqual(setup_taxonomy.parse_md5("4DA5A5C11F0CABFB6F5E659DB5E90C50  new_taxdump.tar.gz\n"),
                         "4da5a5c11f0cabfb6f5e659db5e90c50")
        with self.assertRaises(ValueError):
            setup_taxonomy.parse_md5("<html>not found</html>")
