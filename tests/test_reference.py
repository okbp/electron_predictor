import tempfile
import unittest
from pathlib import Path

from ko_detector.reference import (
    is_ko_list,
    parse_ko_field,
    read_config,
    read_configs,
    read_reference,
    searchable_kos,
    write_config,
)

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "config" / "chemolithoautotroph_donor_acceptor_donor_ko_reference.tsv"
ACCEPTOR_REFERENCE = ROOT / "config" / "chemolithoautotroph_donor_acceptor_acceptor_ko_reference.tsv"

HEADER = "Functional group\tElectron donor\tGene / enzyme\tKO\tPriority\tDetectable via KEGG\tDegeneracy and caveats\tRecommended tool\n"


class ParseKoFieldTest(unittest.TestCase):
    def test_single(self):
        self.assertEqual(parse_ko_field("K17218"), ["K17218"])

    def test_comma_list(self):
        self.assertEqual(parse_ko_field("K10944, K10945, K10946"), ["K10944", "K10945", "K10946"])

    def test_range(self):
        self.assertEqual(parse_ko_field("K17222-K17227"), [f"K1722{i}" for i in range(2, 8)])

    def test_range_with_en_dash_and_mixed(self):
        self.assertEqual(parse_ko_field("K00001–K00002, K00005"), ["K00001", "K00002", "K00005"])

    def test_no_ko(self):
        for text in ("Poor KEGG coverage", "Not in KEGG", "See KEGG", ""):
            self.assertEqual(parse_ko_field(text), [])

    def test_slash_separated_groups(self):
        self.assertEqual(parse_ko_field("K03385, K15876 / K00362, K00363"), ["K03385", "K15876", "K00362", "K00363"])

    def test_is_ko_list(self):
        for text in ("K17218", "K17229, K17230", "K17222-K17227", "K00368 / K15864", "K00001; K00002 and K00003"):
            self.assertTrue(is_ko_list(text), text)
        for text in ("Closely related to psrABC (K08352); poor KO resolution", "Poor KEGG coverage", "See KEGG", ""):
            self.assertFalse(is_ko_list(text), text)

    def test_reversed_range(self):
        with self.assertRaises(ValueError):
            parse_ko_field("K00005-K00001")


class ReferenceAndConfigTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write_reference(self, body):
        path = self.dir / "ref.tsv"
        path.write_text(HEADER + body, encoding="utf-8")
        return path

    def test_read_reference_skips_notes_and_keeps_rows_without_ko(self):
        path = self.write_reference(
            "SOB\tHS-\tfccAB\tK17229, K17230\tAlternative\tYes\tsee \"literature\"\t\n"
            "SOB\tS0\tdsrEFH\tPoor KEGG coverage\tEssential\tDifficult\t\t\n"
            "\t\t\t\t\t\t\t\n"
            "Scope: sulfur donors only.\t\t\t\t\t\t\t\n"
        )
        with self.assertLogs("ko_detector.reference", "WARNING"):
            entries = read_reference(path)
        self.assertEqual([e.ko for e in entries], ["K17229", "K17230", ""])
        self.assertEqual(entries[0].source, "ref.tsv:2")
        self.assertEqual(entries[2].ko_text, "Poor KEGG coverage")
        self.assertEqual({e.role for e in entries}, {"donor"})
        self.assertEqual(searchable_kos(entries), ["K17229", "K17230"])

    def test_config_round_trip(self):
        path = self.write_reference("SOB\tHS-\tfccAB\tK17229\tAlternative\tYes\tsee \"literature\"; a\tb\n")
        entries = read_reference(path)
        config = self.dir / "config.tsv"
        write_config(entries, config, comments=["test"])
        self.assertEqual(read_config(config), entries)

    def test_acceptor_sheet_role_and_description_ko_field(self):
        path = self.dir / "acceptor.tsv"
        path.write_text(
            "Functional group\tElectron acceptor\tGene / enzyme\tKO\n"
            "Nitrate-reducing\tNO2-\tnirK / nirS\tK00368 / K15864\n"
            "SRB\tS2O3 2-\tphsABC\tClosely related to psrABC (K08352); poor KO resolution\n",
            encoding="utf-8",
        )
        with self.assertLogs("ko_detector.reference", "WARNING") as logs:
            entries = read_reference(path)
        self.assertEqual([(e.ko, e.role) for e in entries], [("K00368", "acceptor"), ("K15864", "acceptor"), ("", "acceptor")])
        self.assertIn("is a description", "\n".join(logs.output))

    def test_read_configs_merges_files(self):
        donor = read_reference(self.write_reference("SOB\tHS-\tsqr\tK17218\n"))
        first, second = self.dir / "a.tsv", self.dir / "b.tsv"
        write_config(donor, first)
        write_config(donor, second)
        self.assertEqual(read_configs([first, second]), donor)

    def test_read_config_rejects_invalid_role(self):
        config = self.dir / "config.tsv"
        config.write_text("ko\trole\tfunctional_group\tgene_enzyme\nK17218\tdonr\tSOB\tsqr\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            read_config(config)

    def test_read_config_rejects_invalid_ko(self):
        config = self.dir / "config.tsv"
        config.write_text("ko\tfunctional_group\tgene_enzyme\nK123\tSOB\tsqr\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            read_config(config)


@unittest.skipUnless(REFERENCE.is_file(), "reference sheet not available")
class RealReferenceTest(unittest.TestCase):
    def test_real_reference(self):
        with self.assertLogs("ko_detector.reference", "WARNING"):
            entries = read_reference(REFERENCE)
        kos = [e.ko for e in entries]
        self.assertIn("K17225", kos)
        self.assertNotIn("K22622", kos)
        self.assertEqual(kos.count("K10944"), 3)
        self.assertFalse(any(e.functional_group.startswith("Scope:") for e in entries))
        self.assertEqual({e.role for e in entries}, {"donor"})


@unittest.skipUnless(ACCEPTOR_REFERENCE.is_file(), "acceptor reference sheet not available")
class RealAcceptorReferenceTest(unittest.TestCase):
    def test_real_acceptor_reference(self):
        with self.assertLogs("ko_detector.reference", "WARNING"):
            entries = read_reference(ACCEPTOR_REFERENCE)
        self.assertEqual({e.role for e in entries}, {"acceptor"})
        by_gene = {}
        for e in entries:
            by_gene.setdefault(e.gene_enzyme, []).append(e.ko)
        self.assertEqual(by_gene["narGHJI (cytoplasmic face)"], ["K00370", "K00371", "K00372", "K00373", "K00374"])
        self.assertEqual(by_gene["nirK (Cu-type) / nirS (cd1-type)"], ["K00368", "K15864"])
        self.assertEqual(by_gene["phsABC (thiosulfate reductase)"], [""])
        self.assertEqual([e.gene_enzyme for e in entries if e.ko == "K08352"], ["psrABC (polysulfide reductase)"])
        self.assertFalse(any(e.functional_group.startswith(("Scope:", "Priority:", "Note ", "In methanogens")) for e in entries))
