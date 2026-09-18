import contextlib
import io
import json
import re
import tempfile
import unittest
from pathlib import Path

from ko_detector.cli import main
from ko_detector.html_report import build_payload
from ko_detector.models import STATUS_OK, GenomeResult, Hit, KoEntry
from ko_detector.phenotype import parse_number, read_phenotype, read_tables, token_sources
from ko_detector.results_io import read_electron_categories, read_phenotypes, write_electron_categories, write_phenotypes

COLUMNS = (
    "organism_id", "organism", "genome_accession", "oxygen", "oxygen_status", "oxygen_source_tier", "oxygen_other_values",
    "optimal_temperature", "optimal_temperature_status", "optimal_temperature_source_tier", "optimal_pH", "optimal_pH_kind",
    "habitat", "habitat_status", "electron_donors", "electron_donors_source", "electron_acceptors", "electron_acceptors_source",
    "electron_donor_acceptor_conflicting_reports", "reported_not_used_as_donor", "reported_not_used_as_acceptor",
)
G1, G2, G3 = "GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1"
ROWS = (
    dict(organism_id="org1", organism="Alpha one", genome_accession=G1, oxygen="Facultative aerobic", oxygen_status="verified",
         oxygen_source_tier="2 (genome BioSample)", optimal_temperature="55-60", optimal_temperature_status="secondary source only",
         optimal_pH="-2-4", optimal_pH_kind="optimum", habitat="<I>Hot</I>  spring", habitat_status="verified",
         electron_donors="H2 hydrogen; S2O3(2-) thiosulfate; Fe(II)",
         electron_donors_source="H2 hydrogen: literature (2 papers; 1 review); S2O3(2-) thiosulfate: survey sheet, audited OK; "
                                "Fe(II): survey sheet, audited ERROR (wrong)",
         electron_acceptors="O2 oxygen; H2 hydrogen", electron_acceptors_source="O2 oxygen: literature",
         electron_donor_acceptor_conflicting_reports="S2O3(2-) thiosulfate",
         reported_not_used_as_acceptor="NO3- nitrate [1 paper; putative only]"),
    dict(organism_id="org2", organism="Alpha two", genome_accession=G1, oxygen="facultative anaerobic", oxygen_source_tier="1 (strain)",
         optimal_temperature="65", optimal_temperature_status="verified", optimal_temperature_source_tier="1 (strain)",
         habitat="brand new place", electron_donors="CH4 methane", electron_acceptors="O2 oxygen"),
    dict(organism_id="org3", organism="Beta \u03bc", genome_accession=G2, oxygen="anaerobic", optimal_pH="6.7-6.9"),
    dict(organism_id="org4", organism="No genome", genome_accession="", oxygen="aerobic"),
    dict(organism_id="org5", organism="Shifted", genome_accession=G3, oxygen="aerobic", oxygen_status="1 (strain record)"),
)
COMPOUNDS = (
    "role\ttoken\tcategory\tcompound\n"
    "donor\tH2 hydrogen\tHydrogen oxidizing bacteria\thydrogen\n"
    "donor\tS2O3(2-) thiosulfate\tSulfur oxidizing bacteria\tthiosulfate\n"
    "donor\tFe(II)\tIron oxidizing bacteria\tiron(ii)\n"
    "acceptor\tO2 oxygen\tAerobes / aerobic respiration\toxygen\n"
    "acceptor\tNO3- nitrate\tNitrate-reducing bacteria (denitrifiers / DNRA)\tnitrate\n"
)
CORRECTIONS = f"# comment\ngenome_accession\torganism\ttoken\tfrom_role\tto_role\n{G1}\tAlpha one\tH2 hydrogen\tacceptor\tdonor\n"
HABITATS = "habitat\tclass\tclass_label\tn_organisms\tbasis\nHot spring\thotspring\t\t1\tsingle\n"


def write_inputs(directory: Path, encoding: str = "utf-8") -> Path:
    tables = directory / "config"
    tables.mkdir()
    (tables / "electron_compound_categories.tsv").write_text(COMPOUNDS, encoding="utf-8")
    (tables / "phenotype_corrections.tsv").write_text(CORRECTIONS, encoding="utf-8")
    (tables / "habitat_classes.tsv").write_text(HABITATS, encoding="utf-8")
    lines = ["\t".join(COLUMNS)] + ["\t".join(row.get(column, "") for column in COLUMNS) for row in ROWS]
    path = directory / "phenotype_data.tsv"
    path.write_bytes(("\n".join(lines) + "\n").encode(encoding))
    return path


class ParsingTest(unittest.TestCase):
    def test_numbers(self):
        self.assertEqual(parse_number("25"), (25.0, 25.0))
        self.assertEqual(parse_number("7.5-8"), (7.5, 8.0))
        self.assertEqual(parse_number("-2-30"), (-2.0, 30.0))
        self.assertIsNone(parse_number("warm"))
        self.assertIsNone(parse_number(""))

    def test_token_sources_keep_semicolons_inside_a_source(self):
        text = "H2 hydrogen: literature (2 papers; 1 review); S0 elemental sulfur: survey sheet"
        self.assertEqual(token_sources(["H2 hydrogen", "S0 elemental sulfur"], text), {
            "H2 hydrogen": "literature (2 papers; 1 review)", "S0 elemental sulfur": "survey sheet"})


class ReadPhenotypeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def read(self, encoding="utf-8"):
        path = write_inputs(self.dir, encoding)
        with self.assertLogs("ko_detector.phenotype", "INFO") as logs:
            phenotypes, electron = read_phenotype(path, read_tables(self.dir / "config"))
        return phenotypes, electron, "\n".join(logs.output)

    def test_traits_and_union_of_rows_sharing_a_genome(self):
        phenotypes, _, logs = self.read()
        self.assertEqual(sorted(phenotypes), [G1, G2])        # no accession / shifted rows are skipped
        self.assertIn("org5 Shifted", logs)
        self.assertIn("brand new place", logs)                 # unknown habitat is reported
        first = phenotypes[G1]
        self.assertEqual(first.organisms, ["Alpha one", "Alpha two"])
        oxygen = first.traits["oxygen"]
        self.assertEqual((oxygen.value, oxygen.klass), ("facultative anaerobic", "facultative anaerobic"))  # tier 1 wins
        self.assertEqual(oxygen.other_values, "Facultative aerobic [Alpha one]")
        temperature = first.traits["optimal_temperature"]
        self.assertEqual((temperature.value, temperature.number, temperature.status), ("65", 65.0, "verified"))
        self.assertEqual(temperature.other_values, "55-60 [Alpha one]")
        self.assertEqual(first.traits["optimal_pH"].number, 1.0)
        self.assertEqual(phenotypes[G2].traits["optimal_pH"].number, 6.8)     # not 6.800000000000001
        habitat = first.traits["habitat"]
        self.assertEqual((habitat.value, habitat.klass), ("Hot spring", "hotspring"))
        self.assertEqual(habitat.other_values, "brand new place [Alpha two]")
        self.assertEqual(list(phenotypes[G2].traits), ["oxygen", "optimal_pH"])

    def test_electron_records(self):
        _, electron, logs = self.read()
        self.assertEqual([(r.role, r.compound, r.category, r.consensus) for r in electron[G1]], [
            ("donor", "H2 hydrogen", "Hydrogen oxidizing bacteria", "used"),
            ("donor", "S2O3(2-) thiosulfate", "Sulfur oxidizing bacteria", "CONFLICTING"),
            ("donor", "Fe(II)", "Iron oxidizing bacteria", "audit_error"),
            ("acceptor", "O2 oxygen", "Aerobes / aerobic respiration", "used"),   # H2 moved to donor: already there
            ("acceptor", "NO3- nitrate", "Nitrate-reducing bacteria (denitrifiers / DNRA)", "not_used"),
            ("donor", "CH4 methane", "Other", "used"),                              # second row, unknown token
        ])
        self.assertEqual(electron[G1][0].source, "literature (2 papers; 1 review)")
        self.assertEqual(electron[G1][4].source, "reported not used (1 paper; putative only)")
        self.assertEqual(electron[G2], [])
        self.assertIn("donor:CH4 methane", logs)

    def test_reads_cp932(self):
        phenotypes, _, logs = self.read("cp932")
        self.assertEqual(sorted(phenotypes), [G1, G2])
        self.assertIn("read as cp932", logs)

    def test_missing_columns(self):
        write_inputs(self.dir)
        path = self.dir / "bad.tsv"
        path.write_text("organism_id\toxygen\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            read_phenotype(path, read_tables(self.dir / "config"))

    def test_tsv_round_trip_and_payload(self):
        phenotypes, electron, _ = self.read()
        results = [GenomeResult(g, STATUS_OK, hits=[Hit(g, "K17218", "WP_1", "1", "2", "1e-5", True)]) for g in (G1, G2, G3)]
        write_phenotypes(self.dir / "p.tsv", results, phenotypes)
        write_electron_categories(self.dir / "e.tsv", results, electron)
        self.assertEqual(read_phenotypes(self.dir / "p.tsv"), phenotypes)
        self.assertEqual(read_electron_categories(self.dir / "e.tsv"), electron)

        entries = [KoEntry(ko="K17218", functional_group="SOB", substrate="HS-", gene_enzyme="sqr", role="donor")]
        payload = build_payload(entries, results, "ja", electron=electron, phenotypes=phenotypes)
        self.assertEqual(payload["phenotype"]["traits"], ["oxygen", "optimal_temperature", "optimal_pH", "habitat"])
        genomes = {g["id"]: g for g in payload["genomes"]}
        self.assertEqual(genomes[G1]["phOrg"], "Alpha one; Alpha two")
        self.assertEqual(genomes[G1]["ph"][1][:4], ["65", "", 65.0, "verified"])
        self.assertEqual((genomes[G2]["ph"][1], genomes[G2]["ph"][3]), (None, None))
        self.assertEqual(genomes[G1]["el"][2][3:], ["audit_error", "survey sheet, audited ERROR (wrong)"])
        self.assertNotIn("ph", genomes[G3])


class ScanWithPhenotypeTest(unittest.TestCase):
    def test_scan_and_render_use_the_phenotype_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = write_inputs(root)
            reference = root / "ref.tsv"
            reference.write_text("Functional group\tElectron donor\tGene / enzyme\tKO\nSOB\tHS-\tsqr\tK17218\n", encoding="utf-8")
            for genome in (G1, G3):
                (root / "in" / genome).mkdir(parents=True)
                (root / "in" / genome / f"{genome}.kolist_gene.tsv").write_text(
                    "KO\tgene name\tthrshld\tscore\tE-value\tsignificant\nK17218\tWP_1\t1\t2\t1e-5\t*\n", encoding="utf-8")
            out = root / "out"
            with contextlib.redirect_stdout(io.StringIO()), self.assertLogs("ko_detector", "WARNING"):
                self.assertEqual(main(["scan", "-i", str(root / "in"), "-o", str(out), "-r", str(reference), "--no-taxonomy",
                                       "--phenotype", str(data), "--phenotype-tables", str(root / "config")]), 0)
            self.assertEqual(sorted(read_phenotypes(out / "genome_phenotype.tsv")), [G1])

            def report(path):
                page = path.read_text(encoding="utf-8")
                payload = json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>', page, re.S).group(1))
                return payload, {g["id"]: g for g in payload["genomes"]}

            scanned, genomes = report(out / "report.html")
            self.assertEqual(genomes[G1]["ph"][3][1], "hotspring")
            self.assertNotIn("ph", genomes[G3])

            # render-html reuses the result files (phenotype_data.tsv itself is gone)
            data.unlink()
            rendered = root / "rendered"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["render-html", "-i", str(out), "-o", str(rendered), "--no-taxonomy"]), 0)
            payload, _ = report(rendered / "report_en.html")
            self.assertEqual(payload["genomes"], report(out / "report_en.html")[0]["genomes"])
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["render-html", "-i", str(out), "-o", str(rendered), "--no-taxonomy", "--no-phenotype"]), 0)
            payload, _ = report(rendered / "report.html")
            self.assertIsNone(payload["phenotype"])
            self.assertIsNone(payload["electron"])
