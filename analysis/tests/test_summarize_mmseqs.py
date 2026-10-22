import io
import tempfile
import unittest
from pathlib import Path

from analysis.summarize_mmseqs import normalize_identity, parse_rows, summarize, taxon_label


class SummarizeMmseqsTests(unittest.TestCase):
    def test_identity_fraction_and_percentage_are_normalized(self):
        self.assertAlmostEqual(normalize_identity(0.953), 95.3)
        self.assertAlmostEqual(normalize_identity(95.0), 95.0)

    def test_taxon_annotation_is_extracted_from_target_header(self):
        self.assertEqual(taxon_label("ref_OS=Example_insect_OX=123_GN=COI"), "Example insect")
        self.assertEqual(taxon_label("ref_without_organism"), "(no OS annotation)")

    def test_best_hit_is_selected_by_bitscore_and_query_coverage_is_reported(self):
        text = (
            "query_contig_length_188\ttarget_a_OS=Example_insect_OX=123\t0.953\t93\t1\t0\t94\t186\t50\t80\t1e-12\t64\n"
            "query_contig_length_188\ttarget_b_OS=Example_worm_OX=456\t0.900\t120\t12\t0\t1\t120\t10\t129\t1e-20\t78\n"
        )
        best, valid, skipped = parse_rows(io.StringIO(text))
        self.assertEqual(valid, 2)
        self.assertEqual(skipped, 0)
        row = best["query_contig_length_188"]
        self.assertEqual(row["target"], "target_b_OS=Example_worm_OX=456")
        self.assertAlmostEqual(row["query_coverage_percent"], round(120 / 188 * 100, 3))

    def test_malformed_rows_are_skipped(self):
        best, valid, skipped = parse_rows(io.StringIO("too\tfew\tcolumns\n"))
        self.assertEqual(best, {})
        self.assertEqual(valid, 0)
        self.assertEqual(skipped, 1)

    def test_summary_outputs_are_written(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "hits.tsv"
            source.write_text("query_length_100\tref_OS=Example_species_OX=1\t0.9\t50\t2\t0\t1\t50\t1\t50\t1e-9\t20\n")
            output = root / "out"
            summary = summarize(source, output, "F734")
            self.assertEqual(summary["queries_with_hits"], 1)
            self.assertTrue((output / "F734.per_query.tsv").exists())
            self.assertTrue((output / "F734.top_hit_annotations.tsv").exists())
            self.assertTrue((output / "F734.run_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
