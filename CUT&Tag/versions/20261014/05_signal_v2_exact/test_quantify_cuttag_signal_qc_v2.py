"""Offline regression tests; no real server data or optional Python packages required."""
import csv
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest
import sys

SCRIPT = Path(__file__).with_name("quantify_cuttag_signal_qc_v2.py")
spec = importlib.util.spec_from_file_location("signal_qc", SCRIPT)
qc = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = qc
spec.loader.exec_module(qc)


class SummaryTests(unittest.TestCase):
    def test_zero_filling_not_conditional_mean(self):
        result = qc.integrate([(20, 30, 10.0)], 0, 100, "CPM")
        self.assertEqual(result["mean_full_interval"], 1.0)
        self.assertEqual(result["mean_encoded_bases"], 10.0)
        self.assertEqual(result["encoded_fraction"], .1)

    def test_clipping_and_weighting(self):
        result = qc.integrate([(0, 10, 2.), (10, 30, 4.)], 5, 20, "CPM")
        self.assertAlmostEqual(result["mean_full_interval"], 50/15)
        self.assertEqual(result["encoded_fraction"], 1.)

    def test_no_data_has_zero_observed_coverage_only_for_CPM(self):
        self.assertEqual(qc.integrate(None, 0, 10, "CPM")["mean_full_interval"], 0.)
        self.assertIsNone(qc.integrate(None, 0, 10, "PPOIS")["mean_full_interval"])

    def test_p_value_gaps_are_missing(self):
        result = qc.integrate([(0, 5, 3.)], 0, 10, "PPOIS")
        self.assertIsNone(result["mean_full_interval"])
        self.assertEqual(result["mean_encoded_bases"], 3.)
        self.assertEqual(result["fraction_bp_score_ge_2"], .5)

    def test_nan_is_flagged_not_zero_filled(self):
        result = qc.integrate([(0, 5, math.nan)], 0, 10, "CPM")
        self.assertIsNone(result["mean_full_interval"])
        self.assertEqual(result["status"], "NONFINITE_ENTRIES")

    def test_overlaps_and_negative_values_fail(self):
        with self.assertRaises(ValueError):
            qc.integrate([(0, 10, 1.), (5, 20, 1.)], 0, 30, "CPM")
        with self.assertRaises(ValueError):
            qc.integrate([(0, 10, -1.)], 0, 30, "PPOIS")

    def test_strand_oriented_TSS_bins(self):
        sizes = {"1": 100000}
        gene_p = qc.Region("1", 1000, 2000, "G", "id", "+")
        gene_m = qc.Region("1", 1000, 2000, "G", "id", "-")
        self.assertEqual(qc.tss_bin(gene_p, -100, 100, sizes).key, ("1", 900, 1000))
        self.assertEqual(qc.tss_bin(gene_m, -100, 100, sizes).key, ("1", 2000, 2100))
        self.assertEqual(qc.tss_bin(gene_m, 0, 100, sizes).key, ("1", 1900, 2000))

    def test_positive_overlap_not_adjacency(self):
        index = qc.overlap_index([qc.Region("1", 10, 20), qc.Region("1", 15, 25), qc.Region("1", 40, 50)])
        self.assertFalse(qc.overlaps(qc.Region("1", 25, 40), index))
        self.assertTrue(qc.overlaps(qc.Region("1", 24, 26), index))
        self.assertFalse(qc.overlaps(qc.Region("chr1", 10, 20), index))

    def test_zero_denominator_ratio_and_pseudocounts(self):
        means = {f"{state}_{mark}_rep{rep}": 0. for state in (119, 137) for mark in ("Tag", "Pol2") for rep in (1, 2)}
        means.update({"119_IgG": 0., "137_IgG": 0.})
        means["119_Tag_rep1"] = 1.
        result = qc.comparisons(means)
        self.assertIsNone(result["119_Tag_rep1_ratio_no_pseudocount"])
        self.assertEqual(result["119_Tag_both_above_IgG_descriptive"], 0)
        self.assertAlmostEqual(result["119_Tag_rep1_log2ratio_pc0.01"], math.log2(101))

    def test_missing_contig_not_silently_zero(self):
        reader = qc.Reader.__new__(qc.Reader)
        reader.chroms = {"label": {"1": 100}}
        result = reader.summary("label", qc.Region("chr1", 0, 10))
        self.assertEqual(result["status"], "MISSING_CONTIG")
        self.assertIsNone(result["mean_full_interval"])

    def test_existing_annotation_and_seven_sites(self):
        screen = SCRIPT.parents[1] / "work/hypothesis_screen_import/hypothesis_screen_current_peaks"
        genes = qc.read_bed(screen / "annotations/candidate_genes.bed")
        self.assertEqual(len(genes), 12)
        polr1d = next(g for g in genes if g.name == "POLR1D")
        self.assertEqual(qc.tss_boundary(polr1d), 5044328)
        tag = qc.read_bed(screen / "peak_sets/macs2_119_Tag_rep1-rep2_shared.bed")
        pol2 = qc.overlap_index(qc.read_bed(screen / "peak_sets/macs2_119_Pol2_rep1-rep2_shared.bed"))
        self.assertEqual(sum(qc.overlaps(r, pol2) for r in tag), 7)

    def test_complete_output_column_alignment(self):
        class FakeReader:
            handles = {f"{state}_{mark}_rep{rep}": None for state in (119, 137) for mark in ("Tag", "Pol2") for rep in (1, 2)}
            handles.update({"119_IgG": None, "137_IgG": None, "ppois": None})
            kinds = {label: "PPOIS" if label == "ppois" else "CPM" for label in handles}

            def summary(self, label, region):
                return qc.integrate([(region.start, region.end, 2.)], region.start, region.end, self.kinds[label])

        with tempfile.TemporaryDirectory(dir=SCRIPT.parents[1] / "work") as folder:
            self.assertTrue(Path(folder).resolve().is_relative_to((SCRIPT.parents[1] / "work").resolve()))
            output = Path(folder) / "signals.tsv"
            qc.write_group(output, [qc.Region("1", 0, 10, "G")], FakeReader())
            with output.open(encoding="utf8") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(len(rows[0]), len(rows[1]))
            row = dict(zip(*rows))
            self.assertEqual(row["119_Tag_rep1.mean_full_interval"], "2")
            self.assertEqual(row["ppois.fraction_bp_score_ge_2"], "1")
            self.assertEqual(row["119_Tag_rep1_ratio_no_pseudocount"], "1")

    def test_matching_ignores_Tag_signal_and_preserves_features(self):
        class FakeReader:
            def summary(self, label, r):
                self_label = label
                if self_label not in ("119_IgG", "119_Pol2_rep1", "119_Pol2_rep2"):
                    raise AssertionError("Control selection must not use Tag signal")
                return {"mean_full_interval": 1.}

        class FakeFasta:
            def fetch(self, chrom, start, end):
                return ("ACGT" * (1+(end-start)//4))[:end-start]

        source = qc.Region("1", 1000, 1100)
        genes = [qc.Region("1", start, start+500, "G", str(start), "+") for start in (1000, 21000, 41000, 61000)]
        index = qc.overlap_index([qc.Region("1", 0, 90000)])
        controls, info = qc.matched_regions([source], genes, {"1": 100000}, FakeReader(), FakeFasta(), index, 2, 100, 4)
        self.assertEqual(len(controls), 2)
        self.assertEqual(info[0]["matched_controls"], 2)
        self.assertEqual(len({r.key for r in controls}), 2)
        self.assertTrue(all(r.end-r.start == 100 and r.chrom == "1" for r in controls))


if __name__ == "__main__":
    unittest.main(verbosity=2)
