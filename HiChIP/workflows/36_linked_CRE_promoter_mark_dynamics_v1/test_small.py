#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("pipeline", HERE / "run.py")
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)


def region(log2fc=0.0, support0=1, support6=1, h0=10, h6=10, a0=10, a6=10):
    return {
        "H3K27ac_log2FC": str(log2fc), "supported_0h": str(support0),
        "supported_6h": str(support6), "H3K27ac_0h_raw_tags": str(h0),
        "H3K27ac_6h_raw_tags": str(h6), "ATAC_0h_raw_tags": str(a0),
        "ATAC_6h_raw_tags": str(a6),
    }


class Tests(unittest.TestCase):
    def test_state_definition(self):
        self.assertEqual(P.derive_state(region(0.7), 1.5, 5), "strength")
        self.assertEqual(P.derive_state(region(-0.7), 1.5, 5), "weaken")
        self.assertEqual(P.derive_state(region(0.2), 1.5, 5), "stable")
        self.assertEqual(P.derive_state(region(0.7, support6=0), 1.5, 5), "stable")

    def test_overlap(self):
        self.assertEqual(P.overlap_bp(10, 20, 15, 30), 5)
        self.assertEqual(P.overlap_bp(10, 20, 20, 30), 0)

    def test_rank_and_correlations(self):
        self.assertEqual(P.ranks([1, 1, 3]), [1.5, 1.5, 3.0])
        self.assertAlmostEqual(P.pearson([1, 2, 3], [2, 4, 6]), 1.0)

    def test_signal_calculation(self):
        rows = [{"feature_id": "x"}]
        labels = [f"{assay}_{time}" for assay in (*P.ASSAYS, "IgG") for time in P.TIMES]
        counts = {"x": {label: 10 for label in labels}}
        totals = {label: 100 for label in labels}
        P.add_signals(rows, counts, totals, 10_000_000, 0.5)
        for assay in P.ASSAYS:
            self.assertAlmostEqual(rows[0][f"{assay}_log2FC"], 0.0)
        for assay in P.CHIP_ASSAYS:
            self.assertAlmostEqual(rows[0][f"{assay}_delta_log2_vs_IgG"], 0.0)


if __name__ == "__main__":
    unittest.main()
