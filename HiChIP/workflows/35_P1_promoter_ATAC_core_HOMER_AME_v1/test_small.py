#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("pipeline", HERE / "run.py")
PIPELINE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PIPELINE)


def promoter(region_id, start, end, log2fc, state, support0=1, support6=1,
             h0=100, h6=100, a0=10, a6=10):
    return {
        "region_id": region_id, "chrom": "1", "start": start, "end": end,
        "H3K27ac_log2FC": str(log2fc), "promoter_state": state,
        "supported_0h": str(support0), "supported_6h": str(support6),
        "H3K27ac_0h_raw_tags": str(h0), "H3K27ac_6h_raw_tags": str(h6),
        "ATAC_0h_raw_tags": str(a0), "ATAC_6h_raw_tags": str(a6),
    }


class Tests(unittest.TestCase):
    def test_state_rule(self):
        self.assertEqual(PIPELINE.derive_state(promoter("a", 10, 20, 0.7, "strength"), 1.5, 5), "strength")
        self.assertEqual(PIPELINE.derive_state(promoter("b", 30, 40, -0.7, "weaken"), 1.5, 5), "weaken")
        self.assertEqual(PIPELINE.derive_state(promoter("c", 50, 60, 0.2, "stable"), 1.5, 5), "stable")
        self.assertEqual(
            PIPELINE.derive_state(promoter("d", 70, 80, 0.7, "stable", support6=0), 1.5, 5), "stable"
        )

    def test_nonoverlap_index(self):
        rows = [promoter("a", 10, 20, 0, "stable"), promoter("b", 30, 40, 0, "stable")]
        index = PIPELINE.NonOverlappingRegionIndex(rows)
        self.assertEqual(index.contains("chr1", 10)["region_id"], "a")
        self.assertIsNone(index.contains("1", 20))
        self.assertEqual(index.contains("1", 39)["region_id"], "b")
        with self.assertRaises(ValueError):
            PIPELINE.NonOverlappingRegionIndex([
                promoter("a", 10, 20, 0, "stable"), promoter("x", 19, 25, 0, "stable")
            ])

    def test_gc_matching_without_replacement(self):
        foreground = [
            {"core_id": "f1", "_gc": 0.40}, {"core_id": "f2", "_gc": 0.61},
        ]
        controls = [
            {"core_id": "c1", "_gc": 0.41}, {"core_id": "c2", "_gc": 0.60},
        ]
        pairs, unmatched = PIPELINE.match_gc(foreground, controls, 0.05, 7)
        self.assertEqual(len(pairs), 2)
        self.assertFalse(unmatched)
        self.assertEqual(len({pair[1]["core_id"] for pair in pairs}), 2)

    def test_tsv_header_for_empty_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.tsv"
            PIPELINE.write_tsv(path, [], ["a", "b"])
            self.assertEqual(path.read_text(encoding="utf-8"), "a\tb\n")


if __name__ == "__main__":
    unittest.main()
