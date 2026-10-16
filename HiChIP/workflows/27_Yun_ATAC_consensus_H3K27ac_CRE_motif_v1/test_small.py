#!/usr/bin/env python3
import unittest
from pathlib import Path
import tempfile

import run


class DesignTests(unittest.TestCase):
    def test_h3_consensus_merges_overlap_not_bookends(self):
        peaks = [
            dict(chrom="1", start=10, end=20, source_time="0h"),
            dict(chrom="1", start=15, end=25, source_time="6h"),
            dict(chrom="1", start=25, end=30, source_time="6h"),
        ]
        regions, membership = run.consensus_h3(peaks)
        self.assertEqual([(x["start"], x["end"]) for x in regions], [(10, 25), (25, 30)])
        self.assertEqual(regions[0]["supported_0h"], 1)
        self.assertEqual(regions[0]["supported_6h"], 1)
        self.assertEqual(len(membership), 3)

    def test_atac_consensus_merges_true_overlap(self):
        peaks = [
            dict(chrom="1", start=900, end=1100, summit=1000,
                 source_time="0h", signal_value=5),
            dict(chrom="1", start=1000, end=1200, summit=1080,
                 source_time="6h", signal_value=7),
            dict(chrom="1", start=1400, end=1600, summit=1500,
                 source_time="6h", signal_value=4),
        ]
        cores, membership = run.consensus_atac(peaks)
        self.assertEqual(len(cores), 2)
        self.assertEqual((cores[0]["start"], cores[0]["end"]), (900, 1200))
        self.assertEqual(cores[0]["member_peak_count"], 2)
        self.assertEqual(cores[0]["supported_0h"], 1)
        self.assertEqual(cores[0]["supported_6h"], 1)
        self.assertEqual(len(membership), 3)

    def test_atac_consensus_does_not_merge_bookends(self):
        peaks = [
            dict(chrom="1", start=900, end=1100, summit=1000,
                 source_time="0h", signal_value=1),
            dict(chrom="1", start=1100, end=1300, summit=1200,
                 source_time="6h", signal_value=1),
        ]
        cores, _ = run.consensus_atac(peaks)
        self.assertEqual(len(cores), 2)

    def test_strength_ignores_atac_direction(self):
        region = {"supported_6h": 1, "H3K27ac_6h_raw_tags": 10,
                  "ATAC_6h_raw_tags": 10, "H3K27ac_log2FC": 0.5,
                  "ATAC_log2FC": -2.0}
        self.assertTrue(run.selected(region, "strength", 1.2, 5))

    def test_representative_family_collapses_ctcf(self):
        base = dict(contrast="x", AME_threshold_adjusted_p=0.01,
                    AME_BH_q_within_contrast=0.02, enrichment_ratio=2.0)
        rows = [dict(base, motif_id="m1", motif_alt_id="CTCF"),
                dict(base, motif_id="m2", motif_alt_id="CTCFL", enrichment_ratio=3.0)]
        selected = run.representative_ame(rows)
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["motif_family"], "CTCF")
        self.assertEqual(selected[0]["motif_id"], "m2")

    def test_multicov_parser_accepts_bed6(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.tsv"
            path.write_text("chr1\t0\t100\tID1\t0\t.\t7\t11\n")
            counts = run.parse_multicov(path, ["sample_1", "sample_2"])
        self.assertEqual(counts["ID1"], {"sample_1": 7, "sample_2": 11})


if __name__ == "__main__":
    unittest.main()
