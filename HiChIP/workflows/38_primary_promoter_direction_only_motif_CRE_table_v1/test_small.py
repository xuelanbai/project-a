#!/usr/bin/env python3
"""Small deterministic tests that do not require cluster inputs."""
import math
import unittest

import run


def feature(log2fc=0.0, support0=1, support6=1, h0=10, h6=10, a0=10, a6=10):
    return {
        "H3K27ac_log2FC": log2fc,
        "supported_0h": support0, "supported_6h": support6,
        "H3K27ac_0h_raw_tags": h0, "H3K27ac_6h_raw_tags": h6,
        "ATAC_0h_raw_tags": a0, "ATAC_6h_raw_tags": a6,
    }


class StateTests(unittest.TestCase):
    def test_direction_only_split(self):
        self.assertEqual(run.state_and_reason(feature(1e-6), 1.0, 5)[0], "strength")
        self.assertEqual(run.state_and_reason(feature(-1e-6), 1.0, 5)[0], "weaken")
        self.assertEqual(run.state_and_reason(feature(0.0), 1.0, 5)[0], "stable")

    def test_active_state_support_and_minimum_are_required(self):
        cutoff = 0.01
        self.assertEqual(run.state_and_reason(feature(cutoff, support6=0), 1.0, 5)[0], "stable")
        self.assertEqual(run.state_and_reason(feature(cutoff, a6=4), 1.0, 5)[0], "stable")
        self.assertEqual(run.state_and_reason(feature(-cutoff, h0=4), 1.0, 5)[0], "stable")

    def test_boundary_is_inclusive(self):
        row = feature(math.log2(1.2), h6=5, a6=5)
        self.assertEqual(run.state_and_reason(row, 1.2, 5)[0], "strength")


class PrimaryPromoterTests(unittest.TestCase):
    def test_one_promoter_per_gene_prefers_TSS_containment(self):
        base = {
            "H3K4me3_0h_tags_per_10M": 10, "H3K4me3_6h_tags_per_10M": 10,
            "ATAC_0h_tags_per_10M": 10, "ATAC_6h_tags_per_10M": 10,
        }
        promoters = {
            "P1": dict(base, region_id="P1", start=90, end=110),
            "P2": dict(base, region_id="P2", start=120, end=140),
        }
        associations = [
            {"gene_id": "G1", "region_id": "P2", "chrom": "1", "gene_tss": 100,
             "gene_strand": "+"},
            {"gene_id": "G1", "region_id": "P1", "chrom": "1", "gene_tss": 100,
             "gene_strand": "+"},
        ]
        selected, candidates, qc = run.select_primary_promoters(associations, promoters)
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["region_id"], "P1")
        self.assertEqual(sum(row["is_primary_promoter_for_gene"] for row in candidates), 1)
        self.assertEqual(qc["genes_with_multiple_strict_candidate_promoters"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
