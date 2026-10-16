#!/usr/bin/env python3
"""Dependency-free tests for P1 construction, consensus and state rules."""

from run import GenePromoter, IntervalIndex, build_consensus, build_p1, classify


def boolean_index(rows):
    index = IntervalIndex()
    for chrom, start, end in rows:
        index.add(chrom, start, end, True)
    return index


def main():
    genes = IntervalIndex()
    gene = GenePromoter("1", 800, 1201, "GENE1", "G1", 1000, "+")
    genes.add("1", gene.window_start, gene.window_end, gene)
    peaks_0h, associations_0h, qc_0h = build_p1(
        "0h", [("1", 900, 1100), ("1", 3000, 3200)],
        boolean_index([("1", 950, 980)]),
        boolean_index([("1", 1000, 1050)]),
        genes,
    )
    peaks_6h, associations_6h, _ = build_p1(
        "6h", [("1", 950, 1150)],
        boolean_index([("1", 960, 990)]),
        boolean_index([("1", 1020, 1060)]),
        genes,
    )
    assert len(peaks_0h) == 1 and len(associations_0h) == 1
    assert qc_0h["excluded_no_TSS_window_overlap"] == 1
    assert len(peaks_6h) == 1 and len(associations_6h) == 1
    regions, membership = build_consensus({"0h": peaks_0h, "6h": peaks_6h})
    assert len(regions) == 1 and len(membership) == 2
    assert regions[0]["start"] == 900 and regions[0]["end"] == 1150
    assert regions[0]["supported_0h"] == 1 and regions[0]["supported_6h"] == 1

    strength = dict(H3K27ac_log2FC=1.0, supported_0h=1, supported_6h=1,
                    H3K27ac_0h_raw_tags=10, H3K27ac_6h_raw_tags=12)
    weaken = dict(H3K27ac_log2FC=-1.0, supported_0h=1, supported_6h=1,
                  H3K27ac_0h_raw_tags=12, H3K27ac_6h_raw_tags=10)
    assert classify(strength, 1.5, 5)[0] == "strength"
    assert classify(weaken, 1.5, 5)[0] == "weaken"
    strength["H3K27ac_6h_raw_tags"] = 4
    assert classify(strength, 1.5, 5)[0] == "stable"
    print("SMALL_TEST_OK")


if __name__ == "__main__":
    main()
