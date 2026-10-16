#!/usr/bin/env python3
"""Focused unit tests for overlap and consensus behavior."""
from export_and_count_definitions import (
    MergedIntervalIndex,
    consensus_components,
    intersect_interval_sets,
    merge_true_overlaps,
)


def main() -> None:
    first = MergedIntervalIndex([("1", 10, 20), ("1", 30, 40)])
    second = MergedIntervalIndex([("1", 15, 35)])
    assert first.overlaps("chr1", 19, 21)
    assert not first.overlaps("1", 20, 30)
    assert intersect_interval_sets(first, second) == [("1", 15, 20), ("1", 30, 35)]
    assert merge_true_overlaps([("1", 10, 20), ("1", 20, 30)]) == [("1", 10, 20), ("1", 20, 30)]
    components = consensus_components({"0h": [("1", 10, 20)], "6h": [("1", 15, 25), ("1", 25, 30)]})
    assert [(row["start"], row["end"], row["state"]) for row in components] == [
        (10, 25, "both"),
        (25, 30, "6h_only"),
    ]
    print("SMALL_TEST_OK")


if __name__ == "__main__":
    main()
