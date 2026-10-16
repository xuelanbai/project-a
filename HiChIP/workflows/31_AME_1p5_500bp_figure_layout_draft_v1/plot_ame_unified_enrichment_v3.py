#!/usr/bin/env python3
"""Unified significant AME figure, ranked by enrichment within each state."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from plot_ame_combined import CONTRASTS, read_rows, render_pdf, render_png, render_svg, select_rows
from plot_ame_unified_significant_v2 import build_scene, write_plot_data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    rows = read_rows(args.input)
    selected = {
        state: select_rows(rows, state, "significant_enrichment_rank", args.top_n)
        for state in CONTRASTS
    }
    if not all(selected.values()):
        raise ValueError("One state has no significant enriched representative motifs")

    height, scene, axis_max = build_scene(selected)
    args.output.mkdir(parents=True, exist_ok=True)
    stem = "AME_1p5_500bp_significant_enrichment_rank_unified_layout_v3"
    render_png(args.output / f"{stem}.png", height, scene)
    render_pdf(args.output / f"{stem}.pdf", height, scene)
    render_svg(args.output / f"{stem}.svg", height, scene)
    write_plot_data(args.output / "plot_data_significant_enrichment_rank.tsv", selected)

    parameters = {
        "title": "AME motif enrichment",
        "selection": "representative motif families with enrichment_ratio > 1 and within-contrast BH q <= 0.05",
        "ordering": "log2 enrichment descending within Strength and Weaken; BH q ascending for ties",
        "top_n_per_state": args.top_n,
        "strength_motifs": len(selected["strength"]),
        "weaken_motifs": len(selected["weaken"]),
        "axis_max": axis_max,
        "axis_policy": "shared across both sections",
    }
    (args.output / "design_parameters.json").write_text(
        json.dumps(parameters, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(parameters, ensure_ascii=False))


if __name__ == "__main__":
    main()
