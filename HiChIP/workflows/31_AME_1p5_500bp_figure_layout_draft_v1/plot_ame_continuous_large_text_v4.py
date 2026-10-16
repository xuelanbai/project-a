#!/usr/bin/env python3
"""Continuous significant AME list with large text for composite Figure 3."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from plot_ame_combined import COLORS, CONTRASTS, clean_label, read_rows, render_pdf, render_png, render_svg, select_rows
from plot_ame_unified_significant_v2 import write_plot_data


WIDTH = 900
HEIGHT = 610
PLOT_LEFT = 145
PLOT_RIGHT = 875
PLOT_WIDTH = PLOT_RIGHT - PLOT_LEFT
FIRST_Y = 94
ROW_STEP = 19
BAR_HEIGHT = 13.5
AXIS_Y = 500
AXIS_MAX = 2.0


def build_scene(selected: dict[str, list[dict[str, object]]]) -> list[dict[str, object]]:
    scene: list[dict[str, object]] = []

    def text_item(x: float, y: float, value: str, size: float, color: str,
                  anchor: str = "start", bold: bool = False) -> None:
        scene.append({"kind": "text", "x": x, "y": y, "text": value, "size": size,
                      "color": color, "anchor": anchor, "bold": bold})

    text_item(WIDTH / 2, 30, "AME motif enrichment", 18.0, COLORS["text"], "middle", True)

    for tick in (0, 0.5, 1.0, 1.5, 2.0):
        x = PLOT_LEFT + tick / AXIS_MAX * PLOT_WIDTH
        scene.append({"kind": "line", "x1": x, "y1": 61, "x2": x, "y2": AXIS_Y,
                      "color": COLORS["grid"], "width": 0.65})
        text_item(x, AXIS_Y + 19, f"{tick:g}", 9.4, COLORS["axis"], "middle")
    scene.append({"kind": "line", "x1": PLOT_LEFT, "y1": AXIS_Y, "x2": PLOT_RIGHT, "y2": AXIS_Y,
                  "color": COLORS["axis"], "width": 0.8})

    ordered: list[tuple[str, dict[str, object]]] = []
    for state in ("strength", "weaken"):
        ordered.extend((state, row) for row in selected[state])
    for index, (state, row) in enumerate(ordered):
        center_y = FIRST_Y + index * ROW_STEP
        value = float(row["log2_enrichment"])
        text_item(PLOT_LEFT - 10, center_y, clean_label(row), 10.2, COLORS["text"], "end")
        scene.append({"kind": "rect", "x": PLOT_LEFT, "y": center_y - BAR_HEIGHT / 2,
                      "width": value / AXIS_MAX * PLOT_WIDTH, "height": BAR_HEIGHT,
                      "color": COLORS[state]})

    # A clean legend in the unused upper-right area; state is encoded by colour only.
    scene.append({"kind": "rect", "x": 675, "y": 61, "width": 202, "height": 71, "color": "#FFFFFF"})
    for state, y, windows in (("strength", 82, 3760), ("weaken", 113, 1398)):
        label = "Strength" if state == "strength" else "Weaken"
        scene.append({"kind": "rect", "x": 690, "y": y - 8, "width": 15, "height": 15, "color": COLORS[state]})
        text_item(716, y, label, 11.8, COLORS[state], "start", True)
        text_item(808, y, f"n={windows:,} windows", 9.0, COLORS["axis"], "start")

    text_item((PLOT_LEFT + PLOT_RIGHT) / 2, AXIS_Y + 49,
              "log2(foreground / GC-matched control)", 12.4, COLORS["text"], "middle", True)
    text_item((PLOT_LEFT + PLOT_RIGHT) / 2, AXIS_Y + 87,
              "BH q <= 0.05; 1.5-fold H3K27ac; 500-bp ATAC-core windows",
              9.2, COLORS["axis"], "middle")
    return scene


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
    maximum = max(float(row["log2_enrichment"]) for values in selected.values() for row in values)
    if maximum > AXIS_MAX:
        raise ValueError(f"A bar exceeds the fixed shared x-axis maximum: {maximum}")

    args.output.mkdir(parents=True, exist_ok=True)
    stem = "AME_1p5_500bp_significant_enrichment_continuous_layout_v4"
    scene = build_scene(selected)
    render_png(args.output / f"{stem}.png", HEIGHT, scene)
    render_pdf(args.output / f"{stem}.pdf", HEIGHT, scene)
    render_svg(args.output / f"{stem}.svg", HEIGHT, scene)
    write_plot_data(args.output / "plot_data_significant_enrichment_rank.tsv", selected)
    parameters = {
        "title": "AME motif enrichment",
        "selection": "representative motif families with enrichment_ratio > 1 and within-contrast BH q <= 0.05",
        "ordering": "log2 enrichment descending within Strength, followed immediately by Weaken",
        "non_significant_displayed": False,
        "strength_motifs": len(selected["strength"]),
        "weaken_motifs": len(selected["weaken"]),
        "axis_max": AXIS_MAX,
        "axis_policy": "shared",
        "legend": "upper right",
        "font_scale": "large-text composite-figure version",
    }
    (args.output / "design_parameters.json").write_text(
        json.dumps(parameters, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(parameters, ensure_ascii=False))


if __name__ == "__main__":
    main()
