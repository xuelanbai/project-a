#!/usr/bin/env python3
"""Put significant Strength and Weaken AME motifs in one compact panel."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from plot_ame_combined import (
    COLORS,
    CONTRASTS,
    clean_label,
    nice_axis,
    read_rows,
    render_pdf,
    render_png,
    render_svg,
    select_rows,
)


WIDTH = 720
ROW_STEP = 18
BAR_HEIGHT = 12
PLOT_LEFT = 190
PLOT_RIGHT = 690
PLOT_WIDTH = PLOT_RIGHT - PLOT_LEFT


def significant_q_rows(rows: list[dict[str, str]], state: str, top_n: int) -> list[dict[str, object]]:
    ordered = select_rows(rows, state, "q_rank", 100000)
    return [row for row in ordered if bool(row["significant"])][:top_n]


def build_scene(selected: dict[str, list[dict[str, object]]]) -> tuple[int, list[dict[str, object]], float]:
    strength_count = len(selected["strength"])
    weaken_count = len(selected["weaken"])
    title_y = 24
    strength_header_y = 57
    strength_first_y = 82
    strength_last_y = strength_first_y + (strength_count - 1) * ROW_STEP
    weaken_header_y = strength_last_y + 32
    weaken_first_y = weaken_header_y + 25
    weaken_last_y = weaken_first_y + (weaken_count - 1) * ROW_STEP
    axis_y = weaken_last_y + 22
    height = axis_y + 58

    raw_max = max(float(row["log2_enrichment"]) for state_rows in selected.values() for row in state_rows)
    axis_max, ticks = nice_axis(raw_max * 1.04)
    scene: list[dict[str, object]] = []

    def text_item(x: float, y: float, value: str, size: float, color: str, anchor: str = "start", bold: bool = False) -> None:
        scene.append({"kind": "text", "x": x, "y": y, "text": value, "size": size,
                      "color": color, "anchor": anchor, "bold": bold})

    text_item(WIDTH / 2, title_y, "AME motif enrichment", 12.0, COLORS["text"], "middle", True)
    grid_top = strength_first_y - ROW_STEP * 0.55
    grid_bottom = axis_y
    for tick in ticks:
        x = PLOT_LEFT + tick / axis_max * PLOT_WIDTH
        scene.append({"kind": "line", "x1": x, "y1": grid_top, "x2": x, "y2": grid_bottom,
                      "color": COLORS["grid"], "width": 0.55})
        text_item(x, axis_y + 17, f"{tick:g}", 7.3, COLORS["axis"], "middle")
    scene.append({"kind": "line", "x1": PLOT_LEFT, "y1": axis_y, "x2": PLOT_RIGHT, "y2": axis_y,
                  "color": COLORS["axis"], "width": 0.7})

    sections = {
        "strength": (strength_header_y, strength_first_y, 3760),
        "weaken": (weaken_header_y, weaken_first_y, 1398),
    }
    for state in ("strength", "weaken"):
        header_y, first_y, windows = sections[state]
        title = "Strength" if state == "strength" else "Weaken"
        color = COLORS[state]
        scene.append({"kind": "rect", "x": 24, "y": header_y - 7, "width": 3.5, "height": 14, "color": color})
        text_item(34, header_y, title, 10.0, color, "start", True)
        text_item(91, header_y, f"n={windows:,} windows", 7.4, COLORS["axis"])
        for index, row in enumerate(selected[state]):
            center_y = first_y + index * ROW_STEP
            value = float(row["log2_enrichment"])
            text_item(PLOT_LEFT - 8, center_y, clean_label(row), 7.8, COLORS["text"], "end")
            scene.append({"kind": "rect", "x": PLOT_LEFT, "y": center_y - BAR_HEIGHT / 2,
                          "width": value / axis_max * PLOT_WIDTH, "height": BAR_HEIGHT, "color": color})

    text_item((PLOT_LEFT + PLOT_RIGHT) / 2, axis_y + 35,
              "log2(foreground / GC-matched control)", 8.2, COLORS["text"], "middle")
    text_item((PLOT_LEFT + PLOT_RIGHT) / 2, axis_y + 52,
              "BH q <= 0.05; 1.5-fold H3K27ac; 500-bp ATAC-core windows",
              7.0, COLORS["axis"], "middle")
    return height, scene, axis_max


def write_plot_data(path: Path, selected: dict[str, list[dict[str, object]]]) -> None:
    fields = [
        "section_order", "display_order", "state", "contrast", "motif_family",
        "motif_id", "motif_alt_id", "enrichment_ratio", "log2_enrichment",
        "AME_BH_q_within_contrast", "AME_BH_q_all_contrasts",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for section_order, state in enumerate(("strength", "weaken"), 1):
            for display_order, row in enumerate(selected[state], 1):
                writer.writerow({
                    "section_order": section_order,
                    "display_order": display_order,
                    "state": state,
                    "contrast": row["contrast"],
                    "motif_family": row["motif_family"],
                    "motif_id": row["motif_id"],
                    "motif_alt_id": row["motif_alt_id"],
                    "enrichment_ratio": row["enrichment_ratio"],
                    "log2_enrichment": row["log2_enrichment"],
                    "AME_BH_q_within_contrast": row["AME_BH_q_within_contrast"],
                    "AME_BH_q_all_contrasts": row["AME_BH_q_all_contrasts"],
                })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()
    rows = read_rows(args.input)
    selected = {state: significant_q_rows(rows, state, args.top_n) for state in CONTRASTS}
    if not all(selected.values()):
        raise ValueError("One state has no significant enriched representative motifs")

    height, scene, axis_max = build_scene(selected)
    args.output.mkdir(parents=True, exist_ok=True)
    stem = "AME_1p5_500bp_significant_q_rank_unified_layout_v2"
    render_png(args.output / f"{stem}.png", height, scene)
    render_pdf(args.output / f"{stem}.pdf", height, scene)
    render_svg(args.output / f"{stem}.svg", height, scene)
    write_plot_data(args.output / "plot_data_significant_q_rank.tsv", selected)
    parameters = {
        "title": "AME motif enrichment",
        "selection": "representative motif families with enrichment_ratio > 1 and within-contrast BH q <= 0.05",
        "ordering": "within-contrast BH q ascending, enrichment ratio descending for ties",
        "top_n_per_state": args.top_n,
        "strength_motifs": len(selected["strength"]),
        "weaken_motifs": len(selected["weaken"]),
        "axis_max": axis_max,
        "axis_policy": "shared across both sections",
        "strength_color": COLORS["strength"],
        "weaken_color": COLORS["weaken"],
    }
    (args.output / "design_parameters.json").write_text(
        json.dumps(parameters, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(parameters, ensure_ascii=False))


if __name__ == "__main__":
    main()
