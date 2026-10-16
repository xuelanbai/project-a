#!/usr/bin/env python3
"""Draft three large-text AME layouts with separate Strength/Weaken x-axes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from plot_ame_combined import (
    COLORS,
    CONTRASTS,
    FONT_BOLD,
    clean_label,
    read_rows,
    render_pdf,
    render_png,
    render_svg,
    select_rows,
)
from plot_ame_unified_significant_v2 import write_plot_data


WIDTH = 960
PLOT_LEFT = 190
PLOT_RIGHT = 920
PLOT_WIDTH = PLOT_RIGHT - PLOT_LEFT
WINDOW_COUNTS = {"strength": 3760, "weaken": 1398}
AXES = {
    "strength": (0.6, [0.0, 0.2, 0.4, 0.6]),
    "weaken": (1.8, [0.0, 0.6, 1.2, 1.8]),
}


VARIANTS = {
    # Same row spacing in both panels: shortest and most compact.
    "v5a_compact": {"strength_step": 24.0, "weaken_step": 24.0, "panel_gap": 82.0},
    # Slightly more space for the five Weaken motifs: recommended compromise.
    "v5b_balanced": {"strength_step": 24.0, "weaken_step": 32.0, "panel_gap": 84.0},
    # Give Weaken a near-equal visual block for comparison.
    "v5c_equal_emphasis": {"strength_step": 24.0, "weaken_step": 58.0, "panel_gap": 86.0},
}


def build_scene(
    selected: dict[str, list[dict[str, object]]],
    variant: str,
) -> tuple[int, list[dict[str, object]], dict[str, object]]:
    settings = VARIANTS[variant]
    strength_step = float(settings["strength_step"])
    weaken_step = float(settings["weaken_step"])
    panel_gap = float(settings["panel_gap"])

    title_y = 28.0
    strength_header_y = 72.0
    strength_first_y = 112.0
    strength_axis_y = strength_first_y + (len(selected["strength"]) - 0.5) * strength_step + 24.0
    weaken_header_y = strength_axis_y + panel_gap
    weaken_first_y = weaken_header_y + 40.0
    weaken_axis_y = weaken_first_y + (len(selected["weaken"]) - 0.5) * weaken_step + 24.0
    common_label_y = weaken_axis_y + 53.0
    note_y = common_label_y + 35.0
    height = int(note_y + 26.0)
    scene: list[dict[str, object]] = []

    def text_item(
        x: float,
        y: float,
        value: str,
        size: float,
        color: str,
        anchor: str = "start",
        bold: bool = False,
    ) -> None:
        scene.append(
            {
                "kind": "text",
                "x": x,
                "y": y,
                "text": value,
                "size": size,
                "color": color,
                "anchor": anchor,
                "bold": bold,
            }
        )

    text_item(WIDTH / 2, title_y, "AME motif enrichment", 22.0, COLORS["text"], "middle", True)

    panel_geometry = {
        "strength": (strength_header_y, strength_first_y, strength_axis_y, strength_step),
        "weaken": (weaken_header_y, weaken_first_y, weaken_axis_y, weaken_step),
    }
    for state in ("strength", "weaken"):
        header_y, first_y, axis_y, row_step = panel_geometry[state]
        label = "Strength" if state == "strength" else "Weaken"
        text_item(PLOT_LEFT, header_y, label, 17.0, COLORS[state], "start", True)
        text_item(
            PLOT_RIGHT,
            header_y,
            f"n={WINDOW_COUNTS[state]:,} windows",
            12.5,
            COLORS["axis"],
            "end",
        )

        axis_max, ticks = AXES[state]
        plot_top = first_y - row_step * 0.55
        for tick in ticks:
            x = PLOT_LEFT + tick / axis_max * PLOT_WIDTH
            scene.append(
                {
                    "kind": "line",
                    "x1": x,
                    "y1": plot_top,
                    "x2": x,
                    "y2": axis_y,
                    "color": COLORS["grid"],
                    "width": 0.75,
                }
            )
            text_item(x, axis_y + 20.0, f"{tick:g}", 12.0, COLORS["axis"], "middle")
        scene.append(
            {
                "kind": "line",
                "x1": PLOT_LEFT,
                "y1": axis_y,
                "x2": PLOT_RIGHT,
                "y2": axis_y,
                "color": COLORS["axis"],
                "width": 0.9,
            }
        )

        bar_height = min(18.0, row_step * 0.70)
        for index, row in enumerate(selected[state]):
            center_y = first_y + index * row_step
            value = float(row["log2_enrichment"])
            text_item(PLOT_LEFT - 13.0, center_y, clean_label(row), 13.2, COLORS["text"], "end")
            scene.append(
                {
                    "kind": "rect",
                    "x": PLOT_LEFT,
                    "y": center_y - bar_height / 2,
                    "width": value / axis_max * PLOT_WIDTH,
                    "height": bar_height,
                    "color": COLORS[state],
                }
            )

    text_item(
        (PLOT_LEFT + PLOT_RIGHT) / 2,
        common_label_y,
        "log2(foreground / GC-matched control)",
        16.0,
        COLORS["text"],
        "middle",
        True,
    )
    text_item(
        (PLOT_LEFT + PLOT_RIGHT) / 2,
        note_y,
        "BH q <= 0.05; 1.5-fold H3K27ac; 500-bp ATAC-core windows",
        11.5,
        COLORS["axis"],
        "middle",
    )

    parameters = {
        "variant": variant,
        "selection": "representative motif families with enrichment_ratio > 1 and within-contrast BH q <= 0.05",
        "ordering": "log2 enrichment descending within each state",
        "non_significant_displayed": False,
        "axis_policy": "separate fixed clean axes",
        "axis_max": {state: AXES[state][0] for state in CONTRASTS},
        "strength_motifs": len(selected["strength"]),
        "weaken_motifs": len(selected["weaken"]),
        **settings,
    }
    return height, scene, parameters


def make_comparison(output: Path, generated: list[Path]) -> None:
    images = [Image.open(path).convert("RGB") for path in generated]
    preview_width = 700
    previews: list[Image.Image] = []
    for source in images:
        height = round(source.height * preview_width / source.width)
        previews.append(source.resize((preview_width, height), Image.Resampling.LANCZOS))
    label_height = 58
    gap = 28
    canvas = Image.new(
        "RGB",
        (preview_width * len(previews) + gap * (len(previews) - 1), max(x.height for x in previews) + label_height),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(str(FONT_BOLD), 25)
    for index, (name, preview) in enumerate(zip(VARIANTS, previews)):
        x = index * (preview_width + gap)
        draw.text((x + preview_width / 2, 25), name, fill=COLORS["text"], font=font, anchor="mm")
        canvas.paste(preview, (x, label_height))
    canvas.save(output / "comparison_three_layouts.png", dpi=(200, 200), optimize=True)


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
    for state in CONTRASTS:
        maximum = max(float(row["log2_enrichment"]) for row in selected[state])
        if maximum > AXES[state][0]:
            raise ValueError(f"{state} bar exceeds fixed axis maximum: {maximum}")

    args.output.mkdir(parents=True, exist_ok=True)
    generated_pngs: list[Path] = []
    all_parameters: dict[str, object] = {}
    for variant in VARIANTS:
        height, scene, parameters = build_scene(selected, variant)
        stem = f"AME_1p5_500bp_separate_axes_{variant}"
        png = args.output / f"{stem}.png"
        render_png(png, height, scene)
        render_pdf(args.output / f"{stem}.pdf", height, scene)
        render_svg(args.output / f"{stem}.svg", height, scene)
        generated_pngs.append(png)
        all_parameters[variant] = parameters

    write_plot_data(args.output / "plot_data_significant_enrichment_rank.tsv", selected)
    (args.output / "design_parameters.json").write_text(
        json.dumps(all_parameters, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    make_comparison(args.output, generated_pngs)
    print(json.dumps(all_parameters, ensure_ascii=False))


if __name__ == "__main__":
    main()
