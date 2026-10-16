#!/usr/bin/env python3
"""Weaken-first AME drafts with clean independent axes and an in-panel legend."""

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


WIDTH = 900
PLOT_LEFT = 178
PLOT_RIGHT = 858
PLOT_WIDTH = PLOT_RIGHT - PLOT_LEFT
WINDOW_COUNTS = {"strength": 3760, "weaken": 1398}
AXES = {
    "weaken": (1.8, [0.0, 0.6, 1.2, 1.8]),
    "strength": (0.9, [0.0, 0.3, 0.6, 0.9]),
}
VARIANTS = {
    "v6a_inline_legend": {"legend_layout": "inline", "weaken_step": 36.0},
    "v6b_stacked_legend": {"legend_layout": "stacked", "weaken_step": 40.0},
    "v7_compact_equal_spacing": {
        "legend_layout": "inline",
        "row_step": 22.0,
        "compact": True,
    },
}


def build_scene(
    selected: dict[str, list[dict[str, object]]],
    variant: str,
) -> tuple[int, list[dict[str, object]], dict[str, object]]:
    settings = VARIANTS[variant]
    compact = bool(settings.get("compact", False))
    if compact:
        # A single row pitch and a single bar height are used across both
        # states. This prevents the five-row Weaken panel from being stretched
        # merely because it contains fewer motifs than Strength.
        row_step = float(settings["row_step"])
        weaken_step = row_step
        strength_step = row_step
        title_y = 26.0
        weaken_header_y = 60.0
        weaken_first_y = 92.0
        axis_tail = 20.0
        weaken_axis_y = weaken_first_y + (len(selected["weaken"]) - 0.5) * row_step + axis_tail
        strength_header_y = weaken_axis_y + 42.0
        strength_first_y = strength_header_y + 32.0
        strength_axis_y = strength_first_y + (len(selected["strength"]) - 0.5) * row_step + axis_tail
        common_label_y = strength_axis_y + 38.0
        note_y = common_label_y + 28.0
        height = int(note_y + 24.0)
    else:
        weaken_step = float(settings["weaken_step"])
        strength_step = 24.0
        title_y = 28.0

        weaken_header_y = 71.0
        weaken_first_y = 112.0
        weaken_axis_y = weaken_first_y + (len(selected["weaken"]) - 0.5) * weaken_step + 25.0
        strength_header_y = weaken_axis_y + 83.0
        strength_first_y = strength_header_y + 41.0
        strength_axis_y = strength_first_y + (len(selected["strength"]) - 0.5) * strength_step + 25.0
        common_label_y = strength_axis_y + 52.0
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

    text_item(WIDTH / 2, title_y, "AME motif enrichment", 21.0, COLORS["text"], "middle", True)

    geometry = {
        "weaken": (weaken_header_y, weaken_first_y, weaken_axis_y, weaken_step),
        "strength": (strength_header_y, strength_first_y, strength_axis_y, strength_step),
    }
    for state in ("weaken", "strength"):
        header_y, first_y, axis_y, row_step = geometry[state]
        label = "Weaken" if state == "weaken" else "Strength"
        text_item(PLOT_LEFT, header_y, label, 17.0, COLORS[state], "start", True)

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
            text_item(PLOT_LEFT - 13.0, center_y, clean_label(row), 13.0, COLORS["text"], "end")
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

    # Legend sits below the long CTCF bar, in the empty space to the right of
    # the shorter SPIB/Spi1 bars.
    legend_x = PLOT_LEFT + 0.94 / AXES["weaken"][0] * PLOT_WIDTH
    if settings["legend_layout"] == "inline":
        legend_y = weaken_first_y + weaken_step
        entries = (("weaken", legend_x, legend_y), ("strength", legend_x + 170.0, legend_y))
    else:
        entries = (
            ("weaken", legend_x + 55.0, weaken_first_y + weaken_step * 0.90),
            ("strength", legend_x + 55.0, weaken_first_y + weaken_step * 1.75),
        )
    for state, x, y in entries:
        scene.append({"kind": "rect", "x": x, "y": y - 8.0, "width": 16.0, "height": 16.0, "color": COLORS[state]})
        label = "Weaken" if state == "weaken" else "Strength"
        text_item(x + 24.0, y, label, 13.0, COLORS[state], "start", True)
        text_item(x + 88.0, y, f"n={WINDOW_COUNTS[state]:,}", 10.5, COLORS["axis"], "start")

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

    ctcf = float(selected["weaken"][0]["log2_enrichment"])
    atf4 = float(selected["strength"][0]["log2_enrichment"])
    visual_ratio = (ctcf / AXES["weaken"][0]) / (atf4 / AXES["strength"][0])
    parameters = {
        "variant": variant,
        "panel_order": ["weaken", "strength"],
        "selection": "representative motif families with enrichment_ratio > 1 and within-contrast BH q <= 0.05",
        "ordering": "log2 enrichment descending within each state",
        "non_significant_displayed": False,
        "axis_policy": "independent fixed clean axes",
        "axis_max": {state: AXES[state][0] for state in ("weaken", "strength")},
        "CTCF_to_ATF4_visual_length_ratio": visual_ratio,
        "legend_position": "below CTCF in right-side empty area",
        "bar_height": min(18.0, strength_step * 0.70),
        "strength_row_step": strength_step,
        "weaken_row_step": weaken_step,
        "canvas_width": WIDTH,
        "canvas_height": height,
        **settings,
    }
    return height, scene, parameters


def comparison(output: Path, generated: list[Path]) -> None:
    images = [Image.open(path).convert("RGB") for path in generated]
    preview_width = 780
    previews: list[Image.Image] = []
    for source in images:
        preview_height = round(source.height * preview_width / source.width)
        previews.append(source.resize((preview_width, preview_height), Image.Resampling.LANCZOS))
    gap = 32
    label_height = 58
    canvas = Image.new(
        "RGB",
        (preview_width * len(previews) + gap, max(image.height for image in previews) + label_height),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(str(FONT_BOLD), 25)
    for index, (name, preview) in enumerate(zip(VARIANTS, previews)):
        x = index * (preview_width + gap)
        draw.text((x + preview_width / 2, 25), name, fill=COLORS["text"], font=font, anchor="mm")
        canvas.paste(preview, (x, label_height))
    canvas.save(output / "comparison_weaken_first_legend_layouts.png", dpi=(200, 200), optimize=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--variant", choices=tuple(VARIANTS), action="append")
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
    pngs: list[Path] = []
    parameters: dict[str, object] = {}
    variants = args.variant or list(VARIANTS)
    for variant in variants:
        height, scene, variant_parameters = build_scene(selected, variant)
        stem = f"AME_1p5_500bp_weaken_first_{variant}"
        png = args.output / f"{stem}.png"
        render_png(png, height, scene)
        render_pdf(args.output / f"{stem}.pdf", height, scene)
        render_svg(args.output / f"{stem}.svg", height, scene)
        pngs.append(png)
        parameters[variant] = variant_parameters

    write_plot_data(args.output / "plot_data_significant_enrichment_rank.tsv", selected)
    (args.output / "design_parameters.json").write_text(
        json.dumps(parameters, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    if len(pngs) > 1:
        comparison(args.output, pngs)
    print(json.dumps(parameters, ensure_ascii=False))


if __name__ == "__main__":
    main()
