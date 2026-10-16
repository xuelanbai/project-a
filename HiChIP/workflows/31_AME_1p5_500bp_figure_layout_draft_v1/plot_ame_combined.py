#!/usr/bin/env python3
"""Create compact combined Strength/Weaken AME panels for Figure 3 review."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


CONTRASTS = {
    "strength": "strength_1p5_500bp_ATAC_core",
    "weaken": "weaken_1p5_500bp_ATAC_core",
}

COLORS = {
    "strength": "#CF4A43",
    "strength_pale": "#EFC7C3",
    "weaken": "#3E82AA",
    "weaken_pale": "#C8DDE9",
    "text": "#222222",
    "grid": "#E6E6E6",
    "axis": "#6B6B6B",
}

WIDTH = 900
BASE_HEIGHT = 150
ROW_HEIGHT = 18
PDF_SCALE = 0.72
PNG_SCALE = 3
FONT_REGULAR = Path("C:/Windows/Fonts/arial.ttf")
FONT_BOLD = Path("C:/Windows/Fonts/arialbd.ttf")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "display_order", "state", "contrast", "motif_family", "motif_id",
        "motif_alt_id", "enrichment_ratio", "log2_enrichment",
        "AME_BH_q_within_contrast", "AME_BH_q_all_contrasts", "significant",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def clean_label(row: dict[str, object]) -> str:
    family = str(row["motif_family"]).strip()
    representative = str(row["motif_alt_id"]).strip() or str(row["motif_id"]).strip()
    if family.casefold() == representative.casefold():
        return family
    return f"{family}  {representative}"


def select_rows(rows: list[dict[str, str]], state: str, ordering: str, top_n: int) -> list[dict[str, object]]:
    contrast = CONTRASTS[state]
    parsed: list[dict[str, object]] = []
    for row in rows:
        if row.get("contrast") != contrast:
            continue
        try:
            ratio = float(row["enrichment_ratio"])
            q_within = float(row["AME_BH_q_within_contrast"])
            q_all = float(row["AME_BH_q_all_contrasts"])
        except (TypeError, ValueError):
            continue
        if not math.isfinite(ratio) or ratio <= 1:
            continue
        item: dict[str, object] = dict(row)
        item["enrichment_ratio"] = ratio
        item["log2_enrichment"] = math.log2(ratio)
        item["q_within"] = q_within
        item["q_all"] = q_all
        item["significant"] = item["q_within"] <= 0.05
        item["state"] = state
        parsed.append(item)

    if ordering == "q_rank":
        parsed.sort(key=lambda x: (float(x["q_within"]), -float(x["log2_enrichment"]), str(x["motif_id"])))
    elif ordering == "significant_enrichment_rank":
        parsed = [row for row in parsed if bool(row["significant"])]
        parsed.sort(key=lambda x: (-float(x["log2_enrichment"]), float(x["q_within"]), str(x["motif_id"])))
    else:
        raise ValueError(ordering)
    return parsed[:top_n]


def nice_axis(maximum: float) -> tuple[float, list[float]]:
    target = maximum / 5
    steps = [0.05, 0.1, 0.2, 0.25, 0.5, 1.0, 2.0]
    step = next((candidate for candidate in steps if candidate >= target), steps[-1])
    limit = math.ceil(maximum / step) * step
    ticks = [index * step for index in range(int(round(limit / step)) + 1)]
    return limit, ticks


def build_scene(selected: dict[str, list[dict[str, object]]], ordering: str) -> tuple[int, list[dict[str, object]], dict[str, float]]:
    max_rows = max(len(selected["strength"]), len(selected["weaken"]))
    height = BASE_HEIGHT + ROW_HEIGHT * max_rows
    plot_top = 91
    plot_bottom = height - 67
    panel_origins = {"strength": 18, "weaken": 462}
    plot_left_offset = 145
    plot_width = 270
    row_step = (plot_bottom - plot_top) / max_rows
    bar_height = min(12.0, row_step * 0.68)
    raw_maxima = {
        state: max(float(row["log2_enrichment"]) for row in selected[state])
        for state in CONTRASTS
    }
    if ordering == "q_rank":
        axes = {state: nice_axis(raw_maxima[state] * 1.04) for state in CONTRASTS}
    else:
        shared = nice_axis(max(raw_maxima.values()) * 1.04)
        axes = {state: shared for state in CONTRASTS}
    scene: list[dict[str, object]] = []

    def text_item(x: float, y: float, value: str, size: float, color: str, anchor: str = "start", bold: bool = False) -> None:
        scene.append({"kind": "text", "x": x, "y": y, "text": value, "size": size,
                      "color": color, "anchor": anchor, "bold": bold})

    text_item(WIDTH / 2, 21, "AME motif enrichment (1.5-fold H3K27ac, 500-bp ATAC core)", 11.2, COLORS["text"], "middle", True)
    for state in ("strength", "weaken"):
        origin = panel_origins[state]
        plot_left = origin + plot_left_offset
        plot_right = plot_left + plot_width
        saturated = COLORS[state]
        pale = COLORS[f"{state}_pale"]
        title = "Strength" if state == "strength" else "Weaken"
        subtitle = "q-ranked representatives" if ordering == "q_rank" else "significant representatives ranked by enrichment"
        n_windows = 3760 if state == "strength" else 1398
        text_item(origin, 48, title, 11.0, saturated, "start", True)
        text_item(origin, 65, f"{subtitle}  (n={n_windows:,} windows)", 7.7, COLORS["axis"])

        axis_max, ticks = axes[state]
        for tick in ticks:
            x = plot_left + tick / axis_max * plot_width
            scene.append({"kind": "line", "x1": x, "y1": plot_top - 3, "x2": x, "y2": plot_bottom,
                          "color": COLORS["grid"], "width": 0.55})
            text_item(x, plot_bottom + 17, f"{tick:g}", 7.2, COLORS["axis"], "middle")
        scene.append({"kind": "line", "x1": plot_left, "y1": plot_bottom, "x2": plot_right, "y2": plot_bottom,
                      "color": COLORS["axis"], "width": 0.7})

        for index, row in enumerate(selected[state]):
            center_y = plot_top + (index + 0.5) * row_step
            value = float(row["log2_enrichment"])
            color = saturated if bool(row["significant"]) else pale
            text_item(plot_left - 7, center_y, clean_label(row), 7.6, COLORS["text"], "end")
            scene.append({"kind": "rect", "x": plot_left, "y": center_y - bar_height / 2,
                          "width": value / axis_max * plot_width, "height": bar_height, "color": color})

    text_item(WIDTH / 2, height - 37, "log2(foreground / GC-matched control)", 8.3, COLORS["text"], "middle")
    note = "Panel-specific x-axes; saturated bars: BH q <= 0.05; pale bars: BH q > 0.05" if ordering == "q_rank" else "Shared x-axis; BH q <= 0.05 within each contrast"
    text_item(WIDTH / 2, height - 17, note, 7.2, COLORS["axis"], "middle")
    return height, scene, {state: axes[state][0] for state in CONTRASTS}


def pil_font(size: float, bold: bool, scale: int) -> ImageFont.FreeTypeFont:
    path = FONT_BOLD if bold else FONT_REGULAR
    return ImageFont.truetype(str(path), max(1, round(size * scale)))


def render_png(path: Path, height: int, scene: list[dict[str, object]]) -> None:
    image = Image.new("RGB", (WIDTH * PNG_SCALE, height * PNG_SCALE), "white")
    draw = ImageDraw.Draw(image)
    for item in scene:
        kind = item["kind"]
        if kind == "rect":
            box = tuple(round(float(item[key]) * PNG_SCALE) for key in ("x", "y", "width", "height"))
            x, y, width, rect_height = box
            draw.rectangle((x, y, x + width, y + rect_height), fill=str(item["color"]))
        elif kind == "line":
            draw.line(tuple(round(float(item[key]) * PNG_SCALE) for key in ("x1", "y1", "x2", "y2")),
                      fill=str(item["color"]), width=max(1, round(float(item["width"]) * PNG_SCALE)))
        else:
            font = pil_font(float(item["size"]), bool(item["bold"]), PNG_SCALE)
            anchor = {"start": "lm", "middle": "mm", "end": "rm"}[str(item["anchor"])]
            draw.text((float(item["x"]) * PNG_SCALE, float(item["y"]) * PNG_SCALE), str(item["text"]),
                      font=font, fill=str(item["color"]), anchor=anchor)
    image.save(path, dpi=(300, 300), optimize=True)


def render_svg(path: Path, height: int, scene: list[dict[str, object]]) -> None:
    chunks = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" viewBox="0 0 {WIDTH} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    for item in scene:
        if item["kind"] == "rect":
            chunks.append(f'<rect x="{item["x"]}" y="{item["y"]}" width="{item["width"]}" height="{item["height"]}" fill="{item["color"]}"/>')
        elif item["kind"] == "line":
            chunks.append(f'<line x1="{item["x1"]}" y1="{item["y1"]}" x2="{item["x2"]}" y2="{item["y2"]}" stroke="{item["color"]}" stroke-width="{item["width"]}"/>')
        else:
            anchor = {"start": "start", "middle": "middle", "end": "end"}[str(item["anchor"])]
            weight = "700" if bool(item["bold"]) else "400"
            chunks.append(f'<text x="{item["x"]}" y="{item["y"]}" fill="{item["color"]}" font-family="Arial, Helvetica, sans-serif" font-size="{item["size"]}" font-weight="{weight}" text-anchor="{anchor}" dominant-baseline="middle">{html.escape(str(item["text"]))}</text>')
    chunks.append("</svg>")
    path.write_text("\n".join(chunks) + "\n", encoding="utf-8")


def render_pdf(path: Path, height: int, scene: list[dict[str, object]]) -> None:
    if FONT_REGULAR.is_file() and FONT_BOLD.is_file():
        pdfmetrics.registerFont(TTFont("FigureArial", str(FONT_REGULAR)))
        pdfmetrics.registerFont(TTFont("FigureArialBold", str(FONT_BOLD)))
    regular = "FigureArial" if "FigureArial" in pdfmetrics.getRegisteredFontNames() else "Helvetica"
    bold = "FigureArialBold" if "FigureArialBold" in pdfmetrics.getRegisteredFontNames() else "Helvetica-Bold"
    pdf = canvas.Canvas(str(path), pagesize=(WIDTH * PDF_SCALE, height * PDF_SCALE))
    pdf.scale(PDF_SCALE, PDF_SCALE)
    for item in scene:
        if item["kind"] == "rect":
            pdf.setFillColor(str(item["color"]))
            pdf.rect(float(item["x"]), height - float(item["y"]) - float(item["height"]),
                     float(item["width"]), float(item["height"]), stroke=0, fill=1)
        elif item["kind"] == "line":
            pdf.setStrokeColor(str(item["color"]))
            pdf.setLineWidth(float(item["width"]))
            pdf.line(float(item["x1"]), height - float(item["y1"]), float(item["x2"]), height - float(item["y2"]))
        else:
            size = float(item["size"])
            pdf.setFillColor(str(item["color"]))
            pdf.setFont(bold if bool(item["bold"]) else regular, size)
            x = float(item["x"])
            y = height - float(item["y"]) - size * 0.34
            value = str(item["text"])
            if item["anchor"] == "middle":
                pdf.drawCentredString(x, y, value)
            elif item["anchor"] == "end":
                pdf.drawRightString(x, y, value)
            else:
                pdf.drawString(x, y, value)
    pdf.showPage()
    pdf.save()


def plot(rows: list[dict[str, str]], output: Path, ordering: str, top_n: int) -> list[dict[str, object]]:
    selected = {state: select_rows(rows, state, ordering, top_n) for state in CONTRASTS}
    if not all(selected.values()):
        raise ValueError(f"No rows selected for one or more states: {ordering}")
    height, scene, axis_maxima = build_scene(selected, ordering)
    stem = f"AME_1p5_500bp_combined_{ordering}_layout_v1"
    output.mkdir(parents=True, exist_ok=True)
    render_png(output / f"{stem}.png", height, scene)
    render_pdf(output / f"{stem}.pdf", height, scene)
    render_svg(output / f"{stem}.svg", height, scene)

    exported: list[dict[str, object]] = []
    for state in ("strength", "weaken"):
        for index, row in enumerate(selected[state], 1):
            exported.append({
                "display_order": index,
                "state": state,
                "contrast": row["contrast"],
                "motif_family": row["motif_family"],
                "motif_id": row["motif_id"],
                "motif_alt_id": row["motif_alt_id"],
                "enrichment_ratio": row["enrichment_ratio"],
                "log2_enrichment": row["log2_enrichment"],
                "AME_BH_q_within_contrast": row["AME_BH_q_within_contrast"],
                "AME_BH_q_all_contrasts": row["AME_BH_q_all_contrasts"],
                "significant": int(bool(row["significant"])),
            })
    write_rows(output / f"plot_data_{ordering}.tsv", exported)
    for row in exported:
        row["axis_max"] = axis_maxima[str(row["state"])]
    return exported


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()
    if args.top_n < 1:
        raise ValueError("--top-n must be positive")
    rows = read_rows(args.input)
    counts = {}
    for ordering in ("q_rank", "significant_enrichment_rank"):
        plotted = plot(rows, args.output, ordering, args.top_n)
        counts[ordering] = dict((state, sum(row["state"] == state for row in plotted)) for state in CONTRASTS)
    parameters = {
        "input": str(args.input.resolve()),
        "fold_change": 1.5,
        "motif_window_bp": 500,
        "top_n_per_state": args.top_n,
        "significance_field": "AME_BH_q_within_contrast",
        "significance_cutoff": 0.05,
        "strength_color": COLORS["strength"],
        "weaken_color": COLORS["weaken"],
        "layout": "two aligned horizontal panels",
        "axis_policy": {
            "q_rank": "panel-specific axes because non-significant high ratios compress Strength",
            "significant_enrichment_rank": "shared axis",
        },
        "counts_plotted": counts,
    }
    (args.output / "design_parameters.json").write_text(
        json.dumps(parameters, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(counts, ensure_ascii=False))


if __name__ == "__main__":
    main()
