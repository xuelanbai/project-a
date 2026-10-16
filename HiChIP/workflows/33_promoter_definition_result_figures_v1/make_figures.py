#!/usr/bin/env python3
"""Create publication-style promoter-definition sensitivity figures."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas
from reportlab.lib.colors import Color, HexColor, white


DEFINITIONS = [
    ("TSS_pm2kb", "TSS +/- 2 kb", "TSS"),
    ("P1_full_H3_peak_support", "P1: full H3K27ac peak support", "P1"),
    ("P2_TSS_local_support_full_H3", "P2: TSS-local support", "P2"),
]
STRENGTH = "#D44C43"
WEAKEN = "#3F88AE"
GREY = "#6F7378"
LIGHT_BLUE = "#7FB4CF"
LIGHT_RED = "#D98279"
MOTIF_ORDER = {
    "strength": [
        "ISR/ATF4", "NF-kB", "CREB1", "TFEC", "ATF2", "TFE3", "C/EBP", "MITF",
        "Nr2f6", "USF1", "TFEB", "AP-1", "Arntl", "BACH1", "Atf3", "BNC2",
    ],
    "weaken": ["CTCF", "SPIB", "Spi1", "ETV7", "SPIC"],
}
GENE_PREFERENCE = [
    "FOXO1", "THRB", "H6PD", "CERS5", "BNIP3L", "LPIN3", "LPL", "PPARD",
    "SREBF2", "SREBF1", "FABP3", "NDRG1",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def candidate_gene_counts(root: Path) -> dict[tuple[str, str], int]:
    counts = {}
    for definition, _, _ in DEFINITIONS:
        for state in ("strength", "weaken"):
            rows = read_tsv(root / "targets" / definition / state / "08_candidate_genes.tsv")
            counts[definition, state] = len({row["gene_id"] for row in rows})
    return counts


def bridge_counts(root: Path) -> dict[tuple[str, str, str], int]:
    return {
        (row["promoter_definition"], row["state"], row["category"]): int(row["genes"])
        for row in read_tsv(root / "04_bridge_matrix_summary.tsv")
    }


def promoter_counts(root: Path) -> dict[tuple[str, str], int]:
    counts = {}
    for row in read_tsv(root / "02_promoter_definition_counts.tsv"):
        if row["time"] in {"not_time_specific", "0h", "6h"}:
            counts[row["promoter_definition"], row["time"]] = int(row["unique_output_intervals"])
    return counts


def load_reference_matrix(root: Path, definition: str) -> dict[tuple[str, str, str], int]:
    values = {}
    for state in ("strength", "weaken"):
        path = root / "targets" / definition / f"{state}_AME_reference_panel_bridge.tsv"
        for row in read_tsv(path):
            gene = row["gene_symbol"]
            for motif in MOTIF_ORDER[state]:
                values[gene, state, motif] = int(row.get(motif, "0") or 0)
    return values


def detected_gene_order(matrices: dict[str, dict]) -> list[str]:
    found = {
        gene for matrix in matrices.values() for (gene, _, _), value in matrix.items() if value > 0
    }
    return [gene for gene in GENE_PREFERENCE if gene in found]


def blend(color: str, fraction: float) -> Color:
    base = HexColor(color)
    f = min(1.0, max(0.0, fraction))
    return Color(1 - f * (1 - base.red), 1 - f * (1 - base.green), 1 - f * (1 - base.blue))


def draw_text(c: canvas.Canvas, x: float, y: float, text: str, size: float = 9,
              bold: bool = False, color: str = "#222222", align: str = "left") -> None:
    font = "Helvetica-Bold" if bold else "Helvetica"
    c.setFont(font, size)
    c.setFillColor(HexColor(color))
    if align == "center":
        c.drawCentredString(x, y, text)
    elif align == "right":
        c.drawRightString(x, y, text)
    else:
        c.drawString(x, y, text)


def nice_ticks(ymax: float, count: int = 4) -> list[float]:
    raw = ymax / count
    power = 10 ** math.floor(math.log10(raw))
    fraction = raw / power
    step = (1 if fraction <= 1 else 2 if fraction <= 2 else 5 if fraction <= 5 else 10) * power
    return [value for value in (index * step for index in range(count + 2)) if value <= ymax]


def grouped_bar(c: canvas.Canvas, box: tuple[float, float, float, float], categories: list[str],
                series: list[tuple[str, list[float | None], str]], ymax: float, title: str,
                ylabel: str, panel: str) -> None:
    x, y, w, h = box
    left, right, bottom, top = 55, 12, 45, 38
    px, py = x + left, y + bottom
    pw, ph = w - left - right, h - bottom - top
    draw_text(c, x, y + h - 5, panel, 13, True)
    draw_text(c, x + 21, y + h - 5, title, 12, True)
    c.setStrokeColor(HexColor("#DADADA"))
    c.setLineWidth(0.6)
    for tick in nice_ticks(ymax):
        ty = py + ph * tick / ymax
        c.line(px, ty, px + pw, ty)
        draw_text(c, px - 7, ty - 3, f"{int(tick):,}", 7.5, align="right", color="#555555")
    group_width = pw / len(categories)
    visible_series = len(series)
    bar_width = min(33, group_width * 0.68 / visible_series)
    for index, category in enumerate(categories):
        center = px + group_width * (index + 0.5)
        draw_text(c, center, py - 17, category, 8.3, align="center")
        for sindex, (_, values, color) in enumerate(series):
            value = values[index]
            if value is None or value <= 0:
                continue
            bx = center + (sindex - (visible_series - 1)/2) * bar_width - bar_width * 0.42
            bh = ph * value / ymax
            c.setFillColor(HexColor(color))
            c.rect(bx, py, bar_width * 0.84, bh, fill=1, stroke=0)
            draw_text(c, bx + bar_width * 0.42, py + bh + 5, f"{int(value):,}", 7.8, align="center")
    c.saveState()
    c.translate(x + 11, py + ph/2)
    c.rotate(90)
    draw_text(c, 0, 0, ylabel, 8.5, align="center", color="#444444")
    c.restoreState()
    legend_x = x + w - 10
    for name, _, color in reversed(series):
        width = stringWidth(name, "Helvetica", 7.5)
        legend_x -= width + 19
        c.setFillColor(HexColor(color))
        c.rect(legend_x, y + h - 25, 9, 9, fill=1, stroke=0)
        draw_text(c, legend_x + 13, y + h - 24, name, 7.5)
    c.setStrokeColor(HexColor("#7A7A7A"))
    c.setLineWidth(0.8)
    c.line(px, py, px + pw, py)
    c.line(px, py, px, py + ph)


def make_comparison(root: Path, out: Path) -> None:
    promoters = promoter_counts(root)
    genes = candidate_gene_counts(root)
    bridges = bridge_counts(root)
    keys = [item[0] for item in DEFINITIONS]
    labels = [item[2] for item in DEFINITIONS]
    path = out / "Figure_1_promoter_definition_comparison.pdf"
    c = canvas.Canvas(str(path), pagesize=(1000, 720))
    c.setTitle("Promoter definition sensitivity")
    draw_text(c, 500, 687, "Promoter definition sensitivity", 20, True, align="center")
    grouped_bar(c, (35, 380, 450, 270), labels, [
        ("time-independent", [promoters[keys[0], "not_time_specific"], None, None], GREY),
        ("0 h", [None, promoters[keys[1], "0h"], promoters[keys[2], "0h"]], LIGHT_BLUE),
        ("6 h", [None, promoters[keys[1], "6h"], promoters[keys[2], "6h"]], LIGHT_RED),
    ], 39_000, "Promoter intervals", "Unique intervals", "a")
    grouped_bar(c, (515, 380, 450, 270), labels, [
        ("Strength", [genes[key, "strength"] for key in keys], STRENGTH),
        ("Weaken", [genes[key, "weaken"] for key in keys], WEAKEN),
    ], 4_000, "Motif-linked target genes", "Unique genes", "b")
    grouped_bar(c, (35, 65, 450, 270), labels, [
        ("Lipid", [bridges[key, "strength", "lipid"] for key in keys], "#E9A39C"),
        ("Inflammation", [bridges[key, "strength", "inflammation"] for key in keys], STRENGTH),
    ], 700, "Strength functional targets", "Linked genes", "c")
    grouped_bar(c, (515, 65, 450, 270), labels, [
        ("Lipid", [bridges[key, "weaken", "lipid"] for key in keys], "#A8CDDF"),
        ("Inflammation", [bridges[key, "weaken", "inflammation"] for key in keys], WEAKEN),
    ], 500, "Weaken functional targets", "Linked genes", "d")
    draw_text(c, 500, 23, "TSS is time-independent; P1/P2 use 6 h for Strength and 0 h for Weaken.",
              8.5, align="center", color="#555555")
    c.showPage()
    c.save()


def matrix_geometry(box: tuple[float, float, float, float], genes: list[str]):
    x, y, w, h = box
    left, right, bottom, top = 85, 15, 100, 72
    gx, gy = x + left, y + bottom
    gw, gh = w - left - right, h - bottom - top
    return gx, gy, gw, gh, gw / 21, gh / len(genes)


def draw_matrix(c: canvas.Canvas, box: tuple[float, float, float, float], matrix: dict,
                genes: list[str], title: str, global_max: int, panel: str | None = None,
                show_xlabel: bool = True) -> None:
    x, y, w, h = box
    gx, gy, gw, gh, cw, ch = matrix_geometry(box, genes)
    motifs = MOTIF_ORDER["strength"] + MOTIF_ORDER["weaken"]
    states = ["strength"] * 16 + ["weaken"] * 5
    if panel:
        draw_text(c, x + 4, y + h - 15, panel, 14, True)
    draw_text(c, x + w/2, y + h - 14, title, 15, True, align="center")
    for row, gene in enumerate(genes):
        cy = gy + gh - (row + 1) * ch
        draw_text(c, gx - 9, cy + ch/2 - 3, gene, 9.5, bold=True, align="right")
        for col, (motif, state) in enumerate(zip(motifs, states)):
            value = matrix.get((gene, state, motif), 0)
            fraction = math.log1p(value) / math.log1p(global_max) if value else 0
            color = STRENGTH if state == "strength" else WEAKEN
            c.setFillColor(blend(color, 0.12 + 0.88 * fraction) if value else white)
            c.setStrokeColor(HexColor("#D8D8D8"))
            c.setLineWidth(0.45)
            c.rect(gx + col*cw, cy, cw, ch, fill=1, stroke=1)
            if value:
                text_color = "#FFFFFF" if fraction > 0.56 else "#222222"
                draw_text(c, gx + (col+0.5)*cw, cy + ch/2 - 3, str(value), 8,
                          bold=value >= 3, color=text_color, align="center")
    split = 16
    strip_y = gy + gh + 9
    c.setFillColor(HexColor(STRENGTH))
    c.rect(gx, strip_y, split*cw, 6, fill=1, stroke=0)
    c.setFillColor(HexColor(WEAKEN))
    c.rect(gx + split*cw, strip_y, 5*cw, 6, fill=1, stroke=0)
    draw_text(c, gx + split*cw/2, strip_y + 13, "Strength motifs", 8.5, True,
              color=STRENGTH, align="center")
    draw_text(c, gx + split*cw + 2.5*cw, strip_y + 13, "Weaken motifs", 8.5, True,
              color=WEAKEN, align="center")
    c.setStrokeColor(HexColor("#666666"))
    c.setLineWidth(1.4)
    c.line(gx + split*cw, gy, gx + split*cw, gy + gh)
    if show_xlabel:
        for col, motif in enumerate(motifs):
            c.saveState()
            c.translate(gx + (col + 0.68)*cw, gy - 8)
            c.rotate(55)
            draw_text(c, 0, 0, motif, 7.6, align="right")
            c.restoreState()


def draw_legend(c: canvas.Canvas, page_width: float, y: float, global_max: int) -> None:
    draw_text(c, page_width/2, y + 26, "Cell value = unique motif-bearing enhancer count linked by HiChIP",
              8.5, align="center", color="#444444")
    examples = sorted({1, min(3, global_max), min(6, global_max), global_max})
    start = page_width/2 - 175
    for state_index, (label, color) in enumerate((("Strength", STRENGTH), ("Weaken", WEAKEN))):
        x = start + state_index * 210
        draw_text(c, x, y + 4, label, 8, bold=True, color=color, align="right")
        for index, value in enumerate(examples):
            fraction = math.log1p(value) / math.log1p(global_max)
            bx = x + 12 + index * 34
            c.setFillColor(blend(color, 0.12 + 0.88*fraction))
            c.setStrokeColor(HexColor("#BBBBBB"))
            c.rect(bx, y, 18, 14, fill=1, stroke=1)
            draw_text(c, bx + 9, y - 9, str(value), 6.8, align="center", color="#555555")


def make_single_matrix(out: Path, definition: str, title: str, matrix: dict,
                       genes: list[str], global_max: int) -> None:
    path = out / f"Figure_2_{definition}_motif_enhancer_gene_matrix.pdf"
    c = canvas.Canvas(str(path), pagesize=(1100, 620))
    c.setTitle(title)
    draw_matrix(c, (30, 65, 1040, 530), matrix, genes, title, global_max, show_xlabel=True)
    draw_legend(c, 1100, 24, global_max)
    c.showPage()
    c.save()


def make_combined(out: Path, matrices: dict[str, dict], genes: list[str], global_max: int) -> None:
    path = out / "Figure_3_combined_three_promoter_matrices.pdf"
    c = canvas.Canvas(str(path), pagesize=(1100, 1500))
    c.setTitle("Promoter-sensitive motif-bearing enhancer-gene links")
    draw_text(c, 550, 1465, "Promoter-sensitive motif-bearing enhancer-gene links", 20, True, align="center")
    boxes = [(30, 990, 1040, 430), (30, 530, 1040, 430), (30, 70, 1040, 430)]
    for index, ((definition, title, _), box) in enumerate(zip(DEFINITIONS, boxes)):
        draw_matrix(c, box, matrices[definition], genes, title, global_max,
                    panel=chr(ord("a") + index), show_xlabel=True)
    draw_legend(c, 1100, 23, global_max)
    c.showPage()
    c.save()


def export_source_data(root: Path, out: Path, matrices: dict[str, dict], genes: list[str]) -> None:
    promoter = promoter_counts(root)
    gene_counts = candidate_gene_counts(root)
    bridges = bridge_counts(root)
    comparison = []
    for definition, title, _ in DEFINITIONS:
        comparison.append({
            "promoter_definition": definition, "display_label": title,
            "promoter_intervals_static": promoter.get((definition, "not_time_specific"), "NA"),
            "promoter_intervals_0h": promoter.get((definition, "0h"), "NA"),
            "promoter_intervals_6h": promoter.get((definition, "6h"), "NA"),
            "strength_motif_linked_genes": gene_counts[definition, "strength"],
            "weaken_motif_linked_genes": gene_counts[definition, "weaken"],
            "strength_lipid_genes": bridges[definition, "strength", "lipid"],
            "strength_inflammation_genes": bridges[definition, "strength", "inflammation"],
            "weaken_lipid_genes": bridges[definition, "weaken", "lipid"],
            "weaken_inflammation_genes": bridges[definition, "weaken", "inflammation"],
        })
    write_tsv(out / "Figure_1_source_data.tsv", comparison, list(comparison[0]))
    matrix_rows = []
    for definition, title, _ in DEFINITIONS:
        for gene in genes:
            for state in ("strength", "weaken"):
                for motif in MOTIF_ORDER[state]:
                    matrix_rows.append({
                        "promoter_definition": definition, "display_label": title,
                        "gene_symbol": gene, "state": state, "motif_family": motif,
                        "unique_CRE_count": matrices[definition].get((gene, state, motif), 0),
                    })
    write_tsv(out / "Figure_2_3_source_data.tsv", matrix_rows,
              ["promoter_definition", "display_label", "gene_symbol", "state", "motif_family", "unique_CRE_count"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, out = args.run.resolve(), args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    matrices = {definition: load_reference_matrix(root, definition) for definition, _, _ in DEFINITIONS}
    genes = detected_gene_order(matrices)
    if not genes:
        raise ValueError("No predefined reference-panel genes were linked")
    global_max = max(value for matrix in matrices.values() for value in matrix.values())
    make_comparison(root, out)
    for definition, title, _ in DEFINITIONS:
        make_single_matrix(out, definition, title, matrices[definition], genes, global_max)
    make_combined(out, matrices, genes, global_max)
    export_source_data(root, out, matrices, genes)
    print(f"Created figures in {out}")


if __name__ == "__main__":
    main()
