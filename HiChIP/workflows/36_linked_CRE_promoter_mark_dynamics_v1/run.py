#!/usr/bin/env python3
"""Quantify multi-mark dynamics in final enhancer, P1 promoter and linked pairs.

The task is descriptive because pooled/merged BAMs are used.  Every mark is
counted in the exact same genomic interval at 0 h and 6 h.  Counts are scaled
to 10 million proper-pair main-chromosome alignments, without division by
interval length and without global median centering.  ChIP assays additionally
receive an IgG-relative signal and its 0 h to 6 h change.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone


MAIN = {str(x) for x in range(1, 19)} | {"X", "Y"}
ASSAYS = ("H3K27ac", "H3K4me3", "Pol2", "ATAC", "H3K27me3")
CHIP_ASSAYS = {"H3K27ac", "H3K4me3", "Pol2", "H3K27me3"}
TIMES = ("0h", "6h")


def canon(chrom: str) -> str:
    value = str(chrom).strip()
    return value[3:] if value.startswith("chr") else value


def chrom_key(chrom: str):
    value = canon(chrom)
    return (0, int(value)) if value.isdigit() else (1, value)


def read_tsv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields and not key.startswith("_"):
                    fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_int(row: dict, field: str) -> int:
    try:
        return int(float(row[field]))
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid integer {field}={row.get(field)!r}") from error


def as_float(row: dict, field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid number {field}={row.get(field)!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"Non-finite number {field}={row.get(field)!r}")
    return value


def derive_state(row: dict, fold: float, minimum: int) -> str:
    cutoff = math.log2(fold)
    value = as_float(row, "H3K27ac_log2FC")
    if (
        value >= cutoff and as_int(row, "supported_6h") == 1
        and as_float(row, "H3K27ac_6h_raw_tags") >= minimum
        and as_float(row, "ATAC_6h_raw_tags") >= minimum
    ):
        return "strength"
    if (
        value <= -cutoff and as_int(row, "supported_0h") == 1
        and as_float(row, "H3K27ac_0h_raw_tags") >= minimum
        and as_float(row, "ATAC_0h_raw_tags") >= minimum
    ):
        return "weaken"
    return "stable"


def read_regions(path: Path, region_type: str, fold: float, minimum: int,
                 validate_archived_state=False) -> list[dict]:
    rows = read_tsv(path)
    required = {"region_id", "chrom", "start", "end", "H3K27ac_log2FC",
                "supported_0h", "supported_6h", "H3K27ac_0h_raw_tags",
                "H3K27ac_6h_raw_tags", "ATAC_0h_raw_tags", "ATAC_6h_raw_tags"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Unexpected region matrix schema: {path}")
    result, seen = [], set()
    for raw in rows:
        original_id = raw["region_id"]
        if original_id in seen:
            raise ValueError(f"Duplicate {region_type} region_id: {original_id}")
        seen.add(original_id)
        chrom, start, end = canon(raw["chrom"]), as_int(raw, "start"), as_int(raw, "end")
        if chrom not in MAIN or not 0 <= start < end:
            raise ValueError(f"Invalid interval {chrom}:{start}-{end}")
        state = derive_state(raw, fold, minimum)
        if validate_archived_state and raw.get("promoter_state") != state:
            raise ValueError(f"Archived promoter state mismatch: {original_id}")
        result.append({
            "feature_id": f"{region_type.upper()}|{original_id}",
            "region_type": region_type,
            "region_id": original_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "length_bp": end - start,
            "state": state,
            "source_times": raw.get("source_times", raw.get("P1_support_times", "")),
            "supported_0h": raw["supported_0h"],
            "supported_6h": raw["supported_6h"],
            "source_H3K27ac_log2FC": raw["H3K27ac_log2FC"],
            "source_ATAC_log2FC": raw.get("ATAC_log2FC", "NA"),
        })
    return result


def read_tss_windows(path: Path, half_width: int, lengths: dict[str, int]) -> tuple[list[dict], dict]:
    rows = read_tsv(path)
    required = {"region_id", "chrom", "gene_id", "gene_symbol", "gene_tss", "gene_strand"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Unexpected promoter-gene association schema: {path}")
    by_key, gene_to_features = {}, defaultdict(set)
    for raw in rows:
        chrom, tss = canon(raw["chrom"]), as_int(raw, "gene_tss")
        if chrom not in lengths or not 0 <= tss < lengths[chrom]:
            raise ValueError(f"TSS outside reference: {chrom}:{tss}")
        start, end = max(0, tss - half_width), min(lengths[chrom], tss + half_width)
        key = (raw["gene_id"], chrom, tss)
        feature_id = f"TSS|{raw['gene_id']}|{chrom}|{tss}"
        if key not in by_key:
            by_key[key] = {
                "feature_id": feature_id, "region_type": "TSS_pm2kb",
                "region_id": feature_id, "chrom": chrom, "start": start, "end": end,
                "length_bp": end - start, "state": "gene_TSS", "gene_id": raw["gene_id"],
                "gene_symbol": raw["gene_symbol"], "gene_tss": tss,
                "gene_strand": raw["gene_strand"], "associated_promoter_ids": set(),
                "associated_promoter_states": set(),
            }
        by_key[key]["associated_promoter_ids"].add(raw["region_id"])
        by_key[key]["associated_promoter_states"].add(raw.get("promoter_state", "NA"))
        gene_to_features[raw["gene_id"]].add(feature_id)
    result = []
    for row in by_key.values():
        row["associated_promoter_ids"] = ";".join(sorted(row["associated_promoter_ids"]))
        row["associated_promoter_states"] = ";".join(sorted(row["associated_promoter_states"]))
        result.append(row)
    return result, {key: sorted(value) for key, value in gene_to_features.items()}


def write_bed(path: Path, rows: list[dict], chrom_map: dict[str, str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        for row in sorted(rows, key=lambda x: (chrom_key(x["chrom"]), x["start"], x["end"], x["feature_id"])):
            handle.write(
                f"{chrom_map[row['chrom']]}\t{row['start']}\t{row['end']}\t{row['feature_id']}\n"
            )


def run_command(command, stdout: Path, stderr: Path, command_log: Path):
    command = [str(item) for item in command]
    with command_log.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"command": command, "time": datetime.now(timezone.utc).isoformat()}) + "\n")
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        subprocess.run(command, check=True, stdout=out, stderr=err)


def multicov(bedtools: Path, bed: Path, bams: list[Path], labels: list[str],
             output: Path, command_log: Path, mapq: int) -> dict[str, dict[str, int]]:
    raw = output.with_suffix(".raw.tsv")
    command = [bedtools, "multicov", "-bed", bed, "-bams", *bams, "-q", str(mapq), "-p"]
    run_command(command, raw, output.with_suffix(".stderr.log"), command_log)
    counts = {}
    with raw.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            cols = line.rstrip().split("\t")
            bed_columns = len(cols) - len(labels)
            if bed_columns < 4:
                raise ValueError(f"{raw}:{line_no}: fewer than four BED columns")
            identifier = cols[3]
            if identifier in counts:
                raise ValueError(f"Duplicate multicov feature ID: {identifier}")
            counts[identifier] = {label: int(value) for label, value in zip(labels, cols[-len(labels):])}
    if not counts:
        raise ValueError(f"Empty multicov output: {raw}")
    return counts


def add_signals(rows: list[dict], counts: dict, totals: dict, scale: float, pseudo: float):
    for row in rows:
        values = counts[row["feature_id"]]
        for label, raw in values.items():
            row[f"{label}_raw_tags"] = raw
            row[f"{label}_tags_per_10M"] = raw / totals[label] * scale
        for assay in ASSAYS:
            zero = row[f"{assay}_0h_tags_per_10M"]
            six = row[f"{assay}_6h_tags_per_10M"]
            row[f"{assay}_log2FC"] = math.log2((six + pseudo) / (zero + pseudo))
            if assay in CHIP_ASSAYS:
                enrich0 = math.log2((zero + pseudo) / (row["IgG_0h_tags_per_10M"] + pseudo))
                enrich6 = math.log2((six + pseudo) / (row["IgG_6h_tags_per_10M"] + pseudo))
                row[f"{assay}_log2_vs_IgG_0h"] = enrich0
                row[f"{assay}_log2_vs_IgG_6h"] = enrich6
                row[f"{assay}_delta_log2_vs_IgG"] = enrich6 - enrich0


def source_recount_qc(rows: list[dict]):
    output = []
    for region_type in ("enhancer", "promoter"):
        subset = [row for row in rows if row["region_type"] == region_type]
        for assay in ("H3K27ac", "ATAC"):
            differences = []
            source_field = f"source_{assay}_log2FC"
            for row in subset:
                source = row.get(source_field, "NA")
                if source not in ("", "NA", None):
                    differences.append(abs(float(row[f"{assay}_log2FC"]) - float(source)))
            output.append({
                "region_type": region_type, "assay": assay,
                "regions_compared": len(differences),
                "max_absolute_log2FC_difference": max(differences, default="NA"),
                "median_absolute_log2FC_difference": median(differences),
            })
    return output


def long_mark_table(rows: list[dict]):
    output = []
    for row in rows:
        for assay in ASSAYS:
            for time in TIMES:
                output.append({
                    "feature_id": row["feature_id"], "region_type": row["region_type"],
                    "region_id": row["region_id"], "state": row["state"],
                    "chrom": row["chrom"], "start": row["start"], "end": row["end"],
                    "length_bp": row["length_bp"], "assay": assay, "time": time,
                    "raw_tags": row[f"{assay}_{time}_raw_tags"],
                    "tags_per_10M": row[f"{assay}_{time}_tags_per_10M"],
                    "log2FC_6h_vs_0h": row[f"{assay}_log2FC"],
                    "log2_vs_IgG": row.get(f"{assay}_log2_vs_IgG_{time}", "NA"),
                    "delta_log2_vs_IgG_6h_minus_0h": row.get(f"{assay}_delta_log2_vs_IgG", "NA"),
                })
    return output


def median(values):
    values = [float(value) for value in values if value not in (None, "NA", "") and math.isfinite(float(value))]
    return statistics.median(values) if values else "NA"


def mean(values):
    values = [float(value) for value in values if value not in (None, "NA", "") and math.isfinite(float(value))]
    return statistics.fmean(values) if values else "NA"


def summarize_marks(rows: list[dict]):
    output = []
    groups = defaultdict(list)
    for row in rows:
        groups[(row["region_type"], row["state"])].append(row)
    for (region_type, state), subset in sorted(groups.items()):
        for assay in ASSAYS:
            changes = [float(row[f"{assay}_log2FC"]) for row in subset]
            record = {
                "region_type": region_type, "state": state, "assay": assay,
                "n_regions": len(subset),
                "median_0h_tags_per_10M": median(row[f"{assay}_0h_tags_per_10M"] for row in subset),
                "median_6h_tags_per_10M": median(row[f"{assay}_6h_tags_per_10M"] for row in subset),
                "median_log2FC": median(changes), "mean_log2FC": mean(changes),
                "fraction_log2FC_gt_0": sum(value > 0 for value in changes) / len(changes),
                "fraction_log2FC_lt_0": sum(value < 0 for value in changes) / len(changes),
                "fraction_at_least_1p5_up": sum(value >= math.log2(1.5) for value in changes) / len(changes),
                "fraction_at_least_1p5_down": sum(value <= -math.log2(1.5) for value in changes) / len(changes),
            }
            if assay in CHIP_ASSAYS:
                record["median_delta_log2_vs_IgG"] = median(
                    row[f"{assay}_delta_log2_vs_IgG"] for row in subset
                )
            else:
                record["median_delta_log2_vs_IgG"] = "NA"
            output.append(record)
    return output


def overlap_bp(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def map_loop_links(loop_paths: dict[str, Path], enhancer_by_id: dict, promoter_assoc_path: Path,
                   promoter_by_id: dict, loop_q_max: float):
    associations = read_tsv(promoter_assoc_path)
    by_gene_chrom = defaultdict(list)
    for row in associations:
        promoter_id = row["region_id"]
        if promoter_id not in promoter_by_id:
            raise ValueError(f"Promoter association references unknown ID: {promoter_id}")
        by_gene_chrom[(row["gene_id"], canon(row["chrom"]))].append(row)

    aggregated, failures = {}, Counter()
    for enhancer_state, path in loop_paths.items():
        for raw in read_tsv(path):
            if as_float(raw, "loop_q") > loop_q_max:
                failures["loop_q_above_threshold"] += 1
                continue
            enhancer_id = raw["region_id"]
            enhancer = enhancer_by_id.get(enhancer_id)
            if enhancer is None:
                failures["enhancer_not_in_final_consensus"] += 1
                continue
            if enhancer["state"] != enhancer_state:
                raise ValueError(f"Loop file/state mismatch: {enhancer_id}")
            chrom = canon(raw["promoter_chrom"])
            old_start, old_end = as_int(raw, "promoter_start"), as_int(raw, "promoter_end")
            candidates = []
            for association in by_gene_chrom.get((raw["gene_id"], chrom), []):
                promoter = promoter_by_id[association["region_id"]]
                overlap = overlap_bp(old_start, old_end, promoter["start"], promoter["end"])
                if overlap > 0:
                    candidates.append((association, overlap))
            if not candidates:
                failures["no_current_P1_promoter_for_loop_gene"] += 1
                continue
            for association, overlap in candidates:
                promoter_id = association["region_id"]
                key = (enhancer_state, enhancer_id, promoter_id, raw["gene_id"])
                if key not in aggregated:
                    aggregated[key] = {
                        "enhancer_state": enhancer_state, "enhancer_id": enhancer_id,
                        "promoter_id": promoter_id,
                        "promoter_state": promoter_by_id[promoter_id]["state"],
                        "gene_id": raw["gene_id"], "gene_symbol": raw["gene_symbol"],
                        "gene_tss": raw["gene_tss"], "gene_strand": raw["gene_strand"],
                        "old_promoter_overlap_bp": 0, "loop_ids": set(), "loop_q_values": [],
                    }
                item = aggregated[key]
                item["old_promoter_overlap_bp"] = max(item["old_promoter_overlap_bp"], overlap)
                item["loop_ids"].add(raw["loop_id"])
                item["loop_q_values"].append(as_float(raw, "loop_q"))
    result = []
    for item in aggregated.values():
        item["unique_loop_count"] = len(item["loop_ids"])
        item["best_loop_q"] = min(item["loop_q_values"])
        item["loop_ids"] = ";".join(sorted(item["loop_ids"]))
        item.pop("loop_q_values")
        result.append(item)
    return result, failures


def add_pair_signals(pairs: list[dict], enhancer_by_id: dict, promoter_by_id: dict,
                     tss_by_feature: dict, gene_to_tss: dict):
    output = []
    for pair in pairs:
        enhancer = enhancer_by_id[pair["enhancer_id"]]
        promoter = promoter_by_id[pair["promoter_id"]]
        tss_ids = gene_to_tss.get(pair["gene_id"], [])
        base = dict(pair, TSS_feature_ids=";".join(tss_ids))
        for assay in ASSAYS:
            base[f"enhancer_{assay}_log2FC"] = enhancer[f"{assay}_log2FC"]
            base[f"promoter_{assay}_log2FC"] = promoter[f"{assay}_log2FC"]
            if assay in CHIP_ASSAYS:
                base[f"enhancer_{assay}_delta_log2_vs_IgG"] = enhancer[f"{assay}_delta_log2_vs_IgG"]
                base[f"promoter_{assay}_delta_log2_vs_IgG"] = promoter[f"{assay}_delta_log2_vs_IgG"]
            tss_values = [tss_by_feature[x][f"{assay}_log2FC"] for x in tss_ids if x in tss_by_feature]
            base[f"TSS_{assay}_median_log2FC"] = median(tss_values)
        output.append(base)
    return output


def state_distribution(pairs: list[dict]):
    output = []
    for enhancer_state in ("strength", "weaken"):
        subset = [row for row in pairs if row["enhancer_state"] == enhancer_state]
        all_promoters = {row["promoter_id"] for row in subset}
        for promoter_state in ("strength", "stable", "weaken"):
            current = [row for row in subset if row["promoter_state"] == promoter_state]
            promoters = {row["promoter_id"] for row in current}
            output.append({
                "enhancer_state": enhancer_state, "promoter_state": promoter_state,
                "unique_promoters": len(promoters),
                "total_unique_promoters": len(all_promoters),
                "proportion_unique_promoters": len(promoters) / len(all_promoters) if all_promoters else "NA",
                "unique_enhancers": len({row["enhancer_id"] for row in current}),
                "unique_genes": len({row["gene_id"] for row in current}),
                "unique_enhancer_promoter_gene_rows": len(current),
            })
    return output


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        while end < len(order) and values[order[end]] == values[order[position]]:
            end += 1
        rank = (position + end - 1) / 2 + 1
        for index in order[position:end]:
            result[index] = rank
        position = end
    return result


def pearson(x, y):
    if len(x) < 3:
        return "NA"
    mean_x, mean_y = statistics.fmean(x), statistics.fmean(y)
    numerator = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - mean_x) ** 2 for a in x) * sum((b - mean_y) ** 2 for b in y))
    return numerator / denominator if denominator else "NA"


def pair_concordance(pair_rows: list[dict]):
    unique = {}
    for row in pair_rows:
        unique[(row["enhancer_state"], row["enhancer_id"], row["promoter_id"])] = row
    output = []
    for enhancer_state in ("strength", "weaken"):
        subset = [row for key, row in unique.items() if key[0] == enhancer_state]
        for assay in ASSAYS:
            x = [float(row[f"enhancer_{assay}_log2FC"]) for row in subset]
            y = [float(row[f"promoter_{assay}_log2FC"]) for row in subset]
            n = len(x)
            same = sum((a > 0 and b > 0) or (a < 0 and b < 0) for a, b in zip(x, y))
            output.append({
                "enhancer_state": enhancer_state, "assay": assay, "unique_enhancer_promoter_pairs": n,
                "median_enhancer_log2FC": median(x), "median_promoter_log2FC": median(y),
                "same_direction_fraction": same / n if n else "NA",
                "both_up_fraction": sum(a > 0 and b > 0 for a, b in zip(x, y)) / n if n else "NA",
                "both_down_fraction": sum(a < 0 and b < 0 for a, b in zip(x, y)) / n if n else "NA",
                "pearson_r": pearson(x, y),
                "spearman_r": pearson(ranks(x), ranks(y)) if n >= 3 else "NA",
            })
    return output


def linked_promoter_mark_summary(pair_rows: list[dict], promoter_by_id: dict):
    output = []
    for enhancer_state in ("strength", "weaken"):
        for promoter_state in ("strength", "stable", "weaken"):
            promoter_ids = {
                row["promoter_id"] for row in pair_rows
                if row["enhancer_state"] == enhancer_state and row["promoter_state"] == promoter_state
            }
            promoters = [promoter_by_id[x] for x in sorted(promoter_ids)]
            for assay in ASSAYS:
                changes = [float(row[f"{assay}_log2FC"]) for row in promoters]
                output.append({
                    "enhancer_state": enhancer_state, "promoter_state": promoter_state,
                    "assay": assay, "unique_promoters": len(promoters),
                    "median_promoter_log2FC": median(changes),
                    "mean_promoter_log2FC": mean(changes),
                    "fraction_promoter_log2FC_gt_0": (
                        sum(value > 0 for value in changes) / len(changes) if changes else "NA"
                    ),
                    "fraction_promoter_log2FC_lt_0": (
                        sum(value < 0 for value in changes) / len(changes) if changes else "NA"
                    ),
                    "median_promoter_delta_log2_vs_IgG": (
                        median(row[f"{assay}_delta_log2_vs_IgG"] for row in promoters)
                        if assay in CHIP_ASSAYS else "NA"
                    ),
                })
    return output


def plot_heatmap(summary: list[dict], output: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups = [
        ("enhancer", "strength"), ("enhancer", "stable"), ("enhancer", "weaken"),
        ("promoter", "strength"), ("promoter", "stable"), ("promoter", "weaken"),
    ]
    lookup = {(row["region_type"], row["state"], row["assay"]): row for row in summary}
    matrix = [[float(lookup[(kind, state, assay)]["median_log2FC"]) for assay in ASSAYS]
              for kind, state in groups]
    limit = max(0.5, max(abs(value) for row in matrix for value in row))
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(ASSAYS)), ASSAYS, rotation=35, ha="right")
    ax.set_yticks(range(len(groups)), [f"{kind.capitalize()} {state.capitalize()}" for kind, state in groups])
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=9,
                    color="white" if abs(value) > limit * 0.55 else "black")
    ax.set_title("Median mark change (6 h / 0 h)", fontsize=14, fontweight="bold")
    cbar = fig.colorbar(image, ax=ax, shrink=0.85)
    cbar.set_label(r"median $\log_2$FC")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"Figure_1_region_mark_log2FC_heatmap.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_state_distribution(rows: list[dict], output: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    states = ("strength", "stable", "weaken")
    colors = {"strength": "#D64B43", "stable": "#BDBDBD", "weaken": "#3E88AE"}
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    for x, enhancer_state in enumerate(("strength", "weaken")):
        bottom = 0.0
        for promoter_state in states:
            row = next(r for r in rows if r["enhancer_state"] == enhancer_state and r["promoter_state"] == promoter_state)
            value = float(row["proportion_unique_promoters"])
            ax.bar(x, value, bottom=bottom, color=colors[promoter_state], width=0.62,
                   label=promoter_state.capitalize() if x == 0 else None)
            if value >= 0.05:
                ax.text(x, bottom + value / 2, f"{100*value:.1f}%", ha="center", va="center", fontsize=10)
            bottom += value
    ax.set_xticks((0, 1), ("Strength enhancers", "Weaken enhancers"))
    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion of unique linked P1 promoters")
    ax.set_title("Promoter states linked by HiChIP", fontsize=14, fontweight="bold")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"Figure_2_linked_promoter_state_distribution.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Submit this real-data analysis as a Linux Slurm compute job.")
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    if float(cfg["fold_threshold"]) != 1.5:
        raise ValueError("This analysis is fixed to the final 1.5-fold state definition")

    project = Path(cfg["project"])
    enhancer_run = Path(cfg["source_enhancer_run"])
    promoter_run = Path(cfg["source_promoter_run"])
    hichip_run = Path(cfg["source_HiChIP_target_run"])
    enhancer_matrix = enhancer_run / "03_CRE_signal_matrix.tsv"
    promoter_matrix = promoter_run / "04_P1_promoter_signal_and_state_matrix.tsv"
    promoter_associations = promoter_run / "05_P1_promoter_state_gene_associations.tsv"
    loop_paths = {
        state: hichip_run / "targets/P1_full_H3_peak_support" / state / "05_all_CRE_loop_gene_links.tsv.gz"
        for state in ("strength", "weaken")
    }
    bedtools = Path(cfg["bedtools_env"]) / "bin/bedtools"
    bam_paths = {label: project / "02.bam" / name for label, name in cfg["bam_names"].items()}
    required = [enhancer_matrix, promoter_matrix, promoter_associations, bedtools,
                *loop_paths.values(), *bam_paths.values()]
    for path in required:
        if not path.is_absolute() or not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    for bam in bam_paths.values():
        if not Path(str(bam) + ".bai").is_file():
            raise FileNotFoundError(Path(str(bam) + ".bai"))

    import pysam
    bam_maps, lengths = {}, None
    for label, path in bam_paths.items():
        with pysam.AlignmentFile(str(path), "rb") as bam:
            if not bam.has_index():
                raise ValueError(f"BAM is not indexed: {path}")
            mapping = {canon(name): name for name in bam.references if canon(name) in MAIN}
            current_lengths = {chrom: bam.get_reference_length(name) for chrom, name in mapping.items()}
            if set(mapping) != MAIN:
                raise ValueError(f"Main chromosome mismatch: {label}")
            if lengths is None:
                lengths = current_lengths
            elif current_lengths != lengths:
                raise ValueError(f"Assembly mismatch: {label}")
            bam_maps[label] = mapping
    first_map = bam_maps[next(iter(bam_maps))]
    if any(mapping != first_map for mapping in bam_maps.values()):
        raise ValueError("BAM chromosome naming differs; multicov needs one shared naming scheme")

    output = Path(cfg["results"]) / f"run_{os.environ['SLURM_JOB_ID']}"
    output.mkdir(parents=True, exist_ok=False)
    for name in ("code", "counts", "logs"):
        (output / name).mkdir()
    command_log = output / "logs/commands.jsonl"
    for name in ("run.py", "config.json", "README.md", "run.slurm", "submit.sh", "check_env.sh", "test_small.py"):
        shutil.copy2(Path(__file__).with_name(name), output / "code" / name)
    print(f"results={output}", flush=True)

    fold, minimum = float(cfg["fold_threshold"]), int(cfg["min_raw_tags_active_state"])
    enhancers = read_regions(enhancer_matrix, "enhancer", fold, minimum)
    promoters = read_regions(promoter_matrix, "promoter", fold, minimum, validate_archived_state=True)
    tss_windows, gene_to_tss = read_tss_windows(
        promoter_associations, int(cfg["promoter_TSS_half_width"]), lengths
    )
    all_features = [*enhancers, *promoters, *tss_windows]
    feature_ids = [row["feature_id"] for row in all_features]
    if len(feature_ids) != len(set(feature_ids)):
        raise ValueError("Combined feature IDs are not unique")
    write_bed(output / "counts/all_regions.bed", all_features, first_map)
    chromosome_rows = [{"feature_id": f"CHR|{chrom}", "chrom": chrom, "start": 0,
                        "end": lengths[chrom]} for chrom in sorted(MAIN, key=chrom_key)]
    write_bed(output / "counts/main_chromosomes.bed", chromosome_rows, first_map)

    labels = list(bam_paths)
    libraries = multicov(bedtools, output / "counts/main_chromosomes.bed",
                         [bam_paths[x] for x in labels], labels,
                         output / "counts/library_sizes.tsv", command_log, int(cfg["min_mapq"]))
    totals = {label: sum(row[label] for row in libraries.values()) for label in labels}
    if any(value <= 0 for value in totals.values()):
        raise ValueError(f"Invalid library totals: {totals}")
    counts = multicov(bedtools, output / "counts/all_regions.bed",
                      [bam_paths[x] for x in labels], labels,
                      output / "counts/all_region_counts.tsv", command_log, int(cfg["min_mapq"]))
    if set(counts) != set(feature_ids):
        raise ValueError("multicov feature IDs differ from the requested region set")
    add_signals(all_features, counts, totals, float(cfg["normalization_total_tags"]),
                float(cfg["normalized_pseudocount"]))

    write_tsv(output / "01_enhancer_mark_signal_matrix.tsv", enhancers)
    write_tsv(output / "02_P1_promoter_mark_signal_matrix.tsv", promoters)
    write_tsv(output / "03_gene_TSS_pm2kb_mark_signal_matrix.tsv", tss_windows)
    recount_qc = source_recount_qc([*enhancers, *promoters])
    write_tsv(output / "03_source_vs_recount_QC.tsv", recount_qc)
    if any(
        row["max_absolute_log2FC_difference"] != "NA"
        and float(row["max_absolute_log2FC_difference"]) > 1e-6
        for row in recount_qc
    ):
        raise ValueError("Recounted H3K27ac/ATAC log2FC differs from the archived state matrices")
    write_tsv(output / "04_region_mark_long.tsv", long_mark_table(all_features))
    mark_summary = summarize_marks([*enhancers, *promoters])
    write_tsv(output / "04_region_mark_summary.tsv", mark_summary)

    enhancer_by_id = {row["region_id"]: row for row in enhancers}
    promoter_by_id = {row["region_id"]: row for row in promoters}
    tss_by_feature = {row["feature_id"]: row for row in tss_windows}
    pairs, mapping_failures = map_loop_links(
        loop_paths, enhancer_by_id, promoter_associations, promoter_by_id,
        float(cfg["loop_q_max"])
    )
    pair_rows = add_pair_signals(pairs, enhancer_by_id, promoter_by_id,
                                 tss_by_feature, gene_to_tss)
    write_tsv(output / "05_HiChIP_enhancer_promoter_gene_mark_matrix.tsv", pair_rows)
    distribution = state_distribution(pair_rows)
    write_tsv(output / "06_linked_promoter_state_distribution.tsv", distribution)
    concordance = pair_concordance(pair_rows)
    write_tsv(output / "07_linked_enhancer_promoter_mark_concordance.tsv", concordance)
    linked_promoter_summary = linked_promoter_mark_summary(pair_rows, promoter_by_id)
    write_tsv(output / "08_linked_promoter_mark_summary.tsv", linked_promoter_summary)

    plot_heatmap(mark_summary, output)
    plot_state_distribution(distribution, output)
    metadata = {
        "completed": datetime.now(timezone.utc).isoformat(),
        "mode": "pooled_descriptive_no_biological_replicate_P_or_FDR",
        "region_definitions": {
            "enhancer": "final_v27_full_H3K27ac_consensus_CRE",
            "promoter": "final_P1_full_H3_peak_support_consensus",
            "TSS": f"gene_TSS_plus_minus_{cfg['promoter_TSS_half_width']}_bp",
        },
        "normalization": "proper-pair main-chromosome alignments per 10M; no length division; no median centering",
        "ChIP_secondary_metric": "log2(mark_per10M/IgG_per10M) at each time and 6h-minus-0h delta",
        "counts": {
            "enhancers": len(enhancers), "promoters": len(promoters), "TSS_windows": len(tss_windows),
            "linked_enhancer_promoter_gene_rows": len(pair_rows),
            "loop_mapping_failures": dict(mapping_failures),
            "enhancer_states": dict(Counter(row["state"] for row in enhancers)),
            "promoter_states": dict(Counter(row["state"] for row in promoters)),
        },
        "library_totals": totals,
        "input_BAM_stats": {label: {"path": str(path), "bytes": path.stat().st_size,
                                      "mtime_ns": path.stat().st_mtime_ns}
                            for label, path in bam_paths.items()},
        "input_sha256": {str(path): sha256(path) for path in
                          [enhancer_matrix, promoter_matrix, promoter_associations, *loop_paths.values()]},
        "config": cfg,
    }
    (output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (output / "completed.json").write_text(
        json.dumps({"completed": True, "results": str(output)}, indent=2), encoding="utf-8"
    )
    print(f"completed=1\nresults={output}", flush=True)


if __name__ == "__main__":
    main()
