#!/usr/bin/env python3
"""Motif-only weaken-definition sensitivity analysis reusing v19 region signals."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


MAIN = {str(x) for x in range(1, 19)} | {"X", "Y"}
LIPID = re.compile(r"PPAR|RXR|SREBF|HNF4|ESRR|NR1H|LXR|KLF15|MLXIPL", re.I)


def canon(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


def number(value, default=math.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def write_tsv(path: Path, rows: list[dict], fields=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_region_table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "chrom", "region_id", "supported_0h", "motif_start", "motif_end",
        "motif_sequence_ok", "motif_length_bp", "H3K27ac_0h_raw_tags",
        "ATAC_0h_raw_tags", "H3K27ac_log2FC", "ATAC_log2FC",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Missing required columns in {path}: {sorted(required - set(rows[0] if rows else []))}")
    ids = [row["region_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate region_id in source matrix")
    return rows


def qualifies(row: dict, definition: dict, min_tags: float) -> bool:
    if int(row["supported_0h"]) != 1 or int(row["motif_sequence_ok"]) != 1:
        return False
    if number(row["H3K27ac_0h_raw_tags"]) < min_tags or number(row["ATAC_0h_raw_tags"]) < min_tags:
        return False
    if number(row["H3K27ac_log2FC"]) > -math.log2(float(definition["h3_fold"])):
        return False
    atac = number(row["ATAC_log2FC"])
    rule = definition["atac_rule"]
    if rule == "down":
        return atac < 0
    if rule == "not_up":
        return atac <= math.log2(float(definition["atac_max_fold"]))
    if rule == "ignore":
        return True
    raise ValueError(f"Unknown ATAC rule: {rule}")


def attach_sequences(rows: list[dict], fasta_path: Path) -> dict[str, str]:
    import pysam
    genome = pysam.FastaFile(str(fasta_path))
    try:
        fa_map = {}
        for name in genome.references:
            key = canon(name)
            if key in MAIN:
                if key in fa_map:
                    raise ValueError(f"Ambiguous FASTA chromosome alias: {name}")
                fa_map[key] = name
        for row in rows:
            if int(row["motif_sequence_ok"]) != 1:
                continue
            chrom = canon(row["chrom"])
            start, end = int(row["motif_start"]), int(row["motif_end"])
            seq = genome.fetch(fa_map[chrom], start, end).upper()
            acgt = sum(seq.count(base) for base in "ACGT")
            if len(seq) != int(row["motif_length_bp"]) or not acgt:
                raise ValueError(f"Invalid motif sequence: {row['region_id']}")
            row["_sequence"] = seq
            row["_gc"] = (seq.count("G") + seq.count("C")) / acgt
    finally:
        genome.close()
    return fa_map


def match_controls(foreground: list[dict], controls: list[dict], max_delta: float, seed: int) -> list[tuple[dict, dict]]:
    """Closest-GC 1:1 matching without replacement for fixed-length windows."""
    rng = random.Random(seed)
    buckets: dict[float, list[dict]] = {}
    for row in controls:
        key = round(float(row["_gc"]), 6)
        buckets.setdefault(key, []).append(row)
    for bucket in buckets.values():
        rng.shuffle(bucket)
    keys = sorted(buckets)
    order = list(foreground)
    rng.shuffle(order)
    pairs = []
    for row in order:
        target = float(row["_gc"])
        pos = bisect.bisect_left(keys, target)
        candidates = []
        left, right = pos - 1, pos
        while left >= 0 and target - keys[left] <= max_delta:
            candidates.append(keys[left]); left -= 1
        while right < len(keys) and keys[right] - target <= max_delta:
            candidates.append(keys[right]); right += 1
        if not candidates:
            continue
        chosen_key = min(candidates, key=lambda value: (abs(value - target), value))
        control = buckets[chosen_key].pop()
        if not buckets[chosen_key]:
            del buckets[chosen_key]
            keys.pop(bisect.bisect_left(keys, chosen_key))
        pairs.append((row, control))
    return pairs


def read_catalog(path: Path) -> list[tuple[str, str]]:
    records, seen = [], set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("MOTIF "):
                continue
            fields = line.strip().split(maxsplit=2)
            motif_id, alt_id = fields[1], fields[2] if len(fields) > 2 else "NA"
            if motif_id in seen:
                raise ValueError(f"Duplicate motif ID: {motif_id}")
            seen.add(motif_id); records.append((motif_id, alt_id))
    if not records:
        raise ValueError("Empty MEME motif catalog")
    return records


def norm_header(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def field(row: dict, *aliases: str, default=""):
    normalized = {norm_header(key): value for key, value in row.items() if key is not None}
    for alias in aliases:
        if norm_header(alias) in normalized:
            return normalized[norm_header(alias)]
    return default


def read_ame(path: Path, definition: dict, catalog: list[tuple[str, str]]) -> tuple[list[dict], dict]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if "# AME (Analysis of Motif Enrichment):" not in text or "\n# ame " not in text:
        raise ValueError(f"AME output lacks completion footer: {path}")
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("#")]
    parsed = []
    for raw in csv.DictReader(lines, delimiter="\t"):
        motif_id = field(raw, "motif_ID", "motif id", "motif")
        alt_id = field(raw, "motif_alt_ID", "motif alt id", "motif_name") or "NA"
        tp, fp = number(field(raw, "tp%", default="nan")), number(field(raw, "fp%", default="nan"))
        parsed.append({
            "definition_id": definition["id"], "definition_role": definition["role"],
            "h3_fold_threshold": definition["h3_fold"], "atac_rule": definition["atac_rule"],
            "motif_id": motif_id, "motif_alt_id": alt_id,
            "lipid_focus": int(bool(LIPID.search(f"{motif_id} {alt_id}"))),
            "AME_raw_optimal_p": number(field(raw, "p-value", "pvalue", default="nan")),
            "AME_threshold_adjusted_p": number(field(raw, "adj_p-value", "q-value", default="nan")),
            "AME_E_value": number(field(raw, "E-value", default="nan")),
            "positive_hit_percent": tp, "control_hit_percent": fp,
            "enrichment_ratio": tp / fp if math.isfinite(tp) and math.isfinite(fp) and fp > 0 else "NA",
            "AME_reported": 1, "AME_missing_reason": "", "AME_source": str(path),
        })
    for row in parsed:
        p_value = row["AME_threshold_adjusted_p"]
        if not isinstance(p_value, float) or not math.isfinite(p_value) or not 0 <= p_value <= 1:
            raise ValueError(f"Invalid AME threshold-adjusted P in {path}: {p_value}")
    expected = {motif_id for motif_id, _ in catalog}
    observed_ids = [row["motif_id"] for row in parsed]
    if len(observed_ids) != len(set(observed_ids)) or set(observed_ids) - expected:
        raise ValueError(f"Duplicate or unexpected motif IDs in {path}")
    if len(observed_ids) < int(0.95 * len(catalog)):
        raise ValueError(f"AME reported only {len(observed_ids)}/{len(catalog)} motifs")
    observed = set(observed_ids)
    for motif_id, alt_id in catalog:
        if motif_id not in observed:
            parsed.append({
                "definition_id": definition["id"], "definition_role": definition["role"],
                "h3_fold_threshold": definition["h3_fold"], "atac_rule": definition["atac_rule"],
                "motif_id": motif_id, "motif_alt_id": alt_id,
                "lipid_focus": int(bool(LIPID.search(f"{motif_id} {alt_id}"))),
                "AME_raw_optimal_p": 1.0, "AME_threshold_adjusted_p": 1.0, "AME_E_value": "NA",
                "positive_hit_percent": "NA", "control_hit_percent": "NA", "enrichment_ratio": "NA",
                "AME_reported": 0,
                "AME_missing_reason": "not_reported_by_completed_AME_run_neutral_P1",
                "AME_source": str(path),
            })
    order = {motif_id: index for index, (motif_id, _) in enumerate(catalog)}
    parsed.sort(key=lambda row: order[row["motif_id"]])
    return parsed, {"expected": len(catalog), "reported": len(observed_ids), "neutral_filled": len(catalog)-len(observed_ids)}


def bh_adjust(rows: list[dict], p_field: str, q_field: str) -> None:
    indexed = sorted(((index, float(row[p_field])) for index, row in enumerate(rows)), key=lambda item: item[1])
    total, running = len(indexed), 1.0
    adjusted = [1.0] * total
    for position in range(total - 1, -1, -1):
        _, p_value = indexed[position]
        running = min(running, p_value * total / (position + 1))
        adjusted[position] = min(1.0, running)
    for position, (row_index, _) in enumerate(indexed):
        rows[row_index][q_field] = adjusted[position]


def run_ame_task(executable: Path, database: Path, directory: Path) -> Path:
    output, log = directory / "AME" / "ame.tsv", directory / "AME" / "ame.stderr.log"
    output.parent.mkdir(parents=True, exist_ok=False)
    command = [str(executable), "--verbose", "1", "--text", "--control", str(directory / "control.fa"),
               "--method", "fisher", "--scoring", "max", "--noseq",
               "--evalue-report-threshold", "1000000000000",
               str(directory / "foreground.fa"), str(database)]
    env = dict(os.environ)
    env["PATH"] = str(executable.parent) + os.pathsep + env.get("PATH", "")
    with output.open("w", encoding="utf-8") as stdout, log.open("w", encoding="utf-8") as stderr:
        subprocess.run(command, stdout=stdout, stderr=stderr, check=True, env=env)
    return output


def write_inputs(directory: Path, foreground: list[dict], pairs: list[tuple[dict, dict]], fa_map: dict[str, str]) -> None:
    directory.mkdir(parents=True, exist_ok=False)
    selected = {"foreground": foreground, "control": [control for _, control in pairs]}
    for label, rows in selected.items():
        with (directory / f"{label}.bed").open("w", encoding="utf-8") as bed, (directory / f"{label}.fa").open("w", encoding="utf-8") as fasta:
            for row in rows:
                chrom = fa_map[canon(row["chrom"])]
                bed.write(f"{chrom}\t{row['motif_start']}\t{row['motif_end']}\t{row['region_id']}\n")
                fasta.write(f">{row['region_id']}\n{row['_sequence']}\n")
