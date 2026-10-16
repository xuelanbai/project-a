#!/usr/bin/env python3
"""P1 promoter motif enrichment using the final enhancer ATAC-core standard.

This downstream task deliberately reuses the archived ATAC consensus cores from
the final enhancer analysis.  That guarantees that promoter and enhancer motif
results use identical ATAC peak grouping and representative summits.  P1
promoter states are not recalculated; their archived state matrix is validated
against the fixed 1.5-fold H3K27ac rule and the explicit >=5 active-state read
filter before motif analysis.
"""
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
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone


MAIN = {str(x) for x in range(1, 19)} | {"X", "Y"}
STATES = ("strength", "weaken")
COLORS = {"strength": "#D64B43", "weaken": "#3E88AE"}


def canon(chrom: str) -> str:
    value = str(chrom).strip()
    return value[3:] if value.startswith("chr") else value


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


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class NonOverlappingRegionIndex:
    """Point lookup for an explicitly non-overlapping interval collection."""

    def __init__(self, rows: list[dict]):
        grouped: dict[str, list[dict]] = defaultdict(list)
        for row in rows:
            grouped[row["chrom"]].append(row)
        self.data = {}
        for chrom, values in grouped.items():
            ordered = sorted(values, key=lambda x: (x["start"], x["end"], x["region_id"]))
            for previous, current in zip(ordered, ordered[1:]):
                if current["start"] < previous["end"]:
                    raise ValueError(
                        f"P1 consensus promoters overlap: {previous['region_id']} and {current['region_id']}"
                    )
            self.data[chrom] = ([row["start"] for row in ordered], ordered)

    def contains(self, chrom: str, point: int) -> dict | None:
        starts, rows = self.data.get(canon(chrom), ([], []))
        index = bisect.bisect_right(starts, int(point)) - 1
        if index >= 0 and rows[index]["start"] <= point < rows[index]["end"]:
            return rows[index]
        return None


class IntervalIndex:
    def __init__(self, intervals):
        grouped = defaultdict(list)
        for chrom, start, end in intervals:
            grouped[canon(chrom)].append((int(start), int(end)))
        self.data = {}
        for chrom, values in grouped.items():
            merged = []
            for start, end in sorted(set(values)):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self.data[chrom] = ([x[0] for x in merged], [x[1] for x in merged])

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        starts, ends = self.data.get(canon(chrom), ([], []))
        index = bisect.bisect_left(starts, end) - 1
        return index >= 0 and ends[index] > start


def as_int(row: dict, field: str) -> int:
    try:
        return int(float(row[field]))
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid integer in {field}: {row.get(field)!r}") from error


def as_float(row: dict, field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid number in {field}: {row.get(field)!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"Non-finite number in {field}: {row.get(field)!r}")
    return value


def derive_state(row: dict, fold: float, minimum: int) -> str:
    value = as_float(row, "H3K27ac_log2FC")
    cutoff = math.log2(fold)
    if (
        value >= cutoff
        and as_int(row, "supported_6h") == 1
        and as_float(row, "H3K27ac_6h_raw_tags") >= minimum
        and as_float(row, "ATAC_6h_raw_tags") >= minimum
    ):
        return "strength"
    if (
        value <= -cutoff
        and as_int(row, "supported_0h") == 1
        and as_float(row, "H3K27ac_0h_raw_tags") >= minimum
        and as_float(row, "ATAC_0h_raw_tags") >= minimum
    ):
        return "weaken"
    return "stable"


def active_eligible(row: dict, state: str, minimum: int) -> bool:
    time = "6h" if state == "strength" else "0h"
    return as_int(row, f"supported_{time}") == 1 and all(
        as_float(row, f"{assay}_{time}_raw_tags") >= minimum for assay in ("H3K27ac", "ATAC")
    )


def read_promoters(path: Path, fold: float, minimum: int) -> tuple[list[dict], list[dict]]:
    required = {
        "region_id", "chrom", "start", "end", "supported_0h", "supported_6h",
        "H3K27ac_0h_raw_tags", "H3K27ac_6h_raw_tags", "ATAC_0h_raw_tags",
        "ATAC_6h_raw_tags", "H3K27ac_log2FC", "ATAC_log2FC", "promoter_state",
    }
    rows = read_tsv(path)
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Unexpected P1 state-matrix schema: {path}")
    result, seen, qc = [], set(), []
    mismatches = []
    for raw in rows:
        region_id = raw["region_id"]
        if region_id in seen:
            raise ValueError(f"Duplicate promoter region_id: {region_id}")
        seen.add(region_id)
        chrom, start, end = canon(raw["chrom"]), as_int(raw, "start"), as_int(raw, "end")
        if chrom not in MAIN or not 0 <= start < end:
            raise ValueError(f"Invalid promoter interval: {chrom}:{start}-{end}")
        row = dict(raw, chrom=chrom, start=start, end=end)
        expected = derive_state(row, fold, minimum)
        if raw["promoter_state"] != expected:
            mismatches.append((region_id, raw["promoter_state"], expected))
        row["_derived_state"] = expected
        result.append(row)
    if mismatches:
        preview = "; ".join(":".join(x) for x in mismatches[:10])
        raise ValueError(f"Archived promoter states disagree with fixed rules ({len(mismatches)}): {preview}")
    for state in ("strength", "weaken", "stable"):
        subset = [row for row in result if row["promoter_state"] == state]
        qc.append({
            "promoter_state": state,
            "promoters": len(subset),
            "min_ATAC_0h_raw_tags": min((as_float(x, "ATAC_0h_raw_tags") for x in subset), default="NA"),
            "min_ATAC_6h_raw_tags": min((as_float(x, "ATAC_6h_raw_tags") for x in subset), default="NA"),
            "ATAC_below_5_either_time": sum(
                as_float(x, "ATAC_0h_raw_tags") < minimum or as_float(x, "ATAC_6h_raw_tags") < minimum
                for x in subset
            ),
        })
    return result, qc


def read_cores(path: Path) -> list[dict]:
    required = {
        "core_id", "chrom", "start", "end", "representative_summit",
        "representative_source_time", "summit_start", "summit_end",
    }
    rows = read_tsv(path)
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Unexpected ATAC-core schema: {path}")
    result, seen = [], set()
    for raw in rows:
        core_id = raw["core_id"]
        if core_id in seen:
            raise ValueError(f"Duplicate ATAC core_id: {core_id}")
        seen.add(core_id)
        chrom = canon(raw["chrom"])
        start, end = as_int(raw, "start"), as_int(raw, "end")
        summit = as_int(raw, "representative_summit")
        if chrom not in MAIN or not 0 <= start <= summit < end:
            raise ValueError(f"Invalid ATAC core/summit: {core_id}")
        result.append(dict(raw, chrom=chrom, start=start, end=end,
                           representative_summit=summit, center=summit))
    return result


def read_catalog(path: Path):
    result, seen = [], set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MOTIF "):
                parts = line.strip().split(maxsplit=2)
                motif_id = parts[1]
                if motif_id in seen:
                    raise ValueError(f"Duplicate motif ID in database: {motif_id}")
                seen.add(motif_id)
                result.append((motif_id, parts[2] if len(parts) > 2 else "NA"))
    if not result:
        raise ValueError(f"No motifs in {path}")
    return result


def norm_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def get_field(row: dict, *names, default=""):
    normalized = {norm_header(key): value for key, value in row.items() if key is not None}
    for name in names:
        if norm_header(name) in normalized:
            return normalized[norm_header(name)]
    return default


def number(value, default=math.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_ame(path: Path, contrast: str, catalog):
    lines = [line for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
             if line.strip() and not line.startswith("#")]
    rows = []
    for raw in csv.DictReader(lines, delimiter="\t"):
        motif_id = get_field(raw, "motif_ID", "motif id", "motif")
        alt_id = get_field(raw, "motif_alt_ID", "motif alt id", "motif_name") or "NA"
        tp, fp = number(get_field(raw, "tp%")), number(get_field(raw, "fp%"))
        rows.append({
            "contrast": contrast,
            "motif_id": motif_id,
            "motif_alt_id": alt_id,
            "AME_raw_optimal_p": number(get_field(raw, "p-value", "pvalue")),
            "AME_threshold_adjusted_p": number(get_field(raw, "adj_p-value", "q-value")),
            "AME_E_value": number(get_field(raw, "E-value")),
            "foreground_hit_percent": tp,
            "control_hit_percent": fp,
            "enrichment_ratio": tp / fp if math.isfinite(tp) and math.isfinite(fp) and fp > 0 else "NA",
        })
    expected = {x[0] for x in catalog}
    observed = {x["motif_id"] for x in rows}
    if not rows or len(observed) != len(rows) or observed - expected:
        raise ValueError(f"Invalid, duplicate, or unexpected AME results: {path}")
    if len(observed) < 0.95 * len(expected):
        raise ValueError(f"AME reported only {len(observed)}/{len(expected)} motifs")
    for motif_id, alt_id in catalog:
        if motif_id not in observed:
            rows.append({
                "contrast": contrast, "motif_id": motif_id, "motif_alt_id": alt_id,
                "AME_raw_optimal_p": 1.0, "AME_threshold_adjusted_p": 1.0,
                "AME_E_value": "NA", "foreground_hit_percent": "NA",
                "control_hit_percent": "NA", "enrichment_ratio": "NA",
            })
    return rows


def read_homer(path: Path, contrast: str):
    with path.open(encoding="utf-8", errors="replace") as handle:
        raw = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for row in raw:
        motif = row.get("Motif Name")
        if not motif:
            raise ValueError(f"Unexpected HOMER table header: {path}")
        target = next(number(v.rstrip("%")) for k, v in row.items() if k.startswith("% of Target"))
        background = next(number(v.rstrip("%")) for k, v in row.items() if k.startswith("% of Background"))
        native_q = next(number(v) for k, v in row.items() if k.lower().startswith("q-value"))
        rows.append({
            "contrast": contrast, "motif": motif, "p": number(row.get("P-value")),
            "native_q": native_q, "foreground_percent": target,
            "background_percent": background,
            "enrichment_ratio": target / background if background > 0 else "NA",
        })
    if not rows:
        raise ValueError(f"Empty HOMER results: {path}")
    return rows


def bh_adjust(rows: list[dict], p_field: str, q_field: str):
    indexed = sorted(((i, float(row[p_field])) for i, row in enumerate(rows)), key=lambda x: x[1])
    total, running = len(indexed), 1.0
    adjusted = [1.0] * total
    for position in range(total - 1, -1, -1):
        running = min(running, indexed[position][1] * total / (position + 1))
        adjusted[position] = min(1.0, running)
    for position, (index, _) in enumerate(indexed):
        rows[index][q_field] = adjusted[position]


def motif_family(text: str) -> str:
    value = text.upper()
    rules = [
        (r"CTCF|BORIS", "CTCF"),
        (r"\b(JUN|FOS|BATF|JDP)", "AP-1"),
        (r"NFKB|NF.KB|\bREL[AB]?\b", "NF-kB"),
        (r"CEBP", "C/EBP"),
        (r"PPAR|RXR", "PPAR/RXR"),
        (r"SREBF", "SREBF"),
        (r"ATF4", "ISR/ATF4"),
        (r"\bKLF", "KLF"),
    ]
    for pattern, family in rules:
        if re.search(pattern, value):
            return family
    token = re.split(r"::|[/_(\s]", text.strip())[0]
    return token or "unclassified"


def representative_ame(rows: list[dict]):
    groups = {}
    for row in rows:
        family = motif_family(f"{row['motif_alt_id']} {row['motif_id']}")
        candidate = dict(row, motif_family=family)
        ratio = candidate["enrichment_ratio"] if isinstance(candidate["enrichment_ratio"], (int, float)) else -math.inf
        key = (candidate["contrast"], family)
        rank = (
            float(candidate["AME_BH_q_within_contrast"]), -ratio,
            float(candidate["AME_threshold_adjusted_p"]), candidate["motif_id"],
        )
        if key not in groups or rank < groups[key][0]:
            groups[key] = (rank, candidate)
    return [value[1] for _, value in sorted(groups.items())]


def representative_homer(rows: list[dict]):
    groups = {}
    for row in rows:
        family = motif_family(row["motif"])
        candidate = dict(row, motif_family=family)
        ratio = candidate["enrichment_ratio"] if isinstance(candidate["enrichment_ratio"], (int, float)) else -math.inf
        key = (candidate["contrast"], family)
        rank = (float(candidate["HOMER_BH_q_within_contrast"]), -ratio,
                float(candidate["p"]), candidate["motif"])
        if key not in groups or rank < groups[key][0]:
            groups[key] = (rank, candidate)
    return [value[1] for _, value in sorted(groups.items())]


def overlap_window_index(rows, width):
    grouped = defaultdict(list)
    half = width // 2
    for row in rows:
        grouped[row["chrom"]].append((row["center"] - half, row["center"] + half))
    return {chrom: IntervalIndex((chrom, start, end) for start, end in values)
            for chrom, values in grouped.items()}


def window_overlaps(index, row, width):
    chrom_index = index.get(row["chrom"])
    if chrom_index is None:
        return False
    half = width // 2
    return chrom_index.overlaps(row["chrom"], row["center"] - half, row["center"] + half)


def match_gc(foreground, controls, max_delta, seed):
    rng = random.Random(seed)
    buckets = defaultdict(list)
    for row in controls:
        buckets[round(row["_gc"], 6)].append(row)
    for values in buckets.values():
        rng.shuffle(values)
    keys = sorted(buckets)
    order = list(foreground)
    rng.shuffle(order)
    pairs, unmatched = [], []
    for row in order:
        target = row["_gc"]
        candidates = [key for key in keys if abs(key - target) <= max_delta + 1e-12]
        if not candidates:
            unmatched.append(row["core_id"])
            continue
        chosen = min(candidates, key=lambda x: (abs(x - target), x))
        control = buckets[chosen].pop()
        if not buckets[chosen]:
            del buckets[chosen]
            keys.remove(chosen)
        pairs.append((row, control))
    return pairs, unmatched


def run_command(command, stdout: Path, stderr: Path, command_log: Path, env=None):
    command = [str(item) for item in command]
    command_log.parent.mkdir(parents=True, exist_ok=True)
    with command_log.open("a", encoding="utf-8") as log:
        log.write(json.dumps({"command": command, "time": datetime.now(timezone.utc).isoformat()}) + "\n")
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        subprocess.run(command, check=True, stdout=out, stderr=err, env=env)


def write_bed(path: Path, rows: list[dict], chrom_map: dict[str, str],
              start_key="start", end_key="end", id_key="region_id"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for row in sorted(rows, key=lambda x: (int(canon(x["chrom"]) not in {str(i) for i in range(1, 19)}),
                                               int(canon(x["chrom"])) if canon(x["chrom"]).isdigit() else canon(x["chrom"]),
                                               int(x[start_key]), int(x[end_key]), str(x[id_key]))):
            chrom = chrom_map[canon(row["chrom"])]
            handle.write(f"{chrom}\t{int(row[start_key])}\t{int(row[end_key])}\t{row[id_key]}\n")


def plot_combined(rows: list[dict], output: Path, tool: str, q_field: str,
                  label_field: str, top_n: int, rank_mode: str, sampling: dict):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = {}
    for state in ("weaken", "strength"):
        contrast = f"promoter_{state}_1p5_500bp_ATAC_core"
        current = [
            row for row in rows
            if row["contrast"] == contrast
            and isinstance(row["enrichment_ratio"], (int, float))
            and float(row["enrichment_ratio"]) > 1
            and float(row[q_field]) <= 0.05
        ]
        if rank_mode == "q_rank":
            current.sort(key=lambda r: (float(r[q_field]), -float(r["enrichment_ratio"]), str(r[label_field])))
        else:
            current.sort(key=lambda r: (-float(r["enrichment_ratio"]), float(r[q_field]), str(r[label_field])))
        selected[state] = current[:top_n]

    heights = [max(2.1, 0.35 * max(1, len(selected[state])) + 0.8) for state in ("weaken", "strength")]
    fig, axes = plt.subplots(2, 1, figsize=(10.5, sum(heights) + 1.2),
                             gridspec_kw={"height_ratios": heights})
    for ax, state in zip(axes, ("weaken", "strength")):
        current = selected[state]
        if current:
            values = [math.log2(float(row["enrichment_ratio"])) for row in current]
            labels = [f"{row['motif_family']}  {row[label_field]}"[:58] for row in current]
            ax.barh(range(len(current)), values, color=COLORS[state], height=0.72)
            ax.set_yticks(range(len(current)), labels, fontsize=10.5)
            ax.invert_yaxis()
        else:
            ax.text(0.02, 0.5, "No significant enriched motif family", transform=ax.transAxes,
                    ha="left", va="center", fontsize=11)
            ax.set_yticks([])
        n = sampling.get(state, {}).get("matched_pairs", "NA")
        ax.text(0.99, 0.98, f"{state.capitalize()}   n={n} windows", color=COLORS[state],
                transform=ax.transAxes, ha="right", va="top", fontsize=11, fontweight="bold")
        ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", labelsize=10)
    axes[-1].set_xlabel(r"$\log_2$(foreground / GC-matched control)", fontsize=12)
    fig.suptitle("TF motif enrichment in P1 promoters", fontsize=17, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    stem = output / f"Figure_{tool}_P1_promoter_{rank_mode}"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Submit this analysis as a Linux Slurm compute job; do not run it on the login node.")
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    if float(cfg["fold_threshold"]) != 1.5 or int(cfg["motif_window_bp"]) != 500:
        raise ValueError("This final version fixes the state threshold at 1.5-fold and motif windows at 500 bp")
    if cfg.get("run_de_novo"):
        raise ValueError("HOMER de novo analysis is intentionally disabled")
    if cfg.get("run_homer") and int(cfg["homer_threads"]) + int(bool(cfg.get("run_ame"))) > int(
        os.environ.get("SLURM_CPUS_PER_TASK", "4")
    ):
        raise ValueError("Allocated CPUs are insufficient for concurrent HOMER and AME")

    core_run = Path(cfg["source_ATAC_core_run"])
    promoter_run = Path(cfg["source_P1_promoter_state_run"])
    core_path = core_run / "02_ATAC_consensus_cores.tsv"
    promoter_path = promoter_run / "04_P1_promoter_signal_and_state_matrix.tsv"
    fasta = Path(cfg["reference_fasta"])
    motif_db = Path(cfg["motif_database"])
    analysis_env, motif_env = Path(cfg["analysis_env"]), Path(cfg["motif_env"])
    executables = {
        "ame": motif_env / "bin/ame",
        "homer": analysis_env / "bin/findMotifs.pl",
    }
    required = [core_path, promoter_path, fasta, Path(str(fasta) + ".fai"), motif_db]
    if cfg.get("run_ame"):
        required.append(executables["ame"])
    if cfg.get("run_homer"):
        required.append(executables["homer"])
    for path in required:
        if not path.is_absolute() or not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    output = Path(cfg["results"]) / f"run_{os.environ['SLURM_JOB_ID']}"
    output.mkdir(parents=True, exist_ok=False)
    for name in ("code", "logs", "regions", "motifs"):
        (output / name).mkdir()
    command_log = output / "logs/commands.jsonl"
    for name in ("run.py", "config.json", "README.md", "run.slurm", "submit.sh", "check_env.sh", "test_small.py"):
        shutil.copy2(Path(__file__).with_name(name), output / "code" / name)
    print(f"results={output}", flush=True)

    fold = float(cfg["fold_threshold"])
    minimum = int(cfg["min_raw_tags_active_state"])
    width = int(cfg["motif_window_bp"])
    promoters, promoter_qc = read_promoters(promoter_path, fold, minimum)
    cores = read_cores(core_path)
    write_tsv(output / "01_P1_state_and_ATAC_filter_QC.tsv", promoter_qc)

    import pysam
    genome = pysam.FastaFile(str(fasta))
    fasta_map = {canon(name): name for name in genome.references if canon(name) in MAIN}
    lengths = {chrom: genome.get_reference_length(name) for chrom, name in fasta_map.items()}
    if set(lengths) != MAIN:
        raise ValueError(f"Reference main chromosomes differ: {sorted(set(lengths) ^ MAIN)}")
    chrom_map = fasta_map

    promoter_index = NonOverlappingRegionIndex(promoters)
    mappings, mapping_excluded = [], Counter()
    for core in cores:
        promoter = promoter_index.contains(core["chrom"], core["representative_summit"])
        if promoter is None:
            mapping_excluded["representative_summit_outside_P1_promoter"] += 1
            continue
        row = {key: value for key, value in core.items() if not key.startswith("_")}
        row.update({
            "region_id": promoter["region_id"],
            "promoter_start": promoter["start"],
            "promoter_end": promoter["end"],
            "promoter_state": promoter["promoter_state"],
            "H3K27ac_log2FC": promoter["H3K27ac_log2FC"],
            "ATAC_log2FC": promoter["ATAC_log2FC"],
            "ATAC_stratum": promoter.get("ATAC_stratum", "NA"),
            "center": core["representative_summit"],
            "cores_per_promoter": 0,
            "_promoter": promoter,
        })
        mappings.append(row)
    per_promoter = Counter(row["region_id"] for row in mappings)
    for row in mappings:
        row["cores_per_promoter"] = per_promoter[row["region_id"]]
    write_tsv(output / "02_ATAC_core_to_P1_promoter.tsv", mappings)

    sequence_exclusions = []
    half = width // 2
    for row in mappings:
        start, end = row["center"] - half, row["center"] + half
        reason = None
        if start < 0 or end > lengths[row["chrom"]]:
            reason = "chromosome_edge"
        else:
            sequence = genome.fetch(fasta_map[row["chrom"]], start, end).upper()
            acgt = sum(sequence.count(base) for base in "ACGT")
            if len(sequence) != width or not acgt or (len(sequence) - acgt) / len(sequence) > float(cfg["max_N_fraction"]):
                reason = "invalid_or_N_rich_sequence"
            else:
                row["_sequence"] = sequence
                row["_gc"] = (sequence.count("G") + sequence.count("C")) / acgt
                row["window_start"], row["window_end"] = start, end
        if reason:
            sequence_exclusions.append({"core_id": row["core_id"], "region_id": row["region_id"], "reason": reason})
    genome.close()
    write_tsv(output / "02_motif_sequence_exclusions.tsv", sequence_exclusions,
              ["core_id", "region_id", "reason"])

    motif_jobs, sampling = [], []
    for state_index, state in enumerate(STATES):
        foreground_promoters = [row for row in promoters if row["promoter_state"] == state]
        foreground_ids = {row["region_id"] for row in foreground_promoters}
        state_mappings = [row for row in mappings if row["region_id"] in foreground_ids]
        valid_fg = [row for row in state_mappings if "_sequence" in row]
        write_bed(output / "regions" / f"{state}_P1_promoters.bed", foreground_promoters, chrom_map)
        write_tsv(output / "regions" / f"{state}_ATAC_core_membership.tsv", state_mappings)
        write_bed(output / "regions" / f"{state}_ATAC_core_500bp.bed", valid_fg, chrom_map,
                  "window_start", "window_end", "core_id")

        foreground_core_ids = {row["core_id"] for row in valid_fg}
        window_index = overlap_window_index(valid_fg, width)
        controls = [
            row for row in mappings
            if "_sequence" in row
            and row["core_id"] not in foreground_core_ids
            and active_eligible(row["_promoter"], state, minimum)
            and not window_overlaps(window_index, row, width)
        ]
        seed = int(cfg["seed"]) + state_index
        pairs, unmatched = match_gc(valid_fg, controls, float(cfg["max_GC_difference"]), seed)
        matched_fraction = len(pairs) / len(valid_fg) if valid_fg else 0.0
        contrast = f"promoter_{state}_1p5_500bp_ATAC_core"
        directory = output / "motifs" / contrast
        directory.mkdir()
        write_tsv(directory / "background_pairs.tsv", [{
            "foreground_id": fg["core_id"], "foreground_promoter_id": fg["region_id"],
            "control_id": bg["core_id"], "control_promoter_id": bg["region_id"],
            "foreground_GC": fg["_gc"], "control_GC": bg["_gc"],
            "delta_GC": abs(fg["_gc"] - bg["_gc"]),
        } for fg, bg in pairs])
        write_tsv(directory / "unmatched_foreground.tsv",
                  [{"core_id": core_id} for core_id in unmatched], ["core_id"])
        for label, position in (("foreground", 0), ("control", 1)):
            with (directory / f"{label}.fa").open("w", encoding="utf-8") as handle:
                for pair in pairs:
                    row = pair[position]
                    handle.write(f">{row['core_id']}|{row['region_id']}\n{row['_sequence']}\n")
        status = "planned"
        if len(pairs) < int(cfg["minimum_motif_pairs"]):
            status = "not_tested_fewer_than_minimum_pairs"
        elif matched_fraction < float(cfg["minimum_matched_fraction"]):
            status = "not_tested_low_matched_fraction"
        sampling_row = {
            "contrast": contrast, "fold": fold, "state": state, "window_bp": width,
            "foreground_promoters": len(foreground_promoters),
            "foreground_promoters_with_core": len({row["region_id"] for row in state_mappings}),
            "foreground_cores": len(state_mappings), "sequence_valid_cores": len(valid_fg),
            "controls_available": len(controls), "matched_pairs": len(pairs),
            "unmatched": len(unmatched), "matched_fraction": matched_fraction,
            "seed": seed, "motif_status": status,
        }
        sampling.append(sampling_row)
        if status == "planned":
            motif_jobs.append((contrast, directory))
    write_tsv(output / "03_motif_sampling.tsv", sampling)

    catalog = read_catalog(motif_db)
    ame_rows, homer_rows, motif_status = [], [], {}
    motif_path = dict(os.environ)
    motif_path["PATH"] = os.pathsep.join([
        str(analysis_env / "bin"), str(motif_env / "bin"), motif_path.get("PATH", "")
    ])
    for contrast, directory in motif_jobs:
        print(f"step=motifs contrast={contrast}", flush=True)
        jobs = []
        if cfg.get("run_ame"):
            ame_dir = directory / "AME"
            ame_dir.mkdir()
            ame_output = ame_dir / "ame.tsv"
            jobs.append(("AME", [
                executables["ame"], "--text", "--control", directory / "control.fa",
                "--method", "fisher", "--scoring", "max", "--noseq",
                "--evalue-report-threshold", "1000000000000", directory / "foreground.fa", motif_db,
            ], ame_output, ame_dir / "ame.stderr.log"))
        if cfg.get("run_homer"):
            jobs.append(("HOMER", [
                executables["homer"], directory / "foreground.fa", "fasta", directory / "HOMER",
                "-fastaBg", directory / "control.fa", "-mset", "vertebrates",
                "-p", str(cfg["homer_threads"]), "-nomotif",
            ], directory / "HOMER.stdout.log", directory / "HOMER.stderr.log"))
        try:
            with ThreadPoolExecutor(max_workers=max(1, len(jobs))) as pool:
                futures = [pool.submit(run_command, command, stdout, stderr, command_log, motif_path)
                           for _, command, stdout, stderr in jobs]
                for future in futures:
                    future.result()
            if cfg.get("run_ame"):
                ame_rows.extend(read_ame(directory / "AME/ame.tsv", contrast, catalog))
            if cfg.get("run_homer"):
                homer_rows.extend(read_homer(directory / "HOMER/knownResults.txt", contrast))
            motif_status[contrast] = "completed"
        except Exception as error:
            motif_status[contrast] = f"failed: {error}"
        write_tsv(output / "04_AME_all.tsv", ame_rows)
        write_tsv(output / "04_HOMER_all.tsv", homer_rows)
        (output / "motif_status.json").write_text(json.dumps(motif_status, indent=2), encoding="utf-8")
    if any(value.startswith("failed") for value in motif_status.values()):
        raise RuntimeError("At least one motif job failed; inspect motif_status.json and per-contrast logs")

    sampling_by_state = {row["state"]: row for row in sampling}
    if ame_rows:
        bh_adjust(ame_rows, "AME_threshold_adjusted_p", "AME_BH_q_all_contrasts")
        for contrast in sorted({row["contrast"] for row in ame_rows}):
            bh_adjust([row for row in ame_rows if row["contrast"] == contrast],
                      "AME_threshold_adjusted_p", "AME_BH_q_within_contrast")
        write_tsv(output / "04_AME_all.tsv", ame_rows)
        representatives = representative_ame(ame_rows)
        write_tsv(output / "04_AME_representative_families.tsv", representatives)
        significant = [row for row in representatives
                       if float(row["AME_BH_q_within_contrast"]) <= 0.05
                       and isinstance(row["enrichment_ratio"], (int, float))
                       and float(row["enrichment_ratio"]) > 1]
        write_tsv(output / "04_AME_significant_representative_families.tsv", significant)
        for rank_mode in ("q_rank", "enrichment_rank"):
            plot_combined(representatives, output, "AME", "AME_BH_q_within_contrast",
                          "motif_alt_id", int(cfg["top_motifs_per_panel"]), rank_mode,
                          sampling_by_state)
    if homer_rows:
        for contrast in sorted({row["contrast"] for row in homer_rows}):
            bh_adjust([row for row in homer_rows if row["contrast"] == contrast],
                      "p", "HOMER_BH_q_within_contrast")
        write_tsv(output / "04_HOMER_all.tsv", homer_rows)
        representatives = representative_homer(homer_rows)
        write_tsv(output / "04_HOMER_representative_families.tsv", representatives)
        significant = [row for row in representatives
                       if float(row["HOMER_BH_q_within_contrast"]) <= 0.05
                       and isinstance(row["enrichment_ratio"], (int, float))
                       and float(row["enrichment_ratio"]) > 1]
        write_tsv(output / "04_HOMER_significant_representative_families.tsv", significant)
        for rank_mode in ("q_rank", "enrichment_rank"):
            plot_combined(representatives, output, "HOMER", "HOMER_BH_q_within_contrast",
                          "motif", int(cfg["top_motifs_per_panel"]), rank_mode,
                          sampling_by_state)

    provenance = {
        "completed": datetime.now(timezone.utc).isoformat(),
        "definition": "P1_promoter_1p5_fold_500bp_ATAC_consensus_core",
        "analysis_unit": "same_final_enhancer_ATAC_consensus_core_representative_summit_inside_P1_promoter",
        "promoter_state_source": "archived_P1_full_H3_peak_support_consensus_state_matrix",
        "ATAC_read_filter": f"active-state raw tags >= {minimum}; explicitly validated",
        "ATAC_direction": "annotation_only_not_foreground_filter",
        "background": "same_P1_promoter_universe_active-state-eligible_nonforeground_ATAC_cores_GC_matched_without_replacement",
        "counts": {
            "P1_promoters": len(promoters), "source_ATAC_cores": len(cores),
            "mapped_ATAC_cores": len(mappings), "mapping_exclusions": dict(mapping_excluded),
            "promoter_states": dict(Counter(row["promoter_state"] for row in promoters)),
        },
        "input_sha256": {
            str(promoter_path): sha256(promoter_path), str(core_path): sha256(core_path),
            str(motif_db): sha256(motif_db),
        },
        "config": cfg,
    }
    (output / "run_metadata.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    (output / "completed.json").write_text(
        json.dumps({"completed": True, "results": str(output)}, indent=2), encoding="utf-8"
    )
    print(f"completed=1\nresults={output}", flush=True)


if __name__ == "__main__":
    main()
