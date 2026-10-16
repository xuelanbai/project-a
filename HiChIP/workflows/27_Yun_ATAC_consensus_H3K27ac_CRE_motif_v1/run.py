#!/usr/bin/env python3
"""Yun-inspired ATAC-consensus cores inside dynamic H3K27ac CREs.

Real BAM processing is intentionally restricted to a Linux Slurm compute job.
Biological-replicate P/FDR values are not inferred from pooled BAM files.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_schema(path: Path, label: str) -> tuple[list[tuple[str, int, int]], list[dict]]:
    rows, counts = [], Counter()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split()
            counts[len(cols)] += 1
            if len(cols) < 3:
                raise ValueError(f"{path}:{line_no}: fewer than three BED columns")
            chrom, start, end = canon(cols[0]), int(cols[1]), int(cols[2])
            if chrom in MAIN:
                if not 0 <= start < end:
                    raise ValueError(f"{path}:{line_no}: invalid BED interval")
                rows.append((chrom, start, end))
    if not rows:
        raise ValueError(f"No main-chromosome intervals in {path}")
    qc = [dict(label=label, path=str(path), columns=n, records=count,
               summit_column_10_available=int(n >= 10)) for n, count in sorted(counts.items())]
    return rows, qc


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


class RegionIndex:
    """Index non-overlapping consensus regions for summit-to-CRE assignment."""
    def __init__(self, rows):
        grouped = defaultdict(list)
        for row in rows:
            grouped[row["chrom"]].append(row)
        self.rows, self.starts = {}, {}
        for chrom, values in grouped.items():
            values.sort(key=lambda x: (x["start"], x["end"]))
            self.rows[chrom] = values
            self.starts[chrom] = [x["start"] for x in values]

    def contains(self, chrom: str, point: int):
        values = self.rows.get(canon(chrom), [])
        index = bisect.bisect_right(self.starts.get(canon(chrom), []), point) - 1
        if index >= 0 and point < values[index]["end"]:
            return values[index]
        return None


def read_promoters(gtf: Path, half_width: int, lengths: dict[str, int]):
    result = []
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) != 9 or cols[2] != "gene":
                continue
            chrom = canon(cols[0])
            if chrom not in lengths:
                continue
            tss = int(cols[3]) - 1 if cols[6] == "+" else int(cols[4]) - 1
            result.append((chrom, max(0, tss - half_width), min(lengths[chrom], tss + half_width + 1)))
    if not result:
        raise ValueError("No gene promoters parsed from GTF")
    return result


def consensus_h3(peaks: list[dict]):
    """Merge true overlaps only; book-ended intervals remain separate."""
    unique = {(r["chrom"], r["start"], r["end"], r["source_time"]): r for r in peaks}
    ordered = sorted(unique.values(), key=lambda r: (r["chrom"], r["start"], r["end"], r["source_time"]))
    regions, membership, current = [], [], None

    def finish(component):
        if component is None:
            return
        members = component.pop("_members")
        number = len(regions) + 1
        component["region_id"] = f"H3K27ac_CRE_{number:06d}"
        component["source_times"] = ";".join(sorted({x["source_time"] for x in members}))
        component["member_peak_count"] = len(members)
        for time in ("0h", "6h"):
            component[f"member_peak_count_{time}"] = sum(x["source_time"] == time for x in members)
            component[f"supported_{time}"] = int(component[f"member_peak_count_{time}"] > 0)
        component["length_bp"] = component["end"] - component["start"]
        regions.append(component)
        for member in members:
            membership.append(dict(region_id=component["region_id"], chrom=member["chrom"],
                                   source_time=member["source_time"], peak_start=member["start"],
                                   peak_end=member["end"]))

    for row in ordered:
        if current is None or row["chrom"] != current["chrom"] or row["start"] >= current["end"]:
            finish(current)
            current = dict(chrom=row["chrom"], start=row["start"], end=row["end"], _members=[row])
        else:
            current["end"] = max(current["end"], row["end"])
            current["_members"].append(row)
    finish(current)
    return regions, membership


def parse_macs3(path: Path, source_time: str) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 10:
                raise ValueError(f"{path}:{line_no}: MACS3 narrowPeak lacks column-10 summit")
            chrom, start, end = canon(cols[0]), int(cols[1]), int(cols[2])
            summit_offset = int(cols[9])
            summit = start + summit_offset
            if chrom not in MAIN:
                continue
            if not (0 <= start <= summit < end):
                raise ValueError(f"{path}:{line_no}: invalid summit")
            rows.append(dict(chrom=chrom, start=start, end=end, peak_name=cols[3],
                             score=float(cols[4]), signal_value=float(cols[6]),
                             minus_log10_p=float(cols[7]), minus_log10_q=float(cols[8]),
                             summit=summit, source_time=source_time, source_line=line_no))
    if not rows:
        raise ValueError(f"No main-chromosome MACS3 peaks in {path}")
    return rows


def consensus_atac(peaks: list[dict]):
    """Build cross-time ATAC consensus loci from overlapping full peak intervals.

    The consensus interval groups the original 0 h and 6 h MACS3 peaks and is
    used for same-coordinate ATAC quantification.  It does not invent a new
    summit: after quantification, the representative summit is taken from an
    original peak in the time point with the larger normalized ATAC signal.
    True overlaps are merged; book-ended intervals remain separate.
    """
    ordered = sorted(peaks, key=lambda r: (
        r["chrom"], r["start"], r["end"], r["source_time"], r["summit"], -r["signal_value"]))
    components, current = [], None
    for row in ordered:
        if current is None or row["chrom"] != current["chrom"] or row["start"] >= current["end"]:
            if current is not None:
                components.append(current)
            current = dict(chrom=row["chrom"], start=row["start"], end=row["end"], _members=[row])
        else:
            current["end"] = max(current["end"], row["end"])
            current["_members"].append(row)
    if current is not None:
        components.append(current)

    result, membership = [], []
    for number, component in enumerate(components, 1):
        members = component["_members"]
        core_id = f"ATAC_CORE_{number:06d}"
        locus = dict(
            core_id=core_id,
            chrom=component["chrom"],
            start=component["start"],
            end=component["end"],
            consensus_length_bp=component["end"] - component["start"],
            source_times=";".join(sorted({x["source_time"] for x in members})),
            member_peak_count=len(members),
            member_peak_count_0h=sum(x["source_time"] == "0h" for x in members),
            member_peak_count_6h=sum(x["source_time"] == "6h" for x in members),
            supported_0h=int(any(x["source_time"] == "0h" for x in members)),
            supported_6h=int(any(x["source_time"] == "6h" for x in members)),
            summit_span_bp=max(x["summit"] for x in members) - min(x["summit"] for x in members),
            _members=members,
        )
        result.append(locus)
        for member in members:
            membership.append(dict(
                core_id=core_id,
                consensus_chrom=locus["chrom"],
                consensus_start=locus["start"],
                consensus_end=locus["end"],
                **member,
            ))
    return result, membership


def write_bed(path: Path, rows: list[dict], chrom_map: dict[str, str],
              start_key="start", end_key="end", id_key="region_id"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in sorted(rows, key=lambda r: (r["chrom"], int(r[start_key]), int(r[end_key]), str(r[id_key]))):
            handle.write(f"{chrom_map[row['chrom']]}\t{int(row[start_key])}\t{int(row[end_key])}\t{row[id_key]}\t0\t.\n")


def run_command(command, stdout: Path, stderr: Path, command_log: Path, env=None):
    command = [str(x) for x in command]
    with command_log.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(command) + "\n")
    stdout.parent.mkdir(parents=True, exist_ok=True)
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        subprocess.run(command, stdout=out, stderr=err, check=True, env=env)


def multicov(bedtools: Path, bed: Path, bams: list[Path], labels: list[str],
             output: Path, command_log: Path, mapq: int) -> dict[str, dict[str, int]]:
    raw = output.with_suffix(".raw.tsv")
    command = [bedtools, "multicov", "-bed", bed, "-bams", *bams, "-q", str(mapq), "-p"]
    run_command(command, raw, output.with_suffix(".stderr.log"), command_log)
    return parse_multicov(raw, labels)


def parse_multicov(path: Path, labels: list[str]) -> dict[str, dict[str, int]]:
    """Parse BEDTools multicov while preserving any valid BED3+ columns."""
    counts = {}
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            cols = line.rstrip().split("\t")
            bed_columns = len(cols) - len(labels)
            if bed_columns < 4:
                raise ValueError(f"{path}:{line_no}: fewer than four input BED columns")
            identifier = cols[3]
            if identifier in counts:
                raise ValueError(f"Duplicate interval ID in multicov: {identifier}")
            counts[identifier] = {label: int(value) for label, value in zip(labels, cols[-len(labels):])}
    if not counts:
        raise ValueError(f"Empty multicov output: {path}")
    return counts


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
        rows.append(dict(contrast=contrast, motif_id=motif_id, motif_alt_id=alt_id,
                         AME_raw_optimal_p=number(get_field(raw, "p-value", "pvalue")),
                         AME_threshold_adjusted_p=number(get_field(raw, "adj_p-value", "q-value")),
                         AME_E_value=number(get_field(raw, "E-value")),
                         foreground_hit_percent=tp, control_hit_percent=fp,
                         enrichment_ratio=tp / fp if math.isfinite(tp) and math.isfinite(fp) and fp > 0 else "NA"))
    expected = {x[0] for x in catalog}
    observed = {x["motif_id"] for x in rows}
    if not rows or len(observed) != len(rows) or observed - expected:
        raise ValueError(f"Invalid, duplicate, or unexpected AME results: {path}")
    if len(observed) < 0.95 * len(expected):
        raise ValueError(f"AME reported only {len(observed)}/{len(expected)} motifs")
    for motif_id, alt_id in catalog:
        if motif_id not in observed:
            rows.append(dict(contrast=contrast, motif_id=motif_id, motif_alt_id=alt_id,
                             AME_raw_optimal_p=1.0, AME_threshold_adjusted_p=1.0,
                             AME_E_value="NA", foreground_hit_percent="NA",
                             control_hit_percent="NA", enrichment_ratio="NA"))
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
        rows.append(dict(contrast=contrast, motif=motif, p=number(row.get("P-value")), native_q=native_q,
                         foreground_percent=target, background_percent=background,
                         enrichment_ratio=target / background if background > 0 else "NA"))
    if not rows:
        raise ValueError(f"Empty HOMER results: {path}")
    return rows


def bh_adjust(rows: list[dict], p_field: str, q_field: str):
    indexed = sorted(((i, float(row[p_field])) for i, row in enumerate(rows)), key=lambda x: x[1])
    total, running = len(indexed), 1.0
    result = [1.0] * total
    for position in range(total - 1, -1, -1):
        running = min(running, indexed[position][1] * total / (position + 1))
        result[position] = min(1.0, running)
    for position, (index, _) in enumerate(indexed):
        rows[index][q_field] = result[position]


def motif_family(text: str) -> str:
    value = text.upper()
    rules = [
        (r"CTCF|BORIS", "CTCF"), (r"\b(JUN|FOS|BATF|JDP)", "AP-1"),
        (r"NFKB|NF.KB|\bREL[AB]?\b", "NF-kB"), (r"CEBP", "C/EBP"),
        (r"PPAR|RXR", "PPAR/RXR"), (r"SREBF", "SREBF"),
        (r"ATF4", "ISR/ATF4"), (r"\bKLF", "KLF"),
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
        rank = (float(candidate["AME_BH_q_within_contrast"]), -ratio,
                float(candidate["AME_threshold_adjusted_p"]), candidate["motif_id"])
        if key not in groups or rank < groups[key][0]:
            groups[key] = (rank, candidate)
    return [value[1] for key, value in sorted(groups.items())]


def overlap_window_index(rows, width):
    grouped = defaultdict(list)
    half = width // 2
    for row in rows:
        grouped[row["chrom"]].append((row["center"] - half, row["center"] + half))
    return {chrom: IntervalIndex([(chrom, start, end) for start, end in values]) for chrom, values in grouped.items()}


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


def selected(region, state: str, fold: float, minimum: int):
    time = "6h" if state == "strength" else "0h"
    if int(region[f"supported_{time}"]) != 1:
        return False
    if any(float(region[f"{assay}_{time}_raw_tags"]) < minimum for assay in ("H3K27ac", "ATAC")):
        return False
    cutoff = math.log2(fold)
    value = float(region["H3K27ac_log2FC"])
    return value >= cutoff if state == "strength" else value <= -cutoff


def active_eligible(region, state: str, minimum: int):
    time = "6h" if state == "strength" else "0h"
    return int(region[f"supported_{time}"]) == 1 and all(
        float(region[f"{assay}_{time}_raw_tags"]) >= minimum for assay in ("H3K27ac", "ATAC"))


def make_plot(rows: list[dict], output: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for contrast in sorted({row["contrast"] for row in rows}):
        current = [row for row in rows if row["contrast"] == contrast and
                   isinstance(row["enrichment_ratio"], (int, float)) and row["enrichment_ratio"] > 1]
        current.sort(key=lambda r: (float(r["AME_BH_q_within_contrast"]), -float(r["enrichment_ratio"])))
        current = current[:15]
        fig, ax = plt.subplots(figsize=(8, max(3, 0.42 * len(current) + 1.5)))
        if current:
            values = [math.log2(float(row["enrichment_ratio"])) for row in current]
            colors = ["#C43C39" if float(row["AME_BH_q_within_contrast"]) <= 0.05 else "#AAAAAA" for row in current]
            ax.barh(range(len(current)), values, color=colors)
            ax.set_yticks(range(len(current)), [f"{row['motif_family']} | {row['motif_alt_id']}"[:65] for row in current])
            ax.invert_yaxis()
            ax.set_xlabel("log2(foreground / GC-matched control)")
        else:
            ax.text(0.5, 0.5, "No motif enrichment ratio > 1", ha="center", transform=ax.transAxes)
        ax.set_title(contrast)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(output / f"Figure_AME_representative_{contrast}.png", dpi=180)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--resume", type=Path,
                        help="Reuse MACS3 peaks in an existing incomplete run directory")
    parser.add_argument("--reuse-peaks-from", type=Path,
                        help="Copy ATAC MACS3 peak files from an earlier run into a new run")
    args = parser.parse_args()
    if args.resume and args.reuse_peaks_from:
        raise ValueError("--resume and --reuse-peaks-from are mutually exclusive")
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Submit this real-data analysis as a Linux Slurm job; do not run it on the login node.")
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    if cfg["fold_thresholds"] != [1.2, 1.5] or cfg["motif_window_sizes"] != [200, 500]:
        raise ValueError("This version fixes fold thresholds [1.2, 1.5] and motif windows [200, 500]")
    if cfg["run_de_novo"]:
        raise ValueError("This exploratory version intentionally disables HOMER de novo motifs")
    if cfg["homer_threads"] + 1 > int(os.environ.get("SLURM_CPUS_PER_TASK", "4")):
        raise ValueError("Allocated CPUs are insufficient for concurrent HOMER + AME")

    project = Path(cfg["project"])
    fasta, gtf = Path(cfg["reference_fasta"]), Path(cfg["gtf"])
    bam_paths = {key: project / "02.bam" / value for key, value in cfg["bam_names"].items()}
    peak_paths = {key: project / "06.peak" / value for key, value in cfg["peak_names"].items()}
    legacy_atac_paths = {key: project / "06.peak" / value
                         for key, value in cfg.get("legacy_ATAC_peak_names_for_QC_only", {}).items()}
    analysis_env, macs3_env, bedtools_env, motif_env = map(
        Path, (cfg["analysis_env"], cfg["macs3_env"], cfg["bedtools_env"], cfg["motif_env"]))
    executables = {
        "macs3": macs3_env / "bin/macs3",
        "bedtools": bedtools_env / "bin/bedtools",
        "homer": analysis_env / "bin/findMotifs.pl",
        "ame": motif_env / "bin/ame",
    }
    required = [fasta, Path(str(fasta) + ".fai"), gtf, Path(cfg["motif_database"]),
                *bam_paths.values(), *peak_paths.values(), *legacy_atac_paths.values(), *executables.values()]
    for path in required:
        if not path.is_absolute() or not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    for path in bam_paths.values():
        index = Path(str(path) + ".bai")
        if not index.is_file():
            raise FileNotFoundError(index)

    import pysam
    genome = pysam.FastaFile(str(fasta))
    fasta_map = {canon(name): name for name in genome.references if canon(name) in MAIN}
    lengths = {chrom: genome.get_reference_length(name) for chrom, name in fasta_map.items()}
    if set(lengths) != MAIN:
        raise ValueError(f"Reference main chromosomes differ: {sorted(set(lengths) ^ MAIN)}")
    bam_maps = {}
    for label, path in bam_paths.items():
        with pysam.AlignmentFile(str(path), "rb") as bam:
            if not bam.has_index():
                raise ValueError(f"BAM is not indexed: {path}")
            bam_maps[label] = {canon(name): name for name in bam.references if canon(name) in MAIN}
            if set(bam_maps[label]) != MAIN:
                raise ValueError(f"BAM main chromosome mismatch: {label}")
            for chrom in MAIN:
                if bam.get_reference_length(bam_maps[label][chrom]) != lengths[chrom]:
                    raise ValueError(f"Assembly mismatch: {label} {chrom}")
    first_map = bam_maps["ATAC_0h"]
    if any(mapping != first_map for mapping in bam_maps.values()):
        raise ValueError("BAM chromosome naming differs; bedtools multicov requires a shared naming scheme")

    if args.resume:
        output = args.resume.resolve()
        if not output.is_absolute() or not output.is_dir():
            raise ValueError("--resume must name an existing absolute run directory")
        expected_root = Path(cfg["results"]).resolve()
        if output.parent != expected_root:
            raise ValueError(f"Resume directory must be directly under {expected_root}")
        if (output / "completed.json").exists():
            raise ValueError("Refusing to resume a completed run")
        print(f"resume={output}", flush=True)
    else:
        output = Path(cfg["results"]) / f"run_{os.environ['SLURM_JOB_ID']}"
        output.mkdir(parents=True, exist_ok=False)
    for name in ("logs", "peaks", "counts", "motifs", "code"):
        (output / name).mkdir(exist_ok=True)
    command_log = output / "logs/commands.jsonl"
    for name in ("run.py", "config.json", "README.md", "run.slurm", "submit.sh", "check_env.sh", "test_small.py"):
        shutil.copy2(Path(__file__).with_name(name), output / "code" / name)
    print(f"results={output}", flush=True)

    schema_qc, beds = [], {}
    for label, path in peak_paths.items():
        beds[label], qc = file_schema(path, label)
        schema_qc.extend(qc)
    for label, path in legacy_atac_paths.items():
        _, qc = file_schema(path, label)
        schema_qc.extend(qc)
    write_tsv(output / "01_input_peak_schema.tsv", schema_qc)

    reuse_peak_dir = None
    if args.reuse_peaks_from:
        reuse_peak_dir = args.reuse_peaks_from.resolve()
        if not reuse_peak_dir.is_absolute() or not reuse_peak_dir.is_dir():
            raise ValueError("--reuse-peaks-from must name an existing absolute run directory")

    # Re-call ATAC peaks/summits from merged BAMs. BAMPE is explicit because
    # MACS3 cannot auto-detect paired-end mode.
    genome_size = sum(lengths.values())
    atac_peaks = []
    for time in ("0h", "6h"):
        label = f"ATAC_{time}"
        narrow_peak = output / f"peaks/{label}_peaks.narrowPeak"
        if reuse_peak_dir is not None:
            source_peak = reuse_peak_dir / f"peaks/{label}_peaks.narrowPeak"
            if not source_peak.is_file() or source_peak.stat().st_size == 0:
                raise FileNotFoundError(source_peak)
            shutil.copy2(source_peak, narrow_peak)
            print(f"step=MACS3_ATAC_{time} action=reuse_external source={source_peak}", flush=True)
        elif args.resume and narrow_peak.is_file() and narrow_peak.stat().st_size > 0:
            print(f"step=MACS3_ATAC_{time} action=reuse_existing", flush=True)
        else:
            command = [executables["macs3"], "callpeak", "-t", bam_paths[label], "-f", "BAMPE",
                       "-g", str(genome_size), "-n", label, "--outdir", output / "peaks",
                       "-q", str(cfg["macs3_q"]), "--keep-dup", "all", "--call-summits"]
            print(f"step=MACS3_ATAC_{time} action=run", flush=True)
            run_command(command, output / f"logs/MACS3_{time}.stdout.log",
                        output / f"logs/MACS3_{time}.stderr.log", command_log)
        atac_peaks.extend(parse_macs3(narrow_peak, time))
    with gzip.open(output / "02_ATAC_called_peaks.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        fields = ["chrom", "start", "end", "peak_name", "score", "signal_value", "minus_log10_p",
                  "minus_log10_q", "summit", "source_time", "source_line"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(atac_peaks)

    atac_indexes = {time: IntervalIndex((r["chrom"], r["start"], r["end"])
                                         for r in atac_peaks if r["source_time"] == time)
                    for time in ("0h", "6h")}
    promoters = IntervalIndex(read_promoters(gtf, int(cfg["promoter_half_width"]), lengths))
    h3k4 = IntervalIndex([*beds["H3K4me3_0h"], *beds["H3K4me3_6h"]])
    candidates, excluded = [], Counter()
    for time in ("0h", "6h"):
        for chrom, start, end in beds[f"H3K27ac_{time}"]:
            if end > lengths.get(chrom, -1):
                raise ValueError(f"H3K27ac interval outside reference: {chrom}:{start}-{end}")
            if not atac_indexes[time].overlaps(chrom, start, end):
                excluded["no_same_time_recalled_ATAC_peak"] += 1; continue
            if promoters.overlaps(chrom, start, end):
                excluded["promoter_overlap"] += 1; continue
            if h3k4.overlaps(chrom, start, end):
                excluded["H3K4me3_overlap_either_time"] += 1; continue
            candidates.append(dict(chrom=chrom, start=start, end=end, source_time=time))
    regions, h3_membership = consensus_h3(candidates)
    if not regions:
        raise ValueError("No H3K27ac consensus CREs after filtering")
    write_tsv(output / "03_H3K27ac_consensus_membership.tsv", h3_membership)
    write_bed(output / "03_H3K27ac_consensus_CRE.bed", regions, first_map)

    # Yun-inspired step: full 0 h and 6 h ATAC peak intervals define a shared
    # consensus locus.  Quantification is performed over that same locus in
    # both ATAC BAMs before selecting an original peak summit.
    cores, core_membership = consensus_atac(atac_peaks)
    write_tsv(output / "02_ATAC_consensus_peak_membership.tsv", core_membership)
    write_bed(output / "02_ATAC_consensus_peaks.bed", cores, first_map,
              "start", "end", "core_id")

    # BEDs for BAM quantification all share exact coordinates across samples.
    write_bed(output / "counts/CREs.bed", regions, first_map)
    write_bed(output / "counts/ATAC_consensus_peaks.bed", cores, first_map,
              "start", "end", "core_id")
    main_rows = [dict(chrom=chrom, start=0, end=lengths[chrom], interval_id=f"CHR_{chrom}") for chrom in sorted(MAIN)]
    write_bed(output / "counts/main_chromosomes.bed", main_rows, first_map, id_key="interval_id")
    labels = list(bam_paths)
    libraries = multicov(executables["bedtools"], output / "counts/main_chromosomes.bed",
                         [bam_paths[x] for x in labels], labels, output / "counts/library_sizes.tsv",
                         command_log, int(cfg["min_mapq"]))
    library_totals = {label: sum(row[label] for row in libraries.values()) for label in labels}
    if any(value <= 0 for value in library_totals.values()):
        raise ValueError(f"Invalid library totals: {library_totals}")
    cre_counts = multicov(executables["bedtools"], output / "counts/CREs.bed",
                          [bam_paths[x] for x in labels], labels, output / "counts/CRE_counts.tsv",
                          command_log, int(cfg["min_mapq"]))
    atac_labels = ["ATAC_0h", "ATAC_6h"]
    core_counts = multicov(executables["bedtools"], output / "counts/ATAC_consensus_peaks.bed",
                           [bam_paths[x] for x in atac_labels], atac_labels,
                           output / "counts/ATAC_core_counts.tsv", command_log, int(cfg["min_mapq"]))
    scale = float(cfg["normalization_total_tags"])
    pseudo = float(cfg["normalized_pseudocount"])
    for region in regions:
        counts = cre_counts[region["region_id"]]
        for label in labels:
            region[f"{label}_raw_tags"] = counts[label]
            region[f"{label}_tags_per_10M"] = counts[label] / library_totals[label] * scale
        for assay in ("H3K27ac", "ATAC"):
            region[f"{assay}_log2FC"] = math.log2(
                (region[f"{assay}_6h_tags_per_10M"] + pseudo) /
                (region[f"{assay}_0h_tags_per_10M"] + pseudo))
        atac = region["ATAC_log2FC"]
        region["ATAC_stratum"] = "ATAC_up" if atac > 0 else "ATAC_down" if atac < 0 else "ATAC_unchanged"
    write_tsv(output / "03_CRE_signal_matrix.tsv", regions)

    for core in cores:
        counts = core_counts[core["core_id"]]
        for label in atac_labels:
            core[f"{label}_raw_tags"] = counts[label]
            core[f"{label}_tags_per_10M"] = counts[label] / library_totals[label] * scale
        available = {x["source_time"] for x in core["_members"]}
        if len(available) == 1:
            selected_time = next(iter(available))
            selection_reason = "only_time_with_called_summit"
        else:
            zero = core["ATAC_0h_tags_per_10M"]
            six = core["ATAC_6h_tags_per_10M"]
            if zero != six:
                selected_time = "0h" if zero > six else "6h"
                selection_reason = "larger_normalized_ATAC_count_in_consensus_peak"
            else:
                max_signal = {time: max(x["signal_value"] for x in core["_members"] if x["source_time"] == time)
                              for time in available}
                selected_time = min(available, key=lambda time: (-max_signal[time], time))
                selection_reason = "ATAC_count_tie_then_larger_MACS3_signal"
        options = [x for x in core["_members"] if x["source_time"] == selected_time]
        chosen = min(options, key=lambda x: (-x["signal_value"], -x["minus_log10_q"], x["summit"]))
        core["representative_summit"] = chosen["summit"]
        core["representative_source_time"] = selected_time
        core["representative_MACS3_signal"] = chosen["signal_value"]
        core["representative_peak_name"] = chosen["peak_name"]
        core["summit_selection_reason"] = selection_reason
        core["summit_start"] = chosen["summit"]
        core["summit_end"] = chosen["summit"] + 1
    write_tsv(output / "02_ATAC_consensus_cores.tsv", cores)
    write_bed(output / "02_ATAC_representative_summits.bed", cores, first_map,
              "summit_start", "summit_end", "core_id")

    region_index = RegionIndex(regions)
    mappings, mapping_excluded = [], Counter()
    for core in cores:
        region = region_index.contains(core["chrom"], core["representative_summit"])
        if region is None:
            mapping_excluded["representative_summit_outside_H3K27ac_CRE"] += 1
            continue
        row = {key: value for key, value in core.items() if not key.startswith("_")}
        row.update(region_id=region["region_id"], CRE_start=region["start"], CRE_end=region["end"],
                   center=core["representative_summit"], cores_per_CRE=0, _region=region)
        mappings.append(row)
    per_cre = Counter(row["region_id"] for row in mappings)
    for row in mappings:
        row["cores_per_CRE"] = per_cre[row["region_id"]]
    write_tsv(output / "04_ATAC_core_to_CRE.tsv", mappings)

    # Attach fixed-length sequences; an ATAC core is the motif analysis unit.
    sequence_exclusions = []
    for row in mappings:
        row["_seqs"], row["_gcs"] = {}, {}
        for width in cfg["motif_window_sizes"]:
            start, end = row["center"] - width // 2, row["center"] + width // 2
            reason = None
            if start < 0 or end > lengths[row["chrom"]]:
                reason = "chromosome_edge"
            else:
                sequence = genome.fetch(fasta_map[row["chrom"]], start, end).upper()
                acgt = sum(sequence.count(base) for base in "ACGT")
                if len(sequence) != width or not acgt or (len(sequence) - acgt) / len(sequence) > cfg["max_N_fraction"]:
                    reason = "invalid_or_N_rich_sequence"
                else:
                    row["_seqs"][width] = sequence
                    row["_gcs"][width] = (sequence.count("G") + sequence.count("C")) / acgt
            if reason:
                sequence_exclusions.append(dict(core_id=row["core_id"], width=width, reason=reason))
    genome.close()
    write_tsv(output / "04_motif_sequence_exclusions.tsv", sequence_exclusions,
              ["core_id", "width", "reason"])

    motif_jobs, sampling = [], []
    seed_counter = 0
    for fold in cfg["fold_thresholds"]:
        fold_token = str(fold).replace(".", "p")
        fold_dir = output / f"fold_{fold_token}"
        (fold_dir / "regions").mkdir(parents=True)
        for state in ("strength", "weaken"):
            foreground_cre = [region for region in regions if selected(region, state, float(fold), int(cfg["min_raw_tags_active_state"]))]
            write_bed(fold_dir / "regions" / f"{state}_CRE.bed", foreground_cre, first_map)
            foreground_ids = {region["region_id"] for region in foreground_cre}
            state_mappings = [row for row in mappings if row["region_id"] in foreground_ids]
            write_tsv(fold_dir / "regions" / f"{state}_ATAC_core_membership.tsv", [
                dict({k: v for k, v in row.items() if not k.startswith("_")},
                     H3K27ac_log2FC=row["_region"]["H3K27ac_log2FC"],
                     ATAC_log2FC=row["_region"]["ATAC_log2FC"],
                     ATAC_stratum=row["_region"]["ATAC_stratum"]) for row in state_mappings])
            for width in cfg["motif_window_sizes"]:
                valid_fg = [row for row in state_mappings if width in row["_seqs"]]
                for row in valid_fg:
                    row["window_start"], row["window_end"] = row["center"] - width // 2, row["center"] + width // 2
                write_bed(fold_dir / "regions" / f"{state}_ATAC_core_{width}bp.bed", valid_fg,
                          first_map, "window_start", "window_end", "core_id")
                foreground_core_ids = {row["core_id"] for row in valid_fg}
                window_index = overlap_window_index(valid_fg, width)
                controls = [row for row in mappings if width in row["_seqs"] and
                            row["core_id"] not in foreground_core_ids and
                            active_eligible(row["_region"], state, int(cfg["min_raw_tags_active_state"])) and
                            not window_overlaps(window_index, row, width)]
                for row in valid_fg + controls:
                    row["_gc"] = row["_gcs"][width]
                seed = int(cfg["seed"]) + seed_counter
                seed_counter += 1
                pairs, unmatched = match_gc(valid_fg, controls, float(cfg["max_GC_difference"]), seed)
                matched_fraction = len(pairs) / len(valid_fg) if valid_fg else 0.0
                contrast = f"{state}_{fold_token}_{width}bp_ATAC_core"
                directory = fold_dir / "motifs" / f"{state}_{width}bp"
                directory.mkdir(parents=True)
                write_tsv(directory / "background_pairs.tsv", [dict(
                    foreground_id=fg["core_id"], control_id=bg["core_id"],
                    foreground_GC=fg["_gc"], control_GC=bg["_gc"], delta_GC=abs(fg["_gc"] - bg["_gc"]))
                    for fg, bg in pairs])
                write_tsv(directory / "unmatched_foreground.tsv", [dict(core_id=x) for x in unmatched], ["core_id"])
                for label, position in (("foreground", 0), ("control", 1)):
                    with (directory / f"{label}.fa").open("w", encoding="utf-8") as fa:
                        for pair in pairs:
                            row = pair[position]
                            fa.write(f">{row['core_id']}\n{row['_seqs'][width]}\n")
                status = "planned"
                if len(pairs) < int(cfg["minimum_motif_pairs"]):
                    status = "not_tested_fewer_than_minimum_pairs"
                elif matched_fraction < float(cfg["minimum_matched_fraction"]):
                    status = "not_tested_low_matched_fraction"
                sampling.append(dict(contrast=contrast, fold=fold, state=state, window_bp=width,
                                     foreground_CREs=len(foreground_cre), foreground_cores=len(state_mappings),
                                     sequence_valid_cores=len(valid_fg), controls_available=len(controls),
                                     matched_pairs=len(pairs), unmatched=len(unmatched),
                                     matched_fraction=matched_fraction, seed=seed, motif_status=status))
                if status == "planned":
                    motif_jobs.append((contrast, directory))
    write_tsv(output / "05_motif_sampling.tsv", sampling)

    catalog = read_catalog(Path(cfg["motif_database"]))
    ame_rows, homer_rows, motif_status = [], [], {}
    motif_env_path = dict(os.environ)
    motif_env_path["PATH"] = os.pathsep.join([str(analysis_env / "bin"), str(motif_env / "bin"),
                                               motif_env_path.get("PATH", "")])
    for contrast, directory in motif_jobs:
        print(f"step=motifs contrast={contrast}", flush=True)
        ame_dir = directory / "AME"
        ame_dir.mkdir()
        ame_output = ame_dir / "ame.tsv"
        ame_command = [executables["ame"], "--text", "--control", directory / "control.fa",
                       "--method", "fisher", "--scoring", "max", "--noseq",
                       "--evalue-report-threshold", "1000000000000", directory / "foreground.fa",
                       cfg["motif_database"]]
        homer_command = [executables["homer"], directory / "foreground.fa", "fasta", directory / "HOMER",
                         "-fastaBg", directory / "control.fa", "-mset", "vertebrates", "-p",
                         str(cfg["homer_threads"]), "-nomotif"]
        try:
            with ThreadPoolExecutor(max_workers=2) as pool:
                futures = [
                    pool.submit(run_command, ame_command, ame_output, ame_dir / "ame.stderr.log", command_log, motif_env_path),
                    pool.submit(run_command, homer_command, directory / "HOMER.stdout.log",
                                directory / "HOMER.stderr.log", command_log, motif_env_path),
                ]
                for future in futures:
                    future.result()
            ame_rows.extend(read_ame(ame_output, contrast, catalog))
            homer_rows.extend(read_homer(directory / "HOMER/knownResults.txt", contrast))
            motif_status[contrast] = "completed"
        except Exception as error:
            motif_status[contrast] = f"failed: {error}"
        write_tsv(output / "05_AME_all.tsv", ame_rows)
        write_tsv(output / "05_HOMER_all.tsv", homer_rows)
        (output / "motif_status.json").write_text(json.dumps(motif_status, indent=2), encoding="utf-8")
    if any(value.startswith("failed") for value in motif_status.values()):
        raise RuntimeError("At least one motif job failed; inspect motif_status.json and per-contrast logs")

    if ame_rows:
        bh_adjust(ame_rows, "AME_threshold_adjusted_p", "AME_BH_q_all_contrasts")
        for contrast in sorted({row["contrast"] for row in ame_rows}):
            subset = [row for row in ame_rows if row["contrast"] == contrast]
            bh_adjust(subset, "AME_threshold_adjusted_p", "AME_BH_q_within_contrast")
        write_tsv(output / "05_AME_all.tsv", ame_rows)
        representatives = representative_ame(ame_rows)
        write_tsv(output / "05_AME_representative_families.tsv", representatives)
        make_plot(representatives, output)

    provenance = {
        "completed": datetime.now(timezone.utc).isoformat(),
        "mode": "pooled_descriptive_no_biological_replicate_P_or_FDR",
        "analysis_unit": "ATAC_summit_core_inside_H3K27ac_consensus_CRE",
        "enhancer_coordinates": "full_H3K27ac_BED_intervals_supported_by_same_time_recalled_ATAC_peak",
        "ATAC_consensus": "true-overlap_union_of_full_0h_6h_MACS3_BAMPE_peak_intervals",
        "representative_summit": "original_MACS3_summit_from_higher_normalized_ATAC_signal_timepoint_within_each_consensus_peak",
        "ATAC_direction": "annotation_only_not_foreground_filter",
        "signal_normalization": "same_coordinate_alignment_counts_per_10M_main_chromosome_alignments_no_length_division_no_median_centering",
        "count_unit": "bedtools_multicov_alignments_with_MAPQ_and_proper_pair_filter_not_fragment_counts",
        "library_totals": library_totals,
        "input_BAM_stats": {label: {"path": str(path), "bytes": path.stat().st_size,
                                      "mtime_ns": path.stat().st_mtime_ns}
                            for label, path in bam_paths.items()},
        "reused_ATAC_peak_run": str(reuse_peak_dir) if reuse_peak_dir else None,
        "H3K27ac_candidate_exclusions": dict(excluded),
        "ATAC_core_mapping_exclusions": dict(mapping_excluded),
        "counts": {"called_ATAC_peak_summits": len(atac_peaks), "ATAC_cores": len(cores),
                   "supported_H3K27ac_peaks": len(candidates), "H3K27ac_CREs": len(regions),
                   "mapped_ATAC_cores": len(mappings)},
        "input_sha256": {str(path): sha256(path) for path in
                          [*peak_paths.values(), *legacy_atac_paths.values(), Path(cfg["motif_database"])]},
        "config": cfg,
    }
    (output / "run_metadata.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    (output / "completed.json").write_text(json.dumps({"completed": True, "results": str(output)}, indent=2), encoding="utf-8")
    print(f"completed=1\nresults={output}", flush=True)


if __name__ == "__main__":
    main()
