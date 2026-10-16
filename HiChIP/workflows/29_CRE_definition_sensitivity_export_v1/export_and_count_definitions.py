#!/usr/bin/env python3
"""Export and count promoter/enhancer candidates under explicit definitions.

No BAM is read and no single biological set is selected. This is a lightweight
definition-sensitivity export using the v27 recalled ATAC peaks, source histone
peak BEDs and Ensembl gene TSS annotations.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import os
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


MAIN = {str(value) for value in range(1, 19)} | {"X", "Y"}


def canon(chrom: str) -> str:
    value = str(chrom).strip()
    return value[3:] if value.startswith("chr") else value


def chrom_order(chrom: str) -> tuple[int, str]:
    value = canon(chrom)
    if value.isdigit():
        return int(value), value
    return {"X": 19, "Y": 20}.get(value, 99), value


def ensure_file(path: Path) -> Path:
    path = path.resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {path}")
    return path


def read_fai(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    names, lengths = {}, {}
    with ensure_file(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            cols = line.rstrip().split("\t")
            if len(cols) < 2:
                raise ValueError(f"{path}:{line_no}: invalid FAI row")
            chrom = canon(cols[0])
            if chrom in MAIN:
                names[chrom] = cols[0]
                lengths[chrom] = int(cols[1])
    if set(lengths) != MAIN:
        raise ValueError(f"Reference main-chromosome mismatch: {sorted(set(lengths) ^ MAIN)}")
    return names, lengths


def read_bed(path: Path, lengths: dict[str, int]) -> list[tuple[str, int, int]]:
    rows = set()
    with ensure_file(path).open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split()
            if len(cols) < 3:
                raise ValueError(f"{path}:{line_no}: fewer than three BED columns")
            chrom, start, end = canon(cols[0]), int(cols[1]), int(cols[2])
            if chrom not in MAIN:
                continue
            if not 0 <= start < end <= lengths[chrom]:
                raise ValueError(f"{path}:{line_no}: invalid interval {chrom}:{start}-{end}")
            rows.add((chrom, start, end))
    if not rows:
        raise ValueError(f"No main-chromosome intervals in {path}")
    return sorted(rows, key=lambda row: (chrom_order(row[0]), row[1], row[2]))


def merge_true_overlaps(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Merge true overlaps only; book-ended half-open intervals remain separate."""
    result: list[list] = []
    for chrom, start, end in sorted(set(intervals), key=lambda row: (chrom_order(row[0]), row[1], row[2])):
        if result and chrom == result[-1][0] and start < result[-1][2]:
            result[-1][2] = max(result[-1][2], end)
        else:
            result.append([chrom, start, end])
    return [(row[0], int(row[1]), int(row[2])) for row in result]


class MergedIntervalIndex:
    def __init__(self, intervals: list[tuple[str, int, int]]):
        grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in intervals:
            grouped[canon(chrom)].append((int(start), int(end)))
        self.data: dict[str, tuple[list[int], list[int]]] = {}
        for chrom, values in grouped.items():
            merged = []
            for start, end in sorted(set(values)):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self.data[chrom] = ([row[0] for row in merged], [row[1] for row in merged])

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        starts, ends = self.data.get(canon(chrom), ([], []))
        index = bisect.bisect_left(starts, end) - 1
        return index >= 0 and ends[index] > start

    def overlap_segments(self, chrom: str, start: int, end: int) -> list[tuple[str, int, int]]:
        chrom = canon(chrom)
        starts, ends = self.data.get(chrom, ([], []))
        index = bisect.bisect_right(ends, start)
        result = []
        while index < len(starts) and starts[index] < end:
            left, right = max(start, starts[index]), min(end, ends[index])
            if left < right:
                result.append((chrom, left, right))
            index += 1
        return result

    def intervals(self) -> list[tuple[str, int, int]]:
        result = []
        for chrom, (starts, ends) in self.data.items():
            result.extend((chrom, start, end) for start, end in zip(starts, ends))
        return sorted(result, key=lambda row: (chrom_order(row[0]), row[1], row[2]))


def intersect_interval_sets(
    left: MergedIntervalIndex,
    right: MergedIntervalIndex,
) -> list[tuple[str, int, int]]:
    result = []
    left_by_chrom = defaultdict(list)
    right_by_chrom = defaultdict(list)
    for chrom, start, end in left.intervals():
        left_by_chrom[chrom].append((start, end))
    for chrom, start, end in right.intervals():
        right_by_chrom[chrom].append((start, end))
    for chrom in MAIN:
        first, second = left_by_chrom[chrom], right_by_chrom[chrom]
        i = j = 0
        while i < len(first) and j < len(second):
            start, end = max(first[i][0], second[j][0]), min(first[i][1], second[j][1])
            if start < end:
                result.append((chrom, start, end))
            if first[i][1] <= second[j][1]:
                i += 1
            else:
                j += 1
    return result


ATTRIBUTE = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTRIBUTE.findall(text)}


def read_genes(gtf: Path, lengths: dict[str, int]) -> list[dict]:
    opener = gzip.open if gtf.suffix == ".gz" else open
    genes = {}
    with opener(ensure_file(gtf), "rt", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) != 9 or cols[2] != "gene":
                continue
            chrom = canon(cols[0])
            if chrom not in MAIN:
                continue
            attrs = parse_attributes(cols[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                raise ValueError(f"{gtf}:{line_no}: gene lacks gene_id")
            strand = cols[6]
            tss = int(cols[3]) - 1 if strand == "+" else int(cols[4]) - 1
            row = {
                "gene_id": gene_id,
                "gene_name": attrs.get("gene_name") or gene_id,
                "chrom": chrom,
                "tss": tss,
                "strand": strand,
            }
            previous = genes.get(gene_id)
            if previous is not None and previous != row:
                raise ValueError(f"Inconsistent duplicate gene: {gene_id}")
            genes[gene_id] = row
    return sorted(genes.values(), key=lambda row: (chrom_order(row["chrom"]), row["tss"], row["gene_id"]))


def promoter_windows(genes: list[dict], half_width: int, lengths: dict[str, int]) -> list[dict]:
    result = []
    for gene in genes:
        row = dict(gene)
        row["start"] = max(0, gene["tss"] - half_width)
        row["end"] = min(lengths[gene["chrom"]], gene["tss"] + half_width + 1)
        result.append(row)
    return result


class PromoterIndex:
    def __init__(self, rows: list[dict]):
        grouped = defaultdict(list)
        for row in rows:
            grouped[row["chrom"]].append(row)
        self.data = {}
        for chrom, values in grouped.items():
            values.sort(key=lambda row: (row["start"], row["end"], row["gene_id"]))
            starts = [row["start"] for row in values]
            ends = [row["end"] for row in values]
            if ends != sorted(ends):
                raise ValueError("Promoter ends are not monotonic")
            self.data[chrom] = (values, starts, ends)

    def query(self, chrom: str, start: int, end: int) -> list[dict]:
        values, starts, ends = self.data.get(canon(chrom), ([], [], []))
        first = bisect.bisect_right(ends, start)
        stop = bisect.bisect_left(starts, end)
        return [row for row in values[first:stop] if row["start"] < end and row["end"] > start]


def percentile(values: list[int], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def interval_summary(intervals: list[tuple[str, int, int]]) -> dict:
    unique = sorted(set(intervals), key=lambda row: (chrom_order(row[0]), row[1], row[2]))
    union = merge_true_overlaps(unique)
    lengths = [end - start for _, start, end in union]
    return {
        "coordinate_intervals_unique": len(unique),
        "regions_after_true_overlap_merge": len(union),
        "union_total_bp": sum(lengths),
        "union_median_bp": percentile(lengths, 0.50),
        "union_p95_bp": percentile(lengths, 0.95),
        "union_max_bp": max(lengths) if lengths else 0,
    }


def consensus_components(by_time: dict[str, list[tuple[str, int, int]]]) -> list[dict]:
    rows = []
    for time, intervals in by_time.items():
        rows.extend((chrom, start, end, time) for chrom, start, end in set(intervals))
    rows.sort(key=lambda row: (chrom_order(row[0]), row[1], row[2], row[3]))
    components, current = [], None

    def finish(value):
        if value is None:
            return
        times = value.pop("times")
        value["supported_0h"] = int("0h" in times)
        value["supported_6h"] = int("6h" in times)
        if value["supported_0h"] and value["supported_6h"]:
            value["state"] = "both"
        elif value["supported_0h"]:
            value["state"] = "0h_only"
        else:
            value["state"] = "6h_only"
        components.append(value)

    for chrom, start, end, time in rows:
        if current is None or chrom != current["chrom"] or start >= current["end"]:
            finish(current)
            current = {"chrom": chrom, "start": start, "end": end, "times": {time}}
        else:
            current["end"] = max(current["end"], end)
            current["times"].add(time)
    finish(current)
    return components


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_bed6(
    path: Path,
    intervals: list[tuple[str, int, int]],
    chrom_names: dict[str, str],
    prefix: str,
) -> int:
    """Write unique source intervals without merging them."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(set(intervals), key=lambda row: (chrom_order(row[0]), row[1], row[2]))
    with path.open("w", encoding="utf-8", newline="") as handle:
        for number, (chrom, start, end) in enumerate(rows, 1):
            handle.write(f"{chrom_names[canon(chrom)]}\t{start}\t{end}\t{prefix}_{number:06d}\t0\t.\n")
    return len(rows)


def write_consensus_bed6(
    path: Path,
    components: list[dict],
    chrom_names: dict[str, str],
    prefix: str,
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for number, row in enumerate(components, 1):
            handle.write(
                f'{chrom_names[canon(row["chrom"])]}\t{row["start"]}\t{row["end"]}\t'
                f'{prefix}_{number:06d}\t0\t.\n'
            )
    return len(components)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with ensure_file(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def promoter_definitions(
    time: str,
    h3_peaks: list[tuple[str, int, int]],
    promoters: PromoterIndex,
    atac: MergedIntervalIndex,
    h3k4: MergedIntervalIndex,
) -> tuple[list[dict], dict[str, dict]]:
    common_index = MergedIntervalIndex(intersect_interval_sets(atac, h3k4))
    tss_h3_count = atac_full_count = 0
    full_peaks, full_genes, full_assoc = set(), set(), set()
    local_peaks, local_genes, local_assoc, clipped = set(), set(), set(), []
    common_peaks, common_genes, common_assoc, common_segments = set(), set(), set(), []

    for peak in h3_peaks:
        chrom, start, end = peak
        linked_promoters = promoters.query(chrom, start, end)
        if not linked_promoters:
            continue
        tss_h3_count += 1
        has_atac_full = atac.overlaps(chrom, start, end)
        if has_atac_full:
            atac_full_count += 1
        if has_atac_full and h3k4.overlaps(chrom, start, end):
            full_peaks.add(peak)
            for promoter in linked_promoters:
                full_genes.add(promoter["gene_id"])
                full_assoc.add((peak, promoter["gene_id"]))

        for promoter in linked_promoters:
            segment = (chrom, max(start, promoter["start"]), min(end, promoter["end"]))
            if segment[1] >= segment[2]:
                continue
            if atac.overlaps(*segment) and h3k4.overlaps(*segment):
                local_peaks.add(peak)
                local_genes.add(promoter["gene_id"])
                local_assoc.add((peak, promoter["gene_id"]))
                clipped.append(segment)
            exact_common = common_index.overlap_segments(*segment)
            if exact_common:
                common_peaks.add(peak)
                common_genes.add(promoter["gene_id"])
                common_assoc.add((peak, promoter["gene_id"]))
                common_segments.extend(exact_common)

    methods = {
        "P1_full_H3_peak_support": {
            "description": "H3 overlaps TSS; ATAC and H3K4me3 each overlap the full H3 peak; retain full H3 peak",
            "source_peaks": full_peaks,
            "genes": full_genes,
            "associations": full_assoc,
            "intervals": list(full_peaks),
        },
        "P2_TSS_local_support_full_H3": {
            "description": "ATAC and H3K4me3 each overlap H3-within-the-linked-TSS-window; retain full H3 peak",
            "source_peaks": local_peaks,
            "genes": local_genes,
            "associations": local_assoc,
            "intervals": list(local_peaks),
        },
        "P3_TSS_local_support_clipped": {
            "description": "same support as P2; retain H3 intersection with the linked TSS window",
            "source_peaks": local_peaks,
            "genes": local_genes,
            "associations": local_assoc,
            "intervals": clipped,
        },
        "P4_four_way_common_intersection": {
            "description": "retain exact common bases of H3, linked TSS window, ATAC and H3K4me3",
            "source_peaks": common_peaks,
            "genes": common_genes,
            "associations": common_assoc,
            "intervals": common_segments,
        },
    }
    steps = [{
        "time": time,
        "H3K27ac_total_unique_peaks": len(h3_peaks),
        "H3K27ac_overlap_any_TSS": tss_h3_count,
        "plus_ATAC_over_full_H3": atac_full_count,
        "plus_H3K4me3_over_full_H3_P1": len(full_peaks),
        "local_ATAC_and_H3K4me3_P2_P3": len(local_peaks),
        "four_way_common_P4": len(common_peaks),
    }]
    return steps, methods


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = ensure_file(args.config)
    cfg = json.loads(config_path.read_text(encoding="utf-8"))

    project = Path(cfg["project"]).resolve()
    source_run = Path(cfg["source_v27_run"]).resolve()
    output_root = Path(cfg["output_root"]).resolve()
    if not source_run.is_dir():
        raise FileNotFoundError(source_run)
    job = os.environ.get("SLURM_JOB_ID") or datetime.now(timezone.utc).strftime("manual_%Y%m%dT%H%M%SZ")
    output = output_root / f"export_{job}"
    output.mkdir(parents=True, exist_ok=False)
    code_dir = output / "code"
    code_dir.mkdir()
    bundle = Path(__file__).resolve().parent
    for name in ("export_and_count_definitions.py", "config.json", "run.slurm", "submit.sh", "README.md", "test_small.py", "test_integration.py"):
        source = bundle / name
        if source.is_file():
            shutil.copy2(source, code_dir / name)
    shutil.copy2(config_path, output / "config.used.json")
    print(f"results={output}", flush=True)

    chrom_names, lengths = read_fai(Path(cfg["reference_fai"]))
    genes = read_genes(Path(cfg["gtf"]), lengths)
    atac_paths = {time: source_run / relative for time, relative in cfg["atac_peak_files"].items()}
    histone_paths = {name: project / relative for name, relative in cfg["histone_peak_files"].items()}
    peaks = {
        "ATAC_0h": read_bed(atac_paths["0h"], lengths),
        "ATAC_6h": read_bed(atac_paths["6h"], lengths),
        **{name: read_bed(path, lengths) for name, path in histone_paths.items()},
    }
    indexes = {name: MergedIntervalIndex(rows) for name, rows in peaks.items()}
    h3k4_either = MergedIntervalIndex([*peaks["H3K4me3_0h"], *peaks["H3K4me3_6h"]])

    input_rows = []
    for label, path in {
        "GTF": Path(cfg["gtf"]),
        "reference_FAI": Path(cfg["reference_fai"]),
        "ATAC_0h": atac_paths["0h"],
        "ATAC_6h": atac_paths["6h"],
        **histone_paths,
    }.items():
        stat = ensure_file(path).stat()
        input_rows.append({"label": label, "path": str(path.resolve()), "size_bytes": stat.st_size})
    write_tsv(output / "00_inputs.tsv", input_rows, ["label", "path", "size_bytes"])

    bed_root = output / "BED"
    promoter_steps, promoter_counts, promoter_consensus = [], [], []
    promoter_interval_sets: dict[tuple[int, str, str], list[tuple[str, int, int]]] = {}
    for half_width in [int(value) for value in cfg["tss_half_widths_bp"]]:
        pindex = PromoterIndex(promoter_windows(genes, half_width, lengths))
        method_by_time = {}
        for time in ("0h", "6h"):
            steps, methods = promoter_definitions(
                time,
                peaks[f"H3K27ac_{time}"],
                pindex,
                indexes[f"ATAC_{time}"],
                indexes[f"H3K4me3_{time}"],
            )
            for row in steps:
                row["TSS_half_width_bp"] = half_width
                row["TSS_total_window_bp"] = half_width * 2 + 1
                promoter_steps.append(row)
            method_by_time[time] = methods
            for method, values in methods.items():
                summary = interval_summary(values["intervals"])
                promoter_counts.append({
                    "TSS_half_width_bp": half_width,
                    "TSS_total_window_bp": half_width * 2 + 1,
                    "time": time,
                    "method": method,
                    "source_H3K27ac_peaks": len(values["source_peaks"]),
                    "peak_gene_associations": len(values["associations"]),
                    "unique_genes": len(values["genes"]),
                    **summary,
                    "description": values["description"],
                })
                promoter_interval_sets[(half_width, method, time)] = values["intervals"]
                method_dir = bed_root / "promoter" / f"TSS_pm{half_width}bp" / method
                written = write_bed6(
                    method_dir / f"promoter_{time}.bed",
                    values["intervals"],
                    chrom_names,
                    f"PROM_{method}_{time}",
                )
                if written != summary["coordinate_intervals_unique"]:
                    raise RuntimeError("Promoter BED row count differs from summary")
        for method in method_by_time["0h"]:
            components = consensus_components({
                "0h": method_by_time["0h"][method]["intervals"],
                "6h": method_by_time["6h"][method]["intervals"],
            })
            component_lengths = [row["end"] - row["start"] for row in components]
            states = Counter(row["state"] for row in components)
            promoter_consensus.append({
                "TSS_half_width_bp": half_width,
                "method": method,
                "consensus_regions": len(components),
                "0h_only": states["0h_only"],
                "6h_only": states["6h_only"],
                "both_supported": states["both"],
                "consensus_total_bp": sum(component_lengths),
                "consensus_median_bp": percentile(component_lengths, 0.50),
                "consensus_p95_bp": percentile(component_lengths, 0.95),
                "consensus_max_bp": max(component_lengths) if component_lengths else 0,
            })
            method_dir = bed_root / "promoter" / f"TSS_pm{half_width}bp" / method
            written = write_consensus_bed6(
                method_dir / "promoter_consensus_0h_6h.bed",
                components,
                chrom_names,
                f"PROM_{method}_CONS",
            )
            if written != len(components):
                raise RuntimeError("Promoter consensus BED row count differs from summary")

    write_tsv(
        output / "01_promoter_stepwise_counts.tsv",
        promoter_steps,
        ["TSS_half_width_bp", "TSS_total_window_bp", "time", "H3K27ac_total_unique_peaks",
         "H3K27ac_overlap_any_TSS", "plus_ATAC_over_full_H3", "plus_H3K4me3_over_full_H3_P1",
         "local_ATAC_and_H3K4me3_P2_P3", "four_way_common_P4"],
    )
    promoter_fields = [
        "TSS_half_width_bp", "TSS_total_window_bp", "time", "method", "source_H3K27ac_peaks",
        "peak_gene_associations", "unique_genes", "coordinate_intervals_unique",
        "regions_after_true_overlap_merge", "union_total_bp", "union_median_bp", "union_p95_bp",
        "union_max_bp", "description",
    ]
    write_tsv(output / "02_promoter_definition_counts.tsv", promoter_counts, promoter_fields)
    write_tsv(
        output / "02_promoter_consensus_counts.tsv",
        promoter_consensus,
        ["TSS_half_width_bp", "method", "consensus_regions", "0h_only", "6h_only", "both_supported",
         "consensus_total_bp", "consensus_median_bp", "consensus_p95_bp", "consensus_max_bp"],
    )

    enhancer_steps, enhancer_counts, enhancer_consensus = [], [], []
    enhancer_interval_sets: dict[tuple[int, str, str], list[tuple[str, int, int]]] = {}
    for half_width in [int(value) for value in cfg["tss_half_widths_bp"]]:
        promoter_index = MergedIntervalIndex([
            (row["chrom"], row["start"], row["end"])
            for row in promoter_windows(genes, half_width, lengths)
        ])
        for h3k4_mode in ("same_time", "either_time_v27"):
            for time in ("0h", "6h"):
                h3 = peaks[f"H3K27ac_{time}"]
                atac_supported = [row for row in h3 if indexes[f"ATAC_{time}"].overlaps(*row)]
                after_tss = [row for row in atac_supported if not promoter_index.overlaps(*row)]
                h3k4_index = indexes[f"H3K4me3_{time}"] if h3k4_mode == "same_time" else h3k4_either
                final = [row for row in after_tss if not h3k4_index.overlaps(*row)]
                summary = interval_summary(final)
                enhancer_steps.append({
                    "TSS_half_width_bp": half_width,
                    "TSS_total_window_bp": half_width * 2 + 1,
                    "H3K4me3_exclusion": h3k4_mode,
                    "time": time,
                    "H3K27ac_total_unique_peaks": len(h3),
                    "plus_same_time_ATAC": len(atac_supported),
                    "excluded_TSS": len(atac_supported) - len(after_tss),
                    "after_TSS_exclusion": len(after_tss),
                    "excluded_H3K4me3": len(after_tss) - len(final),
                    "final_H3K27ac_peaks": len(final),
                })
                enhancer_counts.append({
                    "TSS_half_width_bp": half_width,
                    "TSS_total_window_bp": half_width * 2 + 1,
                    "H3K4me3_exclusion": h3k4_mode,
                    "time": time,
                    "source_H3K27ac_peaks": len(final),
                    **summary,
                })
                enhancer_interval_sets[(half_width, h3k4_mode, time)] = final
                definition_dir = bed_root / "enhancer" / f"TSS_pm{half_width}bp" / f"H3K4me3_{h3k4_mode}"
                written = write_bed6(
                    definition_dir / f"enhancer_{time}.bed",
                    final,
                    chrom_names,
                    f"ENH_{h3k4_mode}_{time}",
                )
                if written != summary["coordinate_intervals_unique"]:
                    raise RuntimeError("Enhancer BED row count differs from summary")
            components = consensus_components({
                "0h": enhancer_interval_sets[(half_width, h3k4_mode, "0h")],
                "6h": enhancer_interval_sets[(half_width, h3k4_mode, "6h")],
            })
            states = Counter(row["state"] for row in components)
            lengths_bp = [row["end"] - row["start"] for row in components]
            enhancer_consensus.append({
                "TSS_half_width_bp": half_width,
                "H3K4me3_exclusion": h3k4_mode,
                "consensus_regions": len(components),
                "0h_only": states["0h_only"],
                "6h_only": states["6h_only"],
                "both_supported": states["both"],
                "consensus_total_bp": sum(lengths_bp),
                "consensus_median_bp": percentile(lengths_bp, 0.50),
                "consensus_p95_bp": percentile(lengths_bp, 0.95),
                "consensus_max_bp": max(lengths_bp) if lengths_bp else 0,
            })
            definition_dir = bed_root / "enhancer" / f"TSS_pm{half_width}bp" / f"H3K4me3_{h3k4_mode}"
            written = write_consensus_bed6(
                definition_dir / "enhancer_consensus_0h_6h.bed",
                components,
                chrom_names,
                f"ENH_{h3k4_mode}_CONS",
            )
            if written != len(components):
                raise RuntimeError("Enhancer consensus BED row count differs from summary")

    write_tsv(
        output / "03_enhancer_stepwise_counts.tsv",
        enhancer_steps,
        ["TSS_half_width_bp", "TSS_total_window_bp", "H3K4me3_exclusion", "time",
         "H3K27ac_total_unique_peaks", "plus_same_time_ATAC", "excluded_TSS", "after_TSS_exclusion",
         "excluded_H3K4me3", "final_H3K27ac_peaks"],
    )
    write_tsv(
        output / "04_enhancer_definition_counts.tsv",
        enhancer_counts,
        ["TSS_half_width_bp", "TSS_total_window_bp", "H3K4me3_exclusion", "time",
         "source_H3K27ac_peaks", "coordinate_intervals_unique", "regions_after_true_overlap_merge",
         "union_total_bp", "union_median_bp", "union_p95_bp", "union_max_bp"],
    )
    write_tsv(
        output / "04_enhancer_consensus_counts.tsv",
        enhancer_consensus,
        ["TSS_half_width_bp", "H3K4me3_exclusion", "consensus_regions", "0h_only", "6h_only",
         "both_supported", "consensus_total_bp", "consensus_median_bp", "consensus_p95_bp",
         "consensus_max_bp"],
    )

    expected_membership = read_tsv(source_run / "03_H3K27ac_consensus_membership.tsv")
    expected_sources = {
        (canon(row["chrom"]), int(row["peak_start"]), int(row["peak_end"]), row["source_time"])
        for row in expected_membership
    }
    reproduced_sources = {
        (*interval, time)
        for time in ("0h", "6h")
        for interval in enhancer_interval_sets[(2000, "either_time_v27", time)]
    }
    matrix = read_tsv(source_run / "03_CRE_signal_matrix.tsv")
    expected_consensus = {(canon(row["chrom"]), int(row["start"]), int(row["end"])) for row in matrix}
    reproduced_components = consensus_components({
        time: enhancer_interval_sets[(2000, "either_time_v27", time)] for time in ("0h", "6h")
    })
    reproduced_consensus = {(row["chrom"], row["start"], row["end"]) for row in reproduced_components}
    qc_rows = [
        {"comparison": "source_time_H3K27ac_membership", "expected": len(expected_sources),
         "reproduced": len(reproduced_sources), "missing": len(expected_sources - reproduced_sources),
         "extra": len(reproduced_sources - expected_sources),
         "status": "PASS" if expected_sources == reproduced_sources else "FAIL"},
        {"comparison": "consensus_CRE_coordinates", "expected": len(expected_consensus),
         "reproduced": len(reproduced_consensus), "missing": len(expected_consensus - reproduced_consensus),
         "extra": len(reproduced_consensus - expected_consensus),
         "status": "PASS" if expected_consensus == reproduced_consensus else "FAIL"},
    ]
    write_tsv(output / "05_v27_reproduction_QC.tsv", qc_rows, ["comparison", "expected", "reproduced", "missing", "extra", "status"])
    if any(row["status"] != "PASS" for row in qc_rows):
        raise RuntimeError("Current v27 enhancer definition was not exactly reproduced; see 05_v27_reproduction_QC.tsv")

    bed_manifest = []
    for bed_path in sorted(bed_root.rglob("*.bed")):
        with bed_path.open(encoding="utf-8") as handle:
            row_count = sum(1 for line in handle if line.strip() and not line.startswith(("#", "track", "browser")))
        bed_manifest.append({
            "relative_path": bed_path.relative_to(output).as_posix(),
            "row_count": row_count,
            "size_bytes": bed_path.stat().st_size,
        })
    write_tsv(output / "00_BED_manifest.tsv", bed_manifest, ["relative_path", "row_count", "size_bytes"])

    notes = f"""候选promoter/enhancer定义敏感性导出与统计说明
==============================================

本任务直接导出不同定义下的0 h、6 h及0 h/6 h consensus BED，然后依据实际导出的区间对象生成统计表；不读取BAM，也不选定唯一的最终生物学定义。

TSS窗口：配置中的500、1500和2000均解释为TSS上下游各500、1500和2000 bp，对应完整窗口长度1001、3001和4001 bp。所有重叠均采用BED半开区间下至少1 bp真实重叠。

Promoter四种比较方法：
P1_full_H3_peak_support：H3K27ac peak与TSS窗口重叠；同时间ATAC和H3K4me3各自在完整H3K27ac peak内有至少1 bp支持；输出坐标仍为完整H3K27ac peak。这是较宽松且坐标规则最一致的方法。
P2_TSS_local_support_full_H3：ATAC和H3K4me3必须各自落在“H3K27ac peak与对应gene TSS窗口的交集”内；通过后仍保留完整H3K27ac peak。这是建议优先评估的方法。
P3_TSS_local_support_clipped：支持判据与P2相同，但坐标裁剪为H3K27ac peak与TSS窗口的交集。P2与P3的来源H3K27ac peak数相同，长度与输出区间数可能不同。
P4_four_way_common_intersection：仅保留H3K27ac、对应TSS窗口、ATAC及H3K4me3四者共同覆盖的碱基。这是最严格的方法，可能明显减少并碎片化区间。

Enhancer比较方法：
共同规则为完整H3K27ac peak与同时间ATAC peak至少重叠1 bp，保留完整H3K27ac peak；与TSS窗口或指定H3K4me3集合有至少1 bp重叠就排除整个peak。
same_time：只用同时间H3K4me3排除。
either_time_v27：用0 h和6 h H3K4me3并集排除；其中TSS上下游2 kb版本应精确复现v27现有增强子定义。

输出：
BED/promoter/TSS_pm*/P*/promoter_0h.bed、promoter_6h.bed和promoter_consensus_0h_6h.bed：四种promoter规则的时间点及共识坐标。
BED/enhancer/TSS_pm*/H3K4me3_*/enhancer_0h.bed、enhancer_6h.bed和enhancer_consensus_0h_6h.bed：enhancer规则的时间点及共识坐标。
00_BED_manifest.tsv：逐个BED记录相对路径、实际行数和文件大小。
01_promoter_stepwise_counts.tsv：promoter逐步过滤数量。
02_promoter_definition_counts.tsv：四种promoter方法在各时间点的来源peak数、gene数、区间数与长度。
02_promoter_consensus_counts.tsv：跨0 h/6 h真实重叠合并后的promoter consensus数量及时间点支持状态。
03_enhancer_stepwise_counts.tsv：enhancer逐步过滤数量。
04_enhancer_definition_counts.tsv：各TSS窗口及H3K4me3排除模式的时间点结果。
04_enhancer_consensus_counts.tsv：对应的consensus数量及0h_only、6h_only、both数量。
05_v27_reproduction_QC.tsv：确认TSS上下游2 kb + either-time H3K4me3的结果与v27来源peak和consensus坐标完全一致。

BED文件均为无表头BED6：第1列染色体；第2列0-based、包含的起点；第3列0-based、不包含的终点；第4列本次生成的区间ID；第5列固定为0；第6列为“.”。0 h和6 h文件保留各定义产生的唯一来源区间，不额外合并相互重叠的peak；consensus文件跨两个时间点仅合并真实重叠，首尾相接不合并。

输入ATAC peaks：{atac_paths['0h']}；{atac_paths['6h']}
输入H3K27ac/H3K4me3：{histone_paths}
输入GTF：{Path(cfg['gtf']).resolve()}
"""
    (output / "说明_定义与输出列.txt").write_text(notes, encoding="utf-8")
    completed = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "results": str(output),
        "bed_root": str(bed_root),
        "tss_half_widths_bp": cfg["tss_half_widths_bp"],
        "v27_reproduction": {row["comparison"]: row["status"] for row in qc_rows},
    }
    (output / "completed.json").write_text(json.dumps(completed, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(completed, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
