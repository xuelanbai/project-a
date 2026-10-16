#!/usr/bin/env python3
"""Classify P1 promoter-associated H3K27ac regions as strength/weaken/stable.

P1 identity is time specific: a complete H3K27ac peak must overlap a gene
TSS +/- configured width, and same-time ATAC and H3K4me3 peaks must each
overlap that complete H3K27ac peak by at least 1 bp.  The complete H3K27ac
peak remains the promoter coordinate.

The 0 h and 6 h P1 coordinates are combined into one true-overlap consensus
set.  H3K27ac and ATAC BAMs are then counted over exactly the same consensus
coordinates.  Strength/weaken is defined only by H3K27ac fold change, active
endpoint P1 support, and the active-endpoint H3K27ac minimum count.  ATAC
fold change is annotation only; no ATAC read-count or concordance filter is
applied.
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
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
import sys


@dataclass(frozen=True)
class GenePromoter:
    chrom: str
    window_start: int
    window_end: int
    gene_id: str
    gene_symbol: str
    tss: int
    strand: str


def canonical_chrom(value: str) -> str:
    value = str(value).strip()
    return value[3:] if value.startswith("chr") else value


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    return (int(chrom), "") if chrom.isdigit() else (10_000, chrom)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix == ".gz" else path.open(
        encoding="utf-8-sig"
    )


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields and not key.startswith("_"):
                    fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class IntervalIndex:
    """Small chromosome/bin index returning each overlapping value once."""

    def __init__(self, bin_size: int = 100_000):
        self.bin_size = bin_size
        self.bins: dict[str, dict[int, list[int]]] = defaultdict(lambda: defaultdict(list))
        self.items: list[tuple[int, int, object]] = []

    def add(self, chrom: str, start: int, end: int, value: object) -> None:
        chrom = canonical_chrom(chrom)
        item_id = len(self.items)
        self.items.append((start, end, value))
        for bin_id in range(start // self.bin_size, (end - 1) // self.bin_size + 1):
            self.bins[chrom][bin_id].append(item_id)

    def overlaps(self, chrom: str, start: int, end: int):
        seen: set[int] = set()
        bins = self.bins.get(canonical_chrom(chrom), {})
        for bin_id in range(start // self.bin_size, (end - 1) // self.bin_size + 1):
            for item_id in bins.get(bin_id, ()): 
                if item_id in seen:
                    continue
                seen.add(item_id)
                item_start, item_end, value = self.items[item_id]
                if item_start < end and start < item_end:
                    yield value

    def has_overlap(self, chrom: str, start: int, end: int) -> bool:
        return next(iter(self.overlaps(chrom, start, end)), None) is not None


def read_fai(path: Path, allowed: set[str]) -> tuple[dict[str, int], dict[str, str]]:
    lengths: dict[str, int] = {}
    names: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: malformed FASTA index")
            chrom = canonical_chrom(fields[0])
            if chrom not in allowed:
                continue
            if chrom in lengths:
                raise ValueError(f"{path}:{line_no}: duplicate chromosome alias {chrom}")
            lengths[chrom] = int(fields[1])
            names[chrom] = fields[0]
    if set(lengths) != allowed:
        raise ValueError(f"FASTA main chromosome mismatch: {sorted(set(lengths) ^ allowed)}")
    return lengths, names


def parse_gtf_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in text.split(";"):
        field = field.strip()
        if not field or " " not in field:
            continue
        key, value = field.split(None, 1)
        result[key] = value.strip().strip('"')
    return result


def read_gene_promoters(gtf: Path, lengths: dict[str, int], half_width: int) -> list[GenePromoter]:
    rows: set[GenePromoter] = set()
    with open_text(gtf) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gtf}:{line_no}: expected nine GTF columns")
            if fields[2] != "gene":
                continue
            chrom = canonical_chrom(fields[0])
            if chrom not in lengths:
                continue
            if fields[6] not in {"+", "-"}:
                raise ValueError(f"{gtf}:{line_no}: invalid gene strand")
            start, end = int(fields[3]), int(fields[4])
            if not (1 <= start <= end <= lengths[chrom]):
                raise ValueError(f"{gtf}:{line_no}: gene outside chromosome")
            attributes = parse_gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id", "")
            if not gene_id:
                raise ValueError(f"{gtf}:{line_no}: gene feature lacks gene_id")
            gene_symbol = attributes.get("gene_name", gene_id)
            tss = start - 1 if fields[6] == "+" else end - 1
            rows.add(GenePromoter(
                chrom=chrom,
                window_start=max(0, tss - half_width),
                window_end=min(lengths[chrom], tss + half_width + 1),
                gene_id=gene_id,
                gene_symbol=gene_symbol,
                tss=tss,
                strand=fields[6],
            ))
    if not rows:
        raise ValueError(f"No main-chromosome gene promoters parsed from {gtf}")
    return sorted(rows, key=lambda x: (chrom_sort_key(x.chrom), x.window_start, x.window_end, x.gene_id))


def read_bed(path: Path, allowed: set[str], label: str) -> tuple[list[tuple[str, int, int]], list[dict]]:
    intervals: set[tuple[str, int, int]] = set()
    column_counts: Counter = Counter()
    with path.open(encoding="utf-8-sig", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            column_counts[len(fields)] += 1
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: fewer than three BED columns")
            chrom = canonical_chrom(fields[0])
            start, end = int(fields[1]), int(fields[2])
            if not 0 <= start < end:
                raise ValueError(f"{path}:{line_no}: invalid BED interval")
            if chrom in allowed:
                intervals.add((chrom, start, end))
    if not intervals:
        raise ValueError(f"No main-chromosome intervals: {path}")
    qc = [dict(label=label, path=str(path), columns=n, records=count) for n, count in sorted(column_counts.items())]
    return sorted(intervals, key=lambda x: (chrom_sort_key(x[0]), x[1], x[2])), qc


def index_intervals(rows: list[tuple[str, int, int]]) -> IntervalIndex:
    index = IntervalIndex()
    for chrom, start, end in rows:
        index.add(chrom, start, end, True)
    return index


def build_p1(
    time: str,
    h3_peaks: list[tuple[str, int, int]],
    atac: IntervalIndex,
    h3k4me3: IntervalIndex,
    genes: IntervalIndex,
) -> tuple[list[dict], list[dict], dict]:
    coordinates: set[tuple[str, int, int]] = set()
    associations: set[tuple] = set()
    qc: Counter = Counter()
    for chrom, start, end in h3_peaks:
        qc["input_H3K27ac_peaks"] += 1
        linked_genes = list(genes.overlaps(chrom, start, end))
        if not linked_genes:
            qc["excluded_no_TSS_window_overlap"] += 1
            continue
        qc["H3K27ac_peaks_overlapping_TSS"] += 1
        if not atac.has_overlap(chrom, start, end):
            qc["excluded_no_same_time_ATAC_peak_overlap"] += 1
            continue
        if not h3k4me3.has_overlap(chrom, start, end):
            qc["excluded_no_same_time_H3K4me3_peak_overlap"] += 1
            continue
        coordinates.add((chrom, start, end))
        for gene in linked_genes:
            associations.add((chrom, start, end, gene.gene_id, gene.gene_symbol, gene.tss, gene.strand))
    ordered_coordinates = [
        dict(chrom=chrom, start=start, end=end, source_time=time, source_peak_id=f"P1_{time}_{number:07d}")
        for number, (chrom, start, end) in enumerate(
            sorted(coordinates, key=lambda x: (chrom_sort_key(x[0]), x[1], x[2])), 1
        )
    ]
    peak_ids = {(row["chrom"], row["start"], row["end"]): row["source_peak_id"] for row in ordered_coordinates}
    ordered_associations = []
    for chrom, start, end, gene_id, gene_symbol, tss, strand in sorted(
        associations, key=lambda x: (chrom_sort_key(x[0]), x[1], x[2], x[3], x[5])
    ):
        ordered_associations.append(dict(
            source_time=time,
            source_peak_id=peak_ids[(chrom, start, end)],
            chrom=chrom,
            start=start,
            end=end,
            gene_id=gene_id,
            gene_symbol=gene_symbol,
            gene_tss=tss,
            gene_strand=strand,
        ))
    qc["P1_unique_intervals"] = len(ordered_coordinates)
    qc["P1_peak_gene_associations"] = len(ordered_associations)
    return ordered_coordinates, ordered_associations, dict(qc)


def build_consensus(peaks_by_time: dict[str, list[dict]]) -> tuple[list[dict], list[dict]]:
    ordered = sorted(
        [row for rows in peaks_by_time.values() for row in rows],
        key=lambda x: (chrom_sort_key(x["chrom"]), x["start"], x["end"], x["source_time"]),
    )
    components: list[list[dict]] = []
    current: list[dict] = []
    current_chrom = ""
    current_end = -1
    for row in ordered:
        if not current or row["chrom"] != current_chrom or row["start"] >= current_end:
            if current:
                components.append(current)
            current = [row]
            current_chrom = row["chrom"]
            current_end = row["end"]
        else:
            current.append(row)
            current_end = max(current_end, row["end"])
    if current:
        components.append(current)

    regions, membership = [], []
    for number, members in enumerate(components, 1):
        region_id = f"P1_PROMOTER_{number:07d}"
        start, end = min(x["start"] for x in members), max(x["end"] for x in members)
        region = dict(
            region_id=region_id,
            chrom=members[0]["chrom"],
            start=start,
            end=end,
            length_bp=end - start,
            source_times=";".join(sorted({x["source_time"] for x in members})),
            member_peak_count=len(members),
            member_peak_count_0h=sum(x["source_time"] == "0h" for x in members),
            member_peak_count_6h=sum(x["source_time"] == "6h" for x in members),
        )
        region["supported_0h"] = int(region["member_peak_count_0h"] > 0)
        region["supported_6h"] = int(region["member_peak_count_6h"] > 0)
        regions.append(region)
        for member in members:
            membership.append(dict(
                region_id=region_id,
                consensus_chrom=region["chrom"],
                consensus_start=start,
                consensus_end=end,
                source_time=member["source_time"],
                source_peak_id=member["source_peak_id"],
                source_peak_start=member["start"],
                source_peak_end=member["end"],
            ))
    return regions, membership


def attach_gene_associations(membership: list[dict], associations_by_time: dict[str, list[dict]]) -> list[dict]:
    source_to_region = {row["source_peak_id"]: row["region_id"] for row in membership}
    grouped: dict[tuple, set[str]] = defaultdict(set)
    for time, rows in associations_by_time.items():
        for row in rows:
            key = (
                source_to_region[row["source_peak_id"]], row["gene_id"], row["gene_symbol"],
                row["gene_tss"], row["gene_strand"], row["chrom"],
            )
            grouped[key].add(time)
    output = []
    for key, times in sorted(grouped.items(), key=lambda x: (x[0][5], x[0][3], x[0][1], x[0][0])):
        region_id, gene_id, symbol, tss, strand, chrom = key
        output.append(dict(
            region_id=region_id,
            chrom=chrom,
            gene_id=gene_id,
            gene_symbol=symbol,
            gene_tss=tss,
            gene_strand=strand,
            P1_support_times=";".join(sorted(times)),
            P1_supported_0h=int("0h" in times),
            P1_supported_6h=int("6h" in times),
        ))
    return output


def write_bed(path: Path, rows: list[dict], bam_chrom_names: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for row in sorted(rows, key=lambda x: (chrom_sort_key(x["chrom"]), x["start"], x["end"], x.get("region_id", x.get("source_peak_id", "")))):
            identifier = row.get("region_id", row.get("source_peak_id", "NA"))
            handle.write(f"{bam_chrom_names[row['chrom']]}\t{row['start']}\t{row['end']}\t{identifier}\t0\t.\n")


def run_command(command: list, stdout: Path, stderr: Path, command_log: Path) -> None:
    command = [str(item) for item in command]
    with command_log.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(command, ensure_ascii=False) + "\n")
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        subprocess.run(command, stdout=out, stderr=err, check=True)


def multicov(bedtools: Path, bed: Path, bam_paths: list[Path], labels: list[str], output: Path,
             command_log: Path, mapq: int) -> dict[str, dict[str, int]]:
    raw = output.with_suffix(".raw.tsv")
    run_command(
        [bedtools, "multicov", "-bed", bed, "-bams", *bam_paths, "-q", str(mapq), "-p"],
        raw,
        output.with_suffix(".stderr.log"),
        command_log,
    )
    counts: dict[str, dict[str, int]] = {}
    with raw.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            bed_columns = len(fields) - len(labels)
            if bed_columns < 4:
                raise ValueError(f"{raw}:{line_no}: fewer than four BED columns")
            identifier = fields[3]
            if identifier in counts:
                raise ValueError(f"{raw}:{line_no}: duplicate interval ID {identifier}")
            counts[identifier] = {label: int(value) for label, value in zip(labels, fields[-len(labels):])}
    if not counts:
        raise ValueError(f"Empty BEDTools multicov result: {raw}")
    return counts


def validate_bams(bam_paths: dict[str, Path], allowed: set[str], lengths: dict[str, int]):
    import pysam

    maps: dict[str, dict[str, str]] = {}
    for label, path in bam_paths.items():
        with pysam.AlignmentFile(str(path), "rb") as bam:
            if not bam.has_index():
                raise ValueError(f"BAM is not indexed: {path}")
            mapping = {canonical_chrom(name): name for name in bam.references if canonical_chrom(name) in allowed}
            if set(mapping) != allowed:
                raise ValueError(f"{label}: BAM main chromosome mismatch: {sorted(set(mapping) ^ allowed)}")
            for chrom in allowed:
                if bam.get_reference_length(mapping[chrom]) != lengths[chrom]:
                    raise ValueError(f"{label}: reference length differs for chromosome {chrom}")
            maps[label] = mapping
    first = maps[next(iter(maps))]
    if any(mapping != first for mapping in maps.values()):
        raise ValueError("BAM chromosome naming differs across input files")
    return first


def classify(region: dict, fold: float, min_h3_tags: int) -> tuple[str, str]:
    cutoff = math.log2(fold)
    value = float(region["H3K27ac_log2FC"])
    if value >= cutoff:
        if not region["supported_6h"]:
            return "stable", "strength_like_but_no_6h_P1_support"
        if region["H3K27ac_6h_raw_tags"] < min_h3_tags:
            return "stable", "strength_like_but_low_6h_H3K27ac_count"
        return "strength", "H3K27ac_at_least_1p5_fold_higher_with_6h_P1_support"
    if value <= -cutoff:
        if not region["supported_0h"]:
            return "stable", "weaken_like_but_no_0h_P1_support"
        if region["H3K27ac_0h_raw_tags"] < min_h3_tags:
            return "stable", "weaken_like_but_low_0h_H3K27ac_count"
        return "weaken", "H3K27ac_at_least_1p5_fold_lower_with_0h_P1_support"
    return "stable", "H3K27ac_within_1p5_fold_threshold"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Submit real BAM processing as a Linux Slurm job; do not run it on the login node.")

    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    if float(cfg["fold_threshold"]) != 1.5:
        raise ValueError("This version fixes the H3K27ac fold threshold at 1.5")
    if cfg.get("apply_ATAC_raw_count_filter", False):
        raise ValueError("This version intentionally disables the ATAC raw-count filter")

    allowed = {canonical_chrom(x) for x in cfg["main_chromosomes"]}
    fasta = Path(cfg["reference_fasta"])
    fai = Path(str(fasta) + ".fai")
    gtf = Path(cfg["gtf"])
    peak_paths = {label: Path(path) for label, path in cfg["peak_inputs"].items()}
    bam_paths = {label: Path(path) for label, path in cfg["bam_inputs"].items()}
    bedtools = Path(cfg["bedtools"])
    required = [fasta, fai, gtf, bedtools, *peak_paths.values(), *bam_paths.values()]
    for path in required:
        if not path.is_absolute() or not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    for path in bam_paths.values():
        if not Path(str(path) + ".bai").is_file():
            raise FileNotFoundError(Path(str(path) + ".bai"))

    lengths, _ = read_fai(fai, allowed)
    bam_chrom_names = validate_bams(bam_paths, allowed, lengths)
    output = Path(cfg["output_root"]) / f"run_{os.environ['SLURM_JOB_ID']}"
    output.mkdir(parents=True, exist_ok=False)
    for directory in ("P1", "counts", "states", "logs", "code"):
        (output / directory).mkdir()
    command_log = output / "logs/commands.jsonl"
    for name in ("run.py", "config.json", "README.md", "run.slurm", "submit.sh", "check_env.sh", "test_small.py"):
        shutil.copy2(Path(__file__).with_name(name), output / "code" / name)
    print(f"results={output}", flush=True)

    gene_promoters = read_gene_promoters(gtf, lengths, int(cfg["promoter_half_width"]))
    gene_index = IntervalIndex()
    for gene in gene_promoters:
        gene_index.add(gene.chrom, gene.window_start, gene.window_end, gene)

    bed_rows: dict[str, list[tuple[str, int, int]]] = {}
    schema_rows: list[dict] = []
    for label, path in peak_paths.items():
        bed_rows[label], qc = read_bed(path, allowed, label)
        schema_rows.extend(qc)
    write_tsv(output / "01_input_peak_schema.tsv", schema_rows)

    peaks_by_time: dict[str, list[dict]] = {}
    associations_by_time: dict[str, list[dict]] = {}
    definition_qc: list[dict] = []
    for time in ("0h", "6h"):
        coordinates, associations, qc = build_p1(
            time=time,
            h3_peaks=bed_rows[f"H3K27ac_{time}"],
            atac=index_intervals(bed_rows[f"ATAC_{time}"]),
            h3k4me3=index_intervals(bed_rows[f"H3K4me3_{time}"]),
            genes=gene_index,
        )
        peaks_by_time[time] = coordinates
        associations_by_time[time] = associations
        write_bed(output / f"P1/02_P1_promoters_{time}.bed", coordinates, bam_chrom_names)
        write_tsv(output / f"P1/02_P1_promoters_{time}_gene_associations.tsv", associations)
        definition_qc.append(dict(time=time, **qc))
    write_tsv(output / "02_P1_definition_QC.tsv", definition_qc)

    regions, membership = build_consensus(peaks_by_time)
    if not regions:
        raise ValueError("No P1 consensus promoter regions")
    consensus_associations = attach_gene_associations(membership, associations_by_time)
    write_bed(output / "P1/03_P1_promoter_consensus_0h_6h.bed", regions, bam_chrom_names)
    write_tsv(output / "P1/03_P1_promoter_consensus_membership.tsv", membership)
    write_tsv(output / "P1/03_P1_promoter_consensus_gene_associations.tsv", consensus_associations)

    write_bed(output / "counts/P1_promoter_consensus.bed", regions, bam_chrom_names)
    main_rows = [
        dict(chrom=chrom, start=0, end=lengths[chrom], region_id=f"CHR_{chrom}")
        for chrom in sorted(allowed, key=chrom_sort_key)
    ]
    write_bed(output / "counts/main_chromosomes.bed", main_rows, bam_chrom_names)
    labels = ["H3K27ac_0h", "H3K27ac_6h", "ATAC_0h", "ATAC_6h"]
    library_counts = multicov(
        bedtools,
        output / "counts/main_chromosomes.bed",
        [bam_paths[label] for label in labels],
        labels,
        output / "counts/library_sizes.tsv",
        command_log,
        int(cfg["min_mapq"]),
    )
    library_totals = {label: sum(row[label] for row in library_counts.values()) for label in labels}
    if any(total <= 0 for total in library_totals.values()):
        raise ValueError(f"Invalid library totals: {library_totals}")
    region_counts = multicov(
        bedtools,
        output / "counts/P1_promoter_consensus.bed",
        [bam_paths[label] for label in labels],
        labels,
        output / "counts/P1_promoter_counts.tsv",
        command_log,
        int(cfg["min_mapq"]),
    )

    scale = float(cfg["normalization_total_tags"])
    pseudocount = float(cfg["normalized_pseudocount"])
    min_h3 = int(cfg["min_H3K27ac_raw_tags_active_state"])
    fold = float(cfg["fold_threshold"])
    for region in regions:
        counts = region_counts[region["region_id"]]
        for label in labels:
            region[f"{label}_raw_tags"] = counts[label]
            region[f"{label}_tags_per_10M"] = counts[label] / library_totals[label] * scale
        for assay in ("H3K27ac", "ATAC"):
            region[f"{assay}_log2FC"] = math.log2(
                (region[f"{assay}_6h_tags_per_10M"] + pseudocount)
                / (region[f"{assay}_0h_tags_per_10M"] + pseudocount)
            )
        atac_fc = region["ATAC_log2FC"]
        region["ATAC_stratum"] = "ATAC_up" if atac_fc > 0 else "ATAC_down" if atac_fc < 0 else "ATAC_unchanged"
        region["promoter_state"], region["classification_reason"] = classify(region, fold, min_h3)
        region["ATAC_raw_count_filter_applied"] = 0
        region["ATAC_direction_filter_applied"] = 0
    write_tsv(output / "04_P1_promoter_signal_and_state_matrix.tsv", regions)

    state_counts = Counter(row["promoter_state"] for row in regions)
    reason_counts = Counter((row["promoter_state"], row["classification_reason"]) for row in regions)
    count_rows = [
        dict(promoter_state=state, classification_reason=reason, count=count)
        for (state, reason), count in sorted(reason_counts.items())
    ]
    write_tsv(output / "05_P1_promoter_state_counts.tsv", count_rows)

    for state in ("strength", "weaken", "stable"):
        selected = [row for row in regions if row["promoter_state"] == state]
        write_bed(output / f"states/05_P1_{state}_promoters.bed", selected, bam_chrom_names)
        write_tsv(output / f"states/05_P1_{state}_promoters.tsv", selected)

    region_lookup = {row["region_id"]: row for row in regions}
    state_gene_rows = []
    for association in consensus_associations:
        region = region_lookup[association["region_id"]]
        state_gene_rows.append(dict(
            promoter_state=region["promoter_state"],
            classification_reason=region["classification_reason"],
            promoter_chrom=region["chrom"],
            promoter_start=region["start"],
            promoter_end=region["end"],
            H3K27ac_log2FC=region["H3K27ac_log2FC"],
            ATAC_log2FC=region["ATAC_log2FC"],
            ATAC_stratum=region["ATAC_stratum"],
            **association,
        ))
    write_tsv(output / "05_P1_promoter_state_gene_associations.tsv", state_gene_rows)

    metadata = {
        "completed": datetime.now(timezone.utc).isoformat(),
        "analysis_mode": "pooled_BAM_descriptive_fold_change_matching_final_enhancer_definition",
        "promoter_definition": "P1_full_H3K27ac_peak_overlaps_TSS_window_with_same_time_ATAC_and_H3K4me3_peak_support",
        "promoter_coordinate": "complete_H3K27ac_peak_then_true_overlap_0h_6h_consensus",
        "state_definition": "H3K27ac_6h_vs_0h_tags_per_10M_fold_change_1p5",
        "ATAC_policy": "same_time_peak_support_for_P1_identity_and_log2FC_annotation_only; no_ATAC_raw_count_or_direction_filter",
        "count_unit": "bedtools_multicov_proper_pair_alignments_not_fragment_counts",
        "normalization": "main_chromosome_alignment_counts_per_10M; no_length_division; no_median_centering",
        "library_totals": library_totals,
        "counts": {
            "P1_0h_intervals": len(peaks_by_time["0h"]),
            "P1_6h_intervals": len(peaks_by_time["6h"]),
            "P1_consensus_regions": len(regions),
            "strength": state_counts["strength"],
            "weaken": state_counts["weaken"],
            "stable": state_counts["stable"],
        },
        "input_sha256": {str(path): sha256(path) for path in [gtf, *peak_paths.values()]},
        "config": cfg,
    }
    (output / "run_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    (output / "completed.json").write_text(
        json.dumps({"completed": True, "results": str(output)}, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata["counts"], ensure_ascii=False), flush=True)
    print(f"completed=1\nresults={output}", flush=True)


if __name__ == "__main__":
    main()
