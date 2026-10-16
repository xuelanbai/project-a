#!/usr/bin/env python3
"""Build three explicit promoter definitions for the HiChIP target analysis.

Coordinates are zero-based, half-open.  P1 and P2 retain complete H3K27ac
peak coordinates, but preserve every qualifying peak-to-gene association so
that a peak near multiple TSSs is never silently assigned to all genes after
coordinate-only deduplication.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import intervals


TIMES = ("0h", "6h")
P1 = "P1_full_H3_peak_support"
P2 = "P2_TSS_local_support_full_H3"
TSS = "TSS_pm2kb"


def _chrom_key(chrom: str) -> tuple[int, str]:
    if chrom.isdigit():
        return int(chrom), ""
    return 10_000, chrom


def load_peak_intervals(path: str | Path, allowed_chromosomes: set[str]) -> list[tuple[str, int, int]]:
    """Read the first three BED columns and return unique sorted intervals."""
    result = set()
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: expected at least three BED columns")
            chrom, start, end = intervals._interval(fields[0], fields[1], fields[2], f"{path}:{line_no}")
            if chrom in allowed_chromosomes:
                result.add((chrom, start, end))
    if not result:
        raise ValueError(f"No usable peak intervals: {path}")
    return sorted(result, key=lambda row: (_chrom_key(row[0]), row[1], row[2]))


def make_boolean_index(rows: list[tuple[str, int, int]]) -> intervals.IntervalIndex:
    index = intervals.IntervalIndex()
    for chrom, start, end in rows:
        index.add(chrom, start, end, True)
    return index


def _has_overlap(index: intervals.IntervalIndex, chrom: str, start: int, end: int) -> bool:
    return next(iter(index.overlaps(chrom, start, end)), None) is not None


def _promoter_values(index: intervals.IntervalIndex) -> list[intervals.Promoter]:
    return sorted(
        {value for _, _, value in index._items},
        key=lambda row: (_chrom_key(row.chrom), row.start, row.end, row.gene_id, row.tss),
    )


def _association_index(values: set[intervals.Promoter]) -> intervals.IntervalIndex:
    index = intervals.IntervalIndex()
    for value in sorted(values, key=lambda row: (_chrom_key(row.chrom), row.start, row.end, row.gene_id, row.tss)):
        index.add(value.chrom, value.start, value.end, value)
    index.qc = {"peak_gene_associations": len(values)}
    return index


def _write_gene_associations(values: list[intervals.Promoter], bed_path: Path, tsv_path: Path,
                             definition: str, time: str) -> None:
    with bed_path.open("w", encoding="utf-8", newline="") as bed:
        for row in values:
            bed.write(f"{row.chrom}\t{row.start}\t{row.end}\t{row.gene_id}\t0\t{row.strand}\n")
    fields = [
        "promoter_definition", "time", "chrom", "start", "end", "gene_id",
        "gene_symbol", "gene_tss", "gene_strand",
    ]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in values:
            writer.writerow({
                "promoter_definition": definition, "time": time, "chrom": row.chrom,
                "start": row.start, "end": row.end, "gene_id": row.gene_id,
                "gene_symbol": row.gene_symbol, "gene_tss": row.tss,
                "gene_strand": row.strand,
            })


def _write_unique_bed(values: list[intervals.Promoter], path: Path, prefix: str) -> int:
    coordinates = sorted(
        {(row.chrom, row.start, row.end) for row in values},
        key=lambda row: (_chrom_key(row[0]), row[1], row[2]),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        for number, (chrom, start, end) in enumerate(coordinates, 1):
            handle.write(f"{chrom}\t{start}\t{end}\t{prefix}_{number:07d}\n")
    return len(coordinates)


def _merge_true_overlaps(rows: set[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Merge intervals only when they truly overlap; touching intervals remain separate."""
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in rows:
        by_chrom[chrom].append((start, end))
    merged = []
    for chrom in sorted(by_chrom, key=_chrom_key):
        current_start = current_end = None
        for start, end in sorted(by_chrom[chrom]):
            if current_start is None:
                current_start, current_end = start, end
            elif start < current_end:
                current_end = max(current_end, end)
            else:
                merged.append((chrom, current_start, current_end))
                current_start, current_end = start, end
        if current_start is not None:
            merged.append((chrom, current_start, current_end))
    return merged


def _write_consensus(values_by_time: dict[str, list[intervals.Promoter]], path: Path,
                     membership_path: Path, prefix: str) -> int:
    source = {
        (row.chrom, row.start, row.end, time)
        for time, values in values_by_time.items() for row in values
    }
    merged = _merge_true_overlaps({item[:3] for item in source})
    source_index = intervals.IntervalIndex()
    for member in sorted(source, key=lambda row: (_chrom_key(row[0]), row[1], row[2], row[3])):
        source_index.add(member[0], member[1], member[2], member)
    fields = ["consensus_id", "chrom", "start", "end", "member_time", "member_start", "member_end"]
    with path.open("w", encoding="utf-8", newline="") as bed, \
            membership_path.open("w", encoding="utf-8", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for number, (chrom, start, end) in enumerate(merged, 1):
            consensus_id = f"{prefix}_{number:07d}"
            bed.write(f"{chrom}\t{start}\t{end}\t{consensus_id}\n")
            members = sorted(
                source_index.overlaps(chrom, start, end),
                key=lambda row: (row[1], row[2], row[3]),
            )
            for _, member_start, member_end, time in members:
                writer.writerow({
                    "consensus_id": consensus_id, "chrom": chrom, "start": start, "end": end,
                    "member_time": time, "member_start": member_start, "member_end": member_end,
                })
    return len(merged)


def build_promoter_definitions(gtf: str | Path, chrom_sizes: dict[str, int], half_width: int,
                               peak_inputs: dict, destination: Path) -> tuple[dict, list[dict]]:
    """Create the TSS, P1 and P2 promoter indexes and all audit-friendly BED/TSV files."""
    destination.mkdir(parents=True, exist_ok=False)
    tss_index, _ = intervals.load_promoters(gtf, chrom_sizes, half_width)
    tss_values = _promoter_values(tss_index)
    tss_dir = destination / TSS
    tss_dir.mkdir()
    _write_gene_associations(
        tss_values, tss_dir / "promoters_gene_TSS_pm2kb.bed",
        tss_dir / "promoters_gene_TSS_pm2kb.tsv", TSS, "not_time_specific",
    )

    indexes = {TSS: {time: tss_index for time in TIMES}, P1: {}, P2: {}}
    values_by_method = {P1: {}, P2: {}}
    qc = [{
        "promoter_definition": TSS, "time": "not_time_specific",
        "input_H3K27ac_peaks": "NA", "input_ATAC_peaks": "NA", "input_H3K4me3_peaks": "NA",
        "H3K27ac_peaks_overlapping_TSS": "NA",
        "H3K27ac_peak_gene_pairs_overlapping_TSS": "NA",
        "peak_gene_associations": len(tss_values), "unique_output_intervals": len(tss_values),
    }]

    for time in TIMES:
        required = [f"H3K27ac_{time}", f"ATAC_{time}", f"H3K4me3_{time}"]
        missing = [key for key in required if key not in peak_inputs]
        if missing:
            raise ValueError(f"Missing promoter peak inputs: {missing}")
        h3 = load_peak_intervals(peak_inputs[f"H3K27ac_{time}"], set(chrom_sizes))
        atac_rows = load_peak_intervals(peak_inputs[f"ATAC_{time}"], set(chrom_sizes))
        h3k4_rows = load_peak_intervals(peak_inputs[f"H3K4me3_{time}"], set(chrom_sizes))
        atac = make_boolean_index(atac_rows)
        h3k4 = make_boolean_index(h3k4_rows)
        p1_values: set[intervals.Promoter] = set()
        p2_values: set[intervals.Promoter] = set()
        counts = Counter()
        for chrom, start, end in h3:
            linked = list(tss_index.overlaps(chrom, start, end))
            if not linked:
                continue
            counts["H3K27ac_peaks_overlapping_TSS"] += 1
            full_atac = _has_overlap(atac, chrom, start, end)
            full_h3k4 = _has_overlap(h3k4, chrom, start, end)
            for gene in linked:
                counts["H3K27ac_peak_gene_pairs_overlapping_TSS"] += 1
                value = intervals.Promoter(chrom, start, end, gene.gene_id, gene.gene_symbol, gene.tss, gene.strand)
                if full_atac and full_h3k4:
                    p1_values.add(value)
                local_start, local_end = max(start, gene.start), min(end, gene.end)
                if (local_start < local_end
                        and _has_overlap(atac, chrom, local_start, local_end)
                        and _has_overlap(h3k4, chrom, local_start, local_end)):
                    p2_values.add(value)
        if not p2_values.issubset(p1_values):
            raise AssertionError(f"{time}: P2 associations must be a subset of P1 associations")
        for method, values in ((P1, p1_values), (P2, p2_values)):
            method_dir = destination / method
            method_dir.mkdir(exist_ok=True)
            ordered = sorted(values, key=lambda row: (_chrom_key(row.chrom), row.start, row.end, row.gene_id, row.tss))
            indexes[method][time] = _association_index(values)
            values_by_method[method][time] = ordered
            unique_count = _write_unique_bed(ordered, method_dir / f"promoter_{time}.bed", f"{method}_{time}")
            _write_gene_associations(
                ordered, method_dir / f"promoter_{time}_gene_associations.bed",
                method_dir / f"promoter_{time}_gene_associations.tsv", method, time,
            )
            qc.append({
                "promoter_definition": method, "time": time,
                "input_H3K27ac_peaks": len(h3), "input_ATAC_peaks": len(atac_rows),
                "input_H3K4me3_peaks": len(h3k4_rows),
                "H3K27ac_peaks_overlapping_TSS": counts["H3K27ac_peaks_overlapping_TSS"],
                "H3K27ac_peak_gene_pairs_overlapping_TSS": counts["H3K27ac_peak_gene_pairs_overlapping_TSS"],
                "peak_gene_associations": len(ordered), "unique_output_intervals": unique_count,
            })

    for method in (P1, P2):
        method_dir = destination / method
        consensus_count = _write_consensus(
            values_by_method[method], method_dir / "promoter_consensus_0h_6h.bed",
            method_dir / "promoter_consensus_0h_6h_membership.tsv", f"{method}_consensus",
        )
        qc.append({
            "promoter_definition": method, "time": "consensus_0h_6h",
            "input_H3K27ac_peaks": "NA", "input_ATAC_peaks": "NA", "input_H3K4me3_peaks": "NA",
            "H3K27ac_peaks_overlapping_TSS": "NA",
            "H3K27ac_peak_gene_pairs_overlapping_TSS": "NA",
            "peak_gene_associations": "NA", "unique_output_intervals": consensus_count,
        })
    return indexes, qc
