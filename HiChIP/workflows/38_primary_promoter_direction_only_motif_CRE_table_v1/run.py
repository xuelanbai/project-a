#!/usr/bin/env python3
"""Primary-promoter log2FC reclassification, HiChIP remapping and tables.

This is a fast reuse analysis.  It never changes the final enhancer coordinates
and it never recounts BAMs.  Instead it reuses the archived all-region mark
matrices, filters the final P1 promoter consensus to promoter--gene pairs for
which the complete consensus interval is contained by that gene's TSS +/-2 kb,
selects at most one primary promoter per gene, and reclassifies enhancer and
promoter H3K27ac changes using absolute log2FC cutoffs.

Both 0 h and 6 h FitHiChIP BEDPE files are remapped to all enhancer states and
the primary promoter--gene associations.  The companion motif task scans
500-bp ATAC-consensus-core windows inside the retained primary promoters.
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
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

import intervals


MAIN = {str(x) for x in range(1, 19)} | {"X", "Y"}
ASSAYS = ("H3K27ac", "H3K4me3", "Pol2", "ATAC", "H3K27me3")
STATE_ORDER = ("strength", "stable", "weaken")


def canon(chrom: str) -> str:
    value = str(chrom).strip()
    return value[3:] if value.startswith("chr") else value


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict]:
    with open_text(path, "r") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields and not key.startswith("_"):
                    fields.append(key)
    with open_text(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_bed4(path: Path, rows) -> None:
    """Write headerless BED4 with canonical pig chromosome names."""
    path.parent.mkdir(parents=True, exist_ok=True)
    chrom_order = {str(value): value for value in range(1, 19)} | {"X": 19, "Y": 20}
    ordered = sorted(rows, key=lambda row: (
        chrom_order[canon(row["chrom"])], int(row["start"]), int(row["end"]), row["region_id"]
    ))
    with path.open("w", encoding="utf-8", newline="") as handle:
        for row in ordered:
            handle.write(
                f"{canon(row['chrom'])}\t{int(row['start'])}\t{int(row['end'])}\t{row['region_id']}\n"
            )


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


def threshold_token(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def state_and_reason(row: dict, linear_fold_threshold: float, minimum: int) -> tuple[str, str]:
    fold = float(linear_fold_threshold)
    if fold < 1:
        raise ValueError("A symmetric linear fold threshold must be >=1")
    cutoff = math.log2(fold)
    value = as_float(row, "H3K27ac_log2FC")
    strength_like = value > cutoff if math.isclose(fold, 1.0) else value >= cutoff
    weaken_like = value < -cutoff if math.isclose(fold, 1.0) else value <= -cutoff
    if strength_like:
        if as_int(row, "supported_6h") != 1:
            return "stable", "strength_like_but_no_6h_peak_support"
        if min(as_float(row, "H3K27ac_6h_raw_tags"), as_float(row, "ATAC_6h_raw_tags")) < minimum:
            return "stable", "strength_like_but_active_state_raw_tags_below_minimum"
        return "strength", f"H3K27ac_linear_FC_above_{fold:g}"
    if weaken_like:
        if as_int(row, "supported_0h") != 1:
            return "stable", "weaken_like_but_no_0h_peak_support"
        if min(as_float(row, "H3K27ac_0h_raw_tags"), as_float(row, "ATAC_0h_raw_tags")) < minimum:
            return "stable", "weaken_like_but_active_state_raw_tags_below_minimum"
        return "weaken", f"H3K27ac_linear_FC_below_{1.0 / fold:g}"
    return "stable", f"H3K27ac_not_directional_at_linear_fold_{fold:g}"


def validate_mark_rows(rows: list[dict], expected_type: str) -> dict[str, dict]:
    required = {
        "region_type", "region_id", "chrom", "start", "end", "source_times",
        "supported_0h", "supported_6h", "H3K27ac_0h_raw_tags", "H3K27ac_6h_raw_tags",
        "ATAC_0h_raw_tags", "ATAC_6h_raw_tags", "H3K27ac_log2FC", "ATAC_log2FC",
    }
    if not required.issubset(rows[0]):
        raise ValueError(f"{expected_type} mark matrix missing columns: {sorted(required - set(rows[0]))}")
    result = {}
    for raw in rows:
        if raw["region_type"] != expected_type:
            raise ValueError(f"Unexpected region_type={raw['region_type']!r}")
        region_id = raw["region_id"]
        if region_id in result:
            raise ValueError(f"Duplicate region_id: {region_id}")
        chrom, start, end = canon(raw["chrom"]), as_int(raw, "start"), as_int(raw, "end")
        if chrom not in MAIN or not 0 <= start < end:
            raise ValueError(f"Invalid {expected_type} interval {chrom}:{start}-{end}")
        result[region_id] = dict(raw, chrom=chrom, start=start, end=end, length_bp=end - start)
    return result


def classify_rows(rows: dict[str, dict], thresholds: list[float], primary: float,
                  minimum: int) -> None:
    for row in rows.values():
        for threshold in thresholds:
            state, reason = state_and_reason(row, threshold, minimum)
            token = threshold_token(threshold)
            row[f"state_fold_{token}"] = state
            row[f"classification_reason_fold_{token}"] = reason
        primary_token = threshold_token(primary)
        row["state"] = row[f"state_fold_{primary_token}"]
        if row.get("region_type") == "promoter":
            row["promoter_state"] = row["state"]
        row["state_linear_fold_threshold"] = primary
        row["classification_reason"] = row[f"classification_reason_fold_{primary_token}"]
        for assay in ASSAYS:
            if f"{assay}_log2FC" in row:
                row[f"{assay}_fold_6h_vs_0h"] = 2.0 ** as_float(row, f"{assay}_log2FC")


def strict_promoter_associations(path: Path, promoters: dict[str, dict], half_width: int):
    rows = read_tsv(path)
    required = {
        "region_id", "gene_id", "gene_symbol", "gene_tss", "gene_strand", "chrom",
        "promoter_start", "promoter_end",
    }
    if not required.issubset(rows[0]):
        raise ValueError(f"Promoter association table missing columns: {sorted(required - set(rows[0]))}")
    strict, seen = [], set()
    qc = Counter(input_associations=len(rows))
    for raw in rows:
        region_id = raw["region_id"]
        promoter = promoters.get(region_id)
        if promoter is None:
            raise ValueError(f"Association references unknown promoter: {region_id}")
        tss = as_int(raw, "gene_tss")
        window_start, window_end = max(0, tss - half_width), tss + half_width
        qc["input_unique_promoters"] += 0
        if promoter["start"] < window_start or promoter["end"] > window_end:
            qc["excluded_not_fully_contained_in_gene_TSS_window"] += 1
            continue
        key = (region_id, raw["gene_id"], canon(raw["chrom"]), tss)
        if key in seen:
            qc["duplicate_identical_associations"] += 1
            continue
        seen.add(key)
        strict.append({
            "region_id": region_id, "gene_id": raw["gene_id"],
            "gene_symbol": raw.get("gene_symbol") or raw["gene_id"], "gene_tss": tss,
            "gene_strand": raw["gene_strand"], "chrom": canon(raw["chrom"]),
            "promoter_start": promoter["start"], "promoter_end": promoter["end"],
            "TSS_window_start": window_start, "TSS_window_end": window_end,
            "complete_promoter_inside_TSS_window": 1,
        })
    retained = {row["region_id"] for row in strict}
    qc["input_unique_promoters"] = len({row["region_id"] for row in rows})
    qc["strict_associations"] = len(strict)
    qc["strict_unique_promoters"] = len(retained)
    qc["strict_unique_genes"] = len({row["gene_id"] for row in strict})
    return strict, retained, dict(qc)


def primary_promoter_rank(association: dict, promoter: dict) -> tuple:
    """Gene-centric primary-promoter rank; signal density is only a tie-breaker."""
    tss, start, end = int(association["gene_tss"]), promoter["start"], promoter["end"]
    contains_tss = start <= tss < end
    if contains_tss:
        distance = 0
    elif tss < start:
        distance = start - tss
    else:
        distance = tss - (end - 1)
    length_kb = (end - start) / 1000.0
    h3k4_density = (
        as_float(promoter, "H3K4me3_0h_tags_per_10M")
        + as_float(promoter, "H3K4me3_6h_tags_per_10M")
    ) / 2.0 / length_kb
    atac_density = (
        as_float(promoter, "ATAC_0h_tags_per_10M")
        + as_float(promoter, "ATAC_6h_tags_per_10M")
    ) / 2.0 / length_kb
    association["candidate_contains_TSS"] = int(contains_tss)
    association["candidate_distance_to_TSS_bp"] = distance
    association["candidate_mean_H3K4me3_per10M_per_kb"] = h3k4_density
    association["candidate_mean_ATAC_per10M_per_kb"] = atac_density
    return (not contains_tss, distance, -h3k4_density, -atac_density, association["region_id"])


def select_primary_promoters(associations: list[dict], promoters: dict[str, dict]):
    """Retain at most one promoter association for each gene."""
    grouped = defaultdict(list)
    for association in associations:
        grouped[association["gene_id"]].append(association)
    selected = []
    qc = Counter(strict_genes=len(grouped), strict_candidate_associations=len(associations))
    for gene_id, candidates in sorted(grouped.items()):
        loci = {(row["chrom"], row["gene_tss"], row["gene_strand"]) for row in candidates}
        if len(loci) != 1:
            raise ValueError(f"Gene has inconsistent loci/TSS values: {gene_id}: {sorted(loci)}")
        ranked = sorted(candidates, key=lambda row: primary_promoter_rank(row, promoters[row["region_id"]]))
        chosen_id = ranked[0]["region_id"]
        if len(ranked) > 1:
            qc["genes_with_multiple_strict_candidate_promoters"] += 1
        qc[f"genes_with_{len(ranked)}_strict_candidates"] += 1
        for rank, row in enumerate(ranked, 1):
            row["candidate_promoters_for_gene"] = len(ranked)
            row["primary_promoter_rank_for_gene"] = rank
            row["is_primary_promoter_for_gene"] = int(rank == 1)
            row["selected_primary_promoter_id"] = chosen_id
            row["primary_selection_rule"] = (
                "contains_TSS_then_nearest_interval_then_mean_H3K4me3_density_"
                "then_mean_ATAC_density_then_region_id"
            )
            row["primary_selection_reason"] = (
                "selected_contains_TSS" if rank == 1 and row["candidate_contains_TSS"]
                else "selected_nearest_supported_peak" if rank == 1
                else "not_selected_higher_ranked_candidate"
            )
        selected.append(ranked[0])
    promoter_counts = Counter(row["region_id"] for row in selected)
    qc["primary_promoter_gene_associations"] = len(selected)
    qc["unique_primary_promoter_regions"] = len(promoter_counts)
    qc["primary_promoter_regions_shared_by_multiple_genes"] = sum(
        count > 1 for count in promoter_counts.values()
    )
    qc["maximum_genes_per_primary_promoter_region"] = max(promoter_counts.values(), default=0)
    return selected, associations, dict(qc)


def add_promoter_gene_columns(promoters: dict[str, dict], associations: list[dict]) -> None:
    grouped = defaultdict(list)
    for row in associations:
        grouped[row["region_id"]].append(row)
    for region_id, promoter in promoters.items():
        rows = grouped[region_id]
        promoter["primary_gene_count"] = len({row["gene_id"] for row in rows})
        promoter["primary_gene_ids"] = ";".join(sorted({row["gene_id"] for row in rows}))
        promoter["primary_gene_symbols"] = ";".join(sorted({row["gene_symbol"] for row in rows}))
        promoter["primary_gene_TSS"] = ";".join(str(x) for x in sorted({row["gene_tss"] for row in rows}))


def interval_index(rows, value_builder):
    index = intervals.IntervalIndex()
    for row in rows:
        index.add(canon(row["chrom"]), int(row["start"]), int(row["end"]), value_builder(row))
    return index


def read_loops(path: Path, q_max: float):
    seen = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 9:
                raise ValueError(f"{path}:{line_no}: expected at least 9 BEDPE columns")
            a1 = intervals._interval(fields[0], fields[1], fields[2], f"{path}:{line_no}:anchor1")
            a2 = intervals._interval(fields[3], fields[4], fields[5], f"{path}:{line_no}:anchor2")
            contact = intervals._number(fields[6], "contact_count")
            loop_p = intervals._number(fields[7], "loop_P", probability=True)
            loop_q = intervals._number(fields[8], "loop_Q", probability=True)
            if loop_q > q_max:
                continue
            loop_id = intervals.stable_loop_id(*a1, *a2)
            if loop_id in seen:
                continue
            seen.add(loop_id)
            yield {
                "loop_id": loop_id, "anchor1": a1, "anchor2": a2,
                "contact_count": contact, "loop_p": loop_p, "loop_q": loop_q,
                "loop_source": str(path), "loop_line": line_no,
            }


def link_policy(enhancer_state: str, loop_time: str) -> str:
    if enhancer_state == "strength" and loop_time == "6h":
        return "time_matched_active_state"
    if enhancer_state == "weaken" and loop_time == "0h":
        return "time_matched_active_state"
    if enhancer_state == "stable":
        return "stable_region_both_timepoints_retained"
    return "opposite_timepoint_sensitivity"


def map_loops(loop_paths: dict[str, Path], enhancers: dict[str, dict],
              promoters: dict[str, dict], associations: list[dict], q_max: float):
    enhancer_index = interval_index(enhancers.values(), lambda row: row["region_id"])
    promoter_index = intervals.IntervalIndex()
    for association in associations:
        promoter_index.add(
            association["chrom"], association["promoter_start"], association["promoter_end"], association
        )
    result, seen = [], set()
    qc = Counter()
    for loop_time, path in sorted(loop_paths.items()):
        for loop in read_loops(path, q_max):
            qc[f"{loop_time}_loops_passing_q"] += 1
            for cre_anchor, promoter_anchor, orientation in (
                (loop["anchor1"], loop["anchor2"], "A_to_B"),
                (loop["anchor2"], loop["anchor1"], "B_to_A"),
            ):
                enhancer_ids = list(enhancer_index.overlaps(*cre_anchor))
                promoter_hits = list(promoter_index.overlaps(*promoter_anchor))
                for enhancer_id in enhancer_ids:
                    enhancer = enhancers[enhancer_id]
                    for association in promoter_hits:
                        promoter_id = association["region_id"]
                        key = (loop_time, loop["loop_id"], enhancer_id, promoter_id,
                               association["gene_id"], orientation)
                        if key in seen:
                            continue
                        seen.add(key)
                        promoter = promoters[promoter_id]
                        row = {
                            "loop_time": loop_time, "loop_policy": link_policy(enhancer["state"], loop_time),
                            "loop_id": loop["loop_id"], "loop_p": loop["loop_p"], "loop_q": loop["loop_q"],
                            "contact_count": loop["contact_count"], "loop_source": loop["loop_source"],
                            "loop_line": loop["loop_line"], "orientation": orientation,
                            "CRE_anchor_chrom": cre_anchor[0], "CRE_anchor_start": cre_anchor[1],
                            "CRE_anchor_end": cre_anchor[2],
                            "promoter_anchor_chrom": promoter_anchor[0],
                            "promoter_anchor_start": promoter_anchor[1],
                            "promoter_anchor_end": promoter_anchor[2],
                            "enhancer_id": enhancer_id, "enhancer_state": enhancer["state"],
                            "enhancer_chrom": enhancer["chrom"], "enhancer_start": enhancer["start"],
                            "enhancer_end": enhancer["end"],
                            "promoter_id": promoter_id, "promoter_state": promoter["state"],
                            "promoter_chrom": promoter["chrom"], "promoter_start": promoter["start"],
                            "promoter_end": promoter["end"],
                            "gene_id": association["gene_id"], "gene_symbol": association["gene_symbol"],
                            "gene_tss": association["gene_tss"], "gene_strand": association["gene_strand"],
                        }
                        for prefix, feature in (("enhancer", enhancer), ("promoter", promoter)):
                            for assay in ASSAYS:
                                for suffix in ("0h_tags_per_10M", "6h_tags_per_10M", "log2FC",
                                               "fold_6h_vs_0h"):
                                    field = f"{assay}_{suffix}"
                                    row[f"{prefix}_{field}"] = feature.get(field, "NA")
                        result.append(row)
    qc["link_rows"] = len(result)
    qc["unique_loops_linked"] = len({row["loop_id"] for row in result})
    qc["unique_enhancers_linked"] = len({row["enhancer_id"] for row in result})
    qc["unique_promoters_linked"] = len({row["promoter_id"] for row in result})
    qc["unique_genes_linked"] = len({row["gene_id"] for row in result})
    return result, dict(qc)


def add_link_summaries(enhancers: dict[str, dict], promoters: dict[str, dict], links: list[dict]) -> None:
    by_enhancer, by_promoter = defaultdict(list), defaultdict(list)
    for row in links:
        by_enhancer[row["enhancer_id"]].append(row)
        by_promoter[row["promoter_id"]].append(row)
    for region_id, row in enhancers.items():
        linked = by_enhancer.get(region_id, [])
        row["linked_promoter_count"] = len({x["promoter_id"] for x in linked})
        row["linked_gene_count"] = len({x["gene_id"] for x in linked})
        row["linked_gene_ids"] = ";".join(sorted({x["gene_id"] for x in linked}))
        row["linked_gene_symbols"] = ";".join(sorted({x["gene_symbol"] for x in linked}))
        row["linked_loop_count"] = len({x["loop_id"] for x in linked})
        row["linked_loop_times"] = ";".join(sorted({x["loop_time"] for x in linked}))
    for region_id, row in promoters.items():
        linked = by_promoter.get(region_id, [])
        row["linked_enhancer_count"] = len({x["enhancer_id"] for x in linked})
        row["linked_loop_count"] = len({x["loop_id"] for x in linked})
        row["linked_loop_times"] = ";".join(sorted({x["loop_time"] for x in linked}))


def promoter_gene_table(promoters: dict[str, dict], associations: list[dict]) -> list[dict]:
    rows = []
    for association in associations:
        promoter = promoters[association["region_id"]]
        row = dict(association)
        row.update({key: value for key, value in promoter.items() if not key.startswith("_")})
        rows.append(row)
    return rows


def state_count_rows(enhancers: dict[str, dict], promoters: dict[str, dict],
                     associations: list[dict], thresholds: list[float]):
    rows = []
    for threshold in thresholds:
        field = f"state_fold_{threshold_token(threshold)}"
        for region_type, data in (("enhancer", enhancers), ("promoter", promoters)):
            for state in STATE_ORDER:
                rows.append({
                    "H3K27ac_linear_fold_threshold": threshold,
                    "region_type": region_type, "state": state,
                    "unique_regions": sum(row[field] == state for row in data.values()),
                    "unique_genes": (
                        len({a["gene_id"] for a in associations if data[a["region_id"]][field] == state})
                        if region_type == "promoter" else "NA"
                    ),
                })
    return rows


def require_file(path: Path) -> Path:
    if not path.is_absolute() or not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Submit this analysis through Slurm; do not run it on the management node.")
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    primary_threshold = float(cfg["primary_linear_fold_threshold"])
    thresholds = sorted({float(value) for value in cfg["sensitivity_linear_fold_thresholds"]})
    if thresholds != [1.0, 1.2] or not math.isclose(primary_threshold, 1.0):
        raise ValueError("This version fixes primary directional FC>1/<1 and sensitivity at 1.2-fold")
    if int(cfg["promoter_TSS_half_width"]) != 2000:
        raise ValueError("This version fixes the promoter window at gene TSS +/-2 kb")
    output = Path(cfg["results"]) / f"run_{os.environ['SLURM_JOB_ID']}"
    output.mkdir(parents=True, exist_ok=False)
    (output / "code").mkdir()
    for name in ("run.py", "motif.py", "intervals.py", "config.json", "README.md", "run.slurm",
                 "submit.sh", "check_env.sh", "test_small.py"):
        source = Path(__file__).with_name(name)
        if source.is_file():
            shutil.copy2(source, output / "code" / name)
    print(f"results={output}", flush=True)

    p1_run = Path(cfg["source_P1_promoter_state_run"])
    mark_run = Path(cfg["source_mark_matrix_run"])
    p1_matrix_path = require_file(p1_run / "04_P1_promoter_signal_and_state_matrix.tsv")
    association_path = require_file(p1_run / "05_P1_promoter_state_gene_associations.tsv")
    enhancer_mark_path = require_file(mark_run / "01_enhancer_mark_signal_matrix.tsv")
    promoter_mark_path = require_file(mark_run / "02_P1_promoter_mark_signal_matrix.tsv")
    core_path = require_file(Path(cfg["source_ATAC_core_run"]) / "02_ATAC_consensus_cores.tsv")
    loop_paths = {time: require_file(Path(path)) for time, path in cfg["loops"].items()}

    source_promoters = validate_mark_rows(read_tsv(promoter_mark_path), "promoter")
    source_p1 = {row["region_id"]: row for row in read_tsv(p1_matrix_path)}
    if set(source_promoters) != set(source_p1):
        raise ValueError("Archived promoter state and mark matrices have different region IDs")
    strict_associations, strict_ids, definition_qc = strict_promoter_associations(
        association_path, source_promoters, int(cfg["promoter_TSS_half_width"])
    )
    strict_candidates = {region_id: source_promoters[region_id] for region_id in sorted(strict_ids)}
    classify_rows(strict_candidates, thresholds, primary_threshold,
                  int(cfg["min_raw_tags_active_state"]))
    associations, candidate_associations, primary_qc = select_primary_promoters(
        strict_associations, strict_candidates
    )
    primary_ids = {row["region_id"] for row in associations}
    promoters = {region_id: strict_candidates[region_id] for region_id in sorted(primary_ids)}
    enhancers = validate_mark_rows(read_tsv(enhancer_mark_path), "enhancer")
    classify_rows(enhancers, thresholds, primary_threshold,
                  int(cfg["min_raw_tags_active_state"]))
    add_promoter_gene_columns(promoters, associations)
    for association in candidate_associations:
        promoter = strict_candidates[association["region_id"]]
        association["promoter_state"] = promoter["state"]
        association["classification_reason"] = promoter["classification_reason"]
        for threshold in thresholds:
            token = threshold_token(threshold)
            association[f"promoter_state_fold_{token}"] = promoter[f"state_fold_{token}"]
    definition_qc.update(primary_qc)

    links, link_qc = map_loops(loop_paths, enhancers, promoters, associations, float(cfg["loop_q_max"]))
    add_link_summaries(enhancers, promoters, links)

    write_tsv(output / "01_strict_promoter_definition_QC.tsv", [definition_qc])
    write_tsv(output / "01_state_counts.tsv",
              state_count_rows(enhancers, promoters, associations, thresholds))
    write_tsv(output / "02_strict_promoter_candidate_gene_ranking.tsv", candidate_associations)
    primary_matrix = output / "02_primary_promoter_region_state_mark.tsv"
    write_tsv(primary_matrix, list(promoters.values()))
    write_tsv(output / "02_primary_promoter_gene_state_mark.tsv",
              promoter_gene_table(promoters, associations))
    write_tsv(output / "03_enhancer_region_state_mark.tsv", list(enhancers.values()))
    write_tsv(output / "04_enhancer_promoter_gene_loop_mark.tsv.gz", links)
    write_tsv(output / "04_HiChIP_mapping_QC.tsv", [link_qc])
    bed_dir = output / "BED"
    write_bed4(bed_dir / "primary_promoter_all.bed", promoters.values())
    write_bed4(bed_dir / "enhancer_all.bed", enhancers.values())
    for state in STATE_ORDER:
        write_bed4(bed_dir / f"primary_promoter_{state}_direction_only.bed",
                   (row for row in promoters.values() if row["state"] == state))
        write_bed4(bed_dir / f"enhancer_{state}_direction_only.bed",
                   (row for row in enhancers.values() if row["state"] == state))
    (bed_dir / "README.txt").write_text(
        "All BED files are headerless BED4.\n"
        "Columns: 1 chromosome (Sscrofa11.1 canonical name), 2 zero-based start, "
        "3 half-open end, 4 stable region_id.\n"
        "Promoter coordinates are complete P1 consensus intervals selected as a primary promoter "
        "for at least one gene; enhancer coordinates are unchanged final-v27 CREs.\n",
        encoding="utf-8",
    )

    motif_output = output / "05_primary_promoter_motif_direction_only"
    if cfg.get("run_motif", True):
        command = [
            sys.executable, str(Path(__file__).with_name("motif.py")), "--config", str(args.config),
            "--promoter-matrix", str(primary_matrix), "--output", str(motif_output),
        ]
        print("step=promoter_motif", flush=True)
        subprocess.run(command, check=True)

    metadata = {
        "completed": datetime.now(timezone.utc).isoformat(),
        "analysis": "one_primary_promoter_per_gene_TSS_pm2kb_direction_only_fast_reuse",
        "primary_promoter_definition": (
            "complete archived P1 consensus interval fully contained in its associated gene TSS +/-2 kb; "
            "same-time ATAC and H3K4me3 support inherited from the validated P1 source; at most one "
            "promoter selected per gene by TSS containment/proximity with signal-density tie-breakers"
        ),
        "enhancer_definition": "unchanged final version-27 H3K27ac consensus CRE",
        "state_definition": "primary direction-only split: linear FC>1 Strength and FC<1 Weaken; 1.2-fold retained as sensitivity columns; active-state support and raw H3K27ac/ATAC tags >=5",
        "signal_source": "archived all-region BAM recount matrices; no BAM is recounted in this fast reuse run",
        "loop_policy": "all enhancer states mapped to both 0h and 6h BEDPE; loop_time and time-matched flag retained",
        "counts": {
            "enhancers": len(enhancers), "primary_promoters": len(promoters),
            "primary_promoter_gene_associations": len(associations), "links": len(links),
            "enhancer_states": dict(Counter(x["state"] for x in enhancers.values())),
            "promoter_states": dict(Counter(x["state"] for x in promoters.values())),
        },
        "definition_QC": definition_qc, "HiChIP_QC": link_qc,
        "input_sha256": {str(path): sha256(path) for path in (
            p1_matrix_path, association_path, enhancer_mark_path, promoter_mark_path, core_path, *loop_paths.values()
        )},
        "config": cfg,
    }
    (output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (output / "completed.json").write_text(
        json.dumps({"completed": True, "results": str(output)}, indent=2), encoding="utf-8"
    )
    print(f"completed=1\nresults={output}", flush=True)


if __name__ == "__main__":
    main()
