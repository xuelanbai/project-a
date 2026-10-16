#!/usr/bin/env python3
"""Link final AME motif-bearing ATAC cores to genes through HiChIP loops.

This is a downstream draft for the finalized v27 1.5-fold / 500-bp result.
It intentionally does not reuse v24 target counts.  Coordinates are zero-based,
half-open.  Full H3K27ac CREs define enhancer-to-anchor overlap; FIMO sites are
located in the associated 500-bp ATAC-core sequences and retained with explicit
site/core-to-anchor flags for stricter sensitivity summaries.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import subprocess

import enrichment
import intervals
import promoter_definitions
from ame_helpers import write_tsv


MAIN = {str(value) for value in range(1, 19)} | {"X", "Y"}
STATES = ("strength", "weaken")
CONTRASTS = {
    "strength": "strength_1p5_500bp_ATAC_core",
    "weaken": "weaken_1p5_500bp_ATAC_core",
}


def canon(chrom: str) -> str:
    value = str(chrom).strip()
    return value[3:] if value.startswith("chr") else value


def ensure_file(path: str | os.PathLike[str]) -> Path:
    result = Path(path).resolve()
    if not result.is_file() or result.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {result}")
    return result


def read_tsv(path: str | os.PathLike[str]) -> list[dict[str, str]]:
    path = ensure_file(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Invalid TSV header: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No data rows: {path}")
    return rows


def sha256(path: str | os.PathLike[str]) -> str:
    digest = hashlib.sha256()
    with ensure_file(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: str | os.PathLike[str], expected_length: int = 500) -> dict[str, int]:
    path = ensure_file(path)
    lengths: dict[str, int] = {}
    current = None
    sequence_length = 0
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if current is not None:
                    lengths[current] = sequence_length
                current = text[1:].split()[0]
                if not current or current in lengths:
                    raise ValueError(f"{path}:{line_no}: duplicate or empty FASTA identifier")
                sequence_length = 0
            else:
                if current is None:
                    raise ValueError(f"{path}:{line_no}: sequence before FASTA header")
                sequence_length += len(text)
    if current is not None:
        lengths[current] = sequence_length
    if not lengths:
        raise ValueError(f"No FASTA records: {path}")
    invalid = {key: value for key, value in lengths.items() if value != expected_length}
    if invalid:
        raise ValueError(f"Expected {expected_length}-bp FASTA records; examples: {list(invalid.items())[:5]}")
    return lengths


def read_fai_sizes(path: str | os.PathLike[str]) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with ensure_file(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: invalid FASTA index row")
            chrom = canon(fields[0])
            if chrom not in MAIN:
                continue
            if chrom in sizes:
                raise ValueError(f"{path}:{line_no}: duplicate chromosome alias {chrom}")
            sizes[chrom] = int(fields[1])
    if set(sizes) != MAIN:
        raise ValueError(f"Main-chromosome FAI mismatch: {sorted(set(sizes) ^ MAIN)}")
    return sizes


def meme_blocks(path: str | os.PathLike[str]) -> tuple[str, dict[str, str]]:
    path = ensure_file(path)
    header: list[str] = []
    blocks: dict[str, list[str]] = {}
    current = None
    for line in path.read_text(encoding="utf-8").splitlines(True):
        if line.startswith("MOTIF "):
            current = line.split()[1]
            if current in blocks:
                raise ValueError(f"Duplicate motif ID in MEME file: {current}")
            blocks[current] = []
        (header if current is None else blocks[current]).append(line)
    if not blocks:
        raise ValueError(f"No motif blocks in {path}")
    return "".join(header), {key: "".join(value) for key, value in blocks.items()}


def select_motifs(rows: list[dict[str, str]], state: str, cfg: dict) -> list[dict]:
    contrast = CONTRASTS[state]
    q_field = cfg["ame_q_field"]
    selected = []
    for row in rows:
        if row.get("contrast") != contrast:
            continue
        try:
            q_value = float(row[q_field])
            ratio = float(row["enrichment_ratio"])
        except (TypeError, ValueError):
            continue
        if math.isfinite(q_value) and q_value <= float(cfg["ame_q"]) and ratio > 1:
            selected.append({
                "scan_id": row["motif_id"],
                "motif_id": row["motif_id"],
                "motif_name": row["motif_alt_id"],
                "motif_family": row["motif_family"],
                "contrast": contrast,
                "AME_q": q_value,
                "q_field": q_field,
                "enrichment_ratio": ratio,
                "log2_enrichment": math.log2(ratio),
            })
    selected.sort(key=lambda row: (-row["log2_enrichment"], row["AME_q"], row["motif_id"]))
    families = [row["motif_family"] for row in selected]
    if not selected:
        raise ValueError(f"No significant enriched representative AME motifs for {contrast}")
    if len(families) != len(set(families)):
        raise ValueError(f"Representative-family table contains duplicate families for {contrast}")
    return selected


def write_selected_library(database: Path, selected: list[dict], destination: Path) -> None:
    header, blocks = meme_blocks(database)
    missing = [row["motif_id"] for row in selected if row["motif_id"] not in blocks]
    if missing:
        raise ValueError(f"Selected motifs absent from MEME database: {missing}")
    destination.write_text(header + "".join(blocks[row["motif_id"]] for row in selected), encoding="utf-8")


def load_cre_bed(path: Path) -> dict[str, dict]:
    regions: dict[str, dict] = {}
    with ensure_file(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_no}: expected BED4 or wider")
            chrom, start, end, region_id = canon(fields[0]), int(fields[1]), int(fields[2]), fields[3]
            if chrom not in MAIN or not 0 <= start < end or region_id in regions:
                raise ValueError(f"{path}:{line_no}: invalid or duplicate CRE")
            regions[region_id] = {"region_id": region_id, "chrom": chrom, "start": start, "end": end}
    if not regions:
        raise ValueError(f"No CREs in {path}")
    return regions


def load_core_membership(path: Path, fasta_ids: set[str], regions: dict[str, dict]) -> tuple[dict, dict]:
    rows = read_tsv(path)
    required = {
        "core_id", "chrom", "window_start", "window_end", "region_id", "CRE_start", "CRE_end",
        "H3K27ac_log2FC", "ATAC_log2FC", "ATAC_stratum",
    }
    if not required.issubset(rows[0]):
        raise ValueError(f"{path}: missing columns {sorted(required - set(rows[0]))}")
    cores: dict[str, dict] = {}
    used_regions: dict[str, dict] = {}
    for row in rows:
        core_id = row["core_id"]
        if core_id not in fasta_ids:
            continue
        chrom = canon(row["chrom"])
        start, end = int(row["window_start"]), int(row["window_end"])
        region_id = row["region_id"]
        if core_id in cores or end - start != 500 or region_id not in regions:
            raise ValueError(f"Invalid, duplicate or unmapped core: {core_id}")
        region = regions[region_id]
        if (chrom, int(row["CRE_start"]), int(row["CRE_end"])) != (
            region["chrom"], region["start"], region["end"]
        ):
            raise ValueError(f"CRE coordinate mismatch for {core_id}")
        core = {
            "core_id": core_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "region_id": region_id,
            "H3K27ac_log2FC": float(row["H3K27ac_log2FC"]),
            "ATAC_log2FC": float(row["ATAC_log2FC"]),
            "ATAC_stratum": row["ATAC_stratum"],
        }
        cores[core_id] = core
        if region_id not in used_regions:
            used_regions[region_id] = dict(region, H3K27ac_log2FC=core["H3K27ac_log2FC"],
                                           ATAC_log2FC=core["ATAC_log2FC"],
                                           ATAC_stratum=core["ATAC_stratum"], core_ids=[])
        else:
            prior = used_regions[region_id]
            if any(prior[key] != core[key] for key in ("H3K27ac_log2FC", "ATAC_log2FC", "ATAC_stratum")):
                raise ValueError(f"Conflicting signal values across cores in {region_id}")
        used_regions[region_id]["core_ids"].append(core_id)
    missing = fasta_ids - set(cores)
    if missing:
        raise ValueError(f"FASTA cores absent from membership table: {sorted(missing)[:10]}")
    return cores, used_regions


def run_fimo(cfg: dict, library: Path, fasta: Path, destination: Path) -> Path:
    executable = ensure_file(Path(cfg["motif_env"]) / "bin" / "fimo")
    output = destination / "03_FIMO_raw.tsv"
    command = [str(executable), "--text", "--thresh", str(cfg["fimo_p"]),
               "--motif-pseudo", "0", str(library), str(fasta)]
    (destination / "03_FIMO_command.json").write_text(
        json.dumps({"command": command}, indent=2) + "\n", encoding="utf-8"
    )
    environment = dict(os.environ)
    environment["PATH"] = str(Path(cfg["motif_env"]) / "bin") + os.pathsep + environment.get("PATH", "")
    with output.open("w", encoding="utf-8", newline="") as stdout, \
            (destination / "03_FIMO_stderr.log").open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, env=environment, check=False)
    if completed.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"FIMO failed with exit code {completed.returncode}; see {destination / '03_FIMO_stderr.log'}")
    return output


def parse_sites(path: Path, cores: dict, motifs: dict, destination: Path) -> tuple[dict, dict]:
    sites = defaultdict(lambda: defaultdict(list))
    counts = defaultdict(int)
    fields = [
        "region_id", "core_id", "motif_id", "motif_family", "chrom", "start", "end", "strand",
        "FIMO_p", "site_inside_full_CRE",
    ]
    with path.open(encoding="utf-8") as source, gzip.open(
        destination / "04_FIMO_site_coordinates.tsv.gz", "wt", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader((line for line in source if line.strip() and not line.startswith("#")), delimiter="\t")
        required = {"motif_id", "sequence_name", "start", "stop", "strand", "p-value"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Invalid FIMO text header")
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            core_id, motif_id = row["sequence_name"], row["motif_id"]
            if core_id not in cores or motif_id not in motifs:
                raise ValueError(f"Unexpected FIMO sequence or motif: {core_id}, {motif_id}")
            core = cores[core_id]
            relative_start, relative_end = int(row["start"]) - 1, int(row["stop"])
            if not 0 <= relative_start < relative_end <= 500:
                raise ValueError(f"Invalid FIMO site coordinates in {core_id}")
            p_value = float(row["p-value"])
            if not math.isfinite(p_value) or p_value < 0 or p_value > float(motifs[motif_id]["fimo_p"]):
                raise ValueError(f"FIMO returned site beyond configured threshold: {p_value}")
            start, end = core["start"] + relative_start, core["start"] + relative_end
            region = core["region"]
            inside = region["start"] <= start < end <= region["end"]
            record = {
                "region_id": region["region_id"], "core_id": core_id, "motif_id": motif_id,
                "motif_family": motifs[motif_id]["motif_family"], "chrom": core["chrom"],
                "start": start, "end": end, "strand": row["strand"], "FIMO_p": p_value,
                "site_inside_full_CRE": int(inside),
            }
            writer.writerow(record)
            sites[region["region_id"]][motif_id].append(record)
            counts["reported_sites"] += 1
            counts["sites_inside_full_CRE"] += int(inside)
    return sites, dict(counts)


def loop_links(regions: dict, promoters, path: Path, cutoff: float,
               promoter_definition: str, promoter_time: str):
    index = intervals.IntervalIndex()
    for region_id, row in regions.items():
        index.add(row["chrom"], row["start"], row["end"], region_id)
    seen_loops = set()
    with ensure_file(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 9:
                raise ValueError(f"{path}:{line_no}: expected BEDPE with at least 9 columns")
            anchor1 = intervals._interval(fields[0], fields[1], fields[2], f"{path}:{line_no}:anchor1")
            anchor2 = intervals._interval(fields[3], fields[4], fields[5], f"{path}:{line_no}:anchor2")
            contact = intervals._number(fields[6], "contact_count")
            loop_p = intervals._number(fields[7], "loop P", probability=True)
            loop_q = intervals._number(fields[8], "loop Q", probability=True)
            if loop_q > cutoff:
                continue
            loop_id = intervals.stable_loop_id(*anchor1, *anchor2)
            if loop_id in seen_loops:
                continue
            seen_loops.add(loop_id)
            written = set()
            for cre_anchor, promoter_anchor, orientation in ((anchor1, anchor2, "A_to_B"), (anchor2, anchor1, "B_to_A")):
                for region_id in index.overlaps(*cre_anchor):
                    for promoter in promoters.overlaps(*promoter_anchor):
                        key = (region_id, promoter.gene_id, promoter.start, promoter.end)
                        if key in written:
                            continue
                        written.add(key)
                        yield {
                            "promoter_definition": promoter_definition,
                            "promoter_time": promoter_time,
                            "region_id": region_id, "gene_id": promoter.gene_id,
                            "gene_symbol": promoter.gene_symbol, "gene_tss": promoter.tss,
                            "gene_strand": promoter.strand, "promoter_chrom": promoter.chrom,
                            "promoter_start": promoter.start, "promoter_end": promoter.end,
                            "loop_id": loop_id, "loop_p": loop_p, "loop_q": loop_q,
                            "contact_count": contact, "loop_source": str(path), "loop_line": line_no,
                            "CRE_anchor_chrom": cre_anchor[0], "CRE_anchor_start": cre_anchor[1],
                            "CRE_anchor_end": cre_anchor[2], "promoter_anchor_chrom": promoter_anchor[0],
                            "promoter_anchor_start": promoter_anchor[1], "promoter_anchor_end": promoter_anchor[2],
                            "orientation": orientation,
                        }


def connect_state(state: str, regions: dict, cores: dict, sites: dict, motifs: dict,
                  promoters, promoter_definition: str, promoter_time: str,
                  cfg: dict, destination: Path) -> tuple[list[dict], set[str], dict]:
    link_fields = [
        "promoter_definition", "promoter_time",
        "region_id", "gene_id", "gene_symbol", "gene_tss", "gene_strand", "promoter_chrom",
        "promoter_start", "promoter_end", "loop_id", "loop_p", "loop_q", "contact_count",
        "loop_source", "loop_line", "CRE_anchor_chrom", "CRE_anchor_start", "CRE_anchor_end",
        "promoter_anchor_chrom", "promoter_anchor_start", "promoter_anchor_end", "orientation",
    ]
    evidence_fields = [
        "contrast", "state", "loop_time", "CRE_chrom", "CRE_start", "CRE_end", "H3K27ac_log2FC",
        "ATAC_log2FC", "ATAC_stratum", *link_fields, "motif_id", "motif_name", "motif_family",
        "AME_q", "q_field", "enrichment_ratio", "core_ids", "site_count", "sites_inside_full_CRE",
        "best_FIMO_p", "best_site_start", "best_site_end", "best_site_strand",
        "any_core_overlaps_CRE_anchor", "any_site_overlaps_CRE_anchor",
    ]
    groups = {}
    universe: set[str] = set()
    link_count = evidence_count = 0
    loop_path = ensure_file(cfg["loops"][state])
    with gzip.open(destination / "05_all_CRE_loop_gene_links.tsv.gz", "wt", encoding="utf-8", newline="") as link_out, \
            gzip.open(destination / "06_motif_CRE_loop_gene_evidence.tsv.gz", "wt", encoding="utf-8", newline="") as evidence_out:
        link_writer = csv.DictWriter(link_out, fieldnames=link_fields, delimiter="\t", lineterminator="\n")
        evidence_writer = csv.DictWriter(evidence_out, fieldnames=evidence_fields, delimiter="\t", lineterminator="\n")
        link_writer.writeheader()
        evidence_writer.writeheader()
        for link in loop_links(
            regions, promoters, loop_path, float(cfg["loop_q"]), promoter_definition, promoter_time
        ):
            link_writer.writerow(link)
            link_count += 1
            universe.add(link["gene_id"])
            region = regions[link["region_id"]]
            for motif_id, found in sites.get(link["region_id"], {}).items():
                motif = motifs[motif_id]
                best = min(found, key=lambda item: (item["FIMO_p"], item["start"], item["strand"]))
                motif_core_ids = sorted({item["core_id"] for item in found})
                core_anchor = any(
                    cores[core_id]["start"] < link["CRE_anchor_end"]
                    and link["CRE_anchor_start"] < cores[core_id]["end"]
                    for core_id in motif_core_ids
                )
                site_anchor = any(
                    item["start"] < link["CRE_anchor_end"] and link["CRE_anchor_start"] < item["end"]
                    for item in found
                )
                evidence_writer.writerow({
                    "contrast": CONTRASTS[state], "state": state, "loop_time": cfg["loop_times"][state],
                    "CRE_chrom": region["chrom"], "CRE_start": region["start"], "CRE_end": region["end"],
                    "H3K27ac_log2FC": region["H3K27ac_log2FC"], "ATAC_log2FC": region["ATAC_log2FC"],
                    "ATAC_stratum": region["ATAC_stratum"], **link,
                    "motif_id": motif_id, "motif_name": motif["motif_name"],
                    "motif_family": motif["motif_family"], "AME_q": motif["AME_q"],
                    "q_field": motif["q_field"], "enrichment_ratio": motif["enrichment_ratio"],
                    "core_ids": ";".join(motif_core_ids), "site_count": len(found),
                    "sites_inside_full_CRE": sum(item["site_inside_full_CRE"] for item in found),
                    "best_FIMO_p": best["FIMO_p"], "best_site_start": best["start"],
                    "best_site_end": best["end"], "best_site_strand": best["strand"],
                    "any_core_overlaps_CRE_anchor": int(core_anchor),
                    "any_site_overlaps_CRE_anchor": int(site_anchor),
                })
                evidence_count += 1
                key = (motif["motif_family"], link["gene_id"])
                if key not in groups:
                    groups[key] = {
                        "gene_symbol": link["gene_symbol"], "CREs": set(), "loops": set(), "cores": set(),
                        "motifs": set(), "strict_CREs": set(), "best_FIMO_p": 1.0, "best_loop_q": 1.0,
                    }
                group = groups[key]
                group["CREs"].add(link["region_id"])
                group["loops"].add(link["loop_id"])
                group["cores"].update(motif_core_ids)
                group["motifs"].add(motif_id)
                if site_anchor:
                    group["strict_CREs"].add(link["region_id"])
                group["best_FIMO_p"] = min(group["best_FIMO_p"], best["FIMO_p"])
                group["best_loop_q"] = min(group["best_loop_q"], link["loop_q"])
    candidates = []
    for (family, gene_id), group in sorted(groups.items()):
        candidates.append({
            "contrast": CONTRASTS[state], "gene_id": gene_id, "gene_symbol": group["gene_symbol"],
            "group_level": "program", "tf_label": "AME|" + family, "motif_source": "AME",
            "motif_program": family, "unique_CRE_count": len(group["CREs"]),
            "strict_site_anchor_unique_CRE_count": len(group["strict_CREs"]),
            "unique_loop_count": len(group["loops"]), "unique_core_count": len(group["cores"]),
            "motif_count": len(group["motifs"]), "region_ids": ";".join(sorted(group["CREs"])),
            "loop_ids": ";".join(sorted(group["loops"])), "core_ids": ";".join(sorted(group["cores"])),
            "motif_ids": ";".join(sorted(group["motifs"])), "best_FIMO_p": group["best_FIMO_p"],
            "best_loop_q": group["best_loop_q"],
        })
    candidate_fields = list(enrichment.CANDIDATE_COLUMNS) + [
        "motif_source", "motif_program", "unique_CRE_count", "strict_site_anchor_unique_CRE_count",
        "unique_loop_count", "unique_core_count", "motif_count", "region_ids", "loop_ids", "core_ids",
        "motif_ids", "best_FIMO_p", "best_loop_q",
    ]
    write_tsv(destination / "08_candidate_genes.tsv", candidates, candidate_fields)
    write_tsv(destination / "07_callable_genes.tsv", [{"gene_id": gene} for gene in sorted(universe)], ["gene_id"])
    return candidates, universe, {
        "CRE_loop_gene_links": link_count, "motif_evidence_rows": evidence_count,
        "candidate_gene_family_rows": len(candidates), "callable_genes": len(universe),
    }


def bridge_matrices(function_path: Path, destination: Path, cfg: dict) -> list[dict]:
    rows = read_tsv(function_path)
    summary = []
    for state in STATES:
        contrast = CONTRASTS[state]
        state_rows = [row for row in rows if row["contrast"] == contrast and row["group_level"] == "program"]
        for category, predicate in (
            ("lipid", lambda row: row["function_lipid_metabolism"] == "1"),
            ("inflammation", lambda row: row["function_inflammation_immune"] == "1"),
            ("reference_panel", lambda row: row["gene_symbol"] in set(cfg["reference_bridge_genes"])),
        ):
            chosen = [row for row in state_rows if predicate(row)]
            families = sorted({row["motif_program"] for row in chosen})
            genes = {}
            for row in chosen:
                gene_id = row["gene_id"]
                if gene_id not in genes:
                    genes[gene_id] = {"gene_id": gene_id, "gene_symbol": row["gene_symbol"],
                                      **{family: 0 for family in families}}
                genes[gene_id][row["motif_program"]] = int(row["unique_CRE_count"])
            matrix = sorted(genes.values(), key=lambda row: (-sum(int(row[family]) for family in families),
                                                               row["gene_symbol"], row["gene_id"]))
            path = destination / f"{state}_AME_{category}_bridge.tsv"
            write_tsv(path, matrix, ["gene_id", "gene_symbol", *families])
            summary.append({"state": state, "category": category, "genes": len(matrix),
                            "motif_families": len(families), "matrix": str(path)})
    return summary


def write_core_membership(path: Path, cores: dict, regions: dict) -> None:
    rows = []
    for core_id, core in sorted(cores.items(), key=lambda item: (item[1]["chrom"], item[1]["start"], item[0])):
        region = regions[core["region_id"]]
        rows.append({
            "core_id": core_id, "chrom": core["chrom"], "window_start": core["start"],
            "window_end": core["end"], "region_id": core["region_id"], "CRE_start": region["start"],
            "CRE_end": region["end"], "H3K27ac_log2FC": core["H3K27ac_log2FC"],
            "ATAC_log2FC": core["ATAC_log2FC"], "ATAC_stratum": core["ATAC_stratum"],
        })
    write_tsv(path, rows, ["core_id", "chrom", "window_start", "window_end", "region_id", "CRE_start",
                            "CRE_end", "H3K27ac_log2FC", "ATAC_log2FC", "ATAC_stratum"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--promoters-only", action="store_true",
        help="Export and count the three promoter definitions, then stop before motif/HiChIP analysis.",
    )
    args = parser.parse_args()
    cfg = json.loads(ensure_file(args.config).read_text(encoding="utf-8"))
    job_id = os.environ.get("SLURM_JOB_ID", "local")
    destination = Path(cfg["output_root"]).resolve() / f"run_{job_id}"
    destination.mkdir(parents=True, exist_ok=False)

    chrom_sizes = read_fai_sizes(cfg["reference_fasta"] + ".fai")
    promoter_indexes, promoter_qc = promoter_definitions.build_promoter_definitions(
        cfg["gtf"], chrom_sizes, int(cfg["promoter_half_width"]),
        cfg["promoter_peak_inputs"], destination / "promoters",
    )
    write_tsv(
        destination / "02_promoter_definition_counts.tsv", promoter_qc,
        ["promoter_definition", "time", "input_H3K27ac_peaks", "input_ATAC_peaks",
         "input_H3K4me3_peaks", "H3K27ac_peaks_overlapping_TSS",
         "H3K27ac_peak_gene_pairs_overlapping_TSS", "peak_gene_associations",
         "unique_output_intervals"],
    )
    if args.promoters_only:
        promoter_inputs = [
            Path(cfg["reference_fasta"] + ".fai"), Path(cfg["gtf"]),
            *(Path(value) for value in cfg["promoter_peak_inputs"].values()),
        ]
        metadata = {
            "completed": datetime.now(timezone.utc).isoformat(),
            "mode": "promoters_only",
            "promoter_half_width": cfg["promoter_half_width"],
            "promoter_definitions": {
                "TSS_pm2kb": "gene-level TSS +/-2 kb",
                "P1_full_H3_peak_support": "same-time ATAC and H3K4me3 overlap the full TSS-overlapping H3K27ac peak",
                "P2_TSS_local_support_full_H3": "same-time ATAC and H3K4me3 overlap the gene-specific H3K27ac-peak/TSS-window intersection",
            },
            "QC": promoter_qc,
            "input_sha256": {str(path): sha256(path) for path in promoter_inputs},
            "config": cfg,
        }
        (destination / "run_metadata.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        (destination / "completed.json").write_text(
            json.dumps({"completed": True, "mode": "promoters_only", "run": str(destination)},
                       indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({"completed": str(destination), "mode": "promoters_only",
                          "promoter_QC": promoter_qc}, ensure_ascii=False), flush=True)
        return

    source = Path(cfg["source_v27_run"]).resolve()
    if not source.is_dir():
        raise FileNotFoundError(f"Missing v27 run: {source}")
    database = ensure_file(cfg["motif_database"])
    ame_rows = read_tsv(source / "05_AME_representative_families.tsv")

    qc = {}
    selected_summary = []
    state_data = {}
    for state in STATES:
        work = destination / "motifs" / state
        work.mkdir(parents=True)
        selected = select_motifs(ame_rows, state, cfg)
        for index, row in enumerate(selected, 1):
            selected_summary.append(dict(state=state, display_order=index, **row))
        motifs = {row["motif_id"]: dict(row, fimo_p=float(cfg["fimo_p"])) for row in selected}
        library = work / "01_selected_AME_motifs.meme"
        write_selected_library(database, selected, library)
        write_tsv(work / "01_selected_AME_motifs.tsv", selected,
                  ["scan_id", "motif_id", "motif_name", "motif_family", "contrast", "AME_q",
                   "q_field", "enrichment_ratio", "log2_enrichment"])

        state_source = source / "fold_1p5"
        fasta = ensure_file(state_source / "motifs" / f"{state}_500bp" / "foreground.fa")
        fasta_ids = set(read_fasta(fasta))
        cre_bed = state_source / "regions" / f"{state}_CRE.bed"
        membership = state_source / "regions" / f"{state}_ATAC_core_membership.tsv"
        all_regions = load_cre_bed(cre_bed)
        cores, regions = load_core_membership(membership, fasta_ids, all_regions)
        for core in cores.values():
            core["region"] = regions[core["region_id"]]
        write_core_membership(work / "02_final_core_CRE_membership.tsv", cores, regions)
        fimo = run_fimo(cfg, library, fasta, work)
        sites, site_qc = parse_sites(fimo, cores, motifs, work)
        state_data[state] = {
            "regions": regions, "cores": cores, "sites": sites, "motifs": motifs,
        }
        qc[state] = {
            "selected_motif_families": len(selected), "foreground_cores": len(cores),
            "foreground_CREs": len(regions), "FIMO": site_qc,
        }

    write_tsv(destination / "01_selected_AME_motifs.tsv", selected_summary,
              ["state", "display_order", "scan_id", "motif_id", "motif_name", "motif_family",
               "contrast", "AME_q", "q_field", "enrichment_ratio", "log2_enrichment"])
    target_summary = []
    bridge_summary = []
    function_qc = {}
    for definition in (promoter_definitions.TSS, promoter_definitions.P1, promoter_definitions.P2):
        definition_root = destination / "targets" / definition
        definition_root.mkdir(parents=True)
        all_candidates = []
        universes = {}
        qc[definition] = {}
        for state in STATES:
            promoter_time = ("not_time_specific" if definition == promoter_definitions.TSS
                             else cfg["loop_times"][state])
            promoters = promoter_indexes[definition][cfg["loop_times"][state]]
            state_root = definition_root / state
            state_root.mkdir()
            data = state_data[state]
            candidates, universe, link_qc = connect_state(
                state, data["regions"], data["cores"], data["sites"], data["motifs"],
                promoters, definition, promoter_time, cfg, state_root,
            )
            all_candidates.extend(candidates)
            universes[CONTRASTS[state]] = universe
            qc[definition][state] = link_qc
            target_summary.append({
                "promoter_definition": definition, "state": state,
                "promoter_time": promoter_time, **link_qc,
            })
        function_qc[definition] = enrichment.run_enrichment(
            all_candidates, universes, cfg["annotation_tsv"], definition_root / "functions"
        )
        rows = bridge_matrices(
            definition_root / "functions" / "04_TF_candidate_functions.tsv", definition_root, cfg
        )
        for row in rows:
            bridge_summary.append({"promoter_definition": definition, **row})

    write_tsv(destination / "03_target_definition_counts.tsv", target_summary,
              ["promoter_definition", "state", "promoter_time", "CRE_loop_gene_links",
               "motif_evidence_rows", "candidate_gene_family_rows", "callable_genes"])
    write_tsv(destination / "04_bridge_matrix_summary.tsv", bridge_summary,
              ["promoter_definition", "state", "category", "genes", "motif_families", "matrix"])
    qc["function_annotation"] = function_qc
    qc["promoter_definitions"] = promoter_qc

    inputs = [
        source / "05_AME_representative_families.tsv",
        source / "fold_1p5" / "motifs" / "strength_500bp" / "foreground.fa",
        source / "fold_1p5" / "motifs" / "weaken_500bp" / "foreground.fa",
        source / "fold_1p5" / "regions" / "strength_ATAC_core_membership.tsv",
        source / "fold_1p5" / "regions" / "weaken_ATAC_core_membership.tsv",
        Path(cfg["loops"]["strength"]), Path(cfg["loops"]["weaken"]),
        Path(cfg["gtf"]), database,
    ]
    inputs.extend(Path(value) for value in cfg["promoter_peak_inputs"].values())
    metadata = {
        "completed": datetime.now(timezone.utc).isoformat(),
        "source_v27_run": str(source),
        "definition": "final_AME_1p5_500bp_ATAC_core",
        "motif_selection": f"representative family; {cfg['ame_q_field']} <= {cfg['ame_q']}; enrichment_ratio > 1",
        "FIMO_threshold": f"P <= {cfg['fimo_p']}",
        "primary_loop_link": "full H3K27ac CRE overlaps one anchor >=1 bp; opposite anchor overlaps one of three explicit promoter definitions >=1 bp",
        "promoter_definitions": {
            "TSS_pm2kb": "gene-level TSS +/-2 kb",
            "P1_full_H3_peak_support": "full H3K27ac peak overlaps TSS +/-2 kb and has same-time ATAC plus H3K4me3 support anywhere in the full peak; full H3K27ac coordinate retained",
            "P2_TSS_local_support_full_H3": "same-time ATAC plus H3K4me3 support inside the gene-specific intersection of H3K27ac peak and TSS +/-2 kb; full H3K27ac coordinate retained",
        },
        "strict_sensitivity_flag": "motif site itself overlaps the CRE-side anchor >=1 bp",
        "loop_time_policy": cfg["loop_times"],
        "loop_q": cfg["loop_q"],
        "QC": qc,
        "input_sha256": {str(path): sha256(path) for path in inputs},
        "config": cfg,
    }
    (destination / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (destination / "completed.json").write_text(
        json.dumps({"completed": True, "run": str(destination), "QC": qc}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"completed": str(destination), "QC": qc}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
