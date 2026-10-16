#!/usr/bin/env python3
"""Build a self-contained plotting bundle for unbiased P1 enhancer hubs.

The representative selection never uses gene names or pathway labels:

1. Rank AME families by AME log2 enrichment (top five per state).
2. Keep only evidence where the FIMO motif site overlaps the CRE-side loop
   anchor.
3. For each family, rank enhancers by linked-gene count, loop count,
   absolute H3K27ac log2FC, loop q, FIMO P and stable region ID.
4. Merge families that select the same enhancer into one hub.
5. Test every eligible GO BP and Reactome term for each unique hub target set
   against the state-specific P1 callable-gene universe already recorded by
   run_42014.  BH correction is reported within hub/source and globally across
   all hub/source/term tests in this bundle.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path


STATES = ("strength", "weaken")
SOURCES = ("GO_BP", "Reactome")


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, population: int, successes: int, draws: int) -> float:
    values = [
        log_comb(successes, x)
        + log_comb(population - successes, draws - x)
        - log_comb(population, draws)
        for x in range(k, min(successes, draws) + 1)
    ]
    maximum = max(values)
    return min(1.0, math.exp(maximum) * sum(math.exp(value - maximum) for value in values))


def bh(rows: list[dict], p_field: str, q_field: str) -> None:
    order = sorted(range(len(rows)), key=lambda index: rows[index][p_field])
    ceiling = 1.0
    for rank, index in reversed(list(enumerate(order, 1))):
        ceiling = min(ceiling, rows[index][p_field] * len(rows) / rank)
        rows[index][q_field] = ceiling


def evidence_path(run: Path, state: str) -> Path:
    return (
        run
        / "targets"
        / "P1_full_H3_peak_support"
        / state
        / "06_motif_CRE_loop_gene_evidence.tsv.gz"
    )


def load_selected_families(run: Path, top_n: int):
    motifs = list(read_tsv(run / "01_selected_AME_motifs.tsv"))
    selected = {}
    motif_meta = {}
    for state in STATES:
        state_rows = [row for row in motifs if row["contrast"].startswith(state)]
        state_rows.sort(key=lambda row: (-float(row["log2_enrichment"]), row["motif_family"]))
        selected[state] = [row["motif_family"] for row in state_rows[:top_n]]
        for rank, row in enumerate(state_rows, 1):
            motif_meta[state, row["motif_family"]] = {
                "AME_rank": rank,
                "motif_id": row["motif_id"],
                "motif_name": row["motif_name"],
                "AME_q": float(row["AME_q"]),
                "enrichment_ratio": float(row["enrichment_ratio"]),
                "log2_enrichment": float(row["log2_enrichment"]),
            }
    return selected, motif_meta


def collect_strict_evidence(run: Path, selected: dict):
    groups = {}
    strict_rows = defaultdict(list)
    for state in STATES:
        with gzip.open(evidence_path(run, state), "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                family = row["motif_family"]
                if family not in selected[state] or row["any_site_overlaps_CRE_anchor"] != "1":
                    continue
                key = (state, family, row["region_id"])
                group = groups.setdefault(
                    key,
                    {
                        "state": state,
                        "family": family,
                        "region_id": row["region_id"],
                        "chrom": row["CRE_chrom"],
                        "start": int(row["CRE_start"]),
                        "end": int(row["CRE_end"]),
                        "H3K27ac_log2FC": float(row["H3K27ac_log2FC"]),
                        "ATAC_log2FC": float(row["ATAC_log2FC"]),
                        "ATAC_stratum": row["ATAC_stratum"],
                        "genes": {},
                        "loops": set(),
                        "best_loop_q": 1.0,
                        "best_FIMO_p": 1.0,
                    },
                )
                group["genes"][row["gene_id"]] = row["gene_symbol"]
                group["loops"].add(row["loop_id"])
                group["best_loop_q"] = min(group["best_loop_q"], float(row["loop_q"]))
                group["best_FIMO_p"] = min(group["best_FIMO_p"], float(row["best_FIMO_p"]))
                strict_rows[key].append(row)
    return groups, strict_rows


def choose_hubs(selected: dict, groups: dict):
    chosen = []
    for state in STATES:
        for family in selected[state]:
            candidates = [
                value for (candidate_state, candidate_family, _), value in groups.items()
                if candidate_state == state and candidate_family == family
            ]
            candidates.sort(
                key=lambda row: (
                    -len(row["genes"]),
                    -len(row["loops"]),
                    -abs(row["H3K27ac_log2FC"]),
                    row["best_loop_q"],
                    row["best_FIMO_p"],
                    row["region_id"],
                )
            )
            if not candidates:
                raise ValueError(f"No strict evidence for {state} {family}")
            chosen.append(candidates[0])

    hubs = []
    for state in STATES:
        by_region = {}
        for row in [item for item in chosen if item["state"] == state]:
            hub = by_region.setdefault(
                row["region_id"],
                {
                    "state": state,
                    "region_id": row["region_id"],
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "H3K27ac_log2FC": row["H3K27ac_log2FC"],
                    "ATAC_log2FC": row["ATAC_log2FC"],
                    "ATAC_stratum": row["ATAC_stratum"],
                    "families": [],
                },
            )
            hub["families"].append(row["family"])
        for number, hub in enumerate(by_region.values(), 1):
            hub["hub_id"] = f"{state}_HUB_{number:02d}"
            hubs.append(hub)
    return hubs


def build_hub_tables(hubs: list[dict], strict_rows: dict, motif_meta: dict):
    hub_rows, target_rows, loop_rows, family_rows = [], [], [], []
    hub_gene_sets = {}
    for hub in hubs:
        state, region_id = hub["state"], hub["region_id"]
        evidence = []
        for family in hub["families"]:
            evidence.extend(strict_rows[state, family, region_id])
            meta = motif_meta[state, family]
            family_rows.append({
                "hub_id": hub["hub_id"], "state": state, "region_id": region_id,
                "motif_family": family, **meta,
            })

        genes = {}
        loops = set()
        for row in evidence:
            genes[row["gene_id"]] = row["gene_symbol"]
            loops.add(row["loop_id"])
        hub_gene_sets[state, region_id] = genes
        symbol_count = sum(symbol != gene_id for gene_id, symbol in genes.items())
        hub_rows.append({
            "hub_id": hub["hub_id"], "state": state, "region_id": region_id,
            "chrom": hub["chrom"], "start": hub["start"], "end": hub["end"],
            "H3K27ac_log2FC": hub["H3K27ac_log2FC"], "ATAC_log2FC": hub["ATAC_log2FC"],
            "ATAC_stratum": hub["ATAC_stratum"],
            "motif_families": ";".join(hub["families"]),
            "motif_family_count": len(hub["families"]), "target_gene_count": len(genes),
            "genes_with_symbols": symbol_count, "genes_without_symbols": len(genes) - symbol_count,
            "unique_loop_count": len(loops),
            "best_loop_q": min(float(row["loop_q"]) for row in evidence),
            "best_FIMO_p": min(float(row["best_FIMO_p"]) for row in evidence),
            "target_gene_symbols": ";".join(sorted(genes.values(), key=lambda x: (x.startswith("ENSSSCG"), x))),
        })

        per_gene = defaultdict(list)
        per_loop = defaultdict(list)
        for row in evidence:
            per_gene[row["gene_id"]].append(row)
            per_loop[(row["gene_id"], row["loop_id"])].append(row)
        for gene_id, records in sorted(per_gene.items(), key=lambda item: (item[1][0]["gene_symbol"], item[0])):
            target_rows.append({
                "hub_id": hub["hub_id"], "state": state, "region_id": region_id,
                "gene_id": gene_id, "gene_symbol": records[0]["gene_symbol"],
                "motif_families": ";".join(sorted({row["motif_family"] for row in records})),
                "unique_loop_count": len({row["loop_id"] for row in records}),
                "best_loop_q": min(float(row["loop_q"]) for row in records),
                "promoter_intervals": ";".join(sorted({
                    f"{row['promoter_chrom']}:{row['promoter_start']}-{row['promoter_end']}"
                    for row in records
                })),
                "loop_ids": ";".join(sorted({row["loop_id"] for row in records})),
            })
        for (gene_id, loop_id), records in sorted(per_loop.items()):
            first = records[0]
            loop_rows.append({
                "hub_id": hub["hub_id"], "state": state, "region_id": region_id,
                "gene_id": gene_id, "gene_symbol": first["gene_symbol"], "loop_id": loop_id,
                "motif_families": ";".join(sorted({row["motif_family"] for row in records})),
                "loop_q": first["loop_q"], "loop_p": first["loop_p"],
                "contact_count": first["contact_count"],
                "CRE_anchor_chrom": first["CRE_anchor_chrom"],
                "CRE_anchor_start": first["CRE_anchor_start"], "CRE_anchor_end": first["CRE_anchor_end"],
                "promoter_anchor_chrom": first["promoter_anchor_chrom"],
                "promoter_anchor_start": first["promoter_anchor_start"],
                "promoter_anchor_end": first["promoter_anchor_end"],
                "promoter_start": first["promoter_start"], "promoter_end": first["promoter_end"],
                "gene_tss": first["gene_tss"], "gene_strand": first["gene_strand"],
            })
    return hub_rows, family_rows, target_rows, loop_rows, hub_gene_sets


def build_site_table(run: Path, hubs: list[dict]):
    hub_lookup = {(hub["state"], hub["region_id"]): hub for hub in hubs}
    rows = []
    for state in STATES:
        path = run / "motifs" / state / "04_FIMO_site_coordinates.tsv.gz"
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                hub = hub_lookup.get((state, row["region_id"]))
                if hub is None or row["motif_family"] not in hub["families"]:
                    continue
                rows.append({"hub_id": hub["hub_id"], "state": state, **row})
    return rows


def build_pathways(run: Path, hubs: list[dict], hub_gene_sets: dict):
    term_meta = {}
    pathway_source = (
        run / "targets" / "P1_full_H3_peak_support" / "functions" / "05_TF_pathway_enrichment.tsv"
    )
    for row in read_tsv(pathway_source):
        if row["group_level"] != "program":
            continue
        state = "strength" if row["contrast"].startswith("strength") else "weaken"
        term_meta.setdefault(
            (state, row["source"], row["term_id"]),
            {
                "term_name": row["term_name"], "category": row["category"],
                "background_K": int(row["K"]), "universe_N": int(row["N"]),
                "annotation_releases": row["annotation_releases"],
            },
        )

    wanted = {gene_id for genes in hub_gene_sets.values() for gene_id in genes}
    gene_terms = {}
    function_source = (
        run / "targets" / "P1_full_H3_peak_support" / "functions" / "04_TF_candidate_functions.tsv"
    )
    for row in read_tsv(function_source):
        if row["group_level"] == "program" and row["gene_id"] in wanted and row["gene_id"] not in gene_terms:
            gene_terms[row["gene_id"]] = json.loads(row["function_terms_json"])

    results = []
    for hub in hubs:
        state, region_id = hub["state"], hub["region_id"]
        genes = hub_gene_sets[state, region_id]
        source_annotated = {source: set() for source in SOURCES}
        members = defaultdict(set)
        for gene_id in genes:
            for term in gene_terms.get(gene_id, []):
                source_annotated[term["source"]].add(gene_id)
                members[term["source"], term["term_id"]].add(gene_id)
        for source in SOURCES:
            foreground_n = len(source_annotated[source])
            if not foreground_n:
                continue
            block = []
            for (term_state, term_source, term_id), meta in term_meta.items():
                if term_state != state or term_source != source:
                    continue
                hits = members.get((source, term_id), set())
                hit_k = len(hits)
                p_value = (
                    hypergeom_sf(hit_k, meta["universe_N"], meta["background_K"], foreground_n)
                    if hit_k else 1.0
                )
                enrichment_fold = (
                    (hit_k / foreground_n) / (meta["background_K"] / meta["universe_N"])
                    if hit_k else 0.0
                )
                block.append({
                    "hub_id": hub["hub_id"], "state": state, "region_id": region_id,
                    "motif_families": ";".join(hub["families"]), "source": source,
                    "term_id": term_id, "term_name": meta["term_name"], "category": meta["category"],
                    "hit_k": hit_k, "foreground_n": foreground_n,
                    "background_K": meta["background_K"], "universe_N": meta["universe_N"],
                    "enrichment_fold": enrichment_fold, "p_value": p_value,
                    "q_BH_within_hub_source": 1.0, "q_BH_global_selected_hubs": 1.0,
                    "hit_gene_ids": ";".join(sorted(hits)),
                    "hit_gene_symbols": ";".join(sorted(genes[gene_id] for gene_id in hits)),
                    "annotation_releases": meta["annotation_releases"],
                })
            bh(block, "p_value", "q_BH_within_hub_source")
            results.extend(block)
    bh(results, "p_value", "q_BH_global_selected_hubs")
    results.sort(key=lambda row: (
        row["state"], row["hub_id"], row["source"],
        row["q_BH_global_selected_hubs"], row["p_value"], row["term_id"],
    ))
    significant = [
        row for row in results
        if row["q_BH_global_selected_hubs"] <= 0.05
        and row["enrichment_fold"] > 1 and row["hit_k"] >= 2
    ]
    return results, significant


def copy_sources(run: Path, output: Path) -> list[dict]:
    source_dir = output / "source_tables"
    source_dir.mkdir()
    sources = [
        (run / "01_selected_AME_motifs.tsv", "01_selected_AME_motifs.tsv"),
        (run / "motifs" / "strength" / "04_FIMO_site_coordinates.tsv.gz", "strength_FIMO_site_coordinates.tsv.gz"),
        (run / "motifs" / "weaken" / "04_FIMO_site_coordinates.tsv.gz", "weaken_FIMO_site_coordinates.tsv.gz"),
        (evidence_path(run, "strength"), "strength_motif_CRE_loop_gene_evidence.tsv.gz"),
        (evidence_path(run, "weaken"), "weaken_motif_CRE_loop_gene_evidence.tsv.gz"),
        (
            run / "targets" / "P1_full_H3_peak_support" / "functions" / "04_TF_candidate_functions.tsv",
            "P1_TF_candidate_functions.tsv",
        ),
        (
            run / "targets" / "P1_full_H3_peak_support" / "functions" / "05_TF_pathway_enrichment.tsv",
            "P1_TF_pathway_enrichment.tsv",
        ),
        (run / "02_promoter_definition_counts.tsv", "02_promoter_definition_counts.tsv"),
        (run / "03_target_definition_counts.tsv", "03_target_definition_counts.tsv"),
        (run / "run_metadata.json", "run_metadata.json"),
    ]
    manifest = []
    for source, name in sources:
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = source_dir / name
        shutil.copy2(source, destination)
        manifest.append({
            "source_path": str(source), "copied_file": str(destination.relative_to(output)),
            "size_bytes": destination.stat().st_size, "sha256": sha256(destination),
        })
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top-families", type=int, default=5)
    args = parser.parse_args()
    run = args.run.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    output.mkdir(parents=True)

    selected, motif_meta = load_selected_families(run, args.top_families)
    groups, strict_rows = collect_strict_evidence(run, selected)
    hubs = choose_hubs(selected, groups)
    hub_rows, family_rows, target_rows, loop_rows, hub_gene_sets = build_hub_tables(
        hubs, strict_rows, motif_meta
    )
    site_rows = build_site_table(run, hubs)
    pathway_rows, significant_rows = build_pathways(run, hubs, hub_gene_sets)

    write_tsv(output / "01_representative_enhancer_hubs.tsv", hub_rows, [
        "hub_id", "state", "region_id", "chrom", "start", "end", "H3K27ac_log2FC",
        "ATAC_log2FC", "ATAC_stratum", "motif_families", "motif_family_count",
        "target_gene_count", "genes_with_symbols", "genes_without_symbols",
        "unique_loop_count", "best_loop_q", "best_FIMO_p", "target_gene_symbols",
    ])
    write_tsv(output / "02_hub_motif_families.tsv", family_rows, [
        "hub_id", "state", "region_id", "motif_family", "AME_rank", "motif_id",
        "motif_name", "AME_q", "enrichment_ratio", "log2_enrichment",
    ])
    write_tsv(output / "03_hub_target_genes.tsv", target_rows, [
        "hub_id", "state", "region_id", "gene_id", "gene_symbol", "motif_families",
        "unique_loop_count", "best_loop_q", "promoter_intervals", "loop_ids",
    ])
    write_tsv(output / "04_hub_loop_links.tsv", loop_rows, [
        "hub_id", "state", "region_id", "gene_id", "gene_symbol", "loop_id",
        "motif_families", "loop_q", "loop_p", "contact_count", "CRE_anchor_chrom",
        "CRE_anchor_start", "CRE_anchor_end", "promoter_anchor_chrom",
        "promoter_anchor_start", "promoter_anchor_end", "promoter_start", "promoter_end",
        "gene_tss", "gene_strand",
    ])
    write_tsv(output / "05_hub_motif_sites.tsv", site_rows, [
        "hub_id", "state", "region_id", "core_id", "motif_id", "motif_family",
        "chrom", "start", "end", "strand", "FIMO_p", "site_inside_full_CRE",
    ])
    pathway_fields = [
        "hub_id", "state", "region_id", "motif_families", "source", "term_id",
        "term_name", "category", "hit_k", "foreground_n", "background_K", "universe_N",
        "enrichment_fold", "p_value", "q_BH_within_hub_source",
        "q_BH_global_selected_hubs", "hit_gene_ids", "hit_gene_symbols", "annotation_releases",
    ]
    write_tsv(output / "06_hub_pathways_all_tested.tsv", pathway_rows, pathway_fields)
    write_tsv(output / "07_hub_pathways_significant.tsv", significant_rows, pathway_fields)
    manifest_rows = copy_sources(run, output)
    write_tsv(output / "08_source_file_manifest.tsv", manifest_rows, [
        "source_path", "copied_file", "size_bytes", "sha256",
    ])

    readme = f"""Unbiased enhancer-motif-gene pathway plotting bundle

Source run: {run}
Promoter definition: P1_full_H3_peak_support
States: Strength uses 6 h promoters/loops; Weaken uses 0 h promoters/loops.

Representative selection
1. Families were ranked only by AME log2 enrichment; top {args.top_families} per state were retained.
2. Evidence required any_site_overlaps_CRE_anchor=1.
3. The representative enhancer for each family was ranked by target-gene count, loop count,
   absolute H3K27ac log2FC, loop q, FIMO P and region ID.
4. Families selecting the same enhancer were merged into one hub.
5. No gene name, manual lipid panel or pathway label was used in selection.

Pathway testing
GO BP and Reactome over-representation used the corresponding state-specific P1 callable-gene
background recorded in run_42014. BH q values are supplied within each hub/source and globally
across every hub/source/term test in this bundle. The significant table requires global q <= 0.05,
enrichment fold > 1 and at least two hit genes.

Plot-ready tables
01_representative_enhancer_hubs.tsv: one row per unique enhancer hub.
02_hub_motif_families.tsv: AME statistics and motif members of each hub.
03_hub_target_genes.tsv: one row per hub-target gene.
04_hub_loop_links.tsv: one row per unique strict hub-gene-loop link, with both anchors.
05_hub_motif_sites.tsv: exact FIMO coordinates for selected hub motifs.
06_hub_pathways_all_tested.tsv: all pathway hypotheses and adjusted P values.
07_hub_pathways_significant.tsv: globally significant plot-ready pathway subset.
08_source_file_manifest.tsv: copied source paths, sizes and SHA-256 values.

Coordinate convention: BED-style zero-based, half-open.
Interpretation: motif plus loop supports a putative physical association, not TF occupancy,
causal regulation or target-gene expression change.
"""
    (output / "README.txt").write_text(readme, encoding="utf-8", newline="\n")
    completion = {
        "source_run": str(run), "promoter_definition": "P1_full_H3_peak_support",
        "top_families_per_state": args.top_families, "hub_count": len(hub_rows),
        "hub_family_rows": len(family_rows), "hub_target_gene_rows": len(target_rows),
        "hub_loop_rows": len(loop_rows), "hub_motif_site_rows": len(site_rows),
        "pathway_tests": len(pathway_rows), "significant_pathways": len(significant_rows),
    }
    (output / "completed.json").write_text(
        json.dumps(completion, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n"
    )


if __name__ == "__main__":
    main()
