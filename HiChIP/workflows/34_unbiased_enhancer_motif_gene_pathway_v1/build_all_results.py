#!/usr/bin/env python3
"""Export every enriched motif family and every strict P1 CRE-loop-gene link.

Unlike build_plot_data.py, this module does not subset motif families and does
not restrict the main output to one enhancer per family.  A separate ranking
table is emitted only as an optional plotting index.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from collections import defaultdict
from pathlib import Path

from build_plot_data import (
    STATES,
    SOURCES,
    bh,
    collect_strict_evidence,
    copy_sources,
    hypergeom_sf,
    load_selected_families,
    read_tsv,
    write_tsv,
)


def all_family_rows(run: Path, motif_meta: dict) -> list[dict]:
    rows = []
    for state in STATES:
        state_rows = [row for row in read_tsv(run / "01_selected_AME_motifs.tsv")
                      if row["contrast"].startswith(state)]
        state_rows.sort(key=lambda row: (-float(row["log2_enrichment"]), row["motif_family"]))
        for rank, row in enumerate(state_rows, 1):
            rows.append({
                "state": state, "AME_rank_within_state": rank,
                "contrast": row["contrast"], "motif_family": row["motif_family"],
                "motif_id": row["motif_id"], "motif_name": row["motif_name"],
                "AME_q": row["AME_q"], "enrichment_ratio": row["enrichment_ratio"],
                "log2_enrichment": row["log2_enrichment"],
            })
    return rows


def strict_link_tables(groups: dict, strict_rows: dict, motif_meta: dict):
    enhancer_rows, gene_rows, loop_rows = [], [], []
    for key in sorted(groups, key=lambda x: (x[0], motif_meta[x[0], x[1]]["AME_rank"], x[2])):
        state, family, region_id = key
        group, evidence = groups[key], strict_rows[key]
        enhancer_rows.append({
            "state": state, "motif_family": family,
            "AME_rank_within_state": motif_meta[state, family]["AME_rank"],
            "region_id": region_id, "chrom": group["chrom"], "start": group["start"],
            "end": group["end"], "H3K27ac_log2FC": group["H3K27ac_log2FC"],
            "ATAC_log2FC": group["ATAC_log2FC"], "ATAC_stratum": group["ATAC_stratum"],
            "target_gene_count": len(group["genes"]), "unique_loop_count": len(group["loops"]),
            "best_loop_q": group["best_loop_q"], "best_FIMO_p": group["best_FIMO_p"],
            "target_gene_symbols": ";".join(sorted(group["genes"].values(),
                                                     key=lambda x: (x.startswith("ENSSSCG"), x))),
        })
        per_gene = defaultdict(list)
        per_loop = defaultdict(list)
        for row in evidence:
            per_gene[row["gene_id"]].append(row)
            per_loop[row["gene_id"], row["loop_id"]].append(row)
        for gene_id, records in sorted(per_gene.items(), key=lambda item: (item[1][0]["gene_symbol"], item[0])):
            gene_rows.append({
                "state": state, "motif_family": family, "region_id": region_id,
                "gene_id": gene_id, "gene_symbol": records[0]["gene_symbol"],
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
                "state": state, "motif_family": family, "region_id": region_id,
                "gene_id": gene_id, "gene_symbol": first["gene_symbol"], "loop_id": loop_id,
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
    return enhancer_rows, gene_rows, loop_rows


def selected_site_rows(run: Path, groups: dict, strict_rows: dict) -> list[dict]:
    linked_pairs = {(state, family, region_id) for state, family, region_id in groups}
    anchors = defaultdict(set)
    for key, evidence in strict_rows.items():
        for row in evidence:
            anchors[key].add((row["CRE_anchor_chrom"], int(row["CRE_anchor_start"]),
                              int(row["CRE_anchor_end"])))
    rows = []
    for state in STATES:
        path = run / "motifs" / state / "04_FIMO_site_coordinates.tsv.gz"
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (state, row["motif_family"], row["region_id"])
                if key not in linked_pairs:
                    continue
                start, end = int(row["start"]), int(row["end"])
                overlaps = any(
                    row["chrom"] == chrom and start < anchor_end and anchor_start < end
                    for chrom, anchor_start, anchor_end in anchors[key]
                )
                rows.append({"state": state, **row,
                             "overlaps_any_linked_CRE_anchor": int(overlaps)})
    return rows


def representative_index(groups: dict, motif_meta: dict) -> list[dict]:
    rows = []
    for state in STATES:
        families = sorted(
            {family for candidate_state, family, _ in groups if candidate_state == state},
            key=lambda family: motif_meta[state, family]["AME_rank"],
        )
        for family in families:
            candidates = [value for (candidate_state, candidate_family, _), value in groups.items()
                          if candidate_state == state and candidate_family == family]
            candidates.sort(key=lambda row: (
                -len(row["genes"]), -len(row["loops"]), -abs(row["H3K27ac_log2FC"]),
                row["best_loop_q"], row["best_FIMO_p"], row["region_id"],
            ))
            winner = candidates[0]
            rows.append({
                "state": state, "motif_family": family,
                "AME_rank_within_state": motif_meta[state, family]["AME_rank"],
                "region_id": winner["region_id"], "chrom": winner["chrom"],
                "start": winner["start"], "end": winner["end"],
                "H3K27ac_log2FC": winner["H3K27ac_log2FC"],
                "ATAC_log2FC": winner["ATAC_log2FC"],
                "target_gene_count": len(winner["genes"]),
                "unique_loop_count": len(winner["loops"]),
                "selection_note": "optional_plotting_index_only;not_a_main-result_filter",
            })
    return rows


def strict_family_pathways(run: Path, selected: dict):
    term_meta = {}
    pathway_source = (
        run / "targets" / "P1_full_H3_peak_support" / "functions" / "05_TF_pathway_enrichment.tsv"
    )
    for row in read_tsv(pathway_source):
        if row["group_level"] != "program":
            continue
        state = "strength" if row["contrast"].startswith("strength") else "weaken"
        term_meta.setdefault((state, row["source"], row["term_id"]), {
            "term_name": row["term_name"], "category": row["category"],
            "background_K": int(row["K"]), "universe_N": int(row["N"]),
            "annotation_releases": row["annotation_releases"],
        })

    family_genes = defaultdict(dict)
    gene_terms = {}
    function_source = (
        run / "targets" / "P1_full_H3_peak_support" / "functions" / "04_TF_candidate_functions.tsv"
    )
    for row in read_tsv(function_source):
        if row["group_level"] != "program" or not row["tf_label"].startswith("AME|"):
            continue
        state = "strength" if row["contrast"].startswith("strength") else "weaken"
        family = row["tf_label"][4:]
        if family not in selected[state] or int(row["strict_site_anchor_unique_CRE_count"]) <= 0:
            continue
        family_genes[state, family][row["gene_id"]] = row["gene_symbol"]
        gene_terms.setdefault(row["gene_id"], json.loads(row["function_terms_json"]))

    summary, results = [], []
    for state in STATES:
        for family in selected[state]:
            genes = family_genes[state, family]
            summary.append({
                "state": state, "motif_family": family,
                "strict_target_gene_count": len(genes),
                "genes_with_symbols": sum(symbol != gene_id for gene_id, symbol in genes.items()),
                "genes_without_symbols": sum(symbol == gene_id for gene_id, symbol in genes.items()),
                "target_gene_symbols": ";".join(sorted(genes.values(),
                                                         key=lambda x: (x.startswith("ENSSSCG"), x))),
            })
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
                    fold = (
                        (hit_k / foreground_n) / (meta["background_K"] / meta["universe_N"])
                        if hit_k else 0.0
                    )
                    block.append({
                        "state": state, "motif_family": family, "source": source,
                        "term_id": term_id, "term_name": meta["term_name"],
                        "category": meta["category"], "hit_k": hit_k,
                        "foreground_n": foreground_n, "background_K": meta["background_K"],
                        "universe_N": meta["universe_N"], "enrichment_fold": fold,
                        "p_value": p_value, "q_BH_within_family_source": 1.0,
                        "q_BH_global_all_families": 1.0,
                        "hit_gene_ids": ";".join(sorted(hits)),
                        "hit_gene_symbols": ";".join(sorted(genes[gene_id] for gene_id in hits)),
                        "annotation_releases": meta["annotation_releases"],
                    })
                bh(block, "p_value", "q_BH_within_family_source")
                results.extend(block)
    bh(results, "p_value", "q_BH_global_all_families")
    results.sort(key=lambda row: (
        row["state"], row["motif_family"], row["source"],
        row["q_BH_global_all_families"], row["p_value"], row["term_id"],
    ))
    significant = [row for row in results
                   if row["q_BH_global_all_families"] <= 0.05
                   and row["enrichment_fold"] > 1 and row["hit_k"] >= 2]
    return summary, results, significant


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    # Keep an ASCII workspace junction un-resolved on Windows. Resolving it to
    # the original path can pass non-ASCII command-line paths through a legacy
    # code page in some local Python launches.
    run = Path(os.path.abspath(args.run))
    output = Path(os.path.abspath(args.output))
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)

    selected, motif_meta = load_selected_families(run, 1_000_000)
    groups, strict_rows = collect_strict_evidence(run, selected)
    family_rows = all_family_rows(run, motif_meta)
    enhancer_rows, gene_rows, loop_rows = strict_link_tables(groups, strict_rows, motif_meta)
    site_rows = selected_site_rows(run, groups, strict_rows)
    summary_rows, pathway_rows, significant_rows = strict_family_pathways(run, selected)
    representative_rows = representative_index(groups, motif_meta)

    write_tsv(output / "01_all_enriched_motif_families.tsv", family_rows, [
        "state", "AME_rank_within_state", "contrast", "motif_family", "motif_id",
        "motif_name", "AME_q", "enrichment_ratio", "log2_enrichment",
    ])
    write_tsv(output / "02_all_strict_enhancer_motif_summary.tsv", enhancer_rows, [
        "state", "motif_family", "AME_rank_within_state", "region_id", "chrom", "start",
        "end", "H3K27ac_log2FC", "ATAC_log2FC", "ATAC_stratum", "target_gene_count",
        "unique_loop_count", "best_loop_q", "best_FIMO_p", "target_gene_symbols",
    ])
    write_tsv(output / "03_all_strict_enhancer_motif_target_genes.tsv", gene_rows, [
        "state", "motif_family", "region_id", "gene_id", "gene_symbol",
        "unique_loop_count", "best_loop_q", "promoter_intervals", "loop_ids",
    ])
    write_tsv(output / "04_all_strict_enhancer_motif_loop_links.tsv", loop_rows, [
        "state", "motif_family", "region_id", "gene_id", "gene_symbol", "loop_id",
        "loop_q", "loop_p", "contact_count", "CRE_anchor_chrom", "CRE_anchor_start",
        "CRE_anchor_end", "promoter_anchor_chrom", "promoter_anchor_start",
        "promoter_anchor_end", "promoter_start", "promoter_end", "gene_tss", "gene_strand",
    ])
    write_tsv(output / "05_all_linked_motif_sites.tsv", site_rows, [
        "state", "region_id", "core_id", "motif_id", "motif_family", "chrom", "start",
        "end", "strand", "FIMO_p", "site_inside_full_CRE", "overlaps_any_linked_CRE_anchor",
    ])
    write_tsv(output / "06_motif_family_strict_target_summary.tsv", summary_rows, [
        "state", "motif_family", "strict_target_gene_count", "genes_with_symbols",
        "genes_without_symbols", "target_gene_symbols",
    ])
    pathway_fields = [
        "state", "motif_family", "source", "term_id", "term_name", "category", "hit_k",
        "foreground_n", "background_K", "universe_N", "enrichment_fold", "p_value",
        "q_BH_within_family_source", "q_BH_global_all_families", "hit_gene_ids",
        "hit_gene_symbols", "annotation_releases",
    ]
    write_tsv(output / "07_motif_family_pathways_all_tested.tsv", pathway_rows, pathway_fields)
    write_tsv(output / "08_motif_family_pathways_significant.tsv", significant_rows, pathway_fields)
    write_tsv(output / "09_optional_representative_enhancer_per_family.tsv", representative_rows, [
        "state", "motif_family", "AME_rank_within_state", "region_id", "chrom", "start",
        "end", "H3K27ac_log2FC", "ATAC_log2FC", "target_gene_count",
        "unique_loop_count", "selection_note",
    ])
    manifest = copy_sources(run, output)
    write_tsv(output / "10_source_file_manifest.tsv", manifest, [
        "source_path", "copied_file", "size_bytes", "sha256",
    ])

    readme = f"""All enriched motif families: strict P1 CRE-loop-gene plotting bundle

This is the corrected no-family-subsetting version.

Main rules
- All AME-significant enriched motif families are retained: {len(family_rows)} total.
- All enhancer-motif-gene rows require the FIMO motif site itself to overlap the CRE-side
  HiChIP anchor (any_site_overlaps_CRE_anchor=1).
- No gene name, pathway label or manual reference panel is used to retain a row.
- Files 01-08 are the main complete results. File 09 is only an optional plotting index;
  it does not filter or define the main result.

Files
01: all significant enriched motif families and AME statistics.
02: every strict enhancer-motif pair, including target/loop counts.
03: every strict enhancer-motif-target-gene association.
04: every strict enhancer-motif-gene-loop association with both anchors.
05: FIMO sites for every linked enhancer-motif pair; the final column identifies sites
    overlapping at least one linked CRE anchor.
06: strict target-gene set for every motif family.
07: all GO BP and Reactome tests for strict motif-family target sets.
08: globally significant terms (global BH q <= 0.05, fold > 1, hit genes >= 2).
09: one evidence-ranked enhancer per family for optional compact figures only.
10: source file checksums.

Source run: {run}
Promoter definition: P1_full_H3_peak_support
Coordinate convention: zero-based, half-open.
"""
    (output / "README.txt").write_text(readme, encoding="utf-8", newline="\n")
    completion = {
        "source_run": str(run), "motif_family_count": len(family_rows),
        "strict_enhancer_motif_rows": len(enhancer_rows),
        "strict_enhancer_motif_gene_rows": len(gene_rows),
        "strict_enhancer_motif_loop_rows": len(loop_rows),
        "linked_motif_site_rows": len(site_rows), "motif_family_summary_rows": len(summary_rows),
        "pathway_tests": len(pathway_rows), "significant_pathways": len(significant_rows),
        "optional_representative_rows": len(representative_rows),
    }
    (output / "completed.json").write_text(
        json.dumps(completion, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n"
    )


if __name__ == "__main__":
    main()
