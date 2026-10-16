#!/usr/bin/env python3
"""Offline exact bigWig summaries and exploratory matched Pol II comparisons.

Existing CPM tracks use bamCoverage --skipNAs: gaps within a valid contig
are zero OBSERVED coverage, not proof of absent binding. Missing contigs stay NA.
P-value track gaps stay missing, and regional mean scores are not region P-values.
No peak calling, no FDR test, no automatic 'specific binding' classification.
"""
import argparse
import bisect
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import random
import statistics
import sys


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str = "."
    gene_id: str = "."
    strand: str = "."
    offset: str = "."

    @property
    def key(self):
        return self.chrom, self.start, self.end


def read_bed(path):
    regions = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line or line.startswith(("#", "track", "browser")):
            continue
        cells = line.split("\t")
        if len(cells) < 3:
            raise ValueError(f"Bad BED at {path}:{line_no}")
        chrom, start, end = cells[0], int(cells[1]), int(cells[2])
        if start < 0 or end <= start:
            raise ValueError(f"Invalid interval at {path}:{line_no}")
        regions.append(Region(chrom, start, end, *(cells[3:6] + ["."] * 3)[:3]))
    return regions


def clip_region(chrom, start, end, sizes, **metadata):
    if chrom not in sizes:
        raise ValueError(f"Reference contig not found: {chrom}")
    start, end = max(0, start), min(sizes[chrom], end)
    if start >= end:
        raise ValueError(f"Empty clipped interval: {chrom}:{start}-{end}")
    return Region(chrom, start, end, **metadata)


def tss_boundary(gene):
    if gene.strand not in ("+", "-"):
        raise ValueError(f"Gene strand not specified: {gene}")
    return gene.start if gene.strand == "+" else gene.end


def tss_bin(gene, offset, width, sizes):
    anchor = tss_boundary(gene)
    if gene.strand == "+":
        start, end = anchor + offset, anchor + offset + width
    else:
        start, end = anchor - offset - width, anchor - offset
    return clip_region(gene.chrom, start, end, sizes, name=gene.name,
                       gene_id=gene.gene_id, strand=gene.strand, offset=str(offset))


def integrate(intervals, start, end, kind):
    """Exact length-weighted integration, no zoom levels, no gap interpolation."""
    if end <= start:
        raise ValueError("Non-positive interval length")
    total = 0.0
    represented = positive = p2 = p5 = 0
    maximum = None
    previous_end = start
    nonfinite = 0
    for left, right, value in intervals or ():
        left, right = max(start, left), min(end, right)
        if right <= left:
            continue
        if left < previous_end:
            raise ValueError("Overlapping or unsorted bigWig entries")
        previous_end = right
        width = right - left
        if not math.isfinite(value):
            nonfinite += width
            continue
        if value < 0:
            raise ValueError(f"Negative score in {kind} bigWig: {value}")
        total += width * value
        represented += width
        positive += width if value > 0 else 0
        p2 += width if value >= 2 else 0
        p5 += width if value >= 5 else 0
        maximum = value if maximum is None else max(value, maximum)
    length = end - start
    if kind == "CPM":
        # Explicit NaN entries are not the documented no-read gaps: flag them.
        mean = None if nonfinite else total / length
        peak = maximum if maximum is not None else (None if nonfinite else 0.0)
    else:
        mean = total / length if represented == length else None
        peak = maximum
    return {"mean_full_interval": mean,
            "mean_encoded_bases": total / represented if represented else None,
            "encoded_fraction": represented / length,
            "positive_fraction": positive / length, "maximum": peak,
            "fraction_bp_score_ge_2": p2 / length if kind != "CPM" else None,
            "fraction_bp_score_ge_5": p5 / length if kind != "CPM" else None,
            "status": "NONFINITE_ENTRIES" if nonfinite else "OK" if represented == length else "GAPS_ZERO_OBSERVED_CPM" if kind == "CPM" else "GAPS_NA_PPOIS"}


def unique_file(root, filename):
    matches = list(root.rglob(filename))
    if len(matches) != 1 or not matches[0].is_file() or matches[0].stat().st_size == 0:
        raise ValueError(f"Expected one nonempty {filename} under {root}; found {len(matches)}")
    return matches[0]


def prepare(project, reference):
    results = project / "results_pig_v2"
    screen = results / "hypothesis_screen_current_peaks"
    fai = reference / "genome.fa.fai"
    sizes = {}
    for line in fai.read_text().splitlines():
        cells = line.split("\t")
        sizes[cells[0]] = int(cells[1])
    genes = read_bed(screen / "annotations/Ensembl115.genes.bed")
    candidates = read_bed(screen / "annotations/candidate_genes.bed")
    groups = {"rRNA": read_bed(screen / "annotations/Ensembl115.rRNA.bed"),
              "5S_rRNA": read_bed(screen / "annotations/Ensembl115.5S_rRNA.bed"),
              "candidate_gene_bodies": candidates,
              "candidate_genes_plusminus2kb": read_bed(screen / "annotations/candidate_genes.plusminus2kb.bed")}
    for flank in (500, 2000):
        groups[f"candidate_TSS_plusminus{flank}bp"] = [clip_region(g.chrom, tss_boundary(g) - flank, tss_boundary(g) + flank, sizes,
                                                                 name=g.name, gene_id=g.gene_id, strand=g.strand) for g in candidates]
    groups["candidate_TSS_100bp_bins"] = [tss_bin(g, offset, 100, sizes) for g in candidates for offset in range(-2000, 2000, 100)]
    rrna_bins = []
    for index, r in enumerate(groups["rRNA"], 1):
        for start in range(max(0, r.start - 1000), min(sizes[r.chrom], r.end + 1000), 100):
            rrna_bins.append(clip_region(r.chrom, start, min(start + 100, r.end + 1000), sizes,
                                        name=f"rRNA_{index}:{r.chrom}:{r.start}-{r.end}", offset=str(start-r.start)))
    groups["rRNA_and_flanks_100bp_bins"] = rrna_bins
    for caller in ("macs2", "seacr"):
        for state in (119, 137):
            groups[f"{caller}_{state}_Tag_shared"] = read_bed(screen / f"peak_sets/{caller}_{state}_Tag_rep1-rep2_shared.bed")
        groups[f"{caller}_119_Pol2_shared"] = read_bed(screen / f"peak_sets/{caller}_119_Pol2_rep1-rep2_shared.bed")
    tracks = []
    for state in (119, 137):
        root = results / "03_peak_calling/03_bed_to_bigwig"
        for mark, rep in (("Tag", 1), ("Tag", 2), ("IgG", 1), ("Pol2", 1), ("Pol2", 2)):
            label = f"{state}_{mark}" if mark == "IgG" else f"{state}_{mark}_rep{rep}"
            tracks.append((label, "CPM", unique_file(root, f"3D4_{state}_{mark}_R{rep}.bigWig")))
        for mark in ("Tag", "Pol2"):
            filename = f"3D4_{state}_{mark}_rep1-rep2_pooled_vs_3D4_{state}_IgG_rep1.MACS2_ppois_Sscrofa11.1.bigWig"
            tracks.append((f"{state}_{mark}_pooled_vs_IgG", "PPOIS", unique_file(results / "custom_tracks/macs2_ppois", filename)))
    for name, rows in groups.items():
        for r in rows:
            if r.chrom not in sizes or r.end > sizes[r.chrom]:
                raise ValueError(f"BED/reference disagreement: {name} {r}")
    return sizes, genes, groups, tracks


class Reader:
    def __init__(self, tracks, sizes, pybigwig):
        self.handles = {}
        self.chroms = {}
        self.kinds = {}
        self.manifest = []
        try:
            for label, kind, path in tracks:
                handle = pybigwig.open(str(path))
                if handle is None or not handle.isBigWig():
                    raise ValueError(f"Not a readable bigWig: {path}")
                self.handles[label] = handle
                chromosomes = handle.chroms()
                for chrom, length in chromosomes.items():
                    if chrom not in sizes or length != sizes[chrom]:
                        raise ValueError(f"Reference/header mismatch: {label} {chrom}:{length}")
                self.chroms[label] = chromosomes
                self.kinds[label] = kind
                self.manifest.append({"label": label, "kind": kind, "path": str(path), "size_bytes": path.stat().st_size,
                                      "mtime_ns": path.stat().st_mtime_ns, "missing_contigs": sorted(set(sizes)-set(chromosomes))})
        except Exception:
            self.close()
            raise

    def close(self):
        for handle in self.handles.values():
            handle.close()

    def summary(self, label, r):
        if r.chrom not in self.chroms[label]:
            return {k: None for k in ("mean_full_interval", "mean_encoded_bases", "encoded_fraction", "positive_fraction", "maximum", "fraction_bp_score_ge_2", "fraction_bp_score_ge_5")} | {"status": "MISSING_CONTIG"}
        return integrate(self.handles[label].intervals(r.chrom, r.start, r.end), r.start, r.end, self.kinds[label])


def fmt(value):
    return "NA" if value is None else format(value, ".9g") if isinstance(value, float) else value


def comparisons(means):
    result = {}
    for state in (119, 137):
        igg = means[f"{state}_IgG"]
        for mark in ("Tag", "Pol2"):
            vals = [means[f"{state}_{mark}_rep{rep}"] for rep in (1, 2)]
            prefix = f"{state}_{mark}"
            result[f"{prefix}_both_above_IgG_descriptive"] = None if igg is None or any(v is None for v in vals) else int(all(v > igg for v in vals))
            for rep, v in zip((1, 2), vals):
                base = f"{prefix}_rep{rep}"
                result[f"{base}_minus_IgG"] = None if v is None or igg is None else v - igg
                result[f"{base}_ratio_no_pseudocount"] = None if v is None or igg is None or igg == 0 else v / igg
                for pseudocount in (.01, .1):
                    result[f"{base}_log2ratio_pc{pseudocount}"] = None if v is None or igg is None else math.log2((v + pseudocount)/(igg + pseudocount))
    return result


def write_group(path, regions, reader):
    labels = list(reader.handles)
    comp_names = list(comparisons({x: 0.0 for x in labels if reader.kinds[x] == "CPM"}))
    fields = ["chrom", "start", "end", "length_bp", "name", "gene_id", "strand", "offset_bp"]
    for label in labels:
        for metric in ("mean_full_interval", "mean_encoded_bases", "encoded_fraction", "positive_fraction", "maximum", "fraction_bp_score_ge_2", "fraction_bp_score_ge_5", "status"):
            fields.append(f"{label}.{metric}")
    fields += comp_names
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for r in regions:
            row = [r.chrom, r.start, r.end, r.end-r.start, r.name, r.gene_id, r.strand, r.offset]
            means = {}
            for label in labels:
                data = reader.summary(label, r)
                row.extend(data[m] for m in ("mean_full_interval", "mean_encoded_bases", "encoded_fraction", "positive_fraction", "maximum", "fraction_bp_score_ge_2", "fraction_bp_score_ge_5", "status"))
                if reader.kinds[label] == "CPM":
                    means[label] = data["mean_full_interval"]
            row.extend(comparisons(means).values())
            writer.writerow([fmt(x) for x in row])


def overlap_index(regions):
    by_chrom = {}
    for r in regions:
        by_chrom.setdefault(r.chrom, []).append((r.start, r.end))
    result = {}
    for chrom, values in by_chrom.items():
        merged = []
        for start, end in sorted(values):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        result[chrom] = ([x[0] for x in merged], [x[1] for x in merged])
    return result


def overlaps(r, index):
    starts, ends = index.get(r.chrom, ([], []))
    i = bisect.bisect_left(starts, r.end) - 1
    return i >= 0 and ends[i] > r.start


def closest_tss(r, genes_index):
    positions, rows = genes_index.get(r.chrom, ([], []))
    if not positions:
        return None
    midpoint = (r.start+r.end)/2
    i = bisect.bisect_left(positions, midpoint)
    possibilities = [j for j in (i-1, i) if 0 <= j < len(rows)]
    return rows[min(possibilities, key=lambda j: (abs(positions[j]-midpoint), j))]


def gc_stats(fasta, r):
    seq = fasta.fetch(r.chrom, r.start, r.end).upper()
    valid = sum(seq.count(base) for base in "ACGT")
    return ((seq.count("G") + seq.count("C"))/valid if valid else None, valid/(r.end-r.start))


def caliper(a, b, pseudocount=.01):
    return a is not None and b is not None and abs(math.log2((a+pseudocount)/(b+pseudocount))) <= 1


def matched_regions(sources, genes, sizes, reader, fasta, pol2_index, controls_per_site, attempts, seed):
    """Other Pol II regions, NOT validated negatives; no selection on Tag signal."""
    rng = random.Random(seed)
    index, pools = {}, {}
    for gene in genes:
        if gene.strand in ("+", "-"):
            index.setdefault(gene.chrom, []).append(gene)
            pools.setdefault((gene.chrom, gene.strand), []).append(gene)
    for chrom, rows in index.items():
        rows.sort(key=lambda g: (tss_boundary(g), g.gene_id))
        index[chrom] = ([tss_boundary(g) for g in rows], rows)
    controls, info = [], []
    match_labels = ("119_IgG", "119_Pol2_rep1", "119_Pol2_rep2")
    for number, source in enumerate(sources, 1):
        source_id = f"cooccupancy_{number}:{source.chrom}:{source.start}-{source.end}"
        gene = closest_tss(source, index)
        record = {"source_id": source_id, "source_interval": list(source.key), "matched_controls": 0,
                  "attempts": 0, "status": "INSUFFICIENT_MATCHES"}
        if gene is None:
            record["status"] = "NO_TSS_ANNOTATION"
            info.append(record)
            continue
        offset = source.start - tss_boundary(gene)
        source_distance = abs((source.start+source.end)/2 - tss_boundary(gene))
        source_gc, source_valid = gc_stats(fasta, source)
        record.update(nearest_TSS_gene=gene.name, distance_to_TSS=source_distance, source_GC=source_gc, source_ACGT_fraction=source_valid)
        source_means = [reader.summary(label, source)["mean_full_interval"] for label in match_labels]
        if source_gc is None or any(x is None for x in source_means):
            record["status"] = "SOURCE_MISSING_SIGNAL_OR_SEQUENCE"
            info.append(record)
            continue
        pool = pools[(source.chrom, gene.strand)]
        seen = set()
        for attempt in range(attempts):
            record["attempts"] = attempt+1
            target = rng.choice(pool)
            start = tss_boundary(target)+offset
            end = start + source.end-source.start
            if start < 0 or end > sizes[source.chrom] or (start < source.end+10000 and end > source.start-10000):
                continue
            key = source.chrom, start, end
            if key in seen:
                continue
            seen.add(key)
            r = Region(*key, name=source_id, gene_id=target.gene_id, strand=target.strand)
            if not overlaps(r, pol2_index):
                continue
            actual_gene = closest_tss(r, index)
            distance = abs((start+end)/2-tss_boundary(actual_gene))
            if actual_gene.strand != gene.strand or abs(distance-source_distance) > max(100, source_distance*.1):
                continue
            gc, valid = gc_stats(fasta, r)
            if gc is None or abs(gc-source_gc) > .05 or abs(valid-source_valid) > .05:
                continue
            values = [reader.summary(label, r)["mean_full_interval"] for label in match_labels]
            if not all(caliper(a, b) for a, b in zip(source_means, values)):
                continue
            controls.append(r)
            record["matched_controls"] += 1
            if record["matched_controls"] == controls_per_site:
                record["status"] = "MATCHED_EXPLORATORY_COMPARATORS"
                break
        info.append(record)
        print(f"{source_id}: {record['matched_controls']}/{controls_per_site} matched comparison regions", flush=True)
    return controls, info


def matched_summary(sources, controls, match_info, reader):
    result = []
    for source, info in zip(sources, match_info):
        row = dict(info)
        matched = [r for r in controls if r.name == info["source_id"]]
        for rep in (1, 2):
            label = f"119_Tag_rep{rep}"
            foreground = reader.summary(label, source)["mean_full_interval"]
            fg_igg = reader.summary("119_IgG", source)["mean_full_interval"]
            values = [(reader.summary(label, r)["mean_full_interval"], reader.summary("119_IgG", r)["mean_full_interval"]) for r in matched]
            for pc in (.01, .1):
                fg = None if foreground is None or fg_igg is None else math.log2((foreground+pc)/(fg_igg+pc))
                bg = [math.log2((v+pc)/(igg+pc)) for v, igg in values if v is not None and igg is not None]
                row[f"rep{rep}_log2ratio_pc{pc}"] = fg
                row[f"rep{rep}_comparator_median_pc{pc}"] = statistics.median(bg) if bg else None
                # Descriptive rank only. Source was selected using Tag peak calls.
                row[f"rep{rep}_descriptive_percentile_pc{pc}"] = (sum(v < fg for v in bg)+.5*sum(v == fg for v in bg))/len(bg) if bg and fg is not None else None
        row["interpretation"] = "POST_SELECTION_DESCRIPTIVE_ONLY_NOT_PVALUE_OR_FDR"
        result.append(row)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--controls-per-site", type=int, default=20)
    parser.add_argument("--match-attempts", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261009)
    args = parser.parse_args()
    if args.controls_per_site < 1 or args.match_attempts < 1:
        parser.error("Control counts and attempts must be positive")
    import pyBigWig
    import pysam
    sizes, genes, groups, tracks = prepare(args.project, args.reference)
    reader = Reader(tracks, sizes, pyBigWig)
    fasta = None
    try:
        fasta = pysam.FastaFile(str(args.reference / "genome.fa"))
        if dict(zip(fasta.references, fasta.lengths)) != sizes:
            raise ValueError("Indexed FASTA and FAI disagree")
        pol2_index = overlap_index(groups["macs2_119_Pol2_shared"])
        sources = [r for r in groups["macs2_119_Tag_shared"] if overlaps(r, pol2_index)]
        groups["macs2_119_Tag_Pol2_cooccupancy_Tag_intervals"] = sources
        print(f"Input check OK: {len(tracks)} bigWigs; {len(sources)} strict MACS2 cooccupancy sites; {sum(map(len, groups.values()))} intervals; pyBigWig={pyBigWig.__version__}, pysam={pysam.__version__}", flush=True)
        if args.check:
            return
        if args.out is None:
            parser.error("--out required unless --check")
        args.out.mkdir(parents=True, exist_ok=False)
        controls, match_info = matched_regions(sources, genes, sizes, reader, fasta, pol2_index, args.controls_per_site, args.match_attempts, args.seed)
        groups["matched_Pol2_comparison_regions"] = controls
        for name, regions in groups.items():
            print(f"Summarizing {name}: {len(regions)} intervals", flush=True)
            write_group(args.out / f"{name}.exact_signal.tsv", regions, reader)
        summary = matched_summary(sources, controls, match_info, reader)
        (args.out / "cooccupancy_matched_comparison.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        manifest = {"seed": args.seed, "tracks": reader.manifest, "group_counts": {name: len(rows) for name, rows in groups.items()},
                    "matching": {"same_chromosome": True, "same_width": True, "same_TSS_strand": True, "distance_to_nearest_TSS_tolerance_bp": "max(100, 0.1*source_distance)",
                                 "GC_tolerance": .05, "ACGT_fraction_tolerance": .05, "IgG_and_each_Pol2_rep_log2_caliper": 1,
                                 "matching_pseudocount_CPM": .01, "minimum_distance_from_source_bp": 10000,
                                 "select_on_Tag_signal": False, "controls_per_site": args.controls_per_site, "max_attempts": args.match_attempts},
                    "versions": {"pyBigWig": pyBigWig.__version__, "pysam": pysam.__version__},
                    "limits": ["Not a formal enrichment test: source sites were selected using Tag peak calls.", "Matched regions may be truly Tag-bound; they are comparators, not validated negative controls.",
                               "No ATAC, KO, no-primary-antibody control or mappability track is supplied.", "GC, TSS and IgG matching cannot fully remove accessibility or repetitive alignment biases.",
                               "Only one IgG library per state is available; biological variability of the control is not estimated.",
                               "CPM contains observed coverage after existing BAM filters; no spike-in and no absolute occupancy inference.", "Existing IP and IgG duplicate handling may differ; this job does not test that effect.",
                               "PPOIS mean/max are descriptive per-base scores, not regional P/Q values.", "Gene-level TSS is not a comprehensive transcript TSS annotation."]}
        (args.out / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        (args.out / "README.txt").write_text(
            "Exact summaries use original bigWig intervals, not zoom means.\n"
            "CPM mean_full_interval=sum(score*length)/full region length: gaps on a valid contig mean zero observed filtered coverage. Missing contig and explicit nonfinite entries stay NA.\n"
            "PPOIS mean_full_interval is NA if any bases lack encoded values; mean_encoded_bases and encoded_fraction are also reported. PPOIS means are not region P-values.\n"
            "Both-above-IgG, differences and ratios are descriptive, NOT statistical calls. A zero IgG denominator yields NA for unregularized ratio; pseudocount .01 and .1 CPM results permit sensitivity checks.\n"
            "Matched comparators use chromosome, width, TSS proximity, GC, ACGT fraction, IgG and both Pol2 replicate CPM. They are not chosen by Tag signal. Regions can be truly Tag-bound.\n"
            "Source selection uses the current strict Tag/Pol2 peak sets. Percentiles are post-selection descriptive ranks, NOT P-values, Q-values or proof of specificity. Insufficient matching is explicitly recorded.\n"
            "137-Tag signal at the same loci is retained; sharing with 137 is not classified as artifact.\n"
            "No main pipeline, peak calls, old outputs or bigWigs are overwritten.\n", encoding="utf-8")
        (args.out / "COMPLETED.txt").write_text("Completed exact summaries and exploratory matched-region comparisons.\n", encoding="utf-8")
    finally:
        reader.close()
        if fasta is not None:
            fasta.close()


if __name__ == "__main__":
    main()
