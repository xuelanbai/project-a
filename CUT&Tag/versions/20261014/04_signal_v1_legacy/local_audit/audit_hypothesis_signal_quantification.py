"""Audit the received region summaries without treating conditional means as enrichment."""
import csv
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "work/hypothesis_signal_import_20261009/hypothesis_signal_quantification"
SCREEN = ROOT / "work/hypothesis_screen_import/hypothesis_screen_current_peaks"
OUT = ROOT / "outputs/hypothesis_signal_audit_20261009"


def table(path):
    with path.open(encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter="\t")
        names = [x.lstrip("#").strip("'") for x in next(reader)]
        result = {}
        for cells in reader:
            key = cells[0], int(cells[1]), int(cells[2])
            if key in result:
                raise ValueError(f"Duplicate interval in {path}: {key}")
            result[key] = {name: float(v) for name, v in zip(names[3:], cells[3:])}
    return result


def bed(path):
    return [line.rstrip().split("\t") for line in path.read_text(encoding="utf-8-sig").splitlines() if line and not line.startswith("#")]


def finite(v):
    return v if math.isfinite(v) else None


def describe(values):
    vals = [v for v in values if math.isfinite(v)]
    return {"n": len(values), "finite": len(vals), "missing": len(values) - len(vals), "median": statistics.median(vals) if vals else None,
            "regions_mean_ge_2": sum(v >= 2 for v in vals), "regions_mean_ge_5": sum(v >= 5 for v in vals)}


def overlap(a, b):
    return a[0] == b[0] and max(a[1], b[1]) < min(a[2], b[2])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dataset_info = {}
    joined = {}
    for name in ("rRNA", "5S_rRNA", "candidate_genes_plusminus2kb", "macs2.119_Pol2_replicated_peaks", "seacr.119_Pol2_replicated_peaks"):
        cpm = table(DATA / f"{name}.CPM.regional_means.tsv")
        pp = table(DATA / f"{name}.pooled_MACS2_ppois.regional_means.tsv")
        assert set(cpm) == set(pp), name
        joined[name] = {key: dict(vals, **pp[key]) for key, vals in cpm.items()}
        dataset_info[name] = {"regions": len(cpm), "columns": {col: describe([r[col] for r in rows.values()]) for rows in (cpm, pp) for col in next(iter(rows.values()))}}

    annotations = {(x[0], int(x[1]), int(x[2])): {"gene": x[3], "gene_id": x[4], "strand": x[5]} for x in bed(SCREEN / "annotations/Ensembl115.genes.bed")}
    candidate_annotations = {(x[0], int(x[1]), int(x[2])): {"gene": x[3], "gene_id": x[4], "strand": x[5]} for x in bed(SCREEN / "annotations/candidate_genes.plusminus2kb.bed")}
    rrna = []
    for key, vals in joined["rRNA"].items():
        rrna.append(dict(chrom=key[0], start=key[1], end=key[2], **annotations.get(key, {}), **{col: finite(v) for col, v in vals.items()}))
    candidates = []
    for key, vals in joined["candidate_genes_plusminus2kb"].items():
        candidates.append(dict(chrom=key[0], start=key[1], end=key[2], length_bp=key[2]-key[1], **candidate_annotations[key], **{col: finite(v) for col, v in vals.items()}))

    overlap_stats = {}
    selected_polii = {}
    for caller in ("macs2", "seacr"):
        polii = joined[f"{caller}.119_Pol2_replicated_peaks"]
        tag = [(x[0], int(x[1]), int(x[2])) for x in bed(SCREEN / f"peak_sets/{caller}_119_Tag_rep1-rep2_shared.bed")]
        other_tag = [(x[0], int(x[1]), int(x[2])) for x in bed(SCREEN / f"peak_sets/{caller}_137_Tag_rep1-rep2_shared.bed")]
        # Raw interval counts, not enrichment estimates or enrichment p-values.
        tag_polii = [x for x in tag if any(overlap(x, p) for p in polii)]
        overlap_stats[caller] = {"119_Tag_shared_total": len(tag), "119_Tag_shared_overlaps_119_Pol2": len(tag_polii),
                                 "119_Tag_shared_overlaps_137_Tag": sum(any(overlap(x, p) for p in other_tag) for x in tag)}
        selected = [dict(chrom=p[0], start=p[1], end=p[2], **{col: finite(v) for col, v in vals.items()}) for p, vals in polii.items() if any(overlap(p, x) for x in tag)]
        selected_polii[caller] = {"count": len(selected), "ppois_descriptive_only": {col: describe([row[col] if row[col] is not None else math.nan for row in selected]) for col in ("119_Tag_pooled_vs_IgG", "119_Pol2_pooled_vs_IgG")},
                                  "selected_intervals": selected if caller == "macs2" else []}

    # Annotate each strict MACS2 Tag interval that overlaps a strict Pol II interval.
    # The signal means are over the Pol II interval, NOT over the Tag interval.
    macs_tag = [(x[0], int(x[1]), int(x[2])) for x in bed(SCREEN / "peak_sets/macs2_119_Tag_rep1-rep2_shared.bed")]
    nearest = {}
    for x in bed(SCREEN / "nearest_genes/macs2_119_Tag.nearest_gene.tsv"):
        key = (x[0], int(x[1]), int(x[2]))
        nearest.setdefault(key, []).append({"gene": x[6], "gene_id": x[7], "distance_to_gene_body_bp": int(x[9])})
    strict_sites = []
    for t in macs_tag:
        for p, vals in joined["macs2.119_Pol2_replicated_peaks"].items():
            if overlap(t, p):
                overlap_bp = min(t[2], p[2]) - max(t[1], p[1])
                strict_sites.append({"Tag_interval": list(t), "Pol2_interval": list(p), "overlap_bp": overlap_bp,
                                     "Tag_fraction_overlap": overlap_bp/(t[2]-t[1]), "nearest_gene": nearest.get(t, []),
                                     "Tag_ppois_mean_over_Pol2_interval": vals["119_Tag_pooled_vs_IgG"],
                                     "Pol2_ppois_mean_over_Pol2_interval": vals["119_Pol2_pooled_vs_IgG"]})

    original_rrna = {(x[0], int(x[1]), int(x[2])) for x in bed(SCREEN / "annotations/Ensembl115.rRNA.bed")}
    trna = bed(SCREEN / "annotations/Ensembl115.tRNA.bed")
    result = {"input_folder": str(DATA), "scope": "Audit of received summaries; not a corrected signal reanalysis", "datasets": dataset_info, "candidate_genes": candidates,
              "rRNA_regions": rrna, "omitted_rRNA_regions": [list(x) for x in sorted(original_rrna - set(joined["rRNA"]))],
              "tRNA_annotation": {"total": len(trna), "chromosomes": sorted(set(x[0] for x in trna)), "nuclear_count": sum(x[0] not in ("MT", "chrM") for x in trna)},
              "previous_shared_peak_overlap": overlap_stats, "119_Pol2_shared_regions_overlapping_Tag_shared_peaks": selected_polii,
              "MACS2_strict_cooccupancy_sites": strict_sites,
              "limitations": ["multiBigwigSummary 3.5.1 uses bigWig.stats() default; uncovered bases are excluded from the mean denominator and zoom summaries may be used.",
                              "Candidate intervals span the whole gene plus 2 kb on both sides, not a TSS-centered window.",
                              "nan is missing signal summary; do not replace a missing chromosome with zero.",
                              "Region-average -log10(P) is not a region-level P-value, Q-value, or peak-calling decision.",
                              "No genome-wide pathway enrichment, rDNA-aware remapping, or protein interaction experiment is included."]}
    (OUT / "audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")

    lines = ["# 现有信号定量审计（2026-10-09）", "", "## 结论", "", "这次文件完成了区域信号汇总，但不能作为Tag/IgG富集的最终判定。原因是旧脚本的CPM均值未把同一有效染色体上的未记录碱基统一补零；候选区间又是整段基因±2 kb，而非TSS。应先修正定量，不必重跑比对或call peak。", "", "已有严格重复共同peak计数仍然有效：119-Tag/MACS2 22个，其中7个与119-Pol2共同peak重叠；SEACR 876个，其中845个重叠。两者都属于共占位线索，不能证明蛋白互作；比例不能跨caller直接比较，也没有随机背景检验。", "", "137与rRNA：AEMK02000489.1:43231–45101值得重点复核。现有数值显示137两个重复的CPM条件均值均高于119，但MACS2共同peak只触及该rRNA起始边界31 bp，不能说整个rDNA单位富集。另一条SEACR共同peak几乎覆盖的5.8S rRNA（AEMK02000695.1:10195–10348）在旧均值下137两个重复均低于IgG，应避免将该重叠直接当成验证成功。", "", "此前tRNA无重叠不能用来否定Pol III：输入tRNA注释只有17个MT条目，没有核tRNA。输入23个rRNA条目中MT两条也未出现在这次21行的定量表。", "", "## 数据完整性", "", "|区域集|行数|备注|", "|---|---:|---|"]
    for name, info in dataset_info.items():
        lines.append(f"|{name}|{info['regions']:,}|CPM和ppois坐标一致，已按坐标匹配，不依赖行顺序|")
    lines += ["", "## 候选基因的旧整段区域ppois均值", "", "以下仅为原文件的描述性值。均值低不代表启动子没有局部peak；也不是区域显著性P值。", "", "|候选基因|区域长度bp|119-Tag|137-Tag|119-Pol2|137-Pol2|", "|---|---:|---:|---:|---:|---:|"]
    for row in sorted(candidates, key=lambda x: x["gene"]):
        lines.append(f"|{row['gene']}|{row['length_bp']:,}|" + "|".join(f"{row[col]:.3f}" for col in ("119_Tag_pooled_vs_IgG", "137_Tag_pooled_vs_IgG", "119_Pol2_pooled_vs_IgG", "137_Pol2_pooled_vs_IgG")) + "|")
    lines += ["", "## 原始CPM条件均值示例（不能当正式fold enrichment）", "", "|rRNA坐标（BED零基半开）|119-Tag rep1/rep2|119-IgG|137-Tag rep1/rep2|137-IgG|137 pooled ppois均值|", "|---|---:|---:|---:|---:|---:|"]
    for row in rrna:
        if row["chrom"] not in ("AEMK02000489.1", "AEMK02000695.1"):
            continue
        fmt = lambda v: "NA" if v is None else f"{v:.3f}"
        lines.append(f"|{row['chrom']}:{row['start']}–{row['end']}|{fmt(row['119_Tag_rep1'])}/{fmt(row['119_Tag_rep2'])}|{fmt(row['119_IgG'])}|{fmt(row['137_Tag_rep1'])}/{fmt(row['137_Tag_rep2'])}|{fmt(row['137_IgG'])}|{fmt(row['137_Tag_pooled_vs_IgG'])}|")
    lines += ["", "## 119–Pol II的7个严格MACS2共占位候选", "", "Tag区间为两个Tag重复的共同区间，Pol II区间也要求两个重复共同出现。最近基因是注释线索，不等于受调控靶基因。ppois均值在对应Pol II区间上汇总，不能解释成共同区间的P值，也不是新的显著性检验。", "", "|Tag共同区间（BED）|最近基因|相交bp|占Tag区间|Tag ppois均值|Pol II ppois均值|", "|---|---|---:|---:|---:|---:|"]
    for row in strict_sites:
        t = row["Tag_interval"]
        lines.append(f"|{t[0]}:{t[1]}–{t[2]}|{','.join(g['gene'] for g in row['nearest_gene'])}|{row['overlap_bp']}|{row['Tag_fraction_overlap']:.1%}|{row['Tag_ppois_mean_over_Pol2_interval']:.3f}|{row['Pol2_ppois_mean_over_Pol2_interval']:.3f}|")
    lines += ["", "全体13,826个119-Pol II/MACS2共同区间上，119-Tag pooled ppois区域均值的中位数为0.484，Pol II为8.994；上述7个共占位候选对应中位数分别为2.744和27.588。这支持优先逐位点复核，不能推广为119在所有Pol II位点普遍富集。区域是按Pol II峰选择的，这个比较本身存在选择偏差，不能直接给两种蛋白比较显著性。", "", "SEACR的845/876是以Tag区间计数，798是以相交Pol II区间计数；一对多重叠可导致二者不同。这798个Pol II区间上的119-Tag ppois均值中位数为0.339，故SEACR高重叠率不等于已通过MACS2 IgG-relative信号验证。"]
    lines += ["", "## 下一步", "", "1. 对已有14条bigWig做原始区间精确积分：CPM在有效染色体内的未记录区间补零；缺失染色体标NA，不猜测chr前缀。", "2. 对12个候选基因分开输出基因体、gene-level TSS±2 kb、TSS±500 bp和±2 kb的100 bp分箱。注意gene-level TSS不是所有转录本TSS。", "3. 原有Tag共同peak和119-Pol2共同peak也做同口径定量；复核rRNA邻域分箱，区分内部与边缘信号。", "4. 跨状态CPM比值只是相对信号比较，没有spike-in，不能推断绝对全局占位变化。Inflammation需预定义基因集、匹配背景与正式富集检验，现有12个PPI候选不足以证明。", "5. 重复支持规则不变；无需通过放松MACS2阈值寻找预期结果。", "", "## 方法依据", "", "- [deepTools 3.5.1 的实际汇总代码](https://deeptools.readthedocs.io/en/3.5.1/_modules/deeptools/getScorePerBigWigBin.html)：调用bwh.stats，未设置exact=True。", "- [pyBigWig官方文档](https://github.com/deeptools/pyBigWig#compute-summary-information-on-a-range)：均值不把缺失碱基计入分母；原始intervals或exact统计可避免zoom近似。", "- 本地nf-core/cutandrun 3.2.2 conf/modules.config 的CPM bamCoverage配置含--skipNAs，生成稀疏轨道。", "", "这份审计没有重新计算任何bigWig、FRiP/TSS score或蛋白互作显著性。"]
    (OUT / "结果解读与定量修正说明.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"datasets": dataset_info, "overlap": overlap_stats, "selected_polii": {k: {c: v for c, v in d.items() if c != 'selected_intervals'} for k, d in selected_polii.items()}, "report": str(OUT / "结果解读与定量修正说明.md")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
