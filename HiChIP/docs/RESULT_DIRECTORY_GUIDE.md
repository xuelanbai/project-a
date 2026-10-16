# HiChIP 项目五个目录的数据流程与文件索引

更新：2026-10-15。本文依据本地保存的 `README`、`config.json`、`run_metadata.json`、QC 表和实际文件清单整理。**本地文件夹是结果副本或作图包，不等同于服务器完整运行目录；服务器文件是否仍在线，需要登录服务器再确认。** 本文不重新计算数据，也不修改既有结果。

## 先看哪一版

| 想找的内容 | 本地入口 | 关键版本 | 使用时注意 |
| --- | --- | --- | --- |
| Promoter / enhancer 的坐标定义与 BED | `01.启动子、增强子区域划分` | 29 的定义敏感性导出、34 的 P1 启动子状态 | 29/34 的 P1 判据是 H3K27ac peak **与** TSS±2 kb **重叠**；38 才增加“完整落入对应 gene 的 TSS 窗口”。 |
| 增强子 H3K27ac 变化与 ATAC-core motif | `02.最终-Strength-Weaken_Enhancer_motif` | 27 `run_41129`；31 是选出的 AME 1.5 倍/500 bp 作图包 | 27 同时提供 1.2/1.5 倍及 200/500 bp；不要将 31 当作独立的原始 BAM 全流程。 |
| P1 启动子 motif、各种 mark 的 0/6 h 变化 | `02.最终-Strength-Weaken_promoter_motif` | 35 `run_42671`；36 `run_43021` | 35 的启动子状态来自 34 的 1.5 倍版本。36 文件夹名含“p1.2”，但代码说明和来源状态仍为 **1.5 倍**。 |
| motif-bearing enhancer 通过 HiChIP 连到基因 | `03.通过loop把增强减弱活性的Enhancer连接到Promoter和gene` | 32 `run_42014`；34 v2 为全部显著 motif 家族整理包 | 32 比较旧版三种 promoter 定义；34 v1 仅选代表性 hub，v2 才保留全部 21 个显著家族。均不是 38 的严格“完整落入 + 每基因一个 primary promoter”结果。 |
| 最新严格启动子、方向分组、mark 和 loop 联合表 | `38_primary_promoter_direction_only_motif_CRE_table_v1` | 38 `run_43558` | 方向阈值 1.0 不是“至少变化一倍”的效应量阈值；Weaken motif **未完成检验**，见下文。 |

## 服务器与本地位置

服务器结果共同前缀，以下记作 `R`：

```text
R=/mnt/nas/swine/home/jxwu/Projects/A/07.HiChIP/NFkB_target_loci_0h_6h_v2/06_chromatin_state_v1
```

服务器代码共同前缀，以下记作 `S`：

```text
S=/mnt/nas/swine/home/jxwu/Projects/A/scripts_new
```

本地五个目录均位于 `D:\Lucas临时文件\07_HiChIP结果\`。下表的 `R/…`、`S/…` 是上述完整服务器前缀加相对路径，并非本地 D 盘路径。服务器的运行编号保留在结果目录名中；本地部分子文件夹被手工改名，查找服务器文件应以 `code/config.json` 或 `run_metadata.json` 为准。

| 版本 | 服务器结果位置 | 服务器代码位置 |
| --- | --- | --- |
| 27 增强子与 motif | `R/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1/run_41129` | `S/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1` |
| 29 坐标定义敏感性 | `R/29_CRE_definition_sensitivity_export_v1/export_41527` | `S/29_CRE_definition_sensitivity_export_v1` |
| 32 motif–HiChIP–gene | `R/32_final_AME_1p5_500bp_HiChIP_targets_draft_v1/run_42014` | 该本地副本未核实对应的服务器脚本目录；以结果元数据与归档代码为准。 |
| 34 P1 启动子状态 | `R/34_P1_promoter_strength_weaken_v1/run_42216` | `S/34_P1_promoter_strength_weaken_v1` |
| 35 P1 启动子 motif | `R/35_P1_promoter_ATAC_core_HOMER_AME_v1/run_42671` | `S/35_P1_promoter_ATAC_core_HOMER_AME_v1` |
| 36 区域/基因多 mark 信号 | `R/36_linked_CRE_promoter_mark_dynamics_v1/run_43021` | `S/36_linked_CRE_promoter_mark_dynamics_v1` |
| 38 最新联合表 | `R/38_primary_promoter_direction_only_motif_CRE_table_v1/run_43558` | `S/38_primary_promoter_direction_only_motif_CRE_table_v1` |

`31_AME_1p5_500bp_final_candidate_v1`、`33_promoter_definition_result_figures_v1`、两个 `34_unbiased_enhancer_motif_gene_pathway` 文件夹是本地筛选/作图整理包；本次在这些副本里**未核实独立的服务器输出路径**。其中 34 v1/v2 的源数据来自服务器的 32 `run_42014`，不是服务器的 34 P1 启动子状态作业。两个名称相同的“34”属于不同分析。

## 整体数据流

```text
原始/已处理 BAM、H3K27ac/H3K4me3 peak BED、GTF
        ├─ 27：ATAC BAM 重新 call peak/summit → H3K27ac 共识增强子 → 两时间点计数 → 1.2/1.5 倍状态 → ATAC-core AME/HOMER
        │    ├─ 29：比较 promoter/enhancer 坐标定义，核验 enhancer 与 27 一致
        │    ├─ 31：复制 27 的 1.5 倍、500 bp AME 候选结果用于作图
        │    └─ 32：AME motif 位点 + 0/6 h HiChIP loop → promoter/gene → 33 与 34 v1/v2 作图包
        └─ 34 P1：P1 启动子共识坐标 → 同坐标 H3K27ac/ATAC 计数 → 1.5 倍状态
             ├─ 35：复用 27 的 ATAC core → P1 启动子 AME/HOMER
             └─ 36：对增强子/P1 启动子/TSS 窗口重新统计五种 mark → loop 联合矩阵
                  └─ 38：严格完整落入 TSS 窗口 → 每 gene 选一个 primary promoter
                           → 方向分组与 1.2 倍敏感性表 → 0/6 h loop 联合表 → primary-promoter motif
```

**复刻依赖：**38 不从 BAM 重新计数，而是读取 27 的 ATAC core、34 的 P1 候选及状态矩阵、36 的全区域 BAM 计数矩阵，再读取原始 HiChIP BEDPE。故只复制 38 的 `code/` 与 `run_43558/`，不能从原始 BAM 独立复刻全部流程；需保留上游版本、输入文件、参考基因组和 motif 库。27 的本次运行也复用了 26 的 ATAC peak 调用结果，`code/README.md` 记录了从 BAM 重新调用 MACS3 peak 的复刻方法。

## 1. `01.启动子、增强子区域划分`

### 29：区域定义比较与 BED 导出

本地：`Promoter、Enhancer位置划分方案/29_CRE_definition_sensitivity_export_v1/export_41527/`。输入是两个时间点 H3K27ac/H3K4me3 peak BED、27 产生的 ATAC peak、Ensembl 115 GTF；**不读取 BAM，不做 motif**。比较 TSS±500、±1500、±2000 bp，分别输出 0 h、6 h、跨时间共识 BED。共识只合并实际重叠的 peak，不合并仅首尾相接的区间。

- Promoter 的 P1：完整 H3K27ac peak 与某 gene 的 TSS 窗口至少重叠 1 bp，同时间 ATAC/H3K4me3 各与完整 H3K27ac peak 重叠至少 1 bp；输出仍是**完整 H3K27ac peak**。P2 要求支持信号位于 peak 与该 gene TSS 窗口的交集，仍输出完整 peak；P3 裁剪到交集；P4 只保留四类区间共同覆盖的碱基。
- Enhancer：完整 H3K27ac peak 与同时间 ATAC peak 至少重叠 1 bp；若与 TSS 窗口或设定的 H3K4me3 集合至少重叠 1 bp，则排除整个 H3K27ac peak。TSS±2 kb、`either_time_v27` 模式排除任一时间 H3K4me3，与 27 完全一致。`05_v27_reproduction_QC.tsv` 显示原始 peak membership 48,607、共识 CRE 28,405，均 `PASS`。
- 本地 `BED/promoter/TSS_pm2000bp-最终用这版/P1_full_H3_peak_support-最终用这版/`：`promoter_0h.bed`、`promoter_6h.bed`、`promoter_consensus_0h_6h.bed`。P1 共识共 10,010 个区间；逐时间点为 9,782、9,520 个来源 peak。
- 本地 `BED/enhancer/TSS_pm2000bp-最终用这版/H3K4me3_either_time_v27-最终用这版/`：对应三个时间层级的 enhancer BED；服务器原始输出名称不带本地人工添加的“最终用这版”后缀，精确相对路径可用 `00_BED_manifest.tsv` 核对。
- `01_promoter_stepwise_counts.tsv`、`02_promoter_definition_counts.tsv`、`02_promoter_consensus_counts.tsv`、`03_enhancer_stepwise_counts.tsv`、`04_enhancer_definition_counts.tsv`、`04_enhancer_consensus_counts.tsv`：各定义的逐步保留量、gene/peak/共识区间数及长度；`说明_定义与输出列.txt` 解释列意义。

29 的 P1 是**重叠 TSS 窗口**规则，不是 38 的完整落入规则。不要把 `promoter_consensus_0h_6h.bed` 直接当作“每 gene 唯一的严格启动子”。

### 34：P1 启动子的旧版 1.5 倍状态

本地：`Promoter_strength-weaken划分/34_P1_promoter_strength_weaken_v1/run_42216/`。按上述 P1 建立 0/6 h 完整 H3K27ac peak 和共识坐标，在同一坐标上对 H3K27ac、ATAC 合并 BAM 计数。按各 BAM 主染色体 proper-pair 总数换算 tags/10M，不除以区域长度、不做中位数中心化。H3K27ac 变化 ≥1.5 倍且 6 h 有 P1 支持、6 h H3K27ac 原始 tags ≥5 为 Strength；反向为 Weaken；其余 Stable。**ATAC 原始 tags ≥5 不作为本作业状态过滤**，ATAC 方向仅作注释。结果：Strength 168、Weaken 527、Stable 9,315 个 P1 共识区域。

本地的 P1 文件夹被改名为 `P1-最终p1版本/`，服务器代码说明中仍称 `P1/`。关键文件：`02_P1_promoters_{0h,6h}.bed`、`03_P1_promoter_consensus_0h_6h.bed`、对应的 `*_gene_associations.tsv`；`04_P1_promoter_signal_and_state_matrix.tsv` 是两时间点计数与状态总表；`states/05_P1_{strength,weaken,stable}_promoters.bed` 是状态 BED；`05_P1_promoter_state_counts.tsv` 和 `run_metadata.json` 用于核对。

## 2. `02.最终-Strength-Weaken_Enhancer_motif`

### 27：增强子与 motif 的主体结果

本地：`27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1/run_41129/`。增强子坐标来自 **完整 H3K27ac peak**：同时间 ATAC 至少 1 bp 支持，排除 TSS±2 kb 和任一时间 H3K4me3 重叠，再合并 0/6 h 真正重叠的 H3K27ac peak，得到 28,405 个共识 CRE。0/6 h 的 H3K27ac 与 ATAC BAM 在相同 CRE 坐标计数，tags/10M，不按长度除，也不中位数中心化。活跃端要求 H3K27ac、ATAC 原始 tags 均 ≥5；ATAC 自身的增减**不是** enhancer 状态判据。H3K27ac 变化分 1.2、1.5 倍两档；1.5 倍下 Strength 3,682、Weaken 1,764 个 CRE。

ATAC peak 从 0/6 h BAM 经 MACS3 `BAMPE, q=0.01, --call-summits` 调用；跨时间真正重叠的 ATAC peak 构成 consensus peak。每个 consensus peak 从两个时间点的原始 summit 中选择归一化 ATAC 信号更高时间点的一个代表 summit。落在 H3K27ac CRE 内的 ATAC core 可多于一个。围绕 summit 取 200 或 500 bp 序列，GC 匹配背景，分别做 AME（JASPAR vertebrate）和 HOMER（vertebrate known motif，未做 de novo）。

关键文件：`02_ATAC_consensus_cores.tsv`（代表 summit 及选择依据）、`03_H3K27ac_consensus_CRE.bed`（全部增强子坐标）、`03_CRE_signal_matrix.tsv`（0/6 h 信号和状态依据）、`04_ATAC_core_to_CRE.tsv`（core 与 CRE 的一对多关系）、`fold_1p5/regions/{strength,weaken}_CRE.bed`（最终 1.5 倍状态区间）、`fold_1p5/regions/*_ATAC_core_500bp.bed`（motif 窗口）、`05_AME_all.tsv`、`05_HOMER_all.tsv` 及代表家族表。图以 `Figure_AME_*`、`Figure_HOMER_*` 命名，q 排序与显著结果的富集倍数排序均有 PNG/PDF。`motif_status.json` 比 `05_motif_sampling.tsv` 中较早写入的 `planned` 字段更适合判断最终完成状态。

### 31：选定 AME 1.5 倍、500 bp 候选作图包

本地：`汇总结果文件/31_AME_1p5_500bp_final_candidate_v1/`。只复制 27 的 `strength_1p5_500bp_ATAC_core` 和 `weaken_1p5_500bp_ATAC_core`。`01_AME_tables/` 是完整及代表家族表；`02_raw_AME/` 是 AME 原生输出；`03_motif_inputs/` 有前景/背景 FASTA 与 GC 配对；`04_regions/` 有 CRE/core BED 和映射；`05_figures_current/` 是当时图；`06_provenance/` 保存源运行元数据、哈希和代码快照。`汇总结果图/31_AME_1p5_500bp_figure_layout_draft_v7_compact_equal_spacing/` 是后续排版图，不是新的一轮 motif 统计。

## 3. `02.最终-Strength-Weaken_promoter_motif`

### 35：旧 P1 启动子 motif

本地：`35_P1_promoter_ATAC_core_HOMER_AME_v1/run_42671/`。直接读取 34 的 P1 **1.5 倍**启动子状态和 27 的 ATAC consensus core，不重新定义启动子或重新 call peak。围绕 core summit 扫描 500 bp，非前景 P1 core 中选 GC 匹配背景；AME 与 HOMER 均运行，两个方向的 `motif_status.json` 都为 `completed`。输入前景 168 个 Strength、527 个 Weaken P1；具有 core 的为 152、385 个；有效 core 窗口为 190、392 个。`01_P1_state_and_ATAC_filter_QC.tsv` 核对 ATAC ≥5；该批数据中此过滤未再删除 promoter。

关键文件：`02_ATAC_core_to_P1_promoter.tsv`、`regions/`（状态 BED/core 与窗口）、`03_motif_sampling.tsv`（GC 匹配）、`04_AME_all.tsv`、`04_HOMER_all.tsv`、对应的 `*_representative_families.tsv` 和显著家族表、`Figure_*_q_rank` 与 `Figure_*_enrichment_rank` PNG/PDF。

### 36：增强子、P1、TSS 的五种 mark 动态

本地：`36_linked_CRE_promoter_mark_dynamics_v1-27ac显著p1.2/run_43021/`。虽然**本地目录名含 p1.2**，作业说明明确使用 27/34 原有 **1.5 倍**状态；结果目录的命名不能用来推断阈值。对 27 增强子、34 P1 启动子和 gene TSS±2 kb 三套坐标，分别在 0/6 h 统计 H3K27ac、H3K4me3、Pol II、ATAC、H3K27me3 合并 BAM 的原始计数和 tags/10M；ChIP 类 mark 还提供相对 IgG 的次级指标。利用 32 的 loop–gene 来源表关联 enhancer–promoter–gene。比例按唯一 promoter 计，mark 一致性按唯一 enhancer–promoter pair 计，避免 loop 行数重复加权。

关键文件：`01_enhancer_mark_signal_matrix.tsv`、`02_P1_promoter_mark_signal_matrix.tsv`、`03_gene_TSS_pm2kb_mark_signal_matrix.tsv`、`03_source_vs_recount_QC.tsv`、`04_region_mark_long.tsv`、`04_region_mark_summary.tsv`、`05_HiChIP_enhancer_promoter_gene_mark_matrix.tsv`、`06_linked_promoter_state_distribution.tsv`、`07_linked_enhancer_promoter_mark_concordance.tsv`、`08_linked_promoter_mark_summary.tsv` 和两张 PNG/PDF 图。**36 仍使用旧 P1 重叠定义及 1.5 倍分类；38 从这些已计数区域里重新筛选/标状态。**

## 4. `03.通过loop把增强减弱活性的Enhancer连接到Promoter和gene`

### 32：主 HiChIP 连接与启动子定义比较

本地：`32_final_AME_1p5_500bp_HiChIP_targets_draft_v1/run_42014/`；旁边 `run_42013/` 是更早一次运行，不应与 42014 混合。输入 27 的 1.5 倍/500 bp AME 显著且富集的代表 motif 家族：Strength 16、Weaken 5。FIMO 在 motif 窗口中定位位点，阈值 `P≤1e-4`。完整 enhancer 与 HiChIP 一端至少重叠 1 bp，另一端与 promoter 重叠至少 1 bp；Strength 使用 6 h loop、Weaken 使用 0 h loop，FitHiChIP 文件阈值 q≤0.01。分别比较裸 TSS±2 kb、P1、P2 三种 promoter 定义；“motif 位点本身覆盖 CRE 端 anchor”是另外保留的严格证据标记。

关键文件：`01_selected_AME_motifs.tsv`（输入 motif）、`02_promoter_definition_counts.tsv`、`03_target_definition_counts.tsv`、`04_bridge_matrix_summary.tsv`、`promoters/`（三种定义的时间点 BED 与 gene 映射）、`targets/`（按定义和 Strength/Weaken 分层的连接、候选基因及功能结果）。本地部分子目录人为加了“最终用这版”后缀；服务器用原始名称 `P1_full_H3_peak_support`。

### 33、34 v1、34 v2：下游展示包

- `1版_先挑选一些脂代谢相关gene/33_promoter_definition_result_figures_v1/`：比较三种 promoter 方案，并展示预先选定脂代谢相关基因的矩阵。用于讨论排版和定义，不是无偏全基因结果。
- `2版_看各显著富集motif通过loop关联到的gene/34_unbiased_enhancer_motif_gene_pathway_v1/`：先按 AME 富集值选每方向前 5 个家族，再选代表 enhancer hub；保留 6 个 hub、10 个 motif 家族。目录中的 `01`–`07` 是 hub、motif、靶基因、loop、FIMO 位点、通路检验及显著通路；`08_source_file_manifest.tsv` 记源文件和 SHA-256。**这不是全部显著家族。**
- 同级的 `34_unbiased_enhancer_motif_gene_pathway_v2_all_motifs/`：修正版保留全部 **21 个**显著富集家族；主结果为 `01_all_enriched_motif_families.tsv`、`02_all_strict_enhancer_motif_summary.tsv`、`03_all_strict_enhancer_motif_target_genes.tsv`、`04_all_strict_enhancer_motif_loop_links.tsv`、`05_all_linked_motif_sites.tsv`、`06_motif_family_strict_target_summary.tsv`、`07_motif_family_pathways_all_tested.tsv`、`08_motif_family_pathways_significant.tsv`。`09_optional_representative_enhancer_per_family.tsv` 仅用于紧凑绘图，**不用于定义主结果**；`10_source_file_manifest.tsv` 记溯源。严格关联要求 motif 位点本身至少 1 bp 覆盖 CRE 端 loop anchor。通路结果是基于该版本靶基因集合的富集检验，不等于基因表达变化。

## 5. `38_primary_promoter_direction_only_motif_CRE_table_v1`

本地：`run_43558/`。这是**最近的严格启动子与联合表版本**，但不应误称“全部 motif 已完成的最终版”。

1. 从 34 的 10,010 个 P1 共识区间及其 promoter–gene 关联出发，只保留**完整 H3K27ac 共识区间完全落入该 gene TSS±2 kb** 的关联。原有 12,264 条 promoter–gene 关联中，4,195 条未完全落入而排除，剩 8,069 条严格候选关联。
2. 在每个 gene 内选**至多一个** primary promoter：先看是否覆盖 TSS，再按距离 TSS 最近，仍并列时用双时间点平均 H3K4me3 密度、ATAC 密度，最后用 region ID。结果 7,228 个 gene–primary-promoter 关联，对应 **6,244 个唯一 promoter 区间**。一个区间仍可服务多个相邻 gene；这是“每 gene 一个”，**不是“一 promoter 只连一 gene”**。完整候选排序保存在 `02_strict_promoter_candidate_gene_ranking.tsv`。
3. 直接复用 36 的全区域 BAM 计数矩阵；不重算 BAM。按 H3K27ac 6 h/0 h 方向分组：比值 `>1` 为 Strength、`<1` 为 Weaken；要求活跃时间点有对应 peak 支持，H3K27ac 与 ATAC 原始 tags 均 ≥5。`1.0` 仅表示正负方向，**不是有意义变化幅度或统计显著性阈值**。另输出 1.2 倍敏感性状态列。增强子坐标仍是 27 的 28,405 个 CRE，但 38 对它们也重新按方向分类，所以 38 的 Strength/Weaken 数量**不能直接与 27 的 1.5 倍 motif 前景比较**。
4. 两时间点 HiChIP BEDPE 均参与连接，保留 loop 时间和是否匹配活跃端的标记，不只用 32 的“Strength→6 h、Weaken→0 h”规则。`04_enhancer_promoter_gene_loop_mark.tsv.gz` 含 54,577 行连接记录；不是 54,577 个独立基因或独立 CRE。
5. Primary-promoter motif 沿用 27 的 ATAC summit ±250 bp（500 bp），对方向分组做 AME/HOMER。然而 **Strength 完成，Weaken 因有效 GC 配对比例仅 71.1%（低于设置的 90%）而标记 `not_tested_low_matched_fraction`；不能将 Weaken 的空表或空白图解释为“没有显著 motif”。**即使子目录有 `completed.json`，也要看 `motif_status.json` 和 `03_motif_sampling.tsv`。

结果文件：`01_strict_promoter_definition_QC.tsv`（完整落入及每 gene 选一的计数）、`01_state_counts.tsv`（方向版与 1.2 倍敏感性计数）、`02_primary_promoter_gene_state_mark.tsv`（每 gene 一条）、`02_primary_promoter_region_state_mark.tsv`（每唯一 promoter 一条，故行数不同）、`03_enhancer_region_state_mark.tsv`（全部增强子及 mark）、`04_HiChIP_mapping_QC.tsv`、上述压缩 loop 联合表；`BED/` 有 enhancer/primary-promoter 的全量与方向状态 BED4，列说明见 `BED/README.txt`；`05_primary_promoter_motif_direction_only/` 有已测试的 motif 表、GC 配对、图、状态和代码。最新方向版计数：增强子 Strength 16,186、Weaken 9,853、Stable 2,366；唯一 primary promoter Strength 2,441、Weaken 3,522、Stable 281。1.2 倍计数另见 `01_state_counts.tsv`。

## 输入文件、文件格式与复刻检查

- 上游原始分析文件主要在服务器 `/mnt/nas/swine/home/jxwu/Projects/A/02.bam/`（0/6 h 各 assay 合并、去重复、排序 BAM及 `.bai`）、`/mnt/nas/swine/home/jxwu/Projects/A/06.peak/`（H3K27ac、H3K4me3 等 peak BED），以及 `R/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1/run_41129/peaks/`（27 的 ATAC MACS3 peak）。原始 0/6 h HiChIP loop 位于 `/mnt/nas/swine/home/jxwu/Projects/A/07.HiChIP/`，分别是 `AH-1-7-13.FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bedpe` 和 `AH-2-8-14.FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bedpe`。参考 FASTA/GTF 在 `/mnt/nas/swine/home/jxwu/Projects/references/pig/Sscrofa11.1/`；AME 使用服务器 JASPAR2026 vertebrate MEME 库。**精确输入与 SHA-256 以各作业的 `run_metadata.json` / `config.json` 为准。**
- `.bed` 无表头，0-based、half-open：第 1–3 列是染色体/起点/终点。第 4 列通常是稳定 region ID；**29 导出的 BED6** 另有第 5 列固定 0、第 6 列 `.`；**38 的 BED4** 格式见 `BED/README.txt`。BED 不含 gene 关系时，不要按最近 gene 自行推断，要与相应的 `*_gene_associations.tsv` 或联合 TSV 连接。
- `.tsv` 一般有表头；`.tsv.gz` 是压缩 TSV，可直接用于脚本读取。`*_matrix.tsv` 是区域/关联的宽表；`*_QC.tsv` 是筛选、比对或映射检查；`*_membership.tsv` 显示合并区间与来源 peak/core 的多对多关系；`*_all.tsv` 保留全部 motif 或通路结果；`*_representative_families.tsv` 是为避免相似 motif 重复展示而挑选的家族代表，不等于删除原始 motif 命中。
- 比较不同版本计数时，先核对**单位**：H3K27ac peak、共识 CRE、ATAC core、500 bp 窗口、promoter–gene 关联、唯一 promoter、loop 行、唯一 gene 是不同对象；还要核对阈值（1.5/1.2/仅方向）、TSS 关系（重叠/完整落入）、所用 loop 时间和是否需要 motif 位点直接覆盖 anchor。
- 以上 pooled/merged BAM 信号与 H3K27ac 状态是描述性 0/6 h 比较；没有生物学重复层面的差异检验 P 值或 FDR。HiChIP loop + motif 位点表示候选物理关联，不自动证明 TF 占位、因果调控或靶基因表达变化。
