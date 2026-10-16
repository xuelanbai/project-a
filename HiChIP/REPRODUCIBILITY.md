# 复刻顺序、输入与边界

## 原始输入

脚本的 `config.json` 给出实际服务器路径。主要依赖包括：

- `/mnt/nas/swine/home/jxwu/Projects/A/02.bam/`：0 h、6 h 的合并、去重复、排序并索引的 ATAC、H3K27ac、H3K4me3、Pol II、H3K27me3 和 IgG BAM。不同环节使用的 assay 不完全相同。
- `/mnt/nas/swine/home/jxwu/Projects/A/06.peak/`：H3K27ac、H3K4me3 等原始 peak BED。文件名带 `narrowPeak` 不意味着实际有第 10 列 summit；原 H3K27ac 文件实测只有三列。
- `/mnt/nas/swine/home/jxwu/Projects/A/07.HiChIP/`：0 h、6 h FitHiChIP BEDPE loop。
- `/mnt/nas/swine/home/jxwu/Projects/references/pig/Sscrofa11.1/`：基因组 FASTA、`.fai`、Ensembl 115 GTF；以及 JASPAR2026 vertebrate motif 库。

所有 BED 使用 0-based、half-open 坐标；TSV 有表头。输入 BAM 的主染色体 proper-pair 数用于 tags/10M 归一化。相同区域在两个时间点按同一坐标比较；没有按区域长度除，也没有全局中位数中心化。

## A. 从 BAM 到当前增强子 motif（27）

1. 在 `workflows/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1/config.json` 配置原始 BAM、H3K27ac/H3K4me3 peak、FASTA/GTF、JASPAR 库及输出根目录。
2. `check_env.sh` 只做轻量环境检查；在有 Slurm 的集群用 `submit.sh` 提交。正常的从原始输入复刻**不要设置 `REUSE_PEAKS_FROM`**，这样 0/6 h ATAC BAM 会重新由 MACS3 `BAMPE, q=0.01, --call-summits` 调 peak/summit。本次存档的 `run_41129` 曾复用上一轮已经调出的 ATAC peak；代码仍支持从 BAM 开始。
3. 同时间 H3K27ac peak 与 ATAC peak 至少 1 bp 重叠，排除 TSS±2 kb 与任一时间 H3K4me3 重叠的完整 H3K27ac peak；跨时间只合并真正重叠的 H3K27ac peak。用统一共识坐标统计两时间点 BAM。状态使用 1.2 和 1.5 倍阈值。
4. 将 0/6 h ATAC peak 建立共识，每个共识 peak 选择一个原始 MACS3 summit，映射到增强子；围绕 summit 扫描 200/500 bp。AME 使用 JASPAR vertebrate；HOMER 使用 vertebrate known motif，未做 de novo。
5. 核对 `run_metadata.json`、`motif_status.json`、`03_H3K27ac_consensus_CRE.bed`、`03_CRE_signal_matrix.tsv`、`04_ATAC_core_to_CRE.tsv` 与 `fold_1p5/regions/`。存档运行 `run_41129` 有 28,405 个共识增强子；1.5 倍下 Strength 3,682、Weaken 1,764 个 CRE。该数字是复核参照，不是新运行的硬编码目标。

## B. 从 P1 候选到最新严格 primary promoter（34-P1 → 36 → 38）

1. `34_P1_promoter_strength_weaken_v1` 使用完整 H3K27ac peak 与 gene TSS±2 kb **重叠**、同时间 ATAC/H3K4me3 peak 支持，建立跨时间共识 P1；在 BAM 上对 H3K27ac/ATAC 计数，按 1.5 倍分状态。它是 38 的**候选坐标及映射上游**，不是最新严格定义本身。
2. `36_linked_CRE_promoter_mark_dynamics_v1` 读取 27 的全部增强子、34-P1 的全部启动子、32 的 loop–gene 来源表，对 H3K27ac、H3K4me3、Pol II、ATAC、H3K27me3 和 IgG 的 BAM 在同一组区间重新计数。`01_enhancer_mark_signal_matrix.tsv` 和 `02_P1_promoter_mark_signal_matrix.tsv` 是后续 38 复用的**全区域 BAM 计数矩阵**。目录别名曾写“p1.2”，但存档方法仍是 1.5 倍状态；以代码/元数据而非别名为准。
3. `38_primary_promoter_direction_only_motif_CRE_table_v1` 只保留 **完整 P1 共识区间完全落入对应 gene TSS±2 kb** 的候选；每 gene 优先选覆盖 TSS 的候选，否则选最近者，并在距离相同处使用 H3K4me3/ATAC 密度和 region ID 决胜。保留全部候选排序表，可复查。之后复用 36 的计数、27 的 ATAC core、原始 0/6 h loop，不再读取 BAM 重新计数。
4. 38 的主分组是线性比值大于/小于 1 的**方向分组**，另有 1.2 倍敏感性状态。存档 `run_43558` 中 7,228 个 gene 对应 6,244 个唯一 primary promoter 区间；同一区间可以服务多个基因。Strength motif 完成；Weaken GC 匹配比例 71.1%，低于 90% 门槛，**未运行富集检验**。复刻时必须检查 `motif_status.json`，不能仅凭 `completed.json` 或生成了空图就认定成功。

## C. 其他版本如何使用

- `29_CRE_definition_sensitivity_export_v1`：只读 peak/GTF，输出多种定义的 BED 与统计，不是 38 的完整落入规则；其 TSS±2 kb、任一时间 H3K4me3 排除的 enhancer 模式应与 27 的 28,405 个坐标精确一致。
- `35_P1_promoter_ATAC_core_HOMER_AME_v1`：从 34 的旧 P1 1.5 倍状态和 27 的 ATAC core 生成 500 bp AME/HOMER，两方向均已完成。用于旧定义的对照，不等于 38 的严格 primary promoter motif。
- `32_final_AME_1p5_500bp_HiChIP_targets_draft_v1`：用 27 的 1.5 倍/500 bp AME motif、FIMO、时间匹配的 HiChIP loop 连接增强子与三种 promoter 定义及 gene。下游 `34_unbiased_enhancer_motif_gene_pathway_v1/build_all_results.py` 整理**全部 21 个显著家族**的 v2；`build_plot_data.py` 仅整理 v1 代表 hub。二者不能互相代替。
- `31_AME_1p5_500bp_figure_layout_draft_v1` 和 `33_promoter_definition_result_figures_v1` 是绘图/定义比较脚本。人工调整的 v7 增强子 AME 大图没有在此包内得到端到端脚本复刻验证，不能仅凭这些脚本承诺逐像素重现。

## 验证范围和未完成事项

- 本代码包只复制和整理源文件，没有在本机重跑大 BAM、MACS3、AME/HOMER、FitHiChIP 或通路分析；这些任务应在有相同输入、环境和 Slurm 的计算节点执行。
- 38 的 Weaken motif 仍未完成。若需要在**严格 primary promoter**定义下的 Strength/Weaken 成对 motif 图，必须先重新设计可匹配背景或调整有理据的分析方案，并重新执行；不能从 35 的旧 P1 Weaken 直接移植。
- 32/34-unbiased 的 motif–gene 结果使用旧 P1 promoter 定义。如果最终选择 38 的严格 primary promoter，增强子–promoter–gene 连接及功能统计应按新定义另行重算。
- 运行完成后应比对相同输入的 SHA-256、工具版本、参数、`run_metadata.json`、QC 表和 `motif_status.json`；源运行编号和参考数字只是审计锚点，不可替代实际验证。
