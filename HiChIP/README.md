# Project A：CRE、motif 与 HiChIP 分析代码

这是**单独整理的代码副本**，不包含 BAM、peak、BEDPE、FASTA、motif 数据库或分析结果。`workflows/` 中保留了原脚本的版本名称与文件内容，便于追溯和复刻；没有把几代分析代码拼成一个未经验证的“总脚本”。本仓库中的 `config.json` 记录原服务器路径，换环境运行前必须核对/修改。所有大数据和实际结果仍以服务器运行目录、`run_metadata.json`、`completed.json`、`motif_status.json` 为准。

## 当前采用哪些版本

| 分析对象 | 采用/最新版本 | 状态与解释 |
| --- | --- | --- |
| 增强子坐标、H3K27ac 状态及 ATAC-core motif | [`27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1`](workflows/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1/)；运行 `run_41129` | 当前采用的增强子主体流程。增强子坐标是完整 H3K27ac peak 合并后的共识区间；motif 扫描中心为 ATAC 原始 summit。作图候选为 **AME、H3K27ac 1.5 倍、500 bp**。同时保留 1.2 倍/200 bp 等敏感性结果。 |
| P1 启动子候选和旧版状态 | [`34_P1_promoter_strength_weaken_v1`](workflows/34_P1_promoter_strength_weaken_v1/)；`run_42216` | 上游已完成，P1 定义为完整 H3K27ac peak 与 TSS±2 kb **重叠**，ATAC/H3K4me3 同时间支持；状态阈值 1.5 倍。**不是**后续“完整落入 + 每 gene 一个 primary promoter”的规则。 |
| 严格 primary promoter、方向分组、五种 mark 与 loop 联合表 | [`38_primary_promoter_direction_only_motif_CRE_table_v1`](workflows/38_primary_promoter_direction_only_motif_CRE_table_v1/)；`run_43558` | **当前最新整合版本，但属于探索性方向阈值版本**。完整 P1 共识区间必须落入关联 gene 的 TSS±2 kb；每 gene 选一个 primary promoter。方向阈值为 1.0（log2FC 正/负），另有 1.2 倍敏感性列。Strength motif 已完成；Weaken 因 GC 背景匹配不足**未测试**。不可把整版 motif 标为“全部最终完成”。 |
| 旧 P1 启动子 motif | [`35_P1_promoter_ATAC_core_HOMER_AME_v1`](workflows/35_P1_promoter_ATAC_core_HOMER_AME_v1/)；`run_42671` | 1.5 倍的 Strength/Weaken AME 与 HOMER 均完成，但使用的是旧 P1 重叠定义，不能直接替代 38 的严格 primary promoter motif。 |
| 增强子 motif–HiChIP–靶基因及通路 | [`32_final_AME_1p5_500bp_HiChIP_targets_draft_v1`](workflows/32_final_AME_1p5_500bp_HiChIP_targets_draft_v1/)；`run_42014`，再经 [`34_unbiased_enhancer_motif_gene_pathway_v1`](workflows/34_unbiased_enhancer_motif_gene_pathway_v1/) 整理 | 32 比较 TSS、P1、P2 三种 promoter 定义；下游 `build_all_results.py` 对应**全部 21 个显著 motif 家族**的 v2 整理包。仍沿用旧 P1 重叠定义，尚未用 38 的严格 primary promoter 重新跑完整 motif–gene 连接。 |

上表中的“最终”是指**当前选用的具体分析环节**，不是说不同阈值、不同启动子定义的文件可以混用。基于合并 BAM 的 0/6 h 信号比较没有生物学重复层面的差异检验 P 值或 FDR。

## 代码目录

| 路径 | 用途 | 归类 |
| --- | --- | --- |
| `workflows/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1/` | 从 H3K27ac/ATAC peak 和 BAM 生成增强子、共识 ATAC summit 与 AME/HOMER | 当前增强子主体 |
| `workflows/29_CRE_definition_sensitivity_export_v1/` | 不读 BAM，比较多个 TSS 窗口及 promoter/enhancer 定义并输出 BED | 方法敏感性 |
| `workflows/31_AME_1p5_500bp_figure_layout_draft_v1/` | 增强子 AME 图的多版排版试验脚本 | 作图草稿；**不保证复现人工调整的 v7 图** |
| `workflows/32_final_AME_1p5_500bp_HiChIP_targets_draft_v1/` | FIMO 位点、HiChIP loop、promoter、gene 与功能注释 | 下游候选连接 |
| `workflows/33_promoter_definition_result_figures_v1/` | 三种 promoter 定义的比较图 | 展示/敏感性 |
| `workflows/34_P1_promoter_strength_weaken_v1/` | P1 共识启动子及 H3K27ac 1.5 倍状态 | 38 的上游 |
| `workflows/34_unbiased_enhancer_motif_gene_pathway_v1/` | `build_plot_data.py` 为 v1 代表 hub；`build_all_results.py` 为 v2 全部 motif 家族 | 下游本地整理；两个脚本不可混称 |
| `workflows/35_P1_promoter_ATAC_core_HOMER_AME_v1/` | 旧 P1 1.5 倍状态的 500 bp 启动子 motif | 历史完成版 |
| `workflows/36_linked_CRE_promoter_mark_dynamics_v1/` | 五种 mark 的 0/6 h 计数与 loop 联合矩阵 | 38 的上游计数 |
| `workflows/38_primary_promoter_direction_only_motif_CRE_table_v1/` | 完整落入 TSS、每 gene 一 primary promoter、方向状态、loop 和 motif | 最新探索性整合版 |

同名数字 **34** 是两件不同的事：`34_P1_promoter_strength_weaken_v1` 是服务器启动子计数作业；`34_unbiased_enhancer_motif_gene_pathway_v1` 是本地下游绘图数据整理。请按**全名**引用。

## 服务器位置与复制来源

结果根目录：

```text
/mnt/nas/swine/home/jxwu/Projects/A/07.HiChIP/NFkB_target_loci_0h_6h_v2/06_chromatin_state_v1/
```

集群脚本根目录：

```text
/mnt/nas/swine/home/jxwu/Projects/A/scripts_new/
```

此处的 27、29、34-P1、35、36、38 脚本复制自相应完成作业目录内的 `code/` 快照；34-P1、35、36、38 与本机 `C:\Users\Luca\Documents\A\scripts_new\` 同名脚本的 SHA-256 已核对一致。31、32、33、34-unbiased 取自本机 `scripts_new/`；**这几套本地工作源码不宣称是服务器作业的不可变代码快照**。详细运行号、结果路径、输入/输出和版本区别见仓库中的 [复刻顺序](REPRODUCIBILITY.md) 与原 `workflows/<版本>/README.md`。

## 使用前先检查

1. 先读 [复刻顺序](REPRODUCIBILITY.md)，确认要跑的是 27 的增强子分析、35 的旧 P1 motif，还是 38 的严格 primary promoter 分析；三者前景不是同一批区域。
2. 检查各版本 `config.json` 中的**服务器绝对路径**、参考基因组、motif 库及源作业目录。脚本不随仓库携带大数据；不要直接在管理节点运行读取 BAM 的 `run.py`。
3. 把原始结果的 `run_metadata.json`、`completed.json`、`motif_status.json` 与新运行的文件核对。`completed.json` 不能单独证明每个 motif 对比都完成。
4. 不要把 38 的方向阈值 `1.0` 写成“变化至少一倍”或统计显著差异；它只划分 H3K27ac log2FC 的正负。

原五个结果文件夹及各文件的含义，见仓库内的 [结果目录与文件索引](docs/RESULT_DIRECTORY_GUIDE.md)。这是说明文档副本；大结果文件仍留在服务器与本地原结果目录，不随代码包上传。
