# HiChIP / CRE 代码版本标识（2026-10-16）

本目录是从现有分析脚本**复制并归档**的代码，不是重新运行后发布的新算法版本。原始结果、运行日志和大数据没有上传。请先读 [分析入口](README.md) 与 [复刻顺序](REPRODUCIBILITY.md)。

| 分析层 | 当前标记 | 代码与存档运行 | 尚需注意 |
| --- | --- | --- | --- |
| 增强子区间与 motif 主结果 | **当前选用** | [`27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1`](workflows/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1/)，服务器 `run_41129` | 最终展示候选是 AME、H3K27ac **1.5 倍**、ATAC summit **500 bp**；27 仍保留其他阈值与窗口作敏感性分析。 |
| 增强子 motif–HiChIP–基因 | **当前下游候选结果** | [`32_final_AME_1p5_500bp_HiChIP_targets_draft_v1`](workflows/32_final_AME_1p5_500bp_HiChIP_targets_draft_v1/)，`run_42014`；全部 motif 家族由 [`34_unbiased_enhancer_motif_gene_pathway_v1/build_all_results.py`](workflows/34_unbiased_enhancer_motif_gene_pathway_v1/build_all_results.py) 整理 | 这套连接沿用旧版 P1“与 TSS 窗口重叠”定义；尚未按后续严格 primary promoter 定义重算。 |
| P1 启动子 1.5 倍 motif | **旧定义下已完成的对照** | [`35_P1_promoter_ATAC_core_HOMER_AME_v1`](workflows/35_P1_promoter_ATAC_core_HOMER_AME_v1/)，`run_42671` | Strength 和 Weaken 均完成，但不是最新严格 primary promoter 前景。 |
| 严格 primary promoter、mark 和 loop 联合表 | **最新探索版，不能标为完整最终 motif 结果** | [`38_primary_promoter_direction_only_motif_CRE_table_v1`](workflows/38_primary_promoter_direction_only_motif_CRE_table_v1/)，`run_43558` | 每 gene 最多一个完整落入 TSS±2 kb 的 promoter；主阈值只是方向 `>1`/`<1`。Strength motif 完成，Weaken 因 GC 背景配对不足**未测试**。 |

从原始 BAM 复刻增强子主体需运行 27 且**不设置** `REUSE_PEAKS_FROM`；复刻最新严格 promoter 联合表需按 **27 → 34-P1 → 36 → 38** 的依赖顺序运行。38 自身复用 36 的 BAM 计数矩阵，不能单独从原始 BAM 开始。脚本中的 `config.json` 含原集群绝对路径；运行前必须按实际环境核对。所有版本的状态与限制以各次作业的 `run_metadata.json`、QC 表和 `motif_status.json` 为准。
