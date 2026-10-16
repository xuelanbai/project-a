# 最终代码整理版

|字段|内容|
|---|---|
|版本标识|`v2-main-exact-qc-20261016`|
|整理日期|2026-10-16，Asia/Shanghai|
|代码目录|[`final/`](final/README.md)|
|nf-core 版本|cutandrun 3.2.2，Nextflow 23.10.1 的历史部署|
|主分析|猪参考 v2，MACS2 narrow q≤0.01、匹配 IgG、两个重复正长度相交|
|信号定量|原始 bigWig interval 长度加权积分版本 v2|
|整理行为|复制、分类、编写说明和校验；没有改写已运行分析算法或参数|

“最终版”指本次代码整理选择的主分析冻结版。它由已有版本组合而成，不是新增的 nf-core release，也不把探索分析升级为已经验证的主要结论。

## 包含的脚本

1. `server/offline/prewarm_cutandtag_envs.nf`：有网节点检查/预热已知 Conda 环境。
2. `server/config/run_cutandrun_pig_v2.slurm`：猪参考主流程，non-N gsize、CPM 1 bp、MACS2/SEACR。
3. `server/config/postprocess_cutandtag_pig_v2.slurm`：四组 pooled 对 IgG ppois 轨道及两 caller 的重复共同峰。
4. `server/config/run_cuttag_tss_v2.slurm`：全基因 gene-boundary TSS±3 kb coverage profile。
5. `server/config/collect_cuttag_v2_stats.slurm`：收集已有 published 结果的统计与 provenance。
6. `server/config/screen_cuttag_interaction_hypotheses.slurm`：严格现有峰的候选注释和共占位筛查。
7. `server/config/quantify_cuttag_signal_qc_v2.py` 与 `.slurm`：精确信号定量和匹配比较区域。
8. 两份现有 Python 测试，与各自被测文件同目录保存。

本次副本共 10 个代码文件。它们与 2026-10-14 归档中的对应原件逐字节一致；来源及 SHA256 见 `final/CODE_MANIFEST.json`。

## 探索与历史版本

q0.05 narrow、broad、nolambda、初版 multiBigwigSummary、旧结果 p-value 轨道、本地 Excel/PPT 工具都保存在 `versions/20261014/`。其中 job 51897 为 q0.05；job 52126 为 nolambda；broad 尚未提供成功结果包。它们不在主分析冻结副本内。

## 已知状态

主流程 job 44536 有成功日志；精确定量 job 51892 已提供结果。精确定量的匹配比较区域不足已经记录，不能当作通过了正式背景检验。TSS 独立任务及全部后处理任务仍需以实际日志和完整结果确认；本次上传未在服务器重新执行。

历史 MACS2 BAMPE 调用同时携带 `--nomodel --shift -75 --extsize 150`。整理保留真实命令；正式重跑前需核对 MACS2 2.2.7.1 的实际处理行为。target 来自 markdup BAM、IgG 来自 dedup BAM，duplicate flags 与工具读取行为仍需核实。non-N gsize=2,472,047,704 是当前 FASTA 的近似值。

服务器参数、样本表和环境配置属于外部依赖；缺失清单见 [历史归档的详细说明](versions/20261014/docs/缺失文件和待核实事项.md)。
