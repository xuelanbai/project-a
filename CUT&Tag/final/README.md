# CUT&Tag 主分析最终代码副本

版本：**`v2-main-exact-qc-20261016`**。整理日期：2026-10-16。

这份目录可以单独复制保存。它包含已经用于项目主分析的猪参考 v2 和精确信号定量 v2 代码，没有在整理时改变其算法、参数或硬编码路径。参数探索脚本另存在 GitHub 的历史版本目录。

## 文件和顺序

|顺序|文件（均相对本目录）|用途|主要输出|
|---|---|---|---|
|0|`server/offline/prewarm_cutandtag_envs.nf`|管理节点检查/预热环境|共享 Conda cache；Nextflow 运行记录|
|1|`server/config/run_cutandrun_pig_v2.slurm`|nf-core 猪参考主流程|`results_pig_v2/`、gsize、trace/report|
|2|`server/config/postprocess_cutandtag_pig_v2.slurm`|pooled ppois、两重复共同峰|`custom_tracks/macs2_ppois/`、`reproducible_peaks_anyOverlap/`|
|3a|`server/config/run_cuttag_tss_v2.slurm`|全基因 TSS±3 kb profiles|`custom_tracks/TSS_allGenes/`|
|3b|`server/config/screen_cuttag_interaction_hypotheses.slurm`|现有峰候选注释/共占位|`hypothesis_screen_current_peaks/`|
|4|`server/config/quantify_cuttag_signal_qc_v2.py` 与 `.slurm`|精确信号与匹配比较|`hypothesis_signal_qc_v2/job-<jobid>/`|
|5|`server/config/collect_cuttag_v2_stats.slurm`|收集已生成统计|`logs/cuttag_v2_table_inputs.json.gz`|
|检查|`server/config/test_collect_cuttag_v2_stats.py`、`test_quantify_cuttag_signal_qc_v2.py`|既有回归测试|测试日志；不读取真实服务器数据|

3a 和 3b 是主流程/轨道完成后的独立分支；4 依赖 3b 的注释、共同峰和步骤 2 的轨道。统计收集范围以原脚本实现为准，不能认为它自动统计了所有后续探索目录。

## 部署

原脚本默认路径如下，包含 `&` 时保持双引号：

```bash
PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
APP="$HOME/apps/cutandtag"
REF="$HOME/Projects/references/pig/Sscrofa11.1"
```

将 offline 下的 `.nf` 放在 `$APP/offline/`，将 config 下的 Slurm/Python 配套文件部署至 `$PROJECT/config/`。部署前备份服务器同名代码并核对来源。Slurm 日志路径、环境 hash、Java module 和数据路径均为原部署值，换服务器时需要同步调整。wrapper 指向 `$PROJECT/config/quantify_cuttag_signal_qc_v2.py`，因此必须核对该实际文件，不能只更新 Slurm。

大多数 Slurm 支持 `bash <script> --check` 后再 `sbatch <script>`；`collect_cuttag_v2_stats.slurm` 没有 `--check` 模式。Conda 预热脚本默认检查，安装需要显式 `--install true` 且在有网管理节点执行。计算节点使用已有共享工具与参考数据。

## 参数和限制

主流程显式采用 MACS2 narrow q≤0.01、同状态 IgG、两个重复支持、CPM 1 bp、IgG factor=1.0、SEACR norm，以及由实际 FASTA 非 N 碱基数估算的 gsize。pooled ppois 展示轨道与单重复调峰后共同区间分开解释。候选基因附近有峰或与 Pol II 共占位不能单独证明蛋白物理互作。

本副本缺少服务器实际 `params.yaml`、`samplesheet.csv`、`cutandtag_env.sh` 和 `local_conda.config`，还需要原始数据、参考基因组和已安装环境。精确定量需要 pyBigWig/pysam。已有历史命令中的 BAMPE/shift 组合、target/control duplicate 处理和有限 rDNA 注释仍需核实。精确定量 job 51892 的匹配比较区域不足，不能据此声称通过背景检验。

这是现有代码的整理冻结版；没有在打包或上传时重新跑服务器流程。检查 `CODE_MANIFEST.json` 获取原件、版本身份和 SHA256，`SHA256SUMS.txt` 校验本目录所有文件。

完整版本说明位于 GitHub：<https://github.com/xuelanbai/project-a/tree/main/CUT%26Tag>。
