# 3D4 猪 CUT&Tag 代码版本归档

整理日期：2026-10-14。对象：119/137 两种状态，Tag 与 Pol2 各两个生物学重复，每种状态一份 IgG；参考为 Sscrofa11.1 / Ensembl 115。

这份代码包保存当前工作区找到的项目脚本，按版本和用途分类。原脚本逐字节复制，未在整理过程中修改分析参数、路径或算法。原位置保留。`manifest.json` 记录来源、分类、大小和 SHA256；`SHA256SUMS.txt` 用于核对归档内容。

这里的 v1/v2 指项目分析版本。两版都使用 nf-core/cutandrun 3.2.2；“信号定量 v2”是后处理代码的修订号，不能误写成 nf-core 软件升级。

## 从哪里开始

先读 [版本对照](docs/版本对照与参数说明.md)，再按需要阅读 [运行顺序与输入输出](docs/运行顺序与输入输出.md)。[缺失文件和待核实事项](docs/缺失文件和待核实事项.md) 列出了当前包没有的服务器配置，以及历史代码中尚待核实的方法问题。

`文件清单.md` 按分类列出全部自定义代码。每个代码目录都有简短说明。服务器计算脚本与本地 Excel/PPT 辅助脚本分别存放。

## 目录

|目录|内容|定位|
|---|---|---|
|`00_environment`|Nextflow Conda 预热及完整性检查|在有网管理节点准备，计算节点离线复用|
|`01_v1_legacy`|旧结果目录上的 MACS2 IgG p-value 轨道脚本|历史版本，gsize 保留 2.7e9|
|`02_v2_main`|猪参考 v2 主流程、pooled 轨道、共同峰、TSS、统计收集|当前严格参数分析的代码记录|
|`03_v2_hypothesis_screen`|现有共同峰注释、候选基因/rRNA、共占位|不重新 call peak|
|`04_signal_v1_legacy`|multiBigwigSummary 初版定量及本地审计|后续被精确定量版本替代|
|`05_signal_v2_exact`|原始 bigWig 区间积分、候选 TSS、匹配比较区域|job 51892 对应版本|
|`06_sensitivity_q05`|局部背景保留，Tag narrow q≤0.05|job 51897 对应探索分析|
|`07_sensitivity_broad_nolambda`|broad 和 nolambda 两个独立模式|nolambda job 52126；broad 成功结果未提供|
|`08_local_tables`|统计表创建、填数、片段长度和 IgG 命名补充|本地工具，包含历史草稿|
|`09_local_reports`|2026-10-13 两版 PPT 生成和检查代码|本地工具，依赖原数据与模板|
|`90_nfcore_source`|本地 nf-core/cutandrun 3.2.2 源码快照|供查阅，不能直接认定等同服务器实际源码|
|`docs`|总说明、参数依据、运行关系和缺失清单|先阅读这些文件|
|`tools`|本次代码打包工具|维护归档，非 CUT&Tag 计算流程|

## 三个容易混淆的结果

1. **pooled bigWig**：把 rep1 和 rep2 两个 BAM 一起交给 MACS2，再以同状态 IgG 建背景，输出 `bdgcmp -m ppois` 的逐碱基信号。它是合并展示轨道，不是“两个 bigWig 相加”，也不是重复一致性的检验。
2. **两个重复共同峰**：分别 call peak 后保留正长度相交的部分。`sharedIntervals` / `shared_any` 是共同覆盖区间；`rep1_supported` / `rep2_supported` 是获得另一重复支持的原峰。区间数和原峰数可能不同。
3. **nf-core CPM bigWig**：coverage 轨道。它的名称和输出目录不能证明已经做了 MACS2/SEACR 调峰或 IgG 局部背景校正。SEACR 的 BED 峰列表也不等同于一种“SEACR p-value bigWig”。

样本名中的 `R1/R2` 在流水线结果里表示生物学重复编号；samplesheet 的 `fastq_1/fastq_2` 才是双端测序的两个 mate。

## 归档范围与使用条件

归档包含找到的自定义代码、既有测试、方法说明和本地上游源码，不包含 FASTQ、BAM、bigWig、完整结果数据、Conda 环境、Node 依赖、Excel/PPT 模板或外部参考基因组。部分代码使用写死的原服务器/本地路径。换机器运行前需要逐项修改路径和核实依赖。

v1 完整提交脚本、服务器 `params.yaml`、`samplesheet.csv`、环境初始化脚本及执行器配置未在本次本地检索中找到，因此这不是一份无需补文件即可从原始 FASTQ 完整复跑的独立环境。没有将聊天截图重写为“原始已运行代码”。

整理和完整性校验不代表所有脚本已在服务器重新执行。已完成作业与待核实版本见版本对照；历史脚本中的参数问题保留在代码里，并在说明中明确标出。
