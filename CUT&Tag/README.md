# 3D4 猪 CUT&Tag：代码和版本归档

本目录整理本项目现有代码。参考为 Sscrofa11.1 / Ensembl 115；样本包含 119、137 两种状态，Tag 和 Pol2 各两个生物学重复，以及同状态 IgG。

**当前整理定稿：`v2-main-exact-qc-20261016`，代码单独存放在 [`final/`](final/)。** 它包含猪参考主分析 v2、pooled 展示轨道、两个重复共同峰、TSS coverage profile、统计收集、现有峰注释，以及原始 bigWig 区间精确定量 v2。这是已有代码的冻结副本，未在上传时改变分析参数。

## 如何选择代码

|位置|内容|用途|
|---|---|---|
|[`final/`](final/README.md)|单独整理的 10 份当前主分析/环境/测试脚本|优先查看；需要配合服务器现有配置与数据|
|[`FINAL_VERSION.md`](FINAL_VERSION.md)|最终版标识、范围和验证状态|确认这里的“最终版”具体包含什么|
|[`versions/20261014/`](versions/20261014/README.md)|31 份项目脚本、打包工具、历史说明和上游源码快照|追溯旧版、参数探索和本地统计/报告代码|
|[`docs/版本索引.md`](docs/版本索引.md)|按分析阶段导航|查找某版代码及作业记录|
|[`docs/部署说明.md`](docs/部署说明.md)|服务器部署、路径和运行顺序|运行前阅读|

## 分析原则

主分析保持 MACS2 narrow q≤0.01；两个生物学重复各自调峰后，有正长度交叠的区域才列入共同区间。MACS2 与 SEACR 的结果分别保留，报告以 MACS2 为主。rep1/rep2 是生物学重复，FASTQ 的 mate1/mate2 是双端测序读段。

pooled 展示轨道将两个重复的 BAM 一起交给 MACS2，并使用同状态 IgG 建立背景，再通过 `bdgcmp -m ppois` 转为 bigWig。pooled 新峰与两个重复共同峰是两类统计结果；共同区间数也可能不同于获得另一重复支持的原峰数。

q0.05、broad 和 nolambda 均位于历史参数探索目录，不列入本次 `final/` 主分析副本。nolambda 仍输入 IgG，但取消局部背景校正；其峰数增加不能直接当作结合可靠性提高。

## 依赖和归档边界

本仓库上传代码和说明，未上传 FASTQ、BAM、bigWig、测序结果、实验统计表或汇报文件。历史代码中的路径和软件环境 hash 原样保留；原工作区文件没有被移动或修改。`final/CODE_MANIFEST.json` 记录原件来源和 SHA256，`final/SHA256SUMS.txt` 校验最终版副本。

服务器 `params.yaml`、`samplesheet.csv`、`cutandtag_env.sh` 和 `local_conda.config` 等实际文件未找到本地副本，需要从已有服务器部署补齐。当前代码的 BAMPE/shift 参数与 duplicate 处理仍有历史核实事项，见 [最终版本说明](FINAL_VERSION.md)。本次“最终版”表示代码整理定稿，不表示补齐了独立复跑环境或完成了全部生物学验证。

原始 nf-core/cutandrun 3.2.2 源码快照位于 `versions/20261014/90_nfcore_source/`，保留其 LICENSE。它未与服务器实际部署逐文件比较。项目自定义代码本次未额外指定开源许可证。
