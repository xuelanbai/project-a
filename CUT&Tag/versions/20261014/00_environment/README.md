# 环境准备与离线预热

prewarm 默认检查字面量 Conda 规格、conda-meta 和 SEACR 执行文件。安装模式应在有网管理节点运行。Nextflow 内部 hash 依赖所用版本。缺失的 cutandtag_env.sh/local_conda.config 见总说明。

总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。

|文件|用途|状态|
|---|---|---|
|[prewarm_cutandtag_envs.nf](prewarm_cutandtag_envs.nf)|离线环境预热/检查|历史修正后版本|
