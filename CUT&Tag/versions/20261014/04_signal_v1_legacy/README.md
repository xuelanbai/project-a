# 初版区域信号定量与审计

multiBigwigSummary 版本保留供历史溯源。区域均值的缺失碱基/分母解释在精确定量 v2 修正。本地审计代码依赖原工作区相对结构，移动后需调整 ROOT/DATA/SCREEN/OUT。

总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。

|文件|用途|状态|
|---|---|---|
|[quantify_cuttag_hypothesis_signals.slurm](quantify_cuttag_hypothesis_signals.slurm)|初版 multiBigwigSummary 区域定量|被精确定量 v2 替代|
|[local_audit/audit_hypothesis_signal_quantification.py](local_audit/audit_hypothesis_signal_quantification.py)|本地区域信号结果审计|相对路径依赖原工作区|
|[local_audit/历史审计说明.md](local_audit/历史审计说明.md)|初版定量修正记录|历史方法记录|
