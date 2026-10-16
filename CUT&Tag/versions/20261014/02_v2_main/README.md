# 猪参考主流程 v2、轨道与统计

run→postprocess→TSS 是输入依赖顺序；collect 在需要统计的结果生成后运行。只有 collect 没有 --check 模式。test_collect 与 collect 保持同目录。参数、缓存和重复解释见 docs。

总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。

|文件|用途|状态|
|---|---|---|
|[run_cutandrun_pig_v2.slurm](run_cutandrun_pig_v2.slurm)|猪参考 v2 nf-core driver|44536 成功|
|[postprocess_cutandtag_pig_v2.slurm](postprocess_cutandtag_pig_v2.slurm)|四条 pooled ppois 与两 caller 共同峰|44538 对应后处理|
|[run_cuttag_tss_v2.slurm](run_cuttag_tss_v2.slurm)|全基因 TSS coverage profiles|脚本存在；完成状态待日志核对|
|[collect_cuttag_v2_stats.slurm](collect_cuttag_v2_stats.slurm)|收集 v2 统计至 JSON.GZ|44539 对应任务|
|[test_collect_cuttag_v2_stats.py](test_collect_cuttag_v2_stats.py)|统计收集器既有测试|与 Slurm 同目录|
