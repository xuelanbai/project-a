# broad 与 nolambda 参数敏感性

同一脚本用 broad/nolambda 参数分别提交。broad 输出 broadPeak/gappedPeak，nolambda 输出 narrowPeak。比较集固定为 q05 job-51897；job 52126 是已提供的 nolambda 结果，broad 完成结果尚未提供。

总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。

|文件|用途|状态|
|---|---|---|
|[cuttag_macs2_broad_nolambda_sensitivity.slurm](cuttag_macs2_broad_nolambda_sensitivity.slurm)|broad/nolambda 两模式独立分析|nolambda 52126 成功；broad 未核验|
