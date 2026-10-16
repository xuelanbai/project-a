# narrow q0.05 参数敏感性

job 51897 对应脚本；Tag 两重复分别调用，保留同状态 IgG 和局部背景。共同区间原则不变，不追加 q0.01 二次过滤，不生成新 bigWig。

总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。

|文件|用途|状态|
|---|---|---|
|[cuttag_macs2_q05_sensitivity.slurm](cuttag_macs2_q05_sensitivity.slurm)|Tag narrow q0.05 与严格集合对照|51897 成功|
