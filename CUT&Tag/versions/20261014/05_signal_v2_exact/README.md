# 精确信号定量 v2

Python 和 Slurm 必须同时部署。wrapper 调用 $PROJECT/config/quantify_cuttag_signal_qc_v2.py；测试与 Python 放同目录。该目录与旧交付包中的同名 Python/Slurm 是相同 SHA256，不额外制造新版本。历史方法说明中的尚未执行等状态，以本归档版本表更新。

总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。

|文件|用途|状态|
|---|---|---|
|[quantify_cuttag_signal_qc_v2.py](quantify_cuttag_signal_qc_v2.py)|原始 interval 积分与匹配比较|51892 对应代码|
|[quantify_cuttag_signal_qc_v2.slurm](quantify_cuttag_signal_qc_v2.slurm)|精确定量 v2 Slurm wrapper|51892 对应代码|
|[test_quantify_cuttag_signal_qc_v2.py](test_quantify_cuttag_signal_qc_v2.py)|精确定量既有回归测试|与 helper 同目录|
|[历史方法说明.md](历史方法说明.md)|保留当时信号 QC 与敏感性说明|历史状态文字以总版本表更新|
