# 本地统计表工具的版本

v1_template 是模板生成；v1_filled 是填表与解析/检查；fragment_draft 使用 tmp/预置资料，为历史草稿；fragment_preserve 是读取完成 TSV 并保留公式缓存的后续版；v2 为新版参数说明与统计填数。inspection 是早期参考表查看工具。多数脚本写死路径，依赖 @oai/artifact-tool 或 openpyxl，运行前读代码中的输入输出。

总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。

|文件|用途|状态|
|---|---|---|
|[v1_template/build_cuttag_macs2_stats.mjs](v1_template/build_cuttag_macs2_stats.mjs)|MACS2 项目统计模板生成|本地历史模板|
|[v1_filled/fill_workbook1_cuttag.mjs](v1_filled/fill_workbook1_cuttag.mjs)|原 Sheet3 填统计|本地已填表对应版本|
|[v1_filled/parse_cuttag_inputs.py](v1_filled/parse_cuttag_inputs.py)|解析旧统计输入文本|本地辅助|
|[v1_filled/inspect_multimqc_blocks.py](v1_filled/inspect_multimqc_blocks.py)|检查 MultiQC 统计字段|本地辅助|
|[v1_filled/inspect_current_workbook.py](v1_filled/inspect_current_workbook.py)|检查工作簿及数值|本地辅助|
|[v1_filled/inspect_workbook1_sheet3.mjs](v1_filled/inspect_workbook1_sheet3.mjs)|Sheet3 片段/IgG 列预览|本地辅助|
|[v1_fragment_draft/update_fragment_igg.mjs](v1_fragment_draft/update_fragment_igg.mjs)|tmp 与预置值的片段/IgG 填表草稿|历史草稿；用 preserve 版溯源正式填表|
|[v1_fragment_preserve/update_fragment_igg_preserve.py](v1_fragment_preserve/update_fragment_igg_preserve.py)|完成 TSV 的片段均值/IgG 命名，保留缓存|后续修订|
|[v1_fragment_preserve/check_fragment_igg_workbook.py](v1_fragment_preserve/check_fragment_igg_workbook.py)|核对片段/IgG 填表前后|本地辅助|
|[v2/update_cuttag_v2_workbook.mjs](v2/update_cuttag_v2_workbook.mjs)|v2 参数说明、--fill-v2 新统计|依赖历史 JSON.GZ 与工作簿|
|[inspection/inspect_reference.mjs](inspection/inspect_reference.mjs)|早期参考工作簿查看|历史辅助；只保存代码|
|[inspection/artifact/inspect_reference.mjs](inspection/artifact/inspect_reference.mjs)|早期 artifact 参考工作簿查看|另一个内容不同的原版本|
