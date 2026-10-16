# Yun启发的ATAC consensus core与H3K27ac动态增强子motif分析

## 这一版回答什么

这一版把“增强子活性”和“motif 定位”拆开：

1. H3K27ac peak 的完整 BED 区间仍是增强子坐标来源；输入文件只要前三列有效即可，不要求 narrowPeak 第 10 列 summit。
2. 从 0 h、6 h 合并 ATAC BAM 分别用 MACS3（BAMPE、q=0.01、`--call-summits`）重新调用 peak/summit。
3. 同时间 H3K27ac peak 必须与重新调用的 ATAC peak 至少重叠 1 bp；排除基因 TSS ±2 kb 和任一时间点 H3K4me3 重叠区。
4. 将 0 h、6 h 通过筛选的 H3K27ac peak 合并成不重叠的 H3K27ac consensus CRE；首尾相接但不重叠的区间不合并。
5. 用同一套 consensus CRE 坐标在四个 BAM 中计数；按各 BAM 主染色体 proper-pair 比对数换算 tags per 10M。输入 BAM 已是 `nodup` 的分析成品，因此本版不再额外设置 MAPQ 截断，使定量读取范围与直接读取 BAM 的 MACS3 峰调用保持一致。这里不按区间长度相除，也不做全局中位数中心化。
6. strength/weaken 只由 H3K27ac 相对变化定义，阈值同时输出 1.2 和 1.5 倍；活跃时间点要求 H3K27ac、ATAC 原始计数均 >=5。ATAC log2FC 只作 up/down/unchanged 注释，不作为前景硬过滤条件。
7. 用0 h、6 h重新调用的完整ATAC peak区间建立跨时间ATAC consensus peak：只合并真正重叠的peak，首尾相接但不重叠的peak保持分开。在每个consensus peak的同一坐标中统计两个ATAC BAM的归一化信号；选择信号较高时间点的原始MACS3 peak summit作为代表中心。consensus peak只负责归组和定量，不从合并区间几何中心生成新summit。
8. 一个 H3K27ac CRE 可保留多个相互分离的 ATAC core；motif 的分析单位是“位于 H3K27ac CRE 内的 ATAC core”。
9. 同时扫描代表summit ±100 bp（200 bp）和±250 bp（500 bp）。HOMER只跑vertebrate known motifs；AME使用JASPAR vertebrate，默认不跑de novo。

最终固定输出8个比较：2个H3K27ac变化阈值（1.2、1.5）× 2个方向（strength、weaken）× 2个motif窗口（200 bp、500 bp）。1.2和1.5表示增强子H3K27ac变化倍数阈值，不是motif富集倍数。

这是 pooled BAM 的描述性分析，不能产生生物学重复层面的差异检验 P/FDR。它适合先验证 motif 模式是否稳定；若需要正式差异可及性统计，仍需保留重复分别计数后用 DiffBind/edgeR/DESeq2。

## 环境

不新建环境，复用：

- MACS3：`/mnt/nas/swine/home/jxwu/Projects/FOODIE/envs/macs3_callpeak`
- BEDTools：`/mnt/nas/swine/home/jxwu/Projects/FOODIE/envs/bedtools_intersect`
- Python、HOMER：`/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus`
- AME：`/mnt/nas/swine/home/jxwu/Projects/A/envs/motif`

上传整个目录后先做轻量环境检查，再提交计算节点：

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1
bash check_env.sh
bash submit.sh
```

`check_env.sh` 只检查软件和 Python 模块，不读取大 BAM；`submit.sh` 才通过 Slurm 运行真实数据。

若作业已完成 MACS3、但在后续步骤失败，可传入原结果目录复用两个时间点的 peak：

```bash
bash submit.sh /absolute/path/to/run_JOBID
```

如果要在新结果目录中复用第26版失败作业已经产生的相同MACS3 peak，可用环境变量指定来源，省去重新call peak的时间：

```bash
REUSE_PEAKS_FROM=/mnt/nas/swine/home/jxwu/Projects/A/07.HiChIP/NFkB_target_loci_0h_6h_v2/06_chromatin_state_v1/26_ATAC_consensus_core_HOMER_AME_v1/run_40695 bash submit.sh
```

脚本会把两个`narrowPeak`复制到新run目录并记录来源；其余consensus、计数、分组和motif步骤全部重新运行。若该目录中的peak文件不存在，请直接运行`bash submit.sh`从BAM重新call peak。

登录节点只做轻量检查：

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/27_Yun_ATAC_consensus_H3K27ac_CRE_motif_v1
bash check_env.sh
bash submit.sh
```

不要在登录节点直接运行 `run.py`。

## 关键输出

- `01_input_peak_schema.tsv`：原 H3K27ac/H3K4me3 BED 的列数检查。
- `02_ATAC_called_peaks.tsv.gz`：重新调用的 ATAC peak/summit。
- `02_ATAC_consensus_peak_membership.tsv`：0 h/6 h原始ATAC peak到consensus peak的逐条对应关系。
- `02_ATAC_consensus_peaks.bed`：用于统一定量的ATAC consensus peak完整区间。
- `02_ATAC_consensus_cores.tsv`：consensus peak统一坐标、两个时间点ATAC信号及代表summit的选择依据。
- `02_ATAC_representative_summits.bed`：每个consensus peak最终选择的原始MACS3 summit（1 bp）。
- `03_H3K27ac_consensus_CRE.bed`：完整 H3K27ac consensus CRE。
- `03_CRE_signal_matrix.tsv`：四个 BAM 在统一 CRE 坐标上的计数、归一化信号和 log2FC。
- `04_ATAC_core_to_CRE.tsv`：ATAC core 与 H3K27ac CRE 的映射；一个 CRE 可有多行。
- `fold_1p2/`、`fold_1p5/`：各自strength/weaken CRE BED、ATAC core BED、200/500 bp motif输入及AME/HOMER结果。
- `05_AME_all.tsv`、`05_HOMER_all.tsv`：8个比较的全部motif汇总结果。
- `run_metadata.json`、`code/`：输入/工具/参数及本次代码快照。
