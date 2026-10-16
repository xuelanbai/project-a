# P1 promoter H3K27ac Strength/Weaken v1

本模块按照最终增强子的同坐标H3K27ac变化框架，建立一版独立的P1启动子状态结果。它不运行motif分析，也不依赖HiChIP连接。

## 本版的核心决定

- P1负责定义“可信的活性启动子相关H3K27ac区域”。
- 完整H3K27ac peak是最终promoter BED坐标，ATAC和H3K4me3只作同时间支持，不截短H3K27ac区间。
- Strength/Weaken只由统一坐标上的H3K27ac 6 h/0 h变化决定，固定阈值为1.5倍。
- 保留活跃端H3K27ac原始tags `>=5`的低信号过滤。
- **不使用ATAC raw tags `>=5`过滤，也不要求ATAC与H3K27ac同向变化。**
- ATAC仍在两个位置保留：P1身份需要同时间ATAC peak支持；统一promoter坐标上的ATAC log2FC作为注释输出。

## 为什么取消ATAC reads >=5

增强子最终代码中的ATAC raw tags `>=5`位于Strength/Weaken前景选择函数中，是一个额外的活跃端质量过滤。ATAC motif窗口实际来源于MACS3调用的ATAC peak及其summit，并不依赖这个数值阈值。候选增强子此前又已经要求同时间ATAC peak重叠，因此两层ATAC过滤存在一定重复。

本模块不扫描promoter motif。P1已经要求同时间ATAC peak与完整H3K27ac peak至少重叠1 bp，所以不再用原始ATAC计数删除promoter。这样可以保留“仍然开放，但H3K27ac发生变化”的启动子。

## 分析步骤

1. 从Ensembl 115 GTF获得每个gene feature的TSS，并建立TSS +/-2 kb窗口。
2. 分别在0 h和6 h构建P1：
   - H3K27ac peak与至少一个gene的TSS窗口重叠；
   - 同时间ATAC peak与完整H3K27ac peak至少重叠1 bp；
   - 同时间H3K4me3 peak与完整H3K27ac peak至少重叠1 bp；
   - 输出坐标仍为完整H3K27ac peak。
3. 取`P1_0h union P1_6h`，只合并真正重叠的完整H3K27ac peaks；首尾相接但不重叠的区间保持分开。
4. 保留每个consensus promoter到原始时间点peak及gene/TSS的对应关系。
5. 使用BEDTools multicov，在相同consensus promoter坐标中统计0 h和6 h的H3K27ac及ATAC合并BAM。
6. 按各BAM主染色体proper-pair alignment总数换算tags per 10M；不除以区间长度，不做中位数中心化。
7. 计算：

   `H3K27ac_log2FC = log2[(6h tags/10M + 0.5)/(0h tags/10M + 0.5)]`

8. 固定1.5倍阈值分类：
   - Strength：`log2FC >= log2(1.5)`、具有6 h P1支持、6 h H3K27ac raw tags >=5；
   - Weaken：`log2FC <= -log2(1.5)`、具有0 h P1支持、0 h H3K27ac raw tags >=5；
   - 其余为Stable，并在`classification_reason`中注明原因。
9. 额外计算ATAC log2FC并标注`ATAC_up/down/unchanged`，但不参与状态筛选。

## 输入文件

所有路径均在`config.json`中明确记录：

- `reference_fasta`及`.fai`：Sscrofa11.1参考基因组；
- `gtf`：Ensembl 115 gene annotation；
- `peak_inputs/H3K27ac_{0h,6h}`：时间点H3K27ac peak BED；
- `peak_inputs/ATAC_{0h,6h}`：v27从合并ATAC BAM重新调用的MACS3 peak；
- `peak_inputs/H3K4me3_{0h,6h}`：时间点H3K4me3 peak BED；
- `bam_inputs/H3K27ac_{0h,6h}`：合并、去重复、排序并建立索引的H3K27ac BAM；
- `bam_inputs/ATAC_{0h,6h}`：合并、去重复、排序并建立索引的ATAC BAM。

本版为了与最终增强子结果直接对应，主分析仍使用合并BAM和描述性1.5倍阈值。两个独立生物学重复暂不进入状态判定；后续可在同一consensus promoter坐标上增加DESeq2/edgeR/DiffBind验证层。

## 运行

上传整个目录后，在管理节点执行轻量检查并提交：

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/34_P1_promoter_strength_weaken_v1
bash check_env.sh
bash submit.sh
```

不要在管理节点直接运行`run.py`。真实BAM读取只在Slurm计算节点发生。

## 主要输出

- `01_input_peak_schema.tsv`：六个peak输入文件的列数和记录数。
- `02_P1_definition_QC.tsv`：0 h/6 h P1逐步过滤计数。
- `P1/02_P1_promoters_0h.bed`、`P1/02_P1_promoters_6h.bed`：时间点P1完整H3K27ac坐标。
- `P1/02_P1_promoters_*_gene_associations.tsv`：时间点peak到gene/TSS的精确关系。
- `P1/03_P1_promoter_consensus_0h_6h.bed`：用于统一定量的P1共识坐标。
- `P1/03_P1_promoter_consensus_membership.tsv`：共识区域到原始peak的关系。
- `P1/03_P1_promoter_consensus_gene_associations.tsv`：共识区域到gene/TSS的关系。
- `04_P1_promoter_signal_and_state_matrix.tsv`：全部计数、归一化信号、H3K27ac/ATAC log2FC、状态和过滤原因。
- `states/05_P1_strength_promoters.bed`：Strength promoter完整共识坐标。
- `states/05_P1_weaken_promoters.bed`：Weaken promoter完整共识坐标。
- `states/05_P1_stable_promoters.bed`：Stable promoter完整共识坐标。
- `05_P1_promoter_state_gene_associations.tsv`：状态、坐标、gene/TSS及信号的完整对应表。
- `05_P1_promoter_state_counts.tsv`：各状态及过滤原因统计。
- `run_metadata.json`：规则、输入哈希、参数、library totals和结果计数。
- `code/`：本次作业使用的完整代码快照。

所有BED均无表头，采用0-based、half-open坐标。BED第4列为稳定的region ID。需要gene symbol、TSS或状态时，应读取带表头TSV，而不是从BED推断。
