# Enhancer–P1 promoter multi-mark dynamics

## Purpose

This task quantifies 0 h and 6 h signal for five assays in the final enhancer
and P1 promoter coordinate systems, then attaches the signal changes to the
HiChIP-linked enhancer–promoter–gene relationships.

The promoter and enhancer definitions are not changed:

- Enhancer: final version-27 H3K27ac consensus CREs.
- Promoter: final P1 full-H3K27ac-peak-support consensus promoters.
- State: the same 1.5-fold H3K27ac rule and active-state H3K27ac/ATAC raw-tag
  requirement used in the archived analyses.
- TSS validation window: fixed gene TSS ±2 kb.

## Signal calculation

Assays: H3K27ac, H3K4me3, Pol II, ATAC and H3K27me3.

All merged BAMs are counted over the exact same interval at both time points
with `bedtools multicov -p -q 0`. Counts are scaled by the number of proper-pair
alignments on pig main chromosomes to tags per 10 million. No interval-length
division and no median centering are applied. Consequently, within-region
0 h/6 h log2 fold changes are valid; absolute signal should not be compared
between regions of very different length as if it were density.

For H3K27ac, H3K4me3, Pol II and H3K27me3 the task also reports:

```text
log2(mark tags-per-10M / IgG tags-per-10M)
```

at each time point and its 6 h minus 0 h change. ATAC has no IgG correction.

This is a pooled/merged-BAM descriptive analysis. It does not generate
biological-replicate P values or FDR values.

## HiChIP linkage

The archived strength and weaken loop–gene tables are remapped to the current
P1 consensus promoter IDs using the same gene ID plus interval overlap. The
output collapses duplicate loops for the same enhancer–promoter–gene relation,
retaining the number of unique loops and best loop q value.

Promoter-state proportions are calculated using unique P1 promoter IDs, not
raw loop rows. Mark concordance uses unique enhancer–promoter pairs so that
genes with multiple loop records do not receive extra weight.

## Run

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/36_linked_CRE_promoter_mark_dynamics_v1
bash check_env.sh
bash submit.sh
```

## Main outputs

- `01_enhancer_mark_signal_matrix.tsv`: every final enhancer and all mark signals.
- `02_P1_promoter_mark_signal_matrix.tsv`: every P1 promoter and all mark signals.
- `03_gene_TSS_pm2kb_mark_signal_matrix.tsv`: gene-level TSS ±2 kb signals.
- `03_source_vs_recount_QC.tsv`: verifies that newly counted H3K27ac/ATAC
  changes reproduce the archived enhancer and promoter state matrices.
- `04_region_mark_long.tsv`: tidy long-format values for downstream plotting.
- `04_region_mark_summary.tsv`: median changes and directional fractions by
  region type, state and assay.
- `05_HiChIP_enhancer_promoter_gene_mark_matrix.tsv`: collapsed linked relations
  with enhancer, promoter and TSS mark changes.
- `06_linked_promoter_state_distribution.tsv`: Strength/Stable/Weaken promoter
  proportions linked from each enhancer state.
- `07_linked_enhancer_promoter_mark_concordance.tsv`: pair-level directional
  agreement plus Pearson and Spearman correlations.
- `08_linked_promoter_mark_summary.tsv`: mark changes in unique linked P1
  promoters, stratified by enhancer state and promoter state.
- `Figure_1_region_mark_log2FC_heatmap.{png,pdf}`.
- `Figure_2_linked_promoter_state_distribution.{png,pdf}`.
- `counts/`, `logs/`, `run_metadata.json` and `code/` for reproducibility.
