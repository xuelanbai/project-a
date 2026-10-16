# Final AME 1.5-fold / 500-bp HiChIP targets (draft)

This downstream module starts from the completed v27 run but recomputes every
motif-site and HiChIP target link. It does not reuse v24 target counts.

## Frozen input definition

- State regions: v27 `fold_1p5/regions/{strength,weaken}_CRE.bed`.
- Motif sequences: v27 500-bp ATAC-core foreground FASTA files.
- Motifs: AME representative families with within-contrast BH q <= 0.05 and
  enrichment ratio > 1. This matches the final displayed AME result.
- FIMO: original JASPAR 2026 vertebrate matrices, P <= 1e-4.
- Loops: 6 h for Strength and 0 h for Weaken, FitHiChIP Q <= 0.01.
- Promoters are evaluated in three parallel definitions (all use TSS +/-2 kb):
  `TSS_pm2kb`, `P1_full_H3_peak_support`, and
  `P2_TSS_local_support_full_H3`. Strength uses 6 h P1/P2 promoters and
  Weaken uses 0 h P1/P2 promoters.

P1 requires a full H3K27ac peak overlapping a gene TSS window and same-time
ATAC plus H3K4me3 support anywhere in that full peak. P2 requires both support
tracks inside the gene-specific intersection of the H3K27ac peak and TSS
window. Both retain the full H3K27ac peak as the promoter coordinate. Exact
peak-to-gene associations are retained separately from coordinate-only BEDs.

The primary link requires the full enhancer H3K27ac CRE to overlap one loop
anchor and the opposite anchor to overlap the selected promoter definition.
The output also records whether the specific ATAC core and FIMO motif site
overlap the CRE-side anchor.

## Run

To export and count only the three promoter definitions first (recommended
for the present sensitivity check):

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/32_final_AME_1p5_500bp_HiChIP_targets_draft_v1
bash submit_promoters.sh
```

After reviewing `02_promoter_definition_counts.tsv`, run the complete motif to
HiChIP target comparison:

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/32_final_AME_1p5_500bp_HiChIP_targets_draft_v1
bash submit.sh
```

The job writes `run_<job-id>/` below the configured output root. The most useful
files are:

- `01_selected_AME_motifs.tsv`: exact motif families used downstream.
- `promoters/TSS_pm2kb/`: gene-level TSS +/-2 kb BED and TSV.
- `promoters/{P1,P2}/promoter_{0h,6h}.bed`: unique coordinate catalogs.
- `promoters/{P1,P2}/*_gene_associations.{bed,tsv}`: the exact gene-linked
  intervals used for HiChIP mapping.
- `promoters/{P1,P2}/promoter_consensus_0h_6h.bed`: coordinate union followed
  by true-overlap merging; it is a reference export, not the time-specific
  HiChIP input.
- `motifs/{strength,weaken}/02_final_core_CRE_membership.tsv`: exact core-to-CRE mapping.
- `motifs/{strength,weaken}/04_FIMO_site_coordinates.tsv.gz`: genomic motif sites.
- `targets/{definition}/{strength,weaken}/05_all_CRE_loop_gene_links.tsv.gz`:
  motif-independent CRE links.
- `targets/{definition}/{strength,weaken}/06_motif_CRE_loop_gene_evidence.tsv.gz`:
  full evidence rows.
- `targets/{definition}/{strength,weaken}/08_candidate_genes.tsv`: gene-by-family summaries with
  primary and strict-anchor independent-CRE counts.
- `targets/{definition}/functions/04_TF_candidate_functions.tsv`: annotations.
- `targets/{definition}/*_AME_{lipid,inflammation,reference_panel}_bridge.tsv`:
  plot-ready matrices.
- `02_promoter_definition_counts.tsv` and `03_target_definition_counts.tsv`:
  direct sensitivity comparisons across the three definitions.
- `run_metadata.json`: definitions, QC and input hashes.

## Input/output note

Promoter inputs are the Ensembl 115 GTF plus the time-matched H3K27ac, ATAC
and H3K4me3 peak BED files listed explicitly in `config.json`. Only their first
three BED columns are required. P1/P2 outputs use complete H3K27ac peak
coordinates; ATAC and H3K4me3 are support filters and never replace or trim
those coordinates. The headerless BED files are zero-based, half-open. BED4
files use column 4 as a stable interval ID; BED6 association files use columns
4/5/6 as gene ID, score 0 and gene strand. The companion TSVs contain headers
and should be used when gene symbols, TSS coordinates or method labels are
needed.

This module is intentionally separate from the all-from-BAM reproduction code.
After the downstream definitions and figures are approved, it can be merged
into that final end-to-end package.
