# Promoter-definition result figures

This module compares all three promoter definitions from run 42014 without
selecting a preferred definition.

```powershell
python make_figures.py `
  --run "D:\Lucas临时文件\32_final_AME_1p5_500bp_HiChIP_targets_draft_v1\run_42014" `
  --output "C:\Users\Luca\Documents\A\outputs\33_promoter_definition_result_figures_v1"

powershell -ExecutionPolicy Bypass -File render_pngs.ps1 `
  -OutputDirectory "C:\Users\Luca\Documents\A\outputs\33_promoter_definition_result_figures_v1"
```

Outputs include one definition-sensitivity summary, one matched-scale
motif-enhancer-gene matrix for each promoter definition, one combined matrix,
and TSV source data for every plotted value.

The matrix is restricted to the predeclared lipid-regulatory reference panel
and displays genes detected in at least one promoter definition. Cell values
are counts of unique motif-bearing H3K27ac CREs linked to a gene by HiChIP.
Strength motifs are red and Weaken motifs are blue. All three figures use the
same gene order, motif order and color scale.
