# AME 1.5-fold / 500-bp figure layout draft v1

This code is intentionally separate from the v27 analysis code. It creates two
combined Strength/Weaken figures from the representative-family AME table:

- q-ranked representatives, capped at 20 motifs per state;
- significant representatives (`BH q <= 0.05` within each contrast), ranked by
  enrichment and capped at 20 motifs per state.

Strength uses a muted warm red and Weaken a muted blue. Both panels use the same
x-axis scale. PNG, vector PDF, SVG and the exact plotted TSV rows are exported.

This draft does not perform HiChIP linking and is not yet merged into the final
raw-BAM-to-figure pipeline.
