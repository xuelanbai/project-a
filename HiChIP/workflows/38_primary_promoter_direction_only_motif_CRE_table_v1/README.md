# Primary promoter and directional CRE analysis

## Definitions

### Promoter candidates

The starting P1 promoter intervals are complete H3K27ac consensus intervals
with the previously validated same-time ATAC and H3K4me3 support. A
promoter--gene candidate is retained only when the complete interval is inside
that gene's `TSS +/- 2 kb` window.

### One primary promoter per gene

Each gene contributes at most one promoter--gene association. Candidates are
ranked independently within each gene:

1. an interval containing the TSS is preferred;
2. otherwise, the interval with the smallest distance to the TSS is preferred;
3. ties use the mean 0 h/6 h H3K4me3 signal density, then ATAC signal density;
4. the stable region ID is the final deterministic tie-breaker.

This is a gene-centric restriction. The same genomic promoter interval may be
the selected primary promoter for more than one nearby/bidirectional gene. The
complete candidate table, ranks and selection flags are retained for audit;
candidates are not merged and peak coordinates are not truncated.

### H3K27ac directional states

The exploratory main analysis uses a linear fold boundary of 1:

- Strength: `6 h / 0 h > 1` (`H3K27ac_log2FC > 0`);
- Weaken: `6 h / 0 h < 1` (`H3K27ac_log2FC < 0`);
- Stable: the ratio is exactly 1 or a directional-support check fails.

Strength and Weaken additionally require H3K27ac peak support at the active
time point and active-time H3K27ac and ATAC raw tags both `>=5`.

The original 1.2-fold rule is exported as a sensitivity state column and count
table. Motif enrichment is run on the direction-only groups to retain enough
primary promoters for the exploratory downstream analysis.

A linear boundary of 1 is not an effect-size cutoff: it separates regions only
by direction. Very small positive and negative changes can enter Strength and
Weaken. Results must therefore be labelled `direction-only`, not significant
or differential, and should be checked against the retained 1.2-fold states.

## Reused inputs

- run 41129: final ATAC consensus-core table;
- run 42216: P1 promoter regions and promoter--gene candidates;
- run 43021: 0 h/6 h H3K27ac, H3K4me3, Pol II, ATAC and H3K27me3 matrices;
- raw 0 h and 6 h FitHiChIP BEDPE files;
- Sscrofa11.1 FASTA and JASPAR 2026 vertebrate motif database.

BAMs are not recounted because run 43021 already counted every source promoter
and final enhancer on identical coordinates. The new definition only filters
and ranks existing promoter candidates and relabels states. Input paths and
SHA-256 checksums are recorded in `run_metadata.json`.

## Motif analysis

The primary direction-only promoter states use 500-bp windows centered on the
representative summit of the final ATAC consensus cores. AME and HOMER use
GC-matched nonforeground controls from the same primary-promoter universe.
HOMER is restricted to known vertebrate motifs; de novo motif discovery is
disabled. ATAC direction is reported but is not a foreground filter.

## Run

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/38_primary_promoter_direction_only_motif_CRE_table_v1
bash check_env.sh
bash submit.sh
```

The job requests one compute node, four CPUs and 32 GB RAM. Do not execute
`run.py` directly on the management node.

## Main outputs

- `01_strict_promoter_definition_QC.tsv`: containment and primary-selection
  counts, including genes with multiple candidates.
- `01_state_counts.tsv`: enhancer and primary-promoter state counts for the
  direction-only split and the 1.2-fold sensitivity rule.
- `02_strict_promoter_candidate_gene_ranking.tsv`: every fully contained
  candidate, its rank and its primary-selection flag.
- `02_primary_promoter_region_state_mark.tsv`: one row per unique selected
  primary-promoter interval, both state thresholds, all mark changes and genes.
- `02_primary_promoter_gene_state_mark.tsv`: exactly one promoter association
  per gene; one interval can still represent multiple genes.
- `03_enhancer_region_state_mark.tsv`: every final enhancer, both state
  thresholds, all mark changes and linked genes.
- `04_enhancer_promoter_gene_loop_mark.tsv.gz`: enhancer--loop--primary
  promoter--gene records with loop time and mark changes for both regions.
- `BED/`: all and main-threshold state-specific BED4 files plus column notes.
- `05_primary_promoter_motif_direction_only/`: AME/HOMER tables, GC matching,
  region membership, logs, and q-ranked/enrichment-ranked PNG and PDF figures.
- `run_metadata.json`: definitions, counts, config and input checksums.

## Interpretation

Signals come from merged BAMs, so changes are descriptive and have no
biological-replicate P value or FDR. Counts are not divided by interval length
for state classification because each region is compared with itself between
time points. Signal density is used only to resolve rare equal-distance
primary-promoter ties across intervals of different lengths.
