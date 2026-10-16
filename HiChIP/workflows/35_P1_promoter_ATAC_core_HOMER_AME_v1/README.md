# P1 promoter ATAC-core motif enrichment

## Fixed analysis definition

This task applies the final enhancer motif standard to the already classified
P1 promoter set. It does **not** redefine P1 promoters or recalculate promoter
states.

- Promoter states: archived P1 consensus state matrix, H3K27ac 1.5-fold rule.
- Active-state H3K27ac and ATAC raw tags: both must be at least 5.
  The archived data contain no promoter below this ATAC threshold, so the
  explicit check excludes zero promoters.
- ATAC direction: annotation only, never a foreground filter.
- Motif center: the exact representative summit in the final enhancer ATAC
  consensus-core library (`run_41129`). Reusing this coordinate library is
  intentional so promoter and enhancer results are directly comparable.
- Motif window: 500 bp centered on the representative ATAC summit.
- Multiple cores: all distinct ATAC cores inside a P1 promoter are retained.
- Background: non-foreground cores from the same P1 promoter universe that
  satisfy the corresponding active-time support/read criteria; foreground
  windows and overlapping control windows are excluded, followed by 1:1 GC
  matching without replacement (maximum GC difference 0.05).
- AME: JASPAR2026 vertebrate, Fisher/max scoring.
- HOMER: vertebrate known motifs only; de novo is disabled.
- Significance: within-contrast BH q <= 0.05 and enrichment ratio > 1.
- Display: motif-family representatives, at most 20 per panel. Both q-ranked
  and enrichment-ranked PNG/PDF figures are written.

## Why the ATAC cores are reused

`run_41129/02_ATAC_consensus_cores.tsv` was generated from the original 0 h and
6 h ATAC BAMs by the archived version-27 code. It contains the shared
cross-time consensus peaks and the original MACS3 summit selected from the
time point with higher normalized ATAC signal. Reusing this table prevents a
second MACS3 call from creating small coordinate differences between promoter
and enhancer motif analyses. The upstream BAM-to-core procedure remains fully
reproducible in the version-27 code snapshot.

## Run on the cluster

Upload this directory to:

```text
/mnt/nas/swine/home/jxwu/Projects/A/scripts_new/35_P1_promoter_ATAC_core_HOMER_AME_v1
```

Then run only the light environment check on the management node and submit
the real analysis through Slurm:

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/35_P1_promoter_ATAC_core_HOMER_AME_v1
bash check_env.sh
bash submit.sh
```

## Main outputs

- `01_P1_state_and_ATAC_filter_QC.tsv`: state counts, minimum ATAC counts and
  the explicit >=5-filter audit.
- `02_ATAC_core_to_P1_promoter.tsv`: all representative ATAC summits mapped to
  P1 promoters.
- `regions/`: Strength/Weaken promoter BEDs, mapped core tables and 500-bp BEDs.
- `03_motif_sampling.tsv`: foreground promoter/core counts and GC-match QC.
- `04_AME_all.tsv`, `04_HOMER_all.tsv`: complete motif results.
- `04_*_representative_families.tsv`: one representative per motif family.
- `04_*_significant_representative_families.tsv`: significant enriched family
  representatives used for the result panels.
- `Figure_*_q_rank.{png,pdf}` and
  `Figure_*_enrichment_rank.{png,pdf}`: paired Weaken/Strength panels.
- `motifs/`: foreground/control FASTA files, matching pairs and native raw
  AME/HOMER output directories.
- `run_metadata.json`, `code/`: input hashes, parameters and a code snapshot.
