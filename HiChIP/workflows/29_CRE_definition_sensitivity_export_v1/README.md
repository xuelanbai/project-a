# CRE promoter/enhancer definition sensitivity export

This lightweight job directly exports 0 h, 6 h and cross-time consensus BEDs under several explicit definitions, then counts the resulting intervals. It reads peak BED/narrowPeak files and the Ensembl GTF only; it does not read BAMs and does not run motif enrichment.

The configured TSS values 500, 1500 and 2000 mean `TSS ±500 bp`, `TSS ±1500 bp` and `TSS ±2000 bp`.

Promoters are compared under four coordinate/support rules. Enhancers are compared across three TSS windows and either same-time or the existing v27 either-time H3K4me3 exclusion. The `±2 kb + either_time_v27` enhancer result must exactly reproduce the completed v27 source peaks and consensus coordinates.

Run on the cluster:

```bash
cd /mnt/nas/swine/home/jxwu/Projects/A/scripts_new/29_CRE_definition_sensitivity_export_v1
bash submit.sh
```

Results are written to `29_CRE_definition_sensitivity_export_v1/export_<job-id>`. BED files are organized under `BED/promoter` and `BED/enhancer`.

Optional lightweight checks before submission:

```bash
python test_small.py
python test_integration.py
```
