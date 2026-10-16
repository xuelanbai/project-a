#!/bin/bash -l
#SBATCH --job-name=cuttag_pig_v2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-pig-v2-%j.out

set -euo pipefail

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
REF="$HOME/Projects/references/pig/Sscrofa11.1"
APP="$HOME/apps/cutandtag"
PIPELINE="$APP/pipelines/cutandrun-3.2.2"
PARAMS="$PROJECT/config/params.yaml"
SAMPLES="$PROJECT/config/samplesheet.csv"
CONDA_CONFIG="$APP/conf/local_conda.config"
NEW_RESULTS="$PROJECT/results_pig_v2"

for f in "$REF/genome.fa" "$REF/genome.fa.fai" "$REF/genes.Ensembl_115.gtf.gz" "$PARAMS" "$SAMPLES" "$CONDA_CONFIG" "$PIPELINE/main.nf" "$APP/cutandtag_env.sh"; do
    [[ -s "$f" ]] || { echo "Missing required file: $f" >&2; exit 1; }
done
[[ -d "$PROJECT/work" ]] || { echo "Existing Nextflow work directory missing" >&2; exit 1; }
[[ -d "$PROJECT/logs" ]] || { echo "Existing logs directory missing" >&2; exit 1; }

if [[ "${1:-}" == --check ]]; then
    echo "Inputs present; reference and caller parameters from current params.yaml:"
    grep -nE '^[[:space:]]*(input|genome|fasta|gtf|bowtie2|chrom_sizes|spikein_genome|spikein_fasta|normalisation_mode|normalisation_binsize|peakcaller|macs_gsize|igg_scale_factor|seacr_norm|replicate_threshold):' "$PARAMS" || true
    echo "New results path: $NEW_RESULTS"
    echo "Read-only check finished; no Nextflow job was started."
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: bash $0 --check  OR  sbatch $0" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit using sbatch; do not calculate on mn-01." >&2; exit 1; }
[[ "$(hostname)" != mn-* ]] || { echo "Refusing to calculate on management node." >&2; exit 1; }

module load Java/17.0.15
source "$APP/cutandtag_env.sh"
export NXF_OFFLINE=true
export NXF_OPTS="-Xms1g -Xmx6g"

# Fast documented approximation: count non-N bases in exactly the FASTA used
# by this run. With MAPQ-filtered reads, uniquely mappable length may be
# smaller, so record this value as an approximation, not a species constant.
read -r MACS_GSIZE FASTA_TOTAL < <(
    awk '!/^>/ { total += length($0); n += gsub(/[Nn]/, "") }
         END { printf "%.0f %.0f\n", total-n, total }' "$REF/genome.fa"
)
FAI_TOTAL=$(awk -F '\t' '{ total += $2 } END { printf "%.0f", total }' "$REF/genome.fa.fai")
[[ "$FASTA_TOTAL" == "$FAI_TOTAL" ]] || { echo "FASTA/FAI length mismatch: $FASTA_TOTAL vs $FAI_TOTAL" >&2; exit 1; }
(( MACS_GSIZE > 2000000000 && MACS_GSIZE <= FASTA_TOTAL )) || { echo "Unexpected MACS gsize: $MACS_GSIZE" >&2; exit 1; }
printf 'non_N_bases=%s\tfasta_total=%s\tmethod=non_N_approximation\n' "$MACS_GSIZE" "$FASTA_TOTAL" \
    > "$PROJECT/logs/pig-v2-gsize-${SLURM_JOB_ID}.txt"
mkdir -p "$NEW_RESULTS/custom_metadata"
printf '%s\n' "$MACS_GSIZE" > "$NEW_RESULTS/custom_metadata/macs_gsize_nonN.txt"

cd "$PROJECT"
"$APP/bin/nextflow" run "$PIPELINE/main.nf" \
    -profile mamba \
    -c "$CONDA_CONFIG" \
    -params-file "$PARAMS" \
    --input "$SAMPLES" \
    --fasta "$REF/genome.fa" \
    --gtf "$REF/genes.Ensembl_115.gtf.gz" \
    --outdir "$NEW_RESULTS" \
    --normalisation_mode CPM \
    --normalisation_binsize 1 \
    --igg_scale_factor 1.0 \
    --seacr_norm norm \
    --peakcaller macs2,seacr \
    --use_control \
    --macs_gsize "$MACS_GSIZE" \
    --macs2_qvalue 0.01 \
    --macs2_narrow_peak \
    --consensus_peak_mode group \
    --replicate_threshold 2 \
    --min_peak_overlap 0.000000001 \
    -work-dir "$PROJECT/work" \
    -with-report "$PROJECT/logs/pig-v2-report-${SLURM_JOB_ID}.html" \
    -with-timeline "$PROJECT/logs/pig-v2-timeline-${SLURM_JOB_ID}.html" \
    -with-trace "$PROJECT/logs/pig-v2-trace-${SLURM_JOB_ID}.txt" \
    -resume
