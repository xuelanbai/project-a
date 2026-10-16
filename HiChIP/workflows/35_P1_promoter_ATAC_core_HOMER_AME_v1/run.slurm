#!/bin/bash -l
#SBATCH --job-name=P1_promoter_motifs
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=P1_promoter_motifs_%j.out
#SBATCH --error=P1_promoter_motifs_%j.err
set -euo pipefail

BUNDLE_DIR="${BUNDLE_DIR:-${SLURM_SUBMIT_DIR:?}}"
CONFIG="${CONFIG:-$BUNDLE_DIR/config.json}"
ANALYSIS_ENV="${ANALYSIS_ENV:-/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus}"

export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "job=${SLURM_JOB_ID:?} host=$(hostname) started=$(date --iso-8601=seconds)"
"$ANALYSIS_ENV/bin/python" -u "$BUNDLE_DIR/run.py" --config "$CONFIG"
echo "completed=$(date --iso-8601=seconds)"
