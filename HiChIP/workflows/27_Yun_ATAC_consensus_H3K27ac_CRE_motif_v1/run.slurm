#!/bin/bash -l
#SBATCH --job-name=Yun_ATAC_CRE_motifs
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=Yun_ATAC_CRE_motifs_%j.out
#SBATCH --error=Yun_ATAC_CRE_motifs_%j.err
set -euo pipefail

BUNDLE_DIR="${BUNDLE_DIR:-${SLURM_SUBMIT_DIR:?}}"
CONFIG="${CONFIG:-$BUNDLE_DIR/config.json}"
ANALYSIS_ENV="${ANALYSIS_ENV:-/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus}"

export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "job=${SLURM_JOB_ID:?} host=$(hostname) started=$(date --iso-8601=seconds)"
args=(--config "$CONFIG")
if [[ -n "${RESUME_FROM:-}" ]]; then
    args+=(--resume "$RESUME_FROM")
fi
if [[ -n "${REUSE_PEAKS_FROM:-}" ]]; then
    args+=(--reuse-peaks-from "$REUSE_PEAKS_FROM")
fi
"$ANALYSIS_ENV/bin/python" -u "$BUNDLE_DIR/run.py" "${args[@]}"
echo "completed=$(date --iso-8601=seconds)"
