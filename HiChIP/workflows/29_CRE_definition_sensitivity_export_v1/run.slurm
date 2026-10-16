#!/bin/bash -l
#SBATCH --job-name=CRE_definition_export
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=CRE_definition_export_%j.out
#SBATCH --error=CRE_definition_export_%j.err
set -euo pipefail

BUNDLE_DIR="${BUNDLE_DIR:-${SLURM_SUBMIT_DIR:?}}"
CONFIG="${CONFIG:-$BUNDLE_DIR/config.json}"
ANALYSIS_ENV="${ANALYSIS_ENV:-/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus}"
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1

echo "job=${SLURM_JOB_ID:?} host=$(hostname) started=$(date --iso-8601=seconds)"
"$ANALYSIS_ENV/bin/python" -u "$BUNDLE_DIR/export_and_count_definitions.py" --config "$CONFIG"
echo "completed=$(date --iso-8601=seconds)"
