#!/bin/bash -l
#SBATCH --job-name=P1_promoter_state
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=P1_promoter_state_%j.out
#SBATCH --error=P1_promoter_state_%j.err

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:?}"
CONFIG="${CONFIG:-${SCRIPT_DIR}/config.json}"
PYTHON="${ANALYSIS_ENV:-/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus}/bin/python"
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1

"${PYTHON}" -u "${SCRIPT_DIR}/run.py" --config "${CONFIG}"
