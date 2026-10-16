#!/bin/bash -l
#SBATCH --job-name=promoter_defs
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=promoter_defs_%j.out
#SBATCH --error=promoter_defs_%j.err

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:?}"
PYTHON="${ANALYSIS_ENV:-/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus}/bin/python"
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1

"${PYTHON}" -u "${SCRIPT_DIR}/run.py" --config "${SCRIPT_DIR}/config.json" --promoters-only
