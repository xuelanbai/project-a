#!/bin/bash -l
#SBATCH --job-name=final_AME_targets
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=final_AME_targets_%j.out
#SBATCH --error=final_AME_targets_%j.err

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:?}"
PYTHON="${ANALYSIS_ENV:-/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus}/bin/python"
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PYTHON}" -u "${SCRIPT_DIR}/run.py" --config "${SCRIPT_DIR}/config.json"
