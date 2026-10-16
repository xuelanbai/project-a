#!/bin/bash -l
#SBATCH --job-name=cuttag_signal_qc_v2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-signal-qc-v2-%j.out

set -euo pipefail
trap 'code=$?; printf "FAILED exit=%s line=%s command=%s\n" "$code" "$LINENO" "$BASH_COMMAND" >&2; exit "$code"' ERR

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
REF="$HOME/Projects/references/pig/Sscrofa11.1"
SCRIPT="$PROJECT/config/quantify_cuttag_signal_qc_v2.py"
PYTHON="$HOME/apps/cutandtag/conda/envs/env-0d01150023c68e439ebff3eb959d98bf/bin/python"
export PYTHONNOUSERSITE=1

[[ -x "$PYTHON" ]] || { echo "Missing deeptools environment Python: $PYTHON" >&2; exit 1; }
[[ -s "$SCRIPT" ]] || { echo "Upload Python helper first: $SCRIPT" >&2; exit 1; }
if [[ "${1:-}" == --check && $# -eq 1 ]]; then
    "$PYTHON" "$SCRIPT" --project "$PROJECT" --reference "$REF" --check
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: bash $0 --check OR sbatch $0" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" && "$(hostname)" != mn-* ]] || {
    echo "Submit with sbatch; computation is not run on the management node." >&2
    exit 1
}
ROOT="$PROJECT/results_pig_v2/hypothesis_signal_qc_v2"
OUT="$ROOT/job-${SLURM_JOB_ID}"
mkdir -p "$ROOT"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1

echo "START $(date -Is); no network required; output=$OUT"
"$PYTHON" "$SCRIPT" --project "$PROJECT" --reference "$REF" --out "$OUT"
[[ -s "$OUT/COMPLETED.txt" ]] || { echo "Completion marker missing" >&2; exit 1; }
ARCHIVE="$ROOT/hypothesis_signal_qc_v2-${SLURM_JOB_ID}.tar.gz"
[[ ! -e "$ARCHIVE" && ! -e "$ARCHIVE.partial" ]] || { echo "Archive already exists; refusing overwrite" >&2; exit 1; }
tar -czf "$ARCHIVE.partial" -C "$ROOT" "$(basename "$OUT")"
mv -- "$ARCHIVE.partial" "$ARCHIVE"
echo "DONE $(date -Is)"
echo "Transfer: $ARCHIVE"
