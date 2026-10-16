#!/bin/bash -l
#SBATCH --job-name=cuttag_signal_qc
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-signal-qc-%j.out

set -euo pipefail

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
RESULTS="$PROJECT/results_pig_v2"
SCREEN="$RESULTS/hypothesis_screen_current_peaks"
BWROOT="$RESULTS/03_peak_calling/03_bed_to_bigwig"
PPROOT="$RESULTS/custom_tracks/macs2_ppois"
OUT="$RESULTS/hypothesis_signal_quantification"
DT_BIN="$HOME/apps/cutandtag/conda/envs/env-0d01150023c68e439ebff3eb959d98bf/bin"
MBW="$DT_BIN/multiBigwigSummary"

[[ -x "$MBW" ]] || { echo "Missing executable: $MBW" >&2; exit 1; }
for bed in \
    "$SCREEN/annotations/Ensembl115.rRNA.bed" \
    "$SCREEN/annotations/Ensembl115.5S_rRNA.bed" \
    "$SCREEN/annotations/candidate_genes.plusminus2kb.bed" \
    "$SCREEN/peak_sets/macs2_119_Pol2_rep1-rep2_shared.bed" \
    "$SCREEN/peak_sets/seacr_119_Pol2_rep1-rep2_shared.bed"
do
    [[ -s "$bed" ]] || { echo "Missing or empty BED: $bed" >&2; exit 1; }
done

find_bw() {
    local pattern="$1"
    local -a matches
    mapfile -t matches < <(find -L "$BWROOT" -type f -name "$pattern" -print)
    (( ${#matches[@]} == 1 )) || {
        echo "Expected one CPM bigWig for $pattern; found ${#matches[@]}" >&2
        printf '%s\n' "${matches[@]}" >&2
        return 1
    }
    [[ -s "${matches[0]}" ]] || return 1
    printf '%s\n' "${matches[0]}"
}

declare -A BW PBW
for state in 119 137; do
    BW["${state}_IgG"]=$(find_bw "3D4_${state}_IgG_R1.bigWig")
    for mark in Tag Pol2; do
        for rep in 1 2; do
            BW["${state}_${mark}_${rep}"]=$(find_bw "3D4_${state}_${mark}_R${rep}.bigWig")
        done
        PBW["${state}_${mark}"]="$PPROOT/3D4_${state}_${mark}_rep1-rep2_pooled_vs_3D4_${state}_IgG_rep1.MACS2_ppois_Sscrofa11.1.bigWig"
        [[ -s "${PBW[${state}_${mark}]}" ]] || {
            echo "Missing pooled IgG-relative bigWig: ${PBW[${state}_${mark}]}" >&2
            exit 1
        }
    done
done

if [[ "${1:-}" == --check ]]; then
    echo "Signal inputs present: 10 CPM bigWigs, 4 pooled IgG-relative MACS2 p-value bigWigs, and hypothesis BED files."
    echo "Outputs will contain regional mean values from multiBigwigSummary; no network is required."
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: bash $0 --check OR sbatch $0" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" && "$(hostname)" != mn-* ]] || {
    echo "Submit this script with sbatch." >&2
    exit 1
}

mkdir -p "$OUT"

all_cpm=(
    "${BW[119_Tag_1]}" "${BW[119_Tag_2]}" "${BW[119_IgG]}"
    "${BW[137_Tag_1]}" "${BW[137_Tag_2]}" "${BW[137_IgG]}"
    "${BW[119_Pol2_1]}" "${BW[119_Pol2_2]}"
    "${BW[137_Pol2_1]}" "${BW[137_Pol2_2]}"
)
all_cpm_labels=(
    119_Tag_rep1 119_Tag_rep2 119_IgG
    137_Tag_rep1 137_Tag_rep2 137_IgG
    119_Pol2_rep1 119_Pol2_rep2
    137_Pol2_rep1 137_Pol2_rep2
)
all_pp=(
    "${PBW[119_Tag]}" "${PBW[137_Tag]}"
    "${PBW[119_Pol2]}" "${PBW[137_Pol2]}"
)
all_pp_labels=(
    119_Tag_pooled_vs_IgG 137_Tag_pooled_vs_IgG
    119_Pol2_pooled_vs_IgG 137_Pol2_pooled_vs_IgG
)

summarize_bed() {
    local name="$1" bed="$2"
    "$MBW" BED-file --BED "$bed" \
        --bwfiles "${all_cpm[@]}" --labels "${all_cpm_labels[@]}" \
        --numberOfProcessors 2 \
        --outFileName "$OUT/${name}.CPM.npz" \
        --outRawCounts "$OUT/${name}.CPM.regional_means.tsv"
    "$MBW" BED-file --BED "$bed" \
        --bwfiles "${all_pp[@]}" --labels "${all_pp_labels[@]}" \
        --numberOfProcessors 2 \
        --outFileName "$OUT/${name}.pooled_MACS2_ppois.npz" \
        --outRawCounts "$OUT/${name}.pooled_MACS2_ppois.regional_means.tsv"
}

summarize_bed rRNA "$SCREEN/annotations/Ensembl115.rRNA.bed"
summarize_bed 5S_rRNA "$SCREEN/annotations/Ensembl115.5S_rRNA.bed"
summarize_bed candidate_genes_plusminus2kb "$SCREEN/annotations/candidate_genes.plusminus2kb.bed"

for caller in macs2 seacr; do
    bed="$SCREEN/peak_sets/${caller}_119_Pol2_rep1-rep2_shared.bed"
    "$MBW" BED-file --BED "$bed" \
        --bwfiles \
            "${BW[119_Tag_1]}" "${BW[119_Tag_2]}" "${BW[119_IgG]}" \
            "${BW[119_Pol2_1]}" "${BW[119_Pol2_2]}" \
        --labels 119_Tag_rep1 119_Tag_rep2 119_IgG 119_Pol2_rep1 119_Pol2_rep2 \
        --numberOfProcessors 2 \
        --outFileName "$OUT/${caller}.119_Pol2_replicated_peaks.CPM.npz" \
        --outRawCounts "$OUT/${caller}.119_Pol2_replicated_peaks.CPM.regional_means.tsv"
    "$MBW" BED-file --BED "$bed" \
        --bwfiles "${PBW[119_Tag]}" "${PBW[119_Pol2]}" \
        --labels 119_Tag_pooled_vs_IgG 119_Pol2_pooled_vs_IgG \
        --numberOfProcessors 2 \
        --outFileName "$OUT/${caller}.119_Pol2_replicated_peaks.pooled_MACS2_ppois.npz" \
        --outRawCounts "$OUT/${caller}.119_Pol2_replicated_peaks.pooled_MACS2_ppois.regional_means.tsv"
done

cat > "$OUT/README.txt" <<'EOF'
Values are mean bigWig signal over each BED interval.
CPM files compare individual biological replicates and matched IgG coverage.
pooled_MACS2_ppois files contain MACS2 -log10(P) signal from pooled IP versus matched IgG.
These values test reproducible IgG-relative chromatin enrichment; they do not prove physical protein-protein interaction.
EOF

tar -czf "$RESULTS/hypothesis_signal_quantification.tar.gz" \
    -C "$RESULTS" "$(basename "$OUT")"

echo "Completed: $OUT"
echo "Transfer: $RESULTS/hypothesis_signal_quantification.tar.gz"
