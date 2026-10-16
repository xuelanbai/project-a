#!/bin/bash -l
#SBATCH --job-name=cuttag_tracks_v2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-tracks-v2-%j.out

set -euo pipefail

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
RESULTS="$PROJECT/results_pig_v2"
APP="$HOME/apps/cutandtag"
OUT="$RESULTS/custom_tracks"
GSIZE_FILE="$RESULTS/custom_metadata/macs_gsize_nonN.txt"

MACS="$APP/conda/envs/env-eb4197d7077a33c4c7b48daeaeb02f19/bin/macs2"
BEDTOOLS="$APP/conda/envs/env-8dafdac3843d8f8edf2f9907029bda58/bin/bedtools"
BEDCLIP="$APP/conda/envs/env-865439183fdfc38f229ad4739ca08e17/bin/bedClip"
BEDGRAPH_TO_BIGWIG="$APP/conda/envs/env-5378deb2578a66435a64a62464238ab2/bin/bedGraphToBigWig"

for executable in "$MACS" "$BEDTOOLS" "$BEDCLIP" "$BEDGRAPH_TO_BIGWIG"; do
    [[ -x "$executable" ]] || { echo "Missing executable: $executable" >&2; exit 1; }
done
[[ -s "$GSIZE_FILE" ]] || { echo "Missing gsize file from first v2 job" >&2; exit 1; }
[[ -d "$RESULTS/03_peak_calling/04_called_peaks/seacr" ]] || { echo "Missing v2 SEACR peaks" >&2; exit 1; }

read -r MACS_GSIZE < "$GSIZE_FILE"
[[ "$MACS_GSIZE" =~ ^[0-9]+$ ]] || { echo "Invalid gsize: $MACS_GSIZE" >&2; exit 1; }

find_one() {
    local directory="$1" filename="$2"
    local -a matches
    mapfile -t matches < <(find -L "$directory" -type f -name "$filename" -print)
    if (( ${#matches[@]} != 1 )); then
        echo "Expected exactly one $filename under $directory; found ${#matches[@]}" >&2
        printf '%s\n' "${matches[@]}" >&2
        return 1
    fi
    if [[ "${3:-}" != allow_empty && ! -s "${matches[0]}" ]]; then
        echo "Empty file: ${matches[0]}" >&2
        return 1
    fi
    printf '%s\n' "${matches[0]}"
}

declare -A BAM SEACR MACS_PEAK
for state in 119 137; do
    BAM["${state}_IgG"]=$(find_one "$RESULTS/02_alignment" "3D4_${state}_IgG_R1.target.dedup.sorted.bam")
    for mark in Pol2 Tag; do
        for rep in 1 2; do
            key="${state}_${mark}_${rep}"
            BAM["$key"]=$(find_one "$RESULTS/02_alignment" "3D4_${state}_${mark}_R${rep}.target.markdup.sorted.bam")
            SEACR["$key"]=$(find_one "$RESULTS/03_peak_calling/04_called_peaks/seacr" "3D4_${state}_${mark}_R${rep}.seacr.peaks.stringent.bed" allow_empty)
            MACS_PEAK["$key"]=$(find_one "$RESULTS/03_peak_calling/04_called_peaks/macs2" "3D4_${state}_${mark}_R${rep}.macs2_peaks.narrowPeak" allow_empty)
            for peak_file in "${SEACR[$key]}" "${MACS_PEAK[$key]}"; do
                awk -F '\t' 'NF<3 || $2!~/^[0-9]+$/ || $3!~/^[0-9]+$/ || $3<=$2 {bad=1} END {exit bad}' "$peak_file" || {
                    echo "Invalid BED coordinates: $peak_file" >&2
                    exit 1
                }
            done
        done
    done
done

echo "Input check passed: 10 BAMs, 8 MACS2 peak files, 8 SEACR peak files (empty peak sets allowed)."
echo "Display output: FOUR pooled rep1+rep2 MACS2 p-value bigWigs, each versus the same-state IgG."
echo "MACS2 effective genome size approximation: $MACS_GSIZE"
if [[ "${1:-}" == --check ]]; then
    echo "Read-only check completed; no analysis was started."
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: bash $0 --check  OR  sbatch $0" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit with sbatch; do not compute on mn-01." >&2; exit 1; }
[[ "$(hostname)" != mn-* ]] || { echo "Refusing to compute on management node." >&2; exit 1; }

mkdir -p "$OUT/reproducible_peaks_anyOverlap" "$OUT/macs2_ppois"
CHROM_SIZES="$OUT/Sscrofa11.1.chrom.sizes"
cut -f1,2 "$HOME/Projects/references/pig/Sscrofa11.1/genome.fa.fai" > "$CHROM_SIZES"
[[ -s "$CHROM_SIZES" ]] || { echo "Failed to create chromosome sizes" >&2; exit 1; }

# Primary reproducibility rule: any positive-length overlap, no fraction cutoff.
# Preserve full original supported peaks AND output their shared genomic bases.
# The primary shared BED is independent of pooled peak calling and its p-values.
COUNTS="$OUT/reproducible_peaks_anyOverlap/peak_counts.txt"
printf 'group\tcaller\trep1_called\trep2_called\trep1_supported\trep2_supported\tshared_intervals\n' > "$COUNTS.partial"
reproducible_peaks() {
    local group="$1" caller="$2" rep1="$3" rep2="$4" control="$5"
    local prefix="$OUT/reproducible_peaks_anyOverlap/${group}_rep1-rep2_vs_${control}.${caller}.2of2_anyOverlap"
    "$BEDTOOLS" intersect -a "$rep1" -b "$rep2" -u > "${prefix}.rep1_supported.bed.partial"
    "$BEDTOOLS" intersect -a "$rep2" -b "$rep1" -u > "${prefix}.rep2_supported.bed.partial"
    # Use only BED3 for an unambiguous intersection, retaining all shared bases.
    cut -f1-3 "$rep1" > "${prefix}.rep1.bed3"
    cut -f1-3 "$rep2" > "${prefix}.rep2.bed3"
    "$BEDTOOLS" intersect -a "${prefix}.rep1.bed3" -b "${prefix}.rep2.bed3" \
        | LC_ALL=C sort -k1,1 -k2,2n \
        | "$BEDTOOLS" merge -i - > "${prefix}.sharedIntervals.bed.partial"
    for suffix in rep1_supported.bed rep2_supported.bed sharedIntervals.bed; do
        mv -- "${prefix}.${suffix}.partial" "${prefix}.${suffix}"
    done
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$group" "$caller" \
        "$(wc -l < "$rep1")" "$(wc -l < "$rep2")" \
        "$(wc -l < "${prefix}.rep1_supported.bed")" \
        "$(wc -l < "${prefix}.rep2_supported.bed")" \
        "$(wc -l < "${prefix}.sharedIntervals.bed")" >> "$COUNTS.partial"
}

for state in 119 137; do
    for mark in Pol2 Tag; do
        group="3D4_${state}_${mark}"
        reproducible_peaks "$group" MACS2 "${MACS_PEAK[${state}_${mark}_1]}" "${MACS_PEAK[${state}_${mark}_2]}" "3D4_${state}_IgG_rep1"
        reproducible_peaks "$group" SEACR "${SEACR[${state}_${mark}_1]}" "${SEACR[${state}_${mark}_2]}" "3D4_${state}_IgG_rep1"
    done
done
mv -- "$COUNTS.partial" "$COUNTS"

run_macs2_ppois() {
    local state="$1" mark="$2" label="$3" target1="$4" target2="$5"
    local control="${BAM[${state}_IgG]}"
    local sample="3D4_${state}_${mark}_${label}_vs_3D4_${state}_IgG_rep1"
    local base="${sample}.MACS2_ppois_Sscrofa11.1"
    local dir="$OUT/macs2_ppois/$sample"
    local bw="$OUT/macs2_ppois/${base}.bigWig"
    local -a treatments=("$target1")
    [[ -z "$target2" ]] || treatments+=("$target2")
    if [[ -s "$bw" ]]; then
        echo "SKIP completed bigWig: $bw"
        return
    fi
    mkdir -p "$dir"
    echo "START MACS2 $sample at $(date -Is)"
    "$MACS" callpeak \
        --nomodel --shift -75 --extsize 150 --keep-dup all \
        -q 0.01 -g "$MACS_GSIZE" -f BAMPE -B \
        -n "$base" --outdir "$dir" \
        -t "${treatments[@]}" -c "$control"
    "$MACS" bdgcmp \
        -t "$dir/${base}_treat_pileup.bdg" \
        -c "$dir/${base}_control_lambda.bdg" \
        -m ppois -o "$dir/${base}.ppois.bedGraph"
    LC_ALL=C sort -S 8G -T "$dir" -k1,1 -k2,2n \
        "$dir/${base}.ppois.bedGraph" > "$dir/${base}.ppois.sorted.bedGraph"
    "$BEDCLIP" "$dir/${base}.ppois.sorted.bedGraph" "$CHROM_SIZES" \
        "$dir/${base}.ppois.sorted.clipped.bedGraph"
    "$BEDGRAPH_TO_BIGWIG" "$dir/${base}.ppois.sorted.clipped.bedGraph" "$CHROM_SIZES" \
        "$dir/${base}.bigWig.partial"
    [[ -s "$dir/${base}.bigWig.partial" ]] || { echo "Empty converted bigWig for $sample" >&2; exit 1; }
    mv -- "$dir/${base}.bigWig.partial" "$bw"
    echo "DONE MACS2 $sample at $(date -Is)"
}

for state in 119 137; do
    for mark in Pol2 Tag; do
        run_macs2_ppois "$state" "$mark" rep1-rep2_pooled \
            "${BAM[${state}_${mark}_1]}" "${BAM[${state}_${mark}_2]}"
    done
done

echo "Completed: FOUR pooled-replicate MACS2 IgG-relative p-value bigWigs, plus MACS2/SEACR 2-of-2 any-overlap peak sets and counts."
echo "Primary reproducible regions: $OUT/reproducible_peaks_anyOverlap/*.sharedIntervals.bed"
echo "Pooled MACS2 peaks are separate exploratory outputs; pooling does not establish replicate support."
