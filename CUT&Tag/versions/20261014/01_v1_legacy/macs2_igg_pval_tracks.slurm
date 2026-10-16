#!/usr/bin/env bash
#SBATCH --job-name=cuttag_macs2_pval
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/macs2-igg-pval-%j.out

set -euo pipefail

# Use --check on mn-01 for read-only input validation. Submit with sbatch for work.
case "${1:-}" in
    --check) CHECK_ONLY=1 ;;
    "")      CHECK_ONLY=0 ;;
    *)       echo "Usage: bash $0 --check  OR  sbatch $0" >&2; exit 2 ;;
esac

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
REF="$HOME/Projects/references/pig/Sscrofa11.1"
APP="$HOME/apps/cutandtag"
OUT="$PROJECT/results/custom_macs2_igg_ppois"

MACS="$APP/conda/envs/env-eb4197d7077a33c4c7b48daeaeb02f19/bin/macs2"
BEDCLIP="$APP/conda/envs/env-865439183fdfc38f229ad4739ca08e17/bin/bedClip"
BEDGRAPH_TO_BIGWIG="$APP/conda/envs/env-5378deb2578a66435a64a62464238ab2/bin/bedGraphToBigWig"

# Intentionally reproduce the original MACS2 callpeak -g 2.7e9 for
# comparability with existing MACS2 peaks. This is NOT a validated pig
# effective genome size and should be reassessed in a separate analysis.
MACS_GSIZE=2700000000

for executable in "$MACS" "$BEDCLIP" "$BEDGRAPH_TO_BIGWIG"; do
    [[ -x "$executable" ]] || { echo "Missing executable: $executable" >&2; exit 1; }
done
[[ -d "$PROJECT/results" ]] || { echo "Missing results directory" >&2; exit 1; }
[[ -s "$REF/genome.fa.fai" ]] || { echo "Missing reference FASTA index" >&2; exit 1; }
[[ -d "$PROJECT/logs" ]] || { echo "Missing logs directory" >&2; exit 1; }

find_one_bam() {
    local filename="$1"
    local -a matches
    mapfile -t matches < <(find -L "$PROJECT/results" -type f -name "$filename" -print)
    if (( ${#matches[@]} != 1 )); then
        echo "Expected exactly one BAM named $filename; found ${#matches[@]}" >&2
        printf '%s\n' "${matches[@]}" >&2
        return 1
    fi
    [[ -s "${matches[0]}" ]] || { echo "Empty BAM: ${matches[0]}" >&2; return 1; }
    printf '%s\n' "${matches[0]}"
}

declare -A BAM
for state in 119 137; do
    BAM["${state}_IgG"]=$(find_one_bam "3D4_${state}_IgG_R1.target.dedup.sorted.bam")
    for mark in Pol2 Tag; do
        for rep in 1 2; do
            key="${state}_${mark}_${rep}"
            BAM["$key"]=$(find_one_bam "3D4_${state}_${mark}_R${rep}.target.markdup.sorted.bam")
        done
    done
done

echo "MACS2: $("$MACS" --version 2>&1)"
echo "Original MACS2 gsize retained: $MACS_GSIZE"
for state in 119 137; do
    echo "3D4_${state}_IgG: ${BAM[${state}_IgG]}"
    for mark in Pol2 Tag; do
        for rep in 1 2; do
            echo "3D4_${state}_${mark}_rep${rep}: ${BAM[${state}_${mark}_${rep}]}"
        done
    done
done

if (( CHECK_ONLY )); then
    echo "CHECK PASS: 10 input BAMs, three executables, and reference index found; no files written."
    exit 0
fi

[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit with sbatch; do not run work on mn-01." >&2; exit 1; }
[[ "$(hostname)" != mn-* ]] || { echo "Refusing to compute on management node." >&2; exit 1; }

mkdir -p "$OUT"
CHROM_SIZES="$OUT/Sscrofa11.1.chrom.sizes"
cut -f1,2 "$REF/genome.fa.fai" > "$CHROM_SIZES"
[[ -s "$CHROM_SIZES" ]] || { echo "Failed to create chrom.sizes" >&2; exit 1; }

MANIFEST="$OUT/manifest.tsv"
if [[ ! -e "$MANIFEST" ]]; then
    printf 'sample\ttreatment_bam_1\ttreatment_bam_2\tcontrol_bam\tgsize\tppois_bigwig\n' > "$MANIFEST"
fi

run_one() {
    local state="$1" mark="$2" label="$3" treatment1="$4" treatment2="$5"
    local control="${BAM[${state}_IgG]}"
    local sample="3D4_${state}_${mark}_${label}_vs_3D4_${state}_IgG_rep1"
    local base="${sample}.MACS2_ppois_Sscrofa11.1"
    local sample_dir="$OUT/$sample"
    local output_bw="$OUT/${base}.bigWig"
    local partial_bw="$sample_dir/${base}.bigWig.partial"
    local raw_bdg="$sample_dir/${base}.ppois.bedGraph"
    local sorted_bdg="$sample_dir/${base}.ppois.sorted.bedGraph"
    local clipped_bdg="$sample_dir/${base}.ppois.sorted.clipped.bedGraph"
    local -a treatments=("$treatment1")
    [[ -z "$treatment2" ]] || treatments+=("$treatment2")

    if [[ -s "$output_bw" ]]; then
        echo "SKIP existing nonempty bigWig: $output_bw"
        return
    fi
    mkdir -p "$sample_dir"
    echo "START $sample at $(date -Is)"

    "$MACS" callpeak \
        --nomodel --shift -75 --extsize 150 --keep-dup all \
        -q 0.01 -g "$MACS_GSIZE" -f BAMPE -B \
        -n "$base" --outdir "$sample_dir" \
        -t "${treatments[@]}" -c "$control"

    "$MACS" bdgcmp \
        -t "$sample_dir/${base}_treat_pileup.bdg" \
        -c "$sample_dir/${base}_control_lambda.bdg" \
        -m ppois -o "$raw_bdg"

    LC_ALL=C sort -S 8G -T "$sample_dir" -k1,1 -k2,2n "$raw_bdg" > "$sorted_bdg"
    "$BEDCLIP" "$sorted_bdg" "$CHROM_SIZES" "$clipped_bdg"
    "$BEDGRAPH_TO_BIGWIG" "$clipped_bdg" "$CHROM_SIZES" "$partial_bw"
    [[ -s "$partial_bw" ]] || { echo "Empty converted bigWig: $partial_bw" >&2; exit 1; }
    mv -- "$partial_bw" "$output_bw"
    [[ -s "$output_bw" ]] || { echo "Empty bigWig: $output_bw" >&2; exit 1; }

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$treatment1" "$treatment2" "$control" \
        "$MACS_GSIZE" "$output_bw" >> "$MANIFEST"
    echo "DONE $sample: $output_bw at $(date -Is)"
}

for state in 119 137; do
    for mark in Pol2 Tag; do
        run_one "$state" "$mark" rep1 "${BAM[${state}_${mark}_1]}" ""
        run_one "$state" "$mark" rep2 "${BAM[${state}_${mark}_2]}" ""
        run_one "$state" "$mark" rep1-rep2_pooled \
            "${BAM[${state}_${mark}_1]}" "${BAM[${state}_${mark}_2]}"
    done
done

echo "All 12 MACS2 IgG-control p-value bigWigs completed."
