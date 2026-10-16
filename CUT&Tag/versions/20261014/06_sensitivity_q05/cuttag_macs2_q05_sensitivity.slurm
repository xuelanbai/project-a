#!/bin/bash -l
#SBATCH --job-name=cuttag_q05_Tag
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-q05-Tag-%j.out

set -euo pipefail
trap 'rc=$?; printf "FAILED exit=%s line=%s command=%s\n" "$rc" "$LINENO" "$BASH_COMMAND" >&2; exit "$rc"' ERR

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
RESULTS="$PROJECT/results_pig_v2"
APP="$HOME/apps/cutandtag"
MACS="$APP/conda/envs/env-eb4197d7077a33c4c7b48daeaeb02f19/bin/macs2"
BEDTOOLS="$APP/conda/envs/env-8dafdac3843d8f8edf2f9907029bda58/bin/bedtools"
GSIZE_FILE="$RESULTS/custom_metadata/macs_gsize_nonN.txt"
SCREEN="$RESULTS/hypothesis_screen_current_peaks/peak_sets"
ALIGN="$RESULTS/02_alignment"

for executable in "$MACS" "$BEDTOOLS"; do
    [[ -x "$executable" ]] || { echo "Missing executable: $executable" >&2; exit 1; }
done
[[ -s "$GSIZE_FILE" ]] || { echo "Missing genome size: $GSIZE_FILE" >&2; exit 1; }
read -r MACS_GSIZE < "$GSIZE_FILE"
[[ "$MACS_GSIZE" =~ ^[0-9]+$ && "$MACS_GSIZE" -gt 2000000000 ]] || { echo "Invalid MACS gsize" >&2; exit 1; }
[[ -d "$PROJECT/logs" && -d "$ALIGN" && -d "$SCREEN" ]] || { echo "Required results/log directory missing" >&2; exit 1; }

find_one() {
    local directory="$1" filename="$2"
    local -a matches=()
    mapfile -t matches < <(find -L "$directory" -type f -name "$filename" -print)
    (( ${#matches[@]} == 1 )) || { echo "Expected one $filename; found ${#matches[@]}" >&2; exit 1; }
    [[ -s "${matches[0]}" ]] || { echo "Empty input: ${matches[0]}" >&2; exit 1; }
    printf '%s\n' "${matches[0]}"
}

declare -A BAM
for state in 119 137; do
    BAM["${state}_IgG"]=$(find_one "$ALIGN" "3D4_${state}_IgG_R1.target.dedup.sorted.bam")
    for rep in 1 2; do
        BAM["${state}_Tag_${rep}"]=$(find_one "$ALIGN" "3D4_${state}_Tag_R${rep}.target.markdup.sorted.bam")
    done
    [[ -f "$SCREEN/macs2_${state}_Tag_rep1-rep2_shared.bed" &&
       -f "$SCREEN/macs2_${state}_Pol2_rep1-rep2_shared.bed" ]] || {
        echo "Missing strict reference peak set for state $state" >&2; exit 1;
    }
done

if [[ "${1:-}" == --check && $# -eq 1 ]]; then
    echo "Input check OK: four Tag BAMs, two matched IgG BAMs, strict peak sets, MACS2 and bedtools."
    echo "Effective genome size: $MACS_GSIZE; sensitivity setting: narrow MACS2 q=0.05, BAMPE, same control and duplicate policy."
    "$MACS" --version
    "$BEDTOOLS" --version
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: bash $0 --check OR sbatch $0" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" && "$(hostname)" != mn-* ]] || {
    echo "Use sbatch on a compute node." >&2; exit 1;
}

ROOT="$RESULTS/sensitivity_macs2_q05_Tag"
OUT="$ROOT/job-${SLURM_JOB_ID}"
[[ ! -e "$OUT" ]] || { echo "Output already exists: $OUT" >&2; exit 1; }
mkdir -p "$ROOT"
mkdir -- "$OUT"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

bed3_merged() {
    local infile="$1" outfile="$2"
    awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=3 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $3>$2 {print $1,$2,$3}' "$infile" \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | "$BEDTOOLS" merge -i - > "$outfile"
}

printf 'state\tq05_rep1\tq05_rep2\tq05_2of2_shared\tstrict_q01_2of2_shared\tq05_shared_overlaps_q01\tq05_shared_disjoint_q01\tq05_shared_overlaps_strict_Pol2\n' \
    > "$OUT/counts.tsv.partial"

for state in 119 137; do
    echo "START state=$state at $(date -Is)"
    for rep in 1 2; do
        name="3D4_${state}_Tag_R${rep}.MACS2_q05_vs_3D4_${state}_IgG_R1"
        dir="$OUT/$name"
        mkdir -- "$dir"
        "$MACS" callpeak \
            --nomodel --shift -75 --extsize 150 --keep-dup all \
            -q 0.05 -g "$MACS_GSIZE" -f BAMPE \
            -n "$name" --outdir "$dir" \
            -t "${BAM[${state}_Tag_${rep}]}" -c "${BAM[${state}_IgG]}"
        peak="$dir/${name}_peaks.narrowPeak"
        [[ -f "$peak" ]] || { echo "No MACS2 peak file: $peak" >&2; exit 1; }
        bed3_merged "$peak" "$OUT/3D4_${state}_Tag_R${rep}.q05.bed"
    done

    a="$OUT/3D4_${state}_Tag_R1.q05.bed"
    b="$OUT/3D4_${state}_Tag_R2.q05.bed"
    strict="$SCREEN/macs2_${state}_Tag_rep1-rep2_shared.bed"
    pol2="$SCREEN/macs2_${state}_Pol2_rep1-rep2_shared.bed"
    shared="$OUT/3D4_${state}_Tag_rep1-rep2.q05.2of2_shared.bed"
    "$BEDTOOLS" intersect -a "$a" -b "$b" \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | "$BEDTOOLS" merge -i - > "$shared.partial"
    mv -- "$shared.partial" "$shared"
    "$BEDTOOLS" intersect -a "$shared" -b "$strict" -v \
        > "$OUT/3D4_${state}_Tag_rep1-rep2.q05.disjoint_from_q01.bed"
    "$BEDTOOLS" intersect -a "$shared" -b "$strict" -u \
        > "$OUT/3D4_${state}_Tag_rep1-rep2.q05.overlapping_q01.bed"
    "$BEDTOOLS" intersect -a "$shared" -b "$pol2" -u \
        > "$OUT/3D4_${state}_Tag_rep1-rep2.q05.overlapping_strict_Pol2.bed"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$state" \
        "$(wc -l < "$a")" "$(wc -l < "$b")" "$(wc -l < "$shared")" \
        "$(wc -l < "$strict")" \
        "$(wc -l < "$OUT/3D4_${state}_Tag_rep1-rep2.q05.overlapping_q01.bed")" \
        "$(wc -l < "$OUT/3D4_${state}_Tag_rep1-rep2.q05.disjoint_from_q01.bed")" \
        "$(wc -l < "$OUT/3D4_${state}_Tag_rep1-rep2.q05.overlapping_strict_Pol2.bed")" \
        >> "$OUT/counts.tsv.partial"
    echo "DONE state=$state at $(date -Is)"
done

mv -- "$OUT/counts.tsv.partial" "$OUT/counts.tsv"
printf '%s\n' \
    'Sensitivity analysis only: each Tag replicate called with MACS2 q=0.05 and its same-state IgG; BAMPE, narrow peaks, same non-N gsize and duplicate policy as v2.' \
    'Positive-length overlap of biological replicates defines 2-of-2 shared BED3 bases; union of adjacent/overlapping intersections is merged.' \
    'Counts overlapping strict q=0.01 and strict Pol2 require only >=1 bp. Disjoint-from-q01 intervals do not overlap any strict shared interval.' \
    'No secondary q<=0.01 filter is applied to q=0.05 calls. More peaks alone do not show higher specificity or stronger biological evidence.' \
    'The matched-background audit from job 51892 was insufficient (0,0,0,1,2,7,0 control regions of 20 per site); evaluate IgG/duplicates before biological claims.' \
    > "$OUT/README.txt"
printf 'Completed q=0.05 sensitivity calls for 119/137 Tag biological replicates.\n' > "$OUT/COMPLETED.txt"
ARCHIVE="$ROOT/cuttag_macs2_q05_Tag-${SLURM_JOB_ID}.tar.gz"
[[ ! -e "$ARCHIVE" && ! -e "$ARCHIVE.partial" ]] || { echo "Archive already exists" >&2; exit 1; }
tar -czf "$ARCHIVE.partial" -C "$ROOT" "$(basename "$OUT")"
mv -- "$ARCHIVE.partial" "$ARCHIVE"
echo "Transfer: $ARCHIVE"
