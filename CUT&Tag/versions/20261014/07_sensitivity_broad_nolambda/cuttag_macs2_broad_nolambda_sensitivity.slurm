#!/bin/bash -l
#SBATCH --job-name=cuttag_macs2_sens
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-macs2-sens-%j.out

set -euo pipefail
trap 'rc=$?; printf "FAILED exit=%s line=%s command=%s\n" "$rc" "$LINENO" "$BASH_COMMAND" >&2; exit "$rc"' ERR

CHECK=false
if [[ "${1:-}" == --check ]]; then
    CHECK=true
    MODE="${2:-}"
    [[ $# -eq 2 ]] || { echo "Usage: bash $0 --check broad|nolambda" >&2; exit 2; }
else
    MODE="${1:-}"
    [[ $# -eq 1 ]] || { echo "Usage: sbatch $0 broad|nolambda" >&2; exit 2; }
fi
case "$MODE" in
    broad)    EXPERIMENT=broad_local_q05_broadq10; PEAK_SUFFIX=broadPeak ;;
    nolambda) EXPERIMENT=narrow_nolambda_q05; PEAK_SUFFIX=narrowPeak ;;
    *) echo "Unknown mode: $MODE (use broad or nolambda)" >&2; exit 2 ;;
esac

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
RESULTS="$PROJECT/results_pig_v2"
APP="$HOME/apps/cutandtag"
MACS="$APP/conda/envs/env-eb4197d7077a33c4c7b48daeaeb02f19/bin/macs2"
BEDTOOLS="$APP/conda/envs/env-8dafdac3843d8f8edf2f9907029bda58/bin/bedtools"
GSIZE_FILE="$RESULTS/custom_metadata/macs_gsize_nonN.txt"
ALIGN="$RESULTS/02_alignment"
SCREEN="$RESULTS/hypothesis_screen_current_peaks/peak_sets"
ANNOT="$RESULTS/hypothesis_screen_current_peaks/annotations"
BASE_Q05="$RESULTS/sensitivity_macs2_q05_Tag/job-51897"

for executable in "$MACS" "$BEDTOOLS"; do
    [[ -x "$executable" ]] || { echo "Missing executable: $executable" >&2; exit 1; }
done
[[ -s "$GSIZE_FILE" ]] || { echo "Missing genome size: $GSIZE_FILE" >&2; exit 1; }
read -r MACS_GSIZE < "$GSIZE_FILE"
[[ "$MACS_GSIZE" =~ ^[0-9]+$ && "$MACS_GSIZE" -gt 2000000000 ]] || {
    echo "Invalid MACS gsize: $MACS_GSIZE" >&2; exit 1;
}
[[ -d "$PROJECT/logs" && -d "$ALIGN" && -d "$SCREEN" && -d "$ANNOT" ]] || {
    echo "Required project directories are missing" >&2; exit 1;
}
for filename in candidate_genes.plusminus2kb.bed Ensembl115.rRNA.bed Ensembl115.5S_rRNA.bed Ensembl115.tRNA.bed; do
    [[ -f "$ANNOT/$filename" ]] || { echo "Missing annotation: $ANNOT/$filename" >&2; exit 1; }
done

find_one() {
    local directory="$1" filename="$2"
    local -a matches=()
    mapfile -t matches < <(find -L "$directory" -type f -name "$filename" -print)
    (( ${#matches[@]} == 1 )) || {
        echo "Expected one $filename, found ${#matches[@]}" >&2; exit 1;
    }
    [[ -s "${matches[0]}" ]] || { echo "Empty input: ${matches[0]}" >&2; exit 1; }
    printf '%s\n' "${matches[0]}"
}

declare -A BAM
for state in 119 137; do
    BAM["${state}_IgG"]=$(find_one "$ALIGN" "3D4_${state}_IgG_R1.target.dedup.sorted.bam")
    for rep in 1 2; do
        BAM["${state}_Tag_${rep}"]=$(find_one "$ALIGN" "3D4_${state}_Tag_R${rep}.target.markdup.sorted.bam")
    done
    for input_file in \
        "$SCREEN/macs2_${state}_Tag_rep1-rep2_shared.bed" \
        "$SCREEN/macs2_${state}_Pol2_rep1-rep2_shared.bed" \
        "$BASE_Q05/3D4_${state}_Tag_rep1-rep2.q05.2of2_shared.bed"; do
        [[ -s "$input_file" ]] || { echo "Missing comparison input: $input_file" >&2; exit 1; }
    done
done

if $CHECK; then
    echo "Read-only input check passed for mode=$MODE ($EXPERIMENT)"
    echo "Four Tag BAMs, two same-state IgG BAMs, non-N pig gsize=$MACS_GSIZE"
    echo "Strict q=0.01 and local-lambda narrow q=0.05 comparator sets are present."
    "$MACS" --version
    "$BEDTOOLS" --version
    exit 0
fi
[[ -n "${SLURM_JOB_ID:-}" && "$(hostname)" != mn-* ]] || {
    echo "Submit with sbatch; do not calculate on the login node" >&2; exit 1;
}

ROOT="$RESULTS/sensitivity_macs2_broad_nolambda"
OUT="$ROOT/${EXPERIMENT}-job-${SLURM_JOB_ID}"
ARCHIVE="$ROOT/${EXPERIMENT}-job-${SLURM_JOB_ID}.tar.gz"
[[ ! -e "$OUT" && ! -e "$ARCHIVE" && ! -e "$ARCHIVE.partial" ]] || {
    echo "Output path already exists; refusing to overwrite" >&2; exit 1;
}
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
line_count() { awk 'END{print NR+0}' "$1"; }
bp_total() { awk '{n+=$3-$2} END{printf "%.0f\n",n+0}' "$1"; }

printf 'mode\tstate\trep1_peaks\trep2_peaks\trep1_supported\trep2_supported\tshared_any_intervals\tshared_any_bp\tshared_reciprocal25_intervals\tshared_overlaps_strict_q01\tshared_overlaps_local_q05\tshared_overlaps_strict_Pol2\tshared_overlaps_candidate_gene_2kb\tshared_overlaps_rRNA\tshared_overlaps_5S_rRNA\tshared_overlaps_tRNA\n' \
    > "$OUT/counts.tsv.partial"

for state in 119 137; do
    echo "START mode=$MODE state=$state at $(date -Is)"
    for rep in 1 2; do
        name="3D4_${state}_Tag_R${rep}.MACS2_${EXPERIMENT}_vs_3D4_${state}_IgG_R1"
        dir="$OUT/$name"
        mkdir -- "$dir"
        args=(--nomodel --shift -75 --extsize 150 --keep-dup all -q 0.05
              -g "$MACS_GSIZE" -f BAMPE -n "$name" --outdir "$dir")
        if [[ "$MODE" == broad ]]; then
            args+=(--broad --broad-cutoff 0.10)
        else
            args+=(--nolambda)
        fi
        "$MACS" callpeak "${args[@]}" \
            -t "${BAM[${state}_Tag_${rep}]}" -c "${BAM[${state}_IgG]}"
        peak="$dir/${name}_peaks.${PEAK_SUFFIX}"
        [[ -f "$peak" ]] || { echo "Missing MACS2 peak output: $peak" >&2; exit 1; }
        bed3_merged "$peak" "$OUT/3D4_${state}_Tag_R${rep}.bed"
    done

    a="$OUT/3D4_${state}_Tag_R1.bed"
    b="$OUT/3D4_${state}_Tag_R2.bed"
    shared="$OUT/3D4_${state}_Tag_rep1-rep2.shared_any.bed"
    shared25="$OUT/3D4_${state}_Tag_rep1-rep2.shared_reciprocal25.bed"
    "$BEDTOOLS" intersect -a "$a" -b "$b" \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | "$BEDTOOLS" merge -i - > "$shared.partial"
    mv -- "$shared.partial" "$shared"
    "$BEDTOOLS" intersect -a "$a" -b "$b" -f 0.25 -F 0.25 \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | "$BEDTOOLS" merge -i - > "$shared25.partial"
    mv -- "$shared25.partial" "$shared25"

    strict="$SCREEN/macs2_${state}_Tag_rep1-rep2_shared.bed"
    local_q05="$BASE_Q05/3D4_${state}_Tag_rep1-rep2.q05.2of2_shared.bed"
    pol2="$SCREEN/macs2_${state}_Pol2_rep1-rep2_shared.bed"
    for pair in \
        "strict_q01|$strict" \
        "local_q05|$local_q05" \
        "strict_Pol2|$pol2" \
        "candidate_gene_2kb|$ANNOT/candidate_genes.plusminus2kb.bed" \
        "rRNA|$ANNOT/Ensembl115.rRNA.bed" \
        "5S_rRNA|$ANNOT/Ensembl115.5S_rRNA.bed" \
        "tRNA|$ANNOT/Ensembl115.tRNA.bed"; do
        label="${pair%%|*}"
        target="${pair#*|}"
        "$BEDTOOLS" intersect -a "$shared" -b "$target" -u \
            > "$OUT/3D4_${state}_Tag.shared_any.overlaps_${label}.bed"
        if [[ "$label" == candidate_gene_2kb || "$label" == rRNA || "$label" == 5S_rRNA || "$label" == tRNA ]]; then
            "$BEDTOOLS" intersect -a "$shared" -b "$target" -wa -wb \
                > "$OUT/3D4_${state}_Tag.shared_any.${label}_hits.tsv"
        fi
    done

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$EXPERIMENT" "$state" \
        "$(line_count "$a")" "$(line_count "$b")" \
        "$("$BEDTOOLS" intersect -a "$a" -b "$b" -u | awk 'END{print NR+0}')" \
        "$("$BEDTOOLS" intersect -a "$b" -b "$a" -u | awk 'END{print NR+0}')" \
        "$(line_count "$shared")" "$(bp_total "$shared")" "$(line_count "$shared25")" \
        "$(line_count "$OUT/3D4_${state}_Tag.shared_any.overlaps_strict_q01.bed")" \
        "$(line_count "$OUT/3D4_${state}_Tag.shared_any.overlaps_local_q05.bed")" \
        "$(line_count "$OUT/3D4_${state}_Tag.shared_any.overlaps_strict_Pol2.bed")" \
        "$(line_count "$OUT/3D4_${state}_Tag.shared_any.overlaps_candidate_gene_2kb.bed")" \
        "$(line_count "$OUT/3D4_${state}_Tag.shared_any.overlaps_rRNA.bed")" \
        "$(line_count "$OUT/3D4_${state}_Tag.shared_any.overlaps_5S_rRNA.bed")" \
        "$(line_count "$OUT/3D4_${state}_Tag.shared_any.overlaps_tRNA.bed")" \
        >> "$OUT/counts.tsv.partial"
    echo "DONE mode=$MODE state=$state at $(date -Is)"
done

mv -- "$OUT/counts.tsv.partial" "$OUT/counts.tsv"
printf '%s\n' \
    "MODE=$MODE; MACS2=$("$MACS" --version 2>&1); pig_nonN_gsize=$MACS_GSIZE" \
    'This is exploratory parameter sensitivity only. Main results remain MACS2 narrow q<=0.01 with local lambda and 2-of-2 replicates.' \
    'The broad mode uses q<=0.05 stronger regions and broad cutoff q<=0.10. Its broadPeak outer intervals are not equivalent to narrowPeak; use gappedPeak for nested detail.' \
    'The nolambda mode keeps same-state IgG as control but disables local-lambda bias correction. Extra peaks cannot be taken as validated binding, especially at accessible/repetitive loci.' \
    'The existing --nomodel/--shift/--extsize flags are retained verbatim from job 51897 for a one-variable comparison. BAMPE uses observed paired fragments; these flags are not an additional relaxation.' \
    'Each state is called on individual biological replicates. shared_any is their positive-length intersection. shared_reciprocal25 additionally requires >=25% of each input peak to overlap; it is a secondary robustness measure.' \
    'The raw Tag BAMs keep marked duplicates, while the matched IgG BAMs are deduplicated, as in v2. The mismatch remains a key limitation. No new BAM or bigWig is created.' \
    'Annotation overlaps are descriptive and limited by the reference annotation; absence of rRNA/tRNA overlap does not rule out repetitive rDNA or unannotated Pol III targets.' \
    > "$OUT/README.txt"
printf 'Completed mode=%s job=%s\n' "$MODE" "$SLURM_JOB_ID" > "$OUT/COMPLETED.txt"
tar -czf "$ARCHIVE.partial" -C "$ROOT" "$(basename "$OUT")"
mv -- "$ARCHIVE.partial" "$ARCHIVE"
echo "Transfer: $ARCHIVE"
