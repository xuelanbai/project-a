#!/bin/bash -l
#SBATCH --job-name=cuttag_tss_v2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-tss-v2-%j.out

set -euo pipefail
PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
RESULTS="$PROJECT/results_pig_v2"
REF="$HOME/Projects/references/pig/Sscrofa11.1"
DT_BIN="$HOME/apps/cutandtag/conda/envs/env-0d01150023c68e439ebff3eb959d98bf/bin"
OUT="$RESULTS/custom_tracks/TSS_allGenes"
export MPLBACKEND=Agg

for tool in computeMatrix plotProfile; do
    [[ -x "$DT_BIN/$tool" ]] || { echo "Missing: $DT_BIN/$tool" >&2; exit 1; }
done
for file in "$REF/genome.fa.fai" "$REF/genes.Ensembl_115.gtf.gz"; do
    [[ -s "$file" ]] || { echo "Missing: $file" >&2; exit 1; }
done
find_bw() {
    local -a matches
    mapfile -t matches < <(find -L "$RESULTS/03_peak_calling/03_bed_to_bigwig" -type f -name "$1" -print)
    (( ${#matches[@]} == 1 )) || { echo "Expected one bigWig: $1; found ${#matches[@]}" >&2; return 1; }
    [[ -s "${matches[0]}" ]] || return 1
    printf '%s\n' "${matches[0]}"
}
declare -A BW PBW
for state in 119 137; do
    BW["${state}_IgG"]=$(find_bw "3D4_${state}_IgG_R1.bigWig")
    for mark in Pol2 Tag; do
        PBW["${state}_${mark}"]="$RESULTS/custom_tracks/macs2_ppois/3D4_${state}_${mark}_rep1-rep2_pooled_vs_3D4_${state}_IgG_rep1.MACS2_ppois_Sscrofa11.1.bigWig"
        [[ -s "${PBW[${state}_${mark}]}" ]] || { echo "Missing pooled p-value bigWig; finish the postprocessing job first: ${PBW[${state}_${mark}]}" >&2; exit 1; }
        for rep in 1 2; do
            BW["${state}_${mark}_${rep}"]=$(find_bw "3D4_${state}_${mark}_R${rep}.bigWig")
        done
    done
done
if [[ "${1:-}" == --check ]]; then
    echo "TSS inputs present: 10 CPM bigWigs, 4 pooled IgG-relative p-value bigWigs, Ensembl 115 GTF, genome FAI, deepTools."
    echo "Design: all annotated gene-boundary TSSs, strand-aware +/-3 kb, 10 bp bins; retain zero-signal genes."
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: bash $0 --check OR sbatch $0" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" && "$(hostname)" != mn-* ]] || { echo "Submit with sbatch to a compute node." >&2; exit 1; }
mkdir -p "$OUT"
export MPLCONFIGDIR="$OUT/.matplotlib"
TSS="$OUT/Ensembl115.allGenes.geneBoundaryTSS.bed"
# GTF positions are 1-based inclusive; BED below is 0-based half-open.
# A gene-boundary TSS is an annotation convention, not an experimentally measured dominant TSS.
gzip -dc "$REF/genes.Ensembl_115.gtf.gz" | awk -F '\t' -v OFS='\t' '
    FNR==NR {sizes[$1]=$2; next}
    /^#/ {next}
    $3=="gene" {
        if (!($1 in sizes) || ($7!="+" && $7!="-")) {bad=1; next}
        id=""; n=split($9, attrs, ";")
        for(i=1;i<=n;i++) {
            sub(/^[[:space:]]+/, "", attrs[i])
            if(attrs[i] ~ /^gene_id[[:space:]]/) {
                id=attrs[i]; sub(/^gene_id[[:space:]]+/, "", id); gsub(/"/, "", id)
                sub(/[[:space:]]+$/, "", id)
            }
        }
        pos=($7=="+" ? $4-1 : $5-1)
        if(id=="" || seen[id]++ || pos<0 || pos>=sizes[$1]) {bad=1; next}
        print $1,pos,pos+1,id,0,$7; count++
    }
    END {if(bad || count==0) exit 1}
' "$REF/genome.fa.fai" - > "$TSS.partial" || {
    echo "GTF gene/TSS validation failed (coordinates, IDs, strand or reference chromosome names)." >&2
    exit 1
}
mv -- "$TSS.partial" "$TSS"
echo "Gene-boundary TSSs: $(wc -l < "$TSS")"

for mark in Pol2 Tag; do
    matrix="$OUT/${mark}.allGenes.CPM.TSS_3kb.matrix.gz"
    # No peak-based selection and no --skipZeros: the TSS QC does not select enriched genes first.
    "$DT_BIN/computeMatrix" reference-point --referencePoint TSS \
        -R "$TSS" \
        -S "${BW[119_${mark}_1]}" "${BW[119_${mark}_2]}" "${BW[119_IgG]}" \
           "${BW[137_${mark}_1]}" "${BW[137_${mark}_2]}" "${BW[137_IgG]}" \
        --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
        --binSize 10 --missingDataAsZero --numberOfProcessors 2 \
        --samplesLabel "119_${mark}_rep1" "119_${mark}_rep2" "119_IgG_rep1" \
                       "137_${mark}_rep1" "137_${mark}_rep2" "137_IgG_rep1" \
        --outFileSortedRegions "$OUT/${mark}.allGenes.TSS.included.bed" \
        --outFileName "$matrix"
    "$DT_BIN/plotProfile" -m "$matrix" --perGroup --averageType mean \
        --yMin 0 --yAxisLabel CPM --refPointLabel TSS \
        --plotTitle "$mark: all annotated gene-boundary TSSs" \
        --outFileNameData "$OUT/${mark}.allGenes.CPM.TSS.profile.txt" \
        --outFileName "$OUT/${mark}.allGenes.CPM.TSS.profile.pdf"

    pooled_matrix="$OUT/${mark}.allGenes.pooled_MACS2_ppois.TSS_3kb.matrix.gz"
    "$DT_BIN/computeMatrix" reference-point --referencePoint TSS \
        -R "$TSS" -S "${PBW[119_${mark}]}" "${PBW[137_${mark}]}" \
        --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
        --binSize 10 --missingDataAsZero --numberOfProcessors 2 \
        --samplesLabel "119_${mark}_pooled_vs_119_IgG" "137_${mark}_pooled_vs_137_IgG" \
        --outFileName "$pooled_matrix"
    "$DT_BIN/plotProfile" -m "$pooled_matrix" --perGroup --averageType mean \
        --yMin 0 --yAxisLabel 'MACS2 -log10(P)' --refPointLabel TSS \
        --plotTitle "$mark: pooled replicates versus matched IgG" \
        --outFileNameData "$OUT/${mark}.allGenes.pooled_MACS2_ppois.TSS.profile.txt" \
        --outFileName "$OUT/${mark}.allGenes.pooled_MACS2_ppois.TSS.profile.pdf"
done
echo "TSS profiles completed: $OUT"
echo "Read as coverage profiles, not a standardized ATAC TSS-enrichment score or a differential-binding test."
