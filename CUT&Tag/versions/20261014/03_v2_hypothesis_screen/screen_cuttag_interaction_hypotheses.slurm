#!/bin/bash -l
#SBATCH --job-name=cuttag_hypothesis
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-hypothesis-%j.out

set -euo pipefail

PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
RESULTS="$PROJECT/results_pig_v2"
REF="$HOME/Projects/references/pig/Sscrofa11.1"
OUT="$RESULTS/hypothesis_screen_current_peaks"
GTF="$REF/genes.Ensembl_115.gtf.gz"
CHROMS="$REF/genome.chrom.sizes"
PEAKROOT="$RESULTS/03_peak_calling/04_called_peaks"

source "$HOME/apps/cutandtag/cutandtag_env.sh"

find_tool() {
    local name="$1" path
    path=$(find "$CUTTAG_ROOT/conda/envs" -mindepth 3 -maxdepth 3 \
        -type f -path "*/bin/$name" -perm -u+x -print 2>/dev/null | sort | tail -n 1)
    [[ -n "$path" ]] || { echo "Missing executable: $name" >&2; exit 1; }
    printf '%s\n' "$path"
}

BEDTOOLS=$(find_tool bedtools)

declare -A PEAKS
for caller in macs2 seacr; do
    for state in 119 137; do
        for mark in Tag Pol2; do
            for rep in 1 2; do
                key="${caller}_${state}_${mark}_R${rep}"
                if [[ "$caller" == macs2 ]]; then
                    PEAKS[$key]="$PEAKROOT/macs2/3D4_${state}_${mark}_R${rep}.macs2_peaks.narrowPeak"
                else
                    PEAKS[$key]="$PEAKROOT/seacr/3D4_${state}_${mark}_R${rep}.seacr.peaks.stringent.bed"
                fi
            done
        done
    done
done

for key in "${!PEAKS[@]}"; do
    [[ -s "${PEAKS[$key]}" ]] || { echo "Missing peak file [$key]: ${PEAKS[$key]}" >&2; exit 1; }
done
for file in "$GTF" "$CHROMS"; do
    [[ -s "$file" ]] || { echo "Missing reference file: $file" >&2; exit 1; }
done

if [[ "${1:-}" == --check ]]; then
    echo "Inputs present: MACS2 and SEACR peaks for Tag/Pol2, states 119/137, replicates 1/2."
    echo "Reference present: Ensembl 115 GTF and Sscrofa11.1 chromosome sizes."
    echo "This screen uses existing peak calls; it does not relax peak-calling thresholds."
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: bash $0 --check OR sbatch $0" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" && "$(hostname)" != mn-* ]] || {
    echo "Submit this script with sbatch." >&2
    exit 1
}

mkdir -p "$OUT/peak_sets" "$OUT/annotations" "$OUT/overlaps" "$OUT/nearest_genes"

make_bed3() {
    local src="$1" dst="$2"
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $3>$2 {print $1,$2,$3}' "$src" |
        sort -k1,1 -k2,2n -k3,3n |
        "$BEDTOOLS" merge -i - > "$dst"
}

make_shared() {
    local a="$1" b="$2" dst="$3"
    "$BEDTOOLS" intersect -a "$a" -b "$b" -wo |
        awk 'BEGIN{OFS="\t"} {s=($2>$5?$2:$5); e=($3<$6?$3:$6); if(e>s) print $1,s,e}' |
        sort -k1,1 -k2,2n -k3,3n |
        "$BEDTOOLS" merge -i - > "$dst"
}

for caller in macs2 seacr; do
    for state in 119 137; do
        for mark in Tag Pol2; do
            for rep in 1 2; do
                key="${caller}_${state}_${mark}_R${rep}"
                make_bed3 "${PEAKS[$key]}" "$OUT/peak_sets/${key}.bed"
            done
            make_shared \
                "$OUT/peak_sets/${caller}_${state}_${mark}_R1.bed" \
                "$OUT/peak_sets/${caller}_${state}_${mark}_R2.bed" \
                "$OUT/peak_sets/${caller}_${state}_${mark}_rep1-rep2_shared.bed"
        done
    done
done

# One row per annotated gene. The gene-boundary coordinates are used only for
# candidate-locus checks and nearest-gene assignment; they do not prove regulation.
gzip -dc "$GTF" |
awk -F '\t' 'BEGIN{OFS="\t"}
    /^#/ || $3!="gene" {next}
    {
        id=""; name=""; n=split($9,a,";")
        for(i=1;i<=n;i++) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", a[i])
            if(a[i] ~ /^gene_id[[:space:]]/) {
                id=a[i]; sub(/^gene_id[[:space:]]+/,"",id); gsub(/"/,"",id)
            }
            if(a[i] ~ /^gene_name[[:space:]]/) {
                name=a[i]; sub(/^gene_name[[:space:]]+/,"",name); gsub(/"/,"",name)
            }
        }
        if(name=="") name=id
        if(id!="" && name!="") print $1,$4-1,$5,name,id,$7
    }' |
sort -k1,1 -k2,2n -k3,3n > "$OUT/annotations/Ensembl115.genes.bed"

# Capture any rRNA/tRNA annotations actually present in this GTF. Tandem 45S
# rDNA arrays can be absent or collapsed in a standard assembly, so a zero count
# here is a reference limitation rather than biological evidence of absence.
gzip -dc "$GTF" |
awk -F '\t' 'BEGIN{OFS="\t"}
    /^#/ {next}
    {
        low=tolower($3 " " $9)
        if(low ~ /ribosomal_rna|gene_biotype "rrna"|transcript_biotype "rrna"|(^|[^a-z])rrna([^a-z]|$)/)
            print $1,$4-1,$5,"rRNA",$3,$7
    }' |
sort -k1,1 -k2,2n -k3,3n | "$BEDTOOLS" merge -i - -c 4 -o distinct \
    > "$OUT/annotations/Ensembl115.rRNA.bed"

gzip -dc "$GTF" |
awk -F '\t' 'BEGIN{OFS="\t"}
    /^#/ {next}
    {
        low=tolower($3 " " $9)
        if(low ~ /gene_biotype "trna"|transcript_biotype "trna"|(^|[^a-z])trna([^a-z]|$)/)
            print $1,$4-1,$5,"tRNA",$3,$7
    }' |
sort -k1,1 -k2,2n -k3,3n | "$BEDTOOLS" merge -i - -c 4 -o distinct \
    > "$OUT/annotations/Ensembl115.tRNA.bed"

gzip -dc "$GTF" |
awk -F '\t' 'BEGIN{OFS="\t"}
    /^#/ {next}
    {
        low=tolower($3 " " $9)
        if(low ~ /5s[_ -]?rrna|5s ribosomal|rn5s|rna5s/)
            print $1,$4-1,$5,"5S_rRNA",$3,$7
    }' |
sort -k1,1 -k2,2n -k3,3n | "$BEDTOOLS" merge -i - -c 4 -o distinct \
    > "$OUT/annotations/Ensembl115.5S_rRNA.bed"

cat > "$OUT/annotations/candidate_genes.txt" <<'CANDIDATES'
POLR1D
EGR2
MYD88
JAK2
HIF1A
TRAF6
TGIF2
ACTB
MDM2
POLR2H
CHUK
DVL3
CANDIDATES

awk 'BEGIN{IGNORECASE=1; OFS="\t"}
    NR==FNR {want[toupper($1)]=1; next}
    toupper($4) in want {print}' \
    "$OUT/annotations/candidate_genes.txt" \
    "$OUT/annotations/Ensembl115.genes.bed" \
    > "$OUT/annotations/candidate_genes.bed"

"$BEDTOOLS" slop -i "$OUT/annotations/candidate_genes.bed" -g "$CHROMS" -b 2000 \
    > "$OUT/annotations/candidate_genes.plusminus2kb.bed"

printf 'caller\tstate\tmark\treplicated_peak_count\trRNA_overlap_peaks\t5S_rRNA_overlap_peaks\ttRNA_overlap_peaks\n' \
    > "$OUT/feature_overlap_summary.tsv"

for caller in macs2 seacr; do
    for state in 119 137; do
        for mark in Tag Pol2; do
            set="$OUT/peak_sets/${caller}_${state}_${mark}_rep1-rep2_shared.bed"
            total=$(wc -l < "$set")
            rrna=$("$BEDTOOLS" intersect -u -a "$set" -b "$OUT/annotations/Ensembl115.rRNA.bed" | wc -l)
            r5s=$("$BEDTOOLS" intersect -u -a "$set" -b "$OUT/annotations/Ensembl115.5S_rRNA.bed" | wc -l)
            trna=$("$BEDTOOLS" intersect -u -a "$set" -b "$OUT/annotations/Ensembl115.tRNA.bed" | wc -l)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$caller" "$state" "$mark" "$total" "$rrna" "$r5s" "$trna" \
                >> "$OUT/feature_overlap_summary.tsv"

            "$BEDTOOLS" intersect -wa -wb \
                -a "$OUT/annotations/candidate_genes.plusminus2kb.bed" -b "$set" \
                > "$OUT/overlaps/${caller}_${state}_${mark}.candidate_gene_window_hits.tsv"

            "$BEDTOOLS" closest -a "$set" -b "$OUT/annotations/Ensembl115.genes.bed" -d -t first \
                > "$OUT/nearest_genes/${caller}_${state}_${mark}.nearest_gene.tsv"
            awk 'BEGIN{OFS="\t"} {print $7,$8,$9}' \
                "$OUT/nearest_genes/${caller}_${state}_${mark}.nearest_gene.tsv" |
                sort -u > "$OUT/nearest_genes/${caller}_${state}_${mark}.gene_list_for_enrichment.tsv"
        done
    done
done

# 119 Tag versus 119 Pol II is the directly testable chromatin co-occupancy
# comparison in the current experiment. Report both caller-specific results.
printf 'caller\tTag119_peaks\tPol2_119_peaks\tTag119_peaks_overlapping_Pol2\tTag119_overlap_percent\tPol2_peaks_overlapping_Tag119\tPol2_overlap_percent\n' \
    > "$OUT/119_Tag_vs_119_Pol2_cooccupancy.tsv"
for caller in macs2 seacr; do
    tag="$OUT/peak_sets/${caller}_119_Tag_rep1-rep2_shared.bed"
    pol2="$OUT/peak_sets/${caller}_119_Pol2_rep1-rep2_shared.bed"
    nt=$(wc -l < "$tag"); np=$(wc -l < "$pol2")
    nto=$("$BEDTOOLS" intersect -u -a "$tag" -b "$pol2" | wc -l)
    npo=$("$BEDTOOLS" intersect -u -a "$pol2" -b "$tag" | wc -l)
    awk -v c="$caller" -v nt="$nt" -v np="$np" -v nto="$nto" -v npo="$npo" \
        'BEGIN{OFS="\t"; printf "%s\t%d\t%d\t%d\t%.2f\t%d\t%.2f\n",c,nt,np,nto,(nt?100*nto/nt:0),npo,(np?100*npo/np:0)}' \
        >> "$OUT/119_Tag_vs_119_Pol2_cooccupancy.tsv"
    "$BEDTOOLS" jaccard -a "$tag" -b "$pol2" > "$OUT/overlaps/${caller}_119_Tag_vs_119_Pol2.jaccard.tsv"
done

{
    printf 'Ensembl115_gene_rows\t%s\n' "$(wc -l < "$OUT/annotations/Ensembl115.genes.bed")"
    printf 'candidate_genes_found\t%s\n' "$(wc -l < "$OUT/annotations/candidate_genes.bed")"
    printf 'GTF_rRNA_intervals\t%s\n' "$(wc -l < "$OUT/annotations/Ensembl115.rRNA.bed")"
    printf 'GTF_5S_rRNA_intervals\t%s\n' "$(wc -l < "$OUT/annotations/Ensembl115.5S_rRNA.bed")"
    printf 'GTF_tRNA_intervals\t%s\n' "$(wc -l < "$OUT/annotations/Ensembl115.tRNA.bed")"
    printf 'interpretation\tCUT&Tag occupancy/co-occupancy is supportive chromatin evidence and does not by itself prove physical protein-protein interaction.\n'
    printf 'rDNA_limit\tAbsence of annotated rRNA overlap is inconclusive when tandem rDNA is absent/collapsed or MAPQ filtering removes multi-mapping rDNA fragments.\n'
} > "$OUT/run_summary.tsv"

tar -czf "$RESULTS/hypothesis_screen_current_peaks.tar.gz" \
    -C "$RESULTS" "$(basename "$OUT")"

echo "Completed: $OUT"
echo "Transfer: $RESULTS/hypothesis_screen_current_peaks.tar.gz"
