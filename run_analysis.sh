#!/usr/bin/env bash
# The whole analysis (docs/methods.md, sections 1 to 7): family, clusters, DEXSeq with stageR,
# DRIMSeq sensitivity, permutation diagnostic, events, transcript catalogue and abundance,
# gene-assigned end usage with the strong-evidence flag, summary and figures.
# Run it through the Makefile (make all), which decompresses the inputs and sets the variables below:
# GTF, FAI, GENOME, NCRNA_BED, BSF_FL_BAM, PCF_FL_BAM, BSF_POLYA_BAM, PCF_POLYA_BAM, SITES_DIR
# (<rep>.end_sites.bed), READ_IDS_DIR (<rep>.polyA_read_ids.txt), OUT, THREADS, FIGURES.
set -euo pipefail
DATA="${DATA:-data}"; OUT="${OUT:-results}"
GTF="${GTF:-$DATA/annotation/merged_fixed_6.gtf}"; FAI="${FAI:-$DATA/genome/tb927_56.fa.fai}"
BSF_FL_BAM="${BSF_FL_BAM:-$DATA/BSF.full_length.bam}"; PCF_FL_BAM="${PCF_FL_BAM:-$DATA/PCF.full_length.bam}"
SITES_DIR="${SITES_DIR:-results/sites}"; THREADS="${THREADS:-4}"
READ_IDS_DIR="${READ_IDS_DIR:-$DATA/read_ids}"; GENOME="${GENOME:-$DATA/genome/tb927_56.fa}"
BSF_POLYA_BAM="${BSF_POLYA_BAM:-$DATA/BSF.polyA.bam}"; PCF_POLYA_BAM="${PCF_POLYA_BAM:-$DATA/PCF.polyA.bam}"
NCRNA_BED="${NCRNA_BED:-annotation/ensembl116_noncoding_genes.bed}"; FIGURES="${FIGURES:-1}"
# Part 7 placement rule: mapping quality threshold (primary 20) and the sensitivity thresholds, fixed in the plan
MIN_MAPQ="${MIN_MAPQ:-20}"; MAPQ_SENSITIVITY="${MAPQ_SENSITIVITY:-1 10 30 60}"
S="$(cd "$(dirname "$0")" && pwd)/scripts"
REPS=(bs_rep_1 bs_rep_2 bs_rep_3 pc_rep_1 pc_rep_2 pc_rep_3)
SITES=(); IDS=(); for r in "${REPS[@]}"; do SITES+=("$r=$SITES_DIR/$r.end_sites.bed"); IDS+=("$r=$READ_IDS_DIR/$r.polyA_read_ids.txt"); done
mkdir -p "$OUT/docs"
python "$S/00_manifest.py" start --out "$OUT/docs/run_manifest.json" --inputs "$GTF" "$FAI" "$GENOME" "$NCRNA_BED" "$BSF_FL_BAM" "$PCF_FL_BAM" \
  "$BSF_POLYA_BAM" "$PCF_POLYA_BAM" "$SITES_DIR"/*.end_sites.bed "$READ_IDS_DIR"/*.polyA_read_ids.txt \
  --options "min_ends=10" "pad=100" "smooth=11" "peak_min=5" "extend=15" "valley=0.5" "always_merge=15" "support=5x2" \
  "primary=merge30_peak0.02" "sensitivity=merge50_peak0.05,merge20_peak0.01" "fdr=0.05" "model=DEXSeq+stageR" \
  "events=min_aligned_length=200 tolerance=60 link_distance=100 min_reads=20 min_informative=10 min_cluster=10 min_fraction=0.15 sharp_start=0 (sensitivity 5)" \
  "catalogue=start_tolerance=100 end_tolerance=60 expression_filter=5 reads in 2 replicates; abundance: edgeR QL, TMM over genes and records" \
  "assumption=the six libraries are exchangeable apart from stage; no batch recorded"

python "$S/16_windows_and_family.py" --gtf "$GTF" --fai "$FAI" --full-length-bams "$BSF_FL_BAM" "$PCF_FL_BAM" --sites "${SITES[@]}" --outdir "$OUT/family"
# Part 7 needs the secondary alignment scores of every poly(A) read: one pass over the BAMs
PLACEMENTS="$OUT/family/read_placements.tsv.gz"
python "$S/31_read_placements.py" --bams "bs=$BSF_POLYA_BAM" "pc=$PCF_POLYA_BAM" --out "$PLACEMENTS"
python "$S/17_end_clusters.py" --self-test
for mode in nanopore intergenic; do
  python "$S/17_end_clusters.py" --family "$OUT/family/family_$mode.tsv" --sites "${SITES[@]}" --outdir "$OUT/clusters_${mode}_m30p2"
  python "$S/17_end_clusters.py" --family "$OUT/family/family_$mode.tsv" --sites "${SITES[@]}" --merge 50 --peak-fraction 0.05 --outdir "$OUT/clusters_${mode}_m50p5"
  python "$S/17_end_clusters.py" --family "$OUT/family/family_$mode.tsv" --sites "${SITES[@]}" --merge 20 --peak-fraction 0.01 --outdir "$OUT/clusters_${mode}_m20p1"
  for setting in m30p2 m50p5 m20p1; do Rscript "$S/18_dexseq_stager.R" "$OUT/clusters_${mode}_$setting" 0.05 "$THREADS"; done
  # DRIMSeq on the primary clusters, observed counts, no fallback
  Rscript "$S/09_drimseq.R" "$OUT/clusters_${mode}_m30p2" 0.05 10 0 0 0
  python "$S/19_permutation_diagnostic.py" --folder "$OUT/clusters_${mode}_m30p2" --threads "$THREADS"
  python "$S/20_support_by_setting.py" --primary "primary=$OUT/clusters_${mode}_m30p2" \
    --settings "merge50_peak5=$OUT/clusters_${mode}_m50p5" "merge20_peak1=$OUT/clusters_${mode}_m20p1" \
    --drimseq "$OUT/clusters_${mode}_m30p2" --out "$OUT/docs/support_by_setting_$mode.tsv"
done
python "$S/00_manifest.py" finish --out "$OUT/docs/run_manifest.json" --outputs "$OUT"/family/*.tsv "$OUT"/family/*.bed "$PLACEMENTS" \
  "$OUT"/clusters_*/sites.tsv "$OUT"/clusters_*/site_counts.tsv "$OUT"/clusters_*/windows.tsv "$OUT"/clusters_*/dexseq_*.tsv \
  "$OUT"/clusters_*/stager.tsv "$OUT"/clusters_*_m30p2/drimseq_*.tsv "$OUT"/clusters_*_m30p2/permutation_calls.tsv "$OUT"/docs/support_by_setting_*.tsv
# Interpretation stages (Parts 4 to 7) get their own manifest, so that a rerun of them alone is recorded
python "$S/00_manifest.py" start --out "$OUT/docs/run_manifest_interpretation.json" --inputs "$GENOME" "$NCRNA_BED" "$BSF_POLYA_BAM" "$PCF_POLYA_BAM" \
  "$BSF_FL_BAM" "$PCF_FL_BAM" "$READ_IDS_DIR"/*.polyA_read_ids.txt "$OUT/family/cds_v2.bed" "$OUT/family/transcript_starts.tsv" \
  "$OUT"/clusters_*_m30p2/dexseq_site.tsv "$OUT"/clusters_*_m30p2/dexseq_gene.tsv "$OUT"/clusters_*_m30p2/sites.tsv \
  --options "events=min_aligned_length=200 tolerance=60 link_distance=100 min_reads=20 min_informative=10 min_cluster=10 min_fraction=0.15 sharp_start=0 (sensitivity 5)" \
  "catalogue=start_tolerance=100 end_tolerance=60 expression_filter=5 reads in 2 replicates; abundance: edgeR QL, TMM over genes and records" \
  "gene_assigned=aligned CDS bases>=50, mapq>=$MIN_MAPQ, no secondary alignment with AS>=primary AS, no other same-strand CDS with >=50 aligned bases; min_ends=10 support=5x2; strong=change>=30 points, replicates separated, same direction in both other settings, >=20 assigned ends; mapq sensitivity $MAPQ_SENSITIVITY with the secondary rule and 20 without it" \
  "ranking=tiers by replicate separation, size of change (30, 20 points), support in all analyses, label; within tier by size of change"
# Part 5: independent transcripts, their own family and abundance test
python "$S/22_transcript_catalogue.py" --full-length-bams "$BSF_FL_BAM" "$PCF_FL_BAM" --cds-bed "$OUT/family/cds_v2.bed" \
  --ncrna "$NCRNA_BED" --read-ids "${IDS[@]}" --outdir "$OUT/transcripts"
Rscript "$S/23_transcript_abundance.R" "$OUT/transcripts" "$OUT/clusters_nanopore_m30p2/samples.tsv"
cp "$OUT/transcripts/transcript_abundance.tsv" "$OUT/docs/"
# Part 4: events, with the catalogue and abundance results joined; a sensitivity run counts reproducible alignment starts
for mode in nanopore intergenic; do
  for ss in 0 5; do
    suffix=""; [[ "$ss" != "0" ]] && suffix="_sharpstart_sensitivity"
    python "$S/21_events.py" --dexseq "$OUT/clusters_${mode}_m30p2" --family "$OUT/family/family_$mode.tsv" \
      --transcript-starts "$OUT/family/transcript_starts.tsv" --polya-bams "$BSF_POLYA_BAM" "$PCF_POLYA_BAM" \
      --full-length-bams "$BSF_FL_BAM" "$PCF_FL_BAM" --read-ids "${IDS[@]}" --catalogue "$OUT/transcripts/catalogue.tsv" \
      --abundance "$OUT/transcripts/transcript_abundance.tsv" --ncrna "$NCRNA_BED" --genome "$GENOME" --sharp-start "$ss" \
      --out "$OUT/docs/events_$mode$suffix.tsv"
  done
done
# Part 7: gene-assigned end usage (same windows and clusters, reads assigned to the gene only) and the strong-evidence flag
for mode in nanopore intergenic; do
  for setting in m30p2 m50p5 m20p1; do
    python "$S/28_gene_assigned_counts.py" --family "$OUT/family/family_$mode.tsv" --clusters "$OUT/clusters_${mode}_$setting" \
      --cds-bed "$OUT/family/cds_v2.bed" --polya-bams "$BSF_POLYA_BAM" "$PCF_POLYA_BAM" --read-ids "${IDS[@]}" \
      --min-mapq "$MIN_MAPQ" --placements "$PLACEMENTS" --outdir "$OUT/gene_assigned_${mode}_$setting"
    Rscript "$S/18_dexseq_stager.R" "$OUT/gene_assigned_${mode}_$setting" 0.05 "$THREADS"
  done
  python "$S/20_support_by_setting.py" --primary "primary=$OUT/gene_assigned_${mode}_m30p2" \
    --settings "merge50_peak5=$OUT/gene_assigned_${mode}_m50p5" "merge20_peak1=$OUT/gene_assigned_${mode}_m20p1" \
    --out "$OUT/docs/support_by_setting_gene_assigned_$mode.tsv"
  python "$S/29_strong_end_usage.py" --primary "$OUT/gene_assigned_${mode}_m30p2" \
    --settings "merge50_peak5=$OUT/gene_assigned_${mode}_m50p5" "merge20_peak1=$OUT/gene_assigned_${mode}_m20p1" \
    --window-dexseq "$OUT/clusters_${mode}_m30p2" --events "$OUT/docs/events_$mode.tsv" --out "$OUT/docs/strong_end_usage_$mode.tsv"
done
# Part 7 sensitivity to the placement rule (nanopore family): the whole chain under other mapping quality thresholds, and without the secondary rule
RUNS=("primary=$MIN_MAPQ=yes=$OUT/docs/strong_end_usage_nanopore.tsv=$OUT/gene_assigned_nanopore_m30p2/windows.tsv")
for variant in $MAPQ_SENSITIVITY "nosec"; do
  if [[ "$variant" == "nosec" ]]; then q="$MIN_MAPQ"; sec="no"; pl=(); label="mapq${q}_no_secondary_rule"; else q="$variant"; sec="yes"; pl=(--placements "$PLACEMENTS"); label="mapq$q"; fi
  V="$OUT/mapq_sensitivity/$label"
  for setting in m30p2 m50p5 m20p1; do
    python "$S/28_gene_assigned_counts.py" --family "$OUT/family/family_nanopore.tsv" --clusters "$OUT/clusters_nanopore_$setting" \
      --cds-bed "$OUT/family/cds_v2.bed" --polya-bams "$BSF_POLYA_BAM" "$PCF_POLYA_BAM" --read-ids "${IDS[@]}" \
      --min-mapq "$q" ${pl[@]+"${pl[@]}"} --outdir "$V/gene_assigned_nanopore_$setting"
    Rscript "$S/18_dexseq_stager.R" "$V/gene_assigned_nanopore_$setting" 0.05 "$THREADS"
  done
  python "$S/29_strong_end_usage.py" --primary "$V/gene_assigned_nanopore_m30p2" \
    --settings "merge50_peak5=$V/gene_assigned_nanopore_m50p5" "merge20_peak1=$V/gene_assigned_nanopore_m20p1" \
    --window-dexseq "$OUT/clusters_nanopore_m30p2" --events "$OUT/docs/events_nanopore.tsv" --out "$V/strong_end_usage_nanopore.tsv"
  RUNS+=("$label=$q=$sec=$V/strong_end_usage_nanopore.tsv=$V/gene_assigned_nanopore_m30p2/windows.tsv")
done
python "$S/32_mapq_sensitivity.py" --runs "${RUNS[@]}" --out "$OUT/docs/mapq_sensitivity_nanopore.tsv"
# Part 6: ranking, summary and browser figures of the supported genes
for mode in nanopore intergenic; do
  python "$S/26_rank_events.py" --events "$OUT/docs/events_$mode.tsv" --dexseq-site "$OUT/clusters_${mode}_m30p2/dexseq_site.tsv" \
    --support "$OUT/docs/support_by_setting_$mode.tsv" --out "$OUT/docs/ranked_$mode.tsv"
done
python "$S/24_summary_v2.py" --results "$OUT" --out "$OUT/docs/RESULTS_v2.md"
if [[ "$FIGURES" == "1" ]]; then
  for mode in nanopore intergenic; do
    mkdir -p "$OUT/docs/figures/browser_$mode" "$OUT/docs/figures/zoom_$mode"
    # supported genes of the window analysis plus the flagged genes of the gene-assigned analysis
    for g in $( (tail -n +2 "$OUT/docs/events_$mode.tsv" | cut -f1; awk -F'\t' 'NR>1 && $2=="True"{print $1}' "$OUT/docs/strong_end_usage_$mode.tsv") | sort -u); do
      python "$S/11_browser_figure.py" --gene "$g" --bsf-bam "$BSF_POLYA_BAM" --pcf-bam "$PCF_POLYA_BAM" --cds-bed "$OUT/family/cds_v2.bed" \
        --drimseq "$OUT/clusters_${mode}_m30p2" --sites "${SITES[@]}" --ncrna "$NCRNA_BED" --catalogue "$OUT/transcripts/catalogue.tsv" \
        --out "$OUT/docs/figures/browser_$mode/$g.png" > /dev/null
      python "$S/27_zoom_figure.py" --gene "$g" --family "$OUT/family/family_$mode.tsv" --dexseq "$OUT/clusters_${mode}_m30p2" \
        --events "$OUT/docs/events_$mode.tsv" --assigned "$OUT/gene_assigned_${mode}_m30p2" --cds-bed "$OUT/family/cds_v2.bed" --min-mapq "$MIN_MAPQ" --placements "$PLACEMENTS" --bsf-bam "$BSF_POLYA_BAM" --pcf-bam "$PCF_POLYA_BAM" \
        --read-ids "${IDS[@]}" --out "$OUT/docs/figures/zoom_$mode/$g.png" > /dev/null
    done
  done
  python "$S/25_figure_index_v2.py" --events "Nanopore windows=$OUT/docs/events_nanopore.tsv=$OUT/docs/figures/browser_nanopore" \
    "Intergenic windows=$OUT/docs/events_intergenic.tsv=$OUT/docs/figures/browser_intergenic" \
    --ranked "Nanopore windows=$OUT/docs/ranked_nanopore.tsv" "Intergenic windows=$OUT/docs/ranked_intergenic.tsv" \
    --zoom "Nanopore windows=$OUT/docs/figures/zoom_nanopore" "Intergenic windows=$OUT/docs/figures/zoom_intergenic" --out "$OUT/docs/figures/index.html"
  python "$S/30_strong_index.py" --strong "Nanopore windows (primary)=$OUT/docs/strong_end_usage_nanopore.tsv=$OUT/docs/figures/zoom_nanopore=$OUT/docs/figures/browser_nanopore" \
    "Intergenic windows (sensitivity)=$OUT/docs/strong_end_usage_intergenic.tsv=$OUT/docs/figures/zoom_intergenic=$OUT/docs/figures/browser_intergenic" \
    --out "$OUT/docs/figures/index_strong.html"
fi
python "$S/00_manifest.py" finish --out "$OUT/docs/run_manifest_interpretation.json" --outputs "$OUT"/docs/events_*.tsv "$OUT"/docs/ranked_*.tsv \
  "$OUT"/docs/support_by_setting_*.tsv "$OUT"/docs/strong_end_usage_*.tsv "$OUT"/gene_assigned_*/*.tsv "$OUT/docs/mapq_sensitivity_nanopore.tsv" "$OUT"/mapq_sensitivity/*/*.tsv "$OUT"/mapq_sensitivity/*/gene_assigned_nanopore_*/*.tsv "$OUT/docs/transcript_abundance.tsv" "$OUT"/transcripts/*.tsv "$OUT/docs/RESULTS_v2.md" "$OUT/docs/figures/index.html" "$OUT/docs/figures/index_strong.html"
echo "done: $OUT"
