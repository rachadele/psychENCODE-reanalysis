#!/bin/bash
# Merged-author arm (2026-11-09, pending Rachel's decision; docs/DE_AGREEMENT_DECISIONS.md).
# Author labels mapped to pipeline subclass labels (assets/author_to_pipeline_map.subclass.json)
# before summing, so the five IT types, Lamp5 + Lamp5 Lhx6 and Sst + Sst Chodl each form one group
# that matches one pipeline label. Builder, QC masks, aggregation, limma-voom model, PMI imputation,
# contrasts and the comparison stage are those of scripts/build_pseudobulks.sh, run_de.sh and
# run_comparison_figs.sh; only --label_map differs. Writes only under this folder.
#   1. per-study pseudobulks (bin/build_local_pseudobulks.py --label_map), aggregated by cell type
#   2. limma-voom DE, author arm only (GEMMA_DE_ANALYSIS, local pseudobulks)
#   3. comparison stage vs the existing pipeline-subclass limma-voom DE (identity label map)
#   4. per-sample F1 vs merged author labels (bin/per_cohort_f1.py --merge_author)
#   5. bin/merged_author_de_agreement.py
set -euo pipefail
REPO=/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis
R=$REPO/results_v2.1.2_limma
S=$R/merged_author_2026-11-09
MAP=$REPO/assets/author_to_pipeline_map.subclass.json
ANN=$REPO/results_v2.1.2/annotation
H5=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/psychencode.txt_author_false_process_samples_false/h5ad
MEX=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/all_homo_sapiens_samples/mex
OLDPB=$REPO/results/author_submitted_false/class/experiment_pseudobulks/gemma
PY=/home/rschwartz/anaconda3/envs/scanpyenv/bin/python
WORK=/space/scratch/rschwartz_psyrerun_work/merged_author_2026-11-09
OUT=$S/pseudobulks/author_merged
export PYTHONNOUSERSITE=1
mkdir -p $OUT/experiment_pseudobulks $OUT/ct_pseudobulks $OUT/logs $S/logs

# 1. pseudobulks
build() {
  set -euo pipefail
  c=$1; d=$OUT/experiment_pseudobulks/$c; mkdir -p $d; cd $d
  [ -s ${c}_cell_counts.tsv ] && { echo "skip $c"; return 0; }
  feat=$MEX/$c/$(ls $MEX/$c | head -1)/features.tsv.gz
  srun -p Normal -C thrd64 --mem=48G -c 2 -J merged_pb_$c \
    $PY $REPO/bin/build_local_pseudobulks.py --h5ad $H5/$c/$c.h5ad --author_labels_from_cell_id \
      --label_map $MAP --mex_features $feat --gemma_rows $OLDPB/$c/${c}_aggregated_gemma.tsv.gz \
      --mask $ANN/$c/qc/${c}_mask.tsv --nnz --study_name $c > $OUT/logs/build_$c.log 2>&1
}
export -f build; export OUT ANN H5 MEX OLDPB PY REPO MAP
cat $REPO/shortnames.txt | xargs -P 8 -I{} bash -c 'build {}'
if [ ! -s $OUT/ct_pseudobulks/cell_counts.tsv ]; then
  cd $OUT/ct_pseudobulks
  srun -p Normal -C thrd64 --mem=96G -c 2 -J merged_agg $PY $REPO/bin/aggregate_celltypes_gemma.py \
    --pseudobulk_matrices $OUT/experiment_pseudobulks/*/*_aggregated_gemma.tsv.gz \
    --metadata_files $REPO/gemma/metadata > $OUT/logs/aggregate_counts.log 2>&1
  srun -p Normal -C thrd64 --mem=96G -c 2 -J merged_agg $PY $REPO/bin/aggregate_celltypes_gemma.py \
    --pseudobulk_matrices $OUT/experiment_pseudobulks/*/*_nnz_gemma.tsv.gz \
    --input_suffix _nnz_gemma.tsv.gz --output_suffix _gene_ncells.tsv.gz \
    --metadata_files $REPO/gemma/metadata > $OUT/logs/aggregate_nnz.log 2>&1
  awk 'FNR==1 && NR!=1 {next} {print}' $OUT/experiment_pseudobulks/*/*_cell_counts.tsv > $OUT/ct_pseudobulks/cell_counts.tsv
fi

# 2. DE (author arm only: the pipeline arm has no local directory and is skipped)
mkdir -p $S/nf_run_de && cd $S/nf_run_de
nextflow run $REPO/main.nf -params-file $REPO/params.subclass.json -profile conda \
  -c $R/scripts/slurm.config --outdir $S --skip_comparison --de_method limma_voom \
  --author_ct_pseudobulks $OUT/ct_pseudobulks -w $WORK -resume > $S/logs/nf_de.log 2>&1

# 3. comparison stage, identity map on pipeline labels
$PY -I -c "import json,sys; m=json.load(open(sys.argv[1])); json.dump({v: v for v in sorted(set(m.values()))}, open(sys.argv[2], 'w'), indent=1)" \
  $MAP $S/identity_map.subclass.json
# GEMMA_COMPARISON stages both arms in one task, so identical file names (pipeline labels in both
# arms) collide; its two scripts (bin/aggregate_pairwise.py, bin/average_corr.py) run here unchanged,
# with the arguments of modules/local/aggregate_pairwise and average_de_correlations.
P=$R/author_submitted_false/subclass/limma_voom/gemma; A=$S/author_submitted_true/subclass/limma_voom/gemma
C=$S/comparison/subclass; mkdir -p $C/all_corr $C/average_corr
for d in $A/*/*/; do basename $d; done | sort -u | while read k; do
  mkdir -p $C/all_corr/$k; cd $C/all_corr/$k
  [ -s ${k}_pairwise_corr.tsv ] && continue
  $PY $REPO/bin/aggregate_pairwise.py --contrast $k --pavlab_paths $P/*/$k/*_results.tsv     --author_paths $A/*/$k/*_results.tsv --params $REPO/params.subclass.json     --label_map $S/identity_map.subclass.json > $S/logs/pairwise_$k.log 2>&1
done
cd $C/average_corr
$PY $REPO/bin/average_corr.py --corr_tables $C/all_corr/*/*_pairwise_corr.tsv --annotation_level subclass   --params $REPO/params.subclass.json --label_map $S/identity_map.subclass.json > $S/logs/average_corr.log 2>&1

# 4. per-sample F1 vs merged author labels
F=$S/per_cohort_f1; mkdir -p $F/logs
f1() {
  c=$1; [ -s $F/${c}_per_sample_f1.tsv ] && return 0
  srun -p Normal -C thrd64 --mem=32G -c 1 -J merged_f1_$c $PY $REPO/bin/per_cohort_f1.py --h5ad $H5/$c/$c.h5ad \
    --cta $ANN/$c/predicted_celltypes/${c}_subclass_cell_type.tsv --mask $ANN/$c/qc/${c}_mask.tsv \
    --label_map $MAP --merge_author --study $c --outdir $F > $F/logs/$c.log 2>&1
}
export -f f1; export F
cat $REPO/shortnames.txt | xargs -P 8 -I{} bash -c 'f1 {}'

# 5. agreement vs F1 and cells per donor, next to the current 27 / strict 18 values
$PY $REPO/bin/merged_author_de_agreement.py \
  --merged_corr $S/comparison/subclass/average_corr/combined_filtered_corr.tsv \
  --current_corr $R/comparison/subclass/average_corr/combined/combined_filtered_corr.tsv \
  --current_correlations $R/de_agreement_vs_annotation/correlations_per_cohort.tsv \
  --current_f1_cells $R/de_agreement_vs_annotation/celltype_f1_cells.tsv \
  --label_map $MAP --per_sample_dir $F --merged_cell_counts $OUT/ct_pseudobulks/cell_counts.tsv \
  --current_author_de $R/author_submitted_true/subclass/limma_voom/gemma \
  --merged_author_de $S/author_submitted_true/subclass/limma_voom/gemma \
  --meta_dir $REPO/gemma/metadata --outdir $S/de_agreement
echo MERGED AUTHOR DONE
