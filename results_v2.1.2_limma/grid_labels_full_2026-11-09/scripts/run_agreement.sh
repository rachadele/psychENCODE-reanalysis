#!/bin/bash
# Step 3. Agreement of each arm's DE with the three-cohort author-label DE (bin/compare_label_arms.py),
# per-sample F1 of each arm's labels against author labels (bin/per_cohort_f1.py), then
# bin/grid_labels_de_agreement.py (mean rho per arm, rho vs F1 plain and partial).
set -euo pipefail
REPO=/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis
R=$REPO/results_v2.1.2_limma
S=$R/grid_labels_full_2026-11-09
G=/cosmos/data/aim2-annotation-grid-results/full_cells_2026-11-09
ANN=$REPO/results_v2.1.2/annotation
H5=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/psychencode.txt_author_false_process_samples_false/h5ad
PY=/home/rschwartz/anaconda3/envs/scanpyenv/bin/python
COH="CMC MultiomeBrain SZBDMulti-Seq"
export PYTHONNOUSERSITE=1
declare -A CTA=([grid_knn_smartseq]="$G/cta/scvi_knn__Human_Multiple_Cortical_Areas_SMART-seq__cutoff0__gapfalse/{c}_subclass_cell_type.tsv"
                [grid_rf_aggregated]="$G/cta/scvi_rf__aggregated__cutoff0__gapfalse/{c}_subclass_cell_type.tsv"
                [v212_3cohort]="$ANN/{c}/predicted_celltypes/{c}_subclass_cell_type.tsv")
declare -A MAP=([grid_knn_smartseq]=$REPO/assets/author_to_grid_map.subclass.json
                [grid_rf_aggregated]=$REPO/assets/author_to_grid_map.subclass.json
                [v212_3cohort]=$REPO/assets/author_to_pipeline_map.subclass.json)
mkdir -p $S/comparison $S/per_cohort_f1/logs
AUTH=$S/de/author_3cohort/author_submitted_true/subclass/limma_voom/gemma
armspec=()
for arm in v212_3cohort grid_rf_aggregated grid_knn_smartseq; do
  $PY $REPO/bin/compare_label_arms.py --pavlab_dir $S/de/$arm/author_submitted_false/subclass/limma_voom/gemma \
    --author_dir $AUTH --label_map ${MAP[$arm]} --label $arm --out $S/comparison/label_arms_$arm.tsv
  o=$S/per_cohort_f1/$arm; mkdir -p $o
  for c in $COH; do
    [ -s $o/${c}_per_cohort_f1.tsv ] && continue
    srun -p Normal -C thrd64 --mem=48G -c 1 -J gridf1_$c $PY $REPO/bin/per_cohort_f1.py --h5ad $H5/$c/$c.h5ad \
      --cta ${CTA[$arm]//\{c\}/$c} --mask $ANN/$c/qc/${c}_mask.tsv --label_map ${MAP[$arm]} \
      --study $c --outdir $o > $S/per_cohort_f1/logs/${arm}_$c.log 2>&1 &
  done
  wait
  armspec+=(--arm "$arm=$S/comparison/label_arms_$arm.tsv=$o")
done
$PY $REPO/bin/grid_labels_de_agreement.py "${armspec[@]}" \
  --author_cell_counts $S/pseudobulks/author_3cohort/ct_pseudobulks/cell_counts.tsv \
  --meta_dir $REPO/gemma/metadata --studies $COH \
  --data_version "grid_labels_full_2026-11-09: nf-annotation-grid scVI configs on all cells + v2.1.2 and author labels; QC-masked local pseudobulks, limma-voom, CMC+MultiomeBrain+SZBDMulti-Seq" \
  --extra_inputs $REPO/assets/author_to_grid_map.subclass.json $REPO/assets/author_to_pipeline_map.subclass.json \
    $(for c in $COH; do echo $G/cta/*/${c}_subclass_cell_type.tsv $ANN/$c/qc/${c}_mask.tsv; done) \
  --outdir $S/de_agreement
echo AGREEMENT DONE
