#!/bin/bash
# Step 2. Subclass limma-voom DE over CMC, MultiomeBrain and SZBDMulti-Seq for four label arms,
# with the pseudobulk builder, QC masks, filters and model of results_v2.1.2_limma
# (scripts/build_pseudobulks.sh, scripts/run_de.sh):
#   grid_knn_smartseq  grid labels, scvi_knn, Human_Multiple_Cortical_Areas_SMART-seq ref, cutoff 0
#   grid_rf_aggregated grid labels, scvi_rf, aggregated ref, cutoff 0
#   v212_3cohort       pipeline v2.1.2 labels (existing per-study pseudobulks, re-aggregated)
#   author_3cohort     author labels (existing per-study pseudobulks, re-aggregated)
# The existing DE pools eight cohorts; these arms pool three, so the author reference is rebuilt
# on the same three cohorts. Nothing under results_v2.1.2_limma outside this folder is written.
set -euo pipefail
REPO=/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis
R=$REPO/results_v2.1.2_limma
S=$R/grid_labels_full_2026-11-09
G=/cosmos/data/aim2-annotation-grid-results/full_cells_2026-11-09
ANN=$REPO/results_v2.1.2/annotation
H5=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/psychencode.txt_author_false_process_samples_false/h5ad
MEX=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/all_homo_sapiens_samples/mex
OLDPB=$REPO/results/author_submitted_false/class/experiment_pseudobulks/gemma
PY=/home/rschwartz/anaconda3/envs/scanpyenv/bin/python
COH="CMC MultiomeBrain SZBDMulti-Seq"
export PYTHONNOUSERSITE=1
declare -A CFG=([grid_knn_smartseq]=scvi_knn__Human_Multiple_Cortical_Areas_SMART-seq__cutoff0__gapfalse
                [grid_rf_aggregated]=scvi_rf__aggregated__cutoff0__gapfalse)

build() {  # $1 arm, $2 cohort: per-study pseudobulks from grid CTA, as build_pseudobulks.sh
  set -euo pipefail
  arm=$1; c=$2; d=$S/pseudobulks/$arm/experiment_pseudobulks/$c; mkdir -p $d $S/pseudobulks/$arm/logs; cd $d
  [ -s ${c}_aggregated_gemma.tsv.gz ] && { echo "skip $arm $c"; return 0; }
  feat=$MEX/$c/$(ls $MEX/$c | head -1)/features.tsv.gz
  srun -p Normal -C thrd64 --mem=48G -c 2 -J gridde_pb_$c \
    $PY $REPO/bin/build_local_pseudobulks.py --h5ad $H5/$c/$c.h5ad --cta $G/cta/$CFGID/${c}_subclass_cell_type.tsv \
      --mex_features $feat --gemma_rows $OLDPB/$c/${c}_aggregated_gemma.tsv.gz \
      --mask $ANN/$c/qc/${c}_mask.tsv --nnz --study_name $c > $S/pseudobulks/$arm/logs/build_$c.log 2>&1
}
export -f build; export REPO S G ANN H5 MEX OLDPB PY

aggregate() {  # $1 arm, $2 dir holding <cohort>/<cohort>_{aggregated,nnz}_gemma.tsv.gz and _cell_counts.tsv
  arm=$1; src=$2; o=$S/pseudobulks/$arm/ct_pseudobulks
  [ -s $o/cell_counts.tsv ] && { echo "skip aggregate $arm"; return 0; }
  mkdir -p $o $S/pseudobulks/$arm/logs; cd $o
  agg=(); nnz=(); cc=()
  for c in $COH; do agg+=($src/$c/${c}_aggregated_gemma.tsv.gz); nnz+=($src/$c/${c}_nnz_gemma.tsv.gz); cc+=($src/$c/${c}_cell_counts.tsv); done
  srun -p Normal -C thrd64 --mem=64G -c 2 -J gridde_agg_$arm $PY $REPO/bin/aggregate_celltypes_gemma.py \
    --pseudobulk_matrices ${agg[@]} --metadata_files $REPO/gemma/metadata > $S/pseudobulks/$arm/logs/aggregate_counts.log 2>&1
  srun -p Normal -C thrd64 --mem=64G -c 2 -J gridde_agg_$arm $PY $REPO/bin/aggregate_celltypes_gemma.py \
    --pseudobulk_matrices ${nnz[@]} --input_suffix _nnz_gemma.tsv.gz --output_suffix _gene_ncells.tsv.gz \
    --metadata_files $REPO/gemma/metadata > $S/pseudobulks/$arm/logs/aggregate_nnz.log 2>&1
  awk 'FNR==1 && NR!=1 {next} {print}' ${cc[@]} > $o/cell_counts.tsv
}

for arm in grid_knn_smartseq grid_rf_aggregated; do
  export CFGID=${CFG[$arm]}
  for c in $COH; do echo "$arm $c"; done | xargs -P 3 -L 1 bash -c 'build $0 $1'
  aggregate $arm $S/pseudobulks/$arm/experiment_pseudobulks
done
aggregate v212_3cohort $R/pseudobulks/pavlab_subclass/experiment_pseudobulks
aggregate author_3cohort $R/pseudobulks/author/experiment_pseudobulks

# DE through the pipeline (GEMMA_DE_ANALYSIS, local cell-type pseudobulks), one arm per run
mkdir -p $S/nf_run && cd $S/nf_run
for arm in author_3cohort v212_3cohort grid_rf_aggregated grid_knn_smartseq; do
  if [ $arm = author_3cohort ]; then in="--author_ct_pseudobulks $S/pseudobulks/$arm/ct_pseudobulks"
  else in="--pavlab_ct_pseudobulks $S/pseudobulks/$arm/ct_pseudobulks"; fi
  nextflow run $REPO/main.nf -params-file $REPO/params.subclass.json -profile conda \
    -c $R/scripts/slurm.config --outdir $S/de/$arm --skip_comparison --de_method limma_voom $in \
    -w /space/scratch/rschwartz_psyrerun_work/grid_labels_full_2026-11-09 -resume > $S/logs/nf_de_$arm.log 2>&1
done
echo "GRID DE DONE"
