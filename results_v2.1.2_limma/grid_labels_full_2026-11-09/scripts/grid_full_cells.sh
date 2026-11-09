#!/bin/bash
# Step 1. Two nf-annotation-grid configs on every cell of CMC, MultiomeBrain and SZBDMulti-Seq
# (no per-sample cap; the grid's 1000-cell cap is subsample_n in nf-cohort-curation, upstream of
# the grid). The grid's own bin scripts run with the module arguments of nf-annotation-grid
# (ANNOTATE_QUERY, SCVI_PREDICT, CLASSIFY_QUERY), its cached scVI model and its cached
# 100-cells-per-type references, outside Nextflow so that only these two configs run:
#   scvi_knn__Human_Multiple_Cortical_Areas_SMART-seq__cutoff0__gapfalse  (worse labels)
#   scvi_rf__aggregated__cutoff0__gapfalse                                (control)
# Input: the psychENCODE source h5ads the DE uses (batch key sample_id, equal to the grid's
# 'sample' column). Then writes one CTA table per cohort and config (bin/grid_cta_from_classified.py).
# Skips any output that exists; never writes to the existing grid results.
set -euo pipefail
C=${1:-"CMC MultiomeBrain SZBDMulti-Seq"}
GRID=/space/grp/rschwartz/rschwartz/aim2-reannotation/nf-annotation-grid
REPO=/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis
G=/cosmos/data/aim2-annotation-grid-results/full_cells_2026-11-09
H5=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/psychencode.txt_author_false_process_samples_false/h5ad
MODEL=/cosmos/data/aim2-annotation-grid-cache/scvi_model/homo_sapiens/2024-07-01/scvi-homo_sapiens-2024-07-01
REFS=/cosmos/data/aim2-annotation-grid-cache/refs/homo_sapiens/2024-07-01/brain/dataset_id/sub_100/refs
RELABEL=/space/grp/rschwartz/rschwartz/annotation-benchmark/meta/census_map_human.tsv
GENES=/space/grp/rschwartz/rschwartz/annotation-benchmark/assets/gemma_genes.tsv
KEYS="subclass class family global"
PY=/home/rschwartz/anaconda3/envs/scanpyenv/bin/python
export PYTHONNOUSERSITE=1
mkdir -p $G/logs

mem() { case $1 in SZBDMulti-Seq) echo 300G;; CMC) echo 200G;; *) echo 96G;; esac; }

one() {
  set -euo pipefail
  c=$1; W=$G/work/$c; mkdir -p $W; cd $W
  m=$(mem $c)
  if [ ! -s $W/${c}_processed.h5ad ]; then
    srun -p Normal -C thrd64 --mem=$m -c 20 -J gridfull_annot_$c \
      $PY $GRID/bin/annotate_query.py --query_path $H5/$c/$c.h5ad --model_path $MODEL \
        --batch_key sample_id --gene_mapping $GENES --nmads 5 --seed 42 > $G/logs/annotate_$c.log 2>&1
  fi
  for ref in Human_Multiple_Cortical_Areas_SMART-seq aggregated; do
    mkdir -p $W/predict_$ref; cd $W/predict_$ref
    if [ ! -s probs/${c}_${ref}.rf.prob.df.tsv.gz ] || [ ! -s probs/${c}_${ref}.knn.prob.df.tsv.gz ]; then
      srun -p Normal -C thrd64 --mem=$m -c 8 -J gridfull_pred_$c \
        $PY $GRID/bin/predict_scvi.py --query_path $W/${c}_processed.h5ad --ref_path $REFS/$ref.h5ad \
          --ref_keys $KEYS --n_neighbors 15 > $G/logs/predict_${c}_$ref.log 2>&1
    fi
  done
  for spec in "scvi_knn Human_Multiple_Cortical_Areas_SMART-seq knn" "scvi_rf aggregated rf"; do
    set -- $spec; tag=$1; ref=$2; clf=$3
    cfg=${tag}__${ref}__cutoff0__gapfalse
    O=$G/classified/$c/$cfg; mkdir -p $O; cd $O
    out=${c}_${cfg}_annotated.h5ad
    if [ ! -s $out ]; then
      srun -p Normal -C thrd64 --mem=$m -c 2 -J gridfull_cls_$c \
        $PY $GRID/bin/classify_query.py --raw_h5ad $W/${c}_raw.h5ad \
          --probs $W/predict_$ref/probs/${c}_${ref}.$clf.prob.df.tsv.gz --ref_keys $KEYS \
          --mapping_file $RELABEL --cutoff 0 --method $tag --config_id $cfg --ref_name $ref \
          --out $out > $G/logs/classify_${c}_$cfg.log 2>&1
    fi
    mkdir -p $G/cta/$cfg
    cta=$G/cta/$cfg/${c}_subclass_cell_type.tsv
    if [ ! -s $cta ]; then
      srun -p Normal -C thrd64 --mem=64G -c 1 -J gridfull_cta_$c \
        $PY $REPO/bin/grid_cta_from_classified.py --classified $O/$out --source_h5ad $H5/$c/$c.h5ad \
          --out $cta > $G/logs/cta_${c}_$cfg.log 2>&1
    fi
  done
  echo "done $c"
}
export -f one mem; export GRID REPO G H5 MODEL REFS RELABEL GENES KEYS PY
echo $C | tr ' ' '\n' | xargs -P 3 -I{} bash -c 'one {}'
echo GRID FULL CELLS DONE
