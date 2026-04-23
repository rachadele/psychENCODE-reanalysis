#!/bin/bash
set -euo pipefail

cd /space/grp/rschwartz/rschwartz/psychENCODE-reanalysis

F1_PATH="/space/grp/rschwartz/rschwartz/evaluation_summary.nf/2024-07-01/homo_sapiens_main_branch/100/dataset_id/SCT/gap_false/grant_summary/outliers_removed/whole_cortex/cutoff_0.0/label_metrics_stats_per_label.tsv"

echo "Running comparison-only: class level"
nextflow run main.nf \
    -params-file params.class.json \
    -profile conda \
    --outdir ./results \
    --f1_path "${F1_PATH}" \
    --skip_de_analysis \
    --pavlab_deseq_results "results/author_submitted_false/class/DESeq2/gemma/**tsv" \
    --author_label_deseq_results "results/author_submitted_true/class/DESeq2/gemma/**tsv" \
    -resume

echo "Running comparison-only: subclass level"
nextflow run main.nf \
    -params-file params.subclass.json \
    -profile conda \
    --outdir ./results \
    --f1_path "${F1_PATH}" \
    --skip_de_analysis \
    --pavlab_deseq_results "results/author_submitted_false/subclass/DESeq2/gemma/**tsv" \
    --author_label_deseq_results "results/author_submitted_true/subclass/DESeq2/gemma/**tsv" \
    -resume
