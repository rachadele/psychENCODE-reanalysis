# DE agreement with grid labels on all cells

## Aim

We widened the range of label F1 behind the Figure 4 DE agreement analysis by running pseudobulk DE with labels of known lower quality.

## Methods

We ran two nf-annotation-grid configurations on every cell of CMC, MultiomeBrain and SZBDMulti-Seq (the grid itself labelled 1000 cells per sample): scVI kNN with the Human_Multiple_Cortical_Areas_SMART-seq reference (`grid_knn_smartseq`) and scVI random forest with the aggregated reference (`grid_rf_aggregated`), both at cutoff 0 with 100 reference cells per type. We called the grid's own scripts with its module arguments, cached model and references, outside Nextflow so that only these two configurations ran (`scripts/grid_full_cells.sh`). Labelled h5ads and CTA tables are in `/cosmos/data/aim2-annotation-grid-results/full_cells_2026-11-09/`.

We built QC-masked pseudobulks and ran subclass limma-voom DE with the results_v2.1.2_limma code and model (`scripts/grid_de.sh`). The earlier DE pooled eight cohorts. Here every arm pooled the same three, so we rebuilt the author-label reference (`author_3cohort`) and added pipeline v2.1.2 labels (`v212_3cohort`) from the existing per-study pseudobulks. Author types map to grid labels through `assets/author_to_grid_map.subclass.json`.

`scripts/run_agreement.sh` computed Spearman rho of log2FC and moderated t between each arm and the author DE per cell type and contrast (`bin/compare_label_arms.py`), per-sample F1 against author labels averaged per cohort (`bin/per_cohort_f1.py`), and plain and partial Spearman of rho with F1, controlling for log10 median author cells per donor (`bin/grid_labels_de_agreement.py`). Each pair took the unweighted mean F1 over the cohorts contributing its contrast.

## Results

Mean rho of log2FC was 0.82 for v2.1.2 labels, 0.80 for `grid_rf_aggregated` and 0.75 for `grid_knn_smartseq`, whose macro F1 was 0.54 (`de_agreement/arm_summary.tsv`). Over the 308 pairs of all arms, rho of log2FC correlated with F1 (Spearman 0.32) and kept a positive partial correlation given cells per donor (0.17, `de_agreement/rho_vs_f1.tsv`). The kNN arm carried this association. Within the two arms with high F1, the partial correlation was negative.
