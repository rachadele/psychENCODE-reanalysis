# Merged-author DE arm (2026-11-09)

Status: pending Rachel's decision. Figure 4 and the manuscript were not changed.

## Aim

The pipeline pools seven author subclasses into three labels (five IT types into L2/3-6 IT, LAMP5 and LAMP5 LHX6 into LAMP5, SST and SST Chodl into SST). For these types the author and pipeline pseudobulks held different cells, so we had excluded them from the agreement analysis. We tested whether merging the author types to the pipeline labels before DE removed the need for that exclusion.

## Methods

We mapped each author label to its pipeline subclass label (`assets/author_to_pipeline_map.subclass.json`) before summing counts, using `bin/build_local_pseudobulks.py --label_map`. This gave 21 author groups, each matching one pipeline label. QC masks, aggregation, the limma-voom model (~ Disorder + Age_death + PMI + Biological_Sex + X1000G_ancestry, one fit per type over all eight cohorts, no cohort term), PMI imputation and contrasts were those of `../scripts/run_de.sh`. We compared the new DE with the existing pipeline-subclass DE using the scripts of the pipeline comparison stage (`bin/aggregate_pairwise.py`, `bin/average_corr.py`) over the eight Figure 4 contrasts. F1 was the plain mean of per-sample F1 against merged author labels (`bin/per_cohort_f1.py --merge_author`) over the donors in each type's DE. Cells per donor was the median over the same donors. `bin/merged_author_de_agreement.py` computed plain and partial Spearman correlations across the 21 types. `scripts/run_merged_author.sh` runs every step (tmux session `merged_author_de`).

## Results

The 18 unpooled types reproduced the current DE, ρ, F1 and cells per donor exactly (`de_agreement/checks.json`). The merged pseudobulks equalled the sums of the current per-type columns. Mean agreement of log2FC was 0.88 ± 0.19 over 21 types and eight contrasts, compared with 0.83 ± 0.21 for the current 27 types (`de_agreement/agreement_summary.tsv`). The three merged groups agreed at 0.99 ± 0.01. Across the 21 types, agreement correlated with F1 (ρ = 0.91, p = 1.7e-8) and with cells per donor (ρ = 0.84). Given cells per donor, the partial ρ with F1 was 0.67 (p = 0.0011). Given F1, the partial ρ with cells per donor was 0.33 (p = 0.15). These values were close to those of the strict 18 one-to-one types (0.90, 0.70 and 0.25; `de_agreement/correlations_merged_vs_current.tsv`).

## Conclusion

Merging the pooled author types gave the same conclusion as the strict 18 set without excluding any type. Agreement followed F1 once cells per donor were controlled for.

`manifest.json` lists sha256 and timestamps for every input and output.
