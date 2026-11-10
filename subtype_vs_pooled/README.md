# Author subtypes vs pooled author types in DE

We asked whether the author subtypes that our pipeline pools (five IT types into L2/3-6 IT, Lamp5 with Lamp5_Lhx6, Sst with Sst_Chodl) gave the same DE as their pooled group. We reused existing limma-voom results: 27 unmerged author types and 21 merged author types (`results_v2.1.2_limma/merged_author_2026-11-09`). For the eight Figure 4 contrasts, `subtype_vs_pooled.py` computed Spearman ρ of log2FC and of t over shared genes between each subtype and its merged group, between subtypes of a group, between each subtype and named groups of author types outside its pooled group (below), and between the merged group and our label. We averaged pairs within each contrast and report mean ± SD over contrasts (`summary.tsv`; all pairs in `pairwise_rho.tsv`). Cells per donor was the median over DE donors (`cells_per_donor.tsv`). The merged-vs-ours values matched `de_agreement/celltype_merged.tsv`.

`manifest.json` holds sha256 and timestamps.

## Comparator groups

IT subtypes were compared with GABAergic types, deep-layer non-IT types, glia, and vascular and immune types; Lamp5 and Sst subtypes with IT types, GABAergic types outside their pooled group, deep-layer non-IT types, glia, and vascular and immune types. IT types: L2.3_IT, L4_IT, L5_IT, L6_IT, L6_IT_Car3. GABAergic types: Lamp5, Lamp5_Lhx6, Pax6, Pvalb, Sncg, Sst, Sst_Chodl, Chandelier, Vip. Deep-layer non-IT types: L5_ET, L6_CT, L6b, L5.6_NP. Glia: Astro, Oligo, OPC, Micro. Vascular and immune types: Endo, PC, SMC, VLMC, Immune. The lists are also in `manifest.json` (`comparator_groups`).

## Published subtype DE

We compared the published DE of each author subtype (Emani et al. 2024, `source_data/DEGs/*_DEGcombined.csv`) with our DE on the same subtype, on the merged author type and on our pipeline label. As in Figure S7 (`psychENCODE-author-de-comparison/bin/compare_de_results.py`), we used the published lfc/lfcSE against our limma t over shared genes, for Schizophrenia, Bipolar, ASD and Aging (vs our Age_death). The published lists lack Sst_Chodl, and ASD lacks L6_IT_Car3 and Lamp5_Lhx6, so IT and Lamp5 used three or four contrasts per subtype and Sst used Sst only (13,052 to 15,322 genes per pair; `summary.tsv`, `manifest_plot_t.json`). Our subtype values matched Figure S7 for all 30 pairs.

## Figure (t statistics)

`subtype_vs_pooled_t.png` (from `plot_subtype_vs_pooled_t.py`; sha256 in `manifest_plot_t.json`). `summary.tsv` keeps every comparison, `subtype_vs_siblings` and cells per donor (also `cells_per_donor.tsv`).

Spearman ρ of t statistics, one point per contrast and black bars for the mean over contrasts, categories sorted by mean ρ. (A) Our DE: each author subtype against its pooled author type, other subtypes, the nearest separate types (deep-layer non-IT for IT, other GABAergic for Lamp5 and Sst) and glia. (B) Published subtype DE against our DE on the same subtype, the pooled author type and our label. Our pooled DE agreed with the published subtype results almost as well as our matched subtype DE (IT 0.53 vs 0.57, Lamp5 0.46 vs 0.50, Sst 0.54 vs 0.54), and pooled author and our label gave the same values.
