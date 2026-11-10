# Author subtypes vs pooled author types in DE

We asked whether the author subtypes that our pipeline pools (five IT types into L2/3-6 IT, Lamp5 with Lamp5_Lhx6, Sst with Sst_Chodl) gave the same DE as their pooled group. We reused existing limma-voom results: 27 unmerged author types and 21 merged author types (`results_v2.1.2_limma/merged_author_2026-11-09`). For the eight Figure 4 contrasts, `subtype_vs_pooled.py` computed Spearman ρ of log2FC and of t over shared genes between each subtype and its merged group, between subtypes of a group, between each subtype and named groups of author types outside its pooled group (below), and between the merged group and our label. We averaged pairs within each contrast and report mean ± SD over contrasts (`summary.tsv`; all pairs in `pairwise_rho.tsv`). Cells per donor was the median over DE donors (`cells_per_donor.tsv`). The merged-vs-ours values matched `de_agreement/celltype_merged.tsv`.

Subtypes agreed with their merged group at ρ = 0.76 ± 0.08 (IT), 0.86 ± 0.06 (Lamp5) and 0.68 ± 0.12 (Sst), and with each other at 0.66, 0.64 and 0.32. Agreement fell with subtype size: Sst_Chodl (4 cells per donor) agreed with Sst at 0.32. `manifest.json` holds sha256 and timestamps.

## Comparator groups

IT subtypes were compared with GABAergic types, deep-layer non-IT types, glia, and vascular and immune types; Lamp5 and Sst subtypes with IT types, GABAergic types outside their pooled group, deep-layer non-IT types, glia, and vascular and immune types. IT types: L2.3_IT, L4_IT, L5_IT, L6_IT, L6_IT_Car3. GABAergic types: Lamp5, Lamp5_Lhx6, Pax6, Pvalb, Sncg, Sst, Sst_Chodl, Chandelier, Vip. Deep-layer non-IT types: L5_ET, L6_CT, L6b, L5.6_NP. Glia: Astro, Oligo, OPC, Micro. Vascular and immune types: Endo, PC, SMC, VLMC, Immune. The lists are also in `manifest.json` (`comparator_groups`).

## Figure (t statistics)

`subtype_vs_pooled_t.png` (from `plot_subtype_vs_pooled_t.py`; sha256 in `manifest_plot_t.json`). The figure leaves out pooled author vs our label (shown in Figure 4) and the other comparator groups; `summary.tsv` keeps every comparison, `subtype_vs_siblings` and cells per donor (also `cells_per_donor.tsv`).

Spearman ρ of limma moderated t between each author subtype and its pooled author type, between subtypes, between each subtype and the nearest separate types (deep-layer non-IT for IT, other GABAergic for Lamp5 and Sst) and between each subtype and glia, sorted by mean ρ, one point per contrast and black bars for the mean of eight contrasts. Subtypes agreed best with their pooled type (ρ = 0.82, 0.86, 0.68 for IT, Lamp5, Sst); IT and Lamp5 subtypes agreed with each other (0.70, 0.61) more than with the nearest separate types (0.58, 0.52), Sst and Sst Chodl did not (0.33 vs 0.45), and all exceeded glia (0.28, 0.31, 0.24).
