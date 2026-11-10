# DE agreement analysis: decision log

Figure 4 of sc-manuscript. Producer: `bin/de_agreement_vs_annotation.py`; outputs in `results_v2.1.2_limma/de_agreement_vs_annotation/`.

## 2026-11-10: Figure 4 uses moderated t, decided by Rachel

For consistency with Figure S7, Figure 4 now correlates limma moderated t statistics (column `stat`) instead of log2 fold changes. All other rules were unchanged. `bin/fig_de_agreement.py --stat t` and `bin/merged_author_fig4_inputs.py --stat t` wrote `fig_de_agreement_subclass_merged_t.png`, `fig4_inputs/fig4b_rho_ci_merged_t.tsv`, `fig4_inputs/fig4_agreement_summary_t.tsv`, `fig4_inputs/null_t/` and `fig4_inputs/manifest_fig4_t.json` (sha256, data version, timestamps). The default `--stat log2FoldChange` reproduces the previous outputs.

| | log2FC | t |
|---|---|---|
| Mean ± SD ρ, 21 merged types | 0.88 ± 0.19 | 0.88 ± 0.17 |
| Mean ± SD ρ, class level (27 types) | 0.78 ± 0.25 | 0.79 ± 0.22 |
| Null of matched mean | 0.36 | 0.38 |
| Adjusted means, matched / mismatched | 0.88 / 0.34 | 0.88 / 0.35 |
| ρ, agreement vs F1 (21) | 0.91 | 0.94 |
| ρ, agreement vs cells per donor (21) | 0.84 | 0.92 |
| Partial ρ, F1 given cells (21) | 0.67 (p = 0.0011) | 0.76 (p = 1.1e-4) |
| Partial ρ, cells given F1 (21) | 0.33 (p = 0.15) | 0.66 (p = 0.0015) |
| Partial ρ, F1 given cells (27) | −0.17 (p = 0.42) | −0.12 (p = 0.57) |
| Partial ρ, cells given F1 (27) | 0.53 (p = 0.0050) | 0.60 (p = 0.0013) |

With t, cells per donor stayed significantly associated with agreement once F1 was controlled for.

Panels D and E (`--pooling_dir subtype_vs_pooled`, approved by Rachel 2026-11-10) add the author subtype vs pooled test from `subtype_vs_pooled/summary.tsv` and `pairwise_rho.tsv`; their inputs are in `manifest_fig4_t.json`.

## 2026-11-09: merged-author arm, decided by Rachel and adopted for Figure 4

We reran the author-label DE with author types mapped to pipeline subclass labels before pseudobulking (`assets/author_to_pipeline_map.subclass.json`), so the five IT types, LAMP5 with LAMP5 LHX6, and SST with SST Chodl each formed one group. All 21 groups matched one pipeline label and none were excluded. QC masks, limma-voom model, donors and contrasts were unchanged. The 18 unpooled types reproduced the current DE exactly. F1 and cells per donor followed the de_donors rule above. Outputs, README and manifest: `results_v2.1.2_limma/merged_author_2026-11-09/`. Rachel approved the merged arm on 2026-11-09. It replaced the exclusion of pooled types: Figure 4 (`bin/fig_de_agreement.py --merged_dir`, inputs from `bin/merged_author_fig4_inputs.py`) now shows the 21 merged types, with the 27-type rows kept in 4B for comparison, and the manuscript text uses the merged values. Choosing one representative author type per pooled label after seeing the results was rejected as data dredging.

| | Current, 27 types | Current, strict 18 | Merged, 21 types |
|---|---|---|---|
| Mean ± SD ρ, log2FC | 0.83 ± 0.21 | 0.86 ± 0.20 | 0.88 ± 0.19 |
| Mean ± SD ρ, t | 0.84 ± 0.19 | 0.86 ± 0.18 | 0.88 ± 0.17 |
| ρ, agreement vs F1 | 0.33 (p = 0.097) | 0.90 (p = 2.9e-7) | 0.91 (p = 1.7e-8) |
| ρ, agreement vs cells per donor | 0.58 (p = 0.0014) | 0.81 (p = 4.5e-5) | 0.84 (p = 1.8e-6) |
| Partial ρ, F1 given cells | −0.17 (p = 0.42) | 0.70 (p = 0.0016) | 0.67 (p = 0.0011) |
| Partial ρ, cells given F1 | 0.53 (p = 0.0050) | 0.25 (p = 0.34) | 0.33 (p = 0.15) |

Source: `merged_author_2026-11-09/de_agreement/agreement_summary.tsv` and `correlations_merged_vs_current.tsv` (agreement = mean ρ of log2FC over eight contrasts).

## 2026-11-09, decided by Rachel. Status: final (2026-11-10). Emani et al. released no DE code; our model pools all donors

1. F1 and cells per donor came from every donor in the DE: all eight cohorts and the same donors for every contrast, because limma fitted one model per cell type across all cohorts without a cohort term. A donor was in a type's model when it had metadata and at least one unmasked cell of the type. F1 was the plain mean of per-sample F1 over those donors (cells never pooled, no mean of cohort means). Cells per donor was the median over the same donors, zeros excluded. Option `--support_rule de_donors` (default).
2. One-to-one types were the strict 18: every author type whose pipeline label includes no other author type. LAMP5 and SST were dropped because they share a label with LAMP5 LHX6 and SST Chodl. Option `--one_to_one strict` (default; `own_label` gives the former 20).
3. One table, `celltype_f1_cells.tsv`, fed Figure 4A (cohort F1 bar), 4B, the supplementary scatter and the cells-per-donor adjustment of 4C and S7B (`benchmark-post-hoc/de_permutation_null/matched_effect_adjusted.py`), including the published arm, which formerly used CMC author counts.
4. The DE model was unchanged: no 50-cell filter and no UMI covariate.

## Alternatives

Partial Spearman ρ of mean agreement with F1 given log10 cells per donor (`rule_alternatives.tsv`):

| Rule | 20 own-label types | 18 strict types |
|---|---|---|
| Case cohorts (previous Figure 4B) | 0.76 | 0.78 |
| CMC only, cells pooled | 0.33 | 0.36 |
| All DE donors (adopted) | 0.61 | 0.70 (p = 0.0016) |

Under the adopted rule, cells per donor given F1 had partial ρ = 0.25 (p = 0.34) over the 18 types and 0.53 over all 27. The case-cohort rule did not match the DE donors (audit, 2026-11-09). Its outputs were kept with the suffix `_case_cohorts` and are reproduced exactly by `--support_rule case_cohorts --one_to_one own_label`.

## Outputs (sha256)

| File | sha256 |
|---|---|
| celltype_f1_cells.tsv | b1a2a14c0d1fafdf96fbd090da8832e80d830e55e6b82fdbea4bf272a76e7b0b |
| rule_alternatives.tsv | b121011d17851f9e50eaefd745cf510859c21073357f0c8ca170fd330f60f80c |
| fig4b_rho_ci.tsv | 8f6c93e23239f7b495337efde275b4ef3b6a769230914db2976d3440c3c4214d |
| correlations_per_cohort.tsv | e983dd8556cd71446304790eaf47d61d48603834cf799184f7dc8f7e20cff5f4 |
| pairs_per_cohort.tsv | 6674e2974cbcd63ccc65ecc835715406a0d1ab39553198c256fac10dae6bf2c4 |
| celltype_summary.tsv | e19017841fd90466a727e94bdfc4964893500f8184c350e495250327718d6a8e |
| fig_de_agreement_subclass.png | bc3dded73ba2a7f8c70a296a910abfe6e634ede194fe077d54c29f3b1205684b |
| fig_de_agreement_scatter_subclass.png | de6265d9da5ee586665ac433d41652bf425d7b3ab498082dc92184424a9a9f1c |
| fig4b_rho_ci_case_cohorts.tsv | ed4d27f12d9ed5a1242ba94cc9ac1c09da7fbdfdc89e40f5d790e389b171d49a |
| de_permutation_null/matched_effect_adjusted.tsv | 6f38967415c5e4f36ca8abca0e6400ab2add094e0e7cb71e12294c7fdf1af396 |

`manifest.json` in the output folder records every input with its sha256 and timestamps; data version: pipeline labels 2.1.2, DE `results_v2.1.2_limma`.
