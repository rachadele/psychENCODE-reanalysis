# Comparison of Reanalysis vs. Published DE Results (Cell Class Level)

## Directory Structure

```
results/comparison/class/
├── all_corr/                          # Full pairwise Spearman correlations per contrast
│   └── {contrast}/
│       ├── figs/                      # Heatmaps of pairwise correlations (log2FC and p-value)
│       └── files/                     # TSV: all pipeline-vs-author cell type pair correlations
├── average_corr/                      # Averaged Spearman correlations
│   ├── overall_average/               # Grand mean across all contrasts and cell types
│   ├── per_celltype/                  # Mean correlation per pipeline cell type label
│   ├── filtered_per_contrast/         # Matching cell type pairs only, per contrast
│   └── figs/                          # Boxplots and stripplots of Spearman distributions
└── de_overlap/                        # Jaccard index of DE gene set overlap
    ├── contrast_overlaps/             # Full pairwise Jaccard per contrast (all cell type pairs)
    │   └── {contrast}/
    └── average_overlaps/
        ├── overall_average/           # Grand mean Jaccard
        ├── per_celltype/              # Mean Jaccard per pipeline cell type label
        ├── filtered_per_contrast/     # Matching cell type pairs only, per contrast
        └── figs/                      # Boxplots and stripplots of Jaccard distributions
```

Each contrast directory has both a `figs/` and `files/` subdirectory. Figures include heatmaps for pairwise log2FoldChange and p-value correlations. The full pairwise TSVs (379 lines each) contain all cross-cell-type comparisons; filtered files contain only matched pipeline-to-author cell type pairs (27 pairs).

## Contrasts Analyzed

- **Psychiatric disorders**: Schizophrenia, Bipolar Disorder, ASD, MDD, PTSD, Williams Syndrome (all vs. Control)
- **Covariates**: Age at death, Biological Sex (male vs. female), PMI
- **Ancestry**: EUR, AMR, EAS, SAS, UNKNOWN (all vs. AFR)
- **Intercept**

## Overall Summary

| Metric | Mean | SD |
|--------|------|----|
| Spearman correlation (log2FC) | 0.761 | 0.218 |
| Spearman correlation (p-value) | 0.657 | 0.263 |
| Jaccard index (% DE gene overlap) | 47.4% | 32.3% |

The reanalysis pipeline produces results that are overall well-correlated with the published analysis at the log2FC level (mean rho = 0.76), with more variability in p-value agreement and DE gene set overlap.

## Per Cell Type Agreement

### Spearman Correlation of log2FoldChange

| Cell Type (Pipeline Label) | Mean Spearman (log2FC) | SD | N Pairs |
|---|---|---|---|
| Oligodendrocyte precursor cell | 0.982 | 0.013 | 8 |
| VIP GABAergic interneuron | 0.976 | 0.014 | 8 |
| Astrocyte | 0.967 | 0.027 | 8 |
| PVALB GABAergic interneuron | 0.965 | 0.030 | 8 |
| Chandelier PVALB interneuron | 0.955 | 0.048 | 8 |
| Oligodendrocyte | 0.942 | 0.057 | 8 |
| SNCG GABAergic interneuron | 0.865 | 0.063 | 8 |
| PAX6 GABAergic interneuron | 0.854 | 0.091 | 8 |
| LAMP5 GABAergic interneuron | 0.832 | 0.079 | 16 |
| L2-6 IT glutamatergic neuron | 0.745 | 0.153 | 40 |
| SST GABAergic interneuron | 0.666 | 0.309 | 16 |
| Deep layer non-IT neuron | 0.655 | 0.182 | 32 |
| Microglial cell | 0.638 | 0.339 | 16 |
| Vascular | 0.605 | 0.189 | 32 |

### Spearman Correlation of P-values

| Cell Type (Pipeline Label) | Mean Spearman (p-value) | SD | N Pairs |
|---|---|---|---|
| Oligodendrocyte precursor cell | 0.972 | 0.016 | 8 |
| Astrocyte | 0.952 | 0.027 | 8 |
| VIP GABAergic interneuron | 0.950 | 0.028 | 8 |
| PVALB GABAergic interneuron | 0.949 | 0.033 | 8 |
| Oligodendrocyte | 0.936 | 0.036 | 8 |
| Chandelier PVALB interneuron | 0.921 | 0.062 | 8 |
| SNCG GABAergic interneuron | 0.754 | 0.102 | 8 |
| PAX6 GABAergic interneuron | 0.735 | 0.113 | 8 |
| LAMP5 GABAergic interneuron | 0.704 | 0.138 | 16 |
| L2-6 IT glutamatergic neuron | 0.632 | 0.217 | 40 |
| SST GABAergic interneuron | 0.558 | 0.383 | 16 |
| Microglial cell | 0.554 | 0.291 | 16 |
| Deep layer non-IT neuron | 0.511 | 0.177 | 32 |
| Vascular | 0.437 | 0.202 | 32 |

### Jaccard Index of DE Gene Overlap (%)

| Cell Type (Pipeline Label) | Mean Jaccard (%) | SD | N Pairs |
|---|---|---|---|
| Oligodendrocyte | 89.9 | 9.0 | 8 |
| Oligodendrocyte precursor cell | 89.8 | 8.4 | 8 |
| VIP GABAergic interneuron | 82.9 | 9.0 | 7 |
| Astrocyte | 81.2 | 9.1 | 8 |
| PVALB GABAergic interneuron | 71.6 | 29.8 | 8 |
| Chandelier PVALB interneuron | 65.8 | 35.3 | 6 |
| SNCG GABAergic interneuron | 58.4 | 19.9 | 7 |
| LAMP5 GABAergic interneuron | 48.2 | 21.3 | 12 |
| PAX6 GABAergic interneuron | 47.2 | 32.2 | 7 |
| SST GABAergic interneuron | 43.7 | 38.9 | 12 |
| L2-6 IT glutamatergic neuron | 42.1 | 25.7 | 40 |
| Microglial cell | 31.8 | 29.5 | 16 |
| Deep layer non-IT neuron | 28.1 | 22.5 | 29 |
| Vascular | 25.5 | 29.0 | 24 |

## Per-Contrast Findings

### Disorder Contrasts (Matched Cell Type Pairs)

#### Spearman Correlation of log2FoldChange

| Cell Type Pair (Pipeline -> Author) | SCZ | BD | ASD | MDD | PTSD | Williams |
|---|---|---|---|---|---|---|
| OPC -> OPC | 0.969 | 0.991 | 0.978 | 0.996 | 0.985 | 0.983 |
| Astro -> Astro | 0.913 | 0.990 | 0.936 | 0.974 | 0.986 | 0.976 |
| VIP -> Vip | 0.970 | 0.991 | 0.983 | 0.972 | 0.988 | 0.963 |
| PVALB -> Pvalb | 0.969 | 0.988 | 0.987 | 0.963 | 0.986 | 0.932 |
| Oligo -> Oligo | 0.916 | 0.948 | 0.975 | 0.940 | 0.974 | 0.977 |
| Chandelier -> Chandelier | 0.972 | 0.986 | 0.988 | 0.909 | 0.983 | 0.952 |
| Micro -> Micro | 0.881 | 0.931 | 0.887 | 0.965 | 0.958 | 0.868 |
| Micro -> Immune | 0.644 | 0.444 | 0.256 | 0.620 | **-0.193** | 0.618 |
| SST -> Sst | 0.924 | 0.982 | 0.981 | 0.964 | 0.976 | 0.943 |
| SST -> Sst_Chodl | 0.275 | 0.494 | 0.466 | 0.239 | 0.416 | 0.420 |
| Deep non-IT -> L5_ET | 0.525 | 0.716 | 0.384 | 0.631 | 0.636 | 0.014 |
| Vascular -> VLMC | 0.426 | 0.247 | 0.462 | 0.290 | 0.186 | 0.413 |

Notable: The Immune cell subtype under microglial cell shows a **negative** correlation for PTSD (rho = -0.19), indicating substantial disagreement. Sst_Chodl and VLMC consistently show low agreement across disorders.

#### Jaccard Index of DE Gene Overlap (%)

| Cell Type Pair | SCZ | BD | ASD | MDD | PTSD | Williams |
|---|---|---|---|---|---|---|
| OPC -> OPC | 84.4 | 93.8 | 79.7 | 100.0 | 85.5 | 80.3 |
| Astro -> Astro | 74.4 | 91.0 | 75.5 | 69.2 | 73.0 | 85.1 |
| VIP -> Vip | 63.6 | 89.1 | 82.7 | -- | 87.1 | 82.1 |
| PVALB -> Pvalb | 67.2 | 90.7 | 84.9 | 0.0 | 75.6 | 83.5 |
| Oligo -> Oligo | 100.0 | 93.2 | 81.6 | 100.0 | 75.7 | 82.4 |
| Micro -> Micro | 40.4 | 70.3 | 67.8 | 43.8 | 23.9 | 19.0 |
| Micro -> Immune | 0.0 | 0.0 | 7.7 | 0.0 | 0.0 | 31.8 |
| SST -> Sst_Chodl | 0.0 | 0.0 | 0.0 | -- | -- | 9.6 |
| Deep non-IT -> L5_ET | 0.6 | 19.3 | 1.5 | 0.0 | 16.1 | 0.0 |
| Vascular -> VLMC | 0.0 | 0.0 | 0.3 | -- | -- | 0.0 |
| Vascular -> SMC | 0.0 | 0.0 | 7.9 | -- | -- | 0.0 |

"--" indicates missing/empty values (likely no significant DE genes in one or both analyses).

MDD has the sparsest DE gene overlap, with most cell types showing 0% or missing values. Bipolar Disorder shows the most consistent overlap across cell types.

### Covariate Contrasts (Selected Matched Pairs)

#### Age at Death

High correlations across most cell types (OPC: 0.997, Astro: 0.982, Pvalb: 0.990). Jaccard indices also high (OPC: 95.1%, Oligo: 92.3%, Astro: 89.6%). Low agreement in vascular subtypes (VLMC log2FC rho = 0.485, Jaccard = 1.4%) and Immune cells (log2FC rho = 0.200, Jaccard = 1.8%).

#### Biological Sex (Male vs. Female)

Generally lower correlations than other contrasts. Weakest agreement in deep layer non-IT neurons (L5_ET: 0.291, L6_CT: 0.365) and vascular subtypes (SMC: 0.415, VLMC: 0.426). Jaccard indices are relatively high for abundant cell types (OPC: 100%, Chandelier: 95.8%, PAX6: 95.2%) but low for L6_CT (5.5%) and L6_IT_Car3 (9.8%).

## Key Observations

1. **Highly reproducible cell types**: OPCs, astrocytes, VIP interneurons, PVALB interneurons, and oligodendrocytes show consistently high agreement (rho > 0.93 for log2FC, Jaccard > 70%) across contrasts.

2. **Poorly reproducible cell types/subtypes**: Vascular subtypes (especially VLMC and SMC), Sst_Chodl, Immune cells, and L5_ET neurons show systematically low agreement. These may reflect differences in cell type annotation granularity, low cell numbers, or methodological differences in pseudobulk aggregation.

3. **Contrast-dependent reproducibility**: Bipolar Disorder shows the strongest overall agreement; MDD shows the weakest DE overlap (many cell types with 0% or missing Jaccard values despite moderate correlation), suggesting the reanalysis detected similar effect directions but fewer significant genes.

4. **log2FC vs. p-value agreement**: log2FC correlations are consistently higher than p-value correlations, indicating that effect size directions and magnitudes are more reproducible than statistical significance calls. This is expected given sensitivity to sample composition and multiple testing procedures.

5. **PTSD Immune cell disagreement**: The negative correlation (rho = -0.19) for Immune cells in the PTSD contrast is a notable outlier warranting investigation.

## Files Skipped

The full pairwise TSVs under `all_corr/{contrast}/files/` and `de_overlap/contrast_overlaps/{contrast}/` (379 lines each, covering all cross-cell-type comparisons) were sampled but not fully summarized here. They contain correlations and Jaccard indices for every pipeline cell type vs. every author cell type, including non-matching pairs.
