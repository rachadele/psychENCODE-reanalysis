# Methods: psychENCODE Reanalysis Pipeline

## Overview

This Nextflow DSL2 pipeline (v2.0.0) performs differential expression (DE) analysis of pseudobulk single-cell RNA-seq data and compares DE results obtained using two independent cell type annotation schemes: (1) standardized annotations from the GEMMA single-cell pipeline (sc-pipeline-1.2.0) and (2) original author-submitted annotations. The pipeline operates in two stages: DE analysis and cross-annotation comparison. Analyses were performed at the cell class level.

## Cohorts

Eight cohorts from the PsychENCODE consortium were analyzed:

1. CMC (CommonMind Consortium)
2. PTSDBrainomics
3. SZBDMulti-Seq
4. MultiomeBrain
5. UCLA-ASD
6. Velmeshev-2019.1
7. Velmeshev-2019.2
8. DevBrain

## Stage 1: Differential Expression Analysis

### Data Acquisition

#### GEMMA Pathway

Study-level pseudobulk expression matrices were retrieved from the GEMMA database using `gemma-cli-staging getSingleCellDataMatrix`, aggregated by assay and cell type assignment. Separate retrievals were performed for standardized pipeline annotations (`sc-pipeline-1.2.0-class`) and author-submitted annotations.

#### Manual H5AD Pathway

Alternatively, AnnData H5AD objects were loaded directly with custom cell type annotation files. Pseudobulk aggregation was performed using `scanpy.get.aggregate()` with sum, count_nonzero, and mean aggregation functions, grouping by sample and cell type. QC metrics (total UMI counts, number of genes detected) were computed with `scanpy.pp.calculate_qc_metrics()`. An optional filter removed samples with fewer than 50 cells per cell type.

### Pseudobulk Aggregation Across Studies

Per-study pseudobulk matrices were combined by cell type across all cohorts using outer joins, filling missing genes with zero counts. Cell types with fewer than 16 total samples were excluded from downstream analysis.

For the manual pathway, per-sample metadata was extracted including: sample ID, 1000G ancestry, age at death, biological sex, cohort, disorder diagnosis, individual ID, post-mortem interval (PMI), RIN, and pH. Average UMI per sample (`avg_UMI_sample`) was computed as an additional QC covariate.

### Differential Expression with DESeq2

Pseudobulk count matrices were analyzed per cell type using DESeq2 (R 4.3 environment). The analysis proceeded as follows:

1. **Preprocessing**: Samples with zero library size were removed. Gene and sample identifiers were matched between the count matrix and metadata.

2. **Gene filtering** (optional): Genes were retained if CPM > 0.5 in at least 30% of samples and total counts >= 20 across all samples. CPM was computed using `edgeR::cpm()`.

3. **DESeq2 model**: A DESeqDataSet was constructed with the design formula:
   - GEMMA mode: `~ Disorder + Age_death + PMI + Biological_Sex + X1000G_ancestry`
   - Manual mode: `~ Disorder + Age_death + PMI + Biological_Sex + X1000G_ancestry + avg_UMI_sample`

   "Control" was set as the reference level for the Disorder variable.

4. **Size factor estimation**: `estimateSizeFactors()` with `type = "poscounts"` to handle sparse count data.

5. **Differential expression testing**: `DESeq()` with default parameters (negative binomial generalized linear model, Wald test).

6. **Results extraction**: For each coefficient in the model (disorder contrasts, covariates, ancestry), results were extracted using `results()`, sorted by p-value, and saved as TSV files. Output columns: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, cell_type.

7. **Diagnostic plots**: PCA (colored by cohort), dispersion estimates, Cook's distance histogram, MA plots, and p-value/log2FC distribution histograms were generated per cell type.

### Contrasts Tested

- **Disorder effects** (vs. Control): Schizophrenia, Bipolar Disorder, ASD, MDD, PTSD, Williams Syndrome
- **Covariates**: Age at death (continuous), Biological Sex (male vs. female), PMI (continuous)
- **Ancestry** (vs. AFR): EUR, AMR, EAS, SAS, UNKNOWN
- **Intercept**

## Stage 2: Comparison of DE Results

### Cell Type Mapping

Author-submitted cell type labels were mapped to standardized pipeline labels using curated lookup tables (`params.class.json`). For example, author label "Astro" maps to pipeline label "astrocyte"; "L2.3_IT", "L4_IT", "L5_IT", "L6_IT", and "L6_IT_Car3" all map to "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron" at the class level. A total of 26 author labels were mapped to 14 pipeline class labels.

### Pairwise Spearman Correlation

For each contrast, DESeq2 results from the two annotation schemes were compared across all pairwise combinations of cell types. For each pair, results were merged on gene name, and Spearman rank correlations were computed for:

- **log2FoldChange**: Measures agreement in effect size direction and magnitude
- **p-value**: Measures agreement in statistical significance

Correlations were visualized as heatmaps (pipeline cell types vs. author cell types) per contrast.

### DE Gene Set Overlap (Jaccard Index)

For each contrast and cell type pair, significant DE genes were identified using:
- Adjusted p-value (padj) < 0.05
- |log2FoldChange| > 0 (no fold-change threshold by default)

The Jaccard index was computed as:

```
J(A, B) = |A ∩ B| / |A ∪ B| × 100%
```

where A and B are the sets of significant DE genes from the two annotation schemes.

### Averaging and Summarization

Spearman correlations and Jaccard indices were filtered to retain only matched cell type pairs (author label mapping to the corresponding pipeline label). Cofactor contrasts (Intercept, PMI, ancestry terms) were excluded from summary statistics.

Summary statistics were computed at three levels:

1. **Overall average**: Mean and standard deviation across all matched cell type pairs and disorder/covariate contrasts.
2. **Per cell type**: Mean, standard deviation, and number of pairs per pipeline cell type label, aggregated across contrasts.
3. **Per contrast**: All matched cell type pairs for each individual contrast.

### Visualization

- **Boxplots**: Distribution of Spearman correlations and Jaccard indices across cell types, with mean markers.
- **Stripplots**: Individual data points per cell type, colored by contrast, with jitter.
- **Heatmaps**: Pairwise correlation matrices per contrast (log2FC: red colormap; p-value: blue colormap).

## Software and Environment

- **Workflow manager**: Nextflow (DSL2, >= 23.04.0)
- **Python environment**: scanpy, pandas, numpy, scipy, seaborn, matplotlib
- **R environment** (v4.3): DESeq2, edgeR, Seurat, ggplot2, tidyr, dplyr, argparse, reticulate
- **Data retrieval**: gemma-cli-staging (GEMMA REST API)
- **Execution**: Local executor, conda-managed environments, queue size 90

## Reproducibility

- Random seed: 42 (set in DESeq2 analysis)
- Nextflow caching enabled for intermediate results
- All parameters configurable via `nextflow.config`, `conf/base.config`, `conf/modules.config`, and JSON parameter files