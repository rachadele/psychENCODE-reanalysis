# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **Nextflow DSL2 pipeline** for differential expression (DE) analysis of single-cell RNA-seq data from 8 PsychENCODE cohorts (CMC, PTSDBrainomics, SZBDMulti-Seq, MultiomeBrain, UCLA-ASD, Velmeshev-2019.1, Velmeshev-2019.2, DevBrain). It compares DE results between two cell type annotation schemes: standardized GEMMA annotations and author-submitted annotations.

## Running the Pipeline

```bash
# Class-level analysis
nextflow run main.nf -params-file params.class.json -profile conda

# Subclass-level analysis
nextflow run main.nf -params-file params.subclass.json -profile conda

# Skip DE analysis, run only comparison (requires existing results)
nextflow run main.nf -params-file params.class.json -profile conda --skip_de_analysis \
    --pavlab_deseq_results "./results/author_submitted_false/class/DESeq2/gemma/**tsv" \
    --author_label_deseq_results "./results/author_submitted_true/class/DESeq2/gemma/**tsv"

# Skip comparison, run only DE analysis
nextflow run main.nf -params-file params.class.json -profile conda --skip_comparison
```

Individual Python scripts in `bin/` can be run directly; they use `argparse` for CLI arguments.

## Architecture

### Two-Stage Pipeline

**Stage 1: DE Analysis** (two parallel branches)
- GEMMA branch (`from_gemma=true`): Fetches pseudobulks via `gemma-cli-staging` API → aggregates by cell type → runs DESeq2
- Manual branch (`from_gemma=false`): Loads H5AD files → aggregates pseudobulks → runs DESeq2
- Runs twice: once with `author_submitted=false` (GEMMA/pavlab annotations) and once with `author_submitted=true`

**Stage 2: Comparison** (`gemma_comparison` subworkflow)
- Computes pairwise log2FC correlations between annotation schemes per study/cell type/contrast
- Computes Jaccard overlap of significant DE genes
- Averages metrics across studies

### File Layout

```
main.nf                         # Entry point; delegates to workflows/
workflows/psychencode_reanalysis.nf  # Main orchestration; calls both subworkflows
subworkflows/local/
    gemma_de_analysis.nf        # Stage 1 GEMMA pathway
    manual_de_analysis.nf       # Stage 1 manual H5AD pathway
    gemma_comparison.nf         # Stage 2 comparison
modules/local/                  # 10 Nextflow modules (one process each)
conf/
    base.config                 # Executor settings (local, queue 90), resource labels
    modules.config              # Conda env assignments per process
bin/                            # R and Python scripts called by modules
    DESeq2_analysis.R           # Core DE analysis (negative binomial GLM)
    aggregate_celltypes_*.py    # Pseudobulk aggregation
    aggregate_pairwise.py       # Pairwise log2FC correlation
    de_overlap.py               # Jaccard overlap calculation
    average_corr.py             # Average correlation metrics
    average_de_overlaps.py      # Average overlap metrics
    utils.py                    # Shared visualization/data utilities
params.class.json               # Cell type mappings + parameters for class-level
params.subclass.json            # Cell type mappings + parameters for subclass-level
```

### Key Parameters

| Parameter | Values | Description |
|-----------|--------|-------------|
| `level` | `"class"` / `"subclass"` | Cell type annotation granularity |
| `from_gemma` | bool | GEMMA API vs. local H5AD files |
| `author_submitted` | bool | Toggle annotation scheme |
| `skip_de_analysis` | flag | Skip Stage 1, use existing results |
| `skip_comparison` | flag | Skip Stage 2 |
| `filter_genes` | bool | Pre-filter genes before DE |
| `de_method` | `"deseq2"` / `"limma_voom"` | DE method; limma_voom follows Gemma's LinearModelAnalyzer (GEMMA branch only) |
| `pavlab_ct_pseudobulks` / `author_ct_pseudobulks` | dir | Start an arm from local cell-type pseudobulks instead of the Gemma download |

### DESeq2 Design Formula

GEMMA mode: `~ Disorder + Age_death + PMI + Biological_Sex + X1000G_ancestry`  
Manual mode: adds `+ avg_UMI_sample`

Contrasts: 6 disorders (ASD, Bipolar Disorder, MDD, SCZ, Rett Syndrome, PTSD), 4 ancestry groups, plus covariates.

### Output Structure

```
results/
├── author_submitted_{true|false}/{class|subclass}/
│   ├── experiment_pseudobulks/
│   ├── aggregated_pseudobulks/
│   └── DESeq2/gemma/           # Per-cell-type TSVs with LFC, p-values
└── comparison/{class|subclass}/
    ├── all_corr/               # Pairwise correlations per study
    ├── de_overlap/             # Jaccard overlaps
    └── average_corr/           # Averaged metrics + heatmaps
```

### Dependencies

- **Nextflow** >= 23.04.0 with DSL2
- **Python env** (`psychencode-env`): scanpy, pandas, numpy, scipy, seaborn, matplotlib, anndata
- **R env** (`r-deseq2`): R 4.3, DESeq2, edgeR, Seurat, ggplot2, dplyr, tidyr, argparse, reticulate
- **External**: `gemma-cli-staging` (GEMMA REST API client)

Conda environments are defined in `conf/modules.config` and activated per-process.
