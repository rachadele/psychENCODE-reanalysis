# psychENCODE-reanalysis

Differential expression (DE) analysis pipeline for single-cell RNA-seq data, comparing GEMMA-derived and author-submitted cell type annotations. Reorganized to follow nf-core best practices.

## Project Structure

```
psychENCODE-reanalysis/
├── main.nf                          # Single entry point
├── nextflow.config                  # Main configuration
├── nextflow_schema.json             # Parameter validation schema
│
├── workflows/
│   └── psychencode_reanalysis.nf    # Main workflow orchestration
│
├── subworkflows/local/
│   ├── gemma_de_analysis.nf         # Stage 1: GEMMA DE analysis
│   ├── manual_de_analysis.nf        # Stage 1 alt: Manual H5AD analysis
│   └── gemma_comparison.nf          # Stage 2: Compare DE results
│
├── modules/local/
│   ├── get_gemma_pseudobulks/main.nf
│   ├── aggregate_celltypes_gemma/main.nf
│   ├── aggregate_data_manual/main.nf
│   ├── aggregate_celltypes_manual/main.nf
│   ├── deseq2_analysis_gemma/main.nf
│   ├── deseq2_analysis_manual/main.nf
│   ├── aggregate_pairwise/main.nf
│   ├── de_overlap/main.nf
│   ├── average_de_overlaps/main.nf
│   └── average_de_correlations/main.nf
│
├── conf/
│   ├── base.config                  # Process resources & executor
│   └── modules.config               # Conda environments per process
│
├── bin/                             # R/Python scripts
│
└── deprecated/                      # Old workflow files (for reference)
```

## Usage

### Run Full Pipeline (Class level)
```bash
nextflow run main.nf -params-file params.class.json -profile conda --outdir ./results-class
```

### Run Full Pipeline (Subclass level)
```bash
nextflow run main.nf -params-file params.subclass.json -profile conda --outdir ./results-subclass
```

### Run Only DE Analysis (Stage 1)
```bash
nextflow run main.nf -params-file params.class.json -profile conda --skip_comparison
```

### Run Only Comparison (Stage 2) with Pre-existing Results
```bash
nextflow run main.nf -params-file params.class.json -profile conda --skip_de_analysis \
    --pavlab_deseq_results "./results/author_submitted_false/class/DESeq2/gemma/**tsv" \
    --author_label_deseq_results "./results/author_submitted_true/class/DESeq2/gemma/**tsv"
```

## Stage Control Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_de_analysis` | `true` | Run Stage 1 (DE analysis) |
| `--run_comparison` | `true` | Run Stage 2 (comparison) |
| `--skip_de_analysis` | `false` | Skip Stage 1 |
| `--skip_comparison` | `false` | Skip Stage 2 |

## Pipeline Stages

### Stage 1: Differential Expression Analysis
- **GEMMA_DE_ANALYSIS**: Fetches GEMMA data → aggregates by cell type → DESeq2
- **MANUAL_DE_ANALYSIS**: H5AD files → aggregates → DESeq2

### Stage 2: Comparison
- **GEMMA_COMPARISON**: Pairwise correlations and Jaccard overlaps between annotation schemes

## Configuration

- `conf/base.config`: Process resources (cpus, memory, time)
- `conf/modules.config`: Conda environment assignments
- `params.class.json` / `params.subclass.json`: Cell type mappings (use with `-params-file`)

## Output Structure

Results are organized by stage and parameters:
```
results/
├── author_submitted_{true|false}/
│   └── {class|subclass}/
│       ├── experiment_pseudobulks/    # Stage 1: Raw pseudobulks
│       ├── aggregated_pseudobulks/    # Stage 1: Aggregated by cell type
│       └── DESeq2/                    # Stage 1: DE results
│
└── comparison/
    └── {class|subclass}/
        ├── all_corr/                  # Stage 2: Pairwise correlations
        ├── de_overlap/                # Stage 2: Jaccard overlaps
        └── average_corr/              # Stage 2: Averaged correlations
```

## Conventions

- **Cell Type Naming**: Mapping between GEMMA and project-specific names in JSON param files
- **Reproducibility**: Random seeds set in R scripts (`set.seed(42)`)
- **Process Names**: SCREAMING_SNAKE_CASE (nf-core standard)

## Requirements

- Nextflow >= 23.04.0
- Conda environments:
  - `scanpyenv` (Python processes)
  - `r4.3` (DESeq2 processes)

## Verification

```bash
# Check pipeline syntax
nextflow run main.nf -preview

# View help
nextflow run main.nf --help
```

---
For questions, check `deprecated/` for old workflow patterns or contact the maintainers.

## Author-to-pipeline label map

`assets/author_to_pipeline_map.class.json` and `.subclass.json` map the 27 author cell types (8 cohorts and `source_data/DEGs`) to sc-annotation-pipeline v2.1.2 labels. `--label_map` selects the file (set in `params.class.json` / `params.subclass.json`); with `--label_map null` the inline `gemma_to_gemma_map` is used, which holds the v2.0.0 labels needed to reproduce `results/`.

Merges at both levels: L2.3_IT, L4_IT, L5_IT, L6_IT, L6_IT_Car3 to the L2/3-6 IT neuron; Lamp5, Lamp5_Lhx6 to lamp5; Sst, Sst_Chodl to sst. Class level only: L5.6_NP, L5_ET, L6_CT, L6b to deep layer non-IT; Endo, PC, SMC, VLMC to brain vascular cell. Immune maps to `immune` (class) and `macrophage` (subclass, the most frequent non-microglial immune subclass in all eight cohorts); Micro maps to microglial cell. Pipeline subclasses B cell, T cell, myeloid cell and erythroid lineage cell have no author counterpart.

## DE method: DESeq2 or limma-voom

`--de_method deseq2` (default) runs DESeq2. `--de_method limma_voom` runs `bin/limma_voom_analysis.R`, which follows Gemma's LinearModelAnalyzer for RNA-seq: log2 CPM with Gemma's offsets, Gemma's gene filters (at least 3 expressing cells, nonzero and at least 0.01 log2 CPM variance), voom weights, eBayes moderated t. The design and sample filters match the DESeq2 path. Output keeps the DESeq2 layout and adds DESeq2-named columns, so the comparison steps run unchanged.

`--pavlab_ct_pseudobulks` / `--author_ct_pseudobulks` start an arm from local cell-type pseudobulks (`bin/build_local_pseudobulks.py` with `--mask --nnz`, then `bin/aggregate_celltypes_gemma.py`) instead of the Gemma download. With those, limma-voom uses the adapter's library sizes (masked cells excluded, as in Gemma) and expressing-cell counts. `--author_labels_from_cell_id` builds author-label pseudobulks from the h5ad cell ids. `results_v2.1.2_limma/scripts/` holds the run used for `results_v2.1.2_limma/`.
