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

Results are organized by parameters:
```
results/
└── author_submitted_{true|false}/
    └── {class|subclass}/
        ├── experiment_pseudobulks/
        ├── aggregated_pseudobulks/
        ├── DESeq2/
        └── average_corr/
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
