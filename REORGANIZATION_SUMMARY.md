# nf-core Reorganization Summary

## Overview

The pipeline has been reorganized to follow nf-core best practices with a single entry point and modular stage execution.

## New Directory Structure

```
psychENCODE-reanalysis/
├── main.nf                          # Single entry point
├── nextflow.config                  # Main config with includeConfig
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
└── bin/                             # Scripts (unchanged)
```

## New Stage Control Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_de_analysis` | `true` | Run Stage 1 (DE analysis) |
| `--run_comparison` | `true` | Run Stage 2 (comparison) |
| `--skip_de_analysis` | `false` | Skip Stage 1 |
| `--skip_comparison` | `false` | Skip Stage 2 |

## Usage Examples

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
nextflow run main.nf -params-file params.class.json -profile conda --skip_de_analysis
```

## Key Changes

### 1. Single Entry Point (`main.nf`)
- Replaces the old multi-workflow approach (`runall.sh`)
- Prints pipeline info and handles completion/error events
- Calls the main workflow from `workflows/psychencode_reanalysis.nf`

### 2. Modular Subworkflows
- **GEMMA_DE_ANALYSIS**: Fetches GEMMA data → aggregates → DESeq2
- **MANUAL_DE_ANALYSIS**: H5AD files → aggregates → DESeq2
- **GEMMA_COMPARISON**: Pairwise correlations and Jaccard overlaps

### 3. Process Naming Convention
- All processes renamed to SCREAMING_SNAKE_CASE (nf-core standard)
- Added `label` directives for resource management
- Example: `GetGemmaPseudobulks` → `GET_GEMMA_PSEUDOBULKS`

### 4. Configuration Organization
- `conf/base.config`: Process resources (cpus, memory, time)
- `conf/modules.config`: Conda environment assignments
- `params.class.json` / `params.subclass.json`: Cell type mappings (use with `-params-file`)

### 5. Parameter Schema
- `nextflow_schema.json` enables validation and help text
- Parameters organized into groups: input/output, stage control, analysis options

## Files Preserved (Not Modified)

- `bin/` - All Python and R scripts unchanged
- `params.class.json` and `params.subclass.json` - Still available for reference
- Old workflow files - Kept for reference:
  - `get-gemma-de-results.nf`
  - `gemma-vs-gemma-pseudobulk-de.nf`
  - `author-de-vs-gemma-pseudobulk-de.nf`
- Old modules in `modules/*.nf` - Kept for reference

## Migration Notes

1. The old workflows are still present and can be run directly if needed
2. The new `main.nf` provides the unified entry point
3. Old `runall.sh` is now obsolete - use `main.nf` with stage control params instead

## Verification Commands

```bash
# Check pipeline syntax
nextflow run main.nf -preview

# View help
nextflow run main.nf --help
```
