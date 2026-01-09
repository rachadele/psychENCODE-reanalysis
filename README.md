# psychENCODE-reanalysis

This repository orchestrates large-scale differential expression (DE) analysis of single-cell and bulk RNA-seq data, enabling comparison and aggregation of DE results across studies, cell types, and analysis methods (e.g., GEMMA vs. author-submitted).

## Project Structure
- **Nextflow Pipelines**: Top-level workflows (`author-de-vs-gemma-pseudobulk-de.nf`, `gemma-vs-gemma-pseudobulk-de.nf`, `get-gemma-de-results.nf`) coordinate data wrangling, aggregation, and DE analysis. Pipelines import modular processes from `modules/`.
- **Modules**: Each `.nf` in `modules/` is a reusable Nextflow process (e.g., aggregation, DESeq2, wrangling). These are included in top-level workflows.
- **R/Python Scripts**: Located in `bin/`, these scripts perform core data processing (e.g., `DESeq2_analysis.R`, `aggregate_celltypes_gemma.py`). Many are invoked by Nextflow processes.
- **Configuration**: `nextflow.config` and JSON parameter files (e.g., `params.class.json`) control pipeline behavior, input locations, and cell type mappings.
- **Results**: Output is organized by analysis type and parameters in results folders (e.g., `results_author_submitted_false_from_gemma_true/`).

## Usage

### Running Pipelines
Use `runall.sh` for typical batch runs, or invoke Nextflow directly. Example:
```bash
nextflow gemma-vs-gemma-pseudobulk-de.nf -params-file params.class.json -resume --outdir ./gemma-vs-gemma-class
```

### Parameterization
- Adjust `params.class.json` or `params.subclass.json` to change cell type mappings or analysis level.

### Module Development
- Add new processes as `.nf` files in `modules/` and include them in top-level workflows.

### R/Python Integration
- Scripts in `bin/` are called by Nextflow; ensure CLI arguments match those in the Nextflow process definitions.

## Conventions
- **Cell Type Naming**: Consistent mapping between GEMMA and project-specific cell type names is maintained in JSON param files.
- **Reproducibility**: Random seeds are set in R scripts (e.g., `set.seed(42)`).
- **Resource Management**: Memory and CPU settings are controlled in `nextflow.config` and process directives.
- **Results Structure**: Output directories encode parameter choices (e.g., `results_author_submitted_false_from_gemma_true`).

## Integration
- **External Data**: Input data (e.g., H5AD, metadata) is referenced in config/params files and must be present at specified paths.
- **Conda Environments**: Nextflow and R scripts expect conda to be enabled (`conda.enabled = true` in config). Some R scripts specify conda envs via `reticulate`.

## Examples
- To add a new aggregation step, create a module in `modules/`, then include it in a top-level workflow.
- To run DE analysis for a new cell type, update the relevant param file and rerun the pipeline.

## References
- Top-level workflows: `author-de-vs-gemma-pseudobulk-de.nf`, `gemma-vs-gemma-pseudobulk-de.nf`, `get-gemma-de-results.nf`
- Modules: `modules/`
- Scripts: `bin/`
- Config: `nextflow.config`, `params.class.json`, `params.subclass.json`
- Example run: `runall.sh`

---
For questions about unclear conventions or missing documentation, please check recent pipeline usage in `runall.sh` or contact the maintainers.
