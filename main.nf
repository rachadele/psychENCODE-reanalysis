#!/usr/bin/env nextflow
/*
 * psychENCODE Reanalysis Pipeline
 *
 * Differential expression analysis of single-cell RNA-seq data
 * comparing GEMMA-derived and author-submitted cell type annotations.
 *
 * Usage:
 *   nextflow run main.nf -profile conda --level class
 *   nextflow run main.nf -profile conda --skip_comparison
 *   nextflow run main.nf -profile conda --skip_de_analysis --pavlab_deseq_results "./results/author_submitted_false/class/DESeq2/gemma/**tsv" --author_label_deseq_results "./results/author_submitted_true/class/DESeq2/gemma/**tsv"
 */

nextflow.enable.dsl = 2

// Pipeline version
version = '2.0.0'

// Print pipeline info
log.info """
=========================================
 psychENCODE Reanalysis Pipeline v${version}
=========================================
 Run DE Analysis    : ${params.run_de_analysis && !params.skip_de_analysis}
 Run Comparison     : ${params.run_comparison && !params.skip_comparison}
 Data Source        : ${params.from_gemma ? 'GEMMA' : 'Manual H5AD'}
 Annotation Level   : ${params.level}
 Output Directory   : ${params.outdir}
=========================================
"""

// Include the main workflow
include { PSYCHENCODE_REANALYSIS } from "${projectDir}/workflows/psychencode_reanalysis"

// Entry point
workflow {
    PSYCHENCODE_REANALYSIS()
}

// Completion handler
workflow.onComplete {
    log.info """
=========================================
 Pipeline Complete
=========================================
 Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
 Duration  : ${workflow.duration}
 Output    : ${params.outdir}
=========================================
"""
}

// Error handler
workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
