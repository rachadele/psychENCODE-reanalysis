/*
 * psychENCODE Reanalysis Main Workflow
 *
 * Orchestrates differential expression analysis and comparison of single-cell RNA-seq data.
 * Supports modular execution of individual stages.
 *
 * Stage 1: DE Analysis (GEMMA or Manual pathway)
 * Stage 2: Comparison of DE results between annotation schemes
 */

// Include subworkflows
include { GEMMA_DE_ANALYSIS } from "${projectDir}/subworkflows/local/gemma_de_analysis"
include { MANUAL_DE_ANALYSIS } from "${projectDir}/subworkflows/local/manual_de_analysis"
include { GEMMA_COMPARISON } from "${projectDir}/subworkflows/local/gemma_comparison"

workflow PSYCHENCODE_REANALYSIS {

    main:
    // Determine which stages to run
    def run_de = params.run_de_analysis && !params.skip_de_analysis
    def run_comparison = params.run_comparison && !params.skip_comparison

    // Initialize result channels
    de_results = Channel.empty()

    //=========================================================================
    // STAGE 1: Differential Expression Analysis
    //=========================================================================
    if (run_de) {
        if (params.from_gemma) {
            // GEMMA pathway: fetch data from GEMMA and run DESeq2
            study_names = Channel
                .fromPath(params.study_names)
                .flatMap { file -> file.readLines().collect { it.trim() } }

            GEMMA_DE_ANALYSIS(study_names)
            de_results = GEMMA_DE_ANALYSIS.out.deseq_results

        } else {
            // Manual pathway: use H5AD files with custom annotations
            h5ad_files_channel = Channel
                .fromPath("${params.h5ad_files}/*.h5ad")
                .map { h5ad_file ->
                    def name = h5ad_file.getBaseName()
                    def annotation_file = params.celltype_annotation_files[name]
                    tuple(name, annotation_file, h5ad_file)
                }

            MANUAL_DE_ANALYSIS(h5ad_files_channel)
            de_results = MANUAL_DE_ANALYSIS.out.deseq_results
        }
    }

    //=========================================================================
    // STAGE 2: Comparison of DE Results
    //=========================================================================
    if (run_comparison) {
        // Load pavlab results (from Stage 1 or pre-existing)
        if (run_de && params.author_submitted == false && params.from_gemma) {
            // Use results from Stage 1 - flatten the tuple to get just files
            // this only works if stage 1 was run with pavlab annotations
            pavlab_results = de_results.map { cell_type, files -> files }.flatten()
        } else if (params.pavlab_deseq_results) {
            // Use pre-existing results from params
            pavlab_results = Channel.fromPath(params.pavlab_deseq_results)
        } else {
            error "No pavlab DESeq2 results available. Either run DE analysis stage or provide --pavlab_deseq_results"
        }

        // Load author-label results
        if (run_de && params.author_submitted == true && params.from_gemma) {
            // Use results from Stage 1 - flatten the tuple to get just files
            // this only works if stage 1 was run with author-submitted annotations
            author_results = de_results.map { cell_type, files -> files }.flatten()
        } else
        if (params.author_label_deseq_results) {
            author_results = Channel.fromPath(params.author_label_deseq_results)
        } else {
            error "No author-label DESeq2 results provided. Please provide --author_label_deseq_results"
        }

        GEMMA_COMPARISON(pavlab_results, author_results)
    }

    emit:
    de_results = de_results
}
