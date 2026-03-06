/*
 * GEMMA DE Analysis Subworkflow
 *
 * Fetches pseudobulk data from GEMMA, aggregates by cell type,
 * and runs DESeq2 differential expression analysis.
 */

include { GET_GEMMA_PSEUDOBULKS } from "${projectDir}/modules/local/get_gemma_pseudobulks/main"
include { AGGREGATE_CELLTYPES_GEMMA } from "${projectDir}/modules/local/aggregate_celltypes_gemma/main"
include { DESEQ2_ANALYSIS_GEMMA } from "${projectDir}/modules/local/deseq2_analysis_gemma/main"

workflow GEMMA_DE_ANALYSIS {
    take:
    study_names      // Channel: study name strings
    author_submitted // val: true or false

    main:
    author_submitted_ch = Channel.value(author_submitted)

    // Step 1: Get pseudobulk matrices for each study from GEMMA
    GET_GEMMA_PSEUDOBULKS(study_names, author_submitted_ch)

    // Step 2: Collect all pseudobulks and aggregate by cell type
    aggregated_data = GET_GEMMA_PSEUDOBULKS.out.aggregated_data.collect()
    AGGREGATE_CELLTYPES_GEMMA(aggregated_data, author_submitted_ch)

    // Step 3: Prepare channel with cell type labels (include author_submitted in tuple)
    celltypes_ch = AGGREGATE_CELLTYPES_GEMMA.out.aggregated_celltypes
        .flatMap()
        .map { file ->
            def cell_type = file.getBaseName().split("_pseudobulk_matrix.tsv")[0]
            tuple(cell_type, file)
        }
        .combine(author_submitted_ch)
        .map { cell_type, file, as_val -> tuple(as_val, cell_type, file) }

    // Step 4: Run DESeq2 for each cell type
    DESEQ2_ANALYSIS_GEMMA(celltypes_ch)

    emit:
    deseq_results = DESEQ2_ANALYSIS_GEMMA.out.all_contrasts_gemma
    pseudobulks = AGGREGATE_CELLTYPES_GEMMA.out.aggregated_celltypes
}
