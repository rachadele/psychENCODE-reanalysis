/*
 * GEMMA DE Analysis Subworkflow
 *
 * Fetches pseudobulk data from GEMMA, aggregates by cell type,
 * and runs DESeq2 (params.de_method = 'deseq2', default) or limma-voom
 * (params.de_method = 'limma_voom', Gemma's method) differential expression analysis.
 */

include { GET_GEMMA_PSEUDOBULKS } from "${projectDir}/modules/local/get_gemma_pseudobulks/main"
include { AGGREGATE_CELLTYPES_GEMMA } from "${projectDir}/modules/local/aggregate_celltypes_gemma/main"
include { DESEQ2_ANALYSIS_GEMMA } from "${projectDir}/modules/local/deseq2_analysis_gemma/main"
include { LIMMA_VOOM_ANALYSIS_GEMMA } from "${projectDir}/modules/local/limma_voom_analysis_gemma/main"

workflow GEMMA_DE_ANALYSIS {
    take:
    study_names      // Channel: study name strings
    author_submitted // val: true or false

    main:
    author_submitted_ch = Channel.value(author_submitted)

    // Local cell-type pseudobulks (bin/build_local_pseudobulks.py + bin/aggregate_celltypes_gemma.py)
    // replace Steps 1-2 when params.pavlab_ct_pseudobulks or params.author_ct_pseudobulks is set: a
    // directory of <ct>/<ct>_pseudobulk_matrix.tsv.gz. An arm without its own directory is then skipped.
    def local_input = params.pavlab_ct_pseudobulks || params.author_ct_pseudobulks
    def ct_dir = author_submitted ? params.author_ct_pseudobulks : params.pavlab_ct_pseudobulks

    if (local_input) {
        ct_matrices = ct_dir ? Channel.fromPath("${ct_dir}/*/*_pseudobulk_matrix.tsv.gz", checkIfExists: true) : Channel.empty()
        pseudobulks = ct_matrices
    } else {
        // Step 1: Get pseudobulk matrices for each study from GEMMA
        GET_GEMMA_PSEUDOBULKS(study_names, author_submitted_ch)

        // Step 2: Collect all pseudobulks and aggregate by cell type
        aggregated_data = GET_GEMMA_PSEUDOBULKS.out.aggregated_data.collect()
        AGGREGATE_CELLTYPES_GEMMA(aggregated_data, author_submitted_ch)
        ct_matrices = AGGREGATE_CELLTYPES_GEMMA.out.aggregated_celltypes.flatMap()
        pseudobulks = AGGREGATE_CELLTYPES_GEMMA.out.aggregated_celltypes
    }

    // Step 3: Prepare channel with cell type labels (include author_submitted in tuple). For local input,
    // limma_voom also picks up <ct>/<ct>_gene_ncells.tsv.gz and <dir>/cell_counts.tsv when present
    // (expressing cells per gene; cells and library size per sample), as absolute paths.
    celltypes_ch = ct_matrices
        .map { mat ->
            def cell_type = mat.getName().split("_pseudobulk_matrix.tsv")[0]
            def ncells = ct_dir ? file("${ct_dir}/${cell_type}/${cell_type}_gene_ncells.tsv.gz") : null
            def counts = ct_dir ? file("${ct_dir}/cell_counts.tsv") : null
            tuple(author_submitted, cell_type, mat,
                  ncells && ncells.exists() ? ncells.toString() : "",
                  counts && counts.exists() ? counts.toString() : "")
        }

    // Step 4: Run DE for each cell type (same output layout for both methods)
    if (params.de_method == 'limma_voom') {
        LIMMA_VOOM_ANALYSIS_GEMMA(celltypes_ch)
        de_out = LIMMA_VOOM_ANALYSIS_GEMMA.out.all_contrasts_gemma
    } else if (params.de_method == 'deseq2') {
        DESEQ2_ANALYSIS_GEMMA(celltypes_ch.map { as_val, ct, f, nc, cc -> tuple(as_val, ct, f) })
        de_out = DESEQ2_ANALYSIS_GEMMA.out.all_contrasts_gemma
    } else {
        error "Unknown de_method '${params.de_method}': use deseq2 or limma_voom"
    }

    emit:
    deseq_results = de_out
    pseudobulks = pseudobulks
}
