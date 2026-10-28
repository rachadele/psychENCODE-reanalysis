process LIMMA_VOOM_ANALYSIS_GEMMA {
    tag "${pavlab_cell_type}"
    label 'process_medium'

    publishDir { "${params.outdir}/author_submitted_${author_submitted}/${params.level}/limma_voom/gemma/${pavlab_cell_type}" }, mode: 'copy'

    input:
    tuple val(author_submitted), val(pavlab_cell_type), path(pseudobulk_matrix), val(gene_ncells), val(cell_counts)

    output:
    tuple val(author_submitted), val(pavlab_cell_type), path("**results.tsv"), emit: all_contrasts_gemma
    path "**png", optional: true

    script:
    // gene_ncells / cell_counts: absolute paths or "" (see GEMMA_DE_ANALYSIS). Gemma getSingleCellDataMatrix
    // output has neither: library size then falls back to column sums and the cell-count filters are skipped.
    def nc_opt = gene_ncells ? "--gene_ncells ${gene_ncells}" : ""
    def cc_opt = cell_counts ? "--cell_counts ${cell_counts}" : ""
    """
    Rscript ${projectDir}/bin/limma_voom_analysis.R \
        --pseudobulk_matrix ${pseudobulk_matrix} \
        --metadata ${params.gemma_meta_dir} \
        --cell_type ${pavlab_cell_type} \
        --min_cells_per_sample ${params.limma_min_cells_per_sample} \
        --min_cells_per_gene ${params.limma_min_cells_per_gene} \
        --min_variance ${params.limma_min_variance} \
        ${cc_opt} ${nc_opt}
    """
}
