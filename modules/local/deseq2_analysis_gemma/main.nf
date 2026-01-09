process DESEQ2_ANALYSIS_GEMMA {
    tag "${pavlab_cell_type}"
    label 'process_medium'

    publishDir "${params.outdir}/${params.level}/DESeq2/gemma/${pavlab_cell_type}", mode: 'copy'

    input:
    tuple val(pavlab_cell_type), path(pseudobulk_matrix)

    output:
    tuple val(pavlab_cell_type), path("**results.tsv"), emit: all_contrasts_gemma
    path "**png", optional: true

    script:
    def filter_opt = params.filter_genes ? "--filter_genes" : ""
    """
    Rscript ${projectDir}/bin/DESeq2_analysis.R \
        --pseudobulk_matrix ${pseudobulk_matrix} \
        --metadata ${params.gemma_meta_dir} \
        --mode gemma \
        --cell_type ${pavlab_cell_type} \
        ${filter_opt}
    """
}
