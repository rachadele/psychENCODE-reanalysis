process DESEQ2_ANALYSIS_MANUAL {
    tag "${cell_type}"
    label 'process_medium'

    publishDir "${params.outdir}/${params.level}/DESeq2/manual/${cell_type}", mode: 'copy'

    input:
    tuple val(cell_type), path(pseudobulk_matrix), path(pseudobulk_metadata)

    output:
    tuple val(cell_type), path("**results.tsv"), emit: all_contrasts_manual
    path "**png", optional: true

    script:
    def filter_opt = params.filter_genes ? "--filter_genes" : ""
    """
    Rscript ${projectDir}/bin/DESeq2_analysis.R \
        --pseudobulk_matrix ${pseudobulk_matrix} \
        --metadata ${pseudobulk_metadata} \
        --cell_type ${cell_type} \
        --mode manual \
        ${filter_opt}
    """
}
