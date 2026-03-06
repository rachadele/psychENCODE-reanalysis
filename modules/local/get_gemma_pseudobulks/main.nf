process GET_GEMMA_PSEUDOBULKS {
    tag "${experiment}"
    label 'process_low'

    publishDir { "${params.outdir}/author_submitted_${author_submitted}/${params.level}/experiment_pseudobulks/gemma/${experiment}" }, mode: 'copy'

    input:
    val experiment
    val author_submitted

    output:
    path("**.tsv.gz"), emit: aggregated_data

    script:
    def cta_option = author_submitted ? "author-submitted" : "${params.cta_name}"
    """
    gemma-cli-staging getSingleCellDataMatrix -e ${experiment} --aggregate-by-assay \
        --aggregate-by-cell-type-assignment ${cta_option} \
        --output-file "${experiment}_aggregated_gemma.tsv.gz"
    """
}
