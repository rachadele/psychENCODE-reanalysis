process AVERAGE_DE_CORRELATIONS {
    tag "average_correlations"
    label 'process_low'

    publishDir "${params.outdir}/${params.level}/average_corr/filtered_per_contrast", mode: 'copy', pattern: "filtered**.tsv"
    publishDir "${params.outdir}/${params.level}/average_corr/per_celltype", mode: 'copy', pattern: "**per_celltype_**.tsv"
    publishDir "${params.outdir}/${params.level}/average_corr/overall_average", mode: 'copy', pattern: "**average_corr*.tsv"
    publishDir "${params.outdir}/${params.level}/average_corr/figs", mode: 'copy', pattern: "**png"

    input:
    path corr_tables

    output:
    path "filtered**.tsv", optional: true
    path "**per_celltype_**.tsv", optional: true
    path "**average_corr*.tsv", emit: average_corr
    path "**png", optional: true

    script:
    """
    python ${projectDir}/bin/average_corr.py \
        --corr_tables ${corr_tables} \
        --annotation_level ${params.level}
    """
}
