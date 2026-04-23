process AVERAGE_DE_CORRELATIONS {
    tag "average_correlations"
    label 'process_low'

    publishDir "${params.outdir}/comparison/${params.level}/average_corr/filtered_per_contrast", mode: 'copy', pattern: "filtered**.tsv"
    publishDir "${params.outdir}/comparison/${params.level}/average_corr/per_celltype", mode: 'copy', pattern: "**per_celltype_**.tsv"
    publishDir "${params.outdir}/comparison/${params.level}/average_corr/overall_average", mode: 'copy', pattern: "**average_corr*.tsv"
    publishDir "${params.outdir}/comparison/${params.level}/average_corr/combined", mode: 'copy', pattern: "combined_filtered_corr.tsv"
    publishDir "${params.outdir}/comparison/${params.level}/average_corr/figs", mode: 'copy', pattern: "**png"
    publishDir "${params.outdir}/comparison/${params.level}/average_corr/figs", mode: 'copy', pattern: "heatmap_**.tsv"

    input:
    path corr_tables

    output:
    path "filtered**.tsv", optional: true
    path "**per_celltype_**.tsv", optional: true
    path "**average_corr*.tsv", emit: average_corr
    path "combined_filtered_corr.tsv", optional: true
    path "heatmap_**.tsv", optional: true
    path "**png", optional: true

    script:
    def f1_arg = params.f1_path ? "--f1_path ${params.f1_path}" : ""
    """
    python ${projectDir}/bin/average_corr.py \
        --corr_tables ${corr_tables} \
        --annotation_level ${params.level} \
        --params ${params.params_file} \
        ${f1_arg}
    """
}
