process AGGREGATE_DATA_MANUAL {
    tag "${experiment}"
    label 'process_medium'

    publishDir "${params.outdir}/${params.level}/experiment_pseudobulks/manual/${experiment}", mode: 'copy'

    input:
    tuple val(experiment), path(annotation_file), path(h5ad_file)

    output:
    path "**pseudobulk.h5ad", emit: aggregated_experiments

    script:
    def filter_opt = params.filter_samples ? '--filter_samples' : ''
    """
    python ${projectDir}/bin/aggregate_data_manual.py \
        --h5ad_file ${h5ad_file} \
        --cell_type_column ${params.cell_type_column} \
        --annotation_file ${annotation_file} \
        ${filter_opt}
    """
}
