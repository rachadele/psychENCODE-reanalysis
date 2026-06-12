process AGGREGATE_PAIRWISE {
    tag "${contrast}"
    label 'process_medium'

    publishDir "${params.outdir}/comparison/${params.level}/all_corr/${contrast}/figs/", mode: 'copy', pattern: "**png"
    publishDir "${params.outdir}/comparison/${params.level}/all_corr/${contrast}/files/", mode: 'copy', pattern: "**tsv"

    input:
    tuple val(contrast), path(pavlab_files), path(author_files)

    output:
    path "**pairwise_corr*.png", optional: true, emit: pairwise_heatmaps
    path "**scatter.png", optional: true, emit: scatter_plots
    path "**pairwise_corr*.tsv", emit: pairwise_corrs

    script:
    """
    python ${projectDir}/bin/aggregate_pairwise.py \
        --contrast ${contrast} \
        --pavlab_paths ${pavlab_files} \
        --author_paths ${author_files} \
        --params ${params.params_file}
    """
}
