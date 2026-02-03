process AGGREGATE_CELLTYPES_MANUAL {
    tag "aggregate_celltypes_manual"
    label 'process_medium'

    publishDir "${params.outdir}/author_submitted_${params.author_submitted}/${params.level}/ct_pseudobulks/manual", mode: 'copy'

    input:
    path h5ad_files

    output:
    path "**pseudobulk_matrix.tsv.gz", emit: aggregated_celltypes
    path "**pseudobulk_metadata.tsv", emit: aggregated_celltypes_meta

    script:
    """
    python ${projectDir}/bin/aggregate_celltypes_manual.py \
        --h5ad_files ${h5ad_files} \
        --cell_type_column ${params.cell_type_column}
    """
}
