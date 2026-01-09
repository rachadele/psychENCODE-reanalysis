process AGGREGATE_CELLTYPES_GEMMA {
    tag "AGGREGATE_CELLTYPES_GEMMA"
    label 'process_medium'

    publishDir "${params.outdir}/${params.level}/ct_pseudobulks/gemma", mode: 'copy'

    input:
    path pseudobulk_matrices

    output:
    path "**pseudobulk_matrix.tsv.gz", emit: aggregated_celltypes

    script:
    """
    python ${projectDir}/bin/aggregate_celltypes_gemma.py \
        --pseudobulk_matrices ${pseudobulk_matrices} \
        --metadata_files ${params.gemma_meta_dir}
    """
}
