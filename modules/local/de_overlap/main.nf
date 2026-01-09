process DE_OVERLAP {
    tag "${contrast}"
    label 'process_low'

    publishDir "${params.outdir}/${params.level}/de_overlap/contrast_overlaps/${contrast}", mode: 'copy'

    input:
    tuple val(contrast), path(pavlab_files), path(author_files)

    output:
    path("**pairwise_jaccard.tsv"), emit: overlap_results

    script:
    """
    python ${projectDir}/bin/DE_overlap.py \
        --contrast ${contrast} \
        --sc_pipeline_paths ${pavlab_files} \
        --author_paths ${author_files}
    """
}
