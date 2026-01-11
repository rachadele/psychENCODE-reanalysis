process AVERAGE_DE_OVERLAPS {
    tag "average_overlaps"
    label 'process_low'

    publishDir "${params.outdir}/author_submitted_${params.author_submitted}/${params.level}/de_overlap/average_overlaps/filtered_per_contrast", mode: 'copy', pattern: "filtered**.tsv"
    publishDir "${params.outdir}/author_submitted_${params.author_submitted}/${params.level}/de_overlap/average_overlaps/per_celltype", mode: 'copy', pattern: "per_celltype_average*.tsv"
    publishDir "${params.outdir}/author_submitted_${params.author_submitted}/${params.level}/de_overlap/average_overlaps/overall_average", mode: 'copy', pattern: "average_jaccard.tsv"
    publishDir "${params.outdir}/author_submitted_${params.author_submitted}/${params.level}/de_overlap/average_overlaps/figs", mode: 'copy', pattern: "**png"

    input:
    path de_overlap_files

    output:
    path "filtered**.tsv", optional: true
    path "**per_celltype_average_jaccard.tsv", optional: true
    path "average_jaccard.tsv", emit: average_jaccard
    path "**png", optional: true

    script:
    """
    python ${projectDir}/bin/average_de_overlaps.py \
        --de_overlap_paths ${de_overlap_files} \
        --annotation_level ${params.level}
    """
}
