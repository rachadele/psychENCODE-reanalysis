/*
 * GEMMA Comparison Subworkflow
 *
 * Compares DE results between two annotation schemes (e.g., pavlab vs author labels).
 * Calculates pairwise correlations and Jaccard overlap metrics.
 */

include { AGGREGATE_PAIRWISE } from "${projectDir}/modules/local/aggregate_pairwise/main"
include { DE_OVERLAP } from "${projectDir}/modules/local/de_overlap/main"
include { AVERAGE_DE_OVERLAPS } from "${projectDir}/modules/local/average_de_overlaps/main"
include { AVERAGE_DE_CORRELATIONS } from "${projectDir}/modules/local/average_de_correlations/main"

workflow GEMMA_COMPARISON {
    take:
    pavlab_results    // Channel: path to pavlab DESeq2 results
    author_results    // Channel: path to author-label DESeq2 results

    main:
    // Step 1: Group pavlab results by contrast
    pavlab_by_contrast = pavlab_results
        .map { file ->
            def parts = file.toString().split("/")
            def cell_type = parts[-3]
            def contrast = parts[-2]
            tuple(contrast, file)
        }
        .groupTuple(by: 0)

    // Step 2: Group author results by contrast (with cell type mapping)
    author_by_contrast = author_results
        .map { file ->
            def parts = file.toString().split("/")
            def cell_type = parts[-3]
            def pavlab_cell_type = params.gemma_to_gemma_map[cell_type] ?: cell_type
            def contrast = parts[-2]
            tuple(contrast, file)
        }
        .groupTuple(by: 0)

    // Step 3: Combine channels by contrast for pairwise analysis
    pairwise_channel = pavlab_by_contrast
        .combine(author_by_contrast, by: 0)

    // Step 4: Calculate pairwise correlations
    AGGREGATE_PAIRWISE(pairwise_channel)

    // Step 5: Average correlations across contrasts
    corr_results = AGGREGATE_PAIRWISE.out.pairwise_corrs.collect()
    AVERAGE_DE_CORRELATIONS(corr_results)

    // Step 6: Calculate DE overlap (Jaccard index)
    DE_OVERLAP(pairwise_channel)

    // Step 7: Average overlaps across contrasts
    overlap_results = DE_OVERLAP.out.overlap_results.collect()
    AVERAGE_DE_OVERLAPS(overlap_results)

    emit:
    avg_correlations = AVERAGE_DE_CORRELATIONS.out.average_corr
    avg_overlaps = AVERAGE_DE_OVERLAPS.out.average_jaccard
}
