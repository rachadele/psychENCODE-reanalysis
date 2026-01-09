/*
 * Manual DE Analysis Subworkflow
 *
 * Aggregates pseudobulk data from H5AD files using custom cell type annotations,
 * and runs DESeq2 differential expression analysis.
 */

include { AGGREGATE_DATA_MANUAL } from '../../modules/local/aggregate_data_manual/main'
include { AGGREGATE_CELLTYPES_MANUAL } from '../../modules/local/aggregate_celltypes_manual/main'
include { DESEQ2_ANALYSIS_MANUAL } from '../../modules/local/deseq2_analysis_manual/main'

workflow MANUAL_DE_ANALYSIS {
    take:
    h5ad_files_channel  // Channel: tuple(experiment, annotation_file, h5ad_file)

    main:
    // Step 1: Aggregate each experiment's data manually
    AGGREGATE_DATA_MANUAL(h5ad_files_channel)

    // Step 2: Collect all aggregated experiments and combine by cell type
    aggregated_experiments = AGGREGATE_DATA_MANUAL.out.aggregated_experiments.collect()
    AGGREGATE_CELLTYPES_MANUAL(aggregated_experiments)

    // Step 3: Prepare cell type channel from pseudobulk matrices
    aggregated_celltypes = AGGREGATE_CELLTYPES_MANUAL.out.aggregated_celltypes
        .flatMap()
        .map { file ->
            def cell_type = file.getBaseName().split("_pseudobulk_matrix.tsv")[0]
            tuple(cell_type, file)
        }

    // Step 4: Prepare metadata channel
    aggregated_celltypes_meta = AGGREGATE_CELLTYPES_MANUAL.out.aggregated_celltypes_meta
        .flatMap()
        .map { file ->
            def cell_type = file.getBaseName().split("_pseudobulk_metadata")[0]
            tuple(cell_type, file)
        }

    // Step 5: Combine pseudobulk matrices with their metadata
    ct_pseudobulks_meta = aggregated_celltypes
        .combine(aggregated_celltypes_meta, by: 0)

    // Step 6: Run DESeq2 for each cell type
    DESEQ2_ANALYSIS_MANUAL(ct_pseudobulks_meta)

    emit:
    deseq_results = DESEQ2_ANALYSIS_MANUAL.out.all_contrasts_manual
    pseudobulks = AGGREGATE_CELLTYPES_MANUAL.out.aggregated_celltypes
    metadata = AGGREGATE_CELLTYPES_MANUAL.out.aggregated_celltypes_meta
}
