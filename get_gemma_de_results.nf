#!/usr/bin/env nextflow

// Import required modules
include { GetGemmaPseudobulks } from "${projectDir}/modules/get_gemma_pseudobulks.nf"
include { AggregateCelltypesGemma } from "${projectDir}/modules/aggregate_celltypes_gemma.nf"
include { AggregateDataManual } from "${projectDir}/modules/aggregate_data_manual.nf"
include { AggregateCelltypesManual } from "${projectDir}/modules/aggregate_celltypes_manual.nf"
include { DESeq2AnalysisGemma } from "${projectDir}/modules/DESeq2_analysis_gemma.nf"
include { DESeq2AnalysisManual } from "${projectDir}/modules/DESeq2_analysis_manual.nf"

workflow {
  if (params.from_gemma) {
    Channel
      .fromPath(params.study_names)
      .flatMap { file -> file.readLines().collect { it.trim() } }
      .set { study_names }
    // Aggregate data from GEMMA
    GetGemmaPseudobulks(study_names)
    GetGemmaPseudobulks.out.aggregated_data.collect()
    .set { aggregated_data }
    AggregateCelltypesGemma(aggregated_data)
    AggregateCelltypesGemma.out.aggregated_celltypes.flatMap()
    .set { aggregated_celltypes }

    // extract cell type from channel
    aggregated_celltypes.map { it ->
      def cell_type = it.getBaseName().split("_pseudobulk_matrix.tsv")[0]
      [cell_type, it]
    }
    .set { aggregated_celltypes_channel }
    // Run DESeq2 analysis
    DESeq2AnalysisGemma(aggregated_celltypes_channel)

    DESeq2AnalysisGemma.out.all_contrasts_gemma.flatMap { it ->
      def cell_type = it[0]
      def files = it[1]
      files.collect { results_file ->
        def contrast = results_file.getParent().getBaseName()
        [contrast, cell_type, results_file]
      }
    }
    .set { all_contrasts_pavlab_ct }

  } else {
    Channel.fromPath("${params.h5ad_files}/*.h5ad").map { h5ad_file ->
        def name = h5ad_file.getBaseName()
        [name, h5ad_file]
    }
    .set { h5ad_files_channel }

    AggregateDataManual(h5ad_files_channel).collect()
    .set { aggregated_experiments_channel }

    AggregateCelltypesManual(aggregated_experiments_channel)
    AggregateCelltypesManual.out.aggregated_celltypes
    .flatMap()
    .set { aggregated_celltypes }

    AggregateCelltypesManual.out.aggregated_celltypes_meta
    .flatMap()
    .set { aggregated_celltypes_meta }

    // extract cell type from channel
    aggregated_celltypes.map { it ->
      def cell_type = it.getBaseName().split("_pseudobulk_matrix.tsv")[0]
      [cell_type, it]
    }
    .set { aggregated_celltypes_channel }

    aggregated_celltypes_meta.map { it ->
      def cell_type = it.getBaseName().split("_pseudobulk_metadata")[0]
      [cell_type, it]
    }
    .set { aggregated_celltypes_meta_channel }

    aggregated_celltypes_channel.combine(aggregated_celltypes_meta_channel, by: 0)
    .set { ct_pseudobulks_meta_channel }

    DESeq2AnalysisManual(ct_pseudobulks_meta_channel)
    // flatMap results
    //DESeq2AnalysisManual.out.all_contrasts_manual.flatMap { it ->
      //def cell_type = it[0]
      //def files = it[1]
      //files.collect { results_file ->
        //def contrast = results_file.getParent().getBaseName()
        //[contrast, cell_type, results_file]
      //}
    //}
    //.set { all_contrasts_pavlab_ct }
  //}
	}
}
workflow.onComplete {
  println "Workflow completed. Results are available in the output directory: ${params.outdir}"
}