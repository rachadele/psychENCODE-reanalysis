#!/usr/bin/env nextflow


 save_params_to_file {
    publishDir (
        "${params.outdir}",
        mode: "copy"
    )

    output:
    file "params.yaml"

    script:
    """
    cat <<EOF > params.yaml
    from_gemma : ${params.from_gemma}
    gemma_meta_dir : ${params.gemma_meta_dir}
	  h5ad_files : ${params.h5ad_files}
    outdir: ${params.outdir}
    EOF
    """
}


    // Import all processes from modules
  include { WrangleAuthor } from "${projectDir}/modules/wrangle_author.nf"
  include { GetRelevantSamples } from "${projectDir}/modules/get_relevant_samples.nf"
  include { GetGemmaPseudobulks } from "${projectDir}/modules/get_gemma_pseudobulks.nf"
  include { AggregateCelltypesGemma } from "${projectDir}/modules/aggregate_celltypes_gemma.nf"
  include { AggregateDataManual } from "${projectDir}/modules/aggregate_data_manual.nf"
  include { AggregateCelltypesManual } from "${projectDir}/modules/aggregate_celltypes_manual.nf"
  include { DESeq2AnalysisGemma } from "${projectDir}/modules/DESeq2_analysis_gemma.nf"
  include { DESeq2AnalysisManual } from "${projectDir}/modules/DESeq2_analysis_manual.nf"
  include { DECorr } from "${projectDir}/modules/DE_corr.nf"
  include { AggregatePairwise } from "${projectDir}/modules/aggregate_pairwise.nf"
  include { ComparePseudobulks } from "${projectDir}/modules/compare_pseudobulks.nf"


workflow {

  author_results.map { file ->
    def contrast = file.getBaseName().split("_")[0] // e.g., "ASD_DEGcombined.csv"
    [contrast, file]
  }.set{ author_contrasts }
  
  wrangle_author(author_contrasts)
 

  WrangleAuthor.out.ct_contrasts.flatMap{it ->
   def contrast = it[0]
   def files = it[1]
    files.collect { file ->
      def author_cell_type = file.getBaseName().split("_")[1] // e.g., "Bipolar_Vip_degs.tsv"
      if (params.author_submitted) {
        pavlab_cell_type = params.author_ct_map[author_cell_type] //?: cell_type.replace(".", "_")
      }
      else {
        pavlab_cell_type = params.gemma_ct_map[author_cell_type] //?: cell_type.replace(".", "_")
      }
      [contrast, pavlab_cell_type, author_cell_type, file]
    }
  }
  .set { all_contrasts_author_ct }


  author_pseudobulks = Channel.fromPath(params.author_pseudobulks)
  author_pseudobulks.map { file ->
    def author_cell_type = file.getBaseName().split(".expr.bed")[0].split("__")[0] // e.g., "Astrocyte_pseudobulk_matrix.tsv.gz"
    //def author_cell_type = author_cell_type.split(".")[0] // e.g., "Astrocyte"
    if (params.author_submitted) {
      pavlab_cell_type = params.author_ct_map[author_cell_type]
    } else {
      pavlab_cell_type = params.gemma_ct_map[author_cell_type]
    }
    [pavlab_cell_type, author_cell_type, file]
  }
  .set { author_pseudobulks_channel }

  GetRelevantSamples(author_pseudobulks_channel)


  GetRelevantSamples.out.relevant_samples.map { it ->
    def author_cell_type = it[0]
    def relevant_samples_file = it[1]
    def pavlab_cell_type = params.gemma_ct_map[author_cell_type] ?: author_cell_type.replace(".", "_")
    [pavlab_cell_type, author_cell_type, relevant_samples_file]
  }
  .set { relevant_samples_channel }

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
      def cell_type = it.getBaseName().split("_pseudobulk_matrix.tsv")[0] // e.g., "Astrocyte_pseudobulk_matrix.tsv.gz"
      [cell_type, it]
    }
    .set { aggregated_celltypes_channel }
    // Run DESeq2 analysis
    DESeq2AnalysisGemma(aggregated_celltypes_channel)

    DESeq2AnalysisGemma.out.all_contrasts_gemma.flatMap { it ->
      def cell_type = it[0]
      def files = it[1]
      files.collect { results_file ->
        def contrast = results_file.getParent().getBaseName() // e.g., Disorder_PTSD_vs_Control
        [contrast, cell_type, results_file]
      }
    }
    .set { all_contrasts_pavlab_ct }

    // combine manual contrasts and author contrasts
    
    all_contrasts_pavlab_ct.map { full_contrast, ct, file ->
        def contrast = full_contrast.replaceAll(/Disorder_|_vs_Control/, '')
        
        [contrast, ct, file]
    }.set { pavlab_contrast_channel }

    pavlab_contrast_channel.combine(all_contrasts_author_ct, by: [0,1])
    .set { all_contrasts_channel }
 
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
      def cell_type = it.getBaseName().split("_pseudobulk_matrix.tsv")[0] // e.g., "Astrocyte_pseudobulk_matrix.tsv.gz"
      [cell_type, it]
    }
    .set { aggregated_celltypes_channel }


    aggregated_celltypes_meta.map { it ->
      def cell_type = it.getBaseName().split("_pseudobulk_metadata")[0] // e.g., "Astrocyte_pseudobulk_metadata.tsv"
      [cell_type, it]
    }
    .set { aggregated_celltypes_meta_channel }

    aggregated_celltypes_channel.combine(aggregated_celltypes_meta_channel, by: 0)
    .set { ct_pseudobulks_meta_channel }

    DESeq2AnalysisManual(ct_pseudobulks_meta_channel) 
    // flatMap results

    DESeq2AnalysisManual.out.all_contrasts_manual.flatMap { it ->
      def cell_type = it[0]
      def files = it[1]
      files.collect { results_file ->
        def contrast = results_file.getParent().getBaseName() // e.g., Disorder_PTSD_vs_Control
        [contrast, cell_type, results_file]
      }
    }
    .set { all_contrasts_pavlab_ct }

    all_contrasts_pavlab_ct.view() 
    // combine manual contrasts and author contrasts
    
    all_contrasts_pavlab_ct.map { full_contrast, pavlab_ct, file ->
        def contrast = full_contrast.replaceAll(/Disorder_|_vs_Control/, '')
        
        [contrast, pavlab_ct, file]
    }.set { pavlab_contrast_channel }

    if (params.author_submitted) {
      all_contrasts_author_ct.map { contrast, pavlab_ct, author_ct, file ->
        [contrast, author_ct, pavlab_ct, file ]}
    }.set { all_contrasts_author_ct }
    
    pavlab_contrast_channel.combine(all_contrasts_author_ct, by: [0,1])
    .set { all_contrasts_channel }
  }

  // Run DE correlation
  DECorr(mode, all_contrasts_channel)
 
  pavlab_contrast_channel.map {contrast, ct, file ->
    [contrast, file]
  }.groupTuple(by: 0)
  .set { pavlab_files }


  all_contrasts_author_ct.map {contrast, pavlab_ct, author_ct, file ->
    [contrast, file]
  }.groupTuple(by: 0)
  .set { author_files }

  pavlab_files.combine(author_files, by: 0)
  .set { pairwise_channel }
  // Aggregate pairwise results
  AggregatePairwise(pairwise_channel)
  //// Compare pseudobulks
  //// only compare for valid cell type mappings
  //// need a mapping dictionary of author cell types to pavlab cell types since strings don't map exactly
  author_pseudobulks_channel.combine(aggregated_celltypes_channel, by: 0)
  .set { pseudobulks_combined }

  ComparePseudobulks(mode, pseudobulks_combined)

}


workflow.onComplete {
    println "Successfully completed"
    println ( workflow.success ? 
    """
    ===============================================================================
    Pipeline execution summary
    -------------------------------------------------------------------------------

    Run as      : ${workflow.commandLine}
    output dir : ${params.outdir}

    --------------------------------------------------------------------------------
    ================================================================================
    """.stripIndent() : """
    Failed: ${workflow.errorReport}
    exit status : ${workflow.exitStatus}
    """.stripIndent()
    )
}

workflow.onError = {
println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}