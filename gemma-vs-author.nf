#!/usr/bin/env nextflow


process save_params_to_file {
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

process wrangle_author {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/author_contrasts/${contrast}/files", mode: "copy", pattern: "**tsv"
  publishDir "${params.outdir}/author_contrasts/${contrast}/figs", mode: "copy", pattern: "**png"

  input:
  tuple val(contrast), path(author_results)

  output:
  tuple val(contrast), path("**tsv"), emit: ct_contrasts
  path "**png"

  script:
  """
  python $projectDir/bin/wrangle_author_de_results.py \\
        --author_degs ${author_results} \\
        --contrast ${contrast}
  """
}


process get_relevant_samples {
  publishDir "${params.outdir}/relevant_samples", mode: "copy"

  input:
  tuple val(pavlab_ct), val(author_ct), path(pseudobulk_matrix)

  output:
  tuple val(author_ct), path("**relevant_samples.tsv"), emit: relevant_samples

  script:
  // extract columns of bed.gz

  """
  zcat ${pseudobulk_matrix} | head -1 |
  cut -f 7- | tr '\\t' '\\n' > ${author_ct}_relevant_samples.tsv
  """
}

process get_gemma_pseudobulks {
  publishDir "${params.outdir}/experiment_pseudobulks/gemma/${experiment}", mode: "copy"

  input:
  val experiment

  output:
  path("**.tsv.gz"), emit: aggregated_data

  script:
  """
  if [ ${params.author_submitted} = true ]; then

    gemma-cli getSingleCellDataMatrix -e $experiment --aggregate-by-assay \\
    --aggregate-by-cell-type-assignment author-submitted \\
    --output-file "${experiment}_aggregated_gemma.tsv.gz"
  
  else

    gemma-cli getSingleCellDataMatrix -e $experiment --aggregate-by-assay \\
    --aggregate-by-preferred-cell-type-assignment \\
    --output-file "${experiment}_aggregated_gemma.tsv.gz"
  fi

  """
}

process aggregate_celltypes_gemma {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/ct_pseudobulks/gemma", mode: "copy"

  input:
  path pseudobulk_matrices
  
  output:
  path "**pseudobulk_matrix.tsv.gz", emit: aggregated_celltypes
  
  script:
  """
  python $projectDir/bin/aggregate_celltypes_gemma.py \\
        --pseudobulk_matrices ${pseudobulk_matrices} \\
        --metadata_files ${params.gemma_meta_dir}
  """
  
}

process aggregate_data_manual {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/experiment_pseudobulks/manual/${experiment}", mode: "copy"

  input:
  tuple val(experiment), path(h5ad_file)

  output:
  path "**pseudobulk.h5ad", emit: aggregated_experiments

  script:
  """
  python $projectDir/bin/aggregate_data_manual.py \\
        --h5ad_file ${h5ad_file} \\
  """
}

process aggregate_celltypes_manual {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/ct_pseudobulks/manual", mode: "copy"

  input:
  path h5ad_files

  output:
  path "**pseudobulk_matrix.tsv.gz", emit: aggregated_celltypes
  path "**pseudobulk_metadata.tsv", emit: aggregated_celltypes_meta

  script:
  """
  python $projectDir/bin/aggregate_celltypes_manual.py \\
        --h5ad_files ${h5ad_files}
  """
}

process DESeq2_analysis_gemma {
  conda "/home/rschwartz/anaconda3/envs/r4.3"
  publishDir "${params.outdir}/DESeq2/gemma/${pavlab_cell_type}", mode: "copy"

  input:
  //tuple val(pavlab_cell_type), path(pseudobulk_matrix), val(author_cell_type), path(revelant_samples_file)
  tuple val(pavlab_cell_type), path(pseudobulk_matrix)


  output:
  tuple val(pavlab_cell_type), path("**results.tsv"), emit: all_contrasts_gemma
  path "**png"

  script:
  """
  Rscript $projectDir/bin/DESeq2_analysis.R --pseudobulk_matrix ${pseudobulk_matrix} \\
        --metadata ${params.gemma_meta_dir} \\
        --mode gemma \\
        --cell_type ${pavlab_cell_type}
  """
}

process DESeq2_analysis_manual {
  conda "/home/rschwartz/anaconda3/envs/r4.3"
  publishDir "${params.outdir}/DESeq2/manual/${cell_type}", mode: "copy"

  input:
  tuple val(cell_type), path(pseudobulk_matrix), path(pseudobulk_metadata)

  output:
  tuple val(cell_type), path("**results.tsv"), emit: all_contrasts_manual
  path "**png"

  script:
  """
  Rscript $projectDir/bin/DESeq2_analysis.R --pseudobulk_matrix ${pseudobulk_matrix} --metadata ${pseudobulk_metadata} \\
        --cell_type ${cell_type} \\
        --mode manual

  """
}

process DE_corr {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/DE_corr/${mode}/${contrast}/figs/${pavlab_ct}/${author_ct}", mode: 'copy', pattern: '**png'
  publishDir "${params.outdir}/DE_corr/${mode}/${contrast}/files/${pavlab_ct}/${author_ct}", mode: 'copy', pattern: '**tsv'


  input:
  val mode 
  tuple val(contrast), val(pavlab_ct), path(pavlab_results), val(author_ct), val(author_results)

  output:
  path "**png"
  path "missing_genes.tsv", emit: missing_genes
  path "**DE_merge**.tsv", emit: merged_results

  script:
  """
  python $projectDir/bin/DE_corr.py \\
        --pavlab_results ${pavlab_results} \\
        --author_results ${author_results} \\
        --contrast ${contrast} \\
        --mode ${mode}
  """
}

process aggregate_pairwise {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/all_corr/${contrast}/figs/", mode: "copy", pattern: "**png"
  publishDir "${params.outdir}/all_corr/${contrast}/files/", mode: "copy", pattern: "**tsv"

  input:
  tuple val(contrast), path(pavlab_files), path(author_files)

  output:
  path "pairwise_corr**png"
  path "pairwise_corr**tsv"

  script:

  """
  python $projectDir/bin/aggregate_pairwise.py \\
      --contrast ${contrast} \\
      --pavlab_paths ${pavlab_files} \\
      --author_paths ${author_files}
  """
}


process compare_pseudobulks {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/pseudobulk_comparisons/${mode}/files/${pavlab_cell_type}/${author_cell_type}", mode: "copy", pattern: "**tsv"
  publishDir "${params.outdir}/pseudobulk_comparisons/${mode}/figs/${pavlab_cell_type}/${author_cell_type}", mode: "copy", pattern: "**png"

  input:
  val(mode)
  tuple val(pavlab_cell_type), val(author_cell_type), path(author_pseudobulk),  path(pavlab_pseudobulk)

  output:
  path "**tsv", emit: comparison_results
  path "**png"

  script:
  """
  python $projectDir/bin/compare_pseudobulks.py \\
        --author_pseudobulk ${author_pseudobulk} \\
        --author_cell_type ${author_cell_type} \\
        --pavlab_pseudobulk ${pavlab_pseudobulk} \\
        --pavlab_cell_type ${pavlab_cell_type} \\
        --gemma_metadata ${params.gemma_meta_dir} \\
        --mode ${mode}
  """
}


workflow {

  def mode
  if (params.from_gemma) {
    mode = "gemma"
  } else {
    mode = "manual"
  }

	// Save parameters to a file
	save_params_to_file()

  author_results = Channel.fromPath(params.author_results)

  author_results.map { file ->
    def contrast = file.getBaseName().split("_")[0] // e.g., "ASD_DEGcombined.csv"
    [contrast, file]
  }.set{ author_contrasts }
  
  wrangle_author(author_contrasts)
 

  wrangle_author.out.ct_contrasts.flatMap{it ->
   def contrast = it[0]
   def files = it[1]
    files.collect { file ->
      def author_cell_type = file.getBaseName().split("_")[1] // e.g., "Bipolar_Vip_degs.tsv
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

  get_relevant_samples(author_pseudobulks_channel)



  get_relevant_samples.out.relevant_samples.map { it ->
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
    get_gemma_pseudobulks(study_names)
    get_gemma_pseudobulks.out.aggregated_data.collect()
    .set { aggregated_data }
    aggregate_celltypes_gemma(aggregated_data)
    aggregate_celltypes_gemma.out.aggregated_celltypes.flatMap()
    .set { aggregated_celltypes } 

    // extract cell type from channel
    aggregated_celltypes.map { it ->
      def cell_type = it.getBaseName().split("_pseudobulk_matrix.tsv")[0] // e.g., "Astrocyte_pseudobulk_matrix.tsv.gz"
      [cell_type, it]
    }
    .set { aggregated_celltypes_channel }
    // Run DESeq2 analysis
    DESeq2_analysis_gemma(aggregated_celltypes_channel)

    DESeq2_analysis_gemma.out.all_contrasts_gemma.flatMap { it ->
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


    aggregate_data_manual(h5ad_files_channel).collect()
    .set { aggregated_experiments_channel }

    aggregate_celltypes_manual(aggregated_experiments_channel)
    
    aggregate_celltypes_manual.out.aggregated_celltypes
    .flatMap()
    .set { aggregated_celltypes }

    aggregate_celltypes_manual.out.aggregated_celltypes_meta
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

    

    DESeq2_analysis_manual(ct_pseudobulks_meta_channel) 
    // flatMap results

    DESeq2_analysis_manual.out.all_contrasts_manual.flatMap { it ->
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
  DE_corr(mode, all_contrasts_channel)
 
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
  aggregate_pairwise(pairwise_channel)
  //// Compare pseudobulks
  //// only compare for valid cell type mappings
  //// need a mapping dictionary of author cell types to pavlab cell types since strings don't map exactly
  author_pseudobulks_channel.combine(aggregated_celltypes_channel, by: 0)
  .set { pseudobulks_combined }

  compare_pseudobulks(mode, pseudobulks_combined)


}
