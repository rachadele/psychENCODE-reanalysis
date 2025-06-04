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
    study_names : ${params.study_names}
	  h5ad_files : ${params.h5ad_files}
    outdir: ${params.outdir}
	  EOF
    """
}

process get_gemma_pseudobulks {
  publishDir "${params.outdir}/aggregated/gemma/${experiment}", mode: "copy"

  input:
  val experiment

  output:
  path("**.tsv.gz"), emit: aggregated_data

  script:
  """
  bash $projectDir/bin/aggregateData.sh $experiment
  """
}

process aggregate_celltypes_gemma {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/pseudobulks/gemma", mode: "copy"

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
  publishDir "${params.outdir}/aggregated/manual/${experiment}", mode: "copy"

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
  publishDir "${params.outdir}/pseudobulks/manual", mode: "copy"

  input:
  path h5ad_files

  output:
  path "**pseudobulk_matrix.tsv.gz", emit: aggregated_celltypes

  script:
  """
  python $projectDir/bin/aggregate_celltypes_manual.py \\
        --h5ad_files ${h5ad_files}
  """
}

process DESeq2_analysis_gemma {
  conda "/home/rschwartz/anaconda3/envs/r4.3"
  publishDir "${params.outdir}/DESeq2/gemma/${cell_type}", mode: "copy"

  input:
  tuple val(cell_type), path(pseudobulk_matrix)

  output:
  path "**tsv"
  path "**png"

  script:
  """
  Rscript $projectDir/bin/DESeq2_analysis.R --pseudobulk_matrix ${pseudobulk_matrix} \\
        --gemma_metadata ${params.gemma_meta_dir}
  """
}

process DESeq2_analysis_manual {
  conda "/home/rschwartz/anaconda3/envs/r4.3"
  publishDir "${params.outdir}/DESeq2/manual/${cell_type}", mode: "copy"

  input:
  tuple val(cell_type), path(pseudobulk_matrix)

  output:
  path "**tsv"
  path "**png"

  script:
  """
  Rscript $projectDir/bin/DESeq2_analysis_manual.R --pseudobulk_matrix ${pseudobulk_matrix}
  """
}

workflow {
	// Save parameters to a file
	save_params_to_file()

  if (params.from_gemma) {
    Channel
      .fromPath(params.study_names)
      .flatMap { file -> file.readLines().collect { it.trim() } }
      .set { study_names }
    // Aggregate data from GEMMA
    get_gemma_pseudobulks(study_names)
    .set { aggregated_data_channel }

    aggregated_data_channel.collect()
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
    aggregated_celltypes_channel.view()
    // Run DESeq2 analysis
    DESeq2_analysis(aggregated_celltypes_channel)
 
  } else {
    Channel.fromPath("${params.h5ad_files}/*.h5ad").map { h5ad_file ->
        def name = h5ad_file.getBaseName()
        [name, h5ad_file]
    }
    .set { h5ad_files_channel }

    h5ad_files_channel.view()

    aggregate_data_manual(h5ad_files_channel).collect()
    .set { aggregated_experiments_channel }

    aggregate_celltypes_manual(aggregated_experiments_channel).flatMap()
    .set { aggregated_celltypes }

    // extract cell type from channel
    aggregated_celltypes.map { it ->
      def cell_type = it.getBaseName().split("_pseudobulk_matrix.tsv")[0] // e.g., "Astrocyte_pseudobulk_matrix.tsv.gz"
      [cell_type, it]
    }
    .set { aggregated_celltypes_channel }

    DESeq2_analysis_manual(aggregated_celltypes_channel)
  }

  // Combine RDS files into a single RDS file
  // combine_rds(rds_files_channel)


}