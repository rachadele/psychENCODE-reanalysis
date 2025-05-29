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
  publishDir "${params.outdir}/aggregated", mode: "copy"

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
  publishDir "${params.outdir}/pseudobulks", mode: "copy"

  input:
  path pseudobulk_matrices
  
  output:
  path "**pseudobulk_matrix.tsv.gz", emit: aggregated_pseudobulks
  
  script:
  """
  python $projectDir/bin/aggregate_celltypes_gemma.py \\
        --pseudobulk_matrices ${pseudobulk_matrices} \\
        --metadata_files ${params.gemma_meta_dir}
  """
  
}

process h5ad_to_rds {
  publishDir "${params.outdir}/rds", mode: "copy"

	conda "/home/rschwartz/anaconda3/envs/r4.3"

	input:
	tuple val(name), path(h5ad_file)

	output:
	path "${name}.rds", emit: rds

	script:

	"""
	Rscript $projectDir/bin/h5ad_to_rds.R --h5ad_file ${h5ad_file}
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
    
  }

  else {
    Channel.fromPath("${params.h5ad_files}/*.h5ad").map { h5ad_file ->
        def name = h5ad_file.getBaseName()
        [name, h5ad_file]
    }
    .set { h5ad_files_channel }

    // Convert h5ad files to RDS format
    h5ad_to_rds(h5ad_files_channel)

    h5ad_to_rds.out.rds.collect()
    .set { rds_files_channel }

  }

  // Combine RDS files into a single RDS file
  // combine_rds(rds_files_channel)


}