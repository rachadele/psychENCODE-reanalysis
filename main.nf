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
	  h5ad_files : ${params.h5ad_files}
    outdir: ${params.outdir}
	  EOF
    """
}

process h5ad_to_rds {
  publishDir "${params.outdir}/rds", mode: "copy"

	conda "/home/rschwartz/anaconda3/envs/r4.3"

	input:
	tuple val(name), path(h5ad_file)

	output:
	tuple val(name), path("${name}.rds"), emit: rds

	script:

	"""
	Rscript $projectDir/bin/h5ad_to_rds.R --h5ad_file ${h5ad_file}
	"""
}

process combine_rds {
  memory = '32 GB'
  publishDir "${params.outdir}/combined", mode: "copy"

  conda "/home/rschwartz/anaconda3/envs/r4.3"

  input:
  path rds_files

  output:
  path "combined.rds", emit: combined_rds

  script:
  """
  Rscript $projectDir/bin/combine_rds.R --rds_files ${rds_files}
  """
}

workflow {
	// Save parameters to a file
	save_params_to_file()

  Channel.fromPath("${params.h5ad_files}/*.h5ad").map { h5ad_file ->
      def name = h5ad_file.getBaseName()
      [name, h5ad_file]
  }
  .set { h5ad_files_channel }

  h5ad_files_channel.view()
	// Convert h5ad files to RDS format
	h5ad_to_rds(h5ad_files_channel)

  h5ad_to_rds.out.rds.collect().toList()
  .set { rds_files_channel }

  rds_files_channel.view()

  // Combine RDS files into a single RDS file
  //combine_rds(rds_files_channel)


}