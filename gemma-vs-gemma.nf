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
  python $projectDir/bin/aggregate_pairwise_gemma.py \\
      --contrast ${contrast} \\
      --sc_pipeline_paths ${pavlab_files} \\
      --author_paths ${author_files}
  """
}

workflow {

  pavlab_deseq_results = Channel
    .fromPath("/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results_author_submitted_false_from_gemma_true/DESeq2/gemma/**tsv")
	
  author_label_deseq_results = Channel
    .fromPath("/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results_author_submitted_true_from_gemma_true/DESeq2/gemma/**tsv") 
  
  // get cell type, contrast from each path

  pavlab_deseq_results.map { it ->
    def parts = it.toString().split("/")
    def cell_type = parts[-3]
    def contrast = parts[-2]
    [contrast, it]
    }.groupTuple(by: [0])
    .set { all_contrasts_pavlab_ct }



  author_label_deseq_results.map { it ->
    def parts = it.toString().split("/")
    def cell_type = parts[-3]
    def pavlab_cell_type = params.gemma_to_gemma_map[cell_type] ?: cell_type
    def contrast = parts[-2]
    [contrast, it]
    }
    .groupTuple(by: [0])
    .set { all_contrasts_author_ct }

  // combine the two channels on contrast and cell type
  all_contrasts_pavlab_ct
    .combine(all_contrasts_author_ct, by: 0)
    .set { pairwise_channel }

  //view
  aggregate_pairwise(pairwise_channel)
}