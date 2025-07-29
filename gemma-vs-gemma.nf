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


process DE_corr {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/DE_corr/${contrast}/figs/${pavlab_ct}/${author_ct}", mode: 'copy', pattern: '**png'
  publishDir "${params.outdir}/DE_corr/${contrast}/files/${pavlab_ct}/${author_ct}", mode: 'copy', pattern: '**tsv'


  input:
  tuple val(contrast), val(pavlab_ct), path(pavlab_results), val(author_ct), path(author_results)

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
  """
}

workflow {

  pavlab_deseq_results = Channel
    .fromPath(params.pavlab_deseq_results)

  author_label_deseq_results = Channel
    .fromPath(params.author_label_deseq_results) 

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


    pavlab_deseq_results.map { it ->
    def parts = it.toString().split("/")
    def cell_type = parts[-3]
    def contrast = parts[-2]
    [contrast, cell_type, it]
    }
    .set { ct_contrasts_pavlab }

      author_label_deseq_results.map { it ->
    def parts = it.toString().split("/")
    def cell_type = parts[-3]
    def pavlab_cell_type = params.gemma_to_gemma_map[cell_type] ?: cell_type
    def contrast = parts[-2]
    [contrast, pavlab_cell_type, cell_type, it]
    }
    .set { ct_contrasts_author }

    // combine the cell type specific contrasts from pavlab and author
    ct_contrasts_pavlab.combine(ct_contrasts_author, by: [0,1])
    .set { ct_specific_contrasts }

    ct_specific_contrasts.first().view()

    DE_corr(ct_specific_contrasts)

}