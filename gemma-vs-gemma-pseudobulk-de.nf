#!/usr/bin/env nextflow

// Import all processes from modules
include { AggregatePairwise } from "${projectDir}/modules/aggregate_pairwise.nf"
include { DECorr } from "${projectDir}/modules/DE_corr.nf"
include { DEOverlap } from "${projectDir}/modules/DE_overlap.nf"
include { AverageDEOverlaps } from "${projectDir}/modules/average_de_overlaps.nf"


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
    AggregatePairwise(pairwise_channel)

    // DE overlap
    DEOverlap(pairwise_channel)

    DE_overlap_results = DEOverlap.out.collect()
    AverageDEOverlaps(DE_overlap_results)

    // get cell type, contrast from each path
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
    DECorr(ct_specific_contrasts)

}