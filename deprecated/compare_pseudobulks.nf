process ComparePseudobulks {
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
  python $projectDir/bin/compare_pseudobulks.py \
        --author_pseudobulk ${author_pseudobulk} \
        --author_cell_type ${author_cell_type} \
        --pavlab_pseudobulk ${pavlab_pseudobulk} \
        --pavlab_cell_type ${pavlab_cell_type} \
        --gemma_metadata ${params.gemma_meta_dir} \
        --mode ${mode}
  """
}
