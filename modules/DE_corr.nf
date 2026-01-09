process DECorr {
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
  python $projectDir/bin/DE_corr.py \
        --pavlab_results ${pavlab_results} \
        --author_results ${author_results} \
        --contrast ${contrast}
  """
}
