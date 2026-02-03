process WrangleAuthor {
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
  python $projectDir/bin/wrangle_author_de_results.py \
        --author_degs ${author_results} \
        --contrast ${contrast}
  """
}
