process AggregatePairwise {
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
  python $projectDir/bin/aggregate_pairwise.py \
      --contrast ${contrast} \
      --pavlab_paths ${pavlab_files} \
      --author_paths ${author_files}
  """
}
