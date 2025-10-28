process DEOverlap {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/de_overlap/contrast_overlaps/${contrast}", mode: "copy"
  input:
  tuple val(contrast), path(pavlab_files), path(author_files)
  output:
  path("**pairwise_overlap.tsv")
  script:
  """
  python $projectDir/bin/DE_overlap.py \
      --contrast ${contrast} \
      --sc_pipeline_paths ${pavlab_files} \
      --author_paths ${author_files}
  """
}
