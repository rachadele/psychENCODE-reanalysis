process AverageDEOverlaps {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/de_overlap/average_overlaps/filtered_per_contrast", mode: "copy", pattern: "filtered**.tsv"
  publishDir "${params.outdir}/de_overlap/average_overlaps/per_celltype", mode: "copy", pattern: "per_celltype_averages.tsv"
  publishDir "${params.outdir}/de_overlap/average_overlaps/overall_average", mode: "copy", pattern: "average_de_overlaps.tsv"
  publishDir "${params.outdir}/de_overlap/average_overlaps/figs", mode: "copy", pattern: "**png"
  input:
  path de_overlap_files
  output:
  path "filtered**.tsv"
  path "per_celltype_averages.tsv"
  path "average_de_overlaps.tsv"
  path "**png"
  script:
  """
  python $projectDir/bin/average_de_overlaps.py \
      --de_overlap_paths ${de_overlap_files}
  """
}
