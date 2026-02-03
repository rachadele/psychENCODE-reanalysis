process AverageDECorrelations {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/average_corr/filtered_per_contrast", mode: "copy", pattern: "filtered**.tsv"
  publishDir "${params.outdir}/average_corr/per_celltype", mode: "copy", pattern: "**per_celltype_**.tsv"
  publishDir "${params.outdir}/average_corr/overall_average", mode: "copy", pattern: "**average_corr*.tsv"
  publishDir "${params.outdir}/average_corr/figs", mode: "copy", pattern: "**png"

  input:
  path corr_tables

  output:
  path "filtered**.tsv"
  path "**per_celltype_**.tsv"
  path "**average_corr*.tsv"
  path "**png"
  script:
  """
  python $projectDir/bin/average_corr.py --corr_tables ${corr_tables} \\
      --annotation_level ${params.level}
  """
}
