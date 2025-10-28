process DESeq2AnalysisManual {
  conda "/home/rschwartz/anaconda3/envs/r4.3"
  publishDir "${params.outdir}/DESeq2/manual/${cell_type}", mode: "copy"
  input:
  tuple val(cell_type), path(pseudobulk_matrix), path(pseudobulk_metadata)
  output:
  tuple val(cell_type), path("**results.tsv"), emit: all_contrasts_manual
  path "**png"
  script:
  """
  Rscript $projectDir/bin/DESeq2_analysis.R --pseudobulk_matrix ${pseudobulk_matrix} --metadata ${pseudobulk_metadata} \
        --cell_type ${cell_type} \
        --mode manual
  """
}
