process DESeq2AnalysisGemma {
  conda "/home/rschwartz/anaconda3/envs/r4.3"
  publishDir "${params.outdir}/DESeq2/gemma/${pavlab_cell_type}", mode: "copy"
  input:
  tuple val(pavlab_cell_type), path(pseudobulk_matrix)
  output:
  tuple val(pavlab_cell_type), path("**results.tsv"), emit: all_contrasts_gemma
  path "**png"
  script:
  """
  Rscript $projectDir/bin/DESeq2_analysis.R \\
    --pseudobulk_matrix ${pseudobulk_matrix} \\
    --metadata ${params.gemma_meta_dir} \\
    --mode gemma \\
    --cell_type ${pavlab_cell_type} \\
    ${params.filter_genes ? "--filter_genes" : ""}
  """
}
