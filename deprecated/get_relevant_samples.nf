process GetRelevantSamples {
  publishDir "${params.outdir}/relevant_samples", mode: "copy"
  input:
  tuple val(pavlab_ct), val(author_ct), path(pseudobulk_matrix)
  output:
  tuple val(author_ct), path("**relevant_samples.tsv"), emit: relevant_samples
  script:
  """
  zcat ${pseudobulk_matrix} | head -1 |
  cut -f 7- | tr '\t' '\n' > ${author_ct}_relevant_samples.tsv
  """
}
