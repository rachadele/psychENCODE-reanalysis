process GetGemmaPseudobulks {
  publishDir "${params.outdir}/experiment_pseudobulks/gemma/${experiment}", mode: "copy"
  input:
  val experiment
  output:
  path("**.tsv.gz"), emit: aggregated_data
  script:
  """
  if [ ${params.author_submitted} = true ]; then
    gemma-cli getSingleCellDataMatrix -e $experiment --aggregate-by-assay \
    --aggregate-by-cell-type-assignment author-submitted \
    --output-file "${experiment}_aggregated_gemma.tsv.gz"
  else
    gemma-cli getSingleCellDataMatrix -e $experiment --aggregate-by-assay \
    --aggregate-by-preferred-cell-type-assignment \
    --output-file "${experiment}_aggregated_gemma.tsv.gz"
  fi
  """
}
