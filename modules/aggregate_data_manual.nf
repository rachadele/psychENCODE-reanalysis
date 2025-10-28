process AggregateDataManual {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/experiment_pseudobulks/manual/${experiment}", mode: "copy"
  input:
  tuple val(experiment), path(h5ad_file)
  output:
  path "**pseudobulk.h5ad", emit: aggregated_experiments
  script:
  """
  python $projectDir/bin/aggregate_data_manual.py \
        --h5ad_file ${h5ad_file} \
        --cell_type_column ${params.cell_type_column} \
        ${params.filter_samples ? '--filter_samples' : ''}
  """
}
