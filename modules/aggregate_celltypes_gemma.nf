process AggregateCelltypesGemma {
  conda "/home/rschwartz/anaconda3/envs/scanpyenv"
  publishDir "${params.outdir}/ct_pseudobulks/gemma", mode: "copy"
  input:
  path pseudobulk_matrices
  output:
  path "**pseudobulk_matrix.tsv.gz", emit: aggregated_celltypes
  script:
  """
  python $projectDir/bin/aggregate_celltypes_gemma.py \
        --pseudobulk_matrices ${pseudobulk_matrices} \
        --metadata_files ${params.gemma_meta_dir}
  """
}
