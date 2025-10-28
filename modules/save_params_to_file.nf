process save_params_to_file {
    publishDir (
        "${params.outdir}",
        mode: "copy"
    )
    output:
    file "params.yaml"
    script:
    """
    cat <<EOF > params.yaml
    from_gemma : ${params.from_gemma}
    gemma_meta_dir : ${params.gemma_meta_dir}
    h5ad_files : ${params.h5ad_files}
    outdir: ${params.outdir}
    EOF
    """
}
