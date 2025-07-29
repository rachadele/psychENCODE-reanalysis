import os
import argparse
import glob
import pandas as pd
from functools import reduce
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import rankdata, spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list
from pathlib import Path

def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--contrast", type=str, default="Schizophrenia", help="Contrast to use for DE analysis")
  parser.add_argument("--author_cell_type", type=str, default="Astrocyte", help="Cell type to use for author DE analysis")
  parser.add_argument("--pavlab_cell_type", type=str, default="astrocyte", help="Cell type to use for pavlab DE analysis")
  parser.add_argument("--pavlab_de_results", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results_author_submitted_false_from_gemma_false/DESeq2/manual/astrocyte/Disorder_Schizophrenia_vs_Control/Disorder_Schizophrenia_vs_Control_astrocyte_results.tsv")
  parser.add_argument("--author_de_results", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results_author_submitted_false_from_gemma_false/author_contrasts/Schizophrenia/files/Schizophrenia_Astro_degs.tsv")
  parser.add_argument("--pavlab_pseudobulk", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results_author_submitted_false_from_gemma_false/ct_pseudobulks/manual/astrocyte/astrocyte_pseudobulk_matrix.tsv.gz")
  parser.add_argument("--author_pseudobulk", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/pseudobulks/pseudobulk_expr/Astro.expr.bed.gz")
  parser.add_argument("--author_metadata", type=str,
      default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/PEC2_sample_metadata.txt",
      help="Author metadata file in text format")
  parser.add_argument("--gemma_metadata", type=str,
      default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata",
      help="Metadata directory for all GEMMA experiments")
  
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
    # make pavlab logCPM
    
    

def plot_top_genes_expression(pavlab_pseudobulk, author_pseudobulk, 
                              top_genes, contrast,
                              pavlab_cell_type, author_cell_type):
    """Plot the expression of top genes from both datasets."""
    # filter h5ad objects for top genes
    # filter for disorder in contrast
    pavlab = pavlab_pseudobulk[pavlab_pseudobulk.obs["Disorder"].isin([contrast,"Control"])]
    author = author_pseudobulk[author_pseudobulk.obs["Disorder"].isin([contrast,"Control"])]
    
    pavlab.obs.index = pavlab.obs["Individual_ID"].astype(str)
    author.obs.index = author.obs["Individual_ID"].astype(str)
    # restrict to the same samples
    common_samples = list(set(pavlab.obs["Individual_ID"]).intersection(set(author.obs["Individual_ID"])))
    pavlab = pavlab[common_samples]
    author = author[common_samples]
    
    # check if genes are in both datasets
    common_genes = list(set(top_genes).intersection(set(pavlab.var_names)).intersection(set(author.var_names)))
    
    sc.pl.heatmap(pavlab, var_names=common_genes, groupby=["Disorder"], show=False,
                  save=f"pavlab_{pavlab_cell_type}_{contrast}_top_genes.png")
    
    sc.pl.heatmap(author, var_names=common_genes, groupby=["Disorder"],
                  show=False, 
                  save=f"author_{author_cell_type}_{contrast}_top_genes.png")
 
    
def log_cpm(counts_df, prior_count=1):
    """Calculate log2 CPM with a small prior."""
    lib_sizes = counts_df.sum(axis=0)
    cpm = (counts_df / lib_sizes) * 1e6
    return np.log2(cpm + prior_count)  # add prior to avoid log(0)
  
def make_anndata(pseudobulk_df, metadata_df):
    """Create an AnnData object from a pseudobulk DataFrame and metadata DataFrame."""
    adata = sc.AnnData(pseudobulk_df.T)
    adata.obs = adata.obs_names.to_frame(name="Individual_ID")
    adata.obs = adata.obs.merge(metadata_df, on="Individual_ID", how="left")
    adata.obs["Individual_ID"] = adata.obs["Individual_ID"].astype(str)
    return adata
  
def main():
    
    args = parse_arguments()
    contrast = args.contrast
    author_cell_type = args.author_cell_type
    pavlab_cell_type = args.pavlab_cell_type
    pavlab_de_results_path = args.pavlab_de_results
    author_de_results_path = args.author_de_results
    pavlab_pseudobulk_path = args.pavlab_pseudobulk
    author_pseudobulk_path = args.author_pseudobulk
    author_metadata_path = args.author_metadata
    gemma_metadata_path = args.gemma_metadata
    
    # Load the data
    pavlab_de_results = pd.read_csv(pavlab_de_results_path, sep="\t")
    author_de_results = pd.read_csv(author_de_results_path, sep="\t")
    pavlab_pseudobulk = pd.read_csv(pavlab_pseudobulk_path, sep="\t", index_col=0)
    author_pseudobulk = (pd.read_csv(args.author_pseudobulk, sep="\t")
                         .drop(columns=["#chr", "start", "end", "length", "strand"], errors="ignore")
                         .set_index("gene"))
    
    
    # Load and combine GEMMA metadata
    meta = []
    for fn in Path(args.gemma_metadata).glob("*.tsv"):
        df = pd.read_csv(fn, sep="\t", dtype=str)
        df["Cohort"] = fn.name.replace("_sample_meta.tsv", "")
        meta.append(df)
    gemma_meta = pd.concat(meta, ignore_index=True)
    author_meta = pd.read_csv(args.author_metadata, sep="\t", dtype=str)
    
    pavlab_pseudobulk = log_cpm(pavlab_pseudobulk) 
    
    # make h5ad object for author pseudobulk
    author_pseudobulk = make_anndata(author_pseudobulk, author_meta)
    # make h5ad object for pavlab pseudobulk
    pavlab_pseudobulk = make_anndata(pavlab_pseudobulk, gemma_meta)
    
    

    # plot the difference in expression for the top 2 genes in each de list
    pavlab_top_genes = pavlab_de_results.nlargest(10, "log2FoldChange")["gene"].tolist()
    author_top_genes = author_de_results.nlargest(10, "log2FoldChange")["gene"].tolist()
  
    top_genes = list(set(pavlab_top_genes).union(set(author_top_genes)))
  
    plot_top_genes_expression(pavlab_pseudobulk, author_pseudobulk, top_genes, contrast, author_cell_type, pavlab_cell_type)    
    
    
    
if __name__ == "__main__":
	main()