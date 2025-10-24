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


def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--contrast", type=str, default="Schizophrenia", help="Contrast to use for DE analysis")
  parser.add_argument("--sc_pipeline_paths", type=str, nargs="+", default=[""])
  parser.add_argument("--author_paths", type=str, nargs="+", default=[""])
  parser.add_argument("--fdr_cutoff", type=float, default=0.05, help="FDR cutoff for significance")
  parser.add_argument("--logFC_threshold", type=float, default=0.0, help="Log fold change threshold for significance")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
  
def main():
  args= parse_arguments()
  sc_pipeline_paths = args.sc_pipeline_paths
  author_paths = args.author_paths
  contrast = args.contrast
  fdr_cutoff = args.fdr_cutoff
  logFC_threshold = args.logFC_threshold
  
  # Load the data
  sc_pipeline_de_results = {}
  for path in sc_pipeline_paths:
    de_results = pd.read_csv(path, sep="\t")
    cell_type = de_results["cell_type"].iloc[0]
    sc_pipeline_de_results[cell_type] = de_results
    
  author_de_results = {}
  for path in author_paths:
    de_results = pd.read_csv(path, sep="\t")
    cell_type = de_results["cell_type"].iloc[0]
    author_de_results[cell_type] = de_results
    
  # all pairwise combinations of cell types
  cell_types_sc_pipeline = list(sc_pipeline_de_results.keys())
  cell_types_author = list(author_de_results.keys())
  
  # only compare same cell types
  
  combinations = [(ct_sc_pipeline, ct_author) for ct_sc_pipeline in cell_types_sc_pipeline for ct_author in cell_types_author]
  
  
  overlap_list = []
  
  for combination in combinations:
    ct_sc_pipeline, ct_author = combination
    sc_pipeline_df = sc_pipeline_de_results[ct_sc_pipeline]
    author_df = author_de_results[ct_author]

    # Merge dataframes on gene column
    merged_df = pd.merge(
        author_df, sc_pipeline_df, on="gene", suffixes=("_author", "_gemma"), how="inner"
    )
    print(merged_df.columns)

    # Identify significant genes in each set
    sig_author = set(merged_df.loc[(merged_df["padj_author"] < fdr_cutoff) & (abs(merged_df["log2FoldChange_author"]) > logFC_threshold), "gene"])
    sig_gemma = set(merged_df.loc[(merged_df["padj_gemma"] < fdr_cutoff) & (abs(merged_df["log2FoldChange_gemma"]) > logFC_threshold), "gene"])
    if len(sig_author | sig_gemma) == 0:
        percent_overlap = np.nan
    else:
        percent_overlap = 100 * len(sig_author & sig_gemma) / len(sig_author | sig_gemma)

    overlap_list.append({
        "ct_sc_pipeline": ct_sc_pipeline,
        "ct_author": ct_author,
        "percent_overlap": percent_overlap
    })

  # Create DataFrame from overlap results
  all_overlaps = pd.DataFrame(overlap_list)
  all_overlaps.to_csv(f"{contrast}_pairwise_overlap.tsv", sep="\t", index=False)


if __name__ == "__main__":
  main()