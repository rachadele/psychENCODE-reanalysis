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
  parser.add_argument("--pavlab_paths", type=str, nargs="+", default=["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/DESeq2/manual/astrocyte/Disorder_Schizophrenia_vs_Control/results.tsv","/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/DESeq2/manual/L2_3-6intratelencephalicprojectingglutamatergicneuron/Disorder_Schizophrenia_vs_Control/results.tsv"])
  parser.add_argument("--author_paths", type=str, nargs="+", default=["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/author_contrasts/Schizophrenia/files/Schizophrenia_Astro_degs.tsv","/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/author_contrasts/Schizophrenia/files/Schizophrenia_L2.3.IT_degs.tsv"])               
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
  
def main():
  args= parse_arguments()
  pavlab_paths = args.pavlab_paths
  author_paths = args.author_paths
  contrast = args.contrast
  
  # Load the data
  pavlab_de_results = {}
  for path in pavlab_paths:
    de_results = pd.read_csv(path, sep="\t")
    cell_type = de_results["cell_type"].iloc[0]
    pavlab_de_results[cell_type] = de_results
    
  author_de_results = {}
  for path in author_paths:
    de_results = pd.read_csv(path, sep="\t")
    cell_type = de_results["cell_type"].iloc[0]
    author_de_results[cell_type] = de_results
    
  # all pairwise combinations of cell types
  cell_types_pavlab = list(pavlab_de_results.keys())
  cell_types_author = list(author_de_results.keys())
  combinations = [(ct_pavlab, ct_author) for ct_pavlab in cell_types_pavlab for ct_author in cell_types_author]
  
  
  corr_list = []
  
  for combination in combinations:
    ct_pavlab, ct_author = combination
    pavlab_df = pavlab_de_results[ct_pavlab]
    author_df = author_de_results[ct_author]
    
    # Merge dataframes on gene column
    merged_df = pd.merge(
        author_df, pavlab_df, on="gene", suffixes=("_author", "_gemma"), how="inner"
    )
    
    
    # Calculate Spearman correlation
    merged_df["log2FoldChange_author"] = merged_df["log2FoldChange_author"].astype(float)
    merged_df["log2FoldChange_gemma"] = merged_df["log2FoldChange_gemma"].astype(float)
    
    correlation = merged_df["log2FoldChange_author"].corr(
        merged_df["log2FoldChange_gemma"], method="spearman") 

    # add results to DataFrame
    new_row = {
        "ct_pavlab": ct_pavlab,
        "ct_author": ct_author,
        "spearman_correlation": correlation
    }
    corr_list.append(new_row)



  # Create DataFrame from correlation results
  all_corrs = pd.DataFrame(corr_list)
  
  # Save results to file
  all_corrs.to_csv(f"pairwise_corr_{contrast}.tsv", sep="\t", index=False)

  # Pivot to wide format
  all_corrs_pivot = all_corrs.pivot(
      index="ct_pavlab", columns="ct_author", values="spearman_correlation"
  )

  # --- Optional: Cluster rows and columns ---
  
  # Create the desired order from your list
  desired_col_order = [
      "Lamp5", "Lamp5.Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst.Chodl", "Chandelier", "Vip",
      "L2.3.IT", "L4.IT", "L5.ET", "L5.IT", "L6.IT", "L6.IT.Car3", "L6b", "L5.6.NP", "L6.CT",
      "Oligo", "OPC", "VLMC", "Astro", "Endo", "Immune", "Micro", "PC", "SMC"
  ]

  desired_row_order = [
      "lamp5GABAergiccorticalinterneuron", "pvalbGABAergiccorticalinterneuron", "sncgGABAergiccorticalinterneuron",
      "sstGABAergiccorticalinterneuron", "chandelierpvalbGABAergiccorticalinterneuron",
      "vipGABAergiccorticalinterneuron", "L2_3-6intratelencephalicprojectingglutamatergicneuron",
      "L6bglutamatergiccorticalneuron", "near-projectingglutamatergiccorticalneuron",
      "corticothalamic-projectingglutamatergiccorticalneuron", "oligodendrocyte",
      "oligodendrocyteprecursorcell", "vascularleptomeningealcell", "astrocyte",
      "cerebralcortexendothelialcell", "microglialcell", "pericyte"
  ]

  ordered = all_corrs_pivot.reindex(index=desired_row_order, columns=desired_col_order)
  # drop NA columns
  ordered = ordered.dropna(axis=1, how='all')
  # drop NA rows
  ordered = ordered.dropna(axis=0, how='all')
  # --- Plot heatmap ---
  plt.figure(figsize=(30, 15))
  sns.heatmap(ordered,
              cmap="coolwarm",
              center=0,
              vmin=-1, vmax=1,
              annot=False,
              linewidths=0)

  plt.title(f"Spearman Correlation of DE Results: {contrast}", fontsize=20)
  plt.xlabel("Author Cell Types", fontsize=20)
  plt.ylabel("Pavlab Cell Types", fontsize=20)
  # make font size larger
  plt.xticks(rotation=90, fontsize=20)
  plt.yticks(fontsize=20)
  plt.tight_layout()
  
  # also sav

  plt.savefig(f"pairwise_corr_{contrast}.png")
  plt.close()
  
  
  
  
if __name__ == "__main__":
  main()