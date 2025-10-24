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
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
  
def main():
  args= parse_arguments()
  sc_pipeline_paths = args.sc_pipeline_paths
  author_paths = args.author_paths
  contrast = args.contrast
  
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
  combinations = [(ct_sc_pipeline, ct_author) for ct_sc_pipeline in cell_types_sc_pipeline for ct_author in cell_types_author]
  
  
  corr_list = []
  
  for combination in combinations:
    ct_sc_pipeline, ct_author = combination
    sc_pipeline_df = sc_pipeline_de_results[ct_sc_pipeline]
    author_df = author_de_results[ct_author]
    
    # Merge dataframes on gene column
    merged_df = pd.merge(
        author_df, sc_pipeline_df, on="gene", suffixes=("_author", "_gemma"), how="inner"
    )
    
    
    # Calculate Spearman correlation
    merged_df["log2FoldChange_author"] = merged_df["log2FoldChange_author"].astype(float)
    merged_df["log2FoldChange_gemma"] = merged_df["log2FoldChange_gemma"].astype(float)
    
    correlation = merged_df["log2FoldChange_author"].corr(
        merged_df["log2FoldChange_gemma"], method="spearman") 

    # add results to DataFrame
    new_row = {
        "ct_sc_pipeline": ct_sc_pipeline,
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
      index="ct_sc_pipeline", columns="ct_author", values="spearman_correlation"
  )

  # --- Optional: Cluster rows and columns ---
  print(all_corrs_pivot.index)
  print(all_corrs_pivot.columns)
  # Create the desired order from your list
  desired_col_order = [
        "Lamp5", "Lamp5_Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Chandelier", "Vip",
        "L2_3_IT", "L4_IT", "L5_IT", "L6_IT", "L6_IT_Car3", "L5_ET", "L6b", "L5_6_NP", "L6_CT",
        "Oligo", "OPC", "VLMC", "Astro", "Endo", "Immune", "Micro", "PC", "SMC"
    ]
     
  if "L2.3_IT" in all_corrs_pivot.index:
    desired_row_order = [
        "Lamp5", "Lamp5_Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Chandelier", "Vip",
        "L2.3_IT", "L4_IT", "L5_ET", "L5_IT", "L6_IT", "L6_IT_Car3", "L6b", "L5.6_NP", "L6_CT",
        "Oligo", "OPC", "VLMC", "Astro", "Endo", "Immune", "Micro", "PC", "SMC"
    ]
    
  elif "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron" in all_corrs_pivot.index:
  # change to use periods
  #
    desired_row_order = [
      "lamp5.GABAergic.cortical.interneuron", "pvalb.GABAergic.cortical.interneuron",
      "sncg.GABAergic.cortical.interneuron", "sst.GABAergic.cortical.interneuron",
      "chandelier.pvalb.GABAergic.cortical.interneuron",
      "vip.GABAergic.cortical.interneuron", "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
      "L5.extratelencephalic.projecting.glutamatergic.cortical.neuron", 
      "L6b.glutamatergic.cortical.neuron", "near.projecting.glutamatergic.cortical.neuron",
      "L6.corticothalamic.projecting.glutamatergic.cortical.neuron", "oligodendrocyte",
      "oligodendrocyte.precursor.cell", "vascular.leptomeningeal.cell", "astrocyte",
      "endothelial.cell", "microglial.cell", "pericyte"
    ]
  

  ordered = all_corrs_pivot.reindex(index=desired_row_order, columns=desired_col_order)
  # drop NA columns
  ordered = ordered.dropna(axis=1, how='all')
  # drop NA rows
  ordered = ordered.dropna(axis=0, how='all')
  # --- Plot heatmap ---
  plt.figure(figsize=(30, 15))
  sns.heatmap(ordered,
              # just red
              cmap="Reds",
              center=0.5,
              vmin=0, vmax=1,
              annot=False,
              linewidths=0)

  plt.title(f"Spearman Correlation of DE Results: {contrast}", fontsize=20)
  plt.xlabel("Author Cell Types", fontsize=20)
  plt.ylabel("sc_pipeline Cell Types", fontsize=20)
  # make font size larger
  plt.xticks(rotation=90, fontsize=20)
  plt.yticks(fontsize=20)
  plt.tight_layout()
  
  # also sav

  plt.savefig(f"pairwise_corr_{contrast}.png")
  plt.close()
  
  
  
  
if __name__ == "__main__":
  main()