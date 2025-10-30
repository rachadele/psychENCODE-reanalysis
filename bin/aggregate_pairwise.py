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
  combinations = [(ct_sc_pipeline, ct_author) for ct_sc_pipeline in cell_types_pavlab for ct_author in cell_types_author]
  
  
  corr_list = []
  for combination in combinations:
    ct_sc_pipeline, ct_author = combination
    pavlab_df = pavlab_de_results[ct_sc_pipeline]
    author_df = author_de_results[ct_author]
    # Merge dataframes on gene column
    merged_df = pd.merge(
        author_df, pavlab_df, on="gene", suffixes=("_author", "_gemma"), how="inner"
    )
    # Calculate Spearman correlation for log2FoldChange
    merged_df["log2FoldChange_author"] = merged_df["log2FoldChange_author"].astype(float)
    merged_df["log2FoldChange_gemma"] = merged_df["log2FoldChange_gemma"].astype(float)
    corr_logfc = merged_df["log2FoldChange_author"].corr(
        merged_df["log2FoldChange_gemma"], method="spearman")
    # Calculate Spearman correlation for pvalue
    merged_df["pvalue_author"] = merged_df["pvalue_author"].astype(float)
    merged_df["pvalue_gemma"] = merged_df["pvalue_gemma"].astype(float)
    corr_pvalue = merged_df["pvalue_author"].corr(
        merged_df["pvalue_gemma"], method="spearman")
    # add results to DataFrame
    new_row = {
        "ct_sc_pipeline": ct_sc_pipeline,
        "ct_author": ct_author,
        "spearman_log2FoldChange": corr_logfc,
        "spearman_pvalue": corr_pvalue
    }
    corr_list.append(new_row)

  # Create DataFrame from correlation results
  all_corrs = pd.DataFrame(corr_list)
  all_corrs.to_csv(f"{contrast}_pairwise_corr.tsv", sep="\t", index=False)

  # Pivot to wide format for log2FoldChange
  all_corrs_logfc = all_corrs.pivot(
      index="ct_sc_pipeline", columns="ct_author", values="spearman_log2FoldChange"
  )
  # Pivot to wide format for pvalue
  all_corrs_pvalue = all_corrs.pivot(
      index="ct_sc_pipeline", columns="ct_author", values="spearman_pvalue"
  )

  # Create the desired order from your list
  desired_col_order = [
      "Lamp5", "Lamp5_Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Chandelier", "Vip",
      "L2_3_IT", "L4_IT",  "L5_IT", "L6_IT", "L6_IT_Car3", "L5_ET", "L6b", "L5_6_NP", "L6_CT",
      "Oligo", "OPC", "Astro", "Micro", "Immune", "VLMC", "Endo", "PC", "SMC",
  ]
  print(all_corrs_logfc.index)

  # Use log2FoldChange pivot for ordering
  if "L2.3_IT" in all_corrs_logfc.index:
    desired_row_order = [
        "Lamp5", "Lamp5_Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Chandelier", "Vip",
        "L2.3_IT", "L4_IT", "L5_ET", "L5_IT", "L6_IT", "L6_IT_Car3", "L6b", "L5.6_NP", "L6_CT",
        "Oligo", "OPC", "Astro", "Micro", "Immune", "Endo", "PC", "SMC", "VLMC"
    ]
  elif "vascular" in all_corrs_logfc.index:
    desired_row_order = [
      "lamp5.GABAergic.cortical.interneuron", "pvalb.GABAergic.cortical.interneuron",
      "sncg.GABAergic.cortical.interneuron", "sst.GABAergic.cortical.interneuron",
      "chandelier.pvalb.GABAergic.cortical.interneuron",
      "vip.GABAergic.cortical.interneuron", "L2.3.6.intratelencephalic.projecting.glutamatergicneuron",
      "deep.layer.non.IT", "oligodendrocyte",
      "oligodendrocyte.precursor.cell", "astrocyte",
      "microglial.cell", "vascular"
    ]
  elif "corticothalamic.projecting.glutamatergic.cortical.neuron" in all_corrs_logfc.index:
    desired_row_order = [
      "lamp5.GABAergic.cortical.interneuron", "pvalb.GABAergic.cortical.interneuron",
      "sncg.GABAergic.cortical.interneuron", "sst.GABAergic.cortical.interneuron",
      "chandelier.pvalb.GABAergic.cortical.interneuron",
      "vip.GABAergic.cortical.interneuron", 
      "L2.3.6.intratelencephalic.projecting.glutamatergicneuron",
      "L5.extratelencephalic.projecting.glutamatergic.cortical.neuron", 
      "L6b.glutamatergic.cortical.neuron", 
      "near.projecting.glutamatergic.cortical.neuron", 
      "corticothalamic.projecting.glutamatergic.cortical.neuron", 
      "oligodendrocyte", "oligodendrocyte.precursor.cell", "astrocyte", "microglial.cell", 
      "endothelial.cell", "pericyte", "smooth.muscle.cell"
    ]
  
  
  # Plot heatmap for log2FoldChange
  ordered_logfc = all_corrs_logfc.reindex(index=desired_row_order, columns=desired_col_order)
  ordered_logfc = ordered_logfc.dropna(axis=1, how='all').dropna(axis=0, how='all')
  plt.figure(figsize=(30, 15))
  sns.heatmap(ordered_logfc, cmap="Reds", center=0.5, vmin=0, vmax=1, annot=False, linewidths=0)
  plt.title(f"Spearman Correlation of log2FoldChange: {contrast}", fontsize=20)
  plt.xlabel("Author Cell Types", fontsize=20)
  plt.ylabel("Pavlab Cell Types", fontsize=20)
  plt.xticks(rotation=90, fontsize=20)
  plt.yticks(fontsize=20)
  plt.tight_layout()
  plt.savefig(f"{contrast}_pairwise_corr_log2FoldChange.png")
  plt.close()

  # Plot heatmap for pvalue
  ordered_pvalue = all_corrs_pvalue.reindex(index=desired_row_order, columns=desired_col_order)
  ordered_pvalue = ordered_pvalue.dropna(axis=1, how='all').dropna(axis=0, how='all')
  plt.figure(figsize=(30, 15))
  sns.heatmap(ordered_pvalue, cmap="Blues", center=0.5, vmin=0, vmax=1, annot=False, linewidths=0)
  plt.title(f"Spearman Correlation of pvalue: {contrast}", fontsize=20)
  plt.xlabel("Author Cell Types", fontsize=20)
  plt.ylabel("Pavlab Cell Types", fontsize=20)
  plt.xticks(rotation=90, fontsize=20)
  plt.yticks(fontsize=20)
  plt.tight_layout()
  plt.savefig(f"{contrast}_pairwise_corr_pvalue.png")
  plt.close()
  
  
  
  
if __name__ == "__main__":
  main()