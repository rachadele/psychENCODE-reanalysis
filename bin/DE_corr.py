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


def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--contrast", type=str, default="Schizophrenia",
					  help="Contrast to use for DE analysis")
  parser.add_argument("--author_results", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/author_contrasts/Schizophrenia/Schizophrenia_L2.3.IT_degs.tsv")
  parser.add_argument("--gemma_results", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/DESeq2/manual/L2_3-6intratelencephalicprojectingglutamatergicneuron/Disorder_Schizophrenia_vs_Control/results.tsv")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def label_significance(row):
  if row["padj_author"] < 0.05 and row["padj_gemma"] < 0.05:
      return "Both"
  elif row["padj_author"] < 0.05:
      return "Author Only"
  elif row["padj_gemma"] < 0.05:
      return "Gemma Only"
  else:
      return "Not Significant"
      
def main():
  args = parse_arguments()
  contrast = args.contrast
  author_results = args.author_results
  gemma_results = args.gemma_results
  # Load
  author_df = pd.read_csv(author_results, sep="\t")
  gemma_df = pd.read_csv(gemma_results, sep="\t")
  
  # extract cell type names for each result
  author_cell_type = author_df["cell_type"].unique()[0]
  gemma_cell_type = gemma_df["cell_type"].unique()[0] 
 
  # sort genes in the same order and make sure they are the same genes
  # merge dataframe on gene column
  
  merged_df = pd.merge(
      author_df, gemma_df, on="gene", suffixes=("_author", "_gemma")
  )
  #record missing genes
  
  only_in_gemma = set(gemma_df["gene"]) - set(author_df["gene"])
  only_in_author = set(author_df["gene"]) - set(gemma_df["gene"])

  # writing missing genes to one file with two columns
  missing_genes_df = pd.DataFrame({
      "gene": list(only_in_author.union(only_in_gemma)),
      "in_author": [gene in only_in_author for gene in only_in_author.union(only_in_gemma)],
      "in_gemma": [gene in only_in_gemma for gene in only_in_author.union(only_in_gemma)]
  })
  
  missing_genes_df.to_csv("missing_genes.tsv", sep="\t", index=False)
  
  
 # check if df is empty
  if merged_df.empty:
    # message and exit
    raise ValueError(
        f"No overlapping genes found between author and Gemma results for {contrast}: {author_cell_type} vs {gemma_cell_type}"
    )
  merged_df["Significance"] = merged_df.apply(label_significance, axis=1)
  
  correlation = merged_df["log2FoldChange_author"].corr(
      merged_df["log2FoldChange_gemma"], method="spearman") 

  hue_order = ["Not Significant","Author Only", "Gemma Only","Both"]
  # make hues static colors
  palette = {
    "Not Significant": "lightgray",
    "Author Only": "royalblue",
    "Gemma Only": "darkorange",
    "Both": "crimson"
}
  plt.figure(figsize=(10, 6))

  # Scatter plot with color based on significance
  sns.scatterplot(
      data=merged_df,
      x="log2FoldChange_author",
      y="log2FoldChange_gemma",
      hue="Significance",
      hue_order=hue_order,
      palette=palette,
      alpha=0.6,
      edgecolor=None
  )

  # Add regression line (on all points)
  sns.regplot(
      data=merged_df,
      x="log2FoldChange_author",
      y="log2FoldChange_gemma",
      scatter=False,  # Don't add scatter again
      line_kws={"color": "black", "linewidth": 2}
  )

  plt.title(f"Contrast: {contrast} - {author_cell_type} vs {gemma_cell_type}")
  plt.text(0.05, 0.95, f"Spearman Correlation: {correlation:.2f}",
           transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
  
  # add missing genes length as text
  plt.text(0.05, 0.90, f"Only in Author: {len(only_in_author)}",
           transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
  plt.text(0.05, 0.85, f"Only in Gemma: {len(only_in_gemma)}",
           transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
  
  plt.xlabel(f"Author: {author_cell_type} Log2 Fold Change")
  plt.ylabel(f"Gemma: {gemma_cell_type} Log2 Fold Change")
  plt.legend(title="Significance")
  plt.grid(True)
  plt.tight_layout()
  plt.savefig("DE_correlation.png")
  plt.close()
  
# need to save the merged df using an "outer" merge
  outer_merge = pd.merge(
      author_df, gemma_df, on="gene", suffixes=("_author", "_gemma"), how="outer"
  )
  outer_merge.to_csv(f"{contrast}_{author_cell_type}_{gemma_cell_type}.tsv", sep="\t", index=False)
  

if __name__ == "__main__":
  main()