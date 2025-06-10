import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Set base theme
sns.set(style="whitegrid", font_scale=1.5)
plt.rcParams["axes.grid"] = False


def parse_args():
  """Parse command line arguments."""
  parser = argparse.ArgumentParser(description="Wrangle author pseudobulks")
  parser.add_argument("--author_pseudobulk", type=str,
      default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/pseudobulks/pseudobulk_expr/Endo__VLMC.expr.bed.gz",
      help="Author pseudobulk file in bed format")
  parser.add_argument("--pavlab_pseudobulk", type=str,
      default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/ct_pseudobulks/gemma/endothelialcell/endothelialcell_pseudobulk_matrix.tsv.gz",
      help="PavLab pseudobulk file in tsv format (either GEMMA or manually created)")
  parser.add_argument("--gemma_metadata", type=str,
      default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata",
      help="Metadata directory for all GEMMA experiments")
  parser.add_argument("--author_metadata", type=str,
      default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/PEC2_sample_metadata.txt",
      help="Author metadata file in text format")
  parser.add_argument("--author_cell_type", type=str,
      default="L5.6.NP", help="Cell type name in author pseudobulk file")
  parser.add_argument("--pavlab_cell_type", type=str,
      default="Endo__VLMC", help="Cell type name in PavLab pseudobulk file")
  parser.add_argument("--mode", type=str, default="manual", help="Mode of comparison: 'manual' or 'gemma'")

  return parser.parse_args()

def log_cpm(counts_df, prior_count=1):
    """Calculate log2 CPM with a small prior."""
    lib_sizes = counts_df.sum(axis=0)
    cpm = (counts_df) / lib_sizes * 1e6
    return np.log2(cpm + prior_count)  # add prior to avoid log(0)
  
def plot_corr(merged, spearman_corr):
    # Plot with correlation line
  plt.figure(figsize=(8, 6))
  sns.scatterplot(
      data=merged,
      x="Author_Expression", y="PavLab_Expression",
      alpha=0.6, edgecolor=None, hue="Cohort"
  )
  # Optional: still draw a linear trend line (for visualization only)
  sns.regplot(
      data=merged,
      x="Author_Expression", y="PavLab_Expression",
      scatter=False, color="red"
  )
  plt.text(0.05, 0.95, f"Spearman r = {spearman_corr:.2f}", transform=plt.gca().transAxes, fontsize=12, color="red")
  plt.title("Gene Expression: Author vs PavLab")
  plt.xlabel("Author Expression")
  plt.ylabel("PavLab Expression")
  plt.legend(title="Cohort", bbox_to_anchor=(1.05, 1), loc='upper left')
  plt.tight_layout()
  plt.savefig("gene_by_sample_corr_author_vs_pavlab.png")
  plt.close()

def plot_altmann(comp, args):
    # Bland-Altman plot
  plt.figure(figsize=(10, 10))
  sns.scatterplot(data=comp, x="avg", y="diff", hue="Cohort", alpha=0.7)
  plt.axhline(0, linestyle="--", color="gray")
  plt.title(f"Bland‑Altman for {args.author_cell_type} vs {args.pavlab_cell_type}")
  plt.xlabel("Average log CPM")
  plt.ylabel("Author − PavLab (log CPM)")
  plt.legend(title="Cohort", bbox_to_anchor=(1.05, 1), loc='upper left')
  plt.tight_layout()
  plt.savefig(f"pseudobulk_comparison_plot_{args.author_cell_type}_{args.pavlab_cell_type}.png")
  plt.close()

def main():
  # Load author pseudobulk
  args = parse_args()
  author = (
      pd.read_csv(args.author_pseudobulk, sep="\t")
        .drop(columns=["#chr", "start", "end", "length", "strand"], errors="ignore")
        .set_index("gene")
  )

  # Load PavLab pseudobulk
  pavlab = (
      pd.read_csv(args.pavlab_pseudobulk, sep="\t")
        .set_index("feature_name")
  )

  # Load and combine GEMMA metadata
  meta = []
  for fn in Path(args.gemma_metadata).glob("*.tsv"):
      df = pd.read_csv(fn, sep="\t", dtype=str)
      df["Cohort"] = fn.name.replace("_sample_meta.tsv", "")
      meta.append(df)
  gemma = pd.concat(meta, ignore_index=True)

  # Map PavLab sample IDs to Individual_ID
  if args.mode == "manual":
    mapping = gemma.set_index("sample_id")["Individual_ID"].to_dict()
    pavlab.columns = [mapping.get(x, np.nan) for x in pavlab.columns]
  
  # Load author metadata
  author_meta = pd.read_csv(args.author_metadata, sep="\t", dtype=str)
  
  # Build matching table
  author_ids = set(author.columns)
  pavlab_ids = set(pavlab.columns)
  # get rid of floats 

  # if either are empty, terminate without error
  if not author_ids or not pavlab_ids:
      print("No matching samples found in either author or PavLab pseudobulks.")
      
    
  all_ids = sorted(author_ids | pavlab_ids)

  match_df = pd.DataFrame({
      "Individual_ID": all_ids,
      "in_author": [i in author_ids for i in all_ids],
      "in_pavlab": [i in pavlab_ids for i in all_ids]
  }).merge(author_meta, on="Individual_ID", how="left")
  
  match_df["Cohort"] = match_df["Cohort"].fillna("Unknown")

  match_df.to_csv(
      f"matching_samples_{args.author_cell_type}_{args.pavlab_cell_type}.tsv",
      sep="\t", index=False
  )
  match_df["Individual_ID"] = match_df["Individual_ID"].astype(str) # get rid of floats
  # Summarize cohort
  summary = (
      match_df
        .assign(
            n_matching=lambda df: df.in_author & df.in_pavlab,
            n_only_author=lambda df: df.in_author & ~df.in_pavlab,
            n_only_pavlab=lambda df: ~df.in_author & df.in_pavlab
        )
        
        .groupby("Cohort")
        .agg(
            n_samples=("Individual_ID", "size"),
            n_matching=("n_matching", "sum"),
            n_only_in_author=("n_only_author", "sum"),
            n_only_in_pavlab=("n_only_pavlab", "sum")
        )
        .reset_index()
  )

  pivot = summary.melt(
      id_vars="Cohort",
      value_vars=["n_matching", "n_only_in_author", "n_only_in_pavlab"],
      var_name="source",
      value_name="count"
  )

  # Plot sample matching
  # write text of counts on plot
  plt.figure(figsize=(12, len(pivot['Cohort'].unique()) * 1))  # dynamically scale height
  barplot = sns.barplot(
      data=pivot, 
      x="count", 
      y="Cohort", 
      hue="source", 
      dodge=True
  )

  # Add value labels on bars
  for container in barplot.containers:
      barplot.bar_label(container, fmt="%d", label_type="edge", padding=3, fontsize=10)

  plt.tight_layout()
  plt.savefig(f"sample_matching_summary_{args.author_cell_type}_{args.pavlab_cell_type}.png")
  plt.close()

  # Compute CPM and sums
  pav_cpm = log_cpm(pavlab, prior_count=1)

  author_sum = author.sum(axis=0)
  pav_sum = pav_cpm.sum(axis=0)

  comp = pd.DataFrame({
      "Individual_ID": all_ids,
      "Author_Counts": [author_sum.get(i, np.nan) for i in all_ids],
      "Pavlab_Counts": [pav_sum.get(i, np.nan) for i in all_ids],
  }).merge(author_meta[["Individual_ID", "Cohort"]], on="Individual_ID", how="left")
  
  # fille NA with 0
  comp["Author_Counts"] = comp["Author_Counts"].fillna(0)
  comp["Pavlab_Counts"] = comp["Pavlab_Counts"].fillna(0)
  comp["Cohort"] = comp["Cohort"].fillna("Unknown")

  comp.to_csv(
      f"pseudobulk_comparison_{args.author_cell_type}_{args.pavlab_cell_type}.tsv",
      sep="\t", index=False
  )

  comp["avg"] = comp[["Author_Counts", "Pavlab_Counts"]].mean(axis=1)
  comp["diff"] = comp["Author_Counts"] - comp["Pavlab_Counts"]

  plot_altmann(comp, args)
  
  # Gene-wise correlations
  common_genes = author.index.intersection(pav_cpm.index)
  common_samples = author.columns.intersection(pav_cpm.columns)
  author_mat = author.loc[common_genes, common_samples]
  pavlab_mat = pav_cpm.loc[common_genes, common_samples]

  # Flatten to long format
  author_long = author_mat.stack().reset_index()
  pavlab_long = pavlab_mat.stack().reset_index()

  # Combine into a single DataFrame
  author_long.columns = ["Gene", "Individual_ID", "Author_Expression"]
  pavlab_long.columns = ["Gene", "Individual_ID", "PavLab_Expression"]

  merged = pd.merge(author_long, pavlab_long, on=["Gene", "Individual_ID"])
  merged = merged.merge(
      author_meta[["Individual_ID", "Cohort"]],
      on="Individual_ID",
      how="left"
  )

  # Calculate Spearman correlation
  # don't group by cohort, just calculate overall
  
  spearman_corr, _ = spearmanr(merged["Author_Expression"], merged["PavLab_Expression"])
  plot_corr(merged, spearman_corr)


if __name__ == "__main__":
  main()