import os
import argparse
import glob
import pandas as pd
from functools import reduce
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")


def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--author_metadata", type=str, nargs="+", default ="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/PEC2_sample_metadata.txt")
  parser.add_argument("--metadata_files", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata",
                      help="Path to the metadata files from Gemma")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

  
def plot_summary(full_meta, grouping_col="Cohort_author"):

  cohort_summary = full_meta.groupby(grouping_col).agg({
      "in_both": "sum",
      "only_author": "sum",
      "only_gemma": "sum"
  }).rename(columns={
      "in_both": "In Both",
      "only_author": "Only Author",
      "only_gemma": "Only Gemma"
  })

  # Step 2: Set the index as Cohort (optional, but handy)
  cohort_summary.index.name = grouping_col.replace("_author","")
  # drop columns with all zero values
  cohort_summary = cohort_summary.loc[:, (cohort_summary != 0).any(axis=0)]

  # Now cohort_summary is a DataFrame with Cohort as index and columns as sources
  # You can plot directly with:
  cohort_summary.plot(
      kind="bar",
      stacked=True,
      figsize=(12, 6),
      colormap="Set2",
      edgecolor="black"

  )
  
  newname= grouping_col.replace("_author", "")
  plt.title(f"Comparison of Author and Gemma Samples by {newname}")
  plt.xlabel(newname)
  plt.ylabel("Number of Samples")
  plt.xticks(rotation=90)
  plt.tight_layout()
  plt.savefig(f"author_gemma_comparison_{newname}.png")
  plt.close()
  
   
  
  
def main():
  args = parse_arguments()
  author_metadata = pd.read_csv(args.author_metadata, sep="\t", dtype=str)
  metadata_files = glob.glob(os.path.join(args.metadata_files, "*_sample_meta.tsv"))
  gemma_meta = pd.DataFrame()
  for file in metadata_files:
    temp = pd.read_csv(file, sep="\t", dtype=str)
    gemma_meta = pd.concat([gemma_meta, temp], ignore_index=True)
     
  # summary statistics of differences in "Individual_ID" grouped by "Cohort"
  author_samples =  author_metadata["Individual_ID"].unique()
  gemma_samples = gemma_meta["Individual_ID"].unique()
  intersection = set(author_samples) & set(gemma_samples)
  only_author = set(author_samples) - set(gemma_samples)
  only_gemma = set(gemma_samples) - set(author_samples)
  
  full_meta = pd.merge(author_metadata, gemma_meta, on="Individual_ID", how="outer", suffixes=("_author", "_gemma"))
  full_meta["in_both"] = full_meta["Individual_ID"].isin(intersection)
  full_meta["only_author"] = full_meta["Individual_ID"].isin(only_author)
  full_meta["only_gemma"] = full_meta["Individual_ID"].isin(only_gemma)
 
	# Group by both Cohort and Disorder
  cohort_disorder_summary = full_meta.groupby(["Cohort_author", "Disorder_author"]).agg({
			"in_both": "sum",
			"only_author": "sum",
			"only_gemma": "sum"
	}).reset_index()

	# Rename columns for clarity
  cohort_disorder_summary.columns = ["Cohort", "Disorder", "In Both", "Only Author", "Only Gemma"]

  # write table
  cohort_disorder_summary.to_csv("author_gemma_differences.tsv", index=False, sep="\t")
 
  plot_summary(full_meta, grouping_col="Cohort_author")
  plot_summary(full_meta, grouping_col="Disorder_author")



if __name__ == "__main__":
  main()
