import os
import json
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
import matplotlib.pyplot as plt
import seaborn as sns
import re
from utils import *

def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--de_overlap_paths", type=str, nargs="+", default=[""])
  parser.add_argument("--annotation_level", type=str, default="class", help="Cell type annotation level: class or subclass")
  parser.add_argument("--f1_path", type=str, default=None, help="Path to label_metrics_stats_per_label.tsv for F1 score annotation")
  parser.add_argument("--params", type=str, required=True, help="Path to params JSON file containing gemma_to_gemma_map")
#  parser.add_argument("--logFC_threshold", type=float, default=0.0, help="Log fold change threshold for significance")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args


def main():
	# Plot box and whisker plot for overlaps per cell type
	args = parse_arguments()
	de_overlap_paths = args.de_overlap_paths
	annotation_level = args.annotation_level
	f1_path = args.f1_path
	with open(args.params) as f:
		params = json.load(f)
	gemma_to_gemma_map = params["gemma_to_gemma_map"]

	contrasts_to_remove = ["Intercept","PMI","ancestry"]

	matched_ct_dfs = []
	for path in de_overlap_paths:
		contrast = os.path.basename(path).split("_pairwise_overlap.tsv")[0]
		# regex match to all
		if any(re.search(rf"{rem}", contrast, re.IGNORECASE) for rem in contrasts_to_remove):
			continue
		de_overlap_df = pd.read_csv(path, sep="\t")
  # need to add annotation level argument
		de_overlap_df["contrast"] = contrast
		filtered_df = de_overlap_df[de_overlap_df.apply(cell_types_match, axis=1, gemma_to_gemma_map=gemma_to_gemma_map)]
		# save filtered df to tsv
		filtered_df.to_csv(f"filtered_{os.path.basename(path)}", index=False, sep="\t")
		matched_ct_dfs.append(filtered_df)
	# Combine all filtered DataFrames
	combined_df = pd.concat(matched_ct_dfs, ignore_index=True)
	combined_df.to_csv("combined_filtered_overlaps.tsv", index=False, sep="\t")
	plot_boxplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="jaccard_index")
	plot_stripplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="jaccard_index")
	plot_correlation_heatmap(combined_df, gemma_to_gemma_map, metric="jaccard_index", f1_path=f1_path, annotation_level=annotation_level)
	average_overlap, sd_overlap = compute_overall_averages(combined_df, metric="jaccard_index")
	# write overall average to file
	pd.DataFrame({
		'average_jaccard_index': [average_overlap],
		'sd_jaccard_index': [sd_overlap]
	}).to_csv("average_jaccard.tsv", index=False, sep="\t")

	# write per-cell-type averages to file
	per_celltype_averages = compute_average_per_celltype(combined_df, metric="jaccard_index")
	pd.DataFrame(per_celltype_averages).to_csv("per_celltype_average_jaccard.tsv", index=False, sep="\t")


if __name__ == "__main__":
  main()