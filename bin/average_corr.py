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
  parser.add_argument("--corr_tables", type=str, nargs="+", default=["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-subclass/all_corr/Age_death/files/Age_death_pairwise_corr.tsv"
                                                                     ,"/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-subclass/all_corr/Biological_Sex_male_vs_female/files/Biological_Sex_male_vs_female_pairwise_corr.tsv"])
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
	corr_tables = args.corr_tables
	annotation_level = args.annotation_level
	f1_path = args.f1_path
	with open(args.params) as f:
		params = json.load(f)
	gemma_to_gemma_map = params["gemma_to_gemma_map"]

	contrasts_to_remove = ["Intercept","PMI","ancestry"]

	matched_ct_dfs = []
	for path in corr_tables:
		print(path)
		contrast = os.path.basename(path).split("_pairwise_corr.tsv")[0]
		# regex match to all
		if any(re.search(rf"{rem}", contrast, re.IGNORECASE) for rem in contrasts_to_remove):
			continue
		corr_df = pd.read_csv(path, sep="\t")
  # need to add annotation level argument
		corr_df["contrast"] = contrast
		filtered_df = corr_df[corr_df.apply(cell_types_match, axis=1, gemma_to_gemma_map=gemma_to_gemma_map)]
		# save filtered df to tsv
		filtered_df.to_csv(f"filtered_{os.path.basename(path)}", index=False, sep="\t")
		matched_ct_dfs.append(filtered_df)
	# Combine all filtered DataFrames
	combined_df = pd.concat(matched_ct_dfs, ignore_index=True)
	combined_df.to_csv("combined_filtered_corr.tsv", index=False, sep="\t")
	plot_boxplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_log2FoldChange")
	plot_boxplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_pvalue")
	plot_stripplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_log2FoldChange")
	plot_stripplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_pvalue")
	plot_correlation_heatmap(combined_df, gemma_to_gemma_map, metric="spearman_log2FoldChange", f1_path=f1_path, annotation_level=annotation_level)
	plot_correlation_heatmap(combined_df, gemma_to_gemma_map, metric="spearman_pvalue", f1_path=f1_path, annotation_level=annotation_level)

	average_logFC, sd_logFC = compute_overall_averages(combined_df, metric="spearman_log2FoldChange")
	average_pval, sd_pval = compute_overall_averages(combined_df, metric="spearman_pvalue")
	# write overall average to file
	pd.DataFrame({
		'average_logFC_corr': [average_logFC],
		'sd_logFC_corr': [sd_logFC],
		'average_pval_corr': [average_pval],
		'sd_pval_corr': [sd_pval]
	}).to_csv("average_corr_pvalue_log2FC.tsv", index=False, sep="\t")

	# write per-cell-type averages to file
	per_celltype_averages_fc = compute_average_per_celltype(combined_df, metric="spearman_log2FoldChange")
	pd.DataFrame(per_celltype_averages_fc).to_csv("per_celltype_average_spearman_log2FoldChange.tsv", index=False, sep="\t")
 
 
	per_celltype_averages_pval = compute_average_per_celltype(combined_df, metric="spearman_pvalue")
	pd.DataFrame(per_celltype_averages_pval).to_csv("per_celltype_average_spearman_pvalue.tsv", index=False, sep="\t")

if __name__ == "__main__":
  main()