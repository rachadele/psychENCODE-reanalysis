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
import matplotlib.pyplot as plt
import seaborn as sns
import re
from utils import *

gemma_to_gemma_map_class = {
		"Astro": "astrocyte",
		"Chandelier": "chandelier.pvalb.GABAergic.cortical.interneuron",
		"Endo": "vascular",
		"Immune": "microglial.cell",
		"L2.3_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L4_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L5.6_NP": "deep.layer.non.IT",
		"L5_ET": "deep.layer.non.IT",
		"L5_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L6_CT": "deep.layer.non.IT",
		"L6_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L6_IT_Car3": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L6b": "deep.layer.non.IT",
		"Lamp5": "lamp5.GABAergic.cortical.interneuron",
		"Lamp5_Lhx6": "lamp5.GABAergic.cortical.interneuron",
		"Micro": "microglial.cell",
		"OPC": "oligodendrocyte.precursor.cell",
		"Oligo": "oligodendrocyte",
		"PC": "vascular",
		"Pax6": "PAX6.GABAergic.cortical.interneuron",
		"Pvalb": "pvalb.GABAergic.cortical.interneuron",
		"SMC": "vascular",
		"Sncg": "sncg.GABAergic.cortical.interneuron",
		"Sst": "sst.GABAergic.cortical.interneuron",
		"Sst_Chodl": "sst.GABAergic.cortical.interneuron",
		"VLMC": "vascular",
		"Vip": "vip.GABAergic.cortical.interneuron"
	}

gemma_to_gemma_map_subclass = {
		"Astro": "astrocyte",
		"Chandelier": "chandelier.pvalb.GABAergic.cortical.interneuron",
		"Endo": "endothelial.cell",
		"Immune": "microglial.cell",
		"L2.3_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L4_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L5.6_NP": "near.projecting.glutamatergic.cortical.neuron",
		"L5_ET": "L5.extratelencephalic.projecting.glutamatergic.cortical.neuron",
		"L5_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L6_CT": "corticothalamic.projecting.glutamatergic.cortical.neuron",
		"L6_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L6_IT_Car3": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
		"L6b": "L6b.glutamatergic.cortical.neuron",
		"Lamp5": "lamp5.GABAergic.cortical.interneuron",
		"Lamp5_Lhx6": "lamp5.GABAergic.cortical.interneuron",
		"Micro": "microglial.cell",
		"OPC": "oligodendrocyte.precursor.cell",
		"Oligo": "oligodendrocyte",
		"PC": "pericyte",
		"Pax6": "PAX6.GABAergic.cortical.interneuron",
		"Pvalb": "pvalb.GABAergic.cortical.interneuron",
		"SMC": "smooth.muscle.cell",
		"Sncg": "sncg.GABAergic.cortical.interneuron",
		"Sst": "sst.GABAergic.cortical.interneuron",
		"Sst_Chodl": "sst.GABAergic.cortical.interneuron",
		"VLMC": "vascular.leptomeningeal.cell",
		"Vip": "vip.GABAergic.cortical.interneuron"
	}

def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--corr_tables", type=str, nargs="+", default=["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-subclass/all_corr/Age_death/files/Age_death_pairwise_corr.tsv"
                                                                     ,"/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-subclass/all_corr/Biological_Sex_male_vs_female/files/Biological_Sex_male_vs_female_pairwise_corr.tsv"])
  parser.add_argument("--annotation_level", type=str, default="class", help="Cell type annotation level: class or subclass")
#  parser.add_argument("--logFC_threshold", type=float, default=0.0, help="Log fold change threshold for significance")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args


def main():
	# Plot box and whisker plot for overlaps per cell type
	args = parse_arguments()
	corr_tables = args.corr_tables
	annotation_level = args.annotation_level
 # dynamically set mapping dict based on annotation level
	if annotation_level == "class":
		gemma_to_gemma_map = gemma_to_gemma_map_class
	else:
		gemma_to_gemma_map = gemma_to_gemma_map_subclass

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
	plot_boxplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_log2FoldChange")
	plot_boxplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_pvalue")
	plot_stripplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_log2FoldChange")
	plot_stripplot(combined_df, gemma_to_gemma_map=gemma_to_gemma_map, x="spearman_pvalue")
	plot_correlation_heatmap(combined_df, gemma_to_gemma_map)

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