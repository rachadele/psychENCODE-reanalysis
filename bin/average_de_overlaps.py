
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

gemma_to_gemma_map = {
	"Chandelier": "chandelier.pvalb.GABAergic.cortical.interneuron",
	"Lamp5": "lamp5.GABAergic.cortical.interneuron",
	"Lamp5_Lhx6": "lamp5.GABAergic.cortical.interneuron",
	"Pax6": "caudal.ganglionic.eminence.derived.GABAergic.cortical.interneuron",
	"Pvalb": "pvalb.GABAergic.cortical.interneuron",
	"Sncg": "sncg.GABAergic.cortical.interneuron",
	"Sst": "sst.GABAergic.cortical.interneuron",
	"Sst_Chodl": "sst.GABAergic.cortical.interneuron",
	"Vip": "vip.GABAergic.cortical.interneuron",
	"L2_3_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
	"L4_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
	"L5_6_NP": "near.projecting.glutamatergic.cortical.neuron",
	"L5_ET": "L5.extratelencephalic.projecting.glutamatergic.cortical.neuron",
	"L5_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
	"L6_CT": "L6.corticothalamic.projecting.glutamatergic.cortical.neuron",
	"L6_IT": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
	"L6_IT_Car3": "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
	"L6b": "L6b.glutamatergic.cortical.neuron",
	"Astro": "astrocyte",
	"Endo": "endothelial.cell",
	"Immune": "no match",
	"Micro": "microglial.cell",
	"Micro_PVM": "microglial.cell",
	"Oligo": "oligodendrocyte",
	"OPC": "oligodendrocyte.precursor.cell",
	"PC": "pericyte",
	"SMC": "no match",
	"VLMC": "vascular.leptomeningeal.cell"
}

def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--de_overlap_paths", type=str, nargs="+", default=[""])
#  parser.add_argument("--logFC_threshold", type=float, default=0.0, help="Log fold change threshold for significance")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

# Filter rows where cell types match according to the mapping
def cell_types_match(row):
    ct1 = row['ct_sc_pipeline']
    ct2 = row['ct_author']
    return gemma_to_gemma_map.get(ct2, "no match") == ct1  # returns True if they match and False otherwise

# Compute average overlap per cell type across all contrasts
def compute_average_per_celltype(combined_df):
	result = combined_df.groupby('ct_sc_pipeline')['percent_overlap'].agg(['mean', 'std', 'count']).reset_index()
	result.rename(columns={'mean': 'average_percent_overlap', 'std': 'sd_percent_overlap', 'count': 'n_pairs'}, inplace=True)
	return result

def compute_overall_averages(combined_df):
	# Calculate average percent overlap for matched cell types
	average_overlap = combined_df['percent_overlap'].mean()
	sd_overlap = combined_df['percent_overlap'].std()
	# write overall average to file
	return average_overlap, sd_overlap

# Plot box and whisker plot of percent overlaps per cell type across all contrasts
def plot_overlap_boxplot(df, output_path="overlap_boxplot_per_celltype.png"):

	plt.figure(figsize=(10, 14))
	# Get cell type order from gemma_to_gemma_map values
	from collections import OrderedDict
	celltype_order = list(OrderedDict.fromkeys(gemma_to_gemma_map.values()))
	# Only keep cell types present in the data
	celltype_order = [ct for ct in celltype_order if ct in df['ct_sc_pipeline'].unique()]
	ax = sns.boxplot(y='ct_sc_pipeline', x='percent_overlap', data=df, hue='ct_sc_pipeline', dodge=False, showmeans=False, order=celltype_order)
	# Overlay mean as a dot
	means = df.groupby('ct_sc_pipeline')['percent_overlap'].mean()
	for i, celltype in enumerate(celltype_order):
		if celltype in means:
			ax.scatter(means[celltype], i, color='black', s=60, zorder=10, label='Mean' if i == 0 else "")
	plt.yticks(fontsize=10)
	plt.xlabel('Jaccard Index')
	plt.ylabel('Cell Type')
	plt.title('DE Reproducibility Across All Contrasts')
	# Only show legend once
	handles, labels = ax.get_legend_handles_labels()
	if 'Mean' in labels:
		plt.legend([handles[labels.index('Mean')]], ['Mean'], loc='lower right')
	else:
		plt.legend([],[], frameon=False)
	plt.tight_layout()
	plt.savefig(output_path)
	plt.close()

def main():
	# Plot box and whisker plot for overlaps per cell type
	args = parse_arguments()
	de_overlap_paths = args.de_overlap_paths

	contrasts_to_remove = ["Intercept","PMI","ancestry"]

	matched_ct_dfs = []
	for path in de_overlap_paths:
		contrast = os.path.basename(path).split("_pairwise_overlap.tsv")[0]
		# regex match to all
		if any(re.search(rf"{rem}", contrast, re.IGNORECASE) for rem in contrasts_to_remove):
			continue
		de_overlap_df = pd.read_csv(path, sep="\t")
		filtered_df = de_overlap_df[de_overlap_df.apply(cell_types_match, axis=1)]
		# save filtered df to tsv
		filtered_df.to_csv(f"filtered_{os.path.basename(path)}", index=False, sep="\t")
		matched_ct_dfs.append(filtered_df)
	# Combine all filtered DataFrames
	combined_df = pd.concat(matched_ct_dfs, ignore_index=True)
	plot_overlap_boxplot(combined_df)

	average_overlap, sd_overlap = compute_overall_averages(combined_df)
	# write overall average to file
	pd.DataFrame({
		'average_percent_overlap': [average_overlap],
		'sd_percent_overlap': [sd_overlap]
	}).to_csv("average_de_overlaps.tsv", index=False, sep="\t")

	# write per-cell-type averages to file
	per_celltype_averages = compute_average_per_celltype(combined_df)
	pd.DataFrame(per_celltype_averages).to_csv("per_celltype_averages.tsv", index=False, sep="\t")


if __name__ == "__main__":
  main()