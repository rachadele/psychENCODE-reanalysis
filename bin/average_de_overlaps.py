
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
  parser.add_argument("--de_overlap_paths", type=str, nargs="+", default=[""])
  parser.add_argument("--annotation_level", type=str, default="class", help="Cell type annotation level: class or subclass")
#  parser.add_argument("--logFC_threshold", type=float, default=0.0, help="Log fold change threshold for significance")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def reverse_lookup(mapping, value):
    # need to return all keys that map to the given value
    
    return [k for k, v in mapping.items() if v == value]

# Filter rows where cell types match according to the mapping
def cell_types_match(row, annotation_level):
    ct1 = row['ct_author']
    ct2 = row['ct_sc_pipeline']
    if annotation_level == "class":
        # fix this to reverse match
        return ct1 in reverse_lookup(gemma_to_gemma_map_class, ct2)
    elif annotation_level == "subclass":
        return ct1 in reverse_lookup(gemma_to_gemma_map_subclass, ct2)
    return False

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
def plot_overlap_boxplot(df, output_path="overlap_boxplot_per_celltype.png", annotation_level="class"):

	if annotation_level == "class":
		gemma_to_gemma_map = gemma_to_gemma_map_class
	elif annotation_level == "subclass":
		gemma_to_gemma_map = gemma_to_gemma_map_subclass
	else:
		raise ValueError("annotation_level must be 'class' or 'subclass'")
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
	annotation_level = args.annotation_level

	contrasts_to_remove = ["Intercept","PMI","ancestry"]

	matched_ct_dfs = []
	for path in de_overlap_paths:
		contrast = os.path.basename(path).split("_pairwise_overlap.tsv")[0]
		# regex match to all
		if any(re.search(rf"{rem}", contrast, re.IGNORECASE) for rem in contrasts_to_remove):
			continue
		de_overlap_df = pd.read_csv(path, sep="\t")
  # need to add annotation level argument
  
		filtered_df = de_overlap_df[de_overlap_df.apply(cell_types_match, axis=1, annotation_level=annotation_level)]
		# save filtered df to tsv
		filtered_df.to_csv(f"filtered_{os.path.basename(path)}", index=False, sep="\t")
		matched_ct_dfs.append(filtered_df)
	# Combine all filtered DataFrames
	combined_df = pd.concat(matched_ct_dfs, ignore_index=True)
	plot_overlap_boxplot(combined_df, annotation_level=annotation_level)

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