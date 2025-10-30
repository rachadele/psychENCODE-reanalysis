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
from collections import OrderedDict
plt.rcParams.fontsize = 14
gemma_to_gemma_map = {
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



def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--corr_tables", type=str, nargs="+", default=["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-class/all_corr/Biological_Sex_male_vs_female/files/pairwise_corr_Biological_Sex_male_vs_female.tsv",
                                                                     "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-class/all_corr/Disorder_ASD_vs_Control/files/pairwise_corr_Disorder_ASD_vs_Control.tsv"])
#  parser.add_argument("--logFC_threshold", type=float, default=0.0, help="Log fold change threshold for significance")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

# Filter rows where cell types match according to the mapping
def cell_types_match(row):
    ct1 = row['ct_sc_pipeline']
    ct2 = row['ct_author']
    return gemma_to_gemma_map.get(ct2, "no match") == ct1  # returns True if they match and False otherwise


# Plot box and whisker plot of percent overlaps per cell type across all contrasts
def plot_average_corr_boxplot(df, output_path="spearman_boxplot_per_celltype.png"):

	plt.figure(figsize=(10, 14))
	# Get cell type order from gemma_to_gemma_map values
	celltype_order = list(OrderedDict.fromkeys(gemma_to_gemma_map.values()))
	# Only keep cell types present in the data
	celltype_order = [ct for ct in celltype_order if ct in df['ct_sc_pipeline'].unique()]
	ax = sns.boxplot(y='ct_sc_pipeline', x='spearman_correlation', data=df, hue='ct_sc_pipeline', dodge=False, showmeans=False, order=celltype_order)
	# Overlay mean as a dot
	means = df.groupby('ct_sc_pipeline')['spearman_correlation'].mean()
	for i, celltype in enumerate(celltype_order):
		if celltype in means:
			ax.scatter(means[celltype], i, color='black', s=60, zorder=10, label='Mean' if i == 0 else "")
	plt.yticks(fontsize=10)
	plt.xlabel('Spearman Correlation')
	plt.ylabel('Cell Type')
	plt.title('DE Reproducibility Across All Contrasts', fontsize=16)
	# Only show legend once
	handles, labels = ax.get_legend_handles_labels()
	if 'Mean' in labels:
		plt.legend([handles[labels.index('Mean')]], ['Mean'], loc='lower right')
	else:
		plt.legend([],[], frameon=False)
	plt.tight_layout()
	plt.savefig(output_path)
	plt.close()

def average_correlation(corr_matrix):
    # Exclude diagonal if present
    mask = ~pd.Series(corr_matrix.index).isin([corr_matrix.columns])
    # Average all off-diagonal correlations
    values = corr_matrix.values
    n = values.shape[0]
    # Exclude diagonal
    avg_corr = (values.sum() - np.trace(values)) / (n * (n - 1))
    return avg_corr

def main():
	# Plot box and whisker plot for overlaps per cell type
	args = parse_arguments()
	corr_tables = args.corr_tables

	contrasts_to_remove = ["Intercept","PMI","ancestry"]

	matched_ct_dfs = []
	for path in corr_tables:
		contrast = os.path.basename(path).split("_pairwise_corr")[0]
		# regex match to all
		if any(re.search(rf"{rem}", contrast, re.IGNORECASE) for rem in contrasts_to_remove):
			continue
		pairwise_corr_df = pd.read_csv(path, sep="\t")
		filtered_df = pairwise_corr_df[pairwise_corr_df.apply(cell_types_match, axis=1)]
		# save filtered df to tsv
		#filtered_df.to_csv(f"filtered_{os.path.basename(path)}", index=False, sep="\t")
		matched_ct_dfs.append(filtered_df)
	# Combine all filtered DataFrames
	combined_df = pd.concat(matched_ct_dfs, ignore_index=True)
	plot_average_corr_boxplot(combined_df)


if __name__ == "__main__":
  main()