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
  parser.add_argument("--corr_tables", type=str, nargs="+", default=["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-subclass/all_corr/Age_death/files/Age_death_pairwise_corr.tsv"
                                                                     ,"/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-subclass/all_corr/Biological_Sex_male_vs_female/files/Biological_Sex_male_vs_female_pairwise_corr.tsv"])
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
def compute_average_per_celltype(combined_df, group_col="spearman_log2FoldChange"):
	result = combined_df.groupby('ct_sc_pipeline')[group_col].agg(['mean', 'std', 'count']).reset_index()
	result.rename(columns={'mean': 'average_' + group_col, 'std': 'sd_' + group_col, 'count': 'n_' + group_col}, inplace=True)
	return result

def compute_overall_averages(combined_df):
	# Calculate average percent overlap for matched cell types
	average_logFC = combined_df['spearman_log2FoldChange'].mean()
	sd_logFC = combined_df['spearman_log2FoldChange'].std()

	average_pval = combined_df['spearman_pvalue'].mean()
	sd_pval = combined_df['spearman_pvalue'].std()
	# write overall average to file
	return average_logFC, sd_logFC, average_pval, sd_pval


def compute_average_corr_per_celltype(combined_df):
    result = combined_df.groupby('ct_sc_pipeline')[['spearman_log2FoldChange', 'spearman_pvalue']].agg(['mean', 'std', 'count']).reset_index()
    # Flatten columns
    result.columns = ['ct_sc_pipeline', 'avg_log2FoldChange', 'std_log2FoldChange', 'n_log2FoldChange', 'avg_pvalue', 'std_pvalue', 'n_pvalue']
    return result

def plot_corr_boxplots(df, annotation_level="class", x="spearman_log2FoldChange"):
    if annotation_level == "class":
        gemma_to_gemma_map = gemma_to_gemma_map_class
    elif annotation_level == "subclass":
        gemma_to_gemma_map = gemma_to_gemma_map_subclass
    else:
        raise ValueError("annotation_level must be 'class' or 'subclass'")
    plt.figure(figsize=(10, 14))
    from collections import OrderedDict
    celltype_order = list(OrderedDict.fromkeys(df['ct_sc_pipeline'].unique()))
    # Boxplot for log2FoldChange
    ax1 = sns.boxplot(y='ct_sc_pipeline', x=x, data=df, order=celltype_order)
    means = df.groupby('ct_sc_pipeline')[x].mean()
    for i, celltype in enumerate(celltype_order):
        if celltype in means:
            ax1.scatter(means[celltype], i, color='black', s=60, zorder=10, label='Mean' if i == 0 else "")
    plt.yticks(fontsize=10)
    plt.xlabel(f'Spearman {x}')
    plt.ylabel('Cell Type')
    plt.title(f'Spearman Correlation ({x}) Across Contrasts')
    plt.tight_layout()
    plt.savefig(f'boxplot_spearman_{x}.png')
    plt.close()
   
def plot_corr_stripplot(df, annotation_level="class", x="spearman_log2FoldChange"):
    if annotation_level == "class":
        gemma_to_gemma_map = gemma_to_gemma_map_class
    elif annotation_level == "subclass":
        gemma_to_gemma_map = gemma_to_gemma_map_subclass
    else:
        raise ValueError("annotation_level must be 'class' or 'subclass'")
    plt.figure(figsize=(15, 15))
    from collections import OrderedDict
    celltype_order = list(OrderedDict.fromkeys(df['ct_sc_pipeline'].unique()))
    # make points larger
    ax = sns.stripplot(y='ct_sc_pipeline', x=x, data=df, order=celltype_order, hue='contrast', dodge=True, jitter=True, size=5)
    plt.yticks(fontsize=10)
    plt.xlabel(f'Spearman {x}')
    plt.ylabel('Cell Type')
    plt.title(f'Spearman Correlation ({x}) Across Contrasts (Strip Plot)')
    plt.legend(title='Contrast', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'stripplot_spearman_{x}.png')
    plt.close()


def main():
	# Plot box and whisker plot for overlaps per cell type
	args = parse_arguments()
	corr_tables = args.corr_tables
	annotation_level = args.annotation_level

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
		filtered_df = corr_df[corr_df.apply(cell_types_match, axis=1, annotation_level=annotation_level)]
		# save filtered df to tsv
		filtered_df.to_csv(f"filtered_{os.path.basename(path)}", index=False, sep="\t")
		matched_ct_dfs.append(filtered_df)
	# Combine all filtered DataFrames
	combined_df = pd.concat(matched_ct_dfs, ignore_index=True)
	plot_corr_boxplots(combined_df, annotation_level=annotation_level, x="spearman_log2FoldChange")
	plot_corr_boxplots(combined_df, annotation_level=annotation_level, x="spearman_pvalue")
	plot_corr_stripplot(combined_df, annotation_level=annotation_level, x="spearman_log2FoldChange")
	plot_corr_stripplot(combined_df, annotation_level=annotation_level, x="spearman_pvalue")

	average_logFC, sd_logFC, average_pval, sd_pval = compute_overall_averages(combined_df)
	# write overall average to file
	pd.DataFrame({
		'average_logFC_corr': [average_logFC],
		'sd_logFC_corr': [sd_logFC],
		'average_pval_corr': [average_pval],
		'sd_pval_corr': [sd_pval]
	}).to_csv("average_corr_pvalue_log2FC.tsv", index=False, sep="\t")

	# write per-cell-type averages to file
	per_celltype_averages_fc = compute_average_per_celltype(combined_df, group_col="spearman_log2FoldChange")
	pd.DataFrame(per_celltype_averages_fc).to_csv("per_celltype_average_spearman_log2FoldChange.tsv", index=False, sep="\t")
 
 
	per_celltype_averages_pval = compute_average_per_celltype(combined_df, group_col="spearman_pvalue")
	pd.DataFrame(per_celltype_averages_pval).to_csv("per_celltype_average_spearman_pvalue.tsv", index=False, sep="\t")


	avg_corrs = compute_average_corr_per_celltype(combined_df)
	avg_corrs.to_csv("average_spearman_corrs.tsv", index=False, sep="\t")


if __name__ == "__main__":
  main()