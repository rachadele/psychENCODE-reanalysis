import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict

def reverse_lookup(mapping, value):
    """Return all keys that map to the given value."""
    return [k for k, v in mapping.items() if v == value]

def cell_types_match(row, gemma_to_gemma_map):
    ct1 = row['ct_author']
    ct2 = row['ct_sc_pipeline']
    return ct1 in reverse_lookup(gemma_to_gemma_map, ct2)
     

def plot_boxplot(df, gemma_to_gemma_map, x="percent_overlap", outpath=None):
    plt.figure(figsize=(20, 20))
    # Create a copy with cleaned cell type names (replace "." with " ")
    df_plot = df.copy()
    df_plot['ct_sc_pipeline_clean'] = df_plot['ct_sc_pipeline'].str.replace('.', ' ', regex=False)
    celltype_order = list(OrderedDict.fromkeys(gemma_to_gemma_map.values()))
    celltype_order = [ct for ct in celltype_order if ct in df['ct_sc_pipeline'].unique()]
    celltype_order_clean = [ct.replace('.', ' ') for ct in celltype_order]
    ax = sns.boxplot(y='ct_sc_pipeline_clean', x=x, data=df_plot, hue='ct_sc_pipeline_clean', dodge=False, showmeans=False, order=celltype_order_clean)
    means = df_plot.groupby('ct_sc_pipeline_clean')[x].mean()
    for i, celltype in enumerate(celltype_order_clean):
        if celltype in means:
            ax.scatter(means[celltype], i, color='black', s=100, zorder=10, label='Mean' if i == 0 else "")
    plt.yticks(fontsize=28)
    plt.xticks(fontsize=28)
    plt.xlabel(x.replace('_', ' ').title(), fontsize=32)
    plt.ylabel('Cell Type', fontsize=32)
    plt.title(f'{x.replace("_", " ").title()} Across All Contrasts', fontsize=32)
    handles, labels = ax.get_legend_handles_labels()
    if 'Mean' in labels:
        plt.legend([handles[labels.index('Mean')]], ['Mean'], loc='lower right', fontsize=28)
    else:
        plt.legend([],[], frameon=False)
    plt.tight_layout()
    if outpath is None:
        outpath = f'boxplot_{x}.png'
    plt.savefig(outpath, dpi=150)
    plt.close()

def plot_stripplot(df, gemma_to_gemma_map, x="percent_overlap", output_path=None):
    plt.figure(figsize=(36, 28))
    # Create a copy with cleaned cell type names (replace "." with " ")
    df_plot = df.copy()
    df_plot['ct_sc_pipeline_clean'] = df_plot['ct_sc_pipeline'].str.replace('.', ' ', regex=False)
    celltype_order = list(OrderedDict.fromkeys(df['ct_sc_pipeline'].unique()))
    celltype_order_clean = [ct.replace('.', ' ') for ct in celltype_order]
    ax = sns.stripplot(y='ct_sc_pipeline_clean', x=x, data=df_plot, order=celltype_order_clean, hue='contrast', dodge=False, jitter=True, size=18)
    plt.yticks(fontsize=28)
    plt.xticks(fontsize=28)
    plt.xlabel(x.replace('_', ' ').title(), fontsize=32)
    plt.ylabel('Cell Type', fontsize=32)
    plt.title(f'{x.replace("_", " ").title()} Across Contrasts', fontsize=32)
    plt.legend(title='Contrast', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=24, title_fontsize=28)
    plt.tight_layout()
    if output_path is None:
        output_path = f'stripplot_{x}.png'
    plt.savefig(output_path, dpi=150)
    plt.close()

# Compute average overlap per cell type across all contrasts
def compute_average_per_celltype(combined_df, metric):
    result = combined_df.groupby('ct_sc_pipeline')[metric].agg(['mean', 'std', 'count']).reset_index()
    result.rename(columns={'mean': f'average_{metric}', 'std': f'sd_{metric}', 'count': 'n_pairs'}, inplace=True)
    # sort from lowest to highest average percent overlap
    result = result.sort_values(by=f'average_{metric}', ascending=False)
    return result

def compute_overall_averages(combined_df, metric):
	# Calculate average percent overlap for matched cell types
	average_overlap = combined_df[metric].mean()
	sd_overlap = combined_df[metric].std()
	# write overall average to file
	return average_overlap, sd_overlap
