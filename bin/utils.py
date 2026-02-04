import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
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


def plot_correlation_heatmap(combined_df, gemma_to_gemma_map, metric="spearman_log2FoldChange", output_path=None):
    """Plot a heatmap of a metric (ct_author x contrast).

    Parameters
    ----------
    metric : str
        Column name in combined_df to plot, e.g.
        "spearman_log2FoldChange", "spearman_pvalue", or "jaccard_index".
    """

    # --- metric-specific display settings ---
    metric_settings = {
        "spearman_log2FoldChange": {
            "cbar_label": "Spearman rho (log2FC)",
            "title": "Spearman log2FC Correlation: Author Cell Types x Contrasts",
            "vmin": 0, "vmax": 1, "fmt": ".2f",
        },
        "spearman_pvalue": {
            "cbar_label": "Spearman rho (p-value)",
            "title": "Spearman P-value Correlation: Author Cell Types x Contrasts",
            "vmin": 0, "vmax": 1, "fmt": ".2f",
        },
        "jaccard_index": {
            "cbar_label": "Jaccard Index (%)",
            "title": "Jaccard Index: Author Cell Types x Contrasts",
            "vmin": 0, "vmax": 100, "fmt": ".1f",
        },
    }
    settings = metric_settings.get(metric, {
        "cbar_label": metric.replace("_", " ").title(),
        "title": f"{metric.replace('_', ' ').title()}: Author Cell Types x Contrasts",
        "vmin": 0, "vmax": 1, "fmt": ".2f",
    })

    # --- contrast display names ---
    contrast_display = {
        "Disorder_Schizophrenia_vs_Control": "SCZ",
        "Disorder_Bipolar.Disorder_vs_Control": "BD",
        "Disorder_ASD_vs_Control": "ASD",
        "Disorder_MDD_vs_Control": "MDD",
        "Disorder_PTSD_vs_Control": "PTSD",
        "Disorder_Williams.Syndrome_vs_Control": "WS",
        "Age_death": "Age",
        "Biological_Sex_male_vs_female": "Sex",
    }
    # Also handle contrast names that include the file suffix (jaccard files)
    contrast_display_extra = {
        k + "_pairwise_jaccard.tsv": v for k, v in contrast_display.items()
    }
    contrast_display_all = {**contrast_display, **contrast_display_extra}

    # Keep only rows where ct_author is in the mapping
    df = combined_df[combined_df["ct_author"].isin(gemma_to_gemma_map.keys())].copy()

    # Clean contrast names
    df["contrast_clean"] = df["contrast"].map(contrast_display_all)
    df = df.dropna(subset=["contrast_clean"])

    # Pivot to wide: rows = ct_author, columns = contrast_clean
    pivot = df.pivot_table(
        index="ct_author",
        columns="contrast_clean",
        values=metric,
        aggfunc="mean",
    )

    # Reorder columns in a logical order
    col_order = [c for c in ["SCZ", "BD", "ASD", "MDD", "PTSD", "WS", "Age", "Sex"] if c in pivot.columns]
    pivot = pivot[col_order]

    # --- assign pipeline class to each author label ---
    author_to_class = {k: v for k, v in gemma_to_gemma_map.items() if k in pivot.index}
    classes_present = list(OrderedDict.fromkeys(gemma_to_gemma_map.values()))
    classes_present = [c for c in classes_present if c in author_to_class.values()]

    # colour palette for pipeline classes
    n_classes = len(classes_present)
    class_palette = dict(zip(classes_present, sns.color_palette("tab20", n_colors=n_classes)))

    # Sort rows: group by pipeline class, then by mean value within group (descending)
    row_mean = pivot.mean(axis=1)
    sort_key = pivot.index.map(
        lambda ct: (classes_present.index(author_to_class.get(ct, classes_present[0])), -row_mean.get(ct, 0))
    )
    pivot = pivot.iloc[sort_key.argsort()]

    # Row colour sidebar
    row_colors = pivot.index.map(lambda ct: class_palette[author_to_class[ct]])
    row_colors = list(row_colors)

    # --- clustermap ---
    g = sns.clustermap(
        pivot,
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        cmap="RdYlBu_r",
        vmin=settings["vmin"],
        vmax=settings["vmax"],
        figsize=(max(10, len(col_order) * 1.6), max(10, len(pivot) * 0.55)),
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": settings["cbar_label"], "shrink": 0.5},
        annot=True,
        fmt=settings["fmt"],
        annot_kws={"size": 9},
    )

    # Axis labels
    g.ax_heatmap.set_ylabel("Author Cell Type", fontsize=14)
    g.ax_heatmap.set_xlabel("Contrast", fontsize=14)
    g.ax_heatmap.tick_params(axis="y", labelsize=11)
    g.ax_heatmap.tick_params(axis="x", labelsize=12)

    # --- legend for pipeline classes ---
    legend_handles = [
        mpatches.Patch(
            facecolor=class_palette[cls],
            label=cls.replace(".", " "),
        )
        for cls in classes_present
    ]
    g.ax_heatmap.legend(
        handles=legend_handles,
        title="Pipeline Cell Type Class",
        bbox_to_anchor=(1.35, 1.0),
        loc="upper left",
        fontsize=8,
        title_fontsize=9,
        frameon=True,
    )

    g.fig.suptitle(settings["title"], fontsize=16, y=1.02)
    g.fig.tight_layout(rect=[0, 0, 1, 1])

    if output_path is None:
        output_path = f"heatmap_{metric}.png"
    g.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
