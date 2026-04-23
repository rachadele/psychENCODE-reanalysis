import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
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


def plot_correlation_heatmap(combined_df, gemma_to_gemma_map, metric="spearman_log2FoldChange", f1_path=None, annotation_level="class", output_path=None):
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

    # Pivot to wide: rows = contrast_clean, columns = ct_author
    pivot = df.pivot_table(
        index="contrast_clean",
        columns="ct_author",
        values=metric,
        aggfunc="mean",
    )

    # Reorder rows in a logical order
    row_order = [c for c in ["SCZ", "BD", "ASD", "MDD", "PTSD", "WS", "Age", "Sex"] if c in pivot.index]
    pivot = pivot.loc[row_order]

    # --- assign pipeline class to each author label (columns) ---
    author_to_class = {k: v for k, v in gemma_to_gemma_map.items() if k in pivot.columns}
    classes_present = list(OrderedDict.fromkeys(gemma_to_gemma_map.values()))
    classes_present = [c for c in classes_present if c in author_to_class.values()]

    # colour palette for pipeline classes
    n_classes = len(classes_present)
    class_palette = dict(zip(classes_present, sns.color_palette("tab20", n_colors=n_classes)))

    # --- F1 score loading ---
    # Maps F1 benchmark labels -> pipeline class names used in gemma_to_gemma_map
    f1_label_to_class_map = {
        "class": {
            "Astrocyte":       "astrocyte",
            "Chandelier":      "chandelier.pvalb.GABAergic.cortical.interneuron",
            "Microglia":       "microglial.cell",
            "OPC":             "oligodendrocyte.precursor.cell",
            "Oligodendrocyte": "oligodendrocyte",
            "PAX6":            "PAX6.GABAergic.cortical.interneuron",
            "PVALB":           "pvalb.GABAergic.cortical.interneuron",
            "SNCG":            "sncg.GABAergic.cortical.interneuron",
            "SST":             "sst.GABAergic.cortical.interneuron",
            "VIP":             "vip.GABAergic.cortical.interneuron",
            "Vascular":        "vascular",
            "deep layer non-IT": "deep.layer.non.IT",
            "LAMP5":           "lamp5.GABAergic.cortical.interneuron",
            "L2/3-6 IT":       "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
        },
        "subclass": {
            "Astrocyte":       "astrocyte",
            "Chandelier":      "chandelier.pvalb.GABAergic.cortical.interneuron",
            "Endothelial":     "endothelial.cell",
            "Immune":          "microglial.cell",
            "Microglia":       "microglial.cell",
            "OPC":             "oligodendrocyte.precursor.cell",
            "Oligodendrocyte": "oligodendrocyte",
            "PAX6":            "PAX6.GABAergic.cortical.interneuron",
            "PVALB":           "pvalb.GABAergic.cortical.interneuron",
            "Pericyte":        "pericyte",
            "SNCG":            "sncg.GABAergic.cortical.interneuron",
            "SST":             "sst.GABAergic.cortical.interneuron",
            "SST Chodl":       "sst.GABAergic.cortical.interneuron",
            "VIP":             "vip.GABAergic.cortical.interneuron",
            "VLMC":            "vascular.leptomeningeal.cell",
            "LAMP5":           "lamp5.GABAergic.cortical.interneuron",
            "L2/3-6 IT":      "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
            "L6 IT Car3":      "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
            "L5 ET":           "L5.extratelencephalic.projecting.glutamatergic.cortical.neuron",
            "L5/6 NP":         "near.projecting.glutamatergic.cortical.neuron",
            "L6 CT":           "corticothalamic.projecting.glutamatergic.cortical.neuron",
            "L6b":             "L6b.glutamatergic.cortical.neuron",
        },
    }
    f1_label_to_class = f1_label_to_class_map.get(annotation_level, f1_label_to_class_map["class"])

    ct_to_f1 = {}
    if f1_path is not None:
        f1_df = pd.read_csv(f1_path, sep="\t")
        f1_df = f1_df[(f1_df["key"] == annotation_level) & (f1_df["metric"] == "f1_score")]
        class_to_f1 = {
            f1_label_to_class[row["label"]]: row["mean"]
            for _, row in f1_df.iterrows()
            if row["label"] in f1_label_to_class
        }
        ct_to_f1 = {
            ct: class_to_f1.get(gemma_to_gemma_map.get(ct, ""), np.nan)
            for ct in pivot.columns
        }

    # Sort columns: by F1 score descending (if available), else by pipeline class + mean
    if ct_to_f1:
        pivot = pivot[sorted(pivot.columns, key=lambda ct: ct_to_f1.get(ct, -1), reverse=True)]
    else:
        col_mean = pivot.mean(axis=0)
        sort_key = pivot.columns.map(
            lambda ct: (classes_present.index(author_to_class.get(ct, classes_present[0])), -col_mean.get(ct, 0))
        )
        pivot = pivot.iloc[:, sort_key.argsort()]

    # --- Column colour sidebars ---
    class_colors = [class_palette[author_to_class[ct]] for ct in pivot.columns]

    if ct_to_f1:
        f1_cmap = plt.cm.YlOrRd
        f1_vals = [ct_to_f1.get(ct, np.nan) for ct in pivot.columns]
        f1_norm = plt.Normalize(vmin=0.7, vmax=1.0)
        f1_colors = [
            f1_cmap(f1_norm(v)) if not np.isnan(v) else (0.85, 0.85, 0.85, 1.0)
            for v in f1_vals
        ]
        col_colors = [f1_colors, class_colors]
    else:
        col_colors = class_colors

    # --- clustermap (heatmap only; cbar drawn separately) ---
    g = sns.clustermap(
        pivot,
        row_cluster=False,
        col_cluster=False,
        col_colors=col_colors,
        cmap="RdYlBu_r",
        vmin=settings["vmin"],
        vmax=settings["vmax"],
        figsize=(max(16, len(pivot.columns) * 0.85), max(11, len(row_order) * 2.0)),
        linewidths=0.5,
        linecolor="white",
        cbar_pos=None,
        annot=False,
    )

    # Axis labels
    g.ax_heatmap.set_xlabel("Author Cell Type", fontsize=56, labelpad=22)
    g.ax_heatmap.set_ylabel("Contrast", fontsize=56, labelpad=22)
    g.ax_heatmap.tick_params(axis="x", labelsize=44, rotation=90)
    g.ax_heatmap.tick_params(axis="y", labelsize=46, rotation=0)

    if output_path is None:
        output_path = f"heatmap_{metric}.png"
    g.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(g.fig)

    tsv_path = output_path.replace(".png", ".tsv")
    pivot.to_csv(tsv_path, sep="\t")

    # --- separate legend/colorbar figure ---
    legend_fig = plt.figure(figsize=(26, 14))

    # Spearman / Jaccard colorbar (tall, thin, left side)
    cbar_ax = legend_fig.add_axes([0.05, 0.15, 0.035, 0.72])
    sm_main = plt.cm.ScalarMappable(
        cmap=plt.cm.RdYlBu_r,
        norm=plt.Normalize(vmin=settings["vmin"], vmax=settings["vmax"]),
    )
    sm_main.set_array([])
    legend_fig.colorbar(sm_main, cax=cbar_ax)
    cbar_ax.tick_params(labelsize=40)
    cbar_ax.set_ylabel(settings["cbar_label"], fontsize=48, labelpad=22)

    # F1 score colorbar
    if ct_to_f1:
        f1_ax = legend_fig.add_axes([0.19, 0.15, 0.035, 0.72])
        sm_f1 = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd, norm=plt.Normalize(vmin=0.7, vmax=1.0))
        sm_f1.set_array([])
        legend_fig.colorbar(sm_f1, cax=f1_ax)
        f1_ax.tick_params(labelsize=40)
        f1_ax.set_ylabel("F1 Score", fontsize=48, labelpad=22)

    # Pipeline Cell Type Class legend (right side)
    legend_handles = [
        mpatches.Patch(
            facecolor=class_palette[cls],
            label=cls.replace(".", " "),
        )
        for cls in classes_present
    ]
    legend_ax = legend_fig.add_axes([0.34, 0.02, 0.64, 0.96])
    legend_ax.axis("off")
    legend_ax.legend(
        handles=legend_handles,
        title="Pipeline Cell Type Class",
        loc="center left",
        fontsize=40,
        title_fontsize=48,
        frameon=True,
    )

    legend_path = output_path.replace(".png", "_legend.png")
    legend_fig.savefig(legend_path, dpi=150, bbox_inches="tight")
    plt.close(legend_fig)
