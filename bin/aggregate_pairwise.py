import os
import json
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata, spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list


def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute pairwise DE correlations between GEMMA and author cell types")
    parser.add_argument("--contrast", type=str, default="Schizophrenia")
    parser.add_argument("--pavlab_paths", type=str, nargs="+", required=True)
    parser.add_argument("--author_paths", type=str, nargs="+", required=True)
    parser.add_argument("--params", type=str, required=True, help="Path to params JSON with gemma_to_gemma_map")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def label_significance(row):
    if row["padj_author"] < 0.05 and row["padj_gemma"] < 0.05:
        return "Both"
    elif row["padj_author"] < 0.05:
        return "Author Only"
    elif row["padj_gemma"] < 0.05:
        return "GEMMA Only"
    else:
        return "Not Significant"


def scatter_plot(merged_df, contrast, author_ct, gemma_ct, corr_logfc):
    palette = {
        "Not Significant": "lightgray",
        "Author Only": "royalblue",
        "GEMMA Only": "darkorange",
        "Both": "crimson",
    }
    hue_order = ["Not Significant", "Author Only", "GEMMA Only", "Both"]

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(
        data=merged_df, x="log2FoldChange_author", y="log2FoldChange_gemma",
        hue="Significance", hue_order=hue_order, palette=palette,
        alpha=0.6, edgecolor=None, s=60, ax=ax,
    )
    sns.regplot(
        data=merged_df, x="log2FoldChange_author", y="log2FoldChange_gemma",
        scatter=False, line_kws={"color": "black", "linewidth": 1.5}, ax=ax,
    )
    ax.text(0.05, 0.95, f"Spearman ρ = {corr_logfc:.3f}",
            transform=ax.transAxes, fontsize=13, va="top")
    ax.set_xlabel(f"Author log2FC  ({author_ct})", fontsize=13)
    ax.set_ylabel(f"GEMMA log2FC  ({gemma_ct})", fontsize=13)
    ax.set_title(f"{contrast}  |  {author_ct} → {gemma_ct}", fontsize=14)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    safe_author = author_ct.replace("/", "_").replace(" ", "_")
    fig.savefig(f"{contrast}_{safe_author}_scatter.png", dpi=150)
    plt.close(fig)


def main():
    args = parse_arguments()

    with open(args.params) as f:
        params = json.load(f)
    gemma_to_gemma_map = params["gemma_to_gemma_map"]

    pavlab_de_results = {}
    for path in args.pavlab_paths:
        df = pd.read_csv(path, sep="\t")
        cell_type = df["cell_type"].iloc[0]
        pavlab_de_results[cell_type] = df

    author_de_results = {}
    for path in args.author_paths:
        df = pd.read_csv(path, sep="\t")
        cell_type = df["cell_type"].iloc[0]
        author_de_results[cell_type] = df

    # --- All-vs-all correlations (existing behaviour) ---
    corr_list = []
    for ct_gemma in pavlab_de_results:
        for ct_author in author_de_results:
            merged = pd.merge(
                author_de_results[ct_author], pavlab_de_results[ct_gemma],
                on="gene", suffixes=("_author", "_gemma"), how="inner"
            )
            for col in ["log2FoldChange_author", "log2FoldChange_gemma", "pvalue_author", "pvalue_gemma"]:
                merged[col] = merged[col].astype(float)
            corr_logfc = merged["log2FoldChange_author"].corr(merged["log2FoldChange_gemma"], method="spearman")
            corr_pvalue = merged["pvalue_author"].corr(merged["pvalue_gemma"], method="spearman")
            corr_list.append({
                "ct_sc_pipeline": ct_gemma,
                "ct_author": ct_author,
                "spearman_log2FoldChange": corr_logfc,
                "spearman_pvalue": corr_pvalue,
            })

    all_corrs = pd.DataFrame(corr_list)
    all_corrs.to_csv(f"{args.contrast}_pairwise_corr.tsv", sep="\t", index=False)

    all_corrs_logfc = all_corrs.pivot(index="ct_sc_pipeline", columns="ct_author", values="spearman_log2FoldChange")
    all_corrs_pvalue = all_corrs.pivot(index="ct_sc_pipeline", columns="ct_author", values="spearman_pvalue")

    desired_col_order = [
        "Lamp5", "Lamp5_Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Chandelier", "Vip",
        "L2_3_IT", "L4_IT", "L5_IT", "L6_IT", "L6_IT_Car3", "L5_ET", "L6b", "L5_6_NP", "L6_CT",
        "Oligo", "OPC", "Astro", "Micro", "Immune", "VLMC", "Endo", "PC", "SMC",
    ]
    print(all_corrs_logfc.index)

    if "L2.3_IT" in all_corrs_logfc.index:
        desired_row_order = [
            "Lamp5", "Lamp5_Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Chandelier", "Vip",
            "L2.3_IT", "L4_IT", "L5_ET", "L5_IT", "L6_IT", "L6_IT_Car3", "L6b", "L5.6_NP", "L6_CT",
            "Oligo", "OPC", "Astro", "Micro", "Immune", "Endo", "PC", "SMC", "VLMC"
        ]
    elif "vascular" in all_corrs_logfc.index:
        desired_row_order = [
            "lamp5.GABAergic.cortical.interneuron", "pvalb.GABAergic.cortical.interneuron",
            "sncg.GABAergic.cortical.interneuron", "sst.GABAergic.cortical.interneuron",
            "chandelier.pvalb.GABAergic.cortical.interneuron",
            "vip.GABAergic.cortical.interneuron", "L2.3.6.intratelencephalic.projecting.glutamatergicneuron",
            "deep.layer.non.IT", "oligodendrocyte",
            "oligodendrocyte.precursor.cell", "astrocyte",
            "microglial.cell", "vascular"
        ]
    elif "corticothalamic.projecting.glutamatergic.cortical.neuron" in all_corrs_logfc.index:
        desired_row_order = [
            "lamp5.GABAergic.cortical.interneuron", "pvalb.GABAergic.cortical.interneuron",
            "sncg.GABAergic.cortical.interneuron", "sst.GABAergic.cortical.interneuron",
            "chandelier.pvalb.GABAergic.cortical.interneuron",
            "vip.GABAergic.cortical.interneuron",
            "L2.3.6.intratelencephalic.projecting.glutamatergicneuron",
            "L5.extratelencephalic.projecting.glutamatergic.cortical.neuron",
            "L6b.glutamatergic.cortical.neuron",
            "near.projecting.glutamatergic.cortical.neuron",
            "corticothalamic.projecting.glutamatergic.cortical.neuron",
            "oligodendrocyte", "oligodendrocyte.precursor.cell", "astrocyte", "microglial.cell",
            "endothelial.cell", "pericyte", "smooth.muscle.cell"
        ]
    else:
        desired_row_order = list(all_corrs_logfc.index)

    ordered_logfc = all_corrs_logfc.reindex(index=desired_row_order, columns=desired_col_order)
    ordered_logfc = ordered_logfc.dropna(axis=1, how="all").dropna(axis=0, how="all")
    plt.figure(figsize=(30, 15))
    sns.heatmap(ordered_logfc, cmap="Reds", center=0.5, vmin=0, vmax=1, annot=False, linewidths=0)
    plt.title(f"Spearman Correlation of log2FoldChange: {args.contrast}", fontsize=20)
    plt.xlabel("Author Cell Types", fontsize=20)
    plt.ylabel("Pavlab Cell Types", fontsize=20)
    plt.xticks(rotation=90, fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f"{args.contrast}_pairwise_corr_log2FoldChange.png")
    plt.close()

    ordered_pvalue = all_corrs_pvalue.reindex(index=desired_row_order, columns=desired_col_order)
    ordered_pvalue = ordered_pvalue.dropna(axis=1, how="all").dropna(axis=0, how="all")
    plt.figure(figsize=(30, 15))
    sns.heatmap(ordered_pvalue, cmap="Blues", center=0.5, vmin=0, vmax=1, annot=False, linewidths=0)
    plt.title(f"Spearman Correlation of pvalue: {args.contrast}", fontsize=20)
    plt.xlabel("Author Cell Types", fontsize=20)
    plt.ylabel("Pavlab Cell Types", fontsize=20)
    plt.xticks(rotation=90, fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f"{args.contrast}_pairwise_corr_pvalue.png")
    plt.close()

    # --- Scatter plots for matched pairs only ---
    for author_ct, author_df in author_de_results.items():
        gemma_ct = gemma_to_gemma_map.get(author_ct)
        if gemma_ct is None or gemma_ct not in pavlab_de_results:
            continue
        pavlab_df = pavlab_de_results[gemma_ct]

        merged = pd.merge(author_df, pavlab_df, on="gene", suffixes=("_author", "_gemma"), how="inner")
        if merged.empty:
            continue

        for col in ["log2FoldChange_author", "log2FoldChange_gemma"]:
            merged[col] = merged[col].astype(float)
        for col in ["padj_author", "padj_gemma"]:
            if col not in merged.columns:
                merged[col] = 1.0
            else:
                merged[col] = merged[col].fillna(1.0).astype(float)

        corr_logfc = merged["log2FoldChange_author"].corr(merged["log2FoldChange_gemma"], method="spearman")
        merged["Significance"] = merged.apply(label_significance, axis=1)
        scatter_plot(merged, args.contrast, author_ct, gemma_ct, corr_logfc)


if __name__ == "__main__":
    main()
