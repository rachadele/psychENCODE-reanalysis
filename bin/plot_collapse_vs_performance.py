#!/usr/bin/env python
"""
Show which author cell types collapse into one GEMMA class, and how that
collapse affects log2FC concordance at class vs subclass level.

Two bars per author cell type = mean Spearman rho of log2FC (GEMMA vs author)
at class level (light) and subclass level (dark). Cell types are grouped by
their class-level GEMMA class. Most groups map identically at both levels, so
the two bars are equal; the deep-layer non-IT and vascular collapses are split
into per-subtype GEMMA classes at subclass level, so their subclass bar extends
past the class bar (PC/SMC have no subclass GEMMA class and show class only).
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# class-level GEMMA term -> author cell types (split = collapse undone at subclass)
GROUPS = [
    ("IT", ["L2.3_IT", "L4_IT", "L5_IT", "L6_IT", "L6_IT_Car3"]),
    ("deep-layer non-IT", ["L5.6_NP", "L6b", "L6_CT", "L5_ET"]),
    ("vascular", ["Endo", "VLMC", "PC", "SMC"]),
    ("SST", ["Sst", "Sst_Chodl"]),
    ("LAMP5", ["Lamp5", "Lamp5_Lhx6"]),
    ("microglia", ["Micro", "Immune"]),
    ("1:1", ["OPC", "Vip", "Astro", "Pvalb", "Chandelier", "Oligo", "Sncg", "Pax6"]),
]

CLASS_C, SUB_C = "#a9cce3", "#21618c"  # same hue, light = class, dark = subclass


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--comparison_dir", default="results/comparison")
    p.add_argument("--outdir", default="results/comparison/figs")
    return p.parse_args()


def per_ct_mean(level_dir):
    path = os.path.join(level_dir, "average_corr", "heatmap_data",
                        "heatmap_spearman_log2FoldChange.tsv")
    df = pd.read_csv(path, sep="\t", index_col=0)  # rows=contrast, cols=cell type
    return df.mean(axis=0)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cls = per_ct_mean(os.path.join(args.comparison_dir, "class"))
    sub = per_ct_mean(os.path.join(args.comparison_dir, "subclass"))

    rows, y, group_spans = [], 0.0, []
    for label, cts in GROUPS:
        members = [(ct, cls.get(ct, np.nan), sub.get(ct, np.nan)) for ct in cts]
        members = [m for m in members if not np.isnan(m[1])]
        members.sort(key=lambda t: t[1])  # worst first
        y_bot = y
        for ct, c, s in members:
            rows.append((y, ct, c, s))
            y += 1
        group_spans.append((label, (y_bot + y - 1) / 2))
        y += 1.4

    fig, ax = plt.subplots(figsize=(9, 11))
    off, h = 0.2, 0.36
    for yi, ct, c, s in rows:
        changed = (not np.isnan(s)) and abs(s - c) >= 0.02
        ax.barh(yi - off, c, height=h, color=CLASS_C, edgecolor="#333", linewidth=0.4)
        if not np.isnan(s):
            ax.barh(yi + off, s, height=h, color=SUB_C, edgecolor="#333", linewidth=0.4)
        if changed:
            ax.text(c + 0.008, yi - off, f"{c:.2f}", va="center", fontsize=7, color="#555")
            ax.text(s + 0.008, yi + off, f"{s:.2f}", va="center", fontsize=7.5,
                    color=SUB_C, fontweight="bold")
        else:
            ax.text(c + 0.008, yi, f"{c:.2f}", va="center", fontsize=7.5, color="#444")

    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([r[1] for r in rows], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.08)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Mean Spearman ρ of log2FC  (GEMMA vs author labels)", fontsize=11)

    for label, yc in group_spans:
        ax.text(1.10, yc, label, va="center", ha="left", fontsize=10, clip_on=False)

    ax.grid(axis="x", ls="-", alpha=0.15)
    ax.set_title("Cell-type collapse vs log2FC concordance", fontsize=13,
                 fontweight="bold")

    from matplotlib.patches import Patch
    legend = [Patch(facecolor=CLASS_C, edgecolor="#333", label="class ρ"),
              Patch(facecolor=SUB_C, edgecolor="#333", label="subclass ρ")]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.06),
              ncol=2, fontsize=10, frameon=False)

    fig.subplots_adjust(left=0.14, right=0.84, top=0.94, bottom=0.09)
    out = os.path.join(args.outdir, "collapse_vs_log2fc.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[wrote] {out}")


if __name__ == "__main__":
    main()
