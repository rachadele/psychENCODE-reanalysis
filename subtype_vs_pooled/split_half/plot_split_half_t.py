"""Figure: split-half reproducibility of DE (A) and standard error ratio (B) per type.

Reads rho_per_contrast.tsv, summary.tsv, se_ratio_per_contrast.tsv (this folder) and writes
split_half_t.png, then adds the figure and its inputs to manifest.json (section "figure").
A: Spearman rho of moderated t between donor halves, one point per contrast (mean over 10 splits),
black bar = mean over the eight contrasts. B: median ratio of standard errors, type over pooled
author type, one point per contrast, from the full-data fits. Same x axis in both panels.
Colours: class level of the shared manuscript palette.
"""
import csv
import datetime
import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PALETTE_TSV = "/space/grp/rschwartz/rschwartz/benchmark-post-hoc/manuscript_palette.tsv"
IN = {k: os.path.join(HERE, f) for k, f in
      [("rho", "rho_per_contrast.tsv"), ("summary", "summary.tsv"), ("se", "se_ratio_per_contrast.tsv")]}
OUT_PNG = os.path.join(HERE, "split_half_t.png")
GROUPS = ["IT", "Lamp5", "Sst"]
GROUP_CLASS = {"IT": "L2/3-6 IT", "Lamp5": "LAMP5", "Sst": "SST"}
GROUP_TITLE = {"IT": "L2/3-6 IT", "Lamp5": "Lamp5", "Sst": "Sst"}
# two-line names so cells per donor sits on the third line of every tick label
SUB_LABEL = {"L2.3_IT": "L2/3\nIT", "L4_IT": "L4\nIT", "L5_IT": "L5\nIT", "L6_IT": "L6\nIT",
             "L6_IT_Car3": "L6 IT\nCar3", "Lamp5": "Lamp5\n", "Lamp5_Lhx6": "Lamp5\nLhx6",
             "Sst": "Sst\n", "Sst_Chodl": "Sst\nChodl"}
KIND = {"author": ("Author subtype", "o"), "merged": ("Pooled author type", "s"), "ours": ("Our label", "D")}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    pal = {r["label"]: r["hex"] for r in csv.DictReader(open(PALETTE_TSV), delimiter="\t") if r["level"] == "class"}
    col = {g: pal[GROUP_CLASS[g]] for g in GROUPS}
    rho = pd.read_csv(IN["rho"], sep="\t")
    summ = pd.read_csv(IN["summary"], sep="\t")
    se = pd.read_csv(IN["se"], sep="\t")

    # x layout: per group, subtypes then pooled author then our label, with a gap between groups
    cols, x, pos = [], 0, {}
    for g in GROUPS:
        s = summ[summ.group == g]
        for sch in ["author", "merged", "ours"]:
            for _, r in s[s.scheme == sch].iterrows():
                name = SUB_LABEL.get(r.celltype, "Pooled\nauthor" if sch == "merged" else "Our\nlabel")
                cols.append(dict(group=g, scheme=sch, celltype=r.celltype, x=x,
                                 label=f"{name}\n{r.median_cells_per_donor:.0f}"))
                x += 1
        pos[g] = (cols[-1]["x"] - len(s) + 1, cols[-1]["x"])
        x += 0.8

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 1.2})
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(14, 9), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    rng = np.random.default_rng(0)
    for c in cols:
        kind_lab, mk = KIND[c["scheme"]]
        fill = col[c["group"]] if c["scheme"] == "author" else "white"
        sel = lambda d: d[(d.scheme == c["scheme"]) & (d.celltype == c["celltype"]) & (d.group == c["group"])]
        for ax, d, ycol in [(axA, sel(rho), "mean_rho_t"), (axB, sel(se), "median_se_ratio_vs_pooled")]:
            v = d[ycol].dropna().to_numpy()
            if len(v) == 0:
                continue
            ax.scatter(c["x"] + rng.uniform(-0.15, 0.15, len(v)), v, s=55, marker=mk, facecolor=fill,
                       edgecolor=col[c["group"]], linewidth=1.6, alpha=0.85, zorder=2)
            ax.plot([c["x"] - 0.3, c["x"] + 0.3], [v.mean()] * 2, color="black", lw=2.8, zorder=3)
        if c["scheme"] != "author":
            axA.axvspan(c["x"] - 0.5, c["x"] + 0.5, color="0.94", zorder=0)
            axB.axvspan(c["x"] - 0.5, c["x"] + 0.5, color="0.94", zorder=0)
    for g in GROUPS:
        a, b = pos[g]
        axA.text((a + b) / 2, 1.0, GROUP_TITLE[g], transform=axA.get_xaxis_transform(), ha="center",
                 va="bottom", fontsize=14, fontweight="bold", color=col[g])
    axB.axhline(1, color="0.4", lw=1.2, ls="--", zorder=1)
    axB.set_xticks([c["x"] for c in cols])
    axB.set_xticklabels([c["label"] for c in cols], fontsize=11)
    axB.set_xlim(-0.6, cols[-1]["x"] + 0.6)
    axA.set_ylabel("Split-half Spearman ρ of t", fontsize=13)
    axB.set_ylabel("Standard error,\ntype / pooled", fontsize=13)
    axA.set_ylim(-0.2, 1.0)
    axA.axhline(0, color="0.7", lw=0.8, zorder=1)
    for ax in (axA, axB):
        ax.tick_params(axis="y", labelsize=11)
    axA.text(-0.07, 1.04, "A", transform=axA.transAxes, fontsize=16, fontweight="bold")
    axB.text(-0.07, 1.04, "B", transform=axB.transAxes, fontsize=16, fontweight="bold")
    grey = "0.35"
    handles = [plt.Line2D([], [], ls="", marker="o", markerfacecolor=grey, markeredgecolor=grey, markersize=9,
                          label="Author subtype"),
               plt.Line2D([], [], ls="", marker="s", markerfacecolor="white", markeredgecolor=grey,
                          markeredgewidth=1.6, markersize=9, label="Pooled author type"),
               plt.Line2D([], [], ls="", marker="D", markerfacecolor="white", markeredgecolor=grey,
                          markeredgewidth=1.6, markersize=8, label="Our label"),
               plt.Line2D([], [], color="black", lw=2.8, label="Mean over 8 contrasts")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=12, frameon=False)
    fig.tight_layout(rect=(0.02, 0, 1, 0.95))
    # row label for the third tick-label line (cells per donor), placed after layout
    fig.canvas.draw()
    bb = axB.get_xticklabels()[0].get_window_extent().transformed(fig.transFigure.inverted())
    x0 = axB.get_position().x0
    fig.text(x0 - 0.008, bb.y0 + bb.height / 6, "Cells per donor:", fontsize=11, color="0.25", ha="right", va="center")
    fig.savefig(OUT_PNG, dpi=300, facecolor="white")

    mp = os.path.join(HERE, "manifest.json")
    m = json.load(open(mp))
    ts = lambda p: datetime.datetime.fromtimestamp(os.path.getmtime(p)).astimezone().isoformat()
    m["figure"] = {"created": datetime.datetime.now().astimezone().isoformat(),
                   "script": {"subtype_vs_pooled/split_half/plot_split_half_t.py": {"sha256": sha(__file__)}},
                   "inputs": {p: {"sha256": sha(p), "mtime": ts(p)} for p in list(IN.values()) + [PALETTE_TSV]},
                   "output": {"subtype_vs_pooled/split_half/split_half_t.png": {"sha256": sha(OUT_PNG), "mtime": ts(OUT_PNG)}}}
    json.dump(m, open(mp, "w"), indent=1)


if __name__ == "__main__":
    main()
