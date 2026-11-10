"""Figure: split-half reproducibility of DE (A) and standard error relative to pooled (B) per type.

Reads rho_per_contrast.tsv, summary.tsv, se_ratio_per_contrast.tsv (this folder) and writes
split_half_t.png, then adds the figure and its inputs to manifest.json (section "figure").
A: Spearman rho of moderated t between donor halves, one dot per contrast (mean over 10 donor
splits), black bar = mean over the eight contrasts. B: median ratio of standard errors, type over
pooled author type, one dot per contrast, from the full-data fits. Same x axis in both panels.
Within each group, subtypes are ordered by median cells per donor (most first), then the pooled
author type and our label. Subtype colours are shades of the group's class colour in the shared
manuscript palette (darker = more cells per donor); pooled author type and our label are greys.
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
SUB_NAME = {"L2.3_IT": "L2/3 IT", "L4_IT": "L4 IT", "L5_IT": "L5 IT", "L6_IT": "L6 IT",
            "L6_IT_Car3": "L6 IT Car3", "Lamp5": "Lamp5", "Lamp5_Lhx6": "Lamp5 Lhx6",
            "Sst": "Sst", "Sst_Chodl": "Sst Chodl"}
# two-line tick names so cells per donor sits on the third line of every tick label
TICK = {"L2.3_IT": "L2/3\nIT", "L4_IT": "L4\nIT", "L5_IT": "L5\nIT", "L6_IT": "L6\nIT",
        "L6_IT_Car3": "L6 IT\nCar3", "Lamp5": "Lamp5\n", "Lamp5_Lhx6": "Lamp5\nLhx6",
        "Sst": "Sst\n", "Sst_Chodl": "Sst\nChodl", "merged": "Pooled\nauthor", "ours": "Our\nlabel"}
NEUTRAL = {"merged": ("#4D4D4D", "s", "Pooled author type"), "ours": ("#A6A6A6", "D", "Our label")}


def shades(hex_, n):
    """n shades of a colour, darkest first (mixed toward black, then toward white)."""
    from matplotlib.colors import to_rgb
    c = np.array(to_rgb(hex_))
    f = np.linspace(-0.45, 0.55, n) if n > 2 else np.array([-0.3, 0.45])[:n]
    return [tuple(c * (1 + v)) if v < 0 else tuple(c + (1 - c) * v) for v in f]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    pal = {r["label"]: r["hex"] for r in csv.DictReader(open(PALETTE_TSV), delimiter="\t") if r["level"] == "class"}
    rho = pd.read_csv(IN["rho"], sep="\t")
    summ = pd.read_csv(IN["summary"], sep="\t")
    se = pd.read_csv(IN["se"], sep="\t")

    cols, x, pos, sub_handles = [], 0, {}, []
    for g in GROUPS:
        s = summ[summ.group == g]
        subs = s[s.scheme == "author"].sort_values("median_cells_per_donor", ascending=False)
        shade = shades(pal[GROUP_CLASS[g]], len(subs))
        start = x
        for (_, r), c in zip(subs.iterrows(), shade):
            cols.append(dict(group=g, scheme="author", celltype=r.celltype, x=x, color=c, marker="o",
                             label=f"{TICK[r.celltype]}\n{r.median_cells_per_donor:.0f}"))
            sub_handles.append(plt.Line2D([], [], ls="", marker="o", color=c, markersize=10, label=SUB_NAME[r.celltype]))
            x += 1
        for sch in ["merged", "ours"]:
            r = s[s.scheme == sch].iloc[0]
            c, mk, _ = NEUTRAL[sch]
            cols.append(dict(group=g, scheme=sch, celltype=r.celltype, x=x, color=c, marker=mk,
                             label=f"{TICK[sch]}\n{r.median_cells_per_donor:.0f}"))
            x += 1
        pos[g] = (start, x - 1)
        x += 0.8

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 1.2})
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(16, 9.5), sharex=True, gridspec_kw={"height_ratios": [2.4, 1]})
    rng = np.random.default_rng(0)
    for c in cols:
        sel = lambda d: d[(d.scheme == c["scheme"]) & (d.celltype == c["celltype"]) & (d.group == c["group"])]
        for ax, d, ycol in [(axA, sel(rho), "mean_rho_t"), (axB, sel(se), "median_se_ratio_vs_pooled")]:
            v = d[ycol].dropna().to_numpy()
            if len(v) == 0:
                continue
            ax.scatter(c["x"] + rng.uniform(-0.15, 0.15, len(v)), v, s=60, marker=c["marker"], color=c["color"],
                       edgecolor="0.2", linewidth=0.6, zorder=2)
            ax.plot([c["x"] - 0.32, c["x"] + 0.32], [v.mean()] * 2, color="black", lw=2.8, zorder=3)
    for g in GROUPS:
        a, b = pos[g]
        axA.text((a + b) / 2, 1.0, GROUP_TITLE[g], transform=axA.get_xaxis_transform(), ha="center",
                 va="bottom", fontsize=14, fontweight="bold")
        for ax in (axA, axB):
            if g != GROUPS[0]:
                ax.axvline(a - 0.9, color="0.8", lw=1, zorder=0)
    axB.axhline(1, color="0.4", lw=1.2, ls="--", zorder=1)
    for c in cols:  # the pooled author column is empty in B (ratio 1 by definition): label the line there
        if c["scheme"] == "merged":
            axB.text(c["x"], 1.06, "= 1", ha="center", va="bottom", fontsize=11, color="0.3")
    axB.set_xticks([c["x"] for c in cols])
    axB.set_xticklabels([c["label"] for c in cols], fontsize=11)
    axB.set_xlim(-0.6, cols[-1]["x"] + 0.6)
    axA.set_ylabel("Agreement between donor halves\n(Spearman ρ of t)", fontsize=13)
    axB.set_ylabel("Standard error\nrelative to pooled", fontsize=13)
    axA.set_ylim(-0.15, 1.0)
    axA.axhline(0, color="0.7", lw=0.8, zorder=1)
    axA.text(0.995, 0.97, "One dot = one contrast (mean over 10 donor splits)\nBar = mean over 8 contrasts",
             transform=axA.transAxes, ha="right", va="top", fontsize=12)
    for ax in (axA, axB):
        ax.tick_params(axis="y", labelsize=11)
    axA.text(-0.075, 1.04, "A", transform=axA.transAxes, fontsize=16, fontweight="bold")
    axB.text(-0.075, 1.04, "B", transform=axB.transAxes, fontsize=16, fontweight="bold")
    handles = sub_handles + [plt.Line2D([], [], ls="", marker=NEUTRAL[k][1], color=NEUTRAL[k][0], markersize=10,
                                        label=NEUTRAL[k][2]) for k in ("merged", "ours")]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.835, 0.6), fontsize=12, frameon=False,
               title="Cell type", title_fontsize=13, alignment="left")
    fig.tight_layout(rect=(0.02, 0, 0.83, 0.97))
    # row label for the third tick-label line (cells per donor), placed after layout
    fig.canvas.draw()
    bb = axB.get_xticklabels()[0].get_window_extent().transformed(fig.transFigure.inverted())
    fig.text(axB.get_position().x0 - 0.008, bb.y0 + bb.height / 6, "Cells per donor:", fontsize=11,
             color="0.25", ha="right", va="center")
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
