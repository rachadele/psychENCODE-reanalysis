#!/usr/bin/env python3
"""Plot DE agreement (Spearman rho of limma moderated t) between author subtypes and pooled types.

Reads pairwise_rho.tsv, summary.tsv and cells_per_donor.tsv (written by subtype_vs_pooled.py) and
writes subtype_vs_pooled_t.png plus manifest_plot_t.json. Panel A: per-contrast rho (mean over
pairs within a contrast) for each comparison and group, with the mean over the eight contrasts.
Panel B: each subtype's mean rho with its pooled author type against median cells per donor.
Colours come from the shared manuscript palette (class level).
"""
import csv
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PALETTE_TSV = os.environ.get("MANUSCRIPT_PALETTE",
                             "/space/grp/rschwartz/rschwartz/benchmark-post-hoc/manuscript_palette.tsv")
IN_PAIR = os.path.join(HERE, "pairwise_rho.tsv")
IN_SUM = os.path.join(HERE, "summary.tsv")
IN_CELLS = os.path.join(HERE, "cells_per_donor.tsv")
OUT_PNG = os.path.join(HERE, "subtype_vs_pooled_t.png")
OUT_MAN = os.path.join(HERE, "manifest_plot_t.json")

COMPARISONS = [
    ("subtype_vs_subtype", "Author subtype\nvs subtype"),
    ("subtype_vs_pooled", "Author subtype\nvs pooled"),
    ("baseline_same_class", "Unrelated types,\nsame class"),
    ("baseline_all_other", "Unrelated types,\nall"),
    ("pooled_vs_ours", "Pooled author\nvs our label"),
]
GROUPS = ["IT", "Lamp5", "Sst"]
GROUP_CLASS = {"IT": "L2/3-6 IT", "Lamp5": "LAMP5", "Sst": "SST"}
GROUP_LABEL = {"IT": "IT (5 subtypes)", "Lamp5": "Lamp5 (2)", "Sst": "Sst (2)"}
GROUP_MARKER = {"IT": "o", "Lamp5": "s", "Sst": "D"}
SUBTYPE_LABEL = {"L2.3_IT": "L2/3 IT", "L4_IT": "L4 IT", "L5_IT": "L5 IT", "L6_IT": "L6 IT",
                 "L6_IT_Car3": "L6 IT Car3", "Lamp5": "Lamp5", "Lamp5_Lhx6": "Lamp5 Lhx6",
                 "Sst": "Sst", "Sst_Chodl": "Sst Chodl"}


def colours():
    pal = {r["label"]: r["hex"] for r in csv.DictReader(open(PALETTE_TSV), delimiter="\t")
           if r["level"] == "class"}
    return {g: pal[GROUP_CLASS[g]] for g in GROUPS}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stamp(p):
    return datetime.fromtimestamp(os.path.getmtime(p)).astimezone().isoformat(timespec="seconds")


def main():
    col = colours()
    pair = pd.read_csv(IN_PAIR, sep="\t")
    summ = pd.read_csv(IN_SUM, sep="\t")
    cells = pd.read_csv(IN_CELLS, sep="\t")

    # per-contrast value = mean over pairs within contrast (same as subtype_vs_pooled.py)
    per_c = pair.groupby(["group", "comparison", "contrast"])["rho_t"].mean().reset_index()
    means = per_c.groupby(["group", "comparison"])["rho_t"].mean()
    # check against summary.tsv
    sg = summ[summ.level == "group"].set_index(["group", "comparison"])["mean_rho_t"]
    for k, v in means.items():
        assert abs(v - sg[k]) < 1e-3, (k, v, sg[k])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 1.2})
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw={"width_ratios": [1.75, 1]})

    rng = np.random.default_rng(0)
    dodge = {"IT": -0.25, "Lamp5": 0.0, "Sst": 0.25}
    for i, (comp, _) in enumerate(COMPARISONS):
        for g in GROUPS:
            v = per_c[(per_c.group == g) & (per_c.comparison == comp)]["rho_t"].to_numpy()
            x = i + dodge[g]
            axA.scatter(x + rng.uniform(-0.06, 0.06, len(v)), v, s=30, color=col[g],
                        marker=GROUP_MARKER[g], alpha=0.55, edgecolor="0.3", linewidth=0.4, zorder=2)
            axA.plot([x - 0.1, x + 0.1], [means[(g, comp)]] * 2, color="black", lw=2.5, zorder=3)
    axA.set_xticks(range(len(COMPARISONS)))
    axA.set_xticklabels([lab for _, lab in COMPARISONS], fontsize=11)
    axA.set_xlim(-0.5, len(COMPARISONS) - 0.5)
    axA.set_ylim(-0.05, 1.03)
    axA.set_ylabel("Spearman ρ of t statistics", fontsize=13)
    axA.tick_params(axis="y", labelsize=11)
    for i in range(1, len(COMPARISONS)):
        axA.axvline(i - 0.5, color="0.88", lw=0.8, zorder=0)
    handles = [plt.Line2D([], [], ls="", marker=GROUP_MARKER[g], color=col[g], markeredgecolor="0.3",
                          markersize=8, label=GROUP_LABEL[g]) for g in GROUPS]
    handles.append(plt.Line2D([], [], color="black", lw=2.5, label="Mean of 8 contrasts"))
    axA.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.01), fontsize=11, frameon=False,
               ncol=4, handletextpad=0.3, columnspacing=1.2)
    axA.text(-0.09, 1.06, "A", transform=axA.transAxes, fontsize=16, fontweight="bold")

    # panel B
    st = summ[(summ.level == "type") & (summ.comparison == "subtype_vs_pooled")].copy()
    cpd = cells[cells.scheme == "author"].set_index("celltype")["median_cells_per_donor"]
    st["cells"] = st["type_a"].map(cpd)
    # label positions in data coordinates (None = offset next to the point), chosen by eye on the PNG
    label_at = {"L2.3_IT": ((1082, 1.0), "center"), "L5_IT": ((750, 0.80), "left"),
                "L6_IT": ((750, 0.72), "left"), "L4_IT": ((750, 0.64), "left"),
                "L6_IT_Car3": (None, (9, 0, "left")), "Lamp5": (None, (-10, 0, "right")),
                "Lamp5_Lhx6": (None, (-10, 0, "right")), "Sst": (None, (-10, 0, "right")),
                "Sst_Chodl": (None, (9, 0, "left"))}
    for _, r in st.iterrows():
        g = r["group"]
        axB.errorbar(r["cells"], r["mean_rho_t"], yerr=r["sd_rho_t"], fmt=GROUP_MARKER[g], ms=9,
                     color=col[g], mec="0.3", mew=0.6, ecolor=col[g], elinewidth=1.4, capsize=3, zorder=3)
        pos, how = label_at[r["type_a"]]
        xy = (r["cells"], r["mean_rho_t"])
        if pos is None:
            dx, dy, ha = how
            axB.annotate(SUBTYPE_LABEL[r["type_a"]], xy, xytext=(dx, dy), textcoords="offset points",
                         ha=ha, va="center", fontsize=11)
        else:
            same = abs(pos[0] - xy[0]) < 1
            axB.annotate(SUBTYPE_LABEL[r["type_a"]], xy, xytext=pos, ha=how, va="center", fontsize=11,
                         arrowprops=None if same else dict(arrowstyle="-", color="0.5", lw=0.8,
                                                           shrinkA=2, shrinkB=5))
    axB.set_xscale("log")
    axB.set_xlim(2, 4000)
    axB.set_ylim(-0.05, 1.05)
    axB.set_xlabel("Median cells per donor (subtype)", fontsize=13)
    axB.set_ylabel("ρ of t, subtype vs pooled", fontsize=13)
    axB.tick_params(labelsize=11)
    axB.text(-0.2, 1.06, "B", transform=axB.transAxes, fontsize=16, fontweight="bold")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300, facecolor="white")

    try:
        commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = None
    up = json.load(open(os.path.join(HERE, "manifest.json")))
    man = {
        "created": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "data_version": up.get("data_version"),
        "upstream_manifest": {"path": "subtype_vs_pooled/manifest.json", "created": up.get("created"),
                              "sha256": sha(os.path.join(HERE, "manifest.json"))},
        "repo_commit_at_run": commit,
        "script": {"subtype_vs_pooled/plot_subtype_vs_pooled_t.py": {"sha256": sha(os.path.abspath(__file__))}},
        "inputs": {os.path.relpath(p, os.path.dirname(HERE)) if p.startswith(os.path.dirname(HERE)) else p:
                   {"sha256": sha(p), "mtime": stamp(p)} for p in [IN_PAIR, IN_SUM, IN_CELLS, PALETTE_TSV]},
        "outputs": {"subtype_vs_pooled/subtype_vs_pooled_t.png": {"sha256": sha(OUT_PNG), "mtime": stamp(OUT_PNG)}},
        "group_mean_rho_t": {f"{g}|{c}": round(float(means[(g, c)]), 4) for g in GROUPS for c, _ in COMPARISONS},
    }
    json.dump(man, open(OUT_MAN, "w"), indent=2)
    print(json.dumps(man["group_mean_rho_t"], indent=1))


if __name__ == "__main__":
    main()
