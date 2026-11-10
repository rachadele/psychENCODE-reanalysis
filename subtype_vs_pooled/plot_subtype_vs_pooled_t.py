#!/usr/bin/env python3
"""Plot DE agreement (Spearman rho of limma moderated t) between author subtypes and pooled types.

Reads pairwise_rho.tsv, summary.tsv and cells_per_donor.tsv (written by subtype_vs_pooled.py) and
writes subtype_vs_pooled_t.png plus manifest_plot_t.json. Panel A: per-contrast rho (mean over
pairs within a contrast) for each comparison and group, with the mean over the eight contrasts.
Panel A shows author subtype vs pooled author type, vs subtype, vs the nearest separate types
(NEAREST) and vs glia. One panel; cells per donor stay in cells_per_donor.tsv and summary.tsv.
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

# Panel A categories. "nearest" is the closest separate group of author types: deep-layer non-IT
# types for IT subtypes, GABAergic types outside the pooled group for Lamp5 and Sst subtypes.
NEAREST = {"IT": "vs_deep_non_it", "Lamp5": "vs_gaba", "Sst": "vs_gaba"}
COMPARISONS = [
    ("subtype_vs_pooled", "Author subtype\nvs pooled"),
    ("subtype_vs_subtype", "Author subtype\nvs subtype"),
    ("nearest", "vs nearest\nseparate types"),
    ("vs_glia", "vs glia"),
]
GROUPS = ["IT", "Lamp5", "Sst"]
GROUP_CLASS = {"IT": "L2/3-6 IT", "Lamp5": "LAMP5", "Sst": "SST"}
GROUP_LABEL = {"IT": "IT (5 subtypes); nearest: deep-layer non-IT",
               "Lamp5": "Lamp5 (2); nearest: other GABAergic",
               "Sst": "Sst (2); nearest: other GABAergic"}
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
    # map each group's nearest separate types onto one "nearest" category
    for g, c in NEAREST.items():
        means[(g, "nearest")] = means[(g, c)]
    per_c = pd.concat([per_c] + [per_c[(per_c.group == g) & (per_c.comparison == c)].assign(comparison="nearest")
                                 for g, c in NEAREST.items()], ignore_index=True)
    # one category order for all groups: highest to lowest mean of the group means
    order_mean = pd.Series({c: np.mean([means[(g, c)] for g in GROUPS]) for c, _ in COMPARISONS})
    comps = sorted(COMPARISONS, key=lambda c: -order_mean[c[0]])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 1.2})
    fig, axA = plt.subplots(figsize=(9, 6.6))

    rng = np.random.default_rng(0)
    dodge = {"IT": -0.25, "Lamp5": 0.0, "Sst": 0.25}
    for i, (comp, _) in enumerate(comps):
        for g in GROUPS:
            v = per_c[(per_c.group == g) & (per_c.comparison == comp)]["rho_t"].to_numpy()
            x = i + dodge[g]
            axA.scatter(x + rng.uniform(-0.06, 0.06, len(v)), v, s=36, color=col[g],
                        marker=GROUP_MARKER[g], alpha=0.6, edgecolor="0.3", linewidth=0.4, zorder=2)
            axA.plot([x - 0.11, x + 0.11], [means[(g, comp)]] * 2, color="black", lw=2.5, zorder=3)
    axA.set_xticks(range(len(comps)))
    axA.set_xticklabels([lab for _, lab in comps], fontsize=11)
    axA.set_xlim(-0.5, len(comps) - 0.5)
    ylo = min(-0.05, per_c[per_c.comparison.isin([c for c, _ in comps])].rho_t.min() - 0.04)
    axA.set_ylim(ylo, 1.03)
    axA.set_ylabel("Spearman ρ of t statistics", fontsize=13)
    axA.tick_params(axis="y", labelsize=11)
    for i in range(1, len(comps)):
        axA.axvline(i - 0.5, color="0.88", lw=0.8, zorder=0)

    handles = [plt.Line2D([], [], ls="", marker=GROUP_MARKER[g], color=col[g], markeredgecolor="0.3",
                          markersize=8, label=GROUP_LABEL[g]) for g in GROUPS]
    handles.append(plt.Line2D([], [], color="black", lw=2.5, label="Mean of 8 contrasts"))
    axA.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.01), fontsize=10.5, frameon=False,
               ncol=2, handletextpad=0.3, columnspacing=1.0)

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
        "category_order": [c for c, _ in comps],
        "category_order_mean_rho_t": {c: round(float(order_mean[c]), 4) for c, _ in comps},
        "nearest_separate_types": NEAREST,
        "group_mean_rho_t": {f"{g}|{c}": round(float(means[(g, c)]), 4) for g in GROUPS for c, _ in comps},
    }
    json.dump(man, open(OUT_MAN, "w"), indent=2)
    print(json.dumps({k: man[k] for k in ["category_order_mean_rho_t", "group_mean_rho_t"]}, indent=1))


if __name__ == "__main__":
    main()
