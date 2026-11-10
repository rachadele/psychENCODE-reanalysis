"""Example figure: gene t statistics in donor half A vs half B for one contrast and one split.

Contrast choice: for the IT group, per contrast, pooled author rho minus the mean rho of the five
IT subtypes (rho_per_contrast.tsv); we took the contrast (Williams Syndrome excluded) whose
difference was closest to the mean difference over those contrasts. Seed choice: for that
contrast, the seed minimising the summed squared deviation of the four shown types' rho
(rho_per_split.tsv) from their 10-seed means. Panels: L6 IT Car3, L2/3 IT, pooled author IT and
our L2/3-6 IT label; genes present in both halves; rho as in summarize_split_half.py.
Writes split_half_example_scatter.png and the "figure_example" section of manifest.json.
"""
import csv
import datetime
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CONTRASTS, WORK  # noqa: E402
from plot_split_half_t import GROUP_CLASS, PALETTE_TSV, sha, shades  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PNG = os.path.join(HERE, "split_half_example_scatter.png")
IT = "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron"
SHORT2LONG = {v: k for k, v in CONTRASTS.items()}
LONG_NAME = {"SCZ": "schizophrenia vs control", "BD": "bipolar disorder vs control", "ASD": "ASD vs control",
             "MDD": "MDD vs control", "PTSD": "PTSD vs control", "Age": "age at death", "Sex": "male vs female"}
PANELS = [("author", "L6_IT_Car3", "L6 IT Car3"), ("author", "L2.3_IT", "L2/3 IT"),
          ("merged", IT, "Pooled author type (IT)"), ("ours", IT, "Our label (L2/3-6 IT)")]
NEUTRAL = {"merged": "#4D4D4D", "ours": "#A6A6A6"}


def main():
    rpc = pd.read_csv(os.path.join(HERE, "rho_per_contrast.tsv"), sep="\t")
    rps = pd.read_csv(os.path.join(HERE, "rho_per_split.tsv"), sep="\t")
    summ = pd.read_csv(os.path.join(HERE, "summary.tsv"), sep="\t")
    it = rpc[rpc.group == "IT"]
    diff = (it[it.scheme == "merged"].set_index("contrast")["mean_rho_t"]
            - it[it.scheme == "author"].groupby("contrast")["mean_rho_t"].mean()).drop("WS")
    contrast = (diff - diff.mean()).abs().idxmin()

    shown = rps[(rps.group == "IT") & (rps.contrast == contrast)
                & rps.set_index(["scheme", "celltype"]).index.isin([(s, c) for s, c, _ in PANELS])]
    dev = shown.assign(d=shown.rho_t - shown.groupby(["scheme", "celltype"])["rho_t"].transform("mean"))
    seed = int((dev.d ** 2).groupby(dev.seed).sum().idxmin())

    # subtype colours exactly as in split_half_t.png (shades by cells per donor, darkest first)
    pal = {r["label"]: r["hex"] for r in csv.DictReader(open(PALETTE_TSV), delimiter="\t") if r["level"] == "class"}
    subs = summ[(summ.group == "IT") & (summ.scheme == "author")].sort_values("median_cells_per_donor", ascending=False)
    color = dict(zip(subs.celltype, shades(pal[GROUP_CLASS["IT"]], len(subs))))
    color.update({IT + "|merged": NEUTRAL["merged"], IT + "|ours": NEUTRAL["ours"]})
    cells = summ.set_index(["scheme", "celltype"])["median_cells_per_donor"]

    col = SHORT2LONG[contrast]
    data, rhos = [], {}
    for sch, ct, _ in PANELS:
        a, b = (pd.read_csv(f"{WORK}/fits/{seed}/{h}/{sch}__{ct}.tsv.gz", sep="\t", index_col=0) for h in "AB")
        sh = a.index.intersection(b.index)
        x, y = a.loc[sh, col], b.loc[sh, col]
        ok = x.notna() & y.notna()
        r = spearmanr(x[ok], y[ok]).statistic
        ref = rps[(rps.seed == seed) & (rps.scheme == sch) & (rps.celltype == ct) & (rps.contrast == contrast)].rho_t.iloc[0]
        assert abs(r - ref) < 1e-3, (sch, ct, r, ref)
        data.append((x[ok].to_numpy(), y[ok].to_numpy()))
        rhos[f"{sch}|{ct}"] = round(float(r), 4)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 1.2})
    fig, axes = plt.subplots(1, 4, figsize=(17, 5.2), sharex=True, sharey=True)
    lim = max(np.percentile(np.abs(np.concatenate([np.r_[x, y] for x, y in data])), 99.9), 1) * 1.05
    for ax, (sch, ct, title), (x, y) in zip(axes, PANELS, data):
        c = color[ct] if sch == "author" else color[f"{ct}|{sch}"]
        ax.scatter(x, y, s=5, color=c, alpha=0.35, linewidth=0, rasterized=True)
        ax.axhline(0, color="0.8", lw=0.8, zorder=0)
        ax.axvline(0, color="0.8", lw=0.8, zorder=0)
        ax.plot([-lim, lim], [-lim, lim], color="0.5", lw=1, ls="--", zorder=0)
        ax.set_title(f"{title}\n{cells[(sch, ct)]:.0f} cells per donor", fontsize=13)
        ax.text(0.04, 0.96, f"ρ = {rhos[f'{sch}|{ct}']:.2f}\n{len(x):,} genes", transform=ax.transAxes,
                ha="left", va="top", fontsize=13)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=11)
        ax.set_xlabel("t, donor half 1", fontsize=13)
    axes[0].set_ylabel("t, donor half 2", fontsize=13)
    fig.suptitle(f"Example: {LONG_NAME[contrast]}, donor split {seed} of 10 (one dot = one gene)", fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300, facecolor="white")

    mp = os.path.join(HERE, "manifest.json")
    m = json.load(open(mp))
    ts = lambda p: datetime.datetime.fromtimestamp(os.path.getmtime(p)).astimezone().isoformat()
    ins = [os.path.join(HERE, f) for f in ("rho_per_contrast.tsv", "rho_per_split.tsv", "summary.tsv")] + [PALETTE_TSV] + \
          [f"{WORK}/fits/{seed}/{h}/{s}__{c}.tsv.gz" for s, c, _ in PANELS for h in "AB"]
    m["figure_example"] = {
        "created": datetime.datetime.now().astimezone().isoformat(),
        "contrast": contrast, "seed": seed,
        "contrast_rule": "IT pooled author rho minus mean IT subtype rho per contrast (WS excluded); closest to mean",
        "pooled_minus_subtype_rho": {k: round(float(v), 4) for k, v in diff.items()},
        "seed_rule": "min summed squared deviation of the four shown types' rho from their 10-seed means",
        "rho": rhos,
        "script": {"subtype_vs_pooled/split_half/plot_split_half_example.py": {"sha256": sha(__file__)}},
        "inputs": {p: {"sha256": sha(p), "mtime": ts(p)} for p in ins},
        "output": {"subtype_vs_pooled/split_half/split_half_example_scatter.png": {"sha256": sha(OUT_PNG), "mtime": ts(OUT_PNG)}}}
    json.dump(m, open(mp, "w"), indent=1)
    print(contrast, seed, rhos)


if __name__ == "__main__":
    main()
