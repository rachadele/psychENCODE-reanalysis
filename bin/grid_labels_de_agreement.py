#!/usr/bin/env python3
"""DE agreement with author-label DE against annotation F1, over several label arms.

For each arm (a label set whose DE was run on the same cohorts as the author-label reference):
  * agreement: Spearman rho of log2 fold change and of the moderated t between the arm's DE and
    the author-label DE, per author cell type x contrast (bin/compare_label_arms.py output, the
    comparison used for results_v2.1.2_limma); contrasts SCZ, BD, ASD, MDD, PTSD, WS, Age, Sex
    when present;
  * F1: per-sample F1 of the arm's labels against author labels (bin/per_cohort_f1.py), averaged
    over samples per cohort and cell type; a cell type x contrast pair takes the unweighted mean
    over the cohorts that contribute the contrast's cases (all cohorts for Age and Sex);
  * support: median author cells per donor per cohort (zeros included), averaged the same way,
    log10(x + 1).
Then the mean rho per arm and, over all pairs of all arms (and per arm), plain Spearman of rho
with F1 and partial Spearman controlling for log10 cells per donor. Cells are never pooled
across samples and nothing is weighted.

--arm takes name=compare_tsv=per_cohort_f1_dir (repeatable).
Writes pairs.tsv, per_cohort_f1.tsv, arm_summary.tsv, rho_vs_f1.tsv and manifest.json.
"""
import argparse
import datetime
import glob
import hashlib
import json
import os
import re
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

CONTRAST_DISPLAY = {
    "Disorder_Schizophrenia_vs_Control": "SCZ", "Disorder_Bipolar.Disorder_vs_Control": "BD",
    "Disorder_ASD_vs_Control": "ASD", "Disorder_MDD_vs_Control": "MDD",
    "Disorder_PTSD_vs_Control": "PTSD", "Disorder_Williams.Syndrome_vs_Control": "WS",
    "Age_death": "Age", "Biological_Sex_male_vs_female": "Sex"}
METRICS = {"log2FC": "spearman_log2FoldChange", "t": "spearman_stat"}


def rname(s):
    return re.sub(r"[^A-Za-z0-9_.]", ".", str(s))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~np.isnan(x) & ~np.isnan(y)
    r, p = stats.spearmanr(x[ok], y[ok])
    return r, p, int(ok.sum())


def partial_spearman(x, y, z):
    """Spearman of x and y controlling for z: Pearson of rank residuals, t test on n-3 df."""
    x, y, z = (np.asarray(v, float) for v in (x, y, z))
    ok = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z)
    rx, ry, rz = (stats.rankdata(v[ok]) for v in (x, y, z))
    res = lambda a: a - np.polyval(np.polyfit(rz, a, 1), rz)
    r = np.corrcoef(res(rx), res(ry))[0, 1]
    n = int(ok.sum())
    t = r * np.sqrt((n - 3) / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), n - 3), n


def contrast_cohorts(meta_dir, studies, contrasts):
    meta = {s: pd.read_csv(os.path.join(meta_dir, f"{s}_sample_meta.tsv"), sep="\t", dtype=str) for s in studies}
    out = {}
    for c in contrasts:
        m = re.match(r"Disorder_(.+)_vs_Control$", c)
        out[c] = (sorted(k for k, d in meta.items() if (d["Disorder"].fillna("").map(rname) == m.group(1)).any())
                  if m else sorted(meta))
        if not out[c]:
            raise SystemExit(f"no cohort for contrast {c}")
    return out


def cohort_support(cell_counts, studies):
    cc = pd.read_csv(cell_counts, sep="\t", dtype={"sample": str})
    rows = []
    for st, g in cc[cc["study"].isin(studies)].groupby("study"):
        wide = g.pivot_table(index="sample", columns="ct", values="n_kept", aggfunc="sum").fillna(0)
        rows.append(pd.DataFrame({"study": st, "ct_author": wide.columns,
                                  "median_cells_per_donor": wide.median().values, "n_donors": len(wide)}))
    return pd.concat(rows, ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--arm", action="append", required=True, help="name=compare_tsv=per_cohort_f1_dir")
    p.add_argument("--author_cell_counts", required=True, help="author-label ct_pseudobulks/cell_counts.tsv")
    p.add_argument("--meta_dir", required=True)
    p.add_argument("--studies", nargs="+", required=True)
    p.add_argument("--data_version", default="")
    p.add_argument("--extra_inputs", nargs="*", default=[], help="files whose sha256 go in the manifest")
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    sup = cohort_support(a.author_cell_counts, a.studies)
    inputs = [a.author_cell_counts] + [os.path.join(a.meta_dir, f"{s}_sample_meta.tsv") for s in a.studies]
    pairs, f1_all = [], []
    for spec in a.arm:
        name, cmp_path, f1_dir = spec.split("=")
        f1_files = [os.path.join(f1_dir, f"{s}_per_cohort_f1.tsv") for s in a.studies]
        inputs += [cmp_path] + f1_files + [os.path.join(f1_dir, f"{s}_per_sample_f1.tsv") for s in a.studies]
        f1 = pd.concat([pd.read_csv(f, sep="\t") for f in f1_files], ignore_index=True)
        f1_all.append(f1.assign(arm=name))
        coh = f1.merge(sup, on=["study", "ct_author"], how="left").set_index(["study", "ct_author"])
        cmp = pd.read_csv(cmp_path, sep="\t")
        cmp = cmp[cmp["contrast"].isin(CONTRAST_DISPLAY) & (cmp["status"] == "ok")]
        cmap = contrast_cohorts(a.meta_dir, a.studies, sorted(cmp["contrast"].unique()))
        for r in cmp.itertuples():
            keys = [(s, r.ct_author) for s in cmap[r.contrast] if (s, r.ct_author) in coh.index]
            v = coh.loc[keys] if keys else None
            pairs.append({"arm": name, "contrast": r.contrast, "contrast_display": CONTRAST_DISPLAY[r.contrast],
                          "ct_author": r.ct_author, "ct_label": r.ct_sc_pipeline, "n_genes": r.n_genes,
                          "rho_log2FC": r.spearman_log2FoldChange, "rho_t": r.spearman_stat,
                          "cohorts": ";".join(cmap[r.contrast]), "n_cohorts_with_type": len(keys),
                          "f1": v["f1"].mean() if keys else np.nan,
                          "median_cells_per_donor": v["median_cells_per_donor"].mean() if keys else np.nan})
    pairs = pd.DataFrame(pairs)
    pairs["log10_cells"] = np.log10(pairs["median_cells_per_donor"] + 1)
    f1_all = pd.concat(f1_all, ignore_index=True)
    pairs.to_csv(os.path.join(a.outdir, "pairs.tsv"), sep="\t", index=False, float_format="%.6g")
    f1_all.to_csv(os.path.join(a.outdir, "per_cohort_f1.tsv"), sep="\t", index=False, float_format="%.6g")

    arms = [s.split("=")[0] for s in a.arm]
    summ = []
    for name in arms:
        d, f = pairs[pairs["arm"] == name], f1_all[f1_all["arm"] == name]
        summ.append({"arm": name, "n_pairs": len(d), "n_types": d["ct_author"].nunique(),
                     "mean_rho_log2FC": d["rho_log2FC"].mean(), "mean_rho_t": d["rho_t"].mean(),
                     "median_rho_log2FC": d["rho_log2FC"].median(), "median_rho_t": d["rho_t"].median(),
                     "mean_f1_pairs": d["f1"].mean(), "min_f1_pairs": d["f1"].min(), "max_f1_pairs": d["f1"].max(),
                     "macro_f1_mean_over_cohorts": f.groupby("study")["f1"].mean().mean(),
                     "n_types_without_de": int(f["ct_author"].nunique() - d["ct_author"].nunique())})
    summ = pd.DataFrame(summ)
    summ.to_csv(os.path.join(a.outdir, "arm_summary.tsv"), sep="\t", index=False, float_format="%.6g")

    corr = []
    grid = [n for n in arms if n.startswith("grid_")]
    sets = [("all arms", arms)] + ([("grid arms", grid)] if len(grid) > 1 else []) + [(n, [n]) for n in arms]
    for set_name, members in sets:
        d = pairs[pairs["arm"].isin(members)]
        for mname in METRICS:
            y = d[f"rho_{mname}"]
            for kind, fn, args in (("plain", spearman, (d["f1"], y)),
                                   ("partial | log10 cells per donor", partial_spearman, (d["f1"], y, d["log10_cells"]))):
                r_, p_, n_ = fn(*args)
                corr.append({"points": set_name, "agreement": mname, "test": f"rho ~ F1 ({kind})",
                             "spearman": r_, "p": p_, "n_pairs": n_})
    corr = pd.DataFrame(corr)
    corr.to_csv(os.path.join(a.outdir, "rho_vs_f1.tsv"), sep="\t", index=False, float_format="%.4g")

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    git = lambda *x: subprocess.check_output(["git", "-C", repo, *x], text=True).strip()
    outs = [os.path.join(a.outdir, f) for f in ("pairs.tsv", "per_cohort_f1.tsv", "arm_summary.tsv", "rho_vs_f1.tsv")]
    man = {"created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
           "data_version": a.data_version, "script": os.path.abspath(__file__), "script_sha256": sha256(__file__),
           "repo_commit": git("rev-parse", "HEAD"), "repo_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
           "args": vars(a),
           "input_sha256": {f: {"sha256": sha256(f), "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(f)).isoformat(timespec="seconds")}
                            for f in inputs + a.extra_inputs},
           "output_sha256": {f: sha256(f) for f in outs}}
    json.dump(man, open(os.path.join(a.outdir, "manifest.json"), "w"), indent=1)
    print(summ.to_string(index=False))
    print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
