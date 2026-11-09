#!/usr/bin/env python3
"""DE agreement (author vs pipeline labels) against annotation quality and support.

Per author subclass of one study (default CMC):
  * support: cells per donor in the QC-masked pseudobulks (median over all donors of the
    study, zeros included, and total), for the author-label pseudobulk of the type and for
    the pipeline-label pseudobulk of the label it maps to;
  * CMC F1: F1 of the pipeline labels against the author labels on the cells used in the
    DE (QC-masked cells excluded). An author type t maps to pipeline label L = map[t];
    recall = share of t cells labelled L; precision = share of L cells whose author type maps
    to L (so author types that share L share its precision);
  * mean Spearman rho of log2 fold changes across contrasts (as in the heatmap);
  * benchmark F1 of L (extract_benchmark_f1.py).
Then Spearman correlations across cell types of mean rho with benchmark F1, CMC F1 and
log10 support, and partial Spearman of rho with each F1 controlling for log10 support. This CMC
analysis pools cells across samples; it is kept for comparison as *_cmc_pooled.tsv.

Current analysis (Figure 4; per_cohort_analysis). F1 and cells per donor follow --support_rule:
  de_donors (default; DE_AGREEMENT_DECISIONS.md, 2026-11-09): the limma model is one fit per cell
    type over all cohorts' donors, so every contrast uses the same donors. F1 is the plain mean of
    per-sample F1 (per_cohort_f1.py) over the donors in the DE of that author type (donors with at
    least one unmasked cell of the type and metadata); cells per donor is the median n_kept over the
    same donors (zeros excluded, as limma drops them). Written to celltype_f1_cells.tsv, one row per
    scheme (author, pipeline subclass, pipeline class) and cell type; it is the single F1 /
    cells-per-donor table for Figure 4A-C, the scatter and matched_effect_adjusted.py.
  case_cohorts (superseded): mean over the cohorts that contribute each contrast's cases of the
    per-cohort mean F1 and per-cohort median cells per donor (zeros included). Outputs get a
    _case_cohorts suffix.
--one_to_one: strict (default; 18 types, every type whose pipeline label pools no other author
type) or own_label (20 types; keeps LAMP5 and SST, which share a label with LAMP5 LHX6 and SST Chodl).
Plain and partial Spearman over pairs and per cell type after averaging over contrasts (Figure 4B),
a mixed model on ranks with a cell type random intercept, and rule_alternatives.tsv (each rule and
the CMC pooled values, for all, own-label and strict one-to-one types).
Writes summary TSVs, the combined figure and a provenance JSON to --outdir.
"""
import argparse
import datetime
import hashlib
import json
import os
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from fig_de_agreement import plot_figure
from utils import (author_display, contrast_pivot, load_benchmark_f1, load_label_map,
                   own_label_types)


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


def cmc_f1(h5ad, cta, mask, label_map):
    obs = ad.read_h5ad(h5ad, backed="r").obs
    df = pd.DataFrame({"sample_id": obs["sample_id"].astype(str).values,
                       "row_id": obs.index.astype(str),
                       "cell_id": obs["cell_id"].astype(str).values,
                       "sample_name": obs["sample_name"].astype(str).values})
    # author label from cell_id (<label>[.<n>]_<sample_name>), as build_local_pseudobulks.py
    df["author"] = [rname(re.sub(r"\.\d+$", "", c[:-len(s) - 1])) for c, s in zip(df["cell_id"], df["sample_name"])]
    c = pd.read_csv(cta, sep="\t", dtype=str)
    c["pipeline"] = c["cell_type"].map(rname)
    df = df.merge(c[["sample_id", "cell_id", "pipeline"]].rename(columns={"cell_id": "row_id"}),
                  on=["sample_id", "row_id"], how="left")
    m = pd.read_csv(mask, sep="\t", dtype=str)
    m = m[m["value"].str.lower() == "true"]
    masked = set(zip(m["sample_id"], m["cell_id"]))
    df["masked"] = [(s, r) in masked for s, r in zip(df["sample_id"], df["row_id"])]
    n_all = len(df)
    df = df[~df["masked"] & df["pipeline"].notna()]
    df["author_mapped"] = df["author"].map(label_map)
    rows = []
    for t, lab in label_map.items():
        is_t = df["author"] == t
        if is_t.sum() == 0:
            continue
        pred_l = df["pipeline"] == lab
        recall = (is_t & pred_l).sum() / is_t.sum()
        precision = ((df["author_mapped"] == lab) & pred_l).sum() / max(pred_l.sum(), 1)
        f1 = 0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
        rows.append({"ct_author": t, "pipeline_label": lab, "n_author_cells": int(is_t.sum()),
                     "n_pipeline_cells": int(pred_l.sum()), "cmc_precision": precision,
                     "cmc_recall": recall, "cmc_f1": f1})
    return pd.DataFrame(rows), {"cells_total": n_all, "cells_excluded_masked_or_unlabelled": int(n_all - len(df)),
                                "cells_used": int(len(df))}


def support(cell_counts, study, key):
    cc = pd.read_csv(cell_counts, sep="\t")
    cc = cc[cc["study"] == study]
    wide = cc.pivot_table(index="sample", columns="ct", values="n_kept", aggfunc="sum").fillna(0)
    return pd.DataFrame({f"{key}_median_cells_per_donor": wide.median(),
                         f"{key}_total_cells": wide.sum(),
                         f"{key}_donors_with_cells": (wide > 0).sum(),
                         f"{key}_n_donors": len(wide)})


def contrast_cohorts(meta_dir, contrasts):
    """contrast -> cohorts it draws its cases from. The DE is one model over all cohorts' donors
    (no cohort term), so every contrast uses every cohort's controls; a disorder contrast takes
    its cases from the cohorts with donors of that disorder, Age and Sex from all cohorts."""
    import glob
    meta = {os.path.basename(f).replace("_sample_meta.tsv", ""): pd.read_csv(f, sep="\t", dtype=str)
            for f in sorted(glob.glob(os.path.join(meta_dir, "*_sample_meta.tsv")))}
    out = {}
    for c in contrasts:
        m = re.match(r"Disorder_(.+)_vs_Control$", c)
        if m:
            out[c] = sorted(k for k, d in meta.items()
                            if (d["Disorder"].fillna("").map(rname) == m.group(1)).any())
        else:
            out[c] = sorted(meta)
        if not out[c]:
            raise SystemExit(f"no cohort for contrast {c}")
    return out


def cohort_support(cell_counts):
    """median author cells per donor (zeros included) per cohort and author type."""
    cc = pd.read_csv(cell_counts, sep="\t")
    rows = []
    for st, g in cc.groupby("study"):
        wide = g.pivot_table(index="sample", columns="ct", values="n_kept", aggfunc="sum").fillna(0)
        rows.append(pd.DataFrame({"study": st, "ct_author": wide.columns,
                                  "median_cells_per_donor": wide.median().values, "n_donors": len(wide)}))
    return pd.concat(rows, ignore_index=True)


def mixed_rank_model(pairs, sub):
    """rank(rho) ~ rank(F1) + rank(log10 cells) + (1 | cell type), ranks standardised over the pairs."""
    import statsmodels.formula.api as smf
    d = pairs[sub].copy()
    for c in ("rho", "f1", "log10_cells"):
        r = stats.rankdata(d[c])
        d["r_" + c] = (r - r.mean()) / r.std()
    fit = smf.mixedlm("r_rho ~ r_f1 + r_log10_cells", d, groups=d["ct_author"]).fit(reml=True)
    ci = fit.conf_int()
    return [{"term": t, "beta": fit.params[k], "ci_low": ci.loc[k, 0], "ci_high": ci.loc[k, 1],
             "p": fit.pvalues[k], "n_pairs": len(d), "n_types": d["ct_author"].nunique(),
             "celltype_var": float(fit.cov_re.iloc[0, 0])}
            for t, k in (("F1", "r_f1"), ("log10 cells per donor", "r_log10_cells"))]


def rho_ci_row(label, x, y, z, controlled_for):
    r0, p0, n0 = spearman(x, y)
    r1, p1, _ = partial_spearman(x, y, z)
    row = {"label": label, "n": n0, "controlled_for": controlled_for}
    for kind, r_, p_, k in (("plain", r0, p0, 0), ("partial", r1, p1, 1)):
        zz, se = np.arctanh(r_), 1 / np.sqrt(n0 - 3 - k)
        row.update({f"{kind}_rho": r_, f"{kind}_p": p_, f"{kind}_ci_low": np.tanh(zz - 1.959964 * se),
                    f"{kind}_ci_high": np.tanh(zz + 1.959964 * se)})
    return row


def de_donor_table(per_cohort_dir, meta_dir, counts, label_map):
    """F1 and cells per donor over the donors in the DE (one limma fit per cell type over all cohorts).
    counts: scheme -> cell_counts.tsv. A donor is in the DE of a cell type when it has metadata and a
    nonzero library for that type (limma_voom_analysis.R; min_cells_per_sample = 0)."""
    import glob
    ids = set(pd.concat([pd.read_csv(f, sep="\t", dtype=str)
                         for f in sorted(glob.glob(os.path.join(meta_dir, "*_sample_meta.tsv")))])["Individual_ID"])
    rows = []
    for scheme, path in counts.items():
        cc = pd.read_csv(path, sep="\t", dtype={"sample": str})
        assert ((cc["n_kept"] > 0) == (cc["lib_size"] > 0)).all()
        cc = cc[(cc["lib_size"] > 0) & cc["sample"].isin(ids)]
        g = cc.groupby("ct")["n_kept"]
        rows.append(pd.DataFrame({"scheme": scheme, "celltype": g.size().index, "n_de_donors": g.size().values,
                                  "median_cells_per_donor": g.median().values,
                                  "n_cohorts": cc.groupby("ct")["study"].nunique().values}))
        if scheme == "author":
            de_pairs = set(zip(cc["sample"], cc["ct"]))
    tab = pd.concat(rows, ignore_index=True)
    ps = pd.concat([pd.read_csv(f, sep="\t", dtype={"sample": str})
                    for f in sorted(glob.glob(os.path.join(per_cohort_dir, "*_per_sample_f1.tsv")))], ignore_index=True)
    ps = ps[ps["ct_author"].isin(label_map)]
    f1_pairs = set(zip(ps["sample"], ps["ct_author"]))
    de_pairs = {p for p in de_pairs if p[1] in label_map}
    if f1_pairs != de_pairs:
        raise SystemExit(f"per-sample F1 donors differ from DE donors: {len(f1_pairs - de_pairs)} only in F1, "
                         f"{len(de_pairs - f1_pairs)} only in DE")
    f1 = ps.groupby("ct_author").agg(f1=("f1", "mean"), precision=("precision", "mean"), recall=("recall", "mean"),
                                     n_f1_samples=("f1", "size"))
    tab = tab.merge(f1, left_on="celltype", right_index=True, how="left")
    au = tab["scheme"] == "author"
    tab = tab[~au | tab["celltype"].isin(label_map)]
    au = tab["scheme"] == "author"
    assert tab.loc[au, "f1"].notna().all() and (tab.loc[au, "n_f1_samples"] == tab.loc[au, "n_de_donors"]).all()
    tab["log10_cells"] = np.log10(tab["median_cells_per_donor"] + 1)
    tab["rule"] = "de_donors"
    return tab


def case_cohort_values(per_cohort_dir, meta_dir, cell_counts, label_map, contrasts):
    """Superseded rule: per (contrast, author type), unweighted mean over the cohorts that contribute the
    contrast's cases of the per-cohort mean F1 and per-cohort median cells per donor (zeros included)."""
    import glob
    f1 = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(per_cohort_dir, "*_per_cohort_f1.tsv")))],
                   ignore_index=True)
    sup = cohort_support(cell_counts)
    coh = f1.merge(sup, on=["study", "ct_author"], how="outer", validate="one_to_one")
    if coh[["f1", "median_cells_per_donor"]].isna().any().any():
        raise SystemExit("per-cohort F1 or cells per donor missing:\n" +
                         coh[coh[["f1", "median_cells_per_donor"]].isna().any(axis=1)].to_string())
    coh = coh[coh["ct_author"].isin(label_map)]
    cmap = contrast_cohorts(meta_dir, contrasts)
    ci = coh.set_index(["study", "ct_author"])
    vals = {}
    for c, cs in cmap.items():
        for t in label_map:
            v = ci.loc[[(s, t) for s in cs if (s, t) in ci.index]]
            vals[(c, t)] = (v["f1"].mean(), v["median_cells_per_donor"].mean(), ";".join(cs), len(cs))
    return coh, cmap, vals


def per_type_tests(pt, subsets, y="mean_rho"):
    out = {}
    for sub, m in subsets.items():
        d = pt[m]
        out[sub] = {"n": int(m.sum()),
                    "rho_f1": spearman(d["f1"], d[y])[0], "rho_cells": spearman(d["log10_cells"], d[y])[0],
                    "partial_f1_given_cells": partial_spearman(d["f1"], d[y], d["log10_cells"])[0],
                    "partial_f1_given_cells_p": partial_spearman(d["f1"], d[y], d["log10_cells"])[1],
                    "partial_cells_given_f1": partial_spearman(d["log10_cells"], d[y], d["f1"])[0],
                    "partial_cells_given_f1_p": partial_spearman(d["log10_cells"], d[y], d["f1"])[1]}
    return out


def per_cohort_analysis(combined, label_map, metric, tab, per_cohort_dir, meta_dir, cell_counts, outdir,
                        support_rule="de_donors", one_to_one="strict", f1_cells=None):
    """Agreement per cell type x contrast against F1 and cells per donor under support_rule (see module
    docstring). Returns rho_ci for Figure 4B (per cell type, averaged over contrasts), the tests, the
    mixed model, the contrast -> case cohorts map and the per-type table."""
    from utils import CONTRAST_DISPLAY
    sfx = "" if support_rule == "de_donors" else f"_{support_rule}"
    df = combined[combined["ct_author"].isin(label_map) & combined["contrast"].isin(CONTRAST_DISPLAY)]
    coh, cmap, cvals = case_cohort_values(per_cohort_dir, meta_dir, cell_counts, label_map,
                                          sorted(df["contrast"].unique()))
    coh.to_csv(os.path.join(outdir, "per_cohort_celltype.tsv"), sep="\t", index=False, float_format="%.6g")
    fa = f1_cells[f1_cells["scheme"] == "author"].set_index("celltype")
    rows = []
    for r in df.itertuples():
        if support_rule == "de_donors":
            f, cells, cs, nc = fa.loc[r.ct_author, "f1"], fa.loc[r.ct_author, "median_cells_per_donor"], "all DE donors", int(fa.loc[r.ct_author, "n_cohorts"])
        else:
            f, cells, cs, nc = cvals[(r.contrast, r.ct_author)]
        rows.append({"contrast": r.contrast, "contrast_display": CONTRAST_DISPLAY[r.contrast],
                     "ct_author": r.ct_author, "rho": getattr(r, metric), "cohorts": cs,
                     "n_cohorts": nc, "f1": f, "median_cells_per_donor": cells})
    pairs = pd.DataFrame(rows).dropna(subset=["rho"])
    pairs["log10_cells"] = np.log10(pairs["median_cells_per_donor"] + 1)
    oto = dict(zip(tab["ct_author"], tab["one_to_one_strict" if one_to_one == "strict" else "benchmark_own_label"]))
    pairs["one_to_one"] = pairs["ct_author"].map(oto)
    pairs.to_csv(os.path.join(outdir, f"pairs_per_cohort{sfx}.tsv"), sep="\t", index=False, float_format="%.6g")

    per_type = (pairs.groupby("ct_author")
                .agg(mean_rho=("rho", "mean"), f1=("f1", "mean"), log10_cells=("log10_cells", "mean"),
                     n_contrasts=("rho", "size"), one_to_one=("one_to_one", "first")).reset_index())
    chk = per_type.set_index("ct_author")["mean_rho"] - tab.set_index("ct_author")["mean_rho"]
    assert np.allclose(chk.dropna(), 0, atol=1e-12) and chk.notna().sum() == len(tab)
    per_type.to_csv(os.path.join(outdir, f"celltype_per_cohort{sfx}.tsv"), sep="\t", index=False, float_format="%.6g")

    tests, mm = [], []
    for level, d in (("cell type x contrast", pairs), ("cell type (mean over contrasts)", per_type)):
        y = d["rho"] if "rho" in d else d["mean_rho"]
        for sub, m in (("all", pd.Series(True, index=d.index)), ("one-to-one types", d["one_to_one"].astype(bool))):
            for name, fn, cols in (("agreement ~ F1", spearman, (d["f1"], y)),
                                   ("agreement ~ log10 cells per donor", spearman, (d["log10_cells"], y)),
                                   ("agreement ~ F1 | log10 cells per donor", partial_spearman, (d["f1"], y, d["log10_cells"])),
                                   ("agreement ~ log10 cells per donor | F1", partial_spearman, (d["log10_cells"], y, d["f1"])),
                                   ("F1 ~ log10 cells per donor", spearman, (d["f1"], d["log10_cells"]))):
                r, pv, n = fn(*[c[m] for c in cols])
                tests.append({"level": level, "subset": sub, "test": name, "spearman_rho": r, "p": pv, "n": n,
                              "n_types": d.loc[m, "ct_author"].nunique()})
            if level.startswith("cell type x"):
                for row in mixed_rank_model(d, m):
                    mm.append({"subset": sub, **row})
    tests = pd.DataFrame(tests)
    tests.insert(0, "support_rule", support_rule)
    tests.insert(1, "one_to_one_rule", one_to_one)
    tests.to_csv(os.path.join(outdir, f"correlations_per_cohort{sfx}.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(mm).to_csv(os.path.join(outdir, f"mixed_model_per_cohort{sfx}.tsv"), sep="\t", index=False, float_format="%.6g")

    # Figure 4B: per cell type, averaged over contrasts (points independent across types)
    pt, o = per_type, per_type["one_to_one"].astype(bool)
    rho_ci = pd.DataFrame([
        rho_ci_row(f"Cells per donor, {len(pt)} types", pt["log10_cells"], pt["mean_rho"], pt["f1"], "F1"),
        rho_ci_row(f"F1, {len(pt)} types", pt["f1"], pt["mean_rho"], pt["log10_cells"], "log10 cells per donor"),
        rho_ci_row(f"F1, {int(o.sum())} one-to-one types", pt["f1"][o], pt["mean_rho"][o], pt["log10_cells"][o],
                   "log10 cells per donor")])
    rho_ci.to_csv(os.path.join(outdir, f"fig4b_rho_ci{sfx}.tsv"), sep="\t", index=False, float_format="%.6g")

    # every rule side by side (decision log): de_donors, case_cohorts, CMC pooled; three type subsets
    t = tab.set_index("ct_author")
    alt_vals = {"de_donors": (fa["f1"], fa["log10_cells"]),
                "case_cohorts": (pd.Series({k: np.mean([cvals[(c, k)][0] for c in df.loc[df["ct_author"] == k, "contrast"]])
                                            for k in t.index}),
                                 pd.Series({k: np.mean([np.log10(cvals[(c, k)][1] + 1) for c in df.loc[df["ct_author"] == k, "contrast"]])
                                            for k in t.index})),
                "cmc_pooled": (t["cmc_f1"], t["log10_author_median_cells"])}
    alt = []
    for rule, (f, c) in alt_vals.items():
        d = pd.DataFrame({"mean_rho": t["mean_rho"], "f1": f.reindex(t.index), "log10_cells": c.reindex(t.index)})
        subsets = {"all": pd.Series(True, index=t.index), "own_label": t["benchmark_own_label"].astype(bool),
                   "strict": t["one_to_one_strict"].astype(bool)}
        for sub, v in per_type_tests(d, subsets).items():
            alt.append({"rule": rule, "types": sub, **v})
    alt = pd.DataFrame(alt)
    alt.to_csv(os.path.join(outdir, "rule_alternatives.tsv"), sep="\t", index=False, float_format="%.6g")
    return rho_ci, tests, pd.DataFrame(mm), cmap, per_type


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--combined_corr", required=True, help="subclass combined_filtered_corr.tsv")
    p.add_argument("--params", required=True)
    p.add_argument("--label_map", required=True, help="author -> pipeline subclass map JSON")
    p.add_argument("--benchmark_f1", required=True)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cta", required=True, help="<study>_subclass_cell_type.tsv (pipeline v2.1.2)")
    p.add_argument("--mask", required=True, help="<study>_mask.tsv")
    p.add_argument("--author_cell_counts", required=True, help="pseudobulks/author/ct_pseudobulks/cell_counts.tsv")
    p.add_argument("--pipeline_cell_counts", required=True, help="pseudobulks/pavlab_subclass/ct_pseudobulks/cell_counts.tsv")
    p.add_argument("--study", default="CMC")
    p.add_argument("--level", default="subclass")
    p.add_argument("--metric", default="spearman_log2FoldChange")
    p.add_argument("--outdir", required=True)
    p.add_argument("--null_dir", required=True,
                   help="benchmark-post-hoc/de_permutation_null (pairs.tsv, null_summary.tsv) for panel E")
    p.add_argument("--per_cohort_dir", default=None,
                   help="per_cohort_f1.py outputs (<study>_per_cohort_f1.tsv); default <outdir>/per_cohort")
    p.add_argument("--meta_dir", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "gemma", "metadata"),
                   help="<cohort>_sample_meta.tsv files used by the DE (contrast -> cohorts); default as nextflow.config gemma_meta_dir")
    p.add_argument("--class_cell_counts", default=None,
                   help="pseudobulks/pavlab_class/ct_pseudobulks/cell_counts.tsv; default: sibling of --pipeline_cell_counts")
    p.add_argument("--support_rule", choices=["de_donors", "case_cohorts"], default="de_donors",
                   help="which donors F1 and cells per donor come from (module docstring)")
    p.add_argument("--one_to_one", choices=["strict", "own_label"], default="strict",
                   help="strict: 18 types whose pipeline label pools no other type; own_label: 20 types named as their label")
    p.add_argument("--qa_dir", default=None, help="write layout-alignment JSON and a QA-only PDF here")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    a.per_cohort_dir = a.per_cohort_dir or os.path.join(a.outdir, "per_cohort")
    a.class_cell_counts = a.class_cell_counts or a.pipeline_cell_counts.replace("pavlab_subclass", "pavlab_class")
    assert a.class_cell_counts != a.pipeline_cell_counts

    label_map = load_label_map(a.params, a.label_map)
    combined = pd.read_csv(a.combined_corr, sep="\t")
    pivot = contrast_pivot(combined, label_map, a.metric)
    tab = pd.DataFrame({"ct_author": pivot.columns, "display_name": [author_display(c) for c in pivot.columns],
                        "pipeline_label": [label_map[c] for c in pivot.columns],
                        "mean_rho": pivot.mean().values, "n_contrasts": pivot.notna().sum().values})
    a2l = dict(zip(tab["ct_author"], tab["pipeline_label"]))
    bf1 = load_benchmark_f1(a.benchmark_f1, a.level, a2l)
    own = own_label_types(a2l, a.level)
    tab["benchmark_f1"] = tab["ct_author"].map(bf1)
    tab["benchmark_own_label"] = tab["ct_author"].isin(own)
    # strict one-to-one: the type's pipeline label pools no other author type
    tab["one_to_one_strict"] = tab["pipeline_label"].map(tab["pipeline_label"].value_counts()) == 1

    f1_df, cell_info = cmc_f1(a.h5ad, a.cta, a.mask, label_map)
    tab = tab.merge(f1_df.drop(columns="pipeline_label"), on="ct_author", how="left")
    sa = support(a.author_cell_counts, a.study, "author")
    sp = support(a.pipeline_cell_counts, a.study, "pipeline")
    tab = tab.merge(sa, left_on="ct_author", right_index=True, how="left")
    tab = tab.merge(sp, left_on="pipeline_label", right_index=True, how="left")
    tab["log10_author_median_cells"] = np.log10(tab["author_median_cells_per_donor"] + 1)
    tab["log10_pipeline_median_cells"] = np.log10(tab["pipeline_median_cells_per_donor"] + 1)
    tab = tab.sort_values("mean_rho", ascending=False)

    y = tab["mean_rho"]
    logs = tab["log10_author_median_cells"]
    # author types whose pipeline label pools several author types
    n_per_label = tab["pipeline_label"].map(tab["pipeline_label"].value_counts())
    tab["pooled_pipeline_label"] = n_per_label > 1
    subsets = {
        "all": pd.Series(True, index=tab.index),
        # 7 types: name differs from their pooled pipeline label (IT subtypes, SST Chodl, LAMP5 LHX6)
        "excluding types named differently from their pooled label": tab["benchmark_own_label"],
        # 9 types: every type in a pooled label (adds SST and LAMP5)
        "excluding all types in pooled labels": ~tab["pooled_pipeline_label"],
    }
    tests = []
    def add(name, sub, fn, *cols):
        m = subsets[sub]
        r, pv, n = fn(*[c[m] for c in cols])
        excl = ";".join(tab.loc[~m, "display_name"])
        tests.append({"test": name, "subset": sub, "spearman_rho": r, "p": pv, "n": n, "excluded": excl})
    bf, cf, lp = tab["benchmark_f1"], tab["cmc_f1"], tab["log10_pipeline_median_cells"]
    for sub in subsets:
        add("rho ~ CMC F1", sub, spearman, cf, y)
        add("rho ~ benchmark F1", sub, spearman, bf, y)
        add("rho ~ log10 median author cells per donor", sub, spearman, logs, y)
        add("rho ~ CMC F1 | log10 support (partial)", sub, partial_spearman, cf, y, logs)
        add("rho ~ benchmark F1 | log10 support (partial)", sub, partial_spearman, bf, y, logs)
        add("CMC F1 ~ log10 median author cells", sub, spearman, cf, logs)
        add("benchmark F1 ~ CMC F1", sub, spearman, bf, cf)
    add("rho ~ log10 median pipeline-label cells per donor", "all", spearman, lp, y)
    tests = pd.DataFrame(tests)
    # CMC-only analysis with cells pooled across samples (superseded; kept for comparison)
    tests.to_csv(os.path.join(a.outdir, "correlations_cmc_pooled.tsv"), sep="\t", index=False)

    # Former Figure 4B (CMC only, pooled cells): plain and partial Spearman rho with 95% CIs (Fisher z,
    # SE = 1/sqrt(n - 3 - k), k = number of controlled variables); kept as fig4b_rho_ci_cmc_pooled.tsv.
    own_sub = subsets["excluding types named differently from their pooled label"]
    b_rows = []
    for label, x, z, m in [(f"Cells per donor, {len(tab)} types", logs, cf, subsets["all"]),
                           (f"{a.study} F1, {len(tab)} types", cf, logs, subsets["all"]),
                           (f"{a.study} F1, {int(own_sub.sum())} one-to-one types", cf, logs, own_sub)]:
        r0, p0, n0 = spearman(x[m], y[m])
        r1, p1, n1 = partial_spearman(x[m], y[m], z[m])
        row = {"label": label, "n": n0, "controlled_for": "CMC F1" if x is logs else "log10 median author cells per donor"}
        for kind, r_, p_, k in (("plain", r0, p0, 0), ("partial", r1, p1, 1)):
            zz, se = np.arctanh(r_), 1 / np.sqrt(n0 - 3 - k)
            row.update({f"{kind}_rho": r_, f"{kind}_p": p_, f"{kind}_ci_low": np.tanh(zz - 1.959964 * se),
                        f"{kind}_ci_high": np.tanh(zz + 1.959964 * se)})
        b_rows.append(row)
    rho_ci = pd.DataFrame(b_rows)
    # rows 2 and 3 must equal correlations.tsv
    for i, sub in ((1, "all"), (2, "excluding types named differently from their pooled label")):
        for kind, test in (("plain", "rho ~ CMC F1"), ("partial", "rho ~ CMC F1 | log10 support (partial)")):
            ref = tests[(tests["subset"] == sub) & (tests["test"] == test)].iloc[0]["spearman_rho"]
            assert abs(rho_ci.loc[i, f"{kind}_rho"] - ref) < 1e-12
    rho_ci.to_csv(os.path.join(a.outdir, "fig4b_rho_ci_cmc_pooled.tsv"), sep="\t", index=False, float_format="%.6g")
    tab.to_csv(os.path.join(a.outdir, "celltype_summary_cmc_pooled.tsv"), sep="\t", index=False)

    # Current Figure 4: one F1 / cells-per-donor table (celltype_f1_cells.tsv) for 4A, 4B, the scatter
    # and matched_effect_adjusted.py (4C, S7B)
    f1_cells = de_donor_table(a.per_cohort_dir, a.meta_dir,
                              {"author": a.author_cell_counts, "pipeline_subclass": a.pipeline_cell_counts,
                               "pipeline_class": a.class_cell_counts}, label_map)
    sfx = "" if a.support_rule == "de_donors" else f"_{a.support_rule}"
    f1_cells.to_csv(os.path.join(a.outdir, "celltype_f1_cells.tsv"), sep="\t", index=False, float_format="%.6g")
    rho_ci, tests_pc, mm, cmap, pt = per_cohort_analysis(combined, label_map, a.metric, tab, a.per_cohort_dir,
                                                         a.meta_dir, a.author_cell_counts, a.outdir,
                                                         a.support_rule, a.one_to_one, f1_cells)
    pt = pt.set_index("ct_author")
    tab["f1"] = tab["ct_author"].map(pt["f1"])
    tab["log10_cells"] = tab["ct_author"].map(pt["log10_cells"])
    tab["one_to_one"] = tab["ct_author"].map(pt["one_to_one"]).astype(bool)
    tab["support_rule"] = a.support_rule
    tab.to_csv(os.path.join(a.outdir, f"celltype_summary{sfx}.tsv"), sep="\t", index=False)
    corr_pt = tests_pc[tests_pc["level"] == "cell type (mean over contrasts)"]

    from utils import PIPELINE_TO_BENCH
    from fig_de_agreement import plot_scatter
    plot_scatter(tab, list(tab["ct_author"]), corr_pt,
                 os.path.join(a.outdir, f"fig_de_agreement_scatter_{a.level}{sfx}.png"))
    plot_figure(pivot, tab, label_map, PIPELINE_TO_BENCH[a.level],
                os.path.join(a.outdir, f"fig_de_agreement_{a.level}{sfx}.png"),
                rho_ci, a.null_dir, study=a.study, qa_dir=a.qa_dir)

    prov = {"created": datetime.datetime.now().astimezone().isoformat(), "study": a.study, "level": a.level,
            "cells": cell_info, "support_rule": a.support_rule, "one_to_one": a.one_to_one,
            "one_to_one_types": sorted(tab.loc[tab["one_to_one"], "ct_author"]),
            "decisions": "psychENCODE-reanalysis/docs/DE_AGREEMENT_DECISIONS.md",
            "contrast_cohorts (case_cohorts rule only)": cmap, "inputs": {}, "outputs": {}}
    import glob
    for key, pat in (("per_cohort", os.path.join(a.per_cohort_dir, "*_per_*_f1.tsv")),
                     ("per_cohort", os.path.join(a.per_cohort_dir, "*_cells.json")),
                     ("meta", os.path.join(a.meta_dir, "*_sample_meta.tsv"))):
        for f in sorted(glob.glob(pat)):
            prov["inputs"][f"{key}/{os.path.basename(f)}"] = {"path": os.path.abspath(f), "sha256": sha256(f)}
    for k in ("combined_corr", "label_map", "benchmark_f1", "cta", "mask", "author_cell_counts", "pipeline_cell_counts",
              "class_cell_counts"):
        f = getattr(a, k)
        prov["inputs"][k] = {"path": os.path.abspath(f), "sha256": sha256(f),
                             "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(f)).astimezone().isoformat()}
    for f in ("pairs.tsv", "null_summary.tsv", "matched_effect_adjusted.tsv", "manifest.json"):
        f = os.path.join(a.null_dir, f)
        prov["inputs"][f"null_{os.path.basename(f)}"] = {"path": os.path.abspath(f), "sha256": sha256(f)}
    from fig_de_agreement import PALETTE_TSV
    prov["inputs"]["palette"] = {"path": PALETTE_TSV, "sha256": sha256(PALETTE_TSV)}
    prov["inputs"]["h5ad"] = {"path": os.path.abspath(a.h5ad),
                              "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(a.h5ad)).astimezone().isoformat()}
    bj = os.path.splitext(a.benchmark_f1)[0] + ".json"
    if os.path.exists(bj):
        prov["benchmark_f1_provenance"] = json.load(open(bj))
    prov["data_version"] = {"pipeline_labels": "sc-annotation-pipeline 2.1.2", "de": "results_v2.1.2_limma"}
    for f in (f"celltype_summary{sfx}.tsv", "celltype_summary_cmc_pooled.tsv", "correlations_cmc_pooled.tsv",
              "fig4b_rho_ci_cmc_pooled.tsv", "per_cohort_celltype.tsv", "celltype_f1_cells.tsv", "rule_alternatives.tsv",
              f"pairs_per_cohort{sfx}.tsv", f"celltype_per_cohort{sfx}.tsv", f"correlations_per_cohort{sfx}.tsv",
              f"mixed_model_per_cohort{sfx}.tsv", f"fig4b_rho_ci{sfx}.tsv", f"fig_de_agreement_{a.level}{sfx}.png",
              f"fig_de_agreement_scatter_{a.level}{sfx}.png"):
        prov["outputs"][f] = {"path": os.path.abspath(os.path.join(a.outdir, f)), "sha256": sha256(os.path.join(a.outdir, f))}
    json.dump(prov, open(os.path.join(a.outdir, f"manifest{sfx}.json"), "w"), indent=2)
    print(tests.drop(columns="excluded").round(3).to_string(index=False))
    print(tests_pc.round(3).to_string(index=False))
    print(mm.round(3).to_string(index=False))
    print(pd.read_csv(os.path.join(a.outdir, "rule_alternatives.tsv"), sep="\t").round(3).to_string(index=False))


if __name__ == "__main__":
    main()
