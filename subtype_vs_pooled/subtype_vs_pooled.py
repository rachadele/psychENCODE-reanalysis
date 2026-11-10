#!/usr/bin/env python3
"""Does splitting a pooled pipeline label into author subtypes change the DE signal?

For the three pooled groups (L2/3-6 IT: five author IT types; Lamp5: Lamp5 + Lamp5_Lhx6;
Sst: Sst + Sst_Chodl) and the eight Figure 4 contrasts, Spearman rho of log2FC and of the
moderated t between:
  subtype_vs_pooled    each author subtype (unmerged author DE) and the merged author DE of its group
  subtype_vs_subtype   pairs of author subtypes within a group
  baseline_same_class  each subtype and author types outside its group of the same class
                       (glutamatergic for IT, GABAergic for Lamp5/Sst)
  baseline_all_other   each subtype and every author type outside its group
  pooled_vs_ours       merged author DE of the group and our pipeline-label DE (Figure 4 value)
Genes: inner join of the two tables, as bin/aggregate_pairwise.py. The DE is one limma fit per
cell type over all eight cohorts (no cohort term), so each contrast is one table per type.
Summaries average the pairs of a comparison within each contrast, then give mean and SD over
contrasts. Cells per donor: median n_kept over donors with metadata and a nonzero library
(the de_donors rule of bin/de_agreement_vs_annotation.py).
"""
import datetime
import glob
import hashlib
import itertools
import json
import os
import subprocess

import pandas as pd

REPO = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis"
R = f"{REPO}/results_v2.1.2_limma"
M = f"{R}/merged_author_2026-11-09"
OUT = f"{REPO}/subtype_vs_pooled"
DE = {"author": f"{R}/author_submitted_true/subclass/limma_voom/gemma",
      "merged": f"{M}/author_submitted_true/subclass/limma_voom/gemma",
      "ours": f"{R}/author_submitted_false/subclass/limma_voom/gemma"}
COUNTS = {"author": f"{R}/pseudobulks/author/ct_pseudobulks/cell_counts.tsv",
          "merged": f"{M}/pseudobulks/author_merged/ct_pseudobulks/cell_counts.tsv",
          "ours": f"{R}/pseudobulks/pavlab_subclass/ct_pseudobulks/cell_counts.tsv"}
CONTRASTS = {"Disorder_Schizophrenia_vs_Control": "SCZ", "Disorder_Bipolar.Disorder_vs_Control": "BD",
             "Disorder_ASD_vs_Control": "ASD", "Disorder_MDD_vs_Control": "MDD",
             "Disorder_PTSD_vs_Control": "PTSD", "Disorder_Williams.Syndrome_vs_Control": "WS",
             "Age_death": "Age", "Biological_Sex_male_vs_female": "Sex"}
GROUPS = {"IT": ("L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
                 ["L2.3_IT", "L4_IT", "L5_IT", "L6_IT", "L6_IT_Car3"]),
          "Lamp5": ("lamp5.GABAergic.cortical.interneuron", ["Lamp5", "Lamp5_Lhx6"]),
          "Sst": ("sst.GABAergic.cortical.interneuron", ["Sst", "Sst_Chodl"])}
GLUT = ["L2.3_IT", "L4_IT", "L5_IT", "L6_IT", "L6_IT_Car3", "L5_ET", "L6b", "L5.6_NP", "L6_CT"]
GABA = ["Lamp5", "Lamp5_Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Chandelier", "Vip"]


def path(scheme, ct, contrast):
    return f"{DE[scheme]}/{ct}/{contrast}/{contrast}_{ct}_results.tsv"


_cache, used = {}, set()


def table(scheme, ct, contrast):
    p = path(scheme, ct, contrast)
    if p not in _cache:
        if not os.path.exists(p):
            raise SystemExit(f"missing DE table: {p}")
        _cache[p] = pd.read_csv(p, sep="\t", usecols=["gene", "log2FoldChange", "stat"]).set_index("gene")
        used.add(p)
    return _cache[p]


def rho(a, b):
    j = a.join(b, how="inner", lsuffix="_a", rsuffix="_b")
    return (j["log2FoldChange_a"].corr(j["log2FoldChange_b"], method="spearman"),
            j["stat_a"].corr(j["stat_b"], method="spearman"), len(j))


def cells_per_donor():
    ids = set(pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in
                         sorted(glob.glob(f"{REPO}/gemma/metadata/*_sample_meta.tsv"))])["Individual_ID"])
    out = {}
    for scheme, p in COUNTS.items():
        cc = pd.read_csv(p, sep="\t", dtype={"sample": str})
        cc = cc[(cc["lib_size"] > 0) & cc["sample"].isin(ids)]
        g = cc.groupby("ct")["n_kept"]
        for ct in g.size().index:
            out[(scheme, ct)] = (float(g.median()[ct]), int(g.size()[ct]))
    return out


def main():
    authors = sorted(os.listdir(DE["author"]))
    authors = [a for a in authors if os.path.isdir(f"{DE['author']}/{a}")]
    rows = []
    for contrast, short in CONTRASTS.items():
        def add(group, comparison, sa, a, sb, b):
            r_fc, r_t, n = rho(table(sa, a, contrast), table(sb, b, contrast))
            rows.append({"contrast": short, "group": group, "comparison": comparison, "type_a": a, "scheme_a": sa,
                         "type_b": b, "scheme_b": sb, "rho_log2FC": r_fc, "rho_t": r_t, "n_genes": n})
        for g, (label, subs) in GROUPS.items():
            for s in subs:
                add(g, "subtype_vs_pooled", "author", s, "merged", label)
            for s1, s2 in itertools.combinations(subs, 2):
                add(g, "subtype_vs_subtype", "author", s1, "author", s2)
            same_class = GLUT if g == "IT" else GABA
            for s in subs:
                for o in authors:
                    if o in subs:
                        continue
                    add(g, "baseline_all_other", "author", s, "author", o)
                    if o in same_class:
                        add(g, "baseline_same_class", "author", s, "author", o)
            add(g, "pooled_vs_ours", "merged", label, "ours", label)
    pairs = pd.DataFrame(rows)
    pairs.to_csv(f"{OUT}/pairwise_rho.tsv", sep="\t", index=False, float_format="%.6g")

    cpd = cells_per_donor()

    def summarise(df, keys):
        per_c = df.groupby(keys + ["contrast"])[["rho_log2FC", "rho_t"]].mean().reset_index()
        s = per_c.groupby(keys).agg(mean_rho_log2FC=("rho_log2FC", "mean"), sd_rho_log2FC=("rho_log2FC", "std"),
                                    mean_rho_t=("rho_t", "mean"), sd_rho_t=("rho_t", "std"),
                                    n_contrasts=("contrast", "nunique")).reset_index()
        n_pairs = df.groupby(keys + ["contrast"]).size().groupby(keys).max().rename("pairs_per_contrast")
        return s.merge(n_pairs.reset_index(), on=keys)

    order = ["subtype_vs_pooled", "subtype_vs_subtype", "baseline_same_class", "baseline_all_other", "pooled_vs_ours"]
    by_group = summarise(pairs, ["group", "comparison"]).assign(level="group", type_a="all")
    by_type = summarise(pairs[pairs["comparison"] != "subtype_vs_subtype"], ["group", "comparison", "type_a"]).assign(level="type")
    by_pair = summarise(pairs[pairs["comparison"] == "subtype_vs_subtype"], ["group", "comparison", "type_a", "type_b"]).assign(level="pair")
    summ = pd.concat([by_group, by_type, by_pair], ignore_index=True)
    summ["type_b"] = summ.get("type_b").fillna("") if "type_b" in summ else ""

    def cells(row, col):
        t = row[col]
        if not t or t == "all":
            t = GROUPS[row["group"]][0] if row["comparison"] in ("subtype_vs_pooled", "pooled_vs_ours") or t == "all" else t
            return cpd.get(("merged", t), (float("nan"), 0))
        scheme = "merged" if t == GROUPS[row["group"]][0] else "author"
        return cpd.get((scheme, t), (float("nan"), 0))

    summ["median_cells_per_donor_a"] = [cells(r, "type_a")[0] if r["level"] != "group" else float("nan") for _, r in summ.iterrows()]
    summ["n_de_donors_a"] = [cells(r, "type_a")[1] if r["level"] != "group" else 0 for _, r in summ.iterrows()]
    summ["median_cells_per_donor_b"] = [cells(r, "type_b")[0] if r["level"] == "pair" else float("nan") for _, r in summ.iterrows()]
    summ["median_cells_per_donor_pooled"] = [cpd[("merged", GROUPS[g][0])][0] for g in summ["group"]]
    summ["median_cells_per_donor_ours"] = [cpd[("ours", GROUPS[g][0])][0] for g in summ["group"]]
    summ["comparison"] = pd.Categorical(summ["comparison"], order, ordered=True)
    summ["level"] = pd.Categorical(summ["level"], ["group", "type", "pair"], ordered=True)
    summ = summ.sort_values(["group", "level", "comparison", "type_a", "type_b"])
    cols = ["group", "level", "comparison", "type_a", "type_b", "mean_rho_log2FC", "sd_rho_log2FC", "mean_rho_t",
            "sd_rho_t", "n_contrasts", "pairs_per_contrast", "median_cells_per_donor_a", "n_de_donors_a",
            "median_cells_per_donor_b", "median_cells_per_donor_pooled", "median_cells_per_donor_ours"]
    summ[cols].to_csv(f"{OUT}/summary.tsv", sep="\t", index=False, float_format="%.4g")

    cells_tab = pd.DataFrame([{"group": g, "scheme": sch, "celltype": t,
                               "median_cells_per_donor": cpd[(sch, t)][0], "n_de_donors": cpd[(sch, t)][1]}
                              for g, (lab, subs) in GROUPS.items()
                              for sch, t in [("author", s) for s in subs] + [("merged", lab), ("ours", lab)]])
    cells_tab.to_csv(f"{OUT}/cells_per_donor.tsv", sep="\t", index=False)

    # check against the Figure 4 merged-arm value for the pooled groups
    ref = pd.read_csv(f"{M}/de_agreement/celltype_merged.tsv", sep="\t").set_index("celltype")
    chk = {g: {"here": float(by_group.query("group == @g and comparison == 'pooled_vs_ours'")["mean_rho_log2FC"].iloc[0]),
               "celltype_merged.tsv": float(ref.loc[lab, "mean_rho_log2FC"])} for g, (lab, _) in GROUPS.items()}
    for g, v in chk.items():
        assert abs(v["here"] - v["celltype_merged.tsv"]) < 1e-4, (g, v)

    def sha(f):
        h = hashlib.sha256()
        with open(f, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                h.update(b)
        return h.hexdigest()
    ts = lambda f: datetime.datetime.fromtimestamp(os.path.getmtime(f)).astimezone().isoformat(timespec="seconds")
    rec = lambda fs: {os.path.relpath(f, REPO): {"sha256": sha(f), "mtime": ts(f)} for f in sorted(fs)}
    git = lambda *a: subprocess.check_output(["git", "-C", REPO, *a], text=True).strip()
    inputs = sorted(used) + list(COUNTS.values()) + [f"{M}/de_agreement/celltype_merged.tsv"] \
        + sorted(glob.glob(f"{REPO}/gemma/metadata/*_sample_meta.tsv"))
    man = {"created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
           "data_version": "pipeline labels sc-annotation-pipeline 2.1.2; DE results_v2.1.2_limma (limma-voom, "
                           "27 author types, our subclass labels) and merged_author_2026-11-09 (21 merged author types)",
           "repo": {"branch": git("rev-parse", "--abbrev-ref", "HEAD"), "commit_at_run": git("rev-parse", "HEAD")},
           "script": {os.path.relpath(__file__, REPO): {"sha256": sha(__file__)}},
           "check_pooled_vs_ours_equals_fig4": chk,
           "inputs": rec(inputs),
           "outputs": rec([f"{OUT}/pairwise_rho.tsv", f"{OUT}/summary.tsv", f"{OUT}/cells_per_donor.tsv"])}
    json.dump(man, open(f"{OUT}/manifest.json", "w"), indent=1)
    pd.set_option("display.width", 250)
    print(summ.loc[summ["level"] != "pair", cols[:9] + ["median_cells_per_donor_a"]].round(3).to_string(index=False))
    print(summ.loc[summ["level"] == "pair", cols[:9] + ["median_cells_per_donor_a", "median_cells_per_donor_b"]].round(3).to_string(index=False))
    print(cells_tab.to_string(index=False))


if __name__ == "__main__":
    main()
