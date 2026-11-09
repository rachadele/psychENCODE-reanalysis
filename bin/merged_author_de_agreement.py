#!/usr/bin/env python3
"""DE agreement of the merged-author arm with the pipeline-subclass arm, against F1 and cells per donor.

Merged-author arm: author labels mapped to pipeline subclass labels (--label_map) before the
pseudobulk sums, so each merged author group matches exactly one pipeline label and no type is
excluded. Inputs come from results_v2.1.2_limma/merged_author_2026-11-09/scripts/run_merged_author.sh.

  * agreement: Spearman rho of log2FC and of the moderated t between the merged-author DE and the
    pipeline-subclass DE, per type x contrast (comparison stage, combined_filtered_corr.tsv), over
    the Figure 4 contrasts (utils.CONTRAST_DISPLAY); mean and SD over type x contrast pairs;
  * F1 and cells per donor under the de_donors rule (DE_AGREEMENT_DECISIONS.md): plain mean of
    per-sample F1 against merged author labels (per_cohort_f1.py --merge_author) over the donors in
    each type's DE, and median kept cells per DE donor (de_agreement_vs_annotation.de_donor_table);
  * across merged types: plain Spearman of mean agreement with F1 and with log10(cells + 1), and
    partial Spearman of each controlling for the other, next to the current 27-type and strict
    18-type values (correlations_per_cohort.tsv);
  * checks: for types that pool nothing, the merged-arm DE, rho, F1 and cells per donor equal the
    current author arm.
Writes TSVs, checks.json and manifest.json (sha256 of every input and output) to --outdir.
"""
import argparse
import datetime
import glob
import json
import os

import numpy as np
import pandas as pd

from de_agreement_vs_annotation import de_donor_table, partial_spearman, sha256, spearman
from utils import CONTRAST_DISPLAY

METRICS = {"log2FC": "spearman_log2FoldChange", "t": "spearman_stat"}


def mtime(f):
    return datetime.datetime.fromtimestamp(os.path.getmtime(f)).astimezone().isoformat()


def agreement_summary(pairs, groups):
    rows = []
    for name, m in groups.items():
        d = pairs[m]
        row = {"set": name, "n_types": d["celltype"].nunique(), "n_pairs": len(d)}
        for k, col in METRICS.items():
            row[f"mean_rho_{k}"] = d[col].mean()
            row[f"sd_rho_{k}"] = d[col].std()
        rows.append(row)
    return pd.DataFrame(rows)


def tests(pt, y):
    out = []
    for name, fn, cols in (("agreement ~ F1", spearman, (pt["f1"], pt[y])),
                           ("agreement ~ log10 cells per donor", spearman, (pt["log10_cells"], pt[y])),
                           ("agreement ~ F1 | log10 cells per donor", partial_spearman, (pt["f1"], pt[y], pt["log10_cells"])),
                           ("agreement ~ log10 cells per donor | F1", partial_spearman, (pt["log10_cells"], pt[y], pt["f1"])),
                           ("F1 ~ log10 cells per donor", spearman, (pt["f1"], pt["log10_cells"]))):
        r, p, n = fn(*cols)
        out.append({"test": name, "spearman_rho": r, "p": p, "n": n})
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--merged_corr", required=True, help="merged arm combined_filtered_corr.tsv")
    p.add_argument("--current_corr", required=True, help="current subclass combined_filtered_corr.tsv")
    p.add_argument("--current_correlations", required=True, help="current correlations_per_cohort.tsv")
    p.add_argument("--current_f1_cells", required=True, help="current celltype_f1_cells.tsv")
    p.add_argument("--label_map", required=True, help="author -> pipeline subclass map JSON")
    p.add_argument("--per_sample_dir", required=True, help="per_cohort_f1.py --merge_author outputs")
    p.add_argument("--merged_cell_counts", required=True)
    p.add_argument("--current_author_de", required=True)
    p.add_argument("--merged_author_de", required=True)
    p.add_argument("--meta_dir", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    amap = json.load(open(a.label_map))
    members = {}
    for t, lab in amap.items():
        members.setdefault(lab, []).append(t)
    ident = {lab: lab for lab in members}

    # agreement per type x contrast
    mc = pd.read_csv(a.merged_corr, sep="\t")
    mc = mc[mc["contrast"].isin(CONTRAST_DISPLAY) & (mc["ct_author"] == mc["ct_sc_pipeline"]) & mc["ct_author"].isin(ident)]
    mc = mc.rename(columns={"ct_author": "celltype"})
    mc["n_members"] = mc["celltype"].map(lambda x: len(members[x]))
    cc = pd.read_csv(a.current_corr, sep="\t")
    cc = cc[cc["contrast"].isin(CONTRAST_DISPLAY) & cc["ct_author"].isin(amap)].rename(columns={"ct_author": "celltype"})
    cc["n_members"] = cc["ct_sc_pipeline"].map(lambda x: len(members[x]))
    keep = ["celltype", "contrast", "n_members"] + list(METRICS.values())
    pairs = pd.concat([mc[keep].assign(arm="merged"), cc[keep].assign(arm="current")], ignore_index=True)
    pairs.to_csv(os.path.join(a.outdir, "pairs_agreement.tsv"), sep="\t", index=False, float_format="%.6g")
    m, c = pairs["arm"] == "merged", pairs["arm"] == "current"
    one = pairs["n_members"] == 1
    agree = agreement_summary(pairs, {
        f"merged, all {mc['celltype'].nunique()} types": m,
        "merged, 3 merged groups": m & ~one,
        f"current, all {cc['celltype'].nunique()} types": c,
        f"current, strict {cc.loc[cc['n_members'] == 1, 'celltype'].nunique()} one-to-one types": c & one,
        "current, 9 pooled types": c & ~one})
    agree.to_csv(os.path.join(a.outdir, "agreement_summary.tsv"), sep="\t", index=False, float_format="%.6g")

    # F1 and cells per donor, de_donors rule
    # scheme "author": de_donor_table checks that the per-sample F1 donors equal the DE donors
    fc = de_donor_table(a.per_sample_dir, a.meta_dir, {"author": a.merged_cell_counts}, ident)
    assert (fc["n_f1_samples"] == fc["n_de_donors"]).all()

    pt = (mc.groupby("celltype").agg(mean_rho_log2FC=(METRICS["log2FC"], "mean"), mean_rho_t=(METRICS["t"], "mean"),
                                     n_contrasts=(METRICS["log2FC"], "size")).reset_index())
    pt = pt.merge(fc[["celltype", "f1", "median_cells_per_donor", "log10_cells", "n_de_donors", "n_cohorts"]],
                  on="celltype", how="left", validate="one_to_one")
    pt.insert(1, "author_types", pt["celltype"].map(lambda x: ";".join(members[x])))
    pt.insert(2, "n_members", pt["celltype"].map(lambda x: len(members[x])))
    assert pt[["f1", "log10_cells"]].notna().all().all()
    pt = pt.sort_values("mean_rho_log2FC", ascending=False)
    pt.to_csv(os.path.join(a.outdir, "celltype_merged.tsv"), sep="\t", index=False, float_format="%.6g")

    # correlations, merged vs current
    cur = pd.read_csv(a.current_correlations, sep="\t")
    cur = cur[cur["level"] == "cell type (mean over contrasts)"]
    rows = []
    for y, metric in (("mean_rho_log2FC", "log2FC"), ("mean_rho_t", "t")):
        for r in tests(pt, y):
            rows.append({"set": f"merged, {len(pt)} types", "agreement": metric, **r})
    for sub, name in (("all", "current, 27 types"), ("one-to-one types", "current, strict 18 types")):
        for r in cur[cur["subset"] == sub].itertuples():
            rows.append({"set": name, "agreement": "log2FC", "test": r.test, "spearman_rho": r.spearman_rho,
                         "p": r.p, "n": r.n})
    corr = pd.DataFrame(rows)
    corr.to_csv(os.path.join(a.outdir, "correlations_merged_vs_current.tsv"), sep="\t", index=False, float_format="%.6g")

    # checks: unpooled types identical to the current author arm
    cf = pd.read_csv(a.current_f1_cells, sep="\t")
    cf = cf[cf["scheme"] == "author"].set_index("celltype")
    single = {lab: ms[0] for lab, ms in members.items() if len(ms) == 1}
    ccx = cc.set_index(["celltype", "contrast"])
    chk = []
    for lab, t in single.items():
        r = pt.set_index("celltype").loc[lab]
        rho_cur = ccx.loc[t, METRICS["log2FC"]] if t in ccx.index.get_level_values(0) else pd.Series(dtype=float)
        rho_new = mc.set_index(["celltype", "contrast"]).loc[lab, METRICS["log2FC"]]
        chk.append({"celltype": lab, "author_type": t,
                    "de_identical": de_identical(a.current_author_de, t, a.merged_author_de, lab),
                    "rho_max_abs_diff": float((rho_new - rho_cur.reindex(rho_new.index)).abs().max()),
                    "f1_abs_diff": float(abs(r["f1"] - cf.loc[t, "f1"])),
                    "cells_equal": bool(r["median_cells_per_donor"] == cf.loc[t, "median_cells_per_donor"]),
                    "donors_equal": bool(r["n_de_donors"] == cf.loc[t, "n_de_donors"])})
    chk = pd.DataFrame(chk)
    chk.to_csv(os.path.join(a.outdir, "checks_unpooled_types.tsv"), sep="\t", index=False, float_format="%.6g")
    summary = {"unpooled_types": len(chk), "de_identical": int(chk["de_identical"].sum()),
               "rho_max_abs_diff": float(chk["rho_max_abs_diff"].max()), "f1_max_abs_diff": float(chk["f1_abs_diff"].max()),
               "cells_equal": int(chk["cells_equal"].sum()), "donors_equal": int(chk["donors_equal"].sum())}
    json.dump(summary, open(os.path.join(a.outdir, "checks.json"), "w"), indent=1)

    prov = {"created": datetime.datetime.now().astimezone().isoformat(),
            "data_version": {"pipeline_labels": "sc-annotation-pipeline 2.1.2",
                             "de": "results_v2.1.2_limma (pipeline arm) + merged_author_2026-11-09 (author arm)"},
            "status": "pending Rachel's decision (DE_AGREEMENT_DECISIONS.md, 2026-11-09)",
            "inputs": {}, "outputs": {}}
    for k in ("merged_corr", "current_corr", "current_correlations", "current_f1_cells", "label_map", "merged_cell_counts"):
        f = getattr(a, k)
        prov["inputs"][k] = {"path": os.path.abspath(f), "sha256": sha256(f), "mtime": mtime(f)}
    for f in sorted(glob.glob(os.path.join(a.per_sample_dir, "*_per_sample_f1.tsv"))) + \
            sorted(glob.glob(os.path.join(a.meta_dir, "*_sample_meta.tsv"))):
        prov["inputs"][os.path.basename(f)] = {"path": os.path.abspath(f), "sha256": sha256(f), "mtime": mtime(f)}
    for f in sorted(os.listdir(a.outdir)):
        if f.endswith((".tsv", ".json")) and f != "manifest.json":
            prov["outputs"][f] = {"path": os.path.abspath(os.path.join(a.outdir, f)),
                                  "sha256": sha256(os.path.join(a.outdir, f))}
    json.dump(prov, open(os.path.join(a.outdir, "manifest.json"), "w"), indent=2)

    pd.set_option("display.width", 200)
    print(agree.round(3).to_string(index=False))
    print(pt.drop(columns="author_types").round(3).to_string(index=False))
    print(corr.round(4).to_string(index=False))
    print(json.dumps(summary))


def de_identical(a_dir, a_ct, b_dir, b_ct):
    """DE of author type a_ct (current arm) equals DE of label b_ct (merged arm), contrast by contrast."""
    fa = sorted(glob.glob(os.path.join(a_dir, a_ct, "*", "*_results.tsv")))
    fb = sorted(glob.glob(os.path.join(b_dir, b_ct, "*", "*_results.tsv")))
    ca = [os.path.basename(os.path.dirname(f)) for f in fa]
    cb = [os.path.basename(os.path.dirname(f)) for f in fb]
    if ca != cb or not fa:
        return False
    cols = ["gene", "log2FoldChange", "stat", "pvalue"]
    return all(pd.read_csv(x, sep="\t")[cols].equals(pd.read_csv(y, sep="\t")[cols]) for x, y in zip(fa, fb))


if __name__ == "__main__":
    main()
