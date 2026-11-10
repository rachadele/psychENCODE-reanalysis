"""Split-half reproducibility: Spearman rho of moderated t between halves A and B.

Per seed, type and contrast: rho over genes present in both halves' tables (rho_per_split.tsv).
Per type and contrast: mean over seeds (rho_per_contrast.tsv). Per type: mean and SD over the
eight contrasts of those means (summary.tsv), with median cells per donor and DE donors from
../cells_per_donor.tsv (full data) and the median number of donors per half.
"""
import os
import sys

import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CONTRASTS, OUT, REPO, SEEDS, WORK, jobs  # noqa: E402

MIN_GENES = 10


def main():
    rows, missing = [], []
    for seed in SEEDS:
        for g, sch, ct in jobs():
            f = {h: f"{WORK}/fits/{seed}/{h}/{sch}__{ct}.tsv.gz" for h in "AB"}
            if not all(os.path.exists(p) for p in f.values()):
                missing.append((seed, sch, ct))
                continue
            a, b = (pd.read_csv(f[h], sep="\t", index_col=0) for h in "AB")
            for contrast, short in CONTRASTS.items():
                r = {"seed": seed, "group": g, "scheme": sch, "celltype": ct, "contrast": short,
                     "n_donors_A": int(a["n_samples"].iloc[0]), "n_donors_B": int(b["n_samples"].iloc[0])}
                if contrast in a and contrast in b:
                    sh = a.index.intersection(b.index)
                    x, y = a.loc[sh, contrast], b.loc[sh, contrast]
                    ok = x.notna() & y.notna()
                    r["n_genes"] = int(ok.sum())
                    if ok.sum() >= MIN_GENES:
                        r["rho_t"] = spearmanr(x[ok], y[ok]).statistic
                rows.append(r)
    if missing:
        print(f"missing fits: {len(missing)} (seed, scheme, type), e.g. {missing[:3]}")
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/rho_per_split.tsv", sep="\t", index=False, float_format="%.4g")
    keys = ["group", "scheme", "celltype"]
    pc = df.groupby(keys + ["contrast"], sort=False).agg(
        mean_rho_t=("rho_t", "mean"), sd_rho_t_over_splits=("rho_t", "std"), n_splits=("rho_t", "count"),
        median_n_genes=("n_genes", "median")).reset_index()
    pc.to_csv(f"{OUT}/rho_per_contrast.tsv", sep="\t", index=False, float_format="%.4g")
    s = pc.groupby(keys, sort=False).agg(
        mean_rho_t=("mean_rho_t", "mean"), sd_rho_t=("mean_rho_t", "std"), n_contrasts=("mean_rho_t", "count"),
        n_splits=("n_splits", "max")).reset_index()
    nd = df.assign(n=(df["n_donors_A"] + df["n_donors_B"]) / 2).groupby(keys)["n"].median().rename("median_donors_per_half")
    cpd = pd.read_csv(f"{REPO}/subtype_vs_pooled/cells_per_donor.tsv", sep="\t")
    s = s.merge(nd.reset_index(), on=keys).merge(cpd, on=keys, how="left")
    s.to_csv(f"{OUT}/summary.tsv", sep="\t", index=False, float_format="%.4g")
    print(s.to_string())
    print(pc.pivot_table(index=["scheme", "celltype"], columns="contrast", values="mean_rho_t").round(3).to_string())


if __name__ == "__main__":
    main()
