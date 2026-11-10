"""Donor split halves, one per seed (SEEDS in common.py).

Universe: every donor with metadata (Individual_ID in gemma/metadata/*_sample_meta.tsv), the donors
limma_voom_analysis.R can use. Each donor gets the same half for every cell type.
Stratification: one split serves all eight contrasts. Within each cohort x Disorder x Sex stratum,
donors are sorted by age at death and taken in consecutive pairs; one donor of each pair goes to
each half at random, and an odd last donor goes to a random half. The halves are therefore balanced
on cohort, on each Disorder group (six disorder contrasts), on sex (Sex contrast) and on age
(Age contrast). Output: splits.tsv (seed, donor, cohort, disorder, sex, age, half).
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import META_DIR, OUT, SEEDS  # noqa: E402


def donors():
    d = pd.concat([pd.read_csv(f, sep="\t", dtype=str).assign(cohort=os.path.basename(f).replace("_sample_meta.tsv", ""))
                   for f in sorted(glob.glob(f"{META_DIR}/*.tsv"))])
    d["Disorder"] = d["Disorder"].replace("control", "Control")
    d["age"] = pd.to_numeric(d["Age_death"].str.replace("+", "", regex=False))
    return d.rename(columns={"Individual_ID": "donor", "Disorder": "disorder", "Biological_Sex": "sex"})[
        ["donor", "cohort", "disorder", "sex", "age"]].sort_values("donor").reset_index(drop=True)


def split(d, seed):
    rng = np.random.default_rng(seed)
    half = {}
    for _, s in d.groupby(["cohort", "disorder", "sex"], sort=True):
        ids = s.sort_values(["age", "donor"])["donor"].tolist()
        for i in range(0, len(ids) - 1, 2):
            a, b = (ids[i], ids[i + 1]) if rng.random() < 0.5 else (ids[i + 1], ids[i])
            half[a], half[b] = "A", "B"
        if len(ids) % 2:
            half[ids[-1]] = "A" if rng.random() < 0.5 else "B"
    return d.assign(seed=seed, half=d["donor"].map(half))


def main():
    d = donors()
    assert d["donor"].is_unique
    out = pd.concat([split(d, s) for s in SEEDS])
    out[["seed", "donor", "cohort", "disorder", "sex", "age", "half"]].to_csv(f"{OUT}/splits.tsv", sep="\t", index=False)
    print(out.groupby(["seed", "half"]).size().unstack())


if __name__ == "__main__":
    main()
