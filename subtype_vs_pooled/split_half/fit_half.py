"""One limma-voom fit on one donor half: python3 fit_half.py <scheme> <celltype> <seed> <half>.

half A or B uses splits.tsv; half "all" keeps every column (check that this wrapper reproduces the
existing fit). The pseudobulk matrix and gene_ncells matrix of the existing fit are restricted to
the half's donors and passed, with the same cell_counts file and the pipeline parameters
(min_cells_per_sample 0, min_cells_per_gene 3, min_variance 0.01; nextflow.config), to
bin/limma_voom_analysis.R unchanged. Moderated t of the eight Figure 4 contrasts is kept in
WORK/fits/<seed>/<half>/<scheme>__<celltype>.tsv.gz; the run directory is then removed.
Resumable: an existing output is skipped.
"""
import os
import shutil
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CONTRASTS, META_DIR, OUT, PB, REPO, WORK  # noqa: E402

RSCRIPT = "/home/rschwartz/anaconda3/envs/r4.3/bin/Rscript"


def main(scheme, ct, seed, half):
    out = f"{WORK}/fits/{seed}/{half}/{scheme}__{ct}.tsv.gz"
    if os.path.exists(out):
        return
    run = f"{WORK}/run/{seed}_{half}_{scheme}__{ct}"
    shutil.rmtree(run, ignore_errors=True)
    os.makedirs(run)
    src = f"{PB[scheme]}/{ct}/{ct}"
    mats = {}
    if half == "all":
        keep = None
    else:
        sp = pd.read_csv(f"{OUT}/splits.tsv", sep="\t", dtype=str)
        keep = set(sp.loc[(sp["seed"] == str(seed)) & (sp["half"] == half), "donor"])
    for kind in ("pseudobulk_matrix", "gene_ncells"):
        m = pd.read_csv(f"{src}_{kind}.tsv.gz", sep="\t", dtype=str)
        if keep is not None:
            m = m[["feature_name"] + [c for c in m.columns[1:] if c in keep]]
        mats[kind] = f"{run}/{ct}_{kind}.tsv"
        m.to_csv(mats[kind], sep="\t", index=False)
    cmd = [RSCRIPT, f"{REPO}/bin/limma_voom_analysis.R", "--pseudobulk_matrix", mats["pseudobulk_matrix"],
           "--metadata", META_DIR, "--cell_type", ct, "--min_cells_per_sample", "0",
           "--min_cells_per_gene", "3", "--min_variance", "0.01",
           "--cell_counts", f"{PB[scheme]}/cell_counts.tsv", "--gene_ncells", mats["gene_ncells"]]
    with open(f"{run}/limma.log", "w") as log:
        env = dict(os.environ, PATH=f"{os.path.dirname(RSCRIPT)}:{os.environ['PATH']}")
        subprocess.run(cmd, cwd=run, stdout=log, stderr=subprocess.STDOUT, check=True, env=env)
    n_samples = next(l for l in open(f"{run}/limma.log") if l.startswith("Matching samples")).split("(")[1].split(")")[0]
    cols = []
    for contrast in CONTRASTS:
        p = f"{run}/{contrast}/{contrast}_{ct}_results.tsv"
        if os.path.exists(p):  # a level absent from the half drops its coefficient
            cols.append(pd.read_csv(p, sep="\t", usecols=["gene", "t"]).set_index("gene")["t"].rename(contrast))
    res = pd.concat(cols, axis=1)
    res.index.name = "gene"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + ".tmp.gz"
    res.assign(n_samples=int(n_samples)).to_csv(tmp, sep="\t")
    os.replace(tmp, out)
    shutil.rmtree(run)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
