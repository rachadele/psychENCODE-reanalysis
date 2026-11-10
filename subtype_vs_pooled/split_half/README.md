# Split-half reproducibility of subtype and pooled DE

We asked if pooling author subtypes gave more stable DE. `se_ratio.py` computed median standard error ratios (type over pooled author type) from the existing limma-voom fits. `make_splits.py` split the 307 donors into halves (seeds 1 to 10), pairing donors by age within cohort × Disorder × Sex strata. `fit_half.py` ran `bin/limma_voom_analysis.R` unchanged on each half; with all donors it reproduced the existing t exactly. `summarize_split_half.py` computed Spearman ρ of t between halves, averaged over seeds, then contrasts (`summary.tsv`).

Subtype standard errors exceeded pooled ones (ratios 1.00 to 1.85). Pooled IT (ρ = 0.51) beat every IT subtype (0.39 to 0.50). Pooled Lamp5 and Sst matched their largest subtype (0.36, 0.41) and exceeded Lamp5 Lhx6 (0.27) and Sst Chodl (0.10).

`split_half_t.png`: agreement between donor halves (A) and relative standard error (B); boxes over the eight contrasts, dots for single contrasts.

`split_half_example_scatter.png`: gene t in each donor half for bipolar disorder, split 9, chosen as closest to the mean pooled-minus-subtype difference and the 10-split mean ρ.

`manifest.json` holds sha256, data version, seeds and timestamps.
