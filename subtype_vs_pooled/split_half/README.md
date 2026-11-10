# Split-half reproducibility of subtype and pooled DE

We asked whether pooling author subtypes gave more stable DE than DE on each subtype. `se_ratio.py` computed, from the existing limma-voom fits, the median standard error ratio (type over pooled author type) per contrast. `make_splits.py` split the 307 donors into halves for seeds 1 to 10, pairing donors by age within cohort × Disorder × Sex strata. `fit_half.py` ran `bin/limma_voom_analysis.R` unchanged on each half (all donors reproduced the existing t exactly). `summarize_split_half.py` computed Spearman ρ of t between halves, averaged over seeds per contrast, then over contrasts (`summary.tsv`). `run_split_half.sh` ran all fits (tmux `split_half_de`).

Subtype standard errors exceeded pooled ones (ratios 1.00 to 1.85). Pooled IT (ρ = 0.51) was more reproducible than every IT subtype (0.39 to 0.50). Pooled Lamp5 and Sst matched their largest subtype (0.36, 0.41) and exceeded Lamp5 Lhx6 (0.27) and Sst Chodl (0.10).

`split_half_t.png` (`plot_split_half_t.py`): agreement of t between donor halves (A) and standard error relative to the pooled author type (B) for each subtype, ordered by cells per donor, the pooled author type and our label; one dot per contrast (mean over 10 splits), bars at the mean.

`manifest.json` holds sha256, data version, seeds and timestamps.
