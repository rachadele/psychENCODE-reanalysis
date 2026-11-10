# Split-half reproducibility of subtype and pooled DE

We asked whether pooling author subtypes gave more stable DE than DE on each subtype. First, from the existing limma-voom fits, `se_ratio.py` computed the median ratio of moderated standard errors (subtype over merged author type) per contrast (`se_ratio_summary.tsv`). Second, `make_splits.py` split the 307 donors into halves for ten seeds (1 to 10), pairing donors by age within cohort × Disorder × Sex strata, so that one split served all eight contrasts. `fit_half.py` ran `bin/limma_voom_analysis.R` unchanged on each half for the nine subtypes, three merged author types and three pipeline labels; with all donors it reproduced the existing t exactly. `summarize_split_half.py` computed Spearman ρ of t between halves, averaged over seeds per contrast, then gave mean ± SD over contrasts (`summary.tsv`; per split in `rho_per_split.tsv`). `run_split_half.sh` ran all fits (tmux session `split_half_de`, 75 s per seed).

Subtype standard errors exceeded pooled ones (ratios 1.00 to 1.85). Pooled IT (ρ = 0.51) was more reproducible than every IT subtype (0.39 to 0.50). Pooled Lamp5 and Sst matched their dominant subtype (0.36 and 0.41) and exceeded Lamp5 Lhx6 (0.27) and Sst Chodl (0.10).

`manifest.json` holds sha256, data version, seeds and timestamps.
