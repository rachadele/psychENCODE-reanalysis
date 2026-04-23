# DE Comparison Summary: Subclass Level
**Generated:** 2026-04-08  
**Annotation level:** subclass  
**Results directory:** `results/comparison/subclass/`

---

## Overall Concordance

| Metric | Mean ρ | SD |
|--------|--------|----|
| Spearman log2FC | 0.807 | 0.209 |
| Spearman p-value | 0.724 | 0.247 |

Overall concordance is meaningfully better at subclass level than class level (class: log2FC ρ = 0.761, p-value ρ = 0.657). The improvement reflects more precise GEMMA class assignments that better match individual author subtypes, particularly for deep-layer glutamatergic neurons and vascular cells that suffered from many-to-one collapse at class level.

---

## Per-Cell-Type Concordance

Sorted by mean log2FC ρ ascending (worst first). F1 scores are at the subclass level from the GEMMA sc-pipeline benchmark. `—` = no F1 entry at this level.

| Author CT | GEMMA class | F1 (subclass) | Mean log2FC ρ | Mean p-val ρ | Notes |
|-----------|-------------|--------------|--------------|-------------|-------|
| Immune | microglial.cell | 0.978† | 0.371 | 0.297 | †F1 applies to `Micro`; `Immune` is heterogeneous |
| Sst_Chodl | sst.GABAergic | 0.932 | 0.376 | 0.200 | Rare subtype; many-to-one collapse with Sst |
| VLMC | vascular.leptomeningeal.cell | 0.819 | 0.475 | 0.377 | Own class at subclass but low F1; rare cell type |
| L6_IT_Car3 | L2.3.6.IT | 0.991 | 0.541 | 0.352 | Rare Car3 subtype; many-to-one IT collapse |
| L5_ET | L5.extratelencephalic | 0.678 | 0.674 | 0.500 | Own class at subclass; low F1 still limits concordance |
| L4_IT | L2.3.6.IT | 0.991 | 0.711 | 0.595 | Many-to-one IT collapse; Sex/ASD contrasts drop below 0.70 |
| Lamp5_Lhx6 | lamp5.GABAergic | 0.917 | 0.788 | 0.626 | Many-to-one with Lamp5; Sex contrast drops to 0.66 |
| L5_IT | L2.3.6.IT | 0.991 | 0.797 | 0.695 | Many-to-one IT collapse; Sex drops to 0.68 |
| Endo | endothelial.cell | 0.826 | 0.806 | 0.694 | Own class at subclass (improved from class); PTSD/Sex lower |
| L6_IT | L2.3.6.IT | 0.991 | 0.774 | 0.678 | Many-to-one IT; Sex drops to 0.42 |
| L2.3_IT | L2.3.6.IT | 0.991 | 0.904 | 0.849 | Best of IT subtypes; most representative of GEMMA class |
| Micro | microglial.cell | 0.978 | 0.904 | 0.811 | Good concordance; contrast with `Immune` |
| Pax6 | PAX6.GABAergic | 0.862 | 0.854 | 0.735 | Low F1; WS drops to 0.69 |
| Sncg | sncg.GABAergic | 0.753 | 0.865 | 0.754 | Low F1 but reasonable concordance |
| Lamp5 | lamp5.GABAergic | 0.917 | 0.877 | 0.781 | Moderate; many-to-one with Lamp5_Lhx6 |
| L6_CT | corticothalamic | 0.849 | 0.843 | 0.768 | Own class at subclass (improved from 0.698); Sex drops to 0.52 |
| L6b | L6b.glutamatergic | 0.843 | 0.896 | 0.836 | Own class at subclass (improved from 0.752); Sex drops to 0.62 |
| Oligo | oligodendrocyte | 0.993 | 0.942 | 0.936 | High concordance; Sex contrast slightly lower (0.81) |
| Chandelier | chandelier.pvalb | 0.967 | 0.955 | 0.921 | High concordance |
| Pvalb | pvalb.GABAergic | 0.899 | 0.965 | 0.950 | High concordance |
| Astro | astrocyte | 0.993 | 0.967 | 0.952 | Very high concordance |
| Sst | sst.GABAergic | 0.932 | 0.957 | 0.916 | Good; contrast with Sst_Chodl |
| L5.6_NP | near.projecting | 0.957 | 0.975 | 0.950 | Large improvement vs class (0.697→0.975); own class now |
| Vip | vip.GABAergic | 0.944 | 0.976 | 0.950 | Very high concordance |
| OPC | oligodendrocyte.precursor.cell | 0.995 | 0.982 | 0.972 | Highest concordance overall |

**Not present in heatmap (absent from pairwise comparisons):**
- `PC` → `pericyte` (F1 = 0.0 — classifier completely fails to recall pericytes; likely absent from GEMMA pseudobulks)
- `SMC` → `smooth.muscle.cell` (no F1 benchmark entry; rare cell type not reliably annotated)

---

## Key Improvements Over Class Level

The subclass mapping resolves several many-to-one collapses that hurt class-level concordance:

| Author CT | Class ρ | Subclass ρ | Reason for improvement |
|-----------|---------|-----------|----------------------|
| L5.6_NP | 0.697 | 0.975 | Gets own class (`near.projecting`) at subclass |
| L6b | 0.752 | 0.896 | Gets own class (`L6b.glutamatergic`) at subclass |
| L6_CT | 0.698 | 0.843 | Gets own class (`corticothalamic`) at subclass |
| L5_ET | 0.472 | 0.674 | Gets own class, but low F1 (0.678) still limits concordance |
| Endo | 0.717 | 0.806 | Gets own class (`endothelial.cell`) at subclass |

---

## Mechanisms of Low Concordance

### 1. No benchmark match / heterogeneous author label
**`Immune`** maps to `microglial.cell` but was not evaluated as a distinct benchmark class. Its concordance is the worst in the dataset (mean ρ = 0.371) and highly variable across contrasts — near-zero or negative for PTSD (ρ = −0.19) and Age (ρ = 0.20). `Micro`, mapping to the same GEMMA class with F1 = 0.978, has mean ρ = 0.904. The contrast between these two author labels using the same GEMMA target is the clearest evidence that `Immune` is compositionally heterogeneous and incompatible with `microglial.cell`.

### 2. Low classifier F1 → inconsistent cell assignment
**`L5_ET`** (F1 = 0.678) is the most affected. Despite having its own GEMMA class at subclass level, concordance remains poor (mean ρ = 0.674), driven by near-zero agreement for WS (ρ = 0.38) and Sex (ρ = 0.22). The classifier's low recall (0.67) means many true L5 ET cells are assigned to other classes, producing divergent pseudobulk compositions. This is a classifier quality problem that subclass mapping alone cannot fix.

**`VLMC`** (F1 = 0.819) now has its own `vascular.leptomeningeal.cell` class at subclass level but concordance is still poor (mean ρ = 0.475). The low F1 and likely small cell numbers across cohorts combine to produce unreliable pseudobulks.

**`Sncg`** (F1 = 0.753) and **`Pax6`** (F1 = 0.862) have below-threshold F1 scores with correspondingly reduced concordance.

**`PC` (pericyte)** has F1 = 0.0 with recall = 0.0 — the GEMMA classifier does not successfully annotate any pericytes. This cell type is effectively absent from GEMMA pseudobulks and cannot be compared. It does not appear in the heatmap.

### 3. Many-to-one collapse — fine author subtypes into coarse GEMMA classes
Three many-to-one mappings remain at subclass level:

- **IT neurons**: `L2.3_IT`, `L4_IT`, `L5_IT`, `L6_IT`, `L6_IT_Car3` → all `L2.3.6.intratelencephalic.projecting.glutamatergic.neuron` (n=40 pairs, mean ρ = 0.745). Layer-specific biology is lost. `L6_IT_Car3` is the most affected (mean ρ = 0.541) as the rarest and most transcriptionally distinct subtype.

- **SST subtypes**: `Sst`, `Sst_Chodl` → both `sst.GABAergic` (n=16 pairs, mean ρ = 0.666). `Sst_Chodl` (mean ρ = 0.376) is a rare, transcriptionally extreme subtype whose signal is diluted by collapsing with the much larger `Sst` population.

- **Lamp5 subtypes**: `Lamp5`, `Lamp5_Lhx6` → both `lamp5.GABAergic` (n=16 pairs, mean ρ = 0.832). Less severe than SST because the two subtypes are more similar, but `Lamp5_Lhx6` still underperforms `Lamp5`.

---

## Contrast-Specific Observations

**Sex** drives the most widespread low concordance, affecting nearly all deep-layer and IT glutamatergic subtypes:

| Cell Type | Sex ρ (log2FC) | Pattern |
|-----------|---------------|---------|
| L6_IT | 0.419 | Contrast-driven; other contrasts ≥ 0.74 |
| L6_CT | 0.516 | Contrast-driven |
| L6b | 0.615 | Contrast-driven |
| L4_IT | 0.609 | Contrast-driven |
| Lamp5_Lhx6 | 0.659 | Contrast-driven |
| L5_IT | 0.677 | Contrast-driven |
| Endo | 0.626 | Contrast-driven |
| L6_IT_Car3 | 0.338 | Compound (cell-type + contrast) |

Sex-DE is highly layer- and subtype-specific in cortex. Small compositional differences between annotation schemes produce outsized DE discordance when the true signal is concentrated in a narrow cell population.

**WS (Williams Syndrome)** uniquely breaks down for `L5_ET` (ρ = 0.38). WS is represented in few cohorts with small sample sizes, amplifying stochastic differences between annotation schemes.

**PTSD** shows the most extreme single value: `Immune` ρ = −0.19. Also `Endo` drops to 0.69 in PTSD, consistent with vascular involvement in PTSD biology being annotation-sensitive.

**MDD** has sparse Jaccard overlap across most cell types (many empty cells in the Jaccard heatmap), suggesting insufficient significant DE genes — concordance metrics may be unreliable for MDD at this resolution.

---

## Recommendations

### Exclude from downstream interpretation
- **`Immune`** — heterogeneous catch-all; use `Micro` for microglial conclusions
- **`Sst_Chodl`** — rare subtype fully diluted in GEMMA's SST class; DE results are essentially noise
- **`PC`** — pericyte F1 = 0.0; GEMMA cannot annotate this cell type at all
- **`SMC`** — no F1 benchmark; smooth muscle cells are not reliably recoverable

### Interpret with caution
- **`L5_ET`** — low F1 (0.678); WS and Sex contrasts unreliable even at subclass level
- **`VLMC`** — low F1 (0.819) and low cell numbers; Jaccard near zero for most contrasts
- **`L6_IT_Car3`** — rare Car3 subtype still diluted in the broad IT class
- **`Immune`** for any contrast — but especially PTSD and Age
- **All cell types for Sex contrast** — systematic annotation-scheme sensitivity
- **All cell types for MDD** — insufficient DE signal; overlap metrics unreliable

### No longer concerns at subclass level (improved from class)
- **`L5.6_NP`** — now well-matched (ρ = 0.975); safe to interpret
- **`L6b`** — good concordance (ρ = 0.896) with own GEMMA class
- **`L6_CT`** — reasonable concordance (ρ = 0.843); Sex contrast still lower
- **`Endo`** — adequate concordance (ρ = 0.806) with own endothelial class
