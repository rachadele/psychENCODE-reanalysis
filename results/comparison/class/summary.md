# DE Comparison Summary: Class Level
**Generated:** 2026-04-06  
**Annotation level:** class  
**Results directory:** `results/comparison/class/`

---

## Overall Concordance

| Metric | Mean ρ | SD |
|--------|--------|----|
| Spearman log2FC | 0.761 | 0.218 |
| Spearman p-value | 0.657 | 0.263 |

Overall concordance is moderate. The SD is large, indicating substantial variation across cell types and contrasts. P-value concordance is consistently lower than log2FC concordance, as expected — effect sizes are more reproducible than significance thresholds across annotation schemes.

---

## Per-Cell-Type Concordance

Sorted by mean log2FC ρ ascending (worst first). F1 scores are at the class level from the GEMMA sc-pipeline benchmark.

| Author CT | GEMMA Class | F1 (class) | Mean log2FC ρ | Mean p-val ρ | Notes |
|-----------|-------------|-----------|--------------|-------------|-------|
| Immune | microglial.cell | 0.978† | 0.371 | 0.313 | †F1 applies to `Micro`, not `Immune`; heterogeneous label |
| Sst_Chodl | sst.GABAergic | 0.932 | 0.376 | 0.120 | Rare subtype collapsed into coarse SST class |
| VLMC | vascular | 0.979 | 0.367 | 0.224 | Rare; many-to-one with Endo/PC/SMC into `vascular` |
| L5_ET | deep.layer.non.IT | 0.921† | 0.472 | 0.316 | †Subclass F1 = 0.678; poor classifier recall at fine level |
| L6_IT_Car3 | L2.3.6.IT | 0.991 | 0.539 | 0.289 | Rare Car3 subtype diluted into broad IT class |
| SMC | vascular | 0.979 | 0.548 | 0.313 | Rare smooth muscle; many-to-one vascular collapse |
| L5.6_NP | deep.layer.non.IT | 0.921 | 0.697 | 0.632 | Many-to-one with L5_ET/L6_CT/L6b |
| L6_CT | deep.layer.non.IT | 0.921 | 0.698 | 0.645 | Many-to-one deep-layer collapse |
| L4_IT | L2.3.6.IT | 0.991 | 0.710 | 0.546 | Layer-specific subtype collapsed into broad IT class |
| Endo | vascular | 0.979 | 0.717 | 0.565 | Many-to-one vascular collapse |
| L5_IT | L2.3.6.IT | 0.991 | 0.797 | 0.669 | Many-to-one IT collapse |
| Lamp5_Lhx6 | lamp5.GABAergic | 0.917 | 0.788 | 0.645 | Minor Lamp5 subtype; many-to-one with Lamp5 |
| L6_IT | L2.3.6.IT | 0.991 | 0.773 | 0.608 | Many-to-one IT collapse; Sex contrast drops to 0.41 |
| PC | vascular | 0.979 | 0.787 | 0.624 | Many-to-one vascular collapse |
| L6b | deep.layer.non.IT | 0.921 | 0.752 | 0.673 | Many-to-one deep-layer collapse |
| Micro | microglial.cell | 0.978 | 0.904 | 0.818 | Good concordance; contrast with `Immune` |
| L2.3_IT | L2.3.6.IT | 0.991 | 0.905 | 0.824 | Best among IT subtypes; most similar to GEMMA class |
| Lamp5 | lamp5.GABAergic | 0.917 | 0.877 | 0.742 | Moderate; lower than other GABAergic types |
| Sncg | sncg.GABAergic | 0.753 | 0.865 | 0.748 | Low F1 but decent concordance; small population |
| Pax6 | PAX6.GABAergic | 0.862 | 0.854 | 0.729 | Low F1; PAX6 is a transitional/rare interneuron |
| Sst | sst.GABAergic | 0.932 | 0.957 | 0.870 | Good; contrast with Sst_Chodl |
| Pvalb | pvalb.GABAergic | 0.899 | 0.965 | 0.945 | High concordance despite modest F1 |
| Chandelier | chandelier.pvalb | 0.967 | 0.955 | 0.925 | High concordance |
| Oligo | oligodendrocyte | 0.993 | 0.942 | 0.913 | High concordance; Sex drops to 0.81 |
| Astro | astrocyte | 0.993 | 0.967 | 0.940 | Very high concordance |
| Vip | vip.GABAergic | 0.944 | 0.976 | 0.937 | Very high concordance |
| OPC | oligodendrocyte.precursor.cell | 0.995 | 0.982 | 0.951 | Highest concordance overall |

---

## Mechanisms of Low Concordance

### 1. No benchmark match / heterogeneous author label
**`Immune`** is mapped to `microglial.cell` but was not evaluated as a distinct class in the F1 benchmark. Unlike `Micro` (a clean microglia label with F1 = 0.978 and mean ρ = 0.904), `Immune` likely captures a mixture of microglia, infiltrating macrophages, T cells, and other immune populations depending on the study. This compositional heterogeneity varies across cohorts, producing widely variable DE signals. The PTSD contrast is extreme (ρ = −0.19), suggesting that cohort-specific immune infiltration in PTSD samples drives near-random disagreement between the two annotation schemes.

### 2. Low classifier F1 → inconsistent cell assignment
**`Sncg`** (F1 = 0.753) and **`Pax6`** (F1 = 0.862) are rare interneuron subtypes where the GEMMA classifier struggles. Misclassified cells contaminate the pseudobulk, decorrelating DE from the author-defined population.

**`L5_ET`** is the most severe case. The class-level F1 is 0.921, but the subclass-level F1 drops to 0.678 — the lowest in the benchmark. Many L5 ET cells are reassigned to other deep-layer types, making the GEMMA `deep.layer.non.IT` pseudobulk compositionally very different from the author `L5_ET` pseudobulk. This explains the near-zero concordance for WS (ρ = 0.014) and Sex (ρ = 0.291).

### 3. Many-to-one collapse — fine author subtypes into coarse GEMMA classes
Several groups of author CTs all collapse into a single GEMMA class. The GEMMA pseudobulk averages over the entire class, while each author subtype captures a specific fraction:

- **IT neurons**: `L2.3_IT`, `L4_IT`, `L5_IT`, `L6_IT`, `L6_IT_Car3` → all `L2.3.6.intratelencephalic.projecting.glutamatergic.neuron`  
  L2.3_IT has the best concordance because it is the most abundant subtype and most representative of the GEMMA class average. L6_IT_Car3 is a sparse, transcriptionally distinct subtype that contributes minimally to the GEMMA pseudobulk (mean ρ = 0.539).

- **Deep-layer non-IT**: `L5.6_NP`, `L5_ET`, `L6_CT`, `L6b` → all `deep.layer.non.IT`  
  Each subtype has a distinct transcriptional identity but they are averaged together in GEMMA. L5_ET suffers additionally from poor classifier recall.

- **Vascular**: `Endo`, `PC`, `SMC`, `VLMC` → all `vascular`  
  These are rare cell types with distinct biology. VLMC (mean ρ = 0.367) is the worst, as vascular leptomeningeal cells are the most compositionally distinct from the aggregate `vascular` class.

- **SST subtypes**: `Sst`, `Sst_Chodl` → both `sst.GABAergic`  
  `Sst_Chodl` is a rare, transcriptionally extreme SST subtype. Its signal is heavily diluted in the GEMMA `sst.GABAergic` pseudobulk, yielding near-zero concordance (mean ρ = 0.376, mean Jaccard ≈ 3%).

---

## Contrast-Specific Observations

**Sex** is the most problematic contrast across multiple cell types:

| Cell Type | Sex ρ (log2FC) | Pattern |
|-----------|---------------|---------|
| L6_IT | 0.409 | Contrast-driven; other contrasts ≥ 0.73 |
| L6_CT | 0.365 | Contrast-driven |
| L5_ET | 0.291 | Compound (cell-type + contrast) |
| L5.6_NP | 0.537 | Contrast-driven |
| L6b | 0.618 | Contrast-driven |
| Lamp5_Lhx6 | 0.659 | Contrast-driven |
| PC | 0.596 | Contrast-driven |
| Endo | 0.538 | Contrast-driven |

Sex-DE appears more sensitive to annotation scheme than disease contrasts. This is consistent with sex-biased gene expression being highly cell-type-specific and layer-specific in cortex — small compositional differences between annotation schemes have outsized effects when the true signal is concentrated in a subset of cells within a class.

**WS (Williams Syndrome)** has a unique vulnerability in `L5_ET` (ρ = 0.014) — essentially no agreement. WS is represented in only one or two cohorts with small cell numbers, so stochastic sampling differences between annotation schemes are amplified.

**PTSD** shows severe breakdown in `Immune` (ρ = −0.193) and `VLMC` (ρ = 0.186), likely reflecting cohort-specific cellular composition in PTSD samples (immune infiltration, vascular remodeling).

**MDD** shows widespread missing Jaccard values (many cells shown as empty in the overlap heatmap), suggesting insufficient significant DE genes in most cell types to compute meaningful overlap — concordance metrics may be unreliable for MDD at this annotation level.

---

## Jaccard Overlap Summary

High overlap (>70% for BD, Age, Sex contrasts) in well-characterized classes: OPC, Astro, Oligo, Micro, Chandelier, Vip, Sst, Pvalb.

Overlap is consistently near 0% for:
- `Sst_Chodl` (all contrasts): rare subtype with almost no shared DE genes
- `VLMC` (most contrasts): rare vascular type
- `SMC` (BD, WS, SCZ): smooth muscle cells rarely reach significance
- `L5_ET` (WS = 0%, ASD ≈ 1.5%): classifier failure produces non-overlapping gene sets

---

## Recommendations

### Exclude from downstream interpretation
- **`Immune`** — heterogeneous label; GEMMA maps it to `microglial.cell` but the biology is incompatible. Replace with `Micro` for microglial-specific conclusions.
- **`Sst_Chodl`** — sparse, rare subtype fully diluted in GEMMA's SST class; DE results from the two schemes are essentially uncorrelated.
- **`VLMC`** — rare vascular subtype, near-zero Jaccard in most contrasts, low ρ.

### Interpret with caution
- **`L5_ET`** — very low F1 at subclass level; WS and Sex contrasts are unreliable.
- **`L6_IT_Car3`** — rare Car3 subtype; results driven by cell-number stochasticity.
- **`SMC`** — low cell numbers and many-to-one vascular collapse.
- **All cell types for Sex contrast** — annotation-scheme sensitivity is broadly elevated; treat Sex-DE findings as preliminary pending subclass-level validation.
- **All cell types for MDD** — sparse DE signal; overlap metrics are unreliable.

### Investigate at subclass level
- **`L5_ET`** — low subclass F1 (0.678) suggests the classifier itself needs improvement; subclass GEMMA annotations will not resolve this until classifier recall improves.
- **`L4_IT`, `L6_IT`, `L6_IT_Car3`** — layer-specific biology is lost in the broad IT class; subclass results will provide better-matched pseudobulks.
- **`Lamp5_Lhx6`** — distinct from `Lamp5` biologically; subclass annotations separate these correctly.
- **`L6b`, `L6_CT`, `L5.6_NP`** — all collapse into `deep.layer.non.IT` at class level; subclass results should be examined for these individually.
