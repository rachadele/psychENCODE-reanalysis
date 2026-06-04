# DE Comparison Summary: Class vs Subclass
**Generated:** 2026-06-03
**Annotation levels:** class + subclass (combined)
**Results directories:** `results/comparison/class/`, `results/comparison/subclass/`

This file compares DE concordance (GEMMA vs author-submitted annotations) across both
annotation levels in one place, replacing the two per-level `summary.md` files.

---

## Overall Concordance — Subclass > Class

| Metric | Class mean ρ (SD) | Subclass mean ρ (SD) | Δ (subclass − class) |
|--------|-------------------|----------------------|----------------------|
| Spearman log2FC | 0.761 (0.218) | 0.807 (0.209) | **+0.046** |
| Spearman p-value | 0.657 (0.263) | 0.724 (0.247) | **+0.067** |

Subclass annotation gives meaningfully higher concordance with the author-submitted DE
results, and the SD shrinks slightly — concordance is both higher and more uniform.

**The entire improvement comes from a small number of cell types.** Class and subclass use
the *same* GEMMA pseudobulks for every cell type whose GEMMA-class mapping is unchanged;
their per-cell-type ρ are identical to ~3 decimals. The overall gain is driven exclusively
by the cell types that escape a **many-to-one collapse** at subclass level — i.e. author
subtypes that share one coarse GEMMA class at `class` level but get their own dedicated
GEMMA class at `subclass` level.

### Where the gain comes from

| Author CT | Class GEMMA target | Subclass GEMMA target | Class ρ | Subclass ρ | Δ log2FC |
|-----------|--------------------|-----------------------|---------|------------|----------|
| L5.6_NP | deep.layer.non.IT | near.projecting | 0.697 | 0.975 | **+0.278** |
| L5_ET | deep.layer.non.IT | L5.extratelencephalic | 0.472 | 0.674 | **+0.202** |
| L6_CT | deep.layer.non.IT | corticothalamic | 0.698 | 0.843 | **+0.145** |
| L6b | deep.layer.non.IT | L6b.glutamatergic | 0.752 | 0.896 | **+0.144** |
| VLMC | vascular | vascular.leptomeningeal.cell | 0.367 | 0.475 | **+0.108** |
| Endo | vascular | endothelial.cell | 0.717 | 0.806 | **+0.089** |

The two coarse GEMMA classes that fragment at subclass level — `deep.layer.non.IT`
(→ NP / ET / CT / L6b) and `vascular` (→ endothelial / VLMC / pericyte / SMC) — were the two
worst-performing classes at class level (aggregate ρ 0.655 and 0.605). Splitting them so each
GEMMA pseudobulk matches a single author subtype is what lifts the overall average.

**Caveat:** two vascular author CTs *drop out* of the subclass comparison rather than improve:
`PC` (pericyte, GEMMA F1 = 0.0 — classifier recalls zero pericytes) and `SMC` (no benchmark
entry). At class level they were folded into `vascular` and contributed ρ ≈ 0.79 / 0.55; at
subclass level they have no comparable GEMMA pseudobulk and are absent from the heatmap. So
part of the subclass "improvement" on vascular is also survivorship — the unrecoverable types
simply leave the comparison.

---

## Per-Cell-Type Concordance (both levels)

Sorted by subclass log2FC ρ ascending (worst first). Δ is subclass − class log2FC ρ.
`—` = not present in that level's heatmap.

| Author CT | Class log2FC ρ | Subclass log2FC ρ | Δ | Class p-val ρ | Subclass p-val ρ | Notes |
|-----------|---------------|-------------------|------|--------------|-----------------|-------|
| Immune | 0.371 | 0.371 | 0.00 | 0.313 | 0.297 | Heterogeneous label → `microglial.cell`; worst at both levels |
| Sst_Chodl | 0.376 | 0.376 | 0.00 | 0.120 | 0.200 | Rare subtype; many-to-one collapse with `Sst` at both levels |
| VLMC | 0.367 | 0.475 | +0.11 | 0.224 | 0.377 | Own class at subclass but low F1 (0.819) |
| L6_IT_Car3 | 0.539 | 0.541 | 0.00 | 0.289 | 0.352 | IT collapse persists at both levels |
| L5_ET | 0.472 | 0.674 | +0.20 | 0.316 | 0.500 | Own class at subclass; low F1 (0.678) still limits it |
| L4_IT | 0.710 | 0.711 | 0.00 | 0.546 | 0.595 | IT collapse persists |
| L6_IT | 0.773 | 0.774 | 0.00 | 0.608 | 0.678 | IT collapse; Sex contrast ρ ≈ 0.42 |
| Lamp5_Lhx6 | 0.788 | 0.788 | 0.00 | 0.645 | 0.626 | Lamp5 collapse persists |
| L5_IT | 0.797 | 0.797 | 0.00 | 0.669 | 0.695 | IT collapse |
| Endo | 0.717 | 0.806 | +0.09 | 0.565 | 0.694 | Own `endothelial.cell` class at subclass |
| L6_CT | 0.698 | 0.843 | +0.15 | 0.645 | 0.768 | Own `corticothalamic` class at subclass |
| Pax6 | 0.854 | 0.854 | 0.00 | 0.729 | 0.735 | Low F1 (0.862); transitional interneuron |
| Sncg | 0.865 | 0.865 | 0.00 | 0.748 | 0.754 | Low F1 (0.753); small population |
| Lamp5 | 0.877 | 0.877 | 0.00 | 0.742 | 0.781 | Many-to-one with Lamp5_Lhx6 |
| L6b | 0.752 | 0.896 | +0.14 | 0.673 | 0.836 | Own `L6b.glutamatergic` class at subclass |
| Micro | 0.904 | 0.904 | 0.00 | 0.818 | 0.811 | Good; contrast with `Immune` |
| L2.3_IT | 0.905 | 0.904 | 0.00 | 0.824 | 0.849 | Best IT subtype; most representative of GEMMA class |
| Oligo | 0.942 | 0.942 | 0.00 | 0.913 | 0.936 | High concordance |
| Chandelier | 0.955 | 0.955 | 0.00 | 0.925 | 0.921 | High concordance |
| Sst | 0.957 | 0.957 | 0.00 | 0.870 | 0.916 | Good; contrast with Sst_Chodl |
| Pvalb | 0.965 | 0.965 | 0.00 | 0.945 | 0.950 | High concordance |
| Astro | 0.967 | 0.967 | 0.00 | 0.940 | 0.952 | Very high concordance |
| L5.6_NP | 0.697 | 0.975 | +0.28 | 0.632 | 0.950 | Largest gain; own `near.projecting` class at subclass |
| Vip | 0.976 | 0.976 | 0.00 | 0.937 | 0.950 | Very high concordance |
| OPC | 0.982 | 0.982 | 0.00 | 0.951 | 0.972 | Highest concordance overall |
| PC | 0.787 | — | — | 0.624 | — | Class-only (in `vascular`); pericyte F1 = 0.0, dropped at subclass |
| SMC | 0.548 | — | — | 0.313 | — | Class-only (in `vascular`); no F1 benchmark, dropped at subclass |

---

## Mechanisms

### Resolved by going to subclass (many-to-one collapse broken)
`deep.layer.non.IT` and `vascular` each averaged several distinct author subtypes into one
GEMMA pseudobulk at class level. At subclass level each subtype gets its own GEMMA class, so
the GEMMA and author pseudobulks describe the same population — hence the +0.09 to +0.28 gains
for L5.6_NP, L5_ET, L6_CT, L6b, Endo, VLMC.

### Unchanged at subclass (collapse persists)
Three many-to-one mappings survive at subclass level, so these cell types show **no** gain:
- **IT neurons** — `L2.3_IT`, `L4_IT`, `L5_IT`, `L6_IT`, `L6_IT_Car3` still all map to one
  `L2.3.6.intratelencephalic` class (n = 40 pairs, mean ρ ≈ 0.745). `L6_IT_Car3` (rare, Car3+)
  is most diluted (ρ ≈ 0.54).
- **SST subtypes** — `Sst`, `Sst_Chodl` → `sst.GABAergic`; `Sst_Chodl` stays at ρ ≈ 0.38.
- **Lamp5 subtypes** — `Lamp5`, `Lamp5_Lhx6` → `lamp5.GABAergic`.

### Unchanged at subclass (classifier-quality / label-quality limited)
Going finer cannot help when the limit is the GEMMA classifier or the author label itself:
- **`Immune`** — heterogeneous catch-all mapped to `microglial.cell`; worst at both levels
  (ρ = 0.371), near-zero/negative for PTSD (ρ = −0.19). `Micro` → same GEMMA class scores
  ρ = 0.904, isolating the label (not the mapping) as the problem.
- **`L5_ET`** (F1 0.678) and **`VLMC`** (F1 0.819) — get their own class at subclass but low
  classifier recall keeps the pseudobulks compositionally divergent.
- **`Sncg`** (F1 0.753), **`Pax6`** (F1 0.862) — below-threshold F1, unchanged across levels.
- **`PC`** (F1 0.0) — pericytes are not recalled at all; comparison impossible at subclass.

---

## Contrast-Specific Observations (both levels)

**Sex** is the most discordant contrast at both levels, concentrated in deep-layer and IT
glutamatergic subtypes (L6_IT ρ ≈ 0.42, L6_CT ≈ 0.52, L6b ≈ 0.62, L4_IT ≈ 0.61). Sex-DE is
narrowly layer/subtype-specific in cortex, so small compositional differences between schemes
have outsized effects. Treat Sex-DE as preliminary at either level.

**WS (Williams Syndrome)** uniquely collapses for `L5_ET` (class ρ = 0.01, subclass ρ = 0.38)
— few cohorts, small cell numbers, amplified stochasticity.

**PTSD** drives the single most extreme value at both levels: `Immune` ρ = −0.19; `Endo` also
drops (~0.69) consistent with annotation-sensitive vascular involvement.

**MDD** has sparse significant DE (many empty Jaccard cells) at both levels — overlap metrics
are unreliable for MDD regardless of resolution.

---

## Recommendations

### Use subclass-level results for
- **L5.6_NP, L6b, L6_CT, Endo** — now well-matched (ρ 0.81–0.98); these were unreliable at
  class level and should be read at subclass.
- Any deep-layer or vascular interpretation generally — class level over-aggregates them.

### Exclude from downstream interpretation (both levels)
- **`Immune`** — heterogeneous; use `Micro` for microglial conclusions.
- **`Sst_Chodl`** — rare, fully diluted in the SST class; DE is essentially noise.
- **`PC`** (pericyte F1 = 0.0) and **`SMC`** (no benchmark) — not recoverable by GEMMA.

### Interpret with caution (both levels)
- **`L5_ET`** — low F1 (0.678); WS and Sex unreliable even at subclass.
- **`VLMC`** — low F1 (0.819), low cell numbers; Jaccard near zero in most contrasts.
- **`L6_IT_Car3`** — rare Car3 subtype still diluted in the broad IT class.
- **All cell types for Sex** and **all cell types for MDD** — systematic at both levels.

### Will not improve by changing annotation level
- **IT subtypes** (`L4_IT`, `L5_IT`, `L6_IT`, `L6_IT_Car3`) and **Lamp5_Lhx6** stay collapsed
  at subclass; resolving them needs a finer GEMMA IT/Lamp5 partition, not subclass mapping.
- Classifier-limited types (`L5_ET`, `Sncg`, `Pax6`, `PC`) need classifier improvement, not
  finer labels.
