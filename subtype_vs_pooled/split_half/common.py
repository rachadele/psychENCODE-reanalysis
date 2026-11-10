"""Shared paths and type lists for the split-half analysis (same sources as ../subtype_vs_pooled.py)."""
REPO = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis"
R = f"{REPO}/results_v2.1.2_limma"
M = f"{R}/merged_author_2026-11-09"
OUT = f"{REPO}/subtype_vs_pooled/split_half"
WORK = "/space/scratch/rschwartz_psyrerun_work/split_half"
META_DIR = f"{REPO}/gemma/metadata"
DE = {"author": f"{R}/author_submitted_true/subclass/limma_voom/gemma",
      "merged": f"{M}/author_submitted_true/subclass/limma_voom/gemma",
      "ours": f"{R}/author_submitted_false/subclass/limma_voom/gemma"}
# ct_pseudobulks directories the existing fits used (see each fit's .command.sh)
PB = {"author": f"{R}/pseudobulks/author/ct_pseudobulks",
      "merged": f"{M}/pseudobulks/author_merged/ct_pseudobulks",
      "ours": f"{R}/pseudobulks/pavlab_subclass/ct_pseudobulks"}
CONTRASTS = {"Disorder_Schizophrenia_vs_Control": "SCZ", "Disorder_Bipolar.Disorder_vs_Control": "BD",
             "Disorder_ASD_vs_Control": "ASD", "Disorder_MDD_vs_Control": "MDD",
             "Disorder_PTSD_vs_Control": "PTSD", "Disorder_Williams.Syndrome_vs_Control": "WS",
             "Age_death": "Age", "Biological_Sex_male_vs_female": "Sex"}
GROUPS = {"IT": ("L2.3.6.intratelencephalic.projecting.glutamatergic.neuron",
                 ["L2.3_IT", "L4_IT", "L5_IT", "L6_IT", "L6_IT_Car3"]),
          "Lamp5": ("lamp5.GABAergic.cortical.interneuron", ["Lamp5", "Lamp5_Lhx6"]),
          "Sst": ("sst.GABAergic.cortical.interneuron", ["Sst", "Sst_Chodl"])}
SEEDS = list(range(1, 11))


def de_path(scheme, ct, contrast):
    return f"{DE[scheme]}/{ct}/{contrast}/{contrast}_{ct}_results.tsv"


def jobs():
    """(group, scheme, celltype) for the 15 fits per half."""
    out = []
    for g, (label, subs) in GROUPS.items():
        out += [(g, "author", s) for s in subs] + [(g, "merged", label), (g, "ours", label)]
    return out
