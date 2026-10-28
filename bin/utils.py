import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np
import pandas as pd
import re
from collections import OrderedDict, Counter

def load_label_map(params_path, label_map_path=None):
    """Author -> pipeline label map: the --label_map JSON if given, else the inline
    gemma_to_gemma_map of the params JSON (backward compatible)."""
    import json
    if label_map_path and label_map_path not in ("null", "None"):
        with open(label_map_path) as f:
            return json.load(f)
    with open(params_path) as f:
        return json.load(f)["gemma_to_gemma_map"]


def reverse_lookup(mapping, value):
    """Return all keys that map to the given value."""
    return [k for k, v in mapping.items() if v == value]

def cell_types_match(row, gemma_to_gemma_map):
    ct1 = row['ct_author']
    ct2 = row['ct_sc_pipeline']
    return ct1 in reverse_lookup(gemma_to_gemma_map, ct2)
     

def plot_boxplot(df, gemma_to_gemma_map, x="percent_overlap", outpath=None):
    plt.figure(figsize=(20, 20))
    # Create a copy with cleaned cell type names (replace "." with " ")
    df_plot = df.copy()
    df_plot['ct_sc_pipeline_clean'] = df_plot['ct_sc_pipeline'].str.replace('.', ' ', regex=False)
    celltype_order = list(OrderedDict.fromkeys(gemma_to_gemma_map.values()))
    celltype_order = [ct for ct in celltype_order if ct in df['ct_sc_pipeline'].unique()]
    celltype_order_clean = [ct.replace('.', ' ') for ct in celltype_order]
    ax = sns.boxplot(y='ct_sc_pipeline_clean', x=x, data=df_plot, hue='ct_sc_pipeline_clean', dodge=False, showmeans=False, order=celltype_order_clean)
    means = df_plot.groupby('ct_sc_pipeline_clean')[x].mean()
    for i, celltype in enumerate(celltype_order_clean):
        if celltype in means:
            ax.scatter(means[celltype], i, color='black', s=100, zorder=10, label='Mean' if i == 0 else "")
    plt.yticks(fontsize=28)
    plt.xticks(fontsize=28)
    plt.xlabel(x.replace('_', ' ').title(), fontsize=32)
    plt.ylabel('Cell Type', fontsize=32)
    plt.title(f'{x.replace("_", " ").title()} Across All Contrasts', fontsize=32)
    handles, labels = ax.get_legend_handles_labels()
    if 'Mean' in labels:
        plt.legend([handles[labels.index('Mean')]], ['Mean'], loc='lower right', fontsize=28)
    else:
        plt.legend([],[], frameon=False)
    plt.tight_layout()
    if outpath is None:
        outpath = f'boxplot_{x}.png'
    plt.savefig(outpath, dpi=150)
    plt.close()

def plot_stripplot(df, gemma_to_gemma_map, x="percent_overlap", output_path=None):
    plt.figure(figsize=(36, 28))
    # Create a copy with cleaned cell type names (replace "." with " ")
    df_plot = df.copy()
    df_plot['ct_sc_pipeline_clean'] = df_plot['ct_sc_pipeline'].str.replace('.', ' ', regex=False)
    celltype_order = list(OrderedDict.fromkeys(df['ct_sc_pipeline'].unique()))
    celltype_order_clean = [ct.replace('.', ' ') for ct in celltype_order]
    ax = sns.stripplot(y='ct_sc_pipeline_clean', x=x, data=df_plot, order=celltype_order_clean, hue='contrast', dodge=False, jitter=True, size=18)
    plt.yticks(fontsize=28)
    plt.xticks(fontsize=28)
    plt.xlabel(x.replace('_', ' ').title(), fontsize=32)
    plt.ylabel('Cell Type', fontsize=32)
    plt.title(f'{x.replace("_", " ").title()} Across Contrasts', fontsize=32)
    plt.legend(title='Contrast', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=24, title_fontsize=28)
    plt.tight_layout()
    if output_path is None:
        output_path = f'stripplot_{x}.png'
    plt.savefig(output_path, dpi=150)
    plt.close()

# Compute average overlap per cell type across all contrasts
def compute_average_per_celltype(combined_df, metric):
    result = combined_df.groupby('ct_sc_pipeline')[metric].agg(['mean', 'std', 'count']).reset_index()
    result.rename(columns={'mean': f'average_{metric}', 'std': f'sd_{metric}', 'count': 'n_pairs'}, inplace=True)
    # sort from lowest to highest average percent overlap
    result = result.sort_values(by=f'average_{metric}', ascending=False)
    return result

def compute_overall_averages(combined_df, metric):
	# Calculate average percent overlap for matched cell types
	average_overlap = combined_df[metric].mean()
	sd_overlap = combined_df[metric].std()
	# write overall average to file
	return average_overlap, sd_overlap


# Pipeline label -> benchmark label, per level. Old and new pipeline names both listed.
PIPELINE_TO_BENCH = {
    "class": {
        "astrocyte": "Astrocyte",
        "chandelier.pvalb.GABAergic.cortical.interneuron": "Chandelier",
        "vascular": "Vascular",
        "brain.vascular.cell": "Vascular",
        "immune": "Immune",
        "microglial.cell": "Microglia",
        "oligodendrocyte.precursor.cell": "OPC",
        "oligodendrocyte": "Oligodendrocyte",
        "PAX6.GABAergic.cortical.interneuron": "PAX6",
        "pvalb.GABAergic.cortical.interneuron": "PVALB",
        "sncg.GABAergic.cortical.interneuron": "SNCG",
        "sst.GABAergic.cortical.interneuron": "SST",
        "vip.GABAergic.cortical.interneuron": "VIP",
        "VIP.GABAergic.interneuron": "VIP",
        "deep.layer.non.IT": "deep layer non-IT",
        "lamp5.GABAergic.cortical.interneuron": "LAMP5",
        "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron": "L2/3-6 IT",
    },
    "subclass": {
        "astrocyte": "Astrocyte",
        "chandelier.pvalb.GABAergic.cortical.interneuron": "Chandelier",
        "endothelial.cell": "Endothelial",
        "macrophage": "Immune",
        "microglial.cell": "Microglia",
        "oligodendrocyte.precursor.cell": "OPC",
        "oligodendrocyte": "Oligodendrocyte",
        "PAX6.GABAergic.cortical.interneuron": "PAX6",
        "pvalb.GABAergic.cortical.interneuron": "PVALB",
        "pericyte": "Pericyte",
        "smooth.muscle.cell": "SMC",
        "sncg.GABAergic.cortical.interneuron": "SNCG",
        "sst.GABAergic.cortical.interneuron": "SST",
        "vip.GABAergic.cortical.interneuron": "VIP",
        "VIP.GABAergic.interneuron": "VIP",
        "vascular.leptomeningeal.cell": "VLMC",
        "lamp5.GABAergic.cortical.interneuron": "LAMP5",
        "L2.3.6.intratelencephalic.projecting.glutamatergic.neuron": "L2/3-6 IT",
        "L5.extratelencephalic.projecting.glutamatergic.cortical.neuron": "L5 ET",
        "near.projecting.glutamatergic.cortical.neuron": "L5/6 NP",
        "corticothalamic.projecting.glutamatergic.cortical.neuron": "L6 CT",
        "L6b.glutamatergic.cortical.neuron": "L6b",
    },
}

# Author (PsychENCODE) cell type -> display name; others fall back to author_display()
AUTHOR_DISPLAY = {
    "Astro": "Astrocyte", "Oligo": "Oligodendrocyte", "Micro": "Microglia",
    "Endo": "Endothelial", "PC": "Pericyte", "Vip": "VIP", "Pvalb": "PVALB",
    "Sst": "SST", "Sst_Chodl": "SST Chodl", "Lamp5": "LAMP5", "Lamp5_Lhx6": "LAMP5 LHX6",
    "Sncg": "SNCG", "Pax6": "PAX6",
}

# F1 colour range per level (class values all lie above 0.7; subclass go down to 0.46)
F1_RANGE = {"class": (0.7, 1.0), "subclass": (0.4, 1.0)}


def author_display(ct):
    """L2.3_IT -> L2/3 IT, L5.6_NP -> L5/6 NP, Sst_Chodl -> SST Chodl."""
    if ct in AUTHOR_DISPLAY:
        return AUTHOR_DISPLAY[ct]
    return re.sub(r"(\d)\.(\d)", r"\1/\2", ct).replace("_", " ")



def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _repel_labels(ax, xs, ys, texts, fontsize=9):
    """Label points without overlaps (no adjustText dependency): for each point, try
    offsets at growing radii in 16 directions and keep the first whose box stays inside
    the axes and clears earlier labels and all points. Leader lines for far labels."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pt2px = fig.dpi / 72.0
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    axbb = ax.get_window_extent(r)
    placed = []
    # densest points first so they get the closest slots
    d = np.array([np.sum(np.hypot(*(pts - q).T) < 40) for q in pts])
    angles = np.deg2rad([45, 0, 90, 315, 135, 180, 270, 225, 22, 67, 112, 157, 202, 247, 292, 337])
    for k in np.argsort(-d):
        (px, py), t = pts[k], texts[k]
        if not t:
            continue  # unlabelled point: still avoided by other labels
        probe = ax.text(0, 0, t, fontsize=fontsize)
        bb = probe.get_window_extent(r)
        probe.remove()
        w, h = bb.width, bb.height
        best = None
        for rad in (5, 10, 16, 24, 34, 46, 60, 76):
            for a in angles:
                ox, oy = rad * np.cos(a) * pt2px, rad * np.sin(a) * pt2px
                x0 = px + ox - (w if np.cos(a) < -0.2 else 0 if np.cos(a) > 0.2 else w / 2)
                y0 = py + oy - (h if np.sin(a) < -0.2 else 0 if np.sin(a) > 0.2 else h / 2)
                box = (x0 - 1, y0 - 1, x0 + w + 1, y0 + h + 1)
                if box[0] < axbb.x0 or box[2] > axbb.x1 or box[1] < axbb.y0 or box[3] > axbb.y1:
                    continue
                if any(not (box[2] < b[0] or box[0] > b[2] or box[3] < b[1] or box[1] > b[3]) for b in placed):
                    continue
                if any(box[0] - 4 < qx < box[2] + 4 and box[1] - 4 < qy < box[3] + 4 for qx, qy in pts):
                    continue
                # the label must sit closer to its own point than to any other point
                dist = [np.hypot(max(box[0] - qx, 0, qx - box[2]), max(box[1] - qy, 0, qy - box[3]))
                        for qx, qy in pts]
                if rad <= 10 and min(dd for m, dd in enumerate(dist) if m != k) <= dist[k] + 2:
                    continue
                best = (box, rad)
                break
            if best:
                break
        if best is None:
            best = ((px + 4, py + 4, px + 4 + w, py + 4 + h), 5)
        box, rad = best
        tx, ty = ax.transData.inverted().transform((box[0] + 1, box[1] + 1))
        ax.annotate(t, (xs[k], ys[k]), xytext=(tx, ty), textcoords="data", fontsize=fontsize,
                    ha="left", va="bottom",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="0.55", shrinkA=0, shrinkB=3)
                    if rad > 10 else None)
        placed.append(box)

CONTRAST_DISPLAY = {
    "Disorder_Schizophrenia_vs_Control": "SCZ",
    "Disorder_Bipolar.Disorder_vs_Control": "BD",
    "Disorder_ASD_vs_Control": "ASD",
    "Disorder_MDD_vs_Control": "MDD",
    "Disorder_PTSD_vs_Control": "PTSD",
    "Disorder_Williams.Syndrome_vs_Control": "WS",
    "Age_death": "Age",
    "Biological_Sex_male_vs_female": "Sex",
}


def contrast_pivot(combined_df, gemma_to_gemma_map, metric):
    """contrast (rows, display names) x author cell type table of a metric."""
    disp = dict(CONTRAST_DISPLAY)
    disp.update({k + "_pairwise_jaccard.tsv": v for k, v in CONTRAST_DISPLAY.items()})
    df = combined_df[combined_df["ct_author"].isin(gemma_to_gemma_map.keys())].copy()
    df["contrast_clean"] = df["contrast"].map(disp)
    df = df.dropna(subset=["contrast_clean"])
    pivot = df.pivot_table(index="contrast_clean", columns="ct_author", values=metric, aggfunc="mean")
    row_order = [c for c in ["SCZ", "BD", "ASD", "MDD", "PTSD", "WS", "Age", "Sex"] if c in pivot.index]
    return pivot.loc[row_order]


def load_benchmark_f1(f1_path, annotation_level, author_to_label):
    """author cell type -> benchmark F1 of its pipeline label at annotation_level."""
    to_bench = PIPELINE_TO_BENCH.get(annotation_level, PIPELINE_TO_BENCH["class"])
    f1_df = pd.read_csv(f1_path, sep="\t")
    if "key" in f1_df.columns:
        f1_df = f1_df[f1_df["key"] == annotation_level]
    if "metric" in f1_df.columns:
        f1_df = f1_df[f1_df["metric"] == "f1_score"]
    f1_col = "mean_f1_score" if "mean_f1_score" in f1_df.columns else "mean"
    bench_f1 = dict(zip(f1_df["label"], f1_df[f1_col]))
    return {ct: bench_f1.get(to_bench.get(cls, ""), np.nan) for ct, cls in author_to_label.items()}


def own_label_types(author_to_label, annotation_level):
    """Author types whose own name is a benchmark label (others share a merged label)."""
    to_bench = PIPELINE_TO_BENCH.get(annotation_level, PIPELINE_TO_BENCH["class"])
    return {ct for ct, cls in author_to_label.items() if _norm(author_display(ct)) == _norm(to_bench.get(cls, ""))}


def plot_correlation_heatmap(combined_df, gemma_to_gemma_map, metric="spearman_log2FoldChange", f1_path=None, annotation_level="class", output_path=None,
                             extra_bars=None, scatter_panels=None, f1_range=None):
    """Heatmap of a metric (contrast x author cell type) with annotation bars for the
    benchmark F1 of each author type's pipeline label and for the pipeline label itself,
    and with colour bars and the pipeline-label legend in the same figure.

    metric : column of combined_df, e.g. "spearman_log2FoldChange", "spearman_stat",
        "spearman_pvalue" or "jaccard_index".
    f1_path : TSV with key, label and mean_f1_score (or mean) columns, e.g. from
        extract_benchmark_f1.py. F1 is looked up for the pipeline label at annotation_level.
    extra_bars : list of (label, {author ct: F1}) drawn under the benchmark F1 bar, same scale.
    scatter_panels : list of dicts (x: {author ct: value}, xlabel, label: "all" or
        "outliers", note, xlim) plotted against the column mean of the metric, side by side
        under the heatmap. Default: one panel of benchmark F1 for types with their own label.
    """
    metric_settings = {
        "spearman_log2FoldChange": {"cbar_label": "Spearman ρ (log2 fold change)", "vmin": -1, "vmax": 1},
        "spearman_stat": {"cbar_label": "Spearman ρ (test statistic)", "vmin": -1, "vmax": 1},
        "spearman_pvalue": {"cbar_label": "Spearman ρ (p-value)", "vmin": -1, "vmax": 1},
        "jaccard_index": {"cbar_label": "Jaccard index (%)", "vmin": 0, "vmax": 100},
    }
    settings = metric_settings.get(metric, {"cbar_label": metric.replace("_", " "), "vmin": 0, "vmax": 1})

    pivot = contrast_pivot(combined_df, gemma_to_gemma_map, metric)
    row_order = list(pivot.index)

    author_to_class = {ct: gemma_to_gemma_map[ct] for ct in pivot.columns}
    to_bench = PIPELINE_TO_BENCH.get(annotation_level, PIPELINE_TO_BENCH["class"])
    class_display = lambda cls: to_bench.get(cls, cls.replace(".", " "))

    # --- benchmark F1 of each author type's pipeline label ---
    ct_to_f1 = load_benchmark_f1(f1_path, annotation_level, author_to_class) if f1_path is not None else {}
    extra_bars = extra_bars or []

    # --- column order: pipeline labels by F1 (or mean metric), then author types by mean metric ---
    col_mean = pivot.mean(axis=0)
    class_key = {}
    for cls in set(author_to_class.values()):
        members = [ct for ct, c in author_to_class.items() if c == cls]
        f1 = ct_to_f1.get(members[0], np.nan) if ct_to_f1 else np.nan
        class_key[cls] = f1 if not np.isnan(f1) else col_mean[members].mean() - 10
    cols = sorted(pivot.columns, key=lambda ct: (-class_key[author_to_class[ct]], author_to_class[ct], -col_mean[ct]))
    pivot = pivot[cols]

    classes_present = list(OrderedDict.fromkeys(author_to_class[ct] for ct in cols))
    class_palette = dict(zip(classes_present, sns.color_palette("tab20", n_colors=len(classes_present))))

    # --- layout (inches) ---
    ncol, nrow = len(cols), len(row_order)
    cell = 0.42
    left, right_pad, bottom, top = 1.9, 0.25, 1.55, 0.15
    bar_h, bar_gap = 0.30, 0.06
    n_bars = (2 if ct_to_f1 else 1) + len(extra_bars)
    cbar_block = 2.3 if ct_to_f1 else 1.25
    legend_w = 2.1
    hm_w, hm_h = ncol * cell, nrow * cell
    bars_h = n_bars * bar_h + n_bars * bar_gap
    if scatter_panels is None and ct_to_f1 and metric.startswith("spearman"):
        own = own_label_types(author_to_class, annotation_level)
        scatter_panels = [dict(x={ct: f for ct, f in ct_to_f1.items() if ct in own and not np.isnan(f)},
                               xlabel=f"Annotation F1 (benchmark, {annotation_level})", label="all",
                               note="types without their own benchmark label omitted", xlim=(0.4, 1.04))]
    scatter_panels = [sp for sp in (scatter_panels or []) if sum(not np.isnan(v) for v in sp["x"].values()) >= 3]
    show_b = bool(scatter_panels)
    b_block = 4.9 if show_b else 0.0
    bottom += b_block
    fig_w = left + hm_w + 0.35 + cbar_block + legend_w + right_pad
    fig_h = bottom + hm_h + bars_h + top
    fig = plt.figure(figsize=(fig_w, fig_h))
    fx = lambda x: x / fig_w
    fy = lambda y: y / fig_h

    cmap = plt.cm.RdYlBu_r.copy()
    cmap.set_bad((0.85, 0.85, 0.85))
    ax = fig.add_axes([fx(left), fy(bottom), fx(hm_w), fy(hm_h)])
    mesh = ax.pcolormesh(np.ma.masked_invalid(pivot.values), cmap=cmap, vmin=settings["vmin"],
                         vmax=settings["vmax"], edgecolors="white", linewidth=0.6)
    ax.set_xlim(0, ncol)
    ax.set_ylim(nrow, 0)
    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_xticklabels([author_display(ct) for ct in cols], rotation=90, fontsize=11)
    ax.set_yticks(np.arange(nrow) + 0.5)
    ax.set_yticklabels(row_order, fontsize=12)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("Author cell type", fontsize=13)
    ax.set_ylabel("Contrast", fontsize=13)

    def add_bar(y0, colors, label):
        bax = fig.add_axes([fx(left), fy(y0), fx(hm_w), fy(bar_h)])
        bax.imshow(np.array([colors]), aspect="auto", extent=(0, ncol, 0, 1), interpolation="nearest")
        for i in range(1, ncol):
            bax.axvline(i, color="white", lw=0.6)
        bax.set_xticks([])
        bax.set_yticks([0.5])
        bax.set_yticklabels([label], fontsize=11)
        bax.tick_params(length=0)
        for s in bax.spines.values():
            s.set_visible(False)

    y = bottom + hm_h + bar_gap
    add_bar(y, [mcolors.to_rgba(class_palette[author_to_class[ct]]) for ct in cols], "Pipeline label")
    f1_cmap = plt.cm.YlOrRd
    f1_lo, f1_hi = f1_range or F1_RANGE.get(annotation_level, (0.4, 1.0))
    f1_norm = plt.Normalize(vmin=f1_lo, vmax=f1_hi, clip=True)
    # bottom-up: pipeline label, extra F1 bars (reversed, so the first sits just under the
    # benchmark bar), benchmark F1 on top
    f1_color = lambda v: f1_cmap(f1_norm(v)) if not np.isnan(v) else (0.85, 0.85, 0.85, 1.0)
    for lab, vals in reversed(extra_bars):
        y += bar_h + bar_gap
        add_bar(y, [f1_color(vals.get(ct, np.nan)) for ct in cols], lab)
    if ct_to_f1:
        y += bar_h + bar_gap
        add_bar(y, [f1_color(ct_to_f1[ct]) for ct in cols], "Benchmark F1" if extra_bars else "Annotation F1")

    # --- colour bars ---
    cb_h = hm_h + bars_h
    cx = left + hm_w + 0.45
    cax = fig.add_axes([fx(cx), fy(bottom), fx(0.16), fy(cb_h)])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label(settings["cbar_label"], fontsize=12)
    cb.ax.tick_params(labelsize=10)
    if ct_to_f1:
        cx += 1.1
        fax = fig.add_axes([fx(cx), fy(bottom), fx(0.16), fy(cb_h)])
        sm = plt.cm.ScalarMappable(cmap=f1_cmap, norm=f1_norm)
        sm.set_array([])
        fb = fig.colorbar(sm, cax=fax)
        fb.set_label("Annotation F1" if extra_bars else f"Annotation F1 (benchmark, {annotation_level})", fontsize=12)
        fb.ax.tick_params(labelsize=10)

    # --- pipeline label legend ---
    lax = fig.add_axes([fx(left + hm_w + 0.35 + cbar_block), fy(bottom), fx(legend_w), fy(cb_h)])
    lax.axis("off")
    handles = [mpatches.Patch(facecolor=class_palette[c], label=class_display(c)) for c in classes_present]
    lax.legend(handles=handles, title="Pipeline label", loc="upper left", fontsize=10, title_fontsize=11,
               frameon=False, borderaxespad=0, handlelength=1.2)

    if show_b:
        from scipy.stats import spearmanr
        n_p = len(scatter_panels)
        gap = 0.9
        pw = min((hm_w - gap * (n_p - 1)) / n_p, 10.5)
        rows_out = []
        letters = "BCDEFG"
        for pi, sp in enumerate(scatter_panels):
            b_cts = [ct for ct in cols if not np.isnan(sp["x"].get(ct, np.nan))]
            bx = np.array([sp["x"][ct] for ct in b_cts])
            by = np.array([col_mean[ct] for ct in b_cts])
            rho_b, p_b = spearmanr(bx, by)
            x0 = left + pi * (pw + gap)
            bax = fig.add_axes([fx(x0), fy(0.75), fx(pw), fy(b_block - 1.25)])
            bax.scatter(bx, by, s=46, c=[class_palette[author_to_class[ct]] for ct in b_cts],
                        edgecolors="0.25", linewidths=0.5, zorder=3)
            if sp.get("label", "all") == "outliers":
                lab_idx = [k for k in range(len(b_cts)) if by[k] < 0.85 or bx[k] <= np.quantile(bx, 0.2)]
            else:
                lab_idx = list(range(len(b_cts)))
            _repel_labels(bax, bx, by, [author_display(ct) if k in lab_idx else "" for k, ct in enumerate(b_cts)], fontsize=9)
            bax.set_xlabel(sp["xlabel"], fontsize=12)
            if pi == 0:
                bax.set_ylabel(f"Mean {settings['cbar_label']}", fontsize=12)
            bax.tick_params(labelsize=10)
            if sp.get("xlim"):
                bax.set_xlim(*sp["xlim"])
            else:
                pad = (bx.max() - bx.min()) * 0.06
                bax.set_xlim(bx.min() - pad, bx.max() + pad)
            bax.set_ylim(np.floor(by.min() * 10) / 10 - 0.05, 1.06)
            bax.set_yticks(np.arange(np.ceil(by.min() * 10) / 10 - 0.1, 1.01, 0.1))
            bax.grid(alpha=0.25, lw=0.5)
            for side in ("top", "right"):
                bax.spines[side].set_visible(False)
            p_txt = f"p = {p_b:.2g}" if p_b >= 0.001 else "p < 0.001"
            note = f"\n({sp['note']})" if sp.get("note") else ""
            bax.text(0.02, 0.97, f"Spearman ρ = {rho_b:.2f}, {p_txt}, n = {len(b_cts)}{note}",
                     transform=bax.transAxes, ha="left", va="top", fontsize=10)
            fig.text(fx(0.15 if pi == 0 else x0 - 0.75), fy(b_block - 0.05), letters[pi], fontsize=18, fontweight="bold", va="top")
            rows_out += [{"panel": letters[pi], "xlabel": sp["xlabel"], "ct_author": ct, "display_name": author_display(ct),
                          "x": xv, f"mean_{metric}": yv} for ct, xv, yv in zip(b_cts, bx, by)]
        pd.DataFrame(rows_out).to_csv(
            (output_path or f"heatmap_{metric}.png").replace(".png", "_scatter.tsv"), sep="\t", index=False)
        fig.text(fx(0.15), fy(fig_h - 0.12), "A", fontsize=18, fontweight="bold", va="top")

    if output_path is None:
        output_path = f"heatmap_{metric}.png"
    fig.savefig(output_path, dpi=300, facecolor="white")
    plt.close(fig)

    pivot.to_csv(output_path.replace(".png", ".tsv"), sep="\t")
    pd.DataFrame({
        "ct_author": cols,
        "display_name": [author_display(ct) for ct in cols],
        "pipeline_label": [author_to_class[ct] for ct in cols],
        "benchmark_label": [to_bench.get(author_to_class[ct], "") for ct in cols],
        "annotation_f1": [ct_to_f1.get(ct, np.nan) for ct in cols],
        **{lab: [vals.get(ct, np.nan) for ct in cols] for lab, vals in extra_bars},
    }).to_csv(output_path.replace(".png", "_columns.tsv"), sep="\t", index=False)
