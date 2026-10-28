"""Manuscript figure: DE agreement vs annotation F1 and cell support (sc-manuscript Figure 4).

Two panels at final size (183 mm, double column). A: heatmap of DE agreement (Spearman rho
of log2 fold changes, author vs pipeline labels) per contrast and author type, with bars for
benchmark F1, CMC F1 and cell class, and brackets over author types pooled into one pipeline
label. B: mean rho across contrasts against log10 median CMC cells per donor, coloured by CMC
F1 on the same scale as the F1 bars in A, with a LOWESS smooth; open markers are pooled types.
Text sizes follow the lab rule: >= 10 pt text, >= 12 pt axis titles, >= 14 pt panel letters.
The deliverable is the PNG. With qa_dir, a layout-alignment JSON and a temporary PDF for the
text-size and collision audits are written there (QA only, not a figure export).
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

from utils import _repel_labels, author_display

TXT, AX, LET = 10, 12, 14
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "font.size": TXT,
    "axes.labelsize": AX,
    "xtick.labelsize": TXT,
    "ytick.labelsize": TXT,
    "legend.fontsize": TXT,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

MM = 1 / 25.4
CLASS_OF = {
    "OPC": "Glia", "Astro": "Glia", "Oligo": "Glia",
    "Micro": "Immune", "Immune": "Immune",
    "Endo": "Vascular", "PC": "Vascular", "SMC": "Vascular", "VLMC": "Vascular",
    "Chandelier": "Inhibitory", "Vip": "Inhibitory", "Pvalb": "Inhibitory", "Sst": "Inhibitory",
    "Sst_Chodl": "Inhibitory", "Lamp5": "Inhibitory", "Lamp5_Lhx6": "Inhibitory",
    "Sncg": "Inhibitory", "Pax6": "Inhibitory",
}
CLASS_COLORS = {  # Okabe-Ito subset (colour-blind safe), one hue per cell class
    "Excitatory": "#0072B2", "Inhibitory": "#E69F00", "Glia": "#009E73",
    "Immune": "#CC79A7", "Vascular": "#7F7F7F",
}


def cell_class(ct):
    return CLASS_OF.get(ct, "Excitatory" if ct.startswith("L") else "Glia")


def plot_figure(pivot, tab, label_map, bench_label, out_png, stats_lines, study="CMC", qa_dir=None):
    """pivot: contrast x author type (log2FC rho); tab: celltype_summary rows;
    stats_lines: correlation text for panel B (computed by the caller)."""
    t = tab.set_index("ct_author")
    lab_f1 = t.groupby("pipeline_label")["benchmark_f1"].first()
    cols = sorted(pivot.columns, key=lambda c: (-np.nan_to_num(lab_f1[label_map[c]], nan=-1),
                                                label_map[c], -t.loc[c, "mean_rho"]))
    pv = pivot[cols]
    ncol, nrow = len(cols), len(pv.index)
    pooled = {c for c in cols if sum(label_map[o] == label_map[c] for o in cols) > 1}

    W = 183 * MM
    left, right = 1.12, 0.08
    top = 0.32  # room for the A letter above the bars
    hm_w = W - left - right
    cell_h, bar_h, gap = 0.21, 0.17, 0.04
    hm_h = nrow * cell_h
    y_b, b_h = 0.62, 3.85                     # panel B
    y_hm = y_b + b_h + 0.35 + 1.55            # gap + A tick labels and axis title
    y_br = y_hm + hm_h + gap
    br_h = 0.30
    bars_y0 = y_br + br_h + gap
    H = bars_y0 + 3 * (bar_h + gap) + top
    fig = plt.figure(figsize=(W, H))
    fx, fy = (lambda v: v / W), (lambda v: v / H)

    # ---------- A: heatmap ----------
    cmap = plt.cm.RdYlBu_r.copy()
    cmap.set_bad("#D8D8D8")
    ax = fig.add_axes([fx(left), fy(y_hm), fx(hm_w), fy(hm_h)])
    mesh = ax.pcolormesh(np.ma.masked_invalid(pv.values), cmap=cmap, vmin=-1, vmax=1,
                         edgecolors="white", linewidth=0.5)
    ax.set_xlim(0, ncol); ax.set_ylim(nrow, 0)
    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_xticklabels([author_display(c) for c in cols], rotation=90, ha="right", va="center",
                       rotation_mode="anchor")
    ax.set_yticks(np.arange(nrow) + 0.5); ax.set_yticklabels(pv.index)
    ax.tick_params(length=0, pad=3)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_ylabel("Contrast")
    ax.set_xlabel("Author cell type", labelpad=4)

    f1_lo = min(0.3, np.floor(np.nanmin(t[["benchmark_f1", "cmc_f1"]].values) * 10) / 10)
    f1_cmap, f1_norm = plt.cm.YlOrRd, mpl.colors.Normalize(f1_lo, 1.0)
    f1c = lambda v: f1_cmap(f1_norm(v)) if not np.isnan(v) else (0.85, 0.85, 0.85, 1)

    def bar(y0, colors, label):
        b = fig.add_axes([fx(left), fy(y0), fx(hm_w), fy(bar_h)])
        b.imshow(np.array([colors]), aspect="auto", extent=(0, ncol, 0, 1), interpolation="nearest")
        for i in range(1, ncol):
            b.axvline(i, color="white", lw=0.5)
        b.set_xticks([]); b.set_yticks([0.5]); b.set_yticklabels([label])
        b.tick_params(length=0, pad=3)
        for s in b.spines.values():
            s.set_visible(False)

    bar(bars_y0, [mpl.colors.to_rgba(CLASS_COLORS[cell_class(c)]) for c in cols], "Cell class")
    bar(bars_y0 + (bar_h + gap), [f1c(t.loc[c, "cmc_f1"]) for c in cols], f"{study} F1")
    bar(bars_y0 + 2 * (bar_h + gap), [f1c(t.loc[c, "benchmark_f1"]) for c in cols], "Benchmark F1")

    # brackets over author types pooled into one pipeline label
    bax = fig.add_axes([fx(left), fy(y_br), fx(hm_w), fy(br_h)])
    bax.set_xlim(0, ncol); bax.set_ylim(0, 1); bax.axis("off")
    i = 0
    while i < ncol:
        j = i
        while j + 1 < ncol and label_map[cols[j + 1]] == label_map[cols[i]]:
            j += 1
        if j > i:
            bax.plot([i + 0.15, i + 0.15, j + 0.85, j + 0.85], [0.05, 0.25, 0.25, 0.05], color="0.25", lw=0.8)
            bax.text((i + j + 1) / 2, 0.32, bench_label.get(label_map[cols[i]], ""), ha="center", va="bottom")
        i = j + 1

    # ---------- B: mean rho vs support, coloured by CMC F1 ----------
    key_w = 2.2
    b_w = W - left - right - key_w - 0.4
    bx_ = fig.add_axes([fx(left), fy(y_b), fx(b_w), fy(b_h)])
    x = t.loc[cols, "log10_author_median_cells"].values.astype(float)
    y = t.loc[cols, "mean_rho"].values.astype(float)
    f = t.loc[cols, "cmc_f1"].values.astype(float)
    fit = lowess(y, x, frac=0.8, it=1, return_sorted=True)
    bx_.plot(fit[:, 0], fit[:, 1], color="0.35", lw=1.2, zorder=2)
    for xi, yi, fi, c in zip(x, y, f, cols):
        if c in pooled:
            bx_.scatter(xi, yi, s=34, facecolors="white", edgecolors=[f1c(fi)], linewidths=1.6, zorder=3)
        else:
            bx_.scatter(xi, yi, s=34, facecolors=[f1c(fi)], edgecolors="0.15", linewidths=0.6, zorder=3)
    bx_.set_ylim(np.floor(np.nanmin(y) * 10) / 10 - 0.05, 1.05)
    pad = (np.nanmax(x) - np.nanmin(x)) * 0.05
    bx_.set_xlim(np.nanmin(x) - pad, np.nanmax(x) + pad)
    xs = np.linspace(fit[0, 0], fit[-1, 0], 300)
    _repel_labels(bx_, x, y, [author_display(c) if yi < 0.8 else "" for yi, c in zip(y, cols)],
                  fontsize=TXT, avoid_xy=np.column_stack([xs, np.interp(xs, fit[:, 0], fit[:, 1])]))
    bx_.set_xlabel(f"log₁₀(median {study} cells per donor + 1)")
    bx_.set_ylabel("Mean Spearman ρ")

    # key column beside B: class legend, marker legend, two colour bars
    kx = left + b_w + 0.35
    handles = [mpatches.Patch(facecolor=CLASS_COLORS[k], label=k) for k in CLASS_COLORS]
    handles += [mlines.Line2D([], [], ls="", marker="o", mfc="0.55", mec="0.15", ms=6, label="Own label"),
                mlines.Line2D([], [], ls="", marker="o", mfc="white", mec="0.35", mew=1.6, ms=6, label="Pooled label")]
    sax = fig.add_axes([fx(kx), fy(y_b + b_h - 0.75), fx(key_w), fy(0.75)])
    sax.axis("off")
    sax.text(0, 1, "\n".join(stats_lines), ha="left", va="top", linespacing=1.3)
    lax = fig.add_axes([fx(kx), fy(y_b + 1.55), fx(key_w), fy(b_h - 0.85 - 1.55)])
    lax.axis("off")
    lax.legend(handles=handles, loc="upper left", borderaxespad=0, handlelength=1.0, handleheight=0.9,
               labelspacing=0.3, handletextpad=0.5, title="Cell class (A)", title_fontsize=TXT,
               alignment="left")
    for k, (sm, lab) in enumerate([(mesh, "Spearman ρ (A)"),
                                   (mpl.cm.ScalarMappable(f1_norm, f1_cmap), "Annotation F1 (A, B)")]):
        cax = fig.add_axes([fx(kx), fy(y_b + 1.05 - k * 0.72), fx(key_w - 0.15), fy(0.10)])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=TXT, labelpad=5)
        cax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=TXT, length=2.5, width=0.6)
        cb.set_ticks([-1, 0, 1] if k == 0 else [f1_lo, 1.0])
        cb.outline.set_linewidth(0.6)

    fig.text(fx(0.03), fy(H - 0.03), "A", fontsize=LET, fontweight="bold", va="top")
    fig.text(fx(0.03), fy(y_b + b_h + 0.30), "B", fontsize=LET, fontweight="bold", va="top")

    if qa_dir:
        import sys
        sys.path.insert(0, os.path.expanduser("~/.claude/skills/nature-figure/scripts"))
        from audit_panel_alignment import require_matplotlib_panel_alignment
        os.makedirs(qa_dir, exist_ok=True)
        # B is narrower than A by design (shared key column); left edges must match
        # A is a full-width hero panel and B is narrower (shared key column); the gate has no
        # comparable group for that layout and reports NOT AUDITABLE, which is recorded, not
        # treated as a pass. Any real misalignment it finds (FIX BEFORE DELIVERY) still raises.
        from audit_panel_alignment import PanelAlignmentError
        try:
            require_matplotlib_panel_alignment(
                fig, json_out=os.path.join(qa_dir, "alignment.json"), tolerance_pt=1.5,
                gutter_tolerance_pt=1.5, strict=True, axes=[ax, bx_], panel_ids=["A", "B"])
        except PanelAlignmentError as e:
            if "NOT AUDITABLE" not in str(e):
                raise
        fig.canvas.draw()
        assert abs(ax.get_position().x0 - bx_.get_position().x0) < 1.5 / 72 / W, "A and B left edges differ"
        fig.savefig(os.path.join(qa_dir, "figure_qa.pdf"))
    fig.savefig(out_png, dpi=600, facecolor="white")
    plt.close(fig)
