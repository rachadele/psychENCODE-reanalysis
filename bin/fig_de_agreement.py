"""Manuscript figure: DE agreement vs annotation F1 and cell support (sc-manuscript Figure 4).

At final size (183 mm, double column). A: heatmap of DE agreement (Spearman rho of log2 fold
changes, author vs pipeline labels) per contrast and author type, with bars for benchmark F1,
CMC F1 and cell class, and brackets over author types pooled into one pipeline label.
B, C: mean rho across contrasts against CMC F1 and against log10 median CMC cells per donor,
points coloured by cell class (open: pooled label), OLS line with 95% confidence band, and
Spearman rho, p and n above each panel. Text: >= 10 pt, axis titles >= 12 pt, letters 14 pt.
The only output is the PNG; check_layout() verifies text sizes, page bounds, text overlaps
and panel alignment on the rendered figure and blocks export on failure (report to qa_dir).
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np

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


def check_layout(fig, ax_a, panels, W):
    """Rendered checks, no vector export: text sizes (>= 10 pt, axis titles >= 12 pt), text
    inside the page, no overlapping text boxes, A spanning B to C, B and C the same size."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    fails, texts = [], []
    # axis titles of data axes (colour-bar labels count as ordinary text, >= 10 pt)
    data_axes = [ax_a, *panels]
    titles = {id(a.xaxis.label) for a in data_axes} | {id(a.yaxis.label) for a in data_axes}
    skip = set()
    for a in fig.axes:
        for axis, (lo, hi) in ((a.xaxis, a.get_xlim()), (a.yaxis, a.get_ylim())):
            lo, hi = min(lo, hi), max(lo, hi)
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                if not a.axison or not (lo - 1e-9 <= tk.get_loc() <= hi + 1e-9):
                    skip.update({id(tk.label1), id(tk.label2)})
    for t in fig.findobj(mpl.text.Text):
        if id(t) in skip or not t.get_visible() or not t.get_text().strip():
            continue
        size = t.get_fontsize()
        floor = AX if id(t) in titles else TXT
        if size < floor - 1e-6:
            fails.append(f"text {t.get_text()!r} is {size} pt (< {floor})")
        bb = t.get_window_extent(r)
        if bb.x0 < fb.x0 - 1 or bb.x1 > fb.x1 + 1 or bb.y0 < fb.y0 - 1 or bb.y1 > fb.y1 + 1:
            fails.append(f"text {t.get_text()!r} extends beyond the page")
        texts.append((t.get_text(), bb))
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a_, b_ = texts[i][1], texts[j][1]
            ov = max(0, min(a_.x1, b_.x1) - max(a_.x0, b_.x0)) * max(0, min(a_.y1, b_.y1) - max(a_.y0, b_.y0))
            if ov > 0.05 * min(a_.width * a_.height, b_.width * b_.height):
                fails.append(f"text {texts[i][0]!r} overlaps {texts[j][0]!r}")
    tol = 1.5 / 72 / W
    pa, pb, pc = (q.get_position() for q in (ax_a, *panels))
    if abs(pa.x0 - pb.x0) > tol or abs(pa.x1 - pc.x1) > tol:
        fails.append("panel A does not span B to C")
    if abs(pb.width - pc.width) > tol or abs(pb.height - pc.height) > tol or abs(pb.y0 - pc.y0) > tol:
        fails.append("panels B and C differ in size or baseline")
    return {"n_text": len(texts), "failures": fails}


def plot_figure(pivot, tab, label_map, bench_label, out_png, study="CMC", qa_dir=None):
    """pivot: contrast x author type (log2FC rho); tab: celltype_summary rows."""
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
    y_b, b_h = 0.62, 2.9                      # panels B and C
    y_key, key_h = y_b + b_h + 0.45, 0.85     # key row: class legend and colour bars
    y_hm = y_key + key_h + 0.05 + 1.45        # A tick labels and axis title
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

    # ---------- key row: class legend, marker legend, colour bars ----------
    handles = [mpatches.Patch(facecolor=CLASS_COLORS[k], label=k) for k in CLASS_COLORS]
    handles += [mlines.Line2D([], [], ls="", marker="o", mfc="0.55", mec="0.15", ms=6, label="Own label"),
                mlines.Line2D([], [], ls="", marker="o", mfc="white", mec="0.35", mew=1.6, ms=6, label="Pooled label")]
    lax = fig.add_axes([fx(left - 0.9), fy(y_key), fx(4.6), fy(key_h)])
    lax.axis("off")
    lax.legend(handles=handles, loc="center left", ncol=4, borderaxespad=0, handlelength=1.0,
               handleheight=0.9, columnspacing=0.9, labelspacing=0.4, handletextpad=0.4)
    cb_w = 0.95
    for k, (sm, lab) in enumerate([(mesh, "Spearman ρ"),
                                   (mpl.cm.ScalarMappable(f1_norm, f1_cmap), "Annotation F1")]):
        cx0 = W - right - 2 * cb_w - 0.55 + k * (cb_w + 0.3)
        cax = fig.add_axes([fx(cx0), fy(y_key + 0.3), fx(cb_w), fy(0.10)])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=TXT, labelpad=5)
        cax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=TXT, length=2.5, width=0.6)
        cb.set_ticks([-1, 0, 1] if k == 0 else [f1_lo, 1.0])
        cb.outline.set_linewidth(0.6)

    # ---------- B, C: mean rho vs CMC F1 and vs cells per donor, OLS line with 95% band ----------
    from scipy import stats as sst
    gutter = 0.45
    pw = (hm_w - gutter) / 2
    y = t.loc[cols, "mean_rho"].values.astype(float)
    panels = [("cmc_f1", f"Annotation F1 ({study})"),
              ("log10_author_median_cells", f"Median {study} cells per donor + 1")]
    sax = []
    for k, (col, xl) in enumerate(panels):
        a = fig.add_axes([fx(left + k * (pw + gutter)), fy(y_b), fx(pw), fy(b_h)])
        x = t.loc[cols, col].values.astype(float)
        # OLS line and 95% confidence band for the mean, on the plotted axes
        n = len(x)
        slope, icpt = np.polyfit(x, y, 1)
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 200)
        resid = y - (icpt + slope * x)
        s_err = np.sqrt(np.sum(resid ** 2) / (n - 2))
        se = s_err * np.sqrt(1 / n + (xs - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
        tq = sst.t.ppf(0.975, n - 2)
        yfit = icpt + slope * xs
        logx = col.startswith("log10")
        tr = (lambda v: 10 ** v) if logx else (lambda v: v)   # C is drawn on a log axis
        xs_p = tr(xs)
        for xi, yi, c in zip(tr(x), y, cols):
            colr = CLASS_COLORS[cell_class(c)]
            if c in pooled:
                a.scatter(xi, yi, s=34, facecolors="white", edgecolors=colr, linewidths=1.6, zorder=3)
            else:
                a.scatter(xi, yi, s=34, facecolors=colr, edgecolors="0.15", linewidths=0.6, zorder=3)
        a.fill_between(xs_p, yfit - tq * se, yfit + tq * se, color="0.85", lw=0, zorder=1)
        a.plot(xs_p, yfit, color="0.3", lw=1.2, zorder=2)
        a.set_ylim(np.floor(np.nanmin(y) * 10) / 10 - 0.05, 1.05)
        pad = (np.nanmax(x) - np.nanmin(x)) * 0.06
        if logx:
            a.set_xscale("log")
            a.set_xlim(10 ** (np.nanmin(x) - pad), 10 ** (np.nanmax(x) + pad))
            ticks = [v for v in (1, 10, 100, 1000, 10000) if v <= 10 ** np.nanmax(x)]
            a.set_xticks(ticks); a.set_xticklabels([f"{v:g}" for v in ticks])
            a.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        else:
            a.set_xlim(np.nanmin(x) - pad, np.nanmax(x) + pad)
        lab = [author_display(c) if (yi < 0.8 or xi <= np.quantile(x, 0.1)) else ""
               for xi, yi, c in zip(x, y, cols)]
        _repel_labels(a, tr(x), y, lab, fontsize=TXT, avoid_xy=np.column_stack([xs_p, yfit]))
        a.set_xlabel(xl)
        if k == 0:
            a.set_ylabel("Mean Spearman ρ")
        else:
            a.set_yticklabels([])
        rho, pval = sst.spearmanr(x, y)
        ptxt = f"p = {pval:.2g}" if pval >= 0.001 else "p < 0.001"
        a.set_title(f"Spearman ρ = {rho:.2f}, {ptxt}, n = {n}", loc="left", fontsize=TXT, pad=4)
        sax.append(a)

    fig.text(fx(0.03), fy(H - 0.03), "A", fontsize=LET, fontweight="bold", va="top")
    for k in range(2):
        fig.text(fx(0.03 if k == 0 else left + pw + gutter - 0.3), fy(y_b + b_h + 0.42), "BC"[k],
                 fontsize=LET, fontweight="bold", va="top")

    qa = check_layout(fig, ax, sax, W)
    if qa_dir:
        import json
        os.makedirs(qa_dir, exist_ok=True)
        json.dump(qa, open(os.path.join(qa_dir, "layout_qa.json"), "w"), indent=2)
    if qa["failures"]:
        raise RuntimeError("Figure QA failed: " + "; ".join(qa["failures"]))
    fig.savefig(out_png, dpi=600, facecolor="white")
    plt.close(fig)
