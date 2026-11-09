"""Manuscript figure: DE agreement vs annotation F1 and cell support (sc-manuscript Figure 4).

At final size (183 mm, double column). A: heatmap of DE agreement (Spearman rho of log2 fold
changes, author vs pipeline labels) per contrast and author type, with bars for benchmark F1,
CMC F1 and cell class, and brackets over author types pooled into one pipeline label.
B: mean rho against log10 median CMC cells per donor, Spearman rho and p above. C, D:
added-variable plots for the partial Spearman correlation of mean rho with CMC F1 given cells per
donor (rank residuals of each after regression on rank log10 median cells per donor), for all 27
types (C) and, with residuals computed within the subset, for the 20 types left after dropping
the seven pooled under another type's name (D); plain -> partial rho above each. Points coloured
by cell class (open: those seven pooled types), OLS line with 95% confidence band. E: density of
per-pair rho for matched author/pipeline cell types (above) and mismatched types in the same contrast
(below), observed means as ticks, means adjusted for cells per donor and contrast as diamonds with 95% CI,
permutation null of the matched mean in red (benchmark-post-hoc/de_permutation_null). Text: >= 10 pt, axis titles >= 12 pt, letters 14 pt.
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
# Cell class = class level of the human hierarchy (annotation-benchmark/meta/census_map_human.tsv),
# read from the author -> pipeline class map, in hierarchy order (excitatory, inhibitory CGE then
# MGE, glia, immune, vascular). Colours come from the shared sc-manuscript palette (one TSV for every
# figure; MANUSCRIPT_PALETTE overrides the path): shades of one hue group per family, checked for
# colour-vision deficiencies (benchmark-post-hoc/palette_cvd_check.py).
PALETTE_TSV = os.environ.get("MANUSCRIPT_PALETTE",
                             "/space/grp/rschwartz/rschwartz/benchmark-post-hoc/manuscript_palette.tsv")
CLASS_ORDER = ["L2/3-6 IT", "Deep layer non-IT", "LAMP5", "VIP", "SNCG", "PAX6", "PVALB", "SST",
               "Chandelier", "OPC", "Oligodendrocyte", "Astrocyte", "Microglia", "Immune", "Vascular"]


def _class_colors():
    import csv
    pal = {r["label"].lower(): r["hex"] for r in csv.DictReader(open(PALETTE_TSV), delimiter="\t")
           if r["level"] == "class"}
    return {k: pal[k.lower()] for k in CLASS_ORDER}


CLASS_COLORS = _class_colors()
EDGE = "0.35"  # outline so the lightest class shades separate from white
_CLASS_MAP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets",
                          "author_to_pipeline_map.class.json")


def _class_of():
    import json
    from utils import PIPELINE_TO_BENCH
    disp = {k: ("Deep layer non-IT" if v == "deep layer non-IT" else v)
            for k, v in PIPELINE_TO_BENCH["class"].items()}
    return {a: disp[lab] for a, lab in json.load(open(_CLASS_MAP)).items()}


CLASS_OF = _class_of()


def cell_class(ct):
    return CLASS_OF[ct]


# Matched vs mismatched cell type pairs (benchmark-post-hoc/de_permutation_null). One row per arm:
# density of per-pair rho for matched pairs above the baseline, for mismatched pairs (other cell
# types, same contrast) below it, a tick at each per-pair mean and a red bar at the permutation
# null of the matched mean (mean +/- SD); means adjusted for cells per donor as diamonds with 95% CI. The same drawing is in the other repo's figure module.
MATCH_C, MISMATCH_C, NULL_C = "0.45", "0.78", "#C0392B"


def load_null(null_dir, test, arms):
    """[(arm, matched rho, mismatched rho, summary row)] from pairs.tsv, null_summary.tsv and
    matched_effect_adjusted.tsv (adjusted means added to the row as adj_<group>, _lo, _hi),
    checking that the per-pair means equal both tables."""
    import pandas as pd
    pairs = pd.read_csv(os.path.join(null_dir, "pairs.tsv"), sep="\t")
    summ = pd.read_csv(os.path.join(null_dir, "null_summary.tsv"), sep="\t")
    # marginal means adjusted for cells per donor and contrast (matched_effect_adjusted.py)
    eff = pd.read_csv(os.path.join(null_dir, "matched_effect_adjusted.tsv"), sep="\t")
    out = []
    for arm in arms:
        p = pairs[(pairs["test"] == test) & (pairs["arm"] == arm)]
        s = summ[(summ["test"] == test) & (summ["arm"] == arm) & (summ["null"] == "mismatched_celltype")].iloc[0]
        m, mm = p.loc[p["matched"], "rho"].values, p.loc[~p["matched"], "rho"].values
        assert len(m) == s["n_pairs"] and len(mm) == s["n_mismatched_pairs"]
        assert abs(m.mean() - s["observed_mean"]) < 1e-5 and abs(mm.mean() - s["mismatched_pair_mean"]) < 1e-5
        e = eff[(eff["test"] == test) & (eff["arm"] == arm)]
        s = s.copy()
        for grp, v in (("matched", m), ("mismatched", mm)):
            ob = e[(e["model"] == "observed_mean") & (e["group"] == grp)].iloc[0]
            assert abs(ob["estimate"] - v.mean()) < 1e-5
            r = e[(e["model"] == "adjusted_mean") & (e["group"] == grp)].iloc[0]
            s[f"adj_{grp}"], s[f"adj_{grp}_lo"], s[f"adj_{grp}_hi"] = r["estimate"], r["ci_low"], r["ci_high"]
        out.append((arm, m, mm, s))
    return out


def null_density_rows(ax, rows, labels, half=0.42):
    """rows from load_null; labels: y tick label per row (top to bottom). Per row: density of
    matched rho above the baseline and mismatched rho below, observed means as black ticks (value at
    the tip), means adjusted for cells per donor as hollow diamonds with 95% CI bars (no value label),
    permutation null of the matched mean as a
    red bar. Sets the x limits to the densities; returns them."""
    from scipy.stats import gaussian_kde
    lo = min(min(r[1].min(), r[2].min()) for r in rows)
    hi = max(max(r[1].max(), r[2].max()) for r in rows)
    xs = np.linspace(lo - 0.05, min(hi + 0.05, 1.0), 400)  # tails past the data, never past rho = 1
    for i, (_, m, mm, s) in enumerate(rows):
        y0 = -i
        for grp, v, sign, c in (("matched", m, 1, MATCH_C), ("mismatched", mm, -1, MISMATCH_C)):
            d = gaussian_kde(v)(xs)
            d = d / d.max() * half
            ax.fill_between(xs, y0, y0 + sign * d, color=c, lw=0, zorder=2)
            mu = v.mean()
            ax.plot([mu, mu], [y0, y0 + sign * half * 1.05], color="black", lw=1.2, zorder=3)
            ax.text(mu, y0 + sign * half * 1.08, f"{mu:.2f}", ha="center",
                    va="bottom" if sign > 0 else "top", fontsize=TXT, zorder=4)
            a, a_lo, a_hi = s[f"adj_{grp}"], s[f"adj_{grp}_lo"], s[f"adj_{grp}_hi"]
            ya = y0 + sign * half * 0.62
            ax.plot([a_lo, a_hi], [ya, ya], color="black", lw=1.2, zorder=6, solid_capstyle="butt")
            ax.plot([a], [ya], marker="D", ms=6, mfc="white", mec="black", mew=1.2, ls="", zorder=7)
        nm, ns = s["null_mean"], s["null_sd"]
        ax.add_patch(mpatches.Rectangle((nm - max(ns, 0.006), y0 - half * 0.32), 2 * max(ns, 0.006), half * 0.64,
                                        color=NULL_C, lw=0, zorder=5))
        ax.plot([xs[0], xs[-1]], [y0, y0], color="0.5", lw=0.6, zorder=1)
    ax.set_yticks([-i for i in range(len(rows))])
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-(len(rows) - 1) - half * 1.6, half * 1.6)
    ax.set_xlim(xs[0], xs[-1])
    return xs[0], xs[-1]


def null_legend_handles(null_label="Null of matched mean"):
    import matplotlib.lines as mlines
    return [mpatches.Patch(color=MATCH_C, label="Matched"),
            mpatches.Patch(color=MISMATCH_C, label="Mismatched"),
            mlines.Line2D([], [], color="black", lw=1.2, marker="|", ms=9, mew=1.2, ls="", label="Observed mean"),
            mlines.Line2D([], [], color="black", lw=1.2, marker="D", ms=6, mfc="white", mew=1.2,
                          label="Adjusted mean (95% CI)"),
            mpatches.Patch(color=NULL_C, label=null_label)]

def null_legend(ax, handles, ax_h_in):
    """Two right-aligned legend rows above the axes: fills and mean markers, then the null."""
    kw = dict(loc="lower right", borderaxespad=0.2, handlelength=1.0, handleheight=0.9,
              columnspacing=0.9, handletextpad=0.4)
    top = ax.legend(handles=handles[:4], ncol=4, bbox_to_anchor=(1.0, 1.0 + 0.27 / ax_h_in), **kw)
    ax.add_artist(top)
    ax.legend(handles=handles[4:], ncol=1, bbox_to_anchor=(1.0, 1.0), **kw)


def check_layout(fig, ax_a, panels, W):
    """Rendered checks, no vector export: text sizes (>= 10 pt, axis titles >= 12 pt), text
    inside the page, no overlapping text boxes, A right-aligned with D, B to D the same size."""
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
    pa, pos = ax_a.get_position(), [q.get_position() for q in panels]
    if abs(pa.x1 - pos[-1].x1) > tol:
        fails.append("panel A and the last scatter panel are not right-aligned")
    if any(abs(q.width - pos[0].width) > tol or abs(q.height - pos[0].height) > tol
           or abs(q.y0 - pos[0].y0) > tol for q in pos[1:]):
        fails.append("scatter panels differ in size or baseline")
    return {"n_text": len(texts), "failures": fails}


def plot_figure(pivot, tab, label_map, bench_label, out_png, corr, null_dir, study="CMC", qa_dir=None):
    """pivot: contrast x author type (log2FC rho); tab: celltype_summary rows;
    corr: correlations.tsv rows (plain and partial Spearman annotated on C and D);
    null_dir: benchmark-post-hoc/de_permutation_null (panel E)."""
    t = tab.set_index("ct_author")
    lab_f1 = t.groupby("pipeline_label")["benchmark_f1"].first()
    cols = sorted(pivot.columns, key=lambda c: (-np.nan_to_num(lab_f1[label_map[c]], nan=-1),
                                                label_map[c], -t.loc[c, "mean_rho"]))
    pv = pivot[cols]
    ncol, nrow = len(cols), len(pv.index)
    # open markers: the seven types pooled under another type's name (IT subtypes, SST Chodl,
    # LAMP5 LHX6), i.e. not benchmark_own_label; SST and LAMP5 keep their own name and stay filled
    pooled = {c for c in cols if not bool(t.loc[c, "benchmark_own_label"])}

    W = 183 * MM
    left, right = 1.12, 0.08
    top = 0.32  # room for the A letter above the bars
    hm_w = W - left - right
    cell_h, bar_h, gap = 0.21, 0.17, 0.04
    hm_h = nrow * cell_h
    y_e, e_h = 0.62, 1.45                     # panel E: matched vs mismatched pairs
    y_b, b_h = y_e + e_h + 1.3, 2.6           # panels B, C and D
    y_key, key_h = y_b + b_h + 0.5, 0.62     # key row: marker legend and colour bars
    y_cls, cls_h = y_key + key_h + 0.05, 0.72  # class legend row
    y_hm = y_cls + cls_h + 0.05 + 1.45        # A tick labels and axis title
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

    def bar(y0, colors, label, outline=False):
        b = fig.add_axes([fx(left), fy(y0), fx(hm_w), fy(bar_h)])
        b.imshow(np.array([colors]), aspect="auto", extent=(0, ncol, 0, 1), interpolation="nearest")
        for i in range(1, ncol):
            b.axvline(i, color=EDGE if outline else "white", lw=0.5)
        if outline:
            b.add_patch(mpatches.Rectangle((0, 0), ncol, 1, fill=False, edgecolor=EDGE, lw=0.5, clip_on=False))
        b.set_xticks([]); b.set_yticks([0.5]); b.set_yticklabels([label])
        b.tick_params(length=0, pad=3)
        for s in b.spines.values():
            s.set_visible(False)

    bar(bars_y0, [mpl.colors.to_rgba(CLASS_COLORS[cell_class(c)]) for c in cols], "Cell class", outline=True)
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
    present = [k for k in CLASS_COLORS if any(cell_class(c) == k for c in cols)]
    cax_ = fig.add_axes([fx(left - 0.9), fy(y_cls), fx(W - left + 0.9 - right), fy(cls_h)])
    cax_.axis("off")
    cax_.legend(handles=[mpatches.Patch(facecolor=CLASS_COLORS[k], edgecolor=EDGE, lw=0.5, label=k) for k in present],
                title="Cell class", title_fontsize=TXT, alignment="left",
                loc="center left", ncol=5, borderaxespad=0, handlelength=1.0, handleheight=0.9,
                columnspacing=0.9, labelspacing=0.3, handletextpad=0.4)
    handles = [mlines.Line2D([], [], ls="", marker="o", mfc="0.55", mec="0.15", ms=6, label="Own label"),
               mlines.Line2D([], [], ls="", marker="o", mfc="white", mec="0.35", mew=1.6, ms=6, label="Pooled label")]
    lax = fig.add_axes([fx(left - 0.9), fy(y_key), fx(4.0), fy(key_h)])
    lax.axis("off")
    lax.legend(handles=handles, loc="center left", ncol=2, borderaxespad=0, handlelength=1.0,
               handleheight=0.9, columnspacing=0.9, labelspacing=0.4, handletextpad=0.4)
    cb_w = 0.95
    for k, (sm, lab) in enumerate([(mesh, "Spearman ρ"),
                                   (mpl.cm.ScalarMappable(f1_norm, f1_cmap), "Annotation F1")]):
        cx0 = W - right - 2 * cb_w - 0.55 + k * (cb_w + 0.3)
        cax = fig.add_axes([fx(cx0), fy(y_key + 0.2), fx(cb_w), fy(0.10)])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=TXT, labelpad=5)
        cax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=TXT, length=2.5, width=0.6)
        cb.set_ticks([-1, 0, 1] if k == 0 else [f1_lo, 1.0])
        cb.outline.set_linewidth(0.6)

    # ---------- B: mean rho vs cells per donor; C, D: added-variable plots (partial Spearman) ----------
    from scipy import stats as sst
    lft, g_bc, g_cd = 0.72, 0.85, 0.42         # B and C carry y-axis titles, D only ticks
    pw = (W - right - lft - g_bc - g_cd) / 3
    xs0 = [lft, lft + pw + g_bc, lft + 2 * pw + g_bc + g_cd]
    y = t.loc[cols, "mean_rho"].values.astype(float)
    xf = t.loc[cols, "cmc_f1"].values.astype(float)
    xz = t.loc[cols, "log10_author_median_cells"].values.astype(float)
    disp = np.array([author_display(c) for c in cols])
    mn = lambda v: f"{v:.2f}".replace("-", "\u2212")

    def ols_band(a, x, yv, tr=lambda v: v):
        """OLS line and 95% confidence band for the mean, drawn on the plotted axes."""
        n = len(x)
        slope, icpt = np.polyfit(x, yv, 1)
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 200)
        resid = yv - (icpt + slope * x)
        s_err = np.sqrt(np.sum(resid ** 2) / (n - 2))
        se = s_err * np.sqrt(1 / n + (xs - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
        tq = sst.t.ppf(0.975, n - 2)
        yfit = icpt + slope * xs
        a.fill_between(tr(xs), yfit - tq * se, yfit + tq * se, color="0.85", lw=0, zorder=1)
        a.plot(tr(xs), yfit, color="0.3", lw=1.2, zorder=2)
        return np.column_stack([tr(xs), yfit])

    def points(a, xv, yv, cs=cols):
        for xi, yi, c in zip(xv, yv, cs):
            colr = CLASS_COLORS[cell_class(c)]
            if c in pooled:
                a.scatter(xi, yi, s=34, facecolors="white", edgecolors=colr, linewidths=1.6, zorder=3)
            else:
                a.scatter(xi, yi, s=34, facecolors=colr, edgecolors="0.15", linewidths=0.6, zorder=3)

    sax = []

    # B: mean rho against median cells per donor (log axis)
    a = fig.add_axes([fx(xs0[0]), fy(y_b), fx(pw), fy(b_h)])
    tr = lambda v: 10 ** v
    line = ols_band(a, xz, y, tr)
    points(a, tr(xz), y)
    a.set_ylim(np.floor(np.nanmin(y) * 10) / 10 - 0.05, 1.05)
    pad = (np.nanmax(xz) - np.nanmin(xz)) * 0.06
    a.set_xscale("log")
    a.set_xlim(10 ** (np.nanmin(xz) - pad), 10 ** (np.nanmax(xz) + pad))
    ticks = [v for v in (1, 10, 100, 1000, 10000) if v <= 10 ** np.nanmax(xz)]
    a.set_xticks(ticks); a.set_xticklabels([f"{v:g}" for v in ticks])
    a.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    lab = [d if yi < 0.65 else "" for d, yi in zip(disp, y)]
    _repel_labels(a, tr(xz), y, lab, fontsize=TXT, avoid_xy=line)
    a.set_xlabel("Median cells per donor + 1")
    a.set_ylabel("Mean Spearman ρ")
    rho, pval = sst.spearmanr(xz, y)
    ptxt = f"p = {pval:.2g}" if pval >= 0.001 else "p < 0.001"
    a.set_title(f"ρ = {rho:.2f}, {ptxt}", loc="left", fontsize=TXT, pad=4)
    sax.append(a)

    # C, D: rank residuals of agreement and of CMC F1, each regressed on rank cells per donor,
    # computed within the subset. Their Pearson r is the partial Spearman rho in correlations.tsv.
    def resid(v, z):
        rz, rv = sst.rankdata(z), sst.rankdata(v)
        return rv - np.polyval(np.polyfit(rz, rv, 1), rz)

    own = np.array([c not in pooled for c in cols])
    for k, (sub, m, n_lab) in enumerate([("all", np.ones(len(cols), bool), 2),
                                         ("excluding types named differently from their pooled label", own, 2)]):
        r_ = corr[corr["subset"] == sub]
        pl = r_[r_["test"] == "rho ~ CMC F1"].iloc[0]
        pa_ = r_[r_["test"] == "rho ~ CMC F1 | log10 support (partial)"].iloc[0]
        bx, by = resid(xf[m], xz[m]), resid(y[m], xz[m])
        assert int(pa_["n"]) == m.sum() and abs(np.corrcoef(bx, by)[0, 1] - pa_["spearman_rho"]) < 1e-9
        a = fig.add_axes([fx(xs0[k + 1]), fy(y_b), fx(pw), fy(b_h)])
        line = ols_band(a, bx, by)
        points(a, bx, by, [c for c, mi in zip(cols, m) if mi])
        for v, setl in ((bx, a.set_xlim), (by, a.set_ylim)):
            pad = (v.max() - v.min()) * 0.06
            setl(v.min() - pad, v.max() + pad)
        a.axhline(0, color="0.75", lw=0.6, zorder=0); a.axvline(0, color="0.75", lw=0.6, zorder=0)
        low = set(np.argsort(by)[:n_lab])          # label the lowest adjusted agreement
        _repel_labels(a, bx, by, [d if i in low else "" for i, d in enumerate(disp[m])],
                      fontsize=TXT, avoid_xy=line)
        a.set_title(f"{int(m.sum())} types: ρ {mn(pl['spearman_rho'])} → {mn(pa_['spearman_rho'])}",
                    loc="left", fontsize=TXT, pad=4)
        a.set_xlabel(f"{study} F1 rank, adj.")
        if k == 0:
            a.set_ylabel("Mean ρ rank, adj.")
        sax.append(a)

    fig.text(fx(0.03), fy(H - 0.03), "A", fontsize=LET, fontweight="bold", va="top")
    for k, x0 in enumerate(xs0):
        fig.text(fx(0.03 if k == 0 else x0 - (0.78 if k == 1 else 0.4)), fy(y_b + b_h + 0.42), "BCD"[k],
                 fontsize=LET, fontweight="bold", va="top")

    # ---------- E: per-pair rho, matched vs mismatched cell types, with the permutation null ----------
    ea = fig.add_axes([fx(left), fy(y_e), fx(hm_w - 0.15), fy(e_h)])  # room for the 1.0 tick label
    rows = load_null(null_dir, "labels", ["subclass"])
    ea.set_xticks(np.arange(-1, 1.01, 0.5))
    ea.set_xlabel("Spearman ρ per cell type and contrast pair")
    lo, hi = null_density_rows(ea, rows, ["Author vs\npipeline"])
    s0 = rows[0][3]
    pe = f"p = {s0['empirical_p']:.2g}" if s0["empirical_p"] >= 0.001 else "p < 0.001"
    h = null_legend_handles(f"Null of matched mean ({s0['null_mean']:.2f} ± {s0['null_sd']:.2f}, {pe})")
    null_legend(ea, h, e_h)
    fig.text(fx(0.03), fy(y_e + e_h + 0.65), "E", fontsize=LET, fontweight="bold", va="top")

    qa = check_layout(fig, ax, sax, W)
    if qa_dir:
        import json
        os.makedirs(qa_dir, exist_ok=True)
        json.dump(qa, open(os.path.join(qa_dir, "layout_qa.json"), "w"), indent=2)
    if qa["failures"]:
        raise RuntimeError("Figure QA failed: " + "; ".join(qa["failures"]))
    fig.savefig(out_png, dpi=600, facecolor="white")
    plt.close(fig)
