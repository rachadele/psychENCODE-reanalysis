"""Manuscript figure: DE agreement vs annotation F1 and cell support (sc-manuscript Figure 4).

At final size (183 mm, double column). A: heatmap of DE agreement (Spearman rho of log2 fold
changes, or of limma moderated t with merged_figure(stat="t"), author vs pipeline labels) per contrast and author type, with bars for cohort F1
(PsychENCODE donors) of the subclass and class labels and cell class, and brackets over author types pooled into one pipeline label.
B: Spearman rho of mean agreement with log10 median cells per donor (27 types) and with F1
(27 types; the one-to-one types, by default the 18 whose pipeline label pools no other type), plain
(filled) and partial controlling for the other variable (open; each row label names the variable
controlled for), 95% CIs by Fisher z, p values on the right (fig4b_rho_ci.tsv). F1 (plain mean of
per-sample F1) and cells per donor (median, zeros excluded) come from every donor in the DE of the
type (celltype_f1_cells.tsv); the cohort F1 bar in A is the same F1. C: density of per-pair rho for matched author/pipeline cell types (above) and
mismatched types in the same contrast (below), observed means as ticks, means adjusted for cells per
donor and contrast as diamonds with 95% CI, permutation null of the matched mean in red
(benchmark-post-hoc/de_permutation_null). plot_scatter() draws the per-cell-type scatters behind B
as a separate supplementary PNG. Text: >= 10 pt, axis titles >= 12 pt, letters 14 pt.
The only output is the PNG; check_layout() verifies text sizes, page bounds, text overlaps
and panel alignment on the rendered figure and blocks export on failure (report to qa_dir).
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

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
    inside the page, no overlapping text boxes."""
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
    return {"n_text": len(texts), "failures": fails}


def _panel_b_groups(ba, groups, compact=False):
    """Panel B as grouped rows: a bold heading per variable, then one row per type set, plain (filled)
    and partial (open, controlling for the other variable) rho with 95% CI, p plain -> partial on the right."""
    rows, ys, heads, y = [], [], [], 0
    for head, items in groups[::-1]:          # bottom-up
        for name, r in items[::-1]:
            rows.append((name, r)); ys.append(y); y += 1
        heads.append((head, y - 0.25)); y += 0.7
    off = 0.15
    tr = mpl.transforms.blended_transform_factory(ba.transAxes, ba.transData)
    for (name, r), yy in zip(rows, ys):
        for kind, dy, mfc in (("plain", off, "0.2"), ("partial", -off, "white")):
            ba.plot([r[f"{kind}_ci_low"], r[f"{kind}_ci_high"]], [yy + dy] * 2, color="0.2", lw=1.2, zorder=2,
                    solid_capstyle="butt")
            ba.plot([r[f"{kind}_rho"]], [yy + dy], marker="o", ms=7, mfc=mfc, mec="0.2", mew=1.4, ls="", zorder=4)
        ba.text(1.032, yy, f"{fmt_p(r['plain_p'])} → {fmt_p(r['partial_p'])}", transform=tr,
                ha="left", va="center", fontsize=TXT)
    for head, yy in heads:
        ba.text(-0.02, yy, head, transform=tr, ha="right", va="bottom", fontsize=TXT, fontweight="bold")
    ba.text(1.032, max(ys) + 0.45, "p", transform=tr, ha="left", va="bottom", fontsize=TXT)
    ba.axvline(0, color="0.5", lw=0.8, ls="--", zorder=1)
    ba.set_xlim(-1, 1); ba.set_xticks([-1, 0, 1] if compact else [-1, -0.5, 0, 0.5, 1])
    ba.set_ylim(-0.55, y - 0.2)
    ba.set_yticks(ys); ba.set_yticklabels([n for n, _ in rows])
    ba.tick_params(axis="y", length=0)
    ba.spines["left"].set_visible(False)
    ba.set_xlabel("ρ with mean agreement" if compact else "Spearman ρ with mean agreement")
    hs = [mlines.Line2D([], [], ls="", marker="o", ms=7, mfc="0.2", mec="0.2", label="Plain"),
          mlines.Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec="0.2", mew=1.4,
                        label="Partial, other variable controlled" if compact else
                        "Partial, controlling for the other variable"),
          mlines.Line2D([], [], color="0.2", lw=1.2, label="95% CI")]
    if compact:  # narrow panel: two rows, left-aligned with the row headings
        ba.legend(handles=[hs[0], hs[1], hs[2]], loc="lower left", bbox_to_anchor=(-0.8, 1.0), ncol=2,
                  borderaxespad=0.2, handlelength=1.4, columnspacing=0.9, handletextpad=0.4, labelspacing=0.3)
    else:
        ba.legend(handles=hs, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, borderaxespad=0.2,
                  handlelength=1.4, columnspacing=0.9, handletextpad=0.4)


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.2g}"


# ---------- D, E: author subtypes vs pooled author types (subtype_vs_pooled/) ----------
POOL_COMPARISONS = [("subtype_vs_pooled", "Subtype vs\npooled"), ("subtype_vs_subtype", "Subtype vs\nsubtype"),
                    ("nearest", "vs nearest\nseparate"), ("vs_glia", "vs glia")]
POOL_PUBLISHED = [("published_vs_author_subtype", "Our\nsubtype"), ("published_vs_merged", "Our pooled\ntype"),
                  ("published_vs_ours", "Our\nlabel")]
POOL_GROUPS = {"IT": ("L2/3-6 IT", "o", "IT (5 subtypes)"), "Lamp5": ("LAMP5", "s", "Lamp5 (2)"),
               "Sst": ("SST", "D", "Sst (2)")}
# nearest separate group of author types: deep-layer non-IT for IT, other GABAergic for Lamp5 and Sst
POOL_NEAREST = {"IT": "vs_deep_non_it", "Lamp5": "vs_gaba", "Sst": "vs_gaba"}


def pooling_data(pooling_dir):
    """Per-contrast rho of t (mean over pairs within a contrast) and the mean over contrasts per group and
    comparison, from pooling_dir/pairwise_rho.tsv, checked against summary.tsv (as
    subtype_vs_pooled/plot_subtype_vs_pooled_t.py). Category order: highest to lowest mean of group means."""
    pair = pd.read_csv(os.path.join(pooling_dir, "pairwise_rho.tsv"), sep="\t")
    summ = pd.read_csv(os.path.join(pooling_dir, "summary.tsv"), sep="\t")
    per_c = pair.groupby(["group", "comparison", "contrast"])["rho_t"].mean().reset_index()
    means = per_c.groupby(["group", "comparison"])["rho_t"].mean()
    sg = summ[summ["level"] == "group"].set_index(["group", "comparison"])["mean_rho_t"]
    for k, v in means.items():
        assert abs(v - sg[k]) < 1e-3, (k, v, sg[k])
    for g, c in POOL_NEAREST.items():
        means[(g, "nearest")] = means[(g, c)]
    per_c = pd.concat([per_c] + [per_c[(per_c["group"] == g) & (per_c["comparison"] == c)].assign(comparison="nearest")
                                 for g, c in POOL_NEAREST.items()], ignore_index=True)
    om = {c: np.mean([means[(g, c)] for g in POOL_GROUPS]) for c, _ in POOL_COMPARISONS + POOL_PUBLISHED}
    return (per_c, means, sorted(POOL_COMPARISONS, key=lambda c: -om[c[0]]),
            sorted(POOL_PUBLISHED, key=lambda c: -om[c[0]]))


def draw_pooling(axD, axE, pooling_dir):
    per_c, means, comps, pubs = pooling_data(pooling_dir)
    dodge = {"IT": -0.25, "Lamp5": 0.0, "Sst": 0.25}
    rng = np.random.default_rng(0)
    for ax, cats in ((axD, comps), (axE, pubs)):
        for i, (comp, _) in enumerate(cats):
            for g, (cls, mk, _) in POOL_GROUPS.items():
                v = per_c[(per_c["group"] == g) & (per_c["comparison"] == comp)]["rho_t"].to_numpy()
                x = i + dodge[g]
                ax.scatter(x + rng.uniform(-0.06, 0.06, len(v)), v, s=16, color=CLASS_COLORS[cls], marker=mk,
                           alpha=0.75, edgecolor="0.3", linewidth=0.3, zorder=2)
                ax.plot([x - 0.11, x + 0.11], [means[(g, comp)]] * 2, color="black", lw=1.8, zorder=3)
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels([lab for _, lab in cats])
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.5, len(cats) - 0.5)
        for i in range(1, len(cats)):
            ax.axvline(i - 0.5, color="0.88", lw=0.6, zorder=0)
        ax.set_ylim(-0.05, 1.03)
        ax.set_yticks([0, 0.5, 1])
    axD.set_ylabel("Spearman ρ of t")
    axE.tick_params(axis="y", labelleft=False)
    axD.set_title("Our DE on author subtypes\n(8 contrasts)", fontsize=TXT, pad=4)
    axE.set_title("Published subtype DE\nvs our DE", fontsize=TXT, pad=4)
    handles = [mlines.Line2D([], [], ls="", marker=mk, color=CLASS_COLORS[cls], mec="0.3", ms=6, label=lab)
               for cls, mk, lab in POOL_GROUPS.values()]
    handles.append(mlines.Line2D([], [], color="black", lw=1.8, label="Mean over contrasts"))
    return handles


def plot_figure(pivot, tab, label_map, bench_label, out_png, rho_ci, null_dir, study="CMC", qa_dir=None,
                class_f1_tsv=None, display=None, class_of=None, brackets=True, b_groups=None, pooling_dir=None):
    """pivot: contrast x author type (log2FC rho); tab: celltype_summary rows;
    class_f1_tsv: class-level cohort F1 per author type (class_f1_de_donors.py output); default: env
    FIG4_CLASS_F1_TSV, else celltype_f1_cells_class.tsv beside out_png;
    rho_ci: plain and partial Spearman rho with CIs, one row per test (panel B; fig4b_rho_ci.tsv);
    null_dir: benchmark-post-hoc/de_permutation_null (panel C).
    Merged-author build (merged_figure()): display / class_of map column keys to tick labels and cell
    classes (default author_display, cell_class); brackets=False drops the pooled-type bracket strip;
    b_groups = [(heading, [(row label, rho_ci row), ...]), ...] draws panel B as grouped rows.
    pooling_dir: subtype_vs_pooled/ folder; adds panels D (author subtypes vs pooled types) and E
    (published subtype DE vs our DE) below C."""
    disp_of = (lambda c: display[c]) if display else author_display
    class_fn = (lambda c: class_of[c]) if class_of else cell_class
    t = tab.set_index("ct_author")
    class_f1_tsv = class_f1_tsv or os.environ.get("FIG4_CLASS_F1_TSV") or \
        os.path.join(os.path.dirname(os.path.abspath(out_png)), "celltype_f1_cells_class.tsv")
    cf = pd.read_csv(class_f1_tsv, sep="\t")
    cf = cf[cf["scheme"] == "author"].set_index("celltype")["f1"]
    t["class_f1"] = cf.reindex(t.index)
    assert t["class_f1"].notna().all(), "author types missing from " + class_f1_tsv
    # order: mean cohort F1 of the author types sharing a pipeline label (descending, so pooled
    # types stay adjacent under their bracket), then pipeline label, then cohort F1 (descending).
    # Same rule as figure4_order() in psychENCODE-author-de-comparison compare_de_results.py (S7).
    lab_f1 = t.groupby("pipeline_label")["f1"].mean()
    cols = sorted(pivot.columns, key=lambda c: (-np.nan_to_num(lab_f1[label_map[c]], nan=-1),
                                                label_map[c], -t.loc[c, "f1"]))
    pv = pivot[cols]
    ncol, nrow = len(cols), len(pv.index)

    W = 183 * MM
    left, right = 1.12, 0.08
    top = 0.32  # room for the A letter above the bars
    hm_w = W - left - right
    cell_h, bar_h, gap = 0.21, 0.17, 0.04
    hm_h = nrow * cell_h
    y_d, d_h = 0.78, 1.9                      # panels D, E: subtypes vs pooled types (optional)
    y_e, e_h = (y_d + d_h + 1.55 if pooling_dir else 0.62), 1.45  # panel C: matched vs mismatched pairs
    y_b, b_h = y_e + e_h + 1.3, 1.45          # panel B: plain and partial rho
    b_left, p_col = 2.2, 1.3                 # B row labels on the left, p values on the right
    side = bool(b_groups)                    # grouped B sits beside C on one row (left B, right C)
    if b_groups:
        b_h = 0.7 + 0.5 * sum(len(rows) for _, rows in b_groups)  # grouped rows need more height
        y_b = y_e                            # same row as C
        b_left, p_col = 1.35, 1.15
        b_w = 1.5                            # B data width; p values to its right
        c_x0 = b_left + b_w + p_col + 0.3     # C starts right of B's p values
        e_h = 1.5                            # C axes at the bottom of the row, legend above
    y_key, key_h = y_b + b_h + (0.75 if side else 0.45), 0.62    # key row: colour bars
    y_cls, cls_h = y_key + key_h + 0.05, 0.72  # class legend row
    y_hm = y_cls + cls_h + 0.05 + 1.45        # A tick labels and axis title
    y_br = y_hm + hm_h + gap
    br_h = 0.30 if brackets else 0.0
    bars_y0 = y_br + br_h + gap if brackets else y_br
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
    ax.set_xticklabels([disp_of(c) for c in cols], rotation=90, ha="right", va="center",
                       rotation_mode="anchor")
    ax.set_yticks(np.arange(nrow) + 0.5); ax.set_yticklabels(pv.index)
    ax.tick_params(length=0, pad=3)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_ylabel("Contrast")
    ax.set_xlabel("Author cell type", labelpad=4)

    f1_lo = min(0.3, np.floor(np.nanmin(t[["f1", "class_f1"]].values) * 10) / 10)
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

    bar(bars_y0, [mpl.colors.to_rgba(CLASS_COLORS[class_fn(c)]) for c in cols], "Cell class", outline=True)
    bar(bars_y0 + (bar_h + gap), [f1c(t.loc[c, "f1"]) for c in cols], "Subclass F1")
    bar(bars_y0 + 2 * (bar_h + gap), [f1c(t.loc[c, "class_f1"]) for c in cols], "Class F1")

    # brackets over author types pooled into one pipeline label
    if brackets:
        bax = fig.add_axes([fx(left), fy(y_br), fx(hm_w), fy(br_h)])
        bax.set_xlim(0, ncol); bax.set_ylim(0, 1); bax.axis("off")
    i = 0
    while brackets and i < ncol:
        j = i
        while j + 1 < ncol and label_map[cols[j + 1]] == label_map[cols[i]]:
            j += 1
        if j > i:
            bax.plot([i + 0.15, i + 0.15, j + 0.85, j + 0.85], [0.05, 0.25, 0.25, 0.05], color="0.25", lw=0.8)
            bax.text((i + j + 1) / 2, 0.32, bench_label.get(label_map[cols[i]], ""), ha="center", va="bottom")
        i = j + 1

    # ---------- key row: class legend, colour bars ----------
    present = [k for k in CLASS_COLORS if any(class_fn(c) == k for c in cols)]
    cax_ = fig.add_axes([fx(left - 0.9), fy(y_cls), fx(W - left + 0.9 - right), fy(cls_h)])
    cax_.axis("off")
    cax_.legend(handles=[mpatches.Patch(facecolor=CLASS_COLORS[k], edgecolor=EDGE, lw=0.5, label=k) for k in present],
                title="Cell class", title_fontsize=TXT, alignment="left",
                loc="center left", ncol=5, borderaxespad=0, handlelength=1.0, handleheight=0.9,
                columnspacing=0.9, labelspacing=0.3, handletextpad=0.4)
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

    # ---------- B: plain and partial Spearman rho of agreement, with 95% CIs ----------
    ba = fig.add_axes([fx(b_left), fy(y_b), fx(b_w if side else W - right - b_left - p_col), fy(b_h)])
    if b_groups:
        _panel_b_groups(ba, b_groups, compact=side)
    else:
        rows_b = list(rho_ci.itertuples())
        n_r = len(rows_b)
        off = 0.14
        for i, r in enumerate(rows_b):
            yy = n_r - 1 - i
            for kind, dy, mfc in (("plain", off, "0.2"), ("partial", -off, "white")):
                rho, lo_, hi_ = getattr(r, f"{kind}_rho"), getattr(r, f"{kind}_ci_low"), getattr(r, f"{kind}_ci_high")
                ba.plot([lo_, hi_], [yy + dy] * 2, color="0.2", lw=1.2, zorder=2, solid_capstyle="butt")
                ba.plot([rho], [yy + dy], marker="o", ms=7, mfc=mfc, mec="0.2", mew=1.4, ls="", zorder=4)
            ba.text(1.0 + 0.02 * 1.6, yy, f"{fmt_p(r.plain_p)} → {fmt_p(r.partial_p)}",
                    transform=mpl.transforms.blended_transform_factory(ba.transAxes, ba.transData),
                    ha="left", va="center", fontsize=TXT)
        ba.axvline(0, color="0.5", lw=0.8, ls="--", zorder=1)
        ba.set_xlim(-1, 1)
        ba.set_xticks([-1, -0.5, 0, 0.5, 1])
        ba.set_ylim(-0.6, n_r - 0.4)
        ba.set_yticks(range(n_r))
        # each row names the variable its open (partial) point controls for
        ctrl = lambda v: "cells per donor" if "cells" in str(v) else str(v)
        ba.set_yticklabels([f"{r.label}\nopen: given {ctrl(r.controlled_for)}" for r in rows_b][::-1])
        ba.tick_params(axis="y", length=0)
        ba.spines["left"].set_visible(False)
        ba.set_xlabel("Spearman ρ with mean agreement")
        ba.text(1.0 + 0.02 * 1.6, n_r - 0.4, "p", transform=mpl.transforms.blended_transform_factory(ba.transAxes, ba.transData),
                ha="left", va="bottom", fontsize=TXT)
        ba.legend(handles=[mlines.Line2D([], [], ls="", marker="o", ms=7, mfc="0.2", mec="0.2", label="Plain"),
                           mlines.Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec="0.2", mew=1.4,
                                         label="Partial (open), given the variable named on the row"),
                           mlines.Line2D([], [], color="0.2", lw=1.2, label="95% CI")],
                  loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, borderaxespad=0.2, handlelength=1.4,
                  columnspacing=0.9, handletextpad=0.4)
    fig.text(fx(0.03), fy(H - 0.03), "A", fontsize=LET, fontweight="bold", va="top")
    fig.text(fx(0.03), fy(y_b + b_h + (0.85 if side else 0.62)), "B", fontsize=LET, fontweight="bold", va="top")

    # ---------- C: per-pair rho, matched vs mismatched cell types, with the permutation null ----------
    if side:  # one row (author vs pipeline), no row label; legend stacked above the axes
        ea = fig.add_axes([fx(c_x0), fy(y_e), fx(W - right - 0.12 - c_x0), fy(e_h)])
    else:
        ea = fig.add_axes([fx(left), fy(y_e), fx(hm_w - 0.15), fy(e_h)])  # room for the 1.0 tick label
    rows = load_null(null_dir, "labels", ["subclass"])
    ea.set_xticks(np.arange(-1, 1.01, 0.5))
    ea.set_xlabel("Spearman ρ per pair" if side else "Spearman ρ per cell type and contrast pair")
    lo, hi = null_density_rows(ea, rows, [""] if side else ["Author vs\npipeline"])
    s0 = rows[0][3]
    pe = f"p = {s0['empirical_p']:.2g}" if s0["empirical_p"] >= 0.001 else "p < 0.001"
    if side:
        h = null_legend_handles(f"Null ({s0['null_mean']:.2f} ± {s0['null_sd']:.2f}, {pe})")
        ea.legend(handles=h, ncol=1, loc="lower left", bbox_to_anchor=(0.0, 1.02),
                  borderaxespad=0, handlelength=1.0, handleheight=0.9, columnspacing=0.8, handletextpad=0.4,
                  labelspacing=0.35)
        fig.text(fx(c_x0 - 0.3), fy(y_b + b_h + 0.85), "C", fontsize=LET, fontweight="bold", va="top")
    else:
        h = null_legend_handles(f"Null of matched mean ({s0['null_mean']:.2f} ± {s0['null_sd']:.2f}, {pe})")
        null_legend(ea, h, e_h)
        fig.text(fx(0.03), fy(y_e + e_h + 0.65), "C", fontsize=LET, fontweight="bold", va="top")

    extra = []
    if pooling_dir:
        d_w = 3.5
        axD = fig.add_axes([fx(left), fy(y_d), fx(d_w), fy(d_h)])
        e_x0 = left + d_w + 0.3
        axE = fig.add_axes([fx(e_x0), fy(y_d), fx(W - right - e_x0), fy(d_h)])
        hd = draw_pooling(axD, axE, pooling_dir)
        fig.legend(handles=hd, loc="lower center", bbox_to_anchor=(fx(W / 2), fy(y_d + d_h + 0.5)), ncol=4,
                   borderaxespad=0, handlelength=1.4, columnspacing=1.2, handletextpad=0.4)
        fig.text(fx(0.03), fy(y_d + d_h + 0.45), "D", fontsize=LET, fontweight="bold", va="top")
        fig.text(fx(e_x0 - 0.25), fy(y_d + d_h + 0.45), "E", fontsize=LET, fontweight="bold", va="top")
        extra = [axD, axE]
    qa = check_layout(fig, ax, [ba, ea, *extra], W)
    if qa_dir:
        import json
        os.makedirs(qa_dir, exist_ok=True)
        json.dump(qa, open(os.path.join(qa_dir, "layout_qa.json"), "w"), indent=2)
    if qa["failures"]:
        raise RuntimeError("Figure QA failed: " + "; ".join(qa["failures"]))
    fig.savefig(out_png, dpi=600, facecolor="white")
    plt.close(fig)


def plot_scatter(tab, cols, corr, out_png):
    """Supplementary (not in the manuscript figure set): per-cell-type scatters behind Figure 4B.
    Left: mean rho against median cells per donor (log axis). Middle, right: added-variable plots
    for the partial Spearman of mean rho with F1 given cells per donor, for all types and for
    the one-to-one types (rank residuals computed within the subset). F1 and cells per donor are the
    tab columns f1 and log10_cells (celltype_f1_cells.tsv); corr is the per-cell-type rows of
    correlations_per_cohort.tsv. Points coloured by cell class, open markers for types that are not
    one-to-one (their pipeline label pools other author types)."""
    t = tab.set_index("ct_author")
    pooled = {c for c in cols if not bool(t.loc[c, "one_to_one"])}
    W, right = 183 * MM, 0.08
    y_b, b_h = 0.62, 2.6
    y_leg, leg_h = y_b + b_h + 0.5, 0.95
    H = y_leg + leg_h + 0.05
    fig = plt.figure(figsize=(W, H))
    fx, fy = (lambda v: v / W), (lambda v: v / H)
    present = [k for k in CLASS_COLORS if any(cell_class(c) == k for c in cols)]
    lg = fig.add_axes([fx(0.05), fy(y_leg), fx(W - 0.1), fy(leg_h)])
    lg.axis("off")
    lg.legend(handles=[mpatches.Patch(facecolor=CLASS_COLORS[k], edgecolor=EDGE, lw=0.5, label=k) for k in present]
              + [mlines.Line2D([], [], ls="", marker="o", mfc="white", mec="0.35", mew=1.6, ms=6,
                               label="Pipeline label pools types")],
              loc="center left", ncol=4, borderaxespad=0, handlelength=1.0, handleheight=0.9,
              columnspacing=0.9, labelspacing=0.3, handletextpad=0.4)
    # ---------- B: mean rho vs cells per donor; C, D: added-variable plots (partial Spearman) ----------
    from scipy import stats as sst
    lft, g_bc, g_cd = 0.72, 0.85, 0.42         # B and C carry y-axis titles, D only ticks
    pw = (W - right - lft - g_bc - g_cd) / 3
    xs0 = [lft, lft + pw + g_bc, lft + 2 * pw + g_bc + g_cd]
    y = t.loc[cols, "mean_rho"].values.astype(float)
    xf = t.loc[cols, "f1"].values.astype(float)
    xz = t.loc[cols, "log10_cells"].values.astype(float)
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

    # C, D: rank residuals of agreement and of F1, each regressed on rank cells per donor,
    # computed within the subset. Their Pearson r is the partial Spearman rho in correlations.tsv.
    def resid(v, z):
        rz, rv = sst.rankdata(z), sst.rankdata(v)
        return rv - np.polyval(np.polyfit(rz, rv, 1), rz)

    own = np.array([c not in pooled for c in cols])
    for k, (sub, m, n_lab) in enumerate([("all", np.ones(len(cols), bool), 2),
                                         ("one-to-one types", own, 2)]):
        r_ = corr[corr["subset"] == sub]
        pl = r_[r_["test"] == "agreement ~ F1"].iloc[0]
        pa_ = r_[r_["test"] == "agreement ~ F1 | log10 cells per donor"].iloc[0]
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
        a.set_xlabel("F1 rank, adj.")
        if k == 0:
            a.set_ylabel("Mean ρ rank, adj.")
        sax.append(a)

    for k, x0 in enumerate(xs0):
        fig.text(fx(0.03 if k == 0 else x0 - (0.78 if k == 1 else 0.4)), fy(y_b + b_h + 0.42), "ABC"[k],
                 fontsize=LET, fontweight="bold", va="top")
    qa = check_layout(fig, sax[0], sax, W)
    if qa["failures"]:
        raise RuntimeError("Scatter figure QA failed: " + "; ".join(qa["failures"]))
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)


# merged_figure(stat=...): per-pair agreement column (combined_filtered_corr.tsv), per-type mean column
# (celltype_merged.tsv), permutation-null folder and output suffix
MERGED_STATS = {"log2FoldChange": ("spearman_log2FoldChange", "mean_rho_log2FC", "null", ""),
                "t": ("spearman_stat", "mean_rho_t", "null_t", "_t")}


def merged_figure(merged_dir, current_rho_ci, out_png, qa_dir=None, stat="log2FoldChange", pooling_dir=None):
    """Figure 4 from the merged-author arm (DE_AGREEMENT_DECISIONS.md, 2026-11-09): author types the
    pipeline pools were summed into one pseudobulk per donor, 21 types each matching one pipeline label.
    Inputs, all under merged_dir (results_v2.1.2_limma/merged_author_2026-11-09):
      comparison/subclass/average_corr/combined_filtered_corr.tsv  rho per type x contrast (A)
      de_agreement/celltype_merged.tsv                              F1, cells per donor, mean rho (A, B)
      fig4_inputs/celltype_f1_cells_class_merged.tsv                class F1 (A)
      fig4_inputs/null/                                             permutation null, adjusted means (C)
    (fig4_inputs from bin/merged_author_fig4_inputs.py). current_rho_ci: fig4b_rho_ci.tsv of the
    27-type analysis, whose two 27-type rows sit beside the 21-type rows in B. stat: log2FoldChange
    (default) or t (limma moderated t, column spearman_stat; null from fig4_inputs/null_t). The 27-type
    rows are recomputed from the current comparison table with the same statistic (for log2FoldChange
    they must equal current_rho_ci). Writes fig4b_rho_ci_merged<suffix>.tsv and
    fig4_agreement_summary<suffix>.tsv (mean +- SD per pair: 21 merged types, 27 author types, class
    level) to fig4_inputs/, and manifest_fig4<suffix>.json beside them."""
    import json
    from utils import CONTRAST_DISPLAY, PIPELINE_TO_BENCH
    from de_agreement_vs_annotation import rho_ci_row
    import datetime, hashlib
    corr_col, mean_col, null_sub, suffix = MERGED_STATS[stat]
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    amap = json.load(open(os.path.join(root, "assets", "author_to_pipeline_map.subclass.json")))
    members = {}
    for t_, lab in amap.items():
        members.setdefault(lab, []).append(t_)
    mc = pd.read_csv(os.path.join(merged_dir, "de_agreement", "celltype_merged.tsv"), sep="\t")
    comb = pd.read_csv(os.path.join(merged_dir, "comparison", "subclass", "average_corr", "combined_filtered_corr.tsv"), sep="\t")
    comb = comb[(comb["ct_author"] == comb["ct_sc_pipeline"]) & comb["contrast"].isin(CONTRAST_DISPLAY)]
    comb = comb.assign(c=comb["contrast"].map(CONTRAST_DISPLAY))
    pivot = comb.pivot_table(index="c", columns="ct_author", values=corr_col)
    pivot = pivot.loc[[c for c in ["SCZ", "BD", "ASD", "MDD", "PTSD", "WS", "Age", "Sex"] if c in pivot.index]]
    tab = pd.DataFrame({"ct_author": mc["celltype"], "pipeline_label": mc["celltype"], "f1": mc["f1"],
                        "mean_rho": mc[mean_col], "log10_cells": mc["log10_cells"]})
    assert set(tab["ct_author"]) == set(pivot.columns) == set(members)
    display = {lab: (author_display(ms[0]) if len(ms) == 1 else PIPELINE_TO_BENCH["subclass"][lab])
               for lab, ms in members.items()}
    assert all(len({CLASS_OF[m] for m in ms}) == 1 for ms in members.values())
    class_of = {lab: CLASS_OF[ms[0]] for lab, ms in members.items()}
    n = len(tab)
    rho_ci = pd.DataFrame([
        rho_ci_row(f"Cells per donor, {n} types", tab["log10_cells"], tab["mean_rho"], tab["f1"], "F1"),
        rho_ci_row(f"F1, {n} types", tab["f1"], tab["mean_rho"], tab["log10_cells"], "log10 cells per donor")])
    assert np.allclose(tab.set_index("ct_author")["mean_rho"], pivot.mean().reindex(tab["ct_author"]).values)
    # 27 author types without merging, same statistic: per-type mean over the eight contrasts
    rdir = os.path.dirname(os.path.dirname(os.path.abspath(merged_dir)) + "/")
    cur_corr = os.path.join(rdir, "comparison", "subclass", "average_corr", "combined", "combined_filtered_corr.tsv")
    cls_corr = os.path.join(rdir, "comparison", "class", "average_corr", "combined", "combined_filtered_corr.tsv")
    cur_sum = os.path.join(os.path.dirname(os.path.abspath(current_rho_ci)), "celltype_summary.tsv")
    cc = pd.read_csv(cur_corr, sep="\t")
    cc = cc[cc["contrast"].isin(CONTRAST_DISPLAY)]
    cs = pd.read_csv(cur_sum, sep="\t").set_index("ct_author")
    cs["mean_rho"] = cc.groupby("ct_author")[corr_col].mean().reindex(cs.index)
    assert cs["mean_rho"].notna().all() and len(cs) == 27
    n27 = len(cs)
    c27 = pd.DataFrame([
        rho_ci_row(f"Cells per donor, {n27} types", cs["log10_cells"], cs["mean_rho"], cs["f1"], "F1"),
        rho_ci_row(f"F1, {n27} types", cs["f1"], cs["mean_rho"], cs["log10_cells"], "log10 cells per donor")])
    if stat == "log2FoldChange":  # must reproduce the 27-type rows of the current analysis
        cur = pd.read_csv(current_rho_ci, sep="\t")
        cur = cur[cur["label"].str.contains("27 types")].set_index("label").loc[c27["label"]]
        num = [k for k in c27.columns if k not in ("label", "controlled_for")]
        assert np.allclose(cur[num].values.astype(float), c27[num].values.astype(float), atol=1e-5)
    fin = os.path.join(merged_dir, "fig4_inputs")
    pd.concat([rho_ci, c27]).to_csv(os.path.join(fin, f"fig4b_rho_ci_merged{suffix}.tsv"), sep="\t",
                                    index=False, float_format="%.6g")
    cl = pd.read_csv(cls_corr, sep="\t")
    cl = cl[cl["contrast"].isin(CONTRAST_DISPLAY)]
    summ = pd.DataFrame([{"set": s_, "level": lv, "n_types": nt, "n_pairs": len(v), "mean_rho": v.mean(),
                          "sd_rho": v.std(ddof=1)} for s_, lv, nt, v in
                         [(f"{n} merged types", "subclass", n, pd.Series(pivot.values.ravel()).dropna()),
                          (f"{n27} author types", "subclass", n27, cc[corr_col]),
                          (f"{n27} author types", "class", cl["ct_author"].nunique(), cl[corr_col])]])
    summ.insert(0, "statistic", stat)
    summ_path = os.path.join(fin, f"fig4_agreement_summary{suffix}.tsv")
    summ.to_csv(summ_path, sep="\t", index=False, float_format="%.6g")
    pick = lambda d, s: d[d["label"].str.startswith(s)].iloc[0]
    groups = [("Cells per donor", [("27 author types", pick(c27, "Cells per donor")),
                                   (f"{n} merged types", pick(rho_ci, "Cells per donor"))]),
              ("F1", [("27 author types", pick(c27, "F1")), (f"{n} merged types", pick(rho_ci, "F1"))])]
    plot_figure(pivot, tab, {lab: lab for lab in members}, PIPELINE_TO_BENCH["subclass"], out_png, rho_ci,
                os.path.join(fin, null_sub), qa_dir=qa_dir,
                class_f1_tsv=os.path.join(fin, "celltype_f1_cells_class_merged.tsv"),
                display=display, class_of=class_of, brackets=False, b_groups=groups, pooling_dir=pooling_dir)

    def sha(f):
        h = hashlib.sha256()
        with open(f, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                h.update(b)
        return h.hexdigest()
    rec = lambda f: {"path": os.path.abspath(f), "sha256": sha(f), "mtime": datetime.datetime.fromtimestamp(
        os.path.getmtime(f)).astimezone().isoformat(timespec="seconds")}
    nd = os.path.join(fin, null_sub)
    ins = [os.path.join(merged_dir, "comparison", "subclass", "average_corr", "combined_filtered_corr.tsv"),
           os.path.join(merged_dir, "de_agreement", "celltype_merged.tsv"),
           os.path.join(fin, "celltype_f1_cells_class_merged.tsv"), current_rho_ci, cur_corr, cls_corr, cur_sum,
           os.path.join(root, "assets", "author_to_pipeline_map.subclass.json"),
           os.path.join(nd, "pairs.tsv"), os.path.join(nd, "null_summary.tsv"),
           os.path.join(nd, "matched_effect_adjusted.tsv"), os.path.join(nd, "manifest.json"), PALETTE_TSV]
    if pooling_dir:  # panels D, E
        ins += [os.path.join(pooling_dir, f) for f in ("pairwise_rho.tsv", "summary.tsv", "manifest.json")]
    json.dump({"created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
               "statistic": f"Spearman rho of {stat} ({corr_col})",
               "data_version": "merged_author_2026-11-09 (author arm, 21 merged types) vs results_v2.1.2_limma "
                               "pipeline-subclass DE; sc-annotation-pipeline 2.1.2 labels; limma-voom",
               "script": rec(os.path.abspath(__file__)),
               "inputs": [rec(f) for f in ins],
               "outputs": [rec(f) for f in (out_png, os.path.join(fin, f"fig4b_rho_ci_merged{suffix}.tsv"), summ_path)]},
              open(os.path.join(fin, f"manifest_fig4{suffix}.json"), "w"), indent=1)
    print(summ.round(4).to_string(index=False))
    return pd.concat([rho_ci, c27])


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Figure 4 from the merged-author arm (merged_figure)")
    ap.add_argument("--merged_dir", required=True)
    ap.add_argument("--current_rho_ci", required=True, help="27-type fig4b_rho_ci.tsv")
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--qa_dir", default=None)
    ap.add_argument("--stat", choices=sorted(MERGED_STATS), default="log2FoldChange",
                    help="DE statistic correlated between author and pipeline DE (t = limma moderated t)")
    ap.add_argument("--pooling_dir", default=None, help="subtype_vs_pooled/ folder: adds panels D and E")
    a = ap.parse_args()
    print(merged_figure(a.merged_dir, a.current_rho_ci, a.out_png, a.qa_dir, a.stat, a.pooling_dir).round(4).to_string(index=False))

