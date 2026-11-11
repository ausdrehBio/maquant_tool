"""Summary figures for the top of answers.md (called by run_all.py).

summary_verdicts.png
    One row per check: the verdict, the key metric of the check on its own scale with the zones
    of its verdict rule (green consistent, amber minor deviation, red inconsistent), our value
    (blue dot) and, where the rule compares with it, the authors' own deposited search (orange
    diamond).
summary_ours_vs_authors.png
    Our result divided by the authors' deposited search for the metrics that both reanalyses
    share; 1 means our workflow does exactly as well as the authors' own MaxQuant output.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

import mq_common as mq

VERDICT_COLORS = {"consistent": "#1baf7a", "minor deviation": "#eda100",
                  "inconsistent": "#c23b3a", "not judged": "#c3c2b7"}
OURS = mq.SERIES[0]
AUTHORS = mq.SERIES[1]
TITLES = ["Search settings", "Sample assignment", "Identifications", "Same-file quantities",
          "Statistics code", "Differential results", "Sample structure", "Paper conclusions"]
PART_LABELS = {"4b": "4b viral", "4c": "4c counts", "4e": "4e KS", "4e_strength": "4e size",
               "4e_tail": "4e tails", "4f_induction": "4f wt up", "4f_similar": "4f RSIR≈wt",
               "4f_RFMG_lt_RSIG": "4f RFMG<wt", "4f_AFMG_lt_RSIG": "4f AFMG<wt",
               "4f_RFMG_lt_RFMR": "4f RFMR>RFMG", "4f_AFMG_lt_AFMR": "4f AFMR>AFMG"}


def _specs(answers):
    m = [a.metrics for a in answers]
    q6 = m[5]["comparisons"]
    p_dep = m[6]["procrustes_deposited"]
    proc_cut = 1 - 2 * (1 - p_dep)
    return [
        {"metric": "settings that affect the result and differ", "xlim": (-0.5, 5),
         "zones": [(-0.5, 0.5, "consistent"), (0.5, 5, "minor deviation")],
         "ours": m[0]["relevant_differences"], "fmt": "{:.0f}"},
        {"metric": "smallest profile margin to another condition (Δr)", "xlim": (-0.3, 0.8),
         "zones": [(-0.3, 0, "inconsistent"), (0, 0.8, "consistent")],
         "ours": m[1]["margin_ours"], "authors": m[1]["margin_deposited"], "fmt": "{:.2f}"},
        {"metric": "groups identified by MS/MS per raw file, ours ÷ authors'",
         "xlim": (0.7, 1.35),
         "zones": [(0.7, 0.8, "inconsistent"), (0.8, 0.9, "minor deviation"),
                   (0.9, 1.1, "consistent"), (1.1, 1.25, "minor deviation"),
                   (1.25, 1.35, "inconsistent")],
         "ours": m[2]["msms_ratio"], "authors": 1.0, "fmt": "{:.3f}"},
        {"metric": "same-file difference ÷ replicate difference (LFQ)", "xlim": (0, 1.4),
         "zones": [(0, 0.7, "consistent"), (0.7, 1.0, "minor deviation"), (1.0, 1.4, "inconsistent")],
         "ours": m[3]["ratio"], "ref": (1.0, "replicate noise"), "fmt": "{:.2f}"},
        {"metric": "largest deviation from the published p values (log10 p)", "xlim": (1e-15, 1e-1),
         "log": True, "zones": [(1e-15, 1e-10, "consistent"), (1e-10, 1e-1, "inconsistent")],
         "ours": m[4]["max_dlog10p_t"], "fmt": "{:.0e}"},
        {"metric": "fold-change r vs published, ours − authors' (9 comparisons)",
         "xlim": (-0.15, 0.05),
         "zones": [(-0.15, -0.05, "minor deviation"), (-0.05, 0.05, "consistent")],
         "dots": [(q6[c]["r_ours"] - q6[c]["r_deposited"], c in mq.STRAINS)
                  for c in mq.COMPARISON_NAMES],
         "authors": 0.0, "fmt": "{:+.3f}"},
        {"metric": "Procrustes similarity of the 42 samples to published", "xlim": (0.95, 1.0),
         "zones": [(0.95, proc_cut, "minor deviation"), (proc_cut, 1.0, "consistent")],
         "ours": m[6]["procrustes"], "authors": p_dep, "fmt": "{:.3f}"},
        {"metric": "claim parts of Fig. 4 (grey: rests on the authors' single random draw)",
         "parts": m[7]["parts"]},
    ]


def _badge(ax, verdict):
    color = VERDICT_COLORS[verdict]
    ax.add_patch(FancyBboxPatch((0.02, 0.12), 0.62, 0.36, boxstyle="round,pad=0.02,rounding_size=0.08",
                                transform=ax.transAxes, facecolor=color, edgecolor="none"))
    ax.text(0.33, 0.30, verdict, transform=ax.transAxes, ha="center", va="center",
            color="white", fontsize=9, fontweight="bold")


def _bullet(ax, spec):
    lo, hi = spec["xlim"]
    for a, b, v in spec["zones"]:
        ax.axvspan(a, b, color=VERDICT_COLORS[v], alpha=0.22, lw=0)
    if spec.get("log"):
        ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", labelsize=7, length=2)
    if "ref" in spec:
        x, label = spec["ref"]
        ax.axvline(x, color=mq.INK2, lw=0.8, ls="--")
        ax.text(x, 0.62, f" {label}", fontsize=7, color=mq.INK2, va="center")
    if spec.get("authors") is not None:
        ax.plot(spec["authors"], 0, marker="D", ms=9, mfc=mq.SURFACE, mec=AUTHORS, mew=1.8,
                zorder=3)
    if "dots" in spec:
        for x, filled in spec["dots"]:
            ax.plot(np.clip(x, lo, hi), 0, "o", ms=6.5, mfc=OURS if filled else mq.SURFACE,
                    mec=OURS, mew=1.3, alpha=0.85, zorder=4)
        vals = [x for x, _ in spec["dots"]]
        text = f"{min(vals):+.3f}\n{max(vals):+.3f}"
    else:
        x = spec["ours"]
        ax.plot(np.clip(x, lo, hi), 0, "o", ms=10, color=OURS, zorder=4)
        text = spec["fmt"].format(x)
    ax.set_title(spec["metric"], fontsize=8, color=mq.INK2, fontweight="normal", loc="left",
                 pad=3)
    ax.text(1.01, 0, text, transform=ax.get_yaxis_transform(), fontsize=9, color=OURS,
            va="center", fontweight="bold")


def _parts(ax, spec):
    parts = spec["parts"]
    order = [p for p in PART_LABELS if p in parts]
    for i, p in enumerate(order):
        v = parts[p]["verdict"]
        ax.add_patch(FancyBboxPatch((i + 0.08, -0.45), 0.84, 0.9,
                                    boxstyle="round,pad=0.0,rounding_size=0.12",
                                    facecolor=VERDICT_COLORS[v], edgecolor="none"))
        ax.text(i + 0.5, -0.82, PART_LABELS[p], ha="right", va="top", fontsize=6.5,
                color=mq.INK2, rotation=35, rotation_mode="anchor")
    ax.set_xlim(0, len(order))
    ax.set_ylim(-1, 1)
    ax.axis("off")
    ax.set_title(spec["metric"], fontsize=8, color=mq.INK2, fontweight="normal", loc="left",
                 pad=3)


def verdict_overview(results):
    answers = [a for _, a in results]
    specs = _specs(answers)
    counts = {v: sum(a.verdict == v for a in answers) for v in mq.VERDICTS}
    n = len(answers)
    fig = plt.figure(figsize=(13.5, 1.05 * n + 2.2))
    gs = fig.add_gridspec(n + 2, 3, width_ratios=[2.3, 5.2, 3.6],
                          height_ratios=[0.55] + [1] * n + [0.55], hspace=0.95, wspace=0.08,
                          left=0.015, right=0.985, top=0.945, bottom=0.025)
    # header: counts as a stacked bar
    hax = fig.add_subplot(gs[0, :])
    hax.axis("off")
    left = 0.0
    for v in mq.VERDICTS:
        if counts[v]:
            hax.barh(0, counts[v], left=left, color=VERDICT_COLORS[v], height=0.9, **mq.bar_kw())
            hax.text(left + counts[v] / 2, 0, f"{counts[v]} {v}", ha="center", va="center",
                     color="white", fontsize=10, fontweight="bold")
            left += counts[v]
    hax.set_xlim(0, n)
    hax.set_ylim(-0.5, 0.5)
    title = (f"{counts['consistent']} of {n} checks consistent  ·  {counts['minor deviation']} "
             f"minor deviation  ·  {counts['inconsistent']} inconsistent")
    fig.suptitle(title, x=0.015, y=0.993, ha="left", fontsize=16, fontweight="bold", color=mq.INK)
    for i, (spec, ans) in enumerate(zip(specs, answers)):
        tax = fig.add_subplot(gs[i + 1, 0])
        tax.axis("off")
        tax.text(0.02, 0.82, f"Q{i + 1}  {TITLES[i]}", transform=tax.transAxes, fontsize=10,
                 fontweight="bold", color=mq.INK, va="center")
        _badge(tax, ans.verdict)
        cax = fig.add_subplot(gs[i + 1, 1])
        if "parts" in spec:
            _parts(cax, spec)
        else:
            _bullet(cax, spec)
        rax = fig.add_subplot(gs[i + 1, 2])
        rax.axis("off")
        rax.text(0.17, 0.5, _wrap(ans.short, 48), transform=rax.transAxes, fontsize=8,
                 color=mq.INK2, va="center")
    lax = fig.add_subplot(gs[n + 1, :])
    lax.axis("off")
    handles = [Line2D([], [], marker="o", ls="", color=OURS, ms=8, label="ours"),
               Line2D([], [], marker="D", ls="", mfc=mq.SURFACE, mec=AUTHORS, mew=1.8, ms=8,
                      label="authors' deposited search (baseline)")]
    handles += [Patch(color=VERDICT_COLORS[v], alpha=0.5 if v != "not judged" else 1, label=v)
                for v in ["consistent", "minor deviation", "inconsistent", "not judged"]]
    lax.legend(handles=handles, ncol=6, loc="center", fontsize=8.5, frameon=False)
    return mq.save_fig(fig, "summary_verdicts")


def _wrap(text, width):
    words, lines, line = text.split(), [], ""
    for w in words:
        if len(line) + len(w) + 1 > width:
            lines.append(line)
            line = w
        else:
            line = f"{line} {w}".strip()
    lines.append(line)
    return "\n".join(lines)


def ours_vs_authors(results):
    m = [a.metrics for _, a in results]
    q6 = m[5]["comparisons"]
    judged = [p for p, d in m[7]["parts"].items() if p != "4b" and d["verdict"] != "not judged"]
    rows = [
        ("Q2", "profile correlation with published, same sample", m[1]["r_same_ours"],
         m[1]["r_same_deposited"], "{:.3f}"),
        ("Q3", "groups identified by MS/MS per raw file", m[2]["msms_ours"],
         m[2]["msms_deposited"], "{:,.0f}"),
        ("Q3", "MS/MS spectra assigned per raw file", m[2]["msms_count_ours"],
         m[2]["msms_count_deposited"], "{:,.0f}"),
        ("Q3", "published protein groups matched", m[2]["coverage_ours"],
         m[2]["coverage_deposited"], "{:.1%}"),
        ("Q4", "deviation from published values, residual SD (lower = closer)",
         m[3]["pub_residual_ours"], m[3]["pub_residual_deposited"], "{:.3f}"),
    ]
    for c in mq.COMPARISON_NAMES:
        rows.append(("Q6", f"fold-change r vs published, {mq.comparison_label(c)}",
                     q6[c]["r_ours"], q6[c]["r_deposited"], "{:.3f}"))
    rows += [
        ("Q6", "direction agreement with published hits (median)",
         float(np.median([q6[c]["direction_ours"] for c in q6])),
         float(np.median([q6[c]["direction_deposited"] for c in q6])), "{:.1%}"),
        ("Q6", "significant proteins, all 9 comparisons",
         sum(q6[c]["up_ours"] + q6[c]["down_ours"] for c in q6),
         sum(q6[c]["up_deposited"] + q6[c]["down_deposited"] for c in q6), "{:,.0f}"),
        ("Q6", "recall of the published calls (median)",
         float(np.median([q6[c]["recall_ours"] for c in q6])),
         float(np.median([q6[c]["recall_deposited"] for c in q6])), "{:.1%}"),
        ("Q7", "Procrustes similarity to published", m[6]["procrustes"],
         m[6]["procrustes_deposited"], "{:.3f}"),
        ("Q7", "centroid-distance correlation with published", m[6]["centroid_r"],
         m[6]["centroid_r_deposited"], "{:.3f}"),
        ("Q8", "viral protein pattern, r with published", m[7]["viral_r"],
         m[7]["viral_r_deposited"], "{:.3f}"),
        ("Q8", "judged claim parts, share of seeds in which they hold (mean)",
         float(np.mean([m[7]["parts"][p]["ours"] for p in judged])),
         float(np.mean([m[7]["parts"][p]["deposited"] for p in judged])), "{:.0%}"),
    ]
    ratios = np.array([o / a for _, _, o, a, _ in rows])
    y = np.arange(len(rows))[::-1]
    fig, ax = plt.subplots(figsize=(12.5, 0.36 * len(rows) + 1.6))
    ax.axvspan(0.95, 1.05, color=VERDICT_COLORS["consistent"], alpha=0.15, lw=0)
    ax.axvline(1.0, color=mq.INK2, lw=1)
    inside = (ratios >= 0.95) & (ratios <= 1.05)
    ax.scatter(ratios, y, s=46, color=np.where(inside, OURS, VERDICT_COLORS["minor deviation"]),
               zorder=3)
    for yi, (q, label, o, a, fmt) in zip(y, rows):
        ax.text(0.003, yi, f"{q}  {label}", transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=8, color=mq.INK)
        ax.text(1.01, yi, f"ours {fmt.format(o)}  ·  authors' {fmt.format(a)}",
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.5,
                color=mq.INK2)
    prev = None
    for yi, (q, *_rest) in zip(y, rows):
        if prev is not None and q != prev:
            ax.axhline(yi + 0.5, color=mq.GRID, lw=0.8)
        prev = q
    ax.set_xscale("log")
    lo = min(0.9, ratios.min() * 0.98)
    hi = max(1.1, ratios.max() * 1.02)
    ax.set_xlim(lo, hi)
    ticks = [t for t in (0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2) if lo <= t <= hi]
    ax.set_xticks(ticks, [f"{t:.2f}" for t in ticks])
    ax.minorticks_off()
    ax.set_yticks([])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("ours ÷ authors' deposited search (1 = the same result; green band: within 5%)")
    n_in = int(inside.sum())
    ax.set_title(f"Ours next to the authors' own reanalysis: {n_in} of {len(rows)} metrics within "
                 "5% of the authors' result", loc="left")
    fig.subplots_adjust(left=0.40, right=0.80)
    return mq.save_fig(fig, "summary_ours_vs_authors")


def make(results):
    return [verdict_overview(results), ours_vs_authors(results)]
