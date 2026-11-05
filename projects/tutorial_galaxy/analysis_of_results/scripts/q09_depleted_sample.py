"""Q9. Which sample was depleted for the top serum proteins?

Clean set. Compares the share of each high-abundance serum protein (the
targets of a top-14 depletion, see mq_common.TOP_SERUM_GENES) between the
samples. The depleted sample is the one with the lower albumin share,
albumin being by far the most abundant serum protein.
"""
import numpy as np
from matplotlib.ticker import FuncFormatter

import mq_common as mq

QUESTION = "Which sample was depleted for the top serum proteins?"
TUTORIAL = "Sample1 (Sample2 is pure serum)"
NOTE = ""


def main():
    pg = mq.load_protein_groups()
    cl = mq.clean(pg)
    s = mq.samples(cl)
    colors = mq.sample_colors(s)
    dep, und = mq.depleted_sample(pg)
    top, alb = mq.is_top_serum(cl), mq.is_albumin(cl)
    sh = {x: mq.share(cl, x) for x in s}
    alb_sh = {x: sh[x][alb].sum() for x in s}
    top_sh = {x: sh[x][top].sum() for x in s}
    idx = sh[und][top].sort_values().index  # largest at the top of the plot
    lower = int((sh[dep][idx] < sh[und][idx]).sum())

    gene = mq.genes(cl).str[0]
    names = [f"{g}  {n}" for g, n in zip(gene[idx], mq.labels(cl, maxlen=26)[idx])]
    fig, ax = mq.plt.subplots(figsize=(7.5, 0.32 * len(idx) + 1.2))
    y = np.arange(len(idx))
    lo = np.minimum(sh[dep][idx], sh[und][idx]).replace(0, np.nan)
    hi = np.maximum(sh[dep][idx], sh[und][idx])
    ax.hlines(y, lo, hi, color=mq.GRID, linewidth=2, zorder=1)
    for x in s:
        v = sh[x][idx].replace(0, np.nan)
        ax.scatter(v, y, s=64, color=colors[x], edgecolor=mq.SURFACE, linewidth=1.5, zorder=3,
                   label=f"{x}{' (depleted)' if x == dep else ''}")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
    ax.set_yticks(y, names)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("% of sample intensity (log scale)")
    ax.set_title(f"Top serum proteins: share of intensity, {dep} vs {und}")
    ax.legend(loc="lower right")
    mq.value_grid(ax, "x")
    figure = mq.save_fig(fig, "q09_depleted_sample")

    rows = [(x, mq.pct(alb_sh[x]), mq.pct(top_sh[x])) for x in s]
    text = (
        f"**{dep} was depleted.** Albumin makes up {mq.pct(alb_sh[dep])} of its intensity, "
        f"against {mq.pct(alb_sh[und])} in {und}.\n\n"
        + mq.md_table(rows, ["Sample", "Albumin share", f"All {int(top.sum())} top serum protein groups"])
        + f"\n\n{lower} of the {int(top.sum())} top serum protein groups have a lower share in {dep}. "
        "Shares are relative: removing albumin alone raises the share of every other protein, "
        "so Q10 compares the other targets on the non-albumin intensity."
    )
    return mq.answer(QUESTION, text, figure, short=f"{dep} ({und} is not depleted)")


if __name__ == "__main__":
    print(main())
