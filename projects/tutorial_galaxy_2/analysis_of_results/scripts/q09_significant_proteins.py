"""Q9. How many proteins have a p-value below 0.05?

MSstats adjusted p-value < 0.05, then split by log2 fold change beyond +-0.58 (1.5-fold)
into proteins higher in either condition, as in the tutorial's filter steps. Proteins
missing in one condition get an infinite fold change and an adjusted p-value of 0;
they count as significant here, as in the tutorial's first filter.
"""
import numpy as np

import mq_common as mq

QUESTION = "How many proteins have a p-value below 0.05?"
TUTORIAL = "198 in total; 178 up in metastasizing cSCC, 15 up in RDEB cSCC"
NOTE = "adjusted p-value; up = log2FC beyond ±0.58"


def main():
    cr = mq.load_msstats("ComparisonResult")
    name, up_cond, down_cond = mq.contrast()
    fc = cr[mq.col(cr, "log2FC")]
    adj = cr[mq.col(cr, "adj.pvalue")]
    sig = adj < mq.P_CUTOFF
    up = sig & (fc > mq.FC_CUTOFF)
    down = sig & (fc < -mq.FC_CUTOFF)
    infinite = sig & np.isinf(fc)
    colors = mq.condition_colors()

    finite = np.isfinite(fc) & adj.notna()
    y = -np.log2(adj[finite].clip(lower=1e-300))
    fig, ax = mq.plt.subplots(figsize=(7, 4.6))
    groups = [(finite & ~up & ~down, mq.AXIS, "not significant or < 1.5-fold"),
              (finite & up, colors[up_cond], f"higher in {up_cond}"),
              (finite & down, colors[down_cond], f"higher in {down_cond}")]
    for mask, color, label in groups:
        ax.scatter(fc[mask], -np.log2(adj[mask].clip(lower=1e-300)), s=18, color=color,
                   alpha=0.8, edgecolor="none", label=f"{label} ({int(mask.sum())})")
    ax.axhline(-np.log2(mq.P_CUTOFF), color=mq.INK2, linewidth=0.8)
    for x in (-mq.FC_CUTOFF, mq.FC_CUTOFF):
        ax.axvline(x, color=mq.INK2, linewidth=0.8)
    ax.set_xlabel(f"log2 fold change ({name})")
    ax.set_ylabel("-log2 adjusted p-value")
    ax.set_title(f"{int(sig.sum())} proteins with adjusted p < 0.05")
    ax.legend(loc="upper left", fontsize=8)
    mq.value_grid(ax, "both")
    figure = mq.save_fig(fig, "q09_significant_proteins")

    text = (
        f"**{int(sig.sum())} proteins** have an adjusted p-value below 0.05. Of these, "
        f"**{int(up.sum())} are higher in {up_cond}** and **{int(down.sum())} higher in "
        f"{down_cond}** (log2FC beyond ±{mq.FC_CUTOFF}, comparison {name}); "
        f"{int(sig.sum() - up.sum() - down.sum())} change less than 1.5-fold. "
        f"{int(infinite.sum())} of the significant proteins were missing in one condition "
        "(infinite fold change, not shown in the volcano plot)."
    )
    short = f"{int(sig.sum())} in total; {int(up.sum())} up in {up_cond}, {int(down.sum())} up in {down_cond}"
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
