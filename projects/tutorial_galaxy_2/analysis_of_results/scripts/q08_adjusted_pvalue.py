"""Q8. Why do we filter for the adjusted p-value?

MSstats tests every protein separately. At a raw p-value cutoff of 0.05, about 5% of
the proteins that do not change at all would still pass. The Benjamini-Hochberg
adjusted p-value controls the false discovery rate over all tested proteins instead.
The figure shows the raw p-value distribution against what pure chance would give.
"""
import numpy as np

import mq_common as mq

QUESTION = "Why do we filter for the adjusted p-value?"
TUTORIAL = ("Adjusted p-values account for testing many proteins: a 0.05 cutoff controls the "
            "false discovery rate over all proteins at 5%")
NOTE = ""


def main():
    cr = mq.load_msstats("ComparisonResult")
    p = cr[mq.col(cr, "pvalue")]
    adj = cr[mq.col(cr, "adj.pvalue")]
    tested = p.notna()
    n_tests = int(tested.sum())
    n_raw = int((p < mq.P_CUTOFF).sum())
    n_adj = int((adj < mq.P_CUTOFF).sum())
    chance = mq.P_CUTOFF * n_tests

    fig, ax = mq.plt.subplots(figsize=(7, 3.6))
    bins = np.arange(0, 1.0001, 0.05)
    ax.hist(p[tested], bins=bins, color=mq.NEUTRAL, edgecolor=mq.SURFACE, linewidth=0.8)
    ax.axhline(n_tests / (len(bins) - 1), color=mq.SERIES[1], linewidth=1.5)
    ax.text(1.0, n_tests / (len(bins) - 1), "  expected if no protein changed",
            color=mq.INK2, fontsize=9, va="bottom", ha="right")
    ax.set_xlabel("Raw p-value (bins of 0.05)")
    ax.set_ylabel("Proteins")
    ax.set_title(f"{mq.n(n_raw)} proteins with p < 0.05, {mq.n(n_adj)} with adjusted p < 0.05")
    mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q08_adjusted_pvalue")

    text = (
        f"Because MSstats tests {mq.n(n_tests)} proteins separately. With the raw p-value, about "
        f"{mq.n(chance)} proteins (5% of {mq.n(n_tests)}) would pass p < 0.05 by chance even if no "
        f"protein changed; {mq.n(n_raw)} pass here. The Benjamini-Hochberg adjusted p-value "
        f"controls the false discovery rate instead: among the {mq.n(n_adj)} proteins with "
        "adjusted p < 0.05, about 5% are expected to be false positives. The adjustment is more "
        "conservative, so fewer proteins pass."
    )
    short = f"Multiple testing: {mq.n(n_raw)} with raw p < 0.05 vs {mq.n(n_adj)} adjusted"
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
