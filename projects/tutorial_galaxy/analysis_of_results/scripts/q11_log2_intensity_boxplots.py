"""Q11. Show the distribution of the log2 intensities via box plots.

Clean set; log2 of "Intensity <sample>" and of "LFQ intensity <sample>",
zeros excluded (log2 of 0 is undefined; a 0 means not quantified).
"""
import numpy as np
from matplotlib.colors import to_rgba

import mq_common as mq

QUESTION = "Show the distribution of the log2 intensities via box plots."
TUTORIAL = "Close to normal; depleted sample shifted to higher values"
NOTE = "tutorial plots log2 LFQ (MaxQuant 1.6.10), here log2 Intensity and log2 LFQ"


def main():
    pg = mq.load_protein_groups()
    cl = mq.clean(pg)
    dep = mq.depleted_sample(pg)[0]
    s = mq.samples(cl)
    colors = mq.sample_colors(s)
    measures = [("Intensity", mq.intensity), ("LFQ intensity", mq.lfq)]
    data = {m: [np.log2(f(cl, x)[f(cl, x) > 0]) for x in s] for m, f in measures}

    fig, axes = mq.plt.subplots(1, 2, figsize=(8, 4), sharey=True)
    for ax, (m, _) in zip(axes, measures):
        bp = ax.boxplot(data[m], widths=0.5, patch_artist=True,
                        medianprops={"color": mq.INK, "linewidth": 1.5},
                        whiskerprops={"color": mq.MUTED}, capprops={"color": mq.MUTED},
                        flierprops={"marker": "o", "markersize": 3, "markerfacecolor": mq.MUTED,
                                    "markeredgecolor": "none", "alpha": 0.7})
        for patch, x in zip(bp["boxes"], s):
            patch.set_facecolor(to_rgba(colors[x], 0.25))
            patch.set_edgecolor(colors[x])
            patch.set_linewidth(1.5)
        ax.set_xticks(range(1, len(s) + 1), [f"{x}\nn = {len(d)}" for x, d in zip(s, data[m])])
        ax.set_title(f"log2 {m}")
        ax.tick_params(axis="x", length=0)
        mq.value_grid(ax, "y")
    axes[0].set_ylabel("log2 intensity")
    fig.tight_layout()
    figure = mq.save_fig(fig, "q11_log2_intensity_boxplots")

    rows = []
    for m, _ in measures:
        for x, d in zip(s, data[m]):
            q1, med, q3 = np.percentile(d, [25, 50, 75])
            rows.append((m, x, mq.n(len(d)), f"{med:.2f}", f"{q1:.2f} to {q3:.2f}"))
    med = {m: {x: float(np.median(d)) for x, d in zip(s, data[m])} for m, _ in measures}
    diff_i = med["Intensity"][s[0]] - med["Intensity"][s[1]]
    diff_l = med["LFQ intensity"][s[0]] - med["LFQ intensity"][s[1]]
    text = (
        mq.md_table(rows, ["Measure", "Sample", "Proteins", "Median log2", "Interquartile range"])
        + f"\n\nThe median log2 Intensity differs by {diff_i:+.2f} between {s[0]} and {s[1]}; "
        f"after LFQ normalisation by {diff_l:+.2f}. LFQ intensities only exist for proteins "
        "quantified in both samples, so n is smaller."
    )
    higher = max(s, key=med["Intensity"].get)
    who = f"depleted sample ({dep})" if higher == dep else higher
    short = ("Median log2 Intensity " + ", ".join(f"{x} {med['Intensity'][x]:.2f}" for x in s)
             + f"; {who} higher")
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
