"""Q12. How much is the log2 fold change of serum albumin?

log2(depleted / undepleted) of albumin, from "Intensity" (the measure used
for the abundance questions) and from "LFQ intensity". The median log2 fold
change of all proteins quantified in both samples shows how much of the
difference is a global shift between the samples, which LFQ normalises away.
"""
import numpy as np

import mq_common as mq

QUESTION = "How much is the log2 fold change of serum albumin?"
TUTORIAL = "About -2.91 (from LFQ intensity)"
NOTE = "tutorial: LFQ of MaxQuant 1.6.10; here Intensity (see below)"


def main():
    pg = mq.load_protein_groups()
    cl = mq.clean(pg)
    dep, und = mq.depleted_sample(pg)
    alb = mq.is_albumin(cl)
    i_dep, i_und = mq.intensity(cl, dep), mq.intensity(cl, und)
    l_dep, l_und = mq.lfq(cl, dep), mq.lfq(cl, und)
    fc_int = float(np.log2(i_dep[alb].iloc[0] / i_und[alb].iloc[0]))
    fc_lfq = float(np.log2(l_dep[alb].iloc[0] / l_und[alb].iloc[0]))
    both = (i_dep > 0) & (i_und > 0)
    fc_all = np.log2(i_dep[both] / i_und[both])
    med_all = float(np.median(fc_all))

    fig, ax = mq.plt.subplots(figsize=(5.6, 5.2))
    x, y = np.log2(i_und[both]), np.log2(i_dep[both])
    rest = both & ~alb
    ax.scatter(np.log2(i_und[rest]), np.log2(i_dep[rest]), s=40, color=mq.MUTED, alpha=0.55,
               edgecolor=mq.SURFACE, linewidth=0.8, label="Other proteins")
    ax.scatter(np.log2(i_und[alb]), np.log2(i_dep[alb]), s=80, color=mq.INK, edgecolor=mq.SURFACE,
               linewidth=1.5, zorder=3, label="Albumin")
    lims = [min(x.min(), y.min()) - 0.5, max(x.max(), y.max()) + 0.5]
    ax.plot(lims, lims, color=mq.AXIS, linewidth=1, zorder=0)
    ax.text(lims[1], lims[1], "no change ", ha="right", va="bottom", fontsize=8, color=mq.INK2)
    ax.annotate(f"Albumin\nlog2FC = {fc_int:.2f}".replace("-", "−"), (np.log2(i_und[alb].iloc[0]), np.log2(i_dep[alb].iloc[0])),
                xytext=(-10, -6), textcoords="offset points", ha="right", va="top", fontsize=9, color=mq.INK)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")
    ax.set_xlabel(f"log2 Intensity {und}")
    ax.set_ylabel(f"log2 Intensity {dep} (depleted)")
    ax.set_title("Proteins quantified in both samples")
    ax.legend(loc="upper left")
    mq.value_grid(ax, "both")
    figure = mq.save_fig(fig, "q12_albumin_log2fc")

    direction = "higher" if med_all > 0 else "lower"
    text = (
        f"The log2 fold change of albumin ({dep} depleted vs {und}) is **{fc_int:.2f}** from "
        f"Intensity, i.e. {2 ** -fc_int:.1f}-fold lower in the depleted sample.\n\n"
        f"From LFQ intensity it is {fc_lfq:.2f}. The difference comes from normalisation: across all "
        f"{mq.n(both.sum())} proteins quantified in both samples, the median log2 fold change of "
        f"Intensity is {med_all:+.2f}, so most proteins have {direction} intensities in the depleted "
        "sample. LFQ normalisation largely removes such a global shift; relative to the median protein, "
        f"albumin changes by {fc_int - med_all:.2f}."
    )
    return mq.answer(QUESTION, text, figure, short=f"{fc_int:.2f} (Intensity); {fc_lfq:.2f} (LFQ)")


if __name__ == "__main__":
    print(main())
