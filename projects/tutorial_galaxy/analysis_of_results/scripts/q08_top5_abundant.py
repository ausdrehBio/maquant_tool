"""Q8. What are the top 5 most abundant proteins in both files?

Clean set; abundance = "Intensity <sample>", shown as % of the sample's
summed intensity. Top 5 per sample.
"""
import mq_common as mq

QUESTION = "What are the top 5 most abundant proteins in both files?"
TUTORIAL = ("Sample1: Complement C4-A, Ceruloplasmin, Hemopexin, Serum albumin, Complement factor B. "
            "Sample2: Serum albumin, Immunoglobulin heavy constant gamma 1, Serotransferrin, "
            "Immunoglobulin kappa constant, Haptoglobin")
NOTE = "tutorial ranks by LFQ intensity, here by Intensity (see below)"
TOP_N = 5


def main():
    cl = mq.clean(mq.load_protein_groups())
    s = mq.samples(cl)
    colors = mq.sample_colors(s)
    lab = mq.labels(cl, maxlen=40)
    full = mq.labels(cl, maxlen=80)
    gene = mq.genes(cl).str[0]

    fig, axes = mq.plt.subplots(1, len(s), figsize=(10, 3.2), sharex=True)
    parts, shorts = [], []
    for ax, x in zip(axes, s):
        sh = mq.share(cl, x)
        top = sh.nlargest(TOP_N).index
        bars = ax.barh(lab[top][::-1], sh[top][::-1], color=colors[x], height=0.6, **mq.bar_kw())
        mq.label_bars(ax, bars, [mq.pct(v) for v in sh[top][::-1]])
        ax.set_title(x)
        ax.set_xlabel("% of sample intensity")
        ax.tick_params(axis="y", length=0)
        mq.value_grid(ax, "x")
        inten = mq.intensity(cl, x)
        rows = [(r, full[i], gene[i], mq.sci(inten[i]), mq.pct(sh[i])) for r, i in enumerate(top, 1)]
        parts.append(f"**{x}**\n\n" + mq.md_table(rows, ["Rank", "Protein", "Gene", "Intensity", "Share"]))
        shorts.append(f"{x}: " + ", ".join(full[top]))
    axes[-1].set_xlim(right=axes[-1].get_xlim()[1] * 1.12)  # room for the value labels
    fig.suptitle(f"Top {TOP_N} proteins by share of sample intensity", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    figure = mq.save_fig(fig, "q08_top5_abundant")

    text = "\n\n".join(parts)
    return mq.answer(QUESTION, text, figure, short=". ".join(shorts))


if __name__ == "__main__":
    print(main())
