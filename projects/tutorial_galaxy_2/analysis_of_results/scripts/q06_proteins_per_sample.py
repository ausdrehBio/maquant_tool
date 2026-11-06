"""Q6. Which sample has the lowest amount of proteins after protein summarization?

Proteins with a summarised intensity per run in MSstats' ProteinLevelData (the
tutorial counts the rows per run with Datamash).
"""
import mq_common as mq

QUESTION = "Which sample has the lowest amount of proteins after protein summarization?"
TUTORIAL = "RDEB cSCC4"
NOTE = ""


def main():
    prot = mq.load_msstats("ProteinLevelData")
    prot = prot[prot[mq.col(prot, "LogIntensities")].notna()]
    per_run = prot[mq.col(prot, "originalRUN")].map(mq.sample_of_run).value_counts()
    order = mq.samples()
    counts = per_run.reindex(order, fill_value=0)
    low = counts.idxmin()
    colors = mq.condition_colors()

    fig, ax = mq.plt.subplots(figsize=(10, 3.8))
    bars = ax.bar(order, counts, color=[colors[mq.condition_of(s)] for s in order], width=0.75,
                  **mq.bar_kw())
    mq.label_bars(ax, bars, [mq.n(v) if s == low else "" for s, v in counts.items()])
    ax.set_ylabel("Proteins after summarisation")
    ax.set_title(f"Fewest proteins: {low} ({mq.n(counts[low])})")
    mq.rotate_xticks(ax)
    handles = [mq.plt.Rectangle((0, 0), 1, 1, color=colors[k]) for k in colors]
    ax.legend(handles, list(colors), loc="upper left", bbox_to_anchor=(1.0, 1.0))
    mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q06_proteins_per_sample")

    text = (
        f"**{low}** has the fewest proteins after summarisation: {mq.n(counts[low])}, against a "
        f"median of {mq.n(counts.median())} and a maximum of {mq.n(counts.max())} "
        f"({counts.idxmax()})."
    )
    return mq.answer(QUESTION, text, figure, short=f"{low} ({mq.n(counts[low])} proteins)")


if __name__ == "__main__":
    print(main())
