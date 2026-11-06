"""Q1. How many proteins and features were identified in total?

Proteins = rows of proteinGroups.txt; features = rows of evidence.txt (one peptide
feature per run), as in the tutorial. The figure splits each run's evidence rows into
features identified by MS/MS and features transferred by match between runs.
"""
import mq_common as mq

QUESTION = "How many proteins and features were identified in total?"
TUTORIAL = "2512 protein groups and ~240,000 features"
NOTE = "features = rows of evidence.txt, as in the tutorial"


def main():
    pg, ev = mq.load_protein_groups(), mq.load_evidence()
    run = ev[mq.col(ev, "Raw file")].map(mq.core_name)
    matched = ev[mq.col(ev, "Type")].fillna("").str.contains("MATCH", regex=False)
    order = mq.samples()
    msms = run[~matched].value_counts().reindex(order, fill_value=0)
    mbr = run[matched].value_counts().reindex(order, fill_value=0)
    colors = mq.condition_colors()

    fig, ax = mq.plt.subplots(figsize=(10, 3.8))
    c = [colors[mq.condition_of(s)] for s in order]
    ax.bar(order, msms, color=c, width=0.75, **mq.bar_kw())
    ax.bar(order, mbr, bottom=msms, color=c, alpha=0.4, width=0.75, **mq.bar_kw())
    ax.set_ylabel("Evidence rows (features)")
    ax.set_title(f"{mq.n(len(pg))} protein groups, {mq.n(len(ev))} features")
    mq.rotate_xticks(ax)
    handles = [mq.plt.Rectangle((0, 0), 1, 1, color=colors[k]) for k in colors]
    handles += [mq.plt.Rectangle((0, 0), 1, 1, color=mq.MUTED),
                mq.plt.Rectangle((0, 0), 1, 1, color=mq.MUTED, alpha=0.4)]
    ax.legend(handles, list(colors) + ["identified by MS/MS", "match between runs"],
              loc="upper left", bbox_to_anchor=(1.0, 1.0))
    mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q01_proteins_features_total")

    text = (
        f"**{mq.n(len(pg))} protein groups** and **{mq.n(len(ev))} features** "
        f"(rows of proteinGroups.txt and evidence.txt). {mq.n(matched.sum())} features "
        f"({matched.mean() * 100:.0f}%) were transferred between runs by match between runs; "
        f"per run between {mq.n((msms + mbr).min())} and {mq.n((msms + mbr).max())} features."
    )
    return mq.answer(QUESTION, text, figure,
                     short=f"{mq.n(len(pg))} protein groups, {mq.n(len(ev))} features")


if __name__ == "__main__":
    print(main())
