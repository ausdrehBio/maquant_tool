"""Q3. How large is the share of contaminants per run?

Share of the summed protein group Intensity of a run (decoys excluded) that comes from
groups marked "Potential contaminant", split into contaminant-only groups (CON__, for
example trypsin or human keratins) and contaminant groups with a mouse protein (for example
albumin, whose peptides are shared with bovine serum albumin).
"""
import mq_common as mq

QUESTION = "How large is the share of contaminants per run?"
PAPER = "Not reported"
NOTE = "share of summed protein group Intensity"


def main():
    pg = mq.load_protein_groups()
    pg = pg[~mq.is_decoy(pg)]
    cols = mq.experiment_columns(pg, "Intensity")
    inten = pg[list(cols.values())].astype(float)
    inten.columns = list(cols.keys())
    total = inten.sum()
    con = mq.is_contaminant(pg)
    target = mq.has_target(pg[mq.col(pg, "Protein IDs")])
    only = inten[con & ~target].sum() / total * 100
    with_mouse = inten[con & target].sum() / total * 100
    share = only + with_mouse

    top = (inten[con].div(total, axis=1) * 100).mean(axis=1).sort_values(ascending=False).head(5)
    top_names = [
        f"{mq.gene_of(pg.loc[i, mq.col(pg, 'Protein IDs')])} ({mq.pct(v)})" for i, v in top.items()
    ]

    exps = mq.experiments()
    fig, ax = mq.plt.subplots(figsize=(10, 3.8))
    x = range(len(exps))
    ax.bar(x, with_mouse[exps], color=mq.SERIES[0], label="contaminant group with a mouse protein", **mq.bar_kw())
    ax.bar(x, only[exps], bottom=with_mouse[exps], color=mq.SERIES[1], label="contaminant only (CON__)", **mq.bar_kw())
    ax.set_xticks(list(x), exps)
    mq.rotate_xticks(ax, size=7, angle=90)
    for t in ax.get_xticklabels():
        t.set_color(mq.condition_colors()[mq.condition_of(t.get_text())])
    ax.set_ylabel("% of protein group intensity")
    ax.set_title("Contaminant share per run (tick labels coloured by group)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    mq.value_grid(ax)
    figure = mq.save_fig(fig, "q03_contaminants")

    by_cond = share.groupby([mq.condition_of(e) for e in share.index]).mean().loc[mq.CONDITIONS]
    rows = [(c, mq.pct(by_cond[c])) for c in mq.CONDITIONS]
    text = (
        f"Contaminant groups make up **{mq.pct(share.min())} to {mq.pct(share.max())}** of a "
        f"run's protein group intensity (median {mq.pct(share.median())}; highest in "
        f"{share.idxmax()}). Contaminant-only groups contribute at most {mq.pct(only.max())}, "
        f"contaminant groups with a mouse protein at most {mq.pct(with_mouse.max())}. Largest "
        f"contaminant groups (mean share): {', '.join(top_names)}.\n\n"
        + mq.md_table(rows, ["Group", "Contaminant share (mean)"])
    )
    short = f"{mq.pct(share.min())} to {mq.pct(share.max())} per run (median {mq.pct(share.median())})"
    return mq.answer(QUESTION, text, figure, short)


if __name__ == "__main__":
    print(main())
