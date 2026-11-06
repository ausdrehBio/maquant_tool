"""Q3. How large is the proportion of potential contaminants?

Share of each sample's protein group intensity that comes from groups flagged as
potential contaminant (decoys excluded), split into contaminant groups that contain a
human protein (e.g. albumin, keratins, hemoglobin) and purely non-human ones. The
tutorial reads this from the PTXQC report.
"""
import mq_common as mq

QUESTION = "How large is the proportion of potential contaminants?"
TUTORIAL = "Up to ~60% of the sample intensities (PTXQC report)"
NOTE = "contaminant share of protein group intensity per sample"


def main():
    pg = mq.load_protein_groups()
    pg = pg[~mq.is_decoy(pg)]
    con = mq.is_contaminant(pg)
    human = mq.has_human(pg[mq.col(pg, "Protein IDs")])
    cols = mq.sample_intensity_columns(pg)
    share_h, share_nh = {}, {}
    for s, c in cols.items():
        inten = pg[c].fillna(0).astype(float)
        share_h[s] = inten[con & human].sum() / inten.sum() * 100
        share_nh[s] = inten[con & ~human].sum() / inten.sum() * 100
    total = {s: share_h[s] + share_nh[s] for s in cols}
    top_s = max(total, key=total.get)
    tot_int = pg[mq.col(pg, "Intensity")].fillna(0).astype(float)
    top_con = tot_int[con].nlargest(3).index
    top_genes = [mq.gene_of(pg.loc[i, mq.col(pg, "Protein IDs")]) for i in top_con]
    top_names = ", ".join(top_genes)

    order = list(cols)
    fig, ax = mq.plt.subplots(figsize=(10, 3.8))
    ax.bar(order, [share_h[s] for s in order], color=mq.NEUTRAL, width=0.75,
           label="contaminant groups with a human protein", **mq.bar_kw())
    ax.bar(order, [share_nh[s] for s in order], bottom=[share_h[s] for s in order],
           color=mq.AXIS, width=0.75, label="non-human contaminant groups", **mq.bar_kw())
    colors = mq.condition_colors()
    for tick, s in zip(ax.get_xticklabels(), order):
        tick.set_color(colors[mq.condition_of(s)])
    ax.set_ylabel("% of sample intensity")
    ax.set_title(f"Potential contaminants: up to {mq.pct(total[top_s], 0)} of a sample's intensity")
    mq.rotate_xticks(ax)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))
    mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q03_contaminant_proportion")

    vals = sorted(total.values())
    keratin_note = (
        " These keratins are made by the squamous tumour cells themselves, so in this tissue "
        "much of the flagged intensity is endogenous rather than contamination."
        if top_genes and all(g.startswith("KRT") for g in top_genes) else ""
    )
    text = (
        f"Potential contaminants make up **up to {mq.pct(total[top_s])}** of a sample's protein "
        f"group intensity ({top_s}); the median sample has {mq.pct(vals[len(vals) // 2])}. "
        f"Almost all of it comes from contaminant groups that contain a human protein "
        f"(largest: {top_names}); purely non-human contaminants contribute at most "
        f"{mq.pct(max(share_nh.values()))}.{keratin_note} Tick labels are coloured by condition."
    )
    short = f"Up to {mq.pct(total[top_s])} ({top_s}), median {mq.pct(vals[len(vals) // 2])}"
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
