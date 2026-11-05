"""Q6. Which sample yielded more protein identifications?

Clean set; a protein group counts as identified in a sample if it has at
least one identified peptide there ("Peptides <sample>" > 0). Match between
runs was off, so every identification comes from that sample's own spectra.
summary.txt adds the MS/MS identification rate per raw file as context.
"""
import mq_common as mq

QUESTION = "Which sample yielded more protein identifications?"
TUTORIAL = "Sample 1"
NOTE = ""


def main():
    cl = mq.clean(mq.load_protein_groups())
    s = mq.samples(cl)
    colors = mq.sample_colors(s)
    ids = {x: int(mq.identified_in(cl, x).sum()) for x in s}

    summ = mq.load_summary()
    raw, exp = mq.col(summ, "Raw file"), mq.col(summ, "Experiment")
    id_pct = mq.col(summ, "MS/MS identified [%]")
    seqs = mq.col(summ, "Peptide sequences identified")

    fig, (a1, a2) = mq.plt.subplots(1, 2, figsize=(8, 3.4))
    b1 = a1.bar(s, [ids[x] for x in s], color=[colors[x] for x in s], width=0.6, **mq.bar_kw())
    mq.label_bars(a1, b1, [mq.n(ids[x]) for x in s])
    a1.set_ylabel("Protein groups identified")
    a1.set_title("Protein identifications")
    b2 = a2.bar(summ[raw], summ[id_pct], color=[colors.get(e, mq.NEUTRAL) for e in summ[exp]],
                width=0.6, **mq.bar_kw())
    mq.label_bars(a2, b2, [mq.pct(v) for v in summ[id_pct]])
    a2.set_ylabel("% of MS/MS spectra identified")
    a2.set_title("MS/MS identification rate")
    for ax in (a1, a2):
        ax.tick_params(axis="x", length=0)
        mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q06_protein_ids_per_sample")

    hi = max(s, key=ids.get)
    lo = min(s, key=ids.get)
    rows = [(r[raw], r[exp], mq.n(ids.get(r[exp], 0)), mq.n(r[seqs]), mq.pct(r[id_pct]))
            for _, r in summ.iterrows()]
    text = (
        f"**{hi}** yielded more protein identifications: {mq.n(ids[hi])} vs {mq.n(ids[lo])} "
        f"protein groups ({mq.n(ids[hi] - ids[lo])} more).\n\n"
        + mq.md_table(rows, ["Raw file", "Sample", "Protein groups", "Peptide sequences (summary.txt)", "MS/MS identified"])
    )
    return mq.answer(QUESTION, text, figure, short=f"{hi} ({mq.n(ids[hi])} vs {mq.n(ids[lo])} protein groups)")


if __name__ == "__main__":
    print(main())
