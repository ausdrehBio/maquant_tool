"""Q1. How many proteins were found in total?

Counts all protein groups in proteinGroups.txt and splits them into targets,
contaminants and decoys. MaxQuant reports protein groups: proteins that the
identified peptides cannot tell apart share one row.
"""
import mq_common as mq

QUESTION = "How many proteins were found in total?"
TUTORIAL = "271 protein groups"
NOTE = "MaxQuant 1.6.10 in the tutorial, 2.9.0 here"


def main():
    pg = mq.load_protein_groups()
    dec, con, db = mq.is_decoy(pg), mq.is_contaminant(pg), mq.has_db_protein(pg)
    cats = [
        ("Target (human FASTA)", int((~dec & ~con).sum())),
        ("Contaminant with human protein", int((con & db & ~dec).sum())),
        ("Contaminant only (CON__)", int((con & ~db & ~dec).sum())),
        ("Decoy (REV__)", int(dec.sum())),
    ]
    total = len(pg)
    assert sum(v for _, v in cats) == total
    n_clean = len(mq.clean(pg))

    fig, ax = mq.plt.subplots(figsize=(7, 2.6))
    names, vals = [c for c, _ in cats][::-1], [v for _, v in cats][::-1]
    bars = ax.barh(names, vals, color=mq.NEUTRAL, height=0.6, **mq.bar_kw())
    mq.label_bars(ax, bars, [mq.n(v) for v in vals])
    ax.set_xlabel("Protein groups")
    ax.set_title(f"{mq.n(total)} protein groups in proteinGroups.txt")
    ax.tick_params(axis="y", length=0)
    mq.value_grid(ax, "x")
    figure = mq.save_fig(fig, "q01_total_proteins")

    rows = [(c, mq.n(v)) for c, v in cats] + [("**Total**", f"**{mq.n(total)}**")]
    text = (
        f"**{mq.n(total)} protein groups** were found (all rows of proteinGroups.txt).\n\n"
        + mq.md_table(rows, ["Category", "Protein groups"])
        + f"\n\nWithout the {cats[3][1]} decoys and the {cats[2][1]} contaminant-only groups, "
        f"**{mq.n(n_clean)}** groups remain. This clean set is used for all abundance questions. "
        "Contaminant groups that contain a human protein, such as albumin, stay in it, "
        "because in serum they are genuine sample proteins."
    )
    return mq.answer(QUESTION, text, figure, short=f"{mq.n(total)} protein groups")


if __name__ == "__main__":
    print(main())
