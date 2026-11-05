"""Q2. How many peptides were found in total?

Counts the peptide sequences in peptides.txt, split into targets,
contaminants and decoys, and the peptides identified per sample
("Experiment <sample>" > 0, i.e. at least one evidence in that sample).
"""
import pandas as pd

import mq_common as mq

QUESTION = "How many peptides were found in total?"
TUTORIAL = "2,410 peptides"
NOTE = "MaxQuant version; the tutorial's own peptides table on Zenodo has 2,387 rows"


def main():
    pep = mq.load_peptides()
    dec, con = mq.is_decoy(pep), mq.is_contaminant(pep)
    total = len(pep)
    s = mq.samples(pep)
    colors = mq.sample_colors(s)
    nd = pep[~dec]
    present = pd.DataFrame({x: nd[mq.col(nd, f"Experiment {x}")].fillna(0).gt(0) for x in s})
    per_sample = {x: int(present[x].sum()) for x in s}
    in_all = int(present.all(axis=1).sum())

    fig, ax = mq.plt.subplots(figsize=(6, 3.4))
    names = s + ["In both samples"]
    vals = [per_sample[x] for x in s] + [in_all]
    bars = ax.bar(names, vals, color=[colors[x] for x in s] + [mq.NEUTRAL], width=0.6, **mq.bar_kw())
    mq.label_bars(ax, bars, [mq.n(v) for v in vals])
    ax.set_ylabel("Peptide sequences")
    ax.set_title(f"{mq.n(total)} peptide sequences in peptides.txt, identified per sample")
    ax.tick_params(axis="x", length=0)
    mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q02_total_peptides")

    rows = [
        ("Target", mq.n((~dec & ~con).sum())),
        ("Contaminant", mq.n((con & ~dec).sum())),
        ("Decoy", mq.n(dec.sum())),
        ("**Total**", f"**{mq.n(total)}**"),
    ]
    sample_txt = ", ".join(f"{x}: {mq.n(per_sample[x])}" for x in s)
    text = (
        f"**{mq.n(total)} peptide sequences** were found (all rows of peptides.txt).\n\n"
        + mq.md_table(rows, ["Category", "Peptides"])
        + f"\n\nIdentified per sample (decoys excluded): {sample_txt}; "
        f"{mq.n(in_all)} in both samples."
    )
    return mq.answer(QUESTION, text, figure, short=f"{mq.n(total)} peptides")


if __name__ == "__main__":
    print(main())
