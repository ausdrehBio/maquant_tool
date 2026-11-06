"""Q5. How many proteins were included into the statistical analysis?

Proteins that MSstats summarised and compared (ProteinLevelData / ComparisonResult).
MSstats works on the leading razor protein of each feature, not on MaxQuant's protein
groups, and drops proteins with only one peptide. The figure follows the count from
MaxQuant's protein groups to the compared proteins.
"""
import numpy as np

import mq_common as mq

QUESTION = "How many proteins were included into the statistical analysis?"
TUTORIAL = "1693 (MSstats log)"
NOTE = "depends on the MaxQuant and MSstats versions"


def main():
    pg = mq.load_protein_groups()
    human = int(mq.has_human(pg[mq.col(pg, "Protein IDs")]).sum())
    feat = mq.load_msstats("FeatureLevelData")
    prot = mq.load_msstats("ProteinLevelData")
    cr = mq.load_msstats("ComparisonResult")
    n_input = feat[mq.col(feat, "PROTEIN")].nunique()
    n_summ = prot[mq.col(prot, "Protein")].nunique()
    fc = cr[mq.col(cr, "log2FC")]
    n_tested = cr.loc[fc.notna(), mq.col(cr, "Protein")].nunique()

    steps = [("Protein groups (MaxQuant)", len(pg)), ("with a human protein", human),
             ("MSstats input (≥ 2 peptides)", n_input), ("summarised by dataProcess", n_summ),
             ("with a fold change", n_tested)]
    fig, ax = mq.plt.subplots(figsize=(7, 2.8))
    names, vals = [s for s, _ in steps][::-1], [v for _, v in steps][::-1]
    b = ax.barh(names, vals, color=mq.NEUTRAL, height=0.6, **mq.bar_kw())
    mq.label_bars(ax, b, [mq.n(v) for v in vals])
    ax.set_xlabel("Proteins")
    ax.set_title(f"{mq.n(n_summ)} proteins in the statistical analysis")
    ax.tick_params(axis="y", length=0)
    mq.value_grid(ax, "x")
    figure = mq.save_fig(fig, "q05_proteins_in_statistics")

    text = (
        f"**{mq.n(n_summ)} proteins** were summarised by MSstats and entered the group "
        f"comparison; {mq.n(n_tested)} of them got a fold change. MSstats counts leading razor "
        f"proteins, not protein groups: from {mq.n(len(pg))} MaxQuant protein groups, "
        f"{mq.n(human)} contain a human protein, and {mq.n(n_input)} leading razor proteins with "
        "at least two peptides went into MSstats."
    )
    return mq.answer(QUESTION, text, figure, short=f"{mq.n(n_summ)}")


if __name__ == "__main__":
    print(main())
