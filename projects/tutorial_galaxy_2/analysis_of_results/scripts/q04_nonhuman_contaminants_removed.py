"""Q4. How many proteins were removed as potential non human contaminants?

Tutorial rule: keep only protein group lines that contain "HUMAN". Groups without any
human protein are removed (they are all on MaxQuant's contaminant list). The count is
cross-checked against run_msstats.R, which applies the same rule.
"""
import numpy as np

import mq_common as mq

QUESTION = "How many proteins were removed as potential non human contaminants?"
TUTORIAL = "27 (2518 lines in the protein group file minus 2491 lines after select)"
NOTE = "same rule: keep lines containing HUMAN"
TOP_N = 15


def main():
    pg = mq.load_protein_groups()
    ids = pg[mq.col(pg, "Protein IDs")].fillna("")
    removed = ~mq.has_human(ids)
    n_removed = int(removed.sum())
    from_r = int(mq.msstats_info()["protein_groups_without_human_protein"])
    assert n_removed == from_r, f"Python {n_removed} vs run_msstats.R {from_r}"
    all_con = bool(mq.is_contaminant(pg)[removed].all())

    tot = pg[mq.col(pg, "Intensity")].fillna(0).astype(float)
    top = tot[removed].nlargest(TOP_N).index
    labels = [ids[i].split(";")[0] for i in top][::-1]
    fig, ax = mq.plt.subplots(figsize=(7, 0.3 * len(top) + 1.2))
    vals = np.log10(tot[top].replace(0, np.nan))[::-1]
    ax.barh(labels, vals, color=mq.NEUTRAL, height=0.6, **mq.bar_kw())
    ax.set_xlabel("log10 total intensity")
    ax.set_title(f"{n_removed} removed groups without a human protein (top {len(top)} by intensity)")
    ax.tick_params(axis="y", length=0, labelsize=8)
    mq.value_grid(ax, "x")
    figure = mq.save_fig(fig, "q04_nonhuman_contaminants_removed")

    text = (
        f"**{n_removed} protein groups** were removed: {mq.n(len(pg))} groups before, "
        f"{mq.n(len(pg) - n_removed)} with a human protein after. "
        + ("All of them are on MaxQuant's contaminant list (CON__). " if all_con else
           "Not all of them are flagged as contaminants. ")
        + "Contaminant groups that contain a human protein stay, and their contaminant flag "
        "is cleared before MSstats, so MSstats does not remove them either."
    )
    return mq.answer(QUESTION, text, figure, short=f"{n_removed}")


if __name__ == "__main__":
    print(main())
