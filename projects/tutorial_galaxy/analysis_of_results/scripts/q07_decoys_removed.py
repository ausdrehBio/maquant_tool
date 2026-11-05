"""Q7. How many decoy proteins were removed?

Decoys are reversed protein sequences (REV__) that MaxQuant searches to
estimate the false discovery rate. A decoy group is marked "+" in the Decoy
column. A group that merely contains a REV__ member next to a target
protein is not a decoy and is kept.
"""
import numpy as np

import mq_common as mq

QUESTION = "How many decoy proteins were removed?"
TUTORIAL = "2"
NOTE = ""


def main():
    pg = mq.load_protein_groups()
    dec = mq.is_decoy(pg)
    ids = pg[mq.col(pg, "Protein IDs")].fillna("")
    mixed = ids.str.contains("REV__", regex=False) & ~dec
    score = pg[mq.col(pg, "Score")]
    qval = pg[mq.col(pg, "Q-value")]
    total = pg[mq.col(pg, "Intensity")].fillna(0).astype(float)
    lab = mq.labels(pg)

    fig, ax = mq.plt.subplots(figsize=(7.5, 4.2))
    ok = total > 0
    target = ok & ~dec
    box = {"boxstyle": "round,pad=0.2", "fc": mq.SURFACE, "ec": "none", "alpha": 0.85}
    ax.scatter(score[target], np.log2(total[target]), s=40, color=mq.MUTED, alpha=0.55,
               edgecolor=mq.SURFACE, linewidth=0.8, label="Target and contaminant groups")
    ax.scatter(score[dec & ok], np.log2(total[dec & ok]), s=70, color=mq.INK,
               edgecolor=mq.SURFACE, linewidth=1.5, label="Decoy groups (removed)", zorder=3)
    for k, i in enumerate(pg.index[dec & ok]):
        ax.annotate(ids[i], (score[i], np.log2(total[i])), xytext=(60, 18 - 30 * k),
                    textcoords="offset points", va="center", fontsize=9, color=mq.INK, bbox=box,
                    arrowprops={"arrowstyle": "-", "color": mq.INK2, "linewidth": 0.8})
    for i in pg.index[mixed & ok]:
        ax.scatter(score[i], np.log2(total[i]), s=90, facecolor="none", edgecolor=mq.INK,
                   linewidth=1.2, zorder=4)
        ax.annotate(f"{lab[i]}: contains a REV__ member, kept", (score[i], np.log2(total[i])),
                    xytext=(-40, 40), textcoords="offset points", ha="right", fontsize=9,
                    color=mq.INK2, bbox=box,
                    arrowprops={"arrowstyle": "-", "color": mq.INK2, "linewidth": 0.8})
    ax.set_xlabel("Protein group score (Andromeda)")
    ax.set_ylabel("log2 total intensity")
    ax.set_title(f"{int(dec.sum())} decoy groups, highest score {score[dec].max():.1f} "
                 f"(median target score {score[~dec].median():.0f})")
    ax.legend(loc="lower right")
    mq.value_grid(ax, "both")
    figure = mq.save_fig(fig, "q07_decoys_removed")

    rows = [(ids[i], f"{score[i]:.2f}", f"{qval[i]:.4f}", mq.sci(total[i])) for i in pg.index[dec]]
    mixed_txt = "; ".join(f"`{ids[i]}` ({lab[i]})" for i in pg.index[mixed])
    text = (
        f"**{int(dec.sum())} decoy protein groups were removed.**\n\n"
        + mq.md_table(rows, ["Protein IDs", "Score", "Q-value", "Total intensity"])
        + f"\n\nNot removed: {int(mixed.sum())} group that contains a REV__ member next to a target "
        f"protein, {mixed_txt}. Its leading protein is a real human protein, and MaxQuant does not "
        "mark it as a decoy. Removing every row with REV__ anywhere in Protein IDs would wrongly "
        "drop it."
    )
    return mq.answer(QUESTION, text, figure, short=str(int(dec.sum())))


if __name__ == "__main__":
    print(main())
