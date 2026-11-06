"""Q5. Do the four groups separate in a principal component analysis?

Mouse level: log2 LFQ intensities averaged over the two technical replicates, protein groups
with a value in all 20 mice, each protein centred and scaled to unit variance (as SIMCA's
default UV scaling in the paper). Components from a singular value decomposition.
A group is called separated on a component when all its scores lie on one side of all
scores of the other group(s).
"""
import itertools

import numpy as np

import mq_common as mq

QUESTION = "Do the four groups separate in a principal component analysis?"
PAPER = ("PCA (PC1 vs PC3): healthy control separated from the three db/db groups, "
         "which overlap slightly")
NOTE = "mouse level, complete protein groups, UV scaling"
N_PC = 3


def separated(scores, a, b):
    return scores[a].max() < scores[b].min() or scores[a].min() > scores[b].max()


def main():
    mice = mq.mouse_matrix(mq.log2_matrix(mq.load_protein_groups(), "LFQ intensity"))
    complete = mice.dropna()
    x = complete.T.values  # mice x proteins
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :N_PC] * s[:N_PC]
    explained = s**2 / np.sum(s**2) * 100
    cond = np.array([mq.condition_of_mouse(m) for m in complete.columns])

    dbdb = cond != "CTRL"
    ctrl_sep = [k + 1 for k in range(N_PC) if separated(scores[:, k], cond == "CTRL", dbdb)]
    pair_sep = {}
    for a, b in itertools.combinations(["VEH", "LIRA", "INS"], 2):
        pair_sep[(a, b)] = [k + 1 for k in range(N_PC) if separated(scores[:, k], cond == a, cond == b)]

    colors = mq.condition_colors()
    fig, axes = mq.plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, (i, j) in zip(axes, [(0, 1), (0, 2)]):
        for c in mq.CONDITIONS:
            m = cond == c
            ax.scatter(scores[m, i], scores[m, j], s=46, color=colors[c], edgecolor=mq.SURFACE,
                       linewidth=0.8, label=f"{c}: {mq.CONDITION_LABELS[c]}", zorder=3)
        ax.axhline(0, color=mq.GRID, linewidth=0.8, zorder=1)
        ax.axvline(0, color=mq.GRID, linewidth=0.8, zorder=1)
        ax.set_xlabel(f"PC{i + 1} ({explained[i]:.1f}%)")
        ax.set_ylabel(f"PC{j + 1} ({explained[j]:.1f}%)")
        ax.set_title(f"PC{i + 1} vs PC{j + 1}")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    fig.tight_layout()
    figure = mq.save_fig(fig, "q05_pca")

    def pcs(ks):
        return "PC" + " and PC".join(str(k) for k in ks) if ks else "none of PC1 to PC3"

    pair_text = "; ".join(f"{a} vs {b}: {pcs(ks)}" for (a, b), ks in pair_sep.items())
    text = (
        f"PCA of {mq.n(len(complete))} protein groups quantified in all 20 mice. PC1 to PC3 "
        f"explain {explained[0]:.1f}%, {explained[1]:.1f}% and {explained[2]:.1f}% of the "
        f"variance. **The healthy control group (CTRL) is separated from the three db/db "
        f"groups on {pcs(ctrl_sep)}.** Separation between the db/db groups: {pair_text}."
    )
    short = f"CTRL separated from db/db on {pcs(ctrl_sep)}; db/db pairs: {pair_text}"
    return mq.answer(QUESTION, text, figure, short)


if __name__ == "__main__":
    print(main())
