"""Q6. How many proteins differ between the groups?

MSstats groupComparison (run_msstats.R) for the contrasts of parameter/comparison_matrix.tsv:
VEH-CTRL (db/db vs db/+, both vehicle), LIRA-VEH and INS-VEH (treatment effects within db/db)
and LIRA-INS. Significant: adjusted p-value < 0.05 (Benjamini-Hochberg within a contrast);
differential: significant and log2FC beyond +-0.58 (1.5-fold). Proteins missing in one group
get an infinite fold change; they are counted but not drawn.
"""
import numpy as np

import mq_common as mq

QUESTION = "How many proteins differ between the groups?"
PAPER = ("Table 1: 23 proteins with significant levels (one-way ANOVA with Tukey, shotgun "
         "and SRM), most of them db/db vs db/+")
NOTE = "MSstats, adj. p < 0.05 and > 1.5-fold"
TOP = 8


def main():
    cr = mq.load_msstats("ComparisonResult")
    label_col, fc_col, adj_col = mq.col(cr, "Label"), mq.col(cr, "log2FC"), mq.col(cr, "adj.pvalue")
    prot_col = mq.col(cr, "Protein")
    colors = mq.condition_colors()

    rows, tops = [], []
    fig, axes = mq.plt.subplots(2, 2, figsize=(10.5, 8.4))
    for ax, (label, num, den) in zip(axes.flat, mq.contrasts()):
        d = cr[cr[label_col] == label]
        fc, adj = d[fc_col].astype(float), d[adj_col].astype(float)
        tested = int(adj.notna().sum())
        sig = adj < mq.P_CUTOFF
        up, down = sig & (fc > mq.FC_CUTOFF), sig & (fc < -mq.FC_CUTOFF)
        infinite = sig & np.isinf(fc)
        rows.append((label, f"{num} vs {den}", mq.n(tested), mq.n(sig.sum()), mq.n(up.sum()),
                     mq.n(down.sum()), mq.n(infinite.sum())))

        def names(mask):
            sub = d[mask].assign(a=adj[mask], f=fc[mask].abs()).sort_values(["a", "f"], ascending=[True, False])
            return [mq.gene_of(p) for p in sub[prot_col].head(TOP)]
        tops.append((label, num, den, names(up), names(down)))

        finite = np.isfinite(fc) & adj.notna()
        y = -np.log10(adj.clip(lower=1e-300))
        for mask, color in [(finite & ~up & ~down, mq.AXIS), (finite & up, colors[num]), (finite & down, colors[den])]:
            ax.scatter(fc[mask], y[mask], s=10, color=color, alpha=0.8, edgecolor="none")
        ax.axhline(-np.log10(mq.P_CUTOFF), color=mq.INK2, linewidth=0.7)
        for x in (-mq.FC_CUTOFF, mq.FC_CUTOFF):
            ax.axvline(x, color=mq.INK2, linewidth=0.7)
        ax.set_title(f"{label}: {int(up.sum())} up, {int(down.sum())} down")
        ax.set_xlabel(f"log2 fold change ({num} / {den})")
        ax.set_ylabel("-log10 adjusted p-value")
        mq.value_grid(ax, "both")
    fig.tight_layout()
    figure = mq.save_fig(fig, "q06_differential")

    top_text = "\n".join(
        f"- **{label}**: higher in {num}: {', '.join(u) or 'none'}; higher in {den}: {', '.join(dn) or 'none'}"
        for label, num, den, u, dn in tops
    )
    text = (
        "Proteins with adjusted p < 0.05 and more than 1.5-fold change (MSstats, technical "
        "replicates modelled per mouse):\n\n"
        + mq.md_table(rows, ["Contrast", "Comparison", "Proteins tested", "adj. p < 0.05",
                             "Up (> 1.5-fold)", "Down (> 1.5-fold)", "Missing in one group"])
        + f"\n\nStrongest differential proteins (lowest adjusted p, up to {TOP} per side):\n\n" + top_text
    )
    short = "; ".join(f"{r[0]}: {r[4]} up, {r[5]} down" for r in rows)
    return mq.answer(QUESTION, text, figure, short)


if __name__ == "__main__":
    print(main())
