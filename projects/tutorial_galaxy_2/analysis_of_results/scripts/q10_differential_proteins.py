"""Q10. How many proteins are differentially abundant?

Tutorial definition: adjusted p-value above 0 and below 0.05 (an adjusted p-value of 0
means the protein is missing in one condition) and log2 fold change beyond +-0.58.
The heatmaps show these proteins' summarised intensities per sample (MSstats sample
quantification), scaled per protein, without clustering, as in the tutorial. Gene and
protein names come from the project's FASTA file (the tutorial retrieves them from
UniProt). The three proteins validated by immunostaining in the original publication
(Foell et al. 2018) are checked at the end.
"""
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

import mq_common as mq

QUESTION = "How many proteins are differentially abundant?"
TUTORIAL = "125 up in metastasizing cSCC, 14 up in RDEB cSCC"
NOTE = "adjusted p-value > 0 and < 0.05, log2FC beyond ±0.58"
# proteins stained in the publication: expected higher in which condition
PUBLICATION = {"Q05707": ("COL14A1", "RDEB"), "P02743": ("APCS", "metast"),
               "P12956": ("XRCC6", "metast")}


def sample_matrix():
    """Sample quantification matrix with columns renamed to sample names."""
    m = mq.load_msstats("SampleQuantificationMatrix")
    a = mq._annotation()
    key = (a[mq.col(a, "Condition")].astype(str) + "_" + a[mq.col(a, "BioReplicate")].astype(str))
    rename = dict(zip(key, a["sample"]))
    m = m.rename(columns=rename).set_index(mq.col(m, "Protein"))
    return m[[s for s in mq.samples() if s in m.columns]]


def main():
    cr = mq.load_msstats("ComparisonResult")
    name, up_cond, down_cond = mq.contrast()
    prot_col = mq.col(cr, "Protein")
    fc = cr[mq.col(cr, "log2FC")]
    adj = cr[mq.col(cr, "adj.pvalue")]
    diff = (adj > 0) & (adj < mq.P_CUTOFF) & np.isfinite(fc)
    up = cr[diff & (fc > mq.FC_CUTOFF)].sort_values(mq.col(cr, "log2FC"), ascending=False)
    down = cr[diff & (fc < -mq.FC_CUTOFF)].sort_values(mq.col(cr, "log2FC"))

    mat = sample_matrix()
    cmap = LinearSegmentedColormap.from_list("div", mq.DIVERGING)
    sets = [(up, f"Higher in {up_cond}"), (down, f"Higher in {down_cond}")]
    heights = [max(len(s), 1) for s, _ in sets]
    fig, axes = mq.plt.subplots(2, 1, figsize=(8, 0.11 * sum(heights) + 3),
                                gridspec_kw={"height_ratios": heights})
    colors = mq.condition_colors()
    for ax, (subset, title) in zip(axes, sets):
        rows = mat.reindex(subset[prot_col])
        z = rows.sub(rows.mean(axis=1), axis=0).div(rows.std(axis=1).replace(0, np.nan), axis=0)
        im = ax.imshow(z.to_numpy(dtype=float), aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5,
                       interpolation="nearest")
        ax.set_yticks(range(len(z)), [mq.gene_of(p) for p in z.index],
                      fontsize=max(4.5, min(8, 400 / max(len(z), 1))))
        ax.set_xticks(range(z.shape[1]), z.columns)
        mq.rotate_xticks(ax, size=7)
        for tick, s in zip(ax.get_xticklabels(), z.columns):
            tick.set_color(colors[mq.condition_of(s)])
        ax.set_title(f"{title}: {len(z)} proteins")
        ax.tick_params(length=0)
        for side in ax.spines.values():
            side.set_visible(False)
    cbar = fig.colorbar(im, ax=axes, shrink=0.4, pad=0.02)
    cbar.set_label("z-score per protein", color=mq.INK2)
    figure = mq.save_fig(fig, "q10_differential_proteins")

    acc = cr[prot_col].map(mq.accession)
    pub_rows = []
    for a, (gene, expected) in PUBLICATION.items():
        hit = cr[acc == a]
        if hit.empty:
            pub_rows.append((gene, f"higher in {expected}", "not quantified", "", ""))
            continue
        f, q = hit[mq.col(cr, "log2FC")].iloc[0], hit[mq.col(cr, "adj.pvalue")].iloc[0]
        found = up_cond if f > 0 else down_cond
        ok = "yes" if (found == expected and q < mq.P_CUTOFF) else "no"
        pub_rows.append((gene, f"higher in {expected}", f"{f:.2f}", f"{q:.2g}", ok))

    def listing(subset, k=10):
        genes = [mq.gene_of(p) for p in subset[prot_col]]
        more = f" and {len(genes) - k} more" if len(genes) > k else ""
        return ", ".join(genes[:k]) + more

    text = (
        f"**{len(up)} proteins are higher in {up_cond}** and **{len(down)} higher in "
        f"{down_cond}** (adjusted p between 0 and 0.05, more than 1.5-fold).\n\n"
        f"Strongest in {up_cond}: {listing(up)}.\n\nStrongest in {down_cond}: {listing(down)}.\n\n"
        "Proteins validated by immunostaining in the original publication:\n\n"
        + mq.md_table(pub_rows, ["Gene", "Publication", f"log2FC ({name})", "adj. p", "Confirmed"])
    )
    short = f"{len(up)} up in {up_cond}, {len(down)} up in {down_cond}"
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
