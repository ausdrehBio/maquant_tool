"""Q4. How similar are the technical replicates, and does the SDRF pairing of runs hold?

Pearson correlation of log2 LFQ intensities between runs (clean set, proteins quantified in
both runs). Each mouse kidney homogenate was split into two samples (technical replicates).
The pairing of runs to mice comes from PRIDE's community SDRF, which Table S1 of the paper
(supplement) cannot confirm here; e.g. runs 1 and 10 of each group are mouse 5. If the
pairing is right, a run should be most similar to its SDRF partner.

Rule used: the pairing is supported if the SDRF partner is the most similar run for at least
80% of the runs. MSstats (run_msstats.R) treats the two runs of a mouse as technical
replicates, so a pairing that is not supported must be resolved before the statistics.
"""
import itertools

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

import mq_common as mq

QUESTION = "How similar are the technical replicates, and does the SDRF pairing of runs hold?"
PAPER = "Two technical replicates per kidney, merged for statistics"
NOTE = "pairing from the PRIDE community SDRF"
SUPPORT_SHARE = 0.8
MIN_SHARED = 100  # proteins quantified in both runs for a correlation


def pairing_check():
    """Correlation matrix and per-run comparison of the SDRF partner with the other runs."""
    lfq = mq.log2_matrix(mq.load_protein_groups(), "LFQ intensity")
    corr = lfq.corr(min_periods=MIN_SHARED)
    rows = []
    for e in corr.columns:
        others = corr[e].drop(e).dropna()
        partner = mq.partner_of(e)
        same_group = [o for o in others.index if mq.condition_of(o) == mq.condition_of(e) and o != partner]
        best = others.idxmax()
        rows.append({
            "run": e, "partner": partner, "r_partner": corr.loc[e, partner],
            "r_group_median": others[same_group].median(), "best": best,
            "r_best": others.max(), "best_is_partner": best == partner,
        })
    res = pd.DataFrame(rows).set_index("run")
    best = res["best"].to_dict()
    mutual = sorted({tuple(sorted((a, b))) for a, b in best.items() if best.get(b) == a})
    supported = res["best_is_partner"].mean() >= SUPPORT_SHARE
    return corr, res, mutual, supported


def main():
    corr, res, mutual, supported = pairing_check()
    exps = mq.experiments()
    corr = corr.loc[exps, exps]
    n_best = int(res["best_is_partner"].sum())
    groups_pairs = [(a, b) for a, b in itertools.combinations(exps, 2) if mq.condition_of(a) == mq.condition_of(b)]
    r_partner = res["r_partner"]
    r_other = [corr.loc[a, b] for a, b in groups_pairs if mq.partner_of(a) != b]
    r_between = [corr.loc[a, b] for a, b in itertools.combinations(exps, 2) if mq.condition_of(a) != mq.condition_of(b)]

    fig, (ax1, ax2) = mq.plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
    vals = corr.values[~np.eye(len(exps), dtype=bool)]
    vmin = np.nanpercentile(vals, 2)
    cmap = LinearSegmentedColormap.from_list("seq", [mq.SURFACE, mq.SERIES[0], "#0d2f5c"])
    im = ax1.imshow(corr.values, cmap=cmap, vmin=vmin, vmax=1, interpolation="nearest")
    ax1.set_xticks(range(len(exps)), exps)
    ax1.set_yticks(range(len(exps)), exps)
    mq.rotate_xticks(ax1, size=5, angle=90)
    ax1.tick_params(axis="y", length=0, labelsize=5)
    for k in range(1, len(mq.CONDITIONS)):
        ax1.axhline(k * 10 - 0.5, color=mq.SURFACE, linewidth=1.5)
        ax1.axvline(k * 10 - 0.5, color=mq.SURFACE, linewidth=1.5)
    for e in exps:  # outline the SDRF partner cells
        i, j = exps.index(e), exps.index(mq.partner_of(e))
        ax1.add_patch(mq.plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=mq.SERIES[1], linewidth=1))
    for spine in ax1.spines.values():
        spine.set_visible(False)
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.02, label="Pearson r (log2 LFQ)")
    ax1.set_title("Run-to-run correlation (orange: SDRF partner)")

    colors = mq.condition_colors()
    for i, e in enumerate(exps):
        others = [corr.loc[e, o] for o in exps if o != e and mq.condition_of(o) == mq.condition_of(e) and o != mq.partner_of(e)]
        ax2.scatter([i] * len(others), others, s=10, color=mq.AXIS, edgecolor="none", zorder=2)
        ax2.scatter([i], [res.loc[e, "r_partner"]], s=26, color=colors[mq.condition_of(e)],
                    edgecolor=mq.INK if not res.loc[e, "best_is_partner"] else "none", linewidth=0.8, zorder=3)
    ax2.set_xticks(range(len(exps)), exps)
    mq.rotate_xticks(ax2, size=6, angle=90)
    ax2.set_ylabel("Pearson r")
    ax2.set_title("SDRF partner (coloured) vs other runs of the group (grey)")
    mq.value_grid(ax2)
    fig.tight_layout()
    figure = mq.save_fig(fig, "q04_technical_replicates")

    not_best = res[~res["best_is_partner"]]
    detail = ", ".join(
        f"{e} (partner {r.partner}, r = {r.r_partner:.3f}; most similar {r.best}, r = {r.r_best:.3f})"
        for e, r in not_best.iterrows()
    )
    verdict = (
        "**The SDRF pairing is supported by the data.**" if supported else
        "**The SDRF pairing is not supported by the data**; resolve it before the statistics."
    )
    text = (
        f"For **{n_best} of {len(exps)} runs** the most similar run is the SDRF partner "
        f"(rule: at least {SUPPORT_SHARE:.0%}). {verdict} Median Pearson r of log2 LFQ: SDRF "
        f"partners {np.median(r_partner):.3f}, other runs of the same group "
        f"{np.median(r_other):.3f}, runs of different groups {np.median(r_between):.3f}. "
        f"Mutually most similar pairs in the data: {len(mutual)}, of which "
        f"{sum(mq.partner_of(a) == b for a, b in mutual)} are SDRF pairs."
        + (f"\n\nRuns whose most similar run is not the SDRF partner: {detail}." if len(not_best) else "")
    )
    short = f"partner most similar for {n_best}/{len(exps)} runs; median r {np.median(r_partner):.3f}"
    return mq.answer(QUESTION, text, figure, short)


if __name__ == "__main__":
    print(main())
