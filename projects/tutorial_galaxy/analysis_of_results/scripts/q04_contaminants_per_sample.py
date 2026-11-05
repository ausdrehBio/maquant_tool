"""Q4. In which sample were more contaminants quantified?

Contaminant = MaxQuant's "Potential contaminant" flag (decoys excluded).
Main measure: contaminant share of the sample's total intensity, which is
what the PTXQC report shows. Also reported: the share without albumin
(on MaxQuant's contaminant list, but the main protein of serum) and the
number of contaminant groups with Intensity > 0.
"""
import numpy as np
from matplotlib.ticker import MaxNLocator

import mq_common as mq

QUESTION = "In which sample were more contaminants quantified?"
TUTORIAL = "Sample 2, especially serum albumin (from the PTXQC report)"
NOTE = "both by contaminant share of intensity, the measure PTXQC reports"


def main():
    pg = mq.load_protein_groups()
    nd = pg[~mq.is_decoy(pg)]
    con, alb = mq.is_contaminant(nd), mq.is_albumin(nd)
    s = mq.samples(nd)
    colors = mq.sample_colors(s)
    count, share_all, share_noalb = {}, {}, {}
    for x in s:
        inten = mq.intensity(nd, x)
        count[x] = int((con & (inten > 0)).sum())
        share_all[x] = inten[con].sum() / inten.sum() * 100
        share_noalb[x] = inten[con & ~alb].sum() / inten.sum() * 100

    fig, (a_share, a_count) = mq.plt.subplots(1, 2, figsize=(8, 3.4))
    groups = ["All contaminants", "Without albumin"]
    pos = np.arange(len(groups))
    width = 0.8 / len(s)
    for i, x in enumerate(s):
        vals = [share_all[x], share_noalb[x]]
        b = a_share.bar(pos + (i - (len(s) - 1) / 2) * width, vals, width=width, color=colors[x],
                        label=x, **mq.bar_kw())
        mq.label_bars(a_share, b, [mq.pct(v) for v in vals])
    a_share.set_xticks(pos, groups)
    a_share.set_ylabel("% of sample intensity")
    a_share.set_title("Contaminant share of intensity (PTXQC measure)")
    a_share.legend(loc="upper right")
    b1 = a_count.bar(s, [count[x] for x in s], color=[colors[x] for x in s], width=0.6, **mq.bar_kw())
    mq.label_bars(a_count, b1, [mq.n(count[x]) for x in s])
    a_count.yaxis.set_major_locator(MaxNLocator(integer=True))
    a_count.set_ylabel("Contaminant groups quantified")
    a_count.set_title("Number of contaminant groups")
    for ax in (a_share, a_count):
        ax.tick_params(axis="x", length=0)
        mq.value_grid(ax, "y")
    fig.tight_layout()
    figure = mq.save_fig(fig, "q04_contaminants_per_sample")

    hi_int = max(s, key=share_all.get)
    hi_cnt = max(s, key=count.get)
    by_share = sorted(s, key=share_all.get, reverse=True)
    share_txt = f"{mq.pct(share_all[hi_int])} of the intensity in {hi_int}, against " + ", ".join(
        f"{mq.pct(share_all[x])} in {x}" for x in by_share[1:])
    count_txt = " vs ".join(f"{count[x]} in {x}" for x in s)
    dep = mq.depleted_sample(pg)[0]
    alb_frac = (1 - share_noalb[hi_int] / share_all[hi_int]) * 100
    count_lead = ("the same sample also has more contaminant groups" if hi_cnt == hi_int else
                  f"{hi_cnt} has more contaminant groups")
    reason = (": the depleted sample reveals more low-abundance proteins, contaminants included"
              if hi_cnt == dep and hi_cnt != hi_int else "")
    rows = [(x, mq.pct(share_all[x]), mq.pct(share_noalb[x], 2), count[x]) for x in s]
    text = (
        f"**{hi_int}**: contaminants make up {share_txt} (the measure PTXQC reports).\n\n"
        + mq.md_table(rows, ["Sample", "% intensity", "% intensity without albumin", "Contaminant groups"])
        + f"\n\nIn {hi_int}, albumin accounts for {alb_frac:.0f}% of the contaminant intensity. MaxQuant "
        "flags it as a contaminant, but it is a genuine serum protein. Without albumin, contaminants "
        "(keratins, trypsin, bovine proteins) make up "
        + " and ".join(f"{mq.pct(share_noalb[x], 2)} in {x}" for x in s)
        + f". By number, {count_lead} ({count_txt}){reason}."
    )
    short = (f"{hi_int} by intensity ({share_txt}); by number of groups {hi_cnt} ({count_txt})")
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
