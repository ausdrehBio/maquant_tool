"""Q3. How many proteins were identified and quantified?

Clean set (see mq_common). Identified = row in proteinGroups.txt.
Quantified = Intensity > 0 in at least one sample. LFQ intensities are
reported alongside: MaxLFQ needs at least one peptide quantified in two
samples, so proteins seen in only one sample get no LFQ value.
"""
import pandas as pd

import mq_common as mq

QUESTION = "How many proteins were identified and quantified?"
TUTORIAL = "Sample1: 237, Sample2: 123 (LFQ > 0)"
NOTE = ("tutorial counts LFQ > 0 in all protein groups; MaxQuant 2.9.0 gives no LFQ value to "
        "proteins seen in one sample only, so Intensity > 0 is counted here")


def main():
    cl = mq.clean(mq.load_protein_groups())
    s = mq.samples(cl)
    colors = mq.sample_colors(s)
    quant = pd.DataFrame({x: mq.intensity(cl, x) > 0 for x in s})
    has_lfq = pd.DataFrame({x: mq.lfq(cl, x) > 0 for x in s})

    identified = len(cl)
    quantified = int(quant.any(axis=1).sum())
    per_sample = {x: int(quant[x].sum()) for x in s}
    per_sample_lfq = {x: int(has_lfq[x].sum()) for x in s}
    lfq_any = int(has_lfq.any(axis=1).sum())
    in_all = quant.all(axis=1)
    only = {x: int((quant[x] & ~quant.drop(columns=x).any(axis=1)).sum()) for x in s}
    single = quant.any(axis=1) & ~in_all
    lfq_single = int((single & has_lfq.any(axis=1)).sum())

    fig, (a1, a2) = mq.plt.subplots(1, 2, figsize=(8, 3.4), gridspec_kw={"width_ratios": [1, 1]})
    names1 = ["Identified", "Quantified\n(Intensity > 0)", "With LFQ value"]
    vals1 = [identified, quantified, lfq_any]
    b1 = a1.bar(names1, vals1, color=mq.NEUTRAL, width=0.6, **mq.bar_kw())
    mq.label_bars(a1, b1, [mq.n(v) for v in vals1])
    a1.set_title("Protein groups (clean set)")
    names2 = ["In both"] + [f"Only {x}" for x in s]
    vals2 = [int(in_all.sum())] + [only[x] for x in s]
    b2 = a2.bar(names2, vals2, color=[mq.NEUTRAL] + [colors[x] for x in s], width=0.6, **mq.bar_kw())
    mq.label_bars(a2, b2, [mq.n(v) for v in vals2])
    a2.set_title("Quantified in which sample")
    for ax in (a1, a2):
        ax.set_ylabel("Protein groups")
        ax.tick_params(axis="x", length=0)
        mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q03_identified_quantified")

    only_txt = ", ".join(f"only in {x}: {mq.n(only[x])}" for x in s)
    rows = [(x, mq.n(per_sample[x]), mq.n(per_sample_lfq[x])) for x in s]
    text = (
        f"**{mq.n(identified)} protein groups were identified** (clean set) and "
        f"**{mq.n(quantified)} of them were quantified** (Intensity > 0 in at least one sample): "
        f"{mq.n(in_all.sum())} in both samples, {only_txt}.\n\n"
        + mq.md_table(rows, ["Sample", "Quantified (Intensity > 0)", "With LFQ intensity"])
        + f"\n\n{mq.n(lfq_any)} groups have an LFQ intensity. MaxLFQ is built from peptide ratios "
        f"between samples, so the {mq.n(single.sum())} groups quantified in only one sample have none "
        f"(groups seen in one sample with an LFQ value: {lfq_single})."
    )
    short = ", ".join(f"{x}: {mq.n(per_sample[x])}" for x in s) + " (Intensity > 0, clean set)"
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
