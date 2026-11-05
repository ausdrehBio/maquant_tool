"""Q5. How good was the tryptic digestion (percentage of zero missed cleavages)?

PTXQC method: identified MS/MS spectra (msms.txt) per raw file, without
decoys and contaminants; share with 0 missed cleavages. peptides.txt gives
the same share over all unique peptide sequences.

Rules of thumb for the verdict, applied to the lowest raw file:
at least 90% fully cleaved = good; at least 70% = not ideal but good enough;
below 70% = incomplete digestion.
"""
import pandas as pd

import mq_common as mq

QUESTION = "How good was the tryptic digestion (percentage of zero missed cleavages)?"
TUTORIAL = "75% (Sample1), around 85% (Sample2); not ideal but good enough"
NOTE = "same measure as the PTXQC report the tutorial reads"
GOOD_MIN_PCT = 90
OK_MIN_PCT = 70


def main():
    ms = mq.load_msms()
    ms = ms[~mq.is_decoy(ms) & ~mq.is_contaminant(ms)]
    raw, exp, mc = mq.col(ms, "Raw file"), mq.col(ms, "Experiment"), mq.col(ms, "Missed cleavages")
    counts = pd.crosstab(ms[raw], ms[mc])
    experiment = ms.groupby(raw)[exp].first()
    psms = counts.sum(axis=1)
    zero_pct = counts.get(0, 0) / psms * 100

    pep = mq.load_peptides()
    pep = pep[~mq.is_decoy(pep) & ~mq.is_contaminant(pep)]
    pep_mc = pep[mq.col(pep, "Missed cleavages")]
    pep_pct = (pep_mc == 0).mean() * 100
    max_mc = mq.mqpar_value("maxMissedCleavages")

    colors = mq.sample_colors(mq.samples(mq.load_protein_groups()))
    fig, ax = mq.plt.subplots(figsize=(6, 3.6))
    files = list(counts.index)
    bars = ax.bar(files, [zero_pct[f] for f in files], width=0.6,
                  color=[colors.get(experiment[f], mq.NEUTRAL) for f in files], **mq.bar_kw())
    ax.bar_label(bars, labels=[mq.pct(zero_pct[f]) for f in files], padding=-16, color="#ffffff",
                 fontsize=9, fontweight="bold")
    ax.axhline(pep_pct, color=mq.INK2, linewidth=1)
    ax.text(1.01, pep_pct, f"all peptides\n{mq.pct(pep_pct)}", transform=ax.get_yaxis_transform(),
            va="center", ha="left", color=mq.INK2, fontsize=9)
    ax.set_ylim(0, 100)
    ax.set_ylabel("% PSMs with 0 missed cleavages")
    ax.set_title("Fully cleaved peptide-spectrum matches per raw file")
    ax.tick_params(axis="x", length=0)
    mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q05_missed_cleavages")

    lowest = zero_pct.min()
    if lowest >= GOOD_MIN_PCT:
        verdict_short = "good"
        verdict = f"The digestion was **good**: every raw file is at least {GOOD_MIN_PCT}% fully cleaved"
    elif lowest >= OK_MIN_PCT:
        verdict_short = "not ideal but good enough"
        verdict = (f"The digestion was **not ideal but good enough**: up to {100 - lowest:.0f}% of the "
                   "identified spectra still contain a missed cleavage")
    else:
        verdict_short = "incomplete"
        verdict = (f"The digestion was **incomplete**: in at least one raw file fewer than "
                   f"{OK_MIN_PCT}% of the spectra are fully cleaved")
    spread = zero_pct.max() - zero_pct.min()
    similar = (f"The raw files are close together ({spread:.1f} percentage points apart)."
               if spread < 10 else
               f"The raw files differ by {spread:.1f} percentage points.")
    mc_word = "cleavage" if max_mc == "1" else "cleavages"
    rows =[(f, experiment[f], mq.n(psms[f]), mq.n(counts.get(0, 0)[f]), mq.pct(zero_pct[f])) for f in files]
    text = (
        mq.md_table(rows, ["Raw file", "Sample", "PSMs", "0 missed cleavages", "%"])
        + f"\n\nOver all {mq.n(len(pep))} target peptide sequences (peptides.txt): "
        f"{mq.pct(pep_pct)} without missed cleavage.\n\n"
        f"{verdict}. {similar} Caveat: this search allowed at most {max_mc} missed {mc_word}, "
        "so peptides with more were not identified at all, and the true share of fully cleaved "
        "peptides is somewhat lower."
    )
    short = ", ".join(f"{mq.pct(zero_pct[f])} ({experiment[f]})" for f in files) + f"; {verdict_short}"
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
