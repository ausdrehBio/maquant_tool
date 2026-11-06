"""Q2. How complete is the quantification across the 40 runs?

LFQ intensity > 0 per run and per protein group (clean set). Share of evidence rows that
match between runs transferred into a run (evidence Type MULTI-MATCH), per run.
"""
import mq_common as mq

QUESTION = "How complete is the quantification across the 40 runs?"
PAPER = "Not reported (Progenesis quantifies every aligned feature in every run)"
NOTE = "LFQ intensity > 0; match between runs on"


def main():
    pg = mq.load_protein_groups()
    lfq = mq.log2_matrix(pg, "LFQ intensity")
    inten = mq.log2_matrix(pg, "Intensity")
    n_clean = len(lfq)
    per_run = lfq.notna().sum()
    runs_per_protein = lfq.notna().sum(axis=1)
    complete = int((runs_per_protein == len(lfq.columns)).sum())
    any_lfq = int((runs_per_protein > 0).sum())
    no_lfq = n_clean - any_lfq
    mice = mq.mouse_matrix(lfq)
    complete_mice = int(mice.notna().all(axis=1).sum())

    ev = mq.load_evidence()
    ev = ev[~mq.is_decoy(ev)]
    lookup = dict(zip(mq.annotation()["Raw.file"], mq.annotation()["Experiment"]))
    ev_exp = ev[mq.col(ev, "Raw file")].map(lookup)
    is_match = ev[mq.col(ev, "Type")].astype(str).str.upper().eq("MULTI-MATCH")
    match_share = (is_match.groupby(ev_exp).mean() * 100).reindex(mq.experiments())

    fig, axes = mq.plt.subplots(2, 1, figsize=(10, 6.6))
    mq.run_bars(axes[0], per_run, "protein groups with LFQ")
    axes[0].set_title("Protein groups with LFQ intensity per run")
    mq.condition_legend(axes[0], loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    counts = runs_per_protein.value_counts().reindex(range(len(lfq.columns) + 1), fill_value=0)
    axes[1].bar(counts.index, counts.values, color=mq.NEUTRAL, **mq.bar_kw())
    axes[1].set_xlabel("number of runs with LFQ intensity")
    axes[1].set_ylabel("protein groups")
    axes[1].set_title("Protein groups by number of runs with LFQ intensity")
    mq.value_grid(axes[1])
    fig.tight_layout()
    figure = mq.save_fig(fig, "q02_quantification")

    by_cond = per_run.groupby([mq.condition_of(e) for e in per_run.index]).mean().loc[mq.CONDITIONS]
    match_cond = match_share.groupby([mq.condition_of(e) for e in match_share.index]).mean().loc[mq.CONDITIONS]
    rows = [(c, mq.n(round(by_cond[c])), mq.pct(match_cond[c])) for c in mq.CONDITIONS]
    text = (
        f"Of the {mq.n(n_clean)} protein groups of the clean set, **{mq.n(any_lfq)} have an LFQ "
        f"intensity** in at least one run and **{mq.n(complete)} in all 40 runs**; "
        f"{mq.n(no_lfq)} have none (MaxLFQ needs peptides quantified in at least two runs). "
        f"Per run, {mq.n(per_run.min())} to {mq.n(per_run.max())} protein groups have an LFQ "
        f"intensity (median {mq.n(per_run.median())}); with Intensity, which needs no ratio, "
        f"the median is {mq.n(inten.notna().sum().median())}. After averaging the technical "
        f"replicates, {mq.n(complete_mice)} protein groups have a value in all 20 mice.\n\n"
        f"Match between runs contributed {mq.pct(match_share.min())} to "
        f"{mq.pct(match_share.max())} of the evidence rows of a run (median "
        f"{mq.pct(match_share.median())}).\n\n"
        + mq.md_table(rows, ["Group", "Protein groups with LFQ (mean per run)", "Evidence from match between runs (mean)"])
    )
    short = f"{mq.n(any_lfq)} protein groups with LFQ, {mq.n(complete)} in all 40 runs"
    return mq.answer(QUESTION, text, figure, short)


if __name__ == "__main__":
    print(main())
