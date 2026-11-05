"""Q10. How much did the serum albumin abundance percentage decrease?
Was the depletion overall successful?

Clean set. Albumin share = albumin intensity / sample intensity * 100.
The other top serum proteins are compared as % of the non-albumin
intensity, because removing albumin alone raises every other share.

Criteria for the verdict:
- albumin clearly reduced: share at least 2-fold lower in the depleted sample
- other targets reduced: lower share of the non-albumin intensity
- depth gained: more protein groups identified in the depleted sample
"""
import mq_common as mq

QUESTION = ("How much did the serum albumin abundance percentage decrease? "
            "Was the depletion overall successful?")
TUTORIAL = ("Albumin share down 58 percentage points; quite successful (about twice as many "
            "proteins), some targets unchanged")
NOTE = "tutorial uses LFQ shares, here Intensity shares (see below)"


def main():
    pg = mq.load_protein_groups()
    cl = mq.clean(pg)
    s = mq.samples(cl)
    colors = mq.sample_colors(s)
    dep, und = mq.depleted_sample(pg)
    alb, top = mq.is_albumin(cl), mq.is_top_serum(cl)

    alb_sh = {x: mq.share(cl, x)[alb].iloc[0] for x in s}
    other_top = {}
    for x in s:
        inten = mq.intensity(cl, x)
        other_top[x] = inten[top & ~alb].sum() / inten[~alb].sum() * 100
    ids = {x: int(mq.identified_in(cl, x).sum()) for x in s}
    rank = {x: int((mq.share(cl, x) > alb_sh[x]).sum()) + 1 for x in s}
    # per target: share of the non-albumin intensity, depleted vs undepleted
    nonalb = {x: mq.intensity(cl, x)[~alb].sum() for x in s}
    target_sh = {x: mq.intensity(cl, x)[top & ~alb] / nonalb[x] * 100 for x in s}
    gene = mq.genes(cl).str[0]
    not_reduced = sorted(gene[i] for i in target_sh[dep].index if target_sh[dep][i] >= target_sh[und][i])
    n_targets = int((top & ~alb).sum())

    drop_pp = alb_sh[und] - alb_sh[dep]
    drop_rel = drop_pp / alb_sh[und] * 100
    fold = alb_sh[und] / alb_sh[dep]

    fig, axes = mq.plt.subplots(1, 3, figsize=(10, 3.4))
    panels = [
        ("Albumin", "% of sample intensity", alb_sh, mq.pct),
        ("Other top serum proteins", "% of non-albumin intensity", other_top, mq.pct),
        ("Protein groups identified", "Protein groups", ids, mq.n),
    ]
    for ax, (title, ylabel, vals, fmt) in zip(axes, panels):
        b = ax.bar(s, [vals[x] for x in s], color=[colors[x] for x in s], width=0.6, **mq.bar_kw())
        mq.label_bars(ax, b, [fmt(vals[x]) for x in s])
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xticks(range(len(s)), [f"{x}\n(depleted)" if x == dep else x for x in s])
        ax.tick_params(axis="x", length=0)
        mq.value_grid(ax, "y")
    fig.tight_layout()
    figure = mq.save_fig(fig, "q10_albumin_depletion")

    alb_ok = fold >= 2
    others_ok = other_top[dep] < other_top[und]
    depth_ok = ids[dep] > ids[und]
    exceptions = (f" Not reduced: {len(not_reduced)} of the {n_targets} other target groups "
                  f"({', '.join(not_reduced)}) have the same or a higher share of the non-albumin "
                  f"intensity in {dep}."
                  if not_reduced else " Every other target group was reduced.")
    verdict_short = ("overall successful" if alb_ok and others_ok else
                     "only partly successful" if alb_ok else "not successful")
    if alb_ok and others_ok:
        verdict = ("**Yes, overall successful**: albumin and, taken together, the other top serum "
                   f"proteins were reduced ({mq.pct(other_top[dep])} vs {mq.pct(other_top[und])} of the "
                   f"non-albumin intensity).{exceptions}")
    elif alb_ok:
        verdict = ("**Only partly successful**: albumin was clearly reduced, but the other top serum "
                   f"proteins were not. They make up {mq.pct(other_top[dep])} of the non-albumin "
                   f"intensity in {dep} against {mq.pct(other_top[und])} in {und}.")
    else:
        verdict = "**Not successful**: albumin was not clearly reduced."
    depth = (f"The depleted sample yielded more protein identifications ({mq.n(ids[dep])} vs "
             f"{mq.n(ids[und])}), so depletion did gain depth."
             if depth_ok else
             f"The depleted sample did not yield more protein identifications ({mq.n(ids[dep])} vs "
             f"{mq.n(ids[und])}).")
    rows = [(x + (" (depleted)" if x == dep else ""), mq.pct(alb_sh[x]), rank[x], mq.pct(other_top[x]), mq.n(ids[x]))
            for x in s]
    text = (
        f"The albumin share fell from **{mq.pct(alb_sh[und])}** in {und} to **{mq.pct(alb_sh[dep])}** "
        f"in {dep}: **{drop_pp:.1f} percentage points**, a relative decrease of {drop_rel:.0f}% "
        f"({fold:.1f}-fold).\n\n"
        + mq.md_table(rows, ["Sample", "Albumin share", "Albumin rank",
                             "Other top serum, % of non-albumin", "Protein groups"])
        + f"\n\n{verdict} {depth}"
    )
    short = (f"Albumin share down {drop_pp:.1f} percentage points; {verdict_short} "
             f"({mq.n(ids[dep])} vs {mq.n(ids[und])} protein groups), "
             f"{len(not_reduced)} other targets not reduced")
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
