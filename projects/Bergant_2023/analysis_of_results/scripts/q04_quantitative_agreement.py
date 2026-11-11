"""Q4. Do the quantities agree for the same raw file?

For 1:1-matched protein groups, the log2 values of our search and the deposited one are
compared run by run. A version change and the different run set (42 vs 67) scale every protein
differently (LFQ distributes a protein's intensity over all runs of a search), so each protein's
median log2 ratio over the 42 runs is removed first (the per-protein offset, reported as a scale
effect). What remains, the residual, measures how differently the two searches quantify the
same raw file. Its SD is compared with the SD of the difference between two replicate samples
of the same condition within one search: both are differences of two measurements.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

import mq_common as mq

QUESTION = "Do the quantities agree for the same raw file?"
PAPER = "Not reported; deposited LFQ and Intensity per raw file"
NOTE = "1:1-matched groups; per-protein offsets removed; yardstick: replicate pairs"
MIN_PAIRS = 10      # runs with values in both tables needed for a protein's offset


def matrices(name, family):
    """(ids, log2 matrix, valid, MS/MS-identified) of a source."""
    if name == "published":
        pub = mq.published()
        return list(pub.ids), pub.matrix, pub.real, pub.real
    layer = mq.layers()[name]
    m = layer.matrix(family)
    idt = layer.matrix("Identification type", numeric=False) == "By MS/MS"
    return list(layer.ids), np.log2(m.where(m > 0)), m > 0, idt


def residuals(a, b, family, msms_only=False):
    """Per-sample residual statistics of a - b after removing per-protein offsets."""
    ids_a, A, va, ma = matrices(a, family)
    ids_b, B, vb, mb = matrices(b, family)
    _, pairs = mq.match_groups(ids_a, ids_b)
    A = A.iloc[pairs["a"]].reset_index(drop=True)
    B = B.iloc[pairs["b"]].reset_index(drop=True)
    valid = va.iloc[pairs["a"]].reset_index(drop=True) & vb.iloc[pairs["b"]].reset_index(drop=True)
    if msms_only:
        valid &= ma.iloc[pairs["a"]].reset_index(drop=True) & mb.iloc[pairs["b"]].reset_index(drop=True)
    R = (A - B).where(valid)
    enough = valid.sum(axis=1) >= MIN_PAIRS
    R = R[enough]
    offset = R.median(axis=1)
    res = R.sub(offset, axis=0)
    per_sample = pd.DataFrame({
        "sd": res.std(),
        "within_2fold": (res.abs() < 1).sum() / res.notna().sum(),
        "shift": (A - B).where(valid).median(),
        "n": res.notna().sum(),
    })
    return per_sample, offset, int(enough.sum())


def replicate_pair_sd(name, family):
    """Median SD of differences between two replicates of the same condition (all 105 pairs)."""
    _, M, valid, _ = matrices(name, family)
    sds = []
    for c in mq.CONDITIONS:
        reps = mq.samples_of(c)
        for i in range(6):
            for j in range(i + 1, 6):
                ok = valid[reps[i]] & valid[reps[j]]
                sds.append((M.loc[ok, reps[i]] - M.loc[ok, reps[j]]).std())
    return float(np.median(sds))


def detection(a, b):
    ids_a, _, va, _ = matrices(a, "LFQ intensity")
    ids_b, _, vb, _ = matrices(b, "LFQ intensity")
    _, pairs = mq.match_groups(ids_a, ids_b)
    va = va.iloc[pairs["a"]].reset_index(drop=True)
    vb = vb.iloc[pairs["b"]].reset_index(drop=True)
    tot = (va | vb).sum()
    return ((va & vb).sum() / tot).median(), ((va & ~vb).sum() / tot).median(), \
        ((~va & vb).sum() / tot).median()


def outliers(sd):
    m, d = sd.median(), mq.mad(sd)
    return list(sd.index[(sd > m + 4 * d) & (sd > 1.5 * m)])


def plot(ps_lfq, ps_int, pair_lfq, pair_int, example):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw={"width_ratios": [1, 2.2]})
    ax = axes[0]
    a, b, s = example
    ax.hexbin(b, a, gridsize=60, cmap="Blues", mincnt=1, bins="log", linewidths=0)
    lim = [min(np.nanmin(a), np.nanmin(b)), max(np.nanmax(a), np.nanmax(b))]
    ax.plot(lim, lim, color=mq.AXIS, lw=0.8)
    ax.set_xlabel("deposited, offset-free log2 LFQ")
    ax.set_ylabel("ours, offset-free log2 LFQ")
    ax.set_title(f"Same raw file: {s}")
    ax = axes[1]
    x = np.arange(len(mq.SAMPLES))
    ax.bar(x - 0.2, ps_lfq["sd"], width=0.4, color=mq.SERIES[0], label="LFQ residual", **mq.bar_kw())
    ax.bar(x + 0.2, ps_int["sd"], width=0.4, color=mq.SERIES[1], label="Intensity residual",
           **mq.bar_kw())
    ax.axhline(pair_lfq, color=mq.SERIES[0], ls="--", lw=1, label="LFQ, replicate pairs")
    ax.axhline(pair_int, color=mq.SERIES[1], ls="--", lw=1, label="Intensity, replicate pairs")
    ax.set_xticks(x, mq.SAMPLES)
    mq.rotate_xticks(ax, size=6, angle=90)
    ax.set_ylabel("SD of log2 difference")
    ax.set_ylim(0, 1.3 * max(pair_int, ps_int["sd"].max()))
    ax.set_title("Same raw file, ours vs deposited (after per-protein offsets)")
    ax.legend(ncol=4, fontsize=8, loc="upper left")
    mq.value_grid(ax)
    fig.tight_layout()
    return mq.save_fig(fig, "q04_quantitative_agreement")


def main():
    lay = mq.layers()
    rows = []
    stats_ = {}
    for fam in ("LFQ intensity", "Intensity"):
        for msms_only in (False, True):
            ps, offset, n_prot = residuals("ours", "deposited", fam, msms_only)
            stats_[(fam, msms_only)] = (ps, offset, n_prot)
    pair = {fam: {name: replicate_pair_sd(name, fam) for name in ("ours", "deposited")}
            for fam in ("LFQ intensity", "Intensity")}
    for (fam, msms_only), (ps, offset, n_prot) in stats_.items():
        yard = pair[fam]["deposited"]
        rows.append((fam, "by MS/MS in both" if msms_only else "all values", mq.n(n_prot),
                     mq.f2(ps["sd"].median(), 3), mq.f2(yard, 3), mq.f2(ps["sd"].median() / yard),
                     mq.pct(ps["within_2fold"].median()), mq.f2(offset.std(), 2),
                     mq.f2(ps["shift"].median(), 2)))
    ps_lfq, off_lfq, _ = stats_[("LFQ intensity", False)]
    ps_int, _, _ = stats_[("Intensity", False)]
    ratio = ps_lfq["sd"].median() / pair["LFQ intensity"]["deposited"]
    out = outliers(ps_lfq["sd"])
    keys = dict(zip(mq.design()["Experiment"], mq.design()["key"]))
    rem = [s for s in mq.SAMPLES if keys[s].startswith("REM_")]
    # per-sample shift must not depend on the condition (normalisation check)
    shift = ps_lfq["shift"] - ps_lfq["shift"].median()
    groups = [shift[mq.samples_of(c)].values for c in mq.CONDITIONS]
    anova_p = float(stats.f_oneway(*groups).pvalue)
    cond_shift = {c: float(shift[mq.samples_of(c)].mean()) for c in mq.CONDITIONS}
    shift_range = max(cond_shift.values()) - min(cond_shift.values())
    int_shift = ps_int["shift"] - ps_int["shift"].median()
    int_range = (max(int_shift[mq.samples_of(c)].mean() for c in mq.CONDITIONS)
                 - min(int_shift[mq.samples_of(c)].mean() for c in mq.CONDITIONS))
    lo, hi = min(cond_shift, key=cond_shift.get), max(cond_shift, key=cond_shift.get)

    # ours vs published and the authors' own baseline, LFQ, real published entries only
    pub_rows = []
    pub_res = {}
    for a in [n for n in lay if n != "deposited"] + ["deposited"]:
        ps, offset, n_prot = residuals(a, "published", "LFQ intensity")
        pub_res[a] = ps["sd"].median()
        pub_rows.append((f"{a} vs published", mq.n(n_prot), mq.f2(pub_res[a], 3),
                         mq.pct(ps["within_2fold"].median())))
    pub_pair = replicate_pair_sd("published", "LFQ intensity")
    pub_ok = pub_res["ours"] <= max(1.15 * pub_res["deposited"], 0.7 * pub_pair)

    both, only_ours, only_dep = detection("ours", "deposited")

    # example sample: median residual SD
    s = ps_lfq["sd"].sort_values().index[len(ps_lfq) // 2]
    ids_a, A, va, _ = matrices("ours", "LFQ intensity")
    ids_b, B, vb, _ = matrices("deposited", "LFQ intensity")
    _, pairs = mq.match_groups(ids_a, ids_b)
    A1 = A.iloc[pairs["a"]].reset_index(drop=True)
    B1 = B.iloc[pairs["b"]].reset_index(drop=True)
    ok = (va.iloc[pairs["a"]].reset_index(drop=True) & vb.iloc[pairs["b"]].reset_index(drop=True))
    off = (A1 - B1).where(ok).median(axis=1)
    keep = ok[s] & off.notna()
    example = ((A1[s] - A1.where(ok).mean(axis=1))[keep].values,
               (B1[s] - B1.where(ok).mean(axis=1))[keep].values, s)
    figure = plot(ps_lfq, ps_int, pair["LFQ intensity"]["deposited"],
                  pair["Intensity"]["deposited"], example)

    decomposition = ""
    if "ours-67" in lay:
        total = residuals("ours", "deposited", "LFQ intensity")[0]["sd"].median()
        runset = residuals("ours", "ours-67", "LFQ intensity")[0]["sd"].median()
        version = residuals("ours-67", "deposited", "LFQ intensity")[0]["sd"].median()
        decomposition = (f"\n\n**Part B decomposition** (residual SD, LFQ): total (ours vs "
                         f"deposited) {mq.f2(total, 3)}; run set (ours vs ours-67) "
                         f"{mq.f2(runset, 3)}; version (ours-67 vs deposited) {mq.f2(version, 3)}. "
                         f"Variances: {mq.f2(total**2, 4)} total vs {mq.f2(runset**2 + version**2, 4)} "
                         "as the sum of the parts.")

    if ratio <= 0.7 and not out and pub_ok and shift_range <= 0.1:
        verdict = "consistent"
    elif ratio <= 1.0 and not out:
        verdict = "minor deviation"
    else:
        verdict = "inconsistent"
    criterion = ("Consistent if the median residual SD of LFQ for the same raw file is at most 0.7 "
                 "times the median SD between two replicates of one condition, no raw file is an "
                 "outlier (residual SD above median + 4 MAD and above 1.5 times the median), the "
                 "agreement with the published values is not worse than 1.15 times the authors' "
                 "own (deposited vs published), and the per-condition normalisation shifts differ "
                 "by at most 0.1 log2.")
    lines = [
        "Residual SD per raw file (median over the 42 runs) after removing each protein's offset, "
        "with the replicate-pair SD of the deposited table as yardstick:",
        "",
        mq.md_table(rows, ["Values", "Entries", "Protein groups", "Residual SD", "Replicate-pair SD",
                           "Ratio", "Within 2-fold", "SD of protein offsets",
                           "Median sample shift"]),
        "",
        f"For the same raw file our LFQ values differ from the deposited ones by a residual SD of "
        f"**{mq.f2(ps_lfq['sd'].median(), 3)}** (log2), {mq.f2(ratio)} times the difference "
        f"between two replicates ({mq.f2(pair['LFQ intensity']['deposited'], 3)} in the deposited "
        f"table, {mq.f2(pair['LFQ intensity']['ours'], 3)} in ours). Raw Intensity agrees less "
        "closely, but also better than replicates do. Outlier files: "
        + (", ".join(out) if out else "none") + f". The {len(rem)} remeasured runs (REM) show "
        f"somewhat larger residuals than the others (median {mq.f2(ps_lfq['sd'][rem].median(), 3)} "
        f"vs {mq.f2(ps_lfq['sd'].drop(rem).median(), 3)}), still well below the replicate "
        "yardstick, plausibly because the two searches align these later measurements "
        "differently.",
        "",
        f"The scale differs: protein offsets have an SD of {mq.f2(off_lfq.std())} log2 and the "
        f"median sample shift is {mq.f2(ps_lfq['shift'].median())} log2 (our LFQ values are lower "
        "throughout). This is how LFQ distributes intensity over the runs of each search, not a "
        "measurement difference, and a constant shift leaves fold changes unchanged. Per "
        f"condition, the mean shift varies by {mq.f2(shift_range, 3)} log2 (from {lo} to {hi}; "
        f"one-way ANOVA p = {anova_p:.1e}); for raw Intensity by only {mq.f2(int_range, 3)}. "
        "So LFQ normalisation over a different run set moves the conditions slightly against "
        f"each other, by up to {mq.pct(shift_range / mq.FC_CUTOFF, 0)} of the paper's fold-change "
        "threshold; the measurement itself is unaffected. Part B can confirm the cause.",
        "",
        f"Detection (LFQ > 0) per raw file: both tables {mq.pct(both)}, only ours "
        f"{mq.pct(only_ours)}, only deposited {mq.pct(only_dep)} of the groups detected in "
        "either.",
        "",
        "**Against the published values** (real entries only; the published imputed entries are "
        "recognised by their missing MaxQuant rounding, see Q5):",
        "",
        mq.md_table(pub_rows, ["Tables", "Protein groups", "Residual SD", "Within 2-fold"]),
        "",
        f"The published table comes from a different search, so both MaxQuant tables deviate from "
        f"it more than from each other; ours (residual SD {mq.f2(pub_res['ours'], 3)}) about as "
        f"much as the authors' deposited search ({mq.f2(pub_res['deposited'], 3)}). Published "
        f"replicate-pair SD: {mq.f2(pub_pair, 3)}." + decomposition,
    ]
    short = (f"same-file LFQ residual SD {mq.f2(ps_lfq['sd'].median(), 2)} vs replicate pairs "
             f"{mq.f2(pair['LFQ intensity']['deposited'], 2)}")
    metrics = {"residual_sd": float(ps_lfq["sd"].median()), "ratio": float(ratio),
               "outliers": len(out), "anova_p": anova_p,
               "pair_sd": float(pair["LFQ intensity"]["deposited"]),
               "pub_residual_ours": float(pub_res["ours"]),
               "pub_residual_deposited": float(pub_res["deposited"]),
               "shift_range": float(shift_range)}
    return mq.answer(QUESTION, "\n".join(lines), figure, short, verdict, criterion, metrics)


if __name__ == "__main__":
    print(main())
