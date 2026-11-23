"""Q2. Are the samples assigned to the right pigs and groups?

The deposited raw files carry run numbers only (..._MS<n>.raw). make_parameters_Laybourn_2024.sh
assigned them: upper trachea MS1-14 and MS19-26 with run number = pig (the authors' Spectronaut
experiment labels MS1 as an upper-trachea control), lower trachea MS43-56 and MS61-68 with run
number = pig + 42 (inferred from the numbering), and the group of each pig from GEO GSE268254.
This question tests the assignment with the data: viral peptides (only swH1N1 tissue should
contain detectable virus; the paper quantified no viral protein after huH1N1), interferon-
stimulated proteins (higher in infected pigs), and whether the two tissues of a pig resemble
each other.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

import lb_common as lb

QUESTION = "Are the samples assigned to the right pigs and groups?"
PAPER = "22 pigs: 6 control, 8 swH1N1, 8 huH1N1; upper and lower trachea"
NOTE = "run -> pig from the run numbers, pig -> group from GEO; tested with viral peptides and ISGs"
MIN_VIRAL = 2           # passing viral precursors for a virus-positive sample
FIG1C_LOWEST = {"UT": [9, 10], "LT": [9, 10]}   # swH1N1 pigs labelled lowest in Fig. 1C


def pairing_test(s):
    """Correlate per-pig residuals (sample minus its group mean, log2 LFQ, proteins with values
    in all 44 runs) between the two tissues: is the same pig the best match?"""
    m = np.log2(s.lfq[s.lfq.notna().all(axis=1)])
    res = m.copy()
    for t in lb.TISSUES:
        for g in lb.GROUPS:
            c = lb.samples_of(t, g)
            res[c] = m[c].sub(m[c].mean(axis=1), axis=0)
    ut, lt = lb.samples_of("UT"), lb.samples_of("LT")
    c = np.corrcoef(res[ut].T.values, res[lt].T.values)[:len(ut), len(ut):]
    best = [lb.pig_of(ut[int(np.argmax(c[:, j]))]) == lb.pig_of(lt[j]) for j in range(len(lt))]
    hits = int(sum(best))
    p = stats.binomtest(hits, len(lt), 1 / len(ut), alternative="greater").pvalue
    return hits, len(lt), p, len(m)


Y_GENES = {"DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y", "ZFY", "EIF2S3Y"}


def y_check(s):
    """Peptides of Y-chromosome proteins only: (genes found, peptides, runs with a passing one)."""
    genes = lb.pig_genes()
    ev = s.evidence
    g = ev["Proteins"].fillna("").str.split(";").map(lambda a: {genes.get(x, "?") for x in a})
    y = ev[g.map(lambda x: bool(x) and x <= Y_GENES)]
    found = sorted(set().union(*g[y.index])) if len(y) else []
    return found, y.Sequence.nunique(), y[y.passes].Experiment.nunique()


def main():
    s = lb.search("ours")
    d = lb.design()
    geo = lb.geo_groups()
    vp = lb.viral_precursors("ours")
    isg = lb.isg_score("ours")

    # assignment table
    assign_rows = []
    for g in lb.GROUPS:
        pigs = sorted(p for p, x in geo.items() if x == g)
        ut = d[(d.tissue == "UT") & (d.group == g)].run
        lt = d[(d.tissue == "LT") & (d.group == g)].run
        assign_rows.append((lb.GROUP_LABELS[g], f"{min(pigs)}-{max(pigs)} ({len(pigs)})",
                            f"MS{ut.min()}-{ut.max()}", f"MS{lt.min()}-{lt.max()}"))
    design_ok = all(geo[p] == g for p, g in zip(d.pig, d.group))

    # viral calls
    pos = vp["total"] >= MIN_VIRAL
    viral_rows, wrong_virus = [], []
    for t in lb.TISSUES:
        for g in lb.GROUPS:
            c = lb.samples_of(t, g)
            x = vp.loc[c]
            viral_rows.append((f"{t} {g}", f"{int(pos[c].sum())}/{len(c)}",
                               f"{int(x.total.min())}-{int(x.total.max())}",
                               int(x.sw_only.sum()), int(x.hu_only.sum()), int(x.shared.sum())))
            for smp in c:
                other = {"sw": "hu_only", "hu": "sw_only"}.get(g)
                if other and vp.loc[smp, other] > 0:
                    wrong_virus.append(smp)
    sw = [x for x in lb.samples() if lb.group_of(x) == "sw"]
    ctr = [x for x in lb.samples() if lb.group_of(x) == "ctr"]
    hu = [x for x in lb.samples() if lb.group_of(x) == "hu"]
    sw_neg = [x for x in sw if not pos[x]]
    ctr_pos = [x for x in ctr if pos[x]]
    ctr_any = [x for x in ctr if vp.loc[x, "total"] > 0]
    hu_pos = [x for x in hu if pos[x]]

    # ISG separation
    isg_rows, overlaps = [], []
    for t in lb.TISSUES:
        cmax = isg[lb.samples_of(t, "ctr")].max()
        for g in ("sw", "hu"):
            c = lb.samples_of(t, g)
            low = [x for x in c if isg[x] <= cmax]
            overlaps += low
            isg_rows.append((f"{t} {g}", lb.f2(isg[c].min()), lb.f2(isg[c].median()),
                             lb.f2(cmax), f"{len(c) - len(low)}/{len(c)}"))

    # MaxDIA transfers: viral groups as reported vs run-filtered
    transfer_rows = []
    lfq, det = s.lfq, s.detected
    for gid, members in sorted(lb.search("ours").viral_groups.items()):
        lab = s.labels[gid]
        cells = []
        for g in lb.GROUPS:
            c = [x for x in lb.samples() if lb.group_of(x) == g]
            cells.append(f"{int(lfq.loc[gid, c].notna().sum())} / {int(det.loc[gid, c].sum())}")
        transfer_rows.append((lab, *cells))

    # viral load order vs Fig. 1C
    order_rows = []
    for t in lb.TISSUES:
        c = lb.samples_of(t, "sw")
        ranked = vp.loc[c, "log2_intensity"].sort_values()
        lowest = [lb.pig_of(x) for x in ranked.index[:2]]
        order_rows.append((t, ", ".join(f"pig {p}" for p in FIG1C_LOWEST[t]),
                           ", ".join(f"pig {p}" for p in lowest),
                           ", ".join(f"{lb.pig_of(x)}" for x in ranked.index)))

    hits, n_lt, p_pair, n_prot = pairing_test(s)
    y_genes, y_peptides, y_runs = y_check(s)
    hu_only_total = int(vp["hu_only"].sum())
    sw_in_hu = int(vp.loc[hu, "sw_only"].sum())
    np_gid = next(gid for gid in s.viral_groups if s.labels[gid].startswith("NP"))
    np_ctr_det = int(det.loc[np_gid, ctr].sum())
    np_gap = float(np.log2(lfq.loc[np_gid, sw]).median() - np.log2(lfq.loc[np_gid, ctr]).median())

    # figure: viral precursors and ISG score per sample
    fig, axes = plt.subplots(2, 1, figsize=(10, 4.8), sharex=True)
    x = np.arange(len(lb.samples()))
    colors = [lb.GROUP_COLORS[lb.group_of(smp)] for smp in lb.samples()]
    axes[0].bar(x, vp.loc[lb.samples(), "total"], color=colors, width=0.8)
    axes[0].axhline(MIN_VIRAL - 0.5, color=lb.INK2, lw=0.8, ls="--")
    axes[0].set_ylabel("viral precursors\n(run-level q <= 0.01)")
    axes[1].bar(x, isg[lb.samples()], color=colors, width=0.8)
    axes[1].axhline(0, color=lb.AXIS, lw=0.6)
    axes[1].set_ylabel("ISG score\n(mean z, 6 proteins)")
    for ax in axes:
        lb.value_grid(ax)
        ax.axvline(21.5, color=lb.AXIS, lw=0.8)
    axes[1].set_xticks(x, [f"{smp[:2]} {lb.pig_of(smp)}" for smp in lb.samples()], rotation=90,
                       fontsize=7)
    axes[0].set_title("Viral peptides and interferon-stimulated proteins per sample")
    handles = [plt.Rectangle((0, 0), 1, 1, color=lb.GROUP_COLORS[g]) for g in lb.GROUPS]
    axes[0].legend(handles, [lb.GROUP_LABELS[g] for g in lb.GROUPS], loc="upper right", ncol=3,
                   fontsize=7)
    figure = lb.save_fig(fig, "q02_viral_isg_per_sample")

    virus_fail = len(sw_neg) + len(ctr_pos) + len(wrong_virus)
    if design_ok and virus_fail == 0 and len(overlaps) <= 2:
        verdict = "consistent"
    elif design_ok and virus_fail <= 2 and len(overlaps) <= 6:
        verdict = "minor deviation"
    else:
        verdict = "inconsistent"
    criterion = (f"Consistent if every swH1N1 sample is virus-positive (>= {MIN_VIRAL} viral "
                 "precursors passing run-level q <= 0.01 and PEP <= 0.1), every control is "
                 "virus-negative, no infected sample has a strain-specific peptide of the other "
                 "strain, and at most 2 of the 32 infected samples have an ISG score at or below "
                 "the highest control of their tissue; minor deviation if at most 2 virus calls "
                 "and at most 6 ISG scores disagree; inconsistent otherwise.")
    lines = [
        "**Assignment** (make_parameters_Laybourn_2024.sh): run number = pig in the upper trachea, "
        "pig + 42 in the lower trachea (inferred), group from GEO GSE268254"
        + (" (all 44 samples agree with GEO)." if design_ok else " (DISAGREES with GEO)."),
        "",
        lb.md_table(assign_rows, ["Group", "Pigs (n)", "Upper trachea", "Lower trachea"]),
        "",
        f"**Viral peptides.** A sample is virus-positive with >= {MIN_VIRAL} viral precursors "
        "passing in that run (sw only: in swH1N1 proteins only; shared: identical in both "
        "strains):",
        "",
        lb.md_table(viral_rows, ["Samples", "Virus-positive", "Viral precursors (range)",
                                 "sw only", "hu only", "shared"]),
        "",
        f"All {len(sw) - len(sw_neg)} of {len(sw)} swH1N1 samples are virus-positive; "
        f"{len(ctr_pos)} of {len(ctr)} controls are"
        + (f" ({len(ctr_any)} controls have a single passing viral precursor: "
           f"{', '.join(ctr_any)}; one false identification per run is what a 1% run-level FDR "
           "allows)" if ctr_any else "")
        + f". {len(hu_pos)} of {len(hu)} huH1N1 samples are virus-positive, in line with the paper "
        "(no viral protein quantified after huH1N1, and a 7-fold (upper) to 37-fold (lower "
        "trachea) lower median viral RNA load than after swH1N1, Fig. 1C). "
        + (f"Peptides specific to the huH1N1 strain pass {hu_only_total} times in all runs "
           "together" if hu_only_total else "No peptide specific to the huH1N1 strain passes in "
           "any run")
        + (f", and swH1N1-specific peptides pass {sw_in_hu} times in huH1N1 samples."
           if sw_in_hu else ", and no swH1N1-specific peptide passes in a huH1N1 sample."),
        "",
        "**Interferon-stimulated proteins** (MX1, MX2, ISG15, OASL, IFIT3, IFI6, the proteins of "
        "Fig. 2; mean z-score of log2 LFQ within each tissue):",
        "",
        lb.md_table(isg_rows, ["Samples", "Lowest", "Median", "Highest control",
                               "Above all controls"]),
        "",
        f"{32 - len(overlaps)} of 32 infected samples score above every control of their tissue"
        + (f"; not: {', '.join(overlaps)}." if overlaps else "."),
        "",
        "**Why run-level detection matters.** Samples with a value as reported / with the "
        "protein group detected in that run (viral protein groups):",
        "",
        lb.md_table(transfer_rows, ["Viral protein group", "controls (12)", "swH1N1 (16)",
                                    "huH1N1 (16)"]),
        "",
        f"MaxDIA reports the viral nucleoprotein in {int(lfq.loc[np_gid, ctr].notna().sum())} of "
        f"12 controls, with a median only {np_gap:.1f} log2 units below the swH1N1 samples, "
        f"although its precursors pass in {np_ctr_det} control run(s): the reported value is a "
        "transferred signal without an identification. The run-filtered layer removes such "
        "values (Q1).",
        "",
        "**Viral load order** (swH1N1; our summed intensity of passing viral precursors vs the "
        "pigs that Fig. 1C labels as lowest by RT-qPCR):",
        "",
        lb.md_table(order_rows, ["Tissue", "Lowest in Fig. 1C", "Lowest in ours",
                                 "Ours, low to high (pig)"]),
        "",
        "Viral protein and viral RNA are different measures from different tissue pieces, so the "
        "order is only a weak check and not part of the verdict.",
        "",
        f"**Same pig in both tissues.** Correlating per-pig residuals (each sample minus its group "
        f"mean, {n_prot:,} proteins with values in all runs) between the tissues, the same pig is "
        f"the best match for {hits} of {n_lt} lower-trachea samples (1 expected by chance; "
        f"binomial p = {p_pair:.3f}). Y-chromosome proteins cannot separate the pigs either: "
        f"{y_peptides} peptide(s) specific to them ({', '.join(y_genes) or 'none'}), passing in "
        f"{y_runs} run(s). The data support "
        "the group of every sample but cannot confirm which pig of a group a lower-trachea "
        "sample comes from; the inferred numbering is kept.",
    ]
    short = (f"virus-positive: swH1N1 {len(sw) - len(sw_neg)}/{len(sw)}, controls "
             f"{len(ctr_pos)}/{len(ctr)}, huH1N1 {len(hu_pos)}/{len(hu)}; ISG score above all "
             f"controls in {32 - len(overlaps)}/32 infected")
    metrics = {"sw_positive": len(sw) - len(sw_neg), "ctr_positive": len(ctr_pos),
               "hu_positive": len(hu_pos), "isg_overlaps": len(overlaps),
               "pairing_hits": hits, "pairing_p": float(p_pair)}
    return lb.answer(QUESTION, "\n".join(lines), figure, short, verdict, criterion, metrics)


if __name__ == "__main__":
    print(main())
