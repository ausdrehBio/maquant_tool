"""Q5. Do our differential results agree with Table S3?

Applies the reconstructed paper pipeline (Q4) to both layers of our search and follows every
published DEP into our results: not quantified (protein not identified, or no value in that
tissue), opposite direction, same direction but not significant, or a DEP in the same
direction. A calibration runs the same pipeline on our data with randomly permuted group labels:
how many DEPs does it report when the groups carry no information?
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

import lb_common as lb

QUESTION = "Do our differential results agree with Table S3?"
PAPER = "623 DEPs in six comparisons (Table S3; Fig. 1A)"
NOTE = "paper pipeline reconstructed (Q4); run-filtered layer approximates Spectronaut's filter"
N_PERM = 100
PERM_SEED = 20240816    # publication date of the paper


def permutation_counts(layer):
    """DEPs per comparison for random group labels (labels permuted within each tissue)."""
    m = lb.search("ours").layer(layer)
    rng = np.random.default_rng(PERM_SEED)
    out = {n: [] for n in lb.COMPARISON_NAMES}
    for _ in range(N_PERM):
        perm = m.copy()
        for t in lb.TISSUES:
            cols = lb.samples_of(t)
            perm[cols] = m[list(rng.permutation(cols))].values
            st = lb.paper_stats(perm, t)
            for n in lb.COMPARISON_NAMES:
                if n.startswith(t):
                    out[n].append(int((st[f"{n}|call"] != 0).sum()))
    return {n: np.array(v) for n, v in out.items()}


def match_type(s, acc, gid):
    if gid is None or pd.isna(gid):
        return "not identified"
    maj = s.pg.loc[int(gid), "Majority protein IDs"].split(";")
    if lb.accession_token(maj[0]) == acc or maj[0] == acc:
        return "leads our group"
    return "member of a larger group"


def main():
    s = lb.search("ours")
    pub = lb.table_s3()
    summ = {layer: lb.agreement_summary("ours", layer) for layer in lb.LAYERS}
    counts = {layer: lb.our_dep_counts("ours", layer) for layer in lb.LAYERS}
    prim = lb.PRIMARY_LAYER
    a = lb.agreement("ours", prim).copy()
    a["match"] = [match_type(s, r.acc, r.group_id) for r in a.itertuples()]
    lb.save_table(a, "q05_published_deps_ours_run-filtered", index=False)
    lb.save_table(lb.agreement("ours", "as reported"), "q05_published_deps_ours_as-reported",
                  index=False)

    # per comparison table
    comp_rows = []
    for name in lb.COMPARISON_NAMES:
        row = [name, int((pub.comparison == name).sum())]
        for layer in lb.LAYERS:
            ag = lb.agreement("ours", layer)
            x = ag[ag.comparison == name]
            q = x[x.outcome != "not quantified"]
            same = q.outcome.isin(["DEP, same direction", "same direction, not DEP"]).mean()
            c = counts[layer][name]
            row += [c["up"] + c["down"], lb.pct(same, 0),
                    f"{int((x.outcome == 'DEP, same direction').sum())} "
                    f"({lb.pct((x.outcome == 'DEP, same direction').mean(), 0)})"]
        comp_rows.append(row)

    # outcomes pooled, by class and by match type (primary layer)
    def outcome_counts(df):
        return [int((df.outcome == o).sum()) for o in lb.OUTCOME_COLORS]
    class_rows = []
    for label, sel in [("all", a.index == a.index),
                       (f"quantitative (|log2 FC| < {lb.PRESENCE_LFC:.0f})", ~a.presence),
                       (f"presence/absence (|log2 FC| >= {lb.PRESENCE_LFC:.0f})", a.presence),
                       ("accession leads our group", a.match == "leads our group"),
                       ("accession in a larger group", a.match == "member of a larger group"),
                       ("adjusted p < 0.001", a.published_padj < 0.001),
                       ("adjusted p >= 0.001", a.published_padj >= 0.001)]:
        x = a[sel]
        q = x[x.outcome != "not quantified"]
        same = q.outcome.isin(["DEP, same direction", "same direction, not DEP"]).mean()
        class_rows.append((label, len(x), *outcome_counts(x), lb.pct(same, 0),
                           lb.pct((x.outcome == "DEP, same direction").mean(), 0)))

    # numbers for the text: (n, same direction, found as DEP) per class
    def cls_stats(sel):
        x = a[sel]
        q = x[x.outcome != "not quantified"]
        same = q.outcome.isin(["DEP, same direction", "same direction, not DEP"]).mean()
        return (len(x), lb.pct(same, 0), lb.pct((x.outcome == "DEP, same direction").mean(), 0))
    cls = {"strong": cls_stats(a.published_padj < 0.001),
           "weak": cls_stats(a.published_padj >= 0.001),
           "quant": cls_stats(~a.presence), "presence": cls_stats(a.presence),
           "member": cls_stats(a.match == "member of a larger group")}
    by_comp = {}
    for name in lb.COMPARISON_NAMES:
        x = a[(a.comparison == name) & (a.outcome != "not quantified")]
        by_comp[name] = x.outcome.isin(["DEP, same direction", "same direction, not DEP"]).mean()
    worst_comp = min(by_comp, key=by_comp.get)
    worst_dir = by_comp[worst_comp]
    best_comp = max(by_comp, key=by_comp.get)

    # fold-change correlation for quantitative published DEPs
    qd = a[(~a.presence) & (a.outcome != "not quantified")]
    r_fc = stats.pearsonr(qd.published_lfc, qd.our_lfc)[0] if len(qd) > 2 else np.nan

    # calibration with random labels
    perm = permutation_counts(prim)
    perm_rows = []
    for name in lb.COMPARISON_NAMES:
        v = perm[name]
        c = counts[prim][name]
        perm_rows.append((name, int((pub.comparison == name).sum()), c["up"] + c["down"],
                          int(np.median(v)),
                          f"{int(np.percentile(v, 5))}-{int(np.percentile(v, 95))}"))
    perm_share = np.median([np.median(perm[n]) / (counts[prim][n]["up"] + counts[prim][n]["down"])
                            for n in lb.COMPARISON_NAMES])

    # figures
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, layer in zip(axes, lb.LAYERS):
        ag = lb.agreement("ours", layer)
        left = np.zeros(len(lb.COMPARISON_NAMES))
        for o, color in lb.OUTCOME_COLORS.items():
            v = np.array([(ag[ag.comparison == n].outcome == o).sum() for n in lb.COMPARISON_NAMES])
            ax.barh(lb.COMPARISON_NAMES, v, left=left, color=color, label=o)
            left += v
        ax.set_title(f"ours, {layer}")
        ax.set_xlabel("published DEPs")
        ax.invert_yaxis() if ax is axes[0] else None
        lb.value_grid(ax, "x")
    axes[1].legend(loc="lower right", fontsize=7)
    fig.suptitle("Where the published DEPs end up in our data (paper pipeline)", fontsize=10,
                 fontweight="bold")
    fig_outcomes = lb.save_fig(fig, "q05_published_dep_outcomes")

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    q_all = a[a.outcome != "not quantified"]
    ax.scatter(q_all.published_lfc, q_all.our_lfc, s=9, alpha=0.6,
               c=[lb.OUTCOME_COLORS[o] for o in q_all.outcome], edgecolor="none")
    ax.axhline(0, color=lb.AXIS, lw=0.6)
    ax.axvline(0, color=lb.AXIS, lw=0.6)
    for v in (-lb.PRESENCE_LFC, lb.PRESENCE_LFC):
        ax.axvline(v, color=lb.AXIS, lw=0.6, ls=":")
    ax.set_xlabel("published log2 FC (Table S3)")
    ax.set_ylabel(f"our log2 FC ({prim}, missing = 0)")
    ax.set_title("Published DEPs: fold changes")
    lb.value_grid(ax, "both")
    fig_scatter = lb.save_fig(fig, "q05_fold_changes")

    direction, recall = summ[prim]["direction"], summ[prim]["recall"]
    v_dir = lb.grade(direction, 0.8, 0.65)
    v_rec = lb.grade(recall, 0.4, 0.2)
    verdict = lb.worst(v_dir, v_rec)
    criterion = (f"On the {prim} layer, pooled over the six comparisons: consistent if at least "
                 "80% of the published DEPs that we quantify change in the same direction and at "
                 "least 40% are DEPs in ours in the same direction; minor deviation at 65% and "
                 "20%; inconsistent below. The recall thresholds are low on purpose: the "
                 "published tests are not corrected across proteins, so part of every list is "
                 "expected to be false positives that no reanalysis reproduces.")
    hdr = ["Comparison", "Published DEPs"]
    for layer in lb.LAYERS:
        hdr += [f"Our DEPs ({layer})", "Same direction", "Published DEPs found"]
    lines = [
        lb.md_table(comp_rows, hdr),
        "",
        f"Pooled ({prim}): of {summ[prim]['published']} published DEPs we quantify "
        f"{summ[prim]['quantified']}; {lb.pct(direction)} of these change in the same direction "
        f"and {lb.pct(recall)} of all published DEPs are DEPs in ours (as reported: "
        f"{lb.pct(summ['as reported']['direction'])} and {lb.pct(summ['as reported']['recall'])}). "
        "The two layers agree about equally well with the paper.",
        "",
        f"**By kind of published DEP** ({prim}; outcome columns as in the figure):",
        "",
        lb.md_table(class_rows, ["Published DEPs", "n", *lb.OUTCOME_COLORS, "Same direction",
                                 "Found as DEP"]),
        "",
        f"- **Strong published DEPs agree; weak ones do not.** Of the {cls['strong'][0]} DEPs "
        f"with a published adjusted p < 0.001, {cls['strong'][1]} change in the same direction "
        f"in ours and {cls['strong'][2]} are DEPs; of the other {cls['weak'][0]}, "
        f"{cls['weak'][1]} and {cls['weak'][2]}. The strong ones include the "
        "interferon-stimulated proteins that carry the paper's biology (Q6).",
        f"- **Quantitative vs presence/absence.** Published DEPs with |log2 FC| < "
        f"{lb.PRESENCE_LFC:.0f} are found as DEPs in {cls['quant'][2]} of cases, those with "
        f"larger values (a protein with values in one group and none in another in the "
        f"authors' data; {int(a.presence.sum())} of {len(a)}) in {cls['presence'][2]}. The "
        "presence/absence entries are often proteins we do not identify (Q3) or see in all "
        "groups alike. For the quantitative ones, the fold changes correlate only weakly with "
        f"ours (r = {lb.f2(r_fc)}, n = {len(qd)}), and {cls['quant'][1]} go in the same "
        "direction.",
        f"- **Protein grouping.** Where the published accession is only a member of a larger "
        f"MaxQuant group ({cls['member'][0]} DEPs), our value belongs to the larger group and "
        f"{cls['member'][2]} are found. ANXA5 (published log2 FC +10 in the lower trachea), "
        "STAT1 (-12, sw vs hu) and DFFA (-14 and -12 in the upper trachea) are such cases: "
        "abundant proteins with values in every sample of ours, which the authors' table lists "
        "as present in one group only, apparently an isoform-level group of Spectronaut.",
        f"- **By comparison**, agreement is lowest for {worst_comp} ({lb.pct(worst_dir, 0)} "
        f"in the same direction, close to the 50% of chance) and highest for {best_comp} "
        f"({lb.pct(by_comp[best_comp], 0)}).",
        "",
        f"**Calibration: random group labels** ({N_PERM} permutations of the labels within each "
        f"tissue, {prim} layer, seed {PERM_SEED}):",
        "",
        lb.md_table(perm_rows, ["Comparison", "Published DEPs", "Our DEPs",
                                "Random labels: median", "5-95%"]),
        "",
        "With labels that carry no information, the paper's pipeline still reports dozens of "
        "DEPs per comparison, and the upper range reaches the size of the published lists "
        f"(the median is {lb.pct(perm_share, 0)} of our real count). The pipeline adjusts only "
        "across the three pairwise tests of a protein, and setting missing values to 0 turns "
        "sporadic detections into large fold changes. A substantial part of each published list "
        "is therefore expected to be noise, which limits the recall any reanalysis can reach. "
        "The same holds for our own DEP lists.",
        "",
        f"Per-protein outcomes: tables/q05_published_deps_ours_run-filtered.tsv (not in git).",
    ]
    short = (f"{lb.pct(direction, 0)} of quantified published DEPs in the same direction, "
             f"{lb.pct(recall, 0)} found as DEPs; random labels give "
             f"{lb.pct(perm_share, 0)} of our DEP count")
    metrics = {"direction": direction, "recall": recall,
               "direction_as_reported": summ["as reported"]["direction"],
               "recall_as_reported": summ["as reported"]["recall"],
               "recall_quantitative": summ[prim]["recall_quantitative"],
               "recall_presence": summ[prim]["recall_presence"], "r_fc_quantitative": r_fc,
               "random_label_share": perm_share}
    return lb.answer(QUESTION, "\n".join(lines), [fig_outcomes, fig_scatter], short, verdict,
                     criterion, metrics)


if __name__ == "__main__":
    print(main())
