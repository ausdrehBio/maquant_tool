"""Q4. Can the published statistics be reconstructed from Table S3?

The paper states (section 2.7): log2 transformation, a check for normality, one-way ANOVA
followed by Tukey's post hoc test "to correct for multiple testing", DEPs at adjusted p < 0.05
and log2 fold change +-1; Spectronaut ran without imputation (section 2.3.3). It does not say
how missing values entered the ANOVA. No protein matrix is deposited, so the reconstruction
rests on the properties of Table S3 itself, checked here: the counts, the cutoffs, the
direction labels, and the fingerprints of the missing-value handling and of Tukey's test.
"""
import numpy as np
import pandas as pd

import lb_common as lb

QUESTION = "Can the published statistics be reconstructed from Table S3?"
PAPER = "one-way ANOVA, Tukey post hoc; padj < 0.05 and |log2 FC| >= 1 (section 2.7)"
NOTE = "missing-value handling not stated; reconstructed: missing = 0 on the log2 scale"
FIG1A = {("UT sw vs ctr",): 60, ("LT sw vs ctr",): 44, ("LT hu vs ctr",): 49,
         ("UT hu vs ctr",): 122,
         ("UT sw vs ctr", "LT sw vs ctr", "LT hu vs ctr", "UT hu vs ctr"): 5}


def venn_region(sets, inside):
    inner = set.intersection(*[sets[c] for c in inside])
    outer = set().union(*[sets[c] for c in sets if c not in inside])
    return inner - outer


def main():
    pub = lb.table_s3()
    titles = pub.attrs["titles"]
    n = len(pub)
    count_rows, count_bad = [], 0
    for name in lb.COMPARISON_NAMES:
        x = pub[pub.comparison == name]
        ok = len(x) == titles[name] == x.acc.nunique()
        count_bad += not ok
        count_rows.append((name, titles[name], len(x), int((x.expr == "Up").sum()),
                           int((x.expr == "Down").sum()), "yes" if ok else "NO"))
    fc_bad = pub[pub.lfc.abs() < lb.FC_CUTOFF]
    p_bad = pub[pub.padj >= lb.P_CUTOFF]
    label_bad = pub[(pub.expr == "Up") != (pub.lfc > 0)]
    violations = len(fc_bad) + len(p_bad) + len(label_bad)

    # Fig. 1A Venn of the four infected-vs-control lists
    sets = {c: set(pub[pub.comparison == c].acc) for c in max(FIG1A, key=len)}
    venn_rows = [(" & ".join(k), v, len(venn_region(sets, k))) for k, v in FIG1A.items()]
    core = sorted(pub[pub.acc.isin(venn_region(sets, max(FIG1A, key=len)))].gene.unique())

    # fingerprints of missing values set to 0
    share4 = (pub.lfc.abs() >= lb.PRESENCE_LFC).mean()
    share8 = (pub.lfc.abs() >= 8).mean()
    ident = []
    for t in lb.TISSUES:
        a = pub[pub.comparison == f"{t} sw vs ctr"].set_index("acc")
        c = pub[pub.comparison == f"{t} sw vs hu"].set_index("acc")
        both = a.index.intersection(c.index)
        same = [acc for acc in both if abs(a.lfc[acc] - c.lfc[acc]) < 1e-9]
        ident += [(t, a.gene[acc], a.lfc[acc]) for acc in same]
    zero_p = pub[pub.padj == 0]
    min_p = pub.loc[pub.padj > 0, "padj"].min()

    # does requiring ANOVA p < 0.05 change our calls? (ours, run-filtered)
    anova_drop, total_dep = 0, 0
    for t in lb.TISSUES:
        st = lb.pipeline("ours", lb.PRIMARY_LAYER)[t]
        for name in lb.COMPARISON_NAMES:
            if name.startswith(t):
                dep = st[f"{name}|call"] != 0
                total_dep += int(dep.sum())
                anova_drop += int((dep & (st["anova_p"] >= 0.05)).sum())

    if count_bad == 0 and violations == 0:
        verdict = "consistent"
    elif violations <= 0.02 * n and count_bad == 0:
        verdict = "minor deviation"
    else:
        verdict = "inconsistent"
    criterion = ("Consistent if every sheet holds the number of DEPs its title states and every "
                 "entry meets the stated cutoffs (adjusted p < 0.05, |log2 FC| >= 1) with an "
                 "Up/Down label matching its sign; minor deviation if at most 2% of the entries "
                 "violate this; inconsistent otherwise.")
    lines = [
        lb.md_table(count_rows, ["Sheet", "DEPs in the title", "Entries", "Up", "Down",
                                 "Title = entries"]),
        "",
        f"**Cutoffs and labels.** {n - violations} of {n} entries meet adjusted p < "
        f"{lb.P_CUTOFF} and |log2 FC| >= {lb.FC_CUTOFF:.0f} with a matching Up/Down label. "
        + (" ".join(f"{r.comparison}: {r.gene} has log2 FC {r.lfc:.3f} (adjusted p "
                    f"{r.padj:.4f})." for r in fc_bad.itertuples()) if len(fc_bad) else "")
        + (f" {len(p_bad)} entries have adjusted p >= {lb.P_CUTOFF}." if len(p_bad) else "")
        + (f" {len(label_bad)} labels contradict the sign." if len(label_bad) else ""),
        "",
        "**Figure 1A** (Venn diagram of the four infected-vs-control lists) against Table S3:",
        "",
        lb.md_table(venn_rows, ["Region (only these lists)", "Fig. 1A", "Table S3"]),
        "",
        f"The DEPs common to all four lists are {', '.join(core)} (paper: OASL, MX1, IFIT3, "
        "ISG15 and DRG2).",
        "",
        "**Missing values were set to 0 on the log2 scale.** Two properties of Table S3 show it:",
        "",
        f"- {lb.pct(share8)} of the entries have |log2 FC| >= 8 and {lb.pct(share4)} >= "
        f"{lb.PRESENCE_LFC:.0f}, up to {pub.lfc.abs().max():.2f} (the viral nucleoprotein). A "
        "2^19-fold change between two groups of tissue is not a ratio of intensities; it is the "
        "mean log2 intensity of a group against a group of zeros. With missing values left out, "
        "no fold change could be computed for a protein absent from one group.",
        f"- In {len(ident)} cases ({len({g for _, g, _ in ident})} proteins: "
        f"{', '.join(sorted({g for _, g, _ in ident}))}) the log2 FC of 'sw vs ctr' equals that "
        "of 'sw vs hu' of the same tissue to the last digit: the protein had no value in either "
        "reference group, so both reference means are 0.",
        "",
        f"**Tukey's test.** {len(zero_p)} entries have an adjusted p of exactly 0 (the "
        f"nucleoprotein, in {', '.join(zero_p.comparison)}), the smallest other value is "
        f"{min_p:.1e}. R's TukeyHSD computes p as 1 minus the studentized-range distribution "
        "function, which becomes exactly 0 once that function rounds to 1; a t-test p value "
        "would stay above 0. The adjustment covers the three pairwise comparisons of a protein, "
        "not the thousands of proteins tested, so it does not control the false discovery "
        "rate across the table (Q5 estimates what this means).",
        "",
        "**Reconstructed pipeline** (lb_common.paper_stats): per tissue, every protein group with "
        "a value in one of its 22 samples; log2; missing values set to 0; one-way ANOVA; Tukey's "
        "HSD (k = 3 groups, 19 residual degrees of freedom); DEP at Tukey-adjusted p < 0.05 and "
        "|log2 FC| >= 1. Whether the authors additionally required an ANOVA p < 0.05 cannot be "
        f"seen; in our data it would remove {anova_drop} of {total_dep} DEPs.",
    ]
    short = (f"{n - violations}/{n} entries meet the stated cutoffs; missing values were set to "
             f"0 ({lb.pct(share8)} of entries with |log2 FC| >= 8)")
    metrics = {"entries": n, "violations": violations, "count_mismatches": count_bad,
               "share_lfc8": share8, "identical_pairs": len(ident),
               "anova_drop": anova_drop, "our_deps": total_dep}
    return lb.answer(QUESTION, "\n".join(lines), None, short, verdict, criterion, metrics)


if __name__ == "__main__":
    print(main())
