"""Q6. Do the differential results agree with the publication?

Runs the paper pipeline (filter, imputation over N seeds, statistics; validated in Q5) on our
proteinGroups and, as the authors' own baseline, on their deposited one, and compares both with
the published fold changes and calls (Supplementary Data 5). Per comparison:
- fold-change correlation on 1:1-matched groups without an imputed value in the 12 samples of
  the two compared conditions, in both tables (imputed values are independent random draws and
  would only measure imputation noise); reference: the deposited baseline and a half-split floor
  (replicates 1-3 vs 4-6 of the published real values, Spearman-Brown corrected to 6 replicates)
- slope (reduced major axis) on the same groups, relative to the baseline's slope
- direction agreement among the published hits
- significance counts against the deposited baseline (the published counts come from a search
  with fewer missing values, Q5)
- precision and recall of our calls against the published calls
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mq_common as mq

QUESTION = "Do the differential results agree with the publication?"
PAPER = ("Supplementary Data 5: 6 strains vs mock and 3 G vs R comparisons; counts e.g. RSIG "
         "124 up/205 down, AFMX 141/404")
NOTE = "paper pipeline on ours and on deposited; groups without imputed values for correlations"


def matched(layer_name, result):
    """1:1 pairs between pipeline rows of a layer and published rows, and the reason table for
    published rows without a usable partner."""
    layer = mq.layers()[layer_name]
    pub = mq.published()
    comp, pairs = mq.match_groups(list(layer.ids), list(pub.ids))
    pos = {r: i for i, r in enumerate(result.rows)}          # layer row -> pipeline index
    pairs = pairs[pairs["a"].isin(pos)]
    pi = pairs["a"].map(pos).values
    bi = pairs["b"].values
    status = pd.Series("not identified", index=range(len(pub.matrix)))
    for _, c in comp.iterrows():
        if c["category"] == "split/merge":
            status[c["b"]] = "split/merge"
        elif c["category"] in ("identical", "1:1 overlap"):
            status[c["b"][0]] = "matched" if c["a"][0] in pos else "filtered"
    return pi, bi, status


def complete(mask, cond_a, cond_b):
    cols = [mq.SAMPLES.index(s) for s in mq.samples_of(cond_a) + mq.samples_of(cond_b)]
    return mask[:, cols].all(axis=1)


def half_split(cond_a, cond_b):
    """Correlation of published fold changes from replicates 1-3 vs 4-6 (real values only)."""
    pub = mq.published()
    ok = complete(pub.real.values, cond_a, cond_b)
    M = pub.matrix[ok]

    def fc(reps):
        return (M[[f"{cond_a}_{r}" for r in reps]].mean(axis=1)
                - M[[f"{cond_b}_{r}" for r in reps]].mean(axis=1))
    r3 = float(np.corrcoef(fc([1, 2, 3]), fc([4, 5, 6]))[0, 1])
    return r3, mq.spearman_brown(r3), int(ok.sum())


def layer_metrics(layer_name):
    res = mq.pipeline(layer_name)
    pub = mq.published()
    hits = mq.published_hits()
    pi, bi, status = matched(layer_name, res)
    out = {}
    for k, (c, a, b, _) in enumerate(mq.COMPARISONS):
        sub = complete(res.observed[pi], a, b) & complete(pub.real.values[bi], a, b)
        x = pub.fc[c].values[bi]
        f0 = res.fc[0, pi, k]
        r_sub = float(np.corrcoef(f0[sub], x[sub])[0, 1])
        slope = mq.rma_slope(x[sub], f0[sub])
        r_all = float(np.median([np.corrcoef(res.fc[s, pi, k], x)[0, 1]
                                 for s in range(res.fc.shape[0])]))
        h = hits[c].values[bi]
        is_hit = h != 0
        direction = float(np.median([(np.sign(res.fc[s, pi, k][is_hit]) == h[is_hit]).mean()
                                     for s in range(res.fc.shape[0])]))
        calls = res.calls[:, pi, k]                                   # seeds x pairs
        tp = ((calls == h) & is_hit).sum(axis=1)
        recall = float(np.median(tp / max(is_hit.sum(), 1)))
        precision = float(np.median(tp / np.maximum((calls != 0).sum(axis=1), 1)))
        union = ((calls != 0) | is_hit).sum(axis=1)
        jaccard = float(np.median(tp / np.maximum(union, 1)))
        counts = res.counts()[:, k, :]
        out[c] = {"r_sub": r_sub, "n_sub": int(sub.sum()), "slope": slope, "r_all": r_all,
                  "direction": direction, "recall": recall, "precision": precision,
                  "jaccard": jaccard, "up": float(np.median(counts[:, 0])),
                  "down": float(np.median(counts[:, 1])),
                  "up_range": (int(counts[:, 0].min()), int(counts[:, 0].max())),
                  "down_range": (int(counts[:, 1].min()), int(counts[:, 1].max())),
                  "hits_comparable": int(is_hit.sum())}
    return out, (pi, bi, status, res)


def count_rule(ours, dep):
    tol = max(0.15 * dep, 2 * np.sqrt(dep))
    d = abs(ours - dep)
    return "consistent" if d <= tol else ("minor deviation" if d <= 2 * tol else "inconsistent")


def comparison_verdict(m_o, m_d, r6):
    rules = {}
    if r6 >= 0.5:
        if m_o["r_sub"] >= m_d["r_sub"] - 0.05:
            rules["fold-change r"] = "consistent"
        elif m_o["r_sub"] >= r6:
            rules["fold-change r"] = "minor deviation"
        else:
            rules["fold-change r"] = "inconsistent"
        ratio = m_o["slope"] / m_d["slope"]
        rules["slope"] = ("consistent" if 0.8 <= ratio <= 1.25 else
                          "minor deviation" if 0.67 <= ratio <= 1.5 else "inconsistent")
    d_o, d_d = m_o["direction"], m_d["direction"]
    if d_o < 0.90 or d_o < d_d - 0.05:
        rules["direction"] = "inconsistent"
    elif d_o >= d_d - 0.02:
        rules["direction"] = "consistent"
    else:
        rules["direction"] = "minor deviation"
    rules["counts"] = mq.worst([count_rule(m_o["up"], m_d["up"]),
                                count_rule(m_o["down"], m_d["down"])])
    return mq.worst(rules.values()), rules


def discordant(pi, bi, status, res):
    """Published hits that our pipeline does not call (in < 5% of seeds), with a reason."""
    pub = mq.published()
    hits = mq.published_hits()
    pos_of = dict(zip(bi, pi))
    rows = []
    for k, (c, a, b, _) in enumerate(mq.COMPARISONS):
        for j in np.flatnonzero(hits[c].values != 0):
            h = hits[c].values[j]
            st = status[j]
            if st == "matched":
                i = pos_of[j]
                frac = float((res.calls[:, i, k] == h).mean())
                if frac >= 0.05:
                    continue
                fc_o = float(np.median(res.fc[:, i, k]))
                imp = not complete(res.observed[[i]], a, b)[0]
                pub_imp = not complete(pub.real.values[[j]], a, b)[0]
                if imp or pub_imp:
                    reason = "imputed values involved"
                elif np.sign(fc_o) != h:
                    reason = "opposite direction"
                elif abs(fc_o) > mq.FC_CUTOFF * 0.75:
                    reason = "borderline (same direction)"
                else:
                    reason = "smaller fold change"
            else:
                frac, fc_o, reason = 0.0, np.nan, st
            rows.append({"comparison": c, "row": j, "published_fc": float(pub.fc[c].values[j]),
                         "published_q": float(pub.q[c].values[j]), "ours_fc": fc_o,
                         "ours_fraction": frac, "reason": reason})
    return pd.DataFrame(rows)


def write_tables(name, res):
    layer = mq.layers()[name]
    ids = layer.pg["Majority protein IDs"].values[res.rows]
    df = pd.DataFrame({"majority_protein_ids": ids,
                       "genes": [";".join(mq.group_genes(mq.id_set(s))) for s in ids],
                       "imputed_share": res.imputed_share()})
    for k, c in enumerate(mq.COMPARISON_NAMES):
        df[f"{c}_log2fc"] = np.median(res.fc[:, :, k], axis=0)
        df[f"{c}_frac_up"] = (res.calls[:, :, k] == 1).mean(axis=0)
        df[f"{c}_frac_down"] = (res.calls[:, :, k] == -1).mean(axis=0)
    return mq.save_table(df, f"differential_{name}", index=False)


def plot(metrics, data):
    pub = mq.published()
    pi, bi, status, res = data["ours"]
    fig, axes = plt.subplots(3, 3, figsize=(11, 10.5), sharex=False, sharey=False)
    for k, (c, a, b, _) in enumerate(mq.COMPARISONS):
        ax = axes.flat[k]
        sub = complete(res.observed[pi], a, b) & complete(pub.real.values[bi], a, b)
        x = pub.fc[c].values[bi]
        y = res.fc[0, pi, k]
        ax.plot(x[~sub], y[~sub], ".", ms=1.5, color=mq.GRID, rasterized=True)
        ax.plot(x[sub], y[sub], ".", ms=1.5, color=mq.SERIES[0], alpha=0.6, rasterized=True)
        lim = np.nanpercentile(np.r_[x, y], [0.2, 99.8])
        ax.plot(lim, lim, color=mq.AXIS, lw=0.8)
        for t in (-mq.FC_CUTOFF, mq.FC_CUTOFF):
            ax.axvline(t, color=mq.GRID, lw=0.6)
            ax.axhline(t, color=mq.GRID, lw=0.6)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        m_o, m_d = metrics["ours"][c], metrics["deposited"][c]
        ax.set_title(f"{mq.comparison_label(c)}", fontsize=10)
        ax.text(0.03, 0.97, f"r = {m_o['r_sub']:.3f} (deposited {m_d['r_sub']:.3f})",
                transform=ax.transAxes, va="top", fontsize=8, color=mq.INK2)
        if k % 3 == 0:
            ax.set_ylabel("ours log2 FC")
        if k >= 6:
            ax.set_xlabel("published log2 FC")
    fig.suptitle("Fold changes, ours vs published (blue: no imputed value in either table)",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    f1 = mq.save_fig(fig, "q06_fold_changes")

    fig, ax = plt.subplots(figsize=(11, 3.8))
    hits = mq.published_hits()
    pos = np.arange(len(mq.COMPARISON_NAMES))
    width = 0.26
    series = [("published", [(hits[c] == 1).sum() for c in mq.COMPARISON_NAMES],
               [(hits[c] == -1).sum() for c in mq.COMPARISON_NAMES], mq.INK)]
    for name in [n for n in ("deposited", "ours", "ours-67") if n in metrics]:
        series.append((name, [metrics[name][c]["up"] for c in mq.COMPARISON_NAMES],
                       [metrics[name][c]["down"] for c in mq.COMPARISON_NAMES],
                       mq.LAYER_COLORS[name]))
    width = 0.8 / len(series)
    for i, (name, up, down, color) in enumerate(series):
        off = (i - (len(series) - 1) / 2) * width
        ax.bar(pos + off, up, width=width, color=color, label=name, **mq.bar_kw())
        ax.bar(pos + off, [-d for d in down], width=width, color=color, alpha=0.55, **mq.bar_kw())
    ax.axhline(0, color=mq.AXIS, lw=0.8)
    ax.set_xticks(pos, mq.COMPARISON_NAMES)
    ax.set_ylabel("significant proteins (up above, down below)")
    ax.set_title("Significance counts (seed medians for the reanalyses)")
    ax.legend(ncol=4, fontsize=8)
    mq.value_grid(ax)
    fig.tight_layout()
    f2 = mq.save_fig(fig, "q06_counts")
    return [f1, f2]


def main():
    names = [n for n in ("ours", "deposited", "ours-67") if n in mq.layers()]
    metrics, data = {}, {}
    for name in names:
        metrics[name], data[name] = layer_metrics(name)
    hs = {c: half_split(a, b) for c, a, b, _ in mq.COMPARISONS}

    rows, verdicts, rule_rows = [], {}, []
    for c in mq.COMPARISON_NAMES:
        m_o, m_d = metrics["ours"][c], metrics["deposited"][c]
        r3, r6, _ = hs[c]
        v, rules = comparison_verdict(m_o, m_d, r6)
        verdicts[c] = v
        rows.append((mq.comparison_label(c), mq.n(m_o["n_sub"]), mq.f2(m_o["r_sub"], 3),
                     mq.f2(m_d["r_sub"], 3), mq.f2(r6, 3) + ("" if r6 >= 0.5 else " (not used)"),
                     f"{mq.f2(m_o['slope'])} / {mq.f2(m_d['slope'])}",
                     f"{mq.pct(m_o['direction'])} / {mq.pct(m_d['direction'])}", v))
        hit = mq.published_hits()[c]
        rule_rows.append((mq.comparison_label(c), f"{int((hit == 1).sum())}/{int((hit == -1).sum())}",
                          f"{int(m_d['up'])}/{int(m_d['down'])}",
                          f"{int(m_o['up'])}/{int(m_o['down'])} ({m_o['up_range'][0]}–"
                          f"{m_o['up_range'][1]} / {m_o['down_range'][0]}–{m_o['down_range'][1]})",
                          f"{mq.pct(m_o['recall'], 0)} / {mq.pct(m_d['recall'], 0)}",
                          f"{mq.pct(m_o['precision'], 0)} / {mq.pct(m_d['precision'], 0)}",
                          rules["counts"]))
    verdict = mq.worst(verdicts.values())

    pi, bi, status, res = data["ours"]
    disc = discordant(pi, bi, status, res)
    disc_dep = discordant(*data["deposited"])
    reasons = pd.DataFrame({"ours": disc["reason"].value_counts(),
                            "deposited": disc_dep["reason"].value_counts()}).fillna(0).astype(int)
    reasons = reasons.sort_values("ours", ascending=False)
    id_rows = [r for r in ("not identified", "split/merge", "filtered") if r in reasons.index]
    id_extra = int(reasons.loc[id_rows, "ours"].sum() - reasons.loc[id_rows, "deposited"].sum())
    n_hits = int((mq.published_hits() != 0).sum().sum())
    pub = mq.published()
    top = disc.assign(score=disc["published_fc"].abs() * -np.log10(disc["published_q"])) \
        .sort_values("score", ascending=False).head(12)
    top_rows = [(mq.comparison_label(r.comparison),
                 ";".join(mq.group_genes(pub.ids.iloc[r.row]))[:28] or pub.genes.iloc[r.row][:28],
                 mq.f2(r.published_fc), mq.f2(r.ours_fc), mq.pct(r.ours_fraction, 0), r.reason)
                for r in top.itertuples()]

    # disagreement by abundance and imputation (ours vs published, pooled over comparisons)
    x = np.stack([pub.fc[c].values[bi] for c in mq.COMPARISON_NAMES], axis=1)
    y = np.median(res.fc[:, pi, :], axis=0)
    diff = np.abs(y - x)
    abundance = pub.matrix.where(pub.real).mean(axis=1).values[bi]
    q = pd.qcut(abundance, 4, labels=["lowest quarter", "second", "third", "highest quarter"])
    imp = res.imputed_share()[pi]
    ib = pd.cut(imp, [-0.01, 0.0, 0.1, 0.3, 1.0], labels=["none", "≤ 10%", "10–30%", "> 30%"])
    by_ab = pd.Series(np.median(diff, axis=1)).groupby(q, observed=True).median()
    by_imp = pd.Series(np.median(diff, axis=1)).groupby(ib, observed=True).median()
    n_imp = pd.Series(ib).value_counts()

    tables = [write_tables(name, data[name][3]) for name in names]
    figures = plot(metrics, data)

    criterion = ("Per comparison, the worst of four rules. Fold-change r (groups without imputed "
                 "values): consistent if ours ≥ deposited − 0.05, minor deviation if ≥ the "
                 "half-split floor; used only where that floor is ≥ 0.5. Slope ratio ours/"
                 "deposited within 0.8–1.25 (minor: 0.67–1.5). Direction agreement among "
                 "published hits ≥ deposited − 0.02 and ≥ 90% (inconsistent below 90% or below "
                 "deposited − 0.05). Counts within max(15%, 2√n) of deposited (minor: twice "
                 "that). The question takes the worst comparison.")
    lines = [
        f"The paper pipeline (Q5) ran on our table and on the deposited one with {mq.N_SEEDS} "
        "imputation seeds. **Agreement of fold changes with the published ones** (1:1-matched "
        "groups without imputed values in the compared conditions; ours / deposited where two "
        "values are given):",
        "",
        mq.md_table(rows, ["Comparison", "Groups", "r ours", "r deposited", "Half-split floor",
                           "Slope", "Direction among hits", "Verdict"]),
        "",
        "**Significance counts and overlap with the published calls** (up/down; seed medians and "
        "ranges; recall: share of published hits we also call; precision: share of our calls "
        "among matched groups that the paper also made):",
        "",
        mq.md_table(rule_rows, ["Comparison", "Published", "Deposited", "Ours (range)",
                                "Recall ours / deposited", "Precision ours / deposited",
                                "Counts rule"]),
        "",
        "Both reanalyses find about 30% fewer significant proteins than the publication, "
        "and to the same extent: the search behind the published table had fewer missing "
        "values (Q5). Against that baseline, our counts, fold changes and directions agree as "
        "closely as the authors' own deposited search does.",
        "",
        f"**Published hits not reproduced** (called in fewer than 5% of seeds): {len(disc)} of "
        f"{n_hits} published calls over the 9 comparisons in our reanalysis, {len(disc_dep)} in "
        "the reanalysis of the authors' own deposited table. Reasons, in the order they are "
        "checked:",
        "",
        mq.md_table([(k, int(r["ours"]), int(r["deposited"])) for k, r in reasons.iterrows()],
                    ["Reason", "Ours", "Deposited"]),
        "",
        "Imputed values: a missing value in one of the 12 compared samples of either table. "
        "Borderline: same direction and more than 75% of the fold-change threshold. Not "
        "identified, split/merge and filtered: no usable partner group in the reanalysed table. "
        "The authors' deposited search loses calls in the same proportions, so most losses come "
        f"from re-searching the data. Ours loses {len(disc) - len(disc_dep)} more, mainly "
        f"through identification and grouping ({id_extra} more not identified, split/merge or "
        "filtered; Q3). Without imputed values, "
        f"{int(reasons.loc['opposite direction', 'ours']) if 'opposite direction' in reasons.index else 0}"
        " published hits change direction in our data.",
        "",
        "Strongest of them (published |log2FC| times −log10 adjusted p):",
        "",
        mq.md_table(top_rows, ["Comparison", "Genes", "Published FC", "Ours FC",
                               "Our seeds significant", "Reason"]),
        "",
        "**Where fold changes disagree** (median |ours − published| per group, pooled over the "
        "comparisons):",
        "",
        mq.md_table([(str(k), mq.f2(v, 3)) for k, v in by_ab.items()],
                    ["Published abundance", "Median |Δ log2FC|"]),
        "",
        mq.md_table([(str(k), int(n_imp.get(k, 0)), mq.f2(v, 3)) for k, v in by_imp.items()],
                    ["Imputed share in our 42 samples", "Groups", "Median |Δ log2FC|"]),
        "",
        "Disagreement grows towards low abundance and with the share of imputed values, as "
        "expected for random imputation draws near the detection limit.",
        "",
        "Per-protein results (median fold change and the share of seeds with an up or down call "
        "per comparison) are written to " + ", ".join(f"`{t}`" for t in tables) + " (not in git).",
    ]
    worst_c = [c for c, v in verdicts.items() if v == verdict]
    short = (f"fold-change r {mq.f2(min(metrics['ours'][c]['r_sub'] for c in mq.STRAINS), 3)}–"
             f"{mq.f2(max(metrics['ours'][c]['r_sub'] for c in mq.STRAINS), 3)} vs mock "
             f"(deposited {mq.f2(min(metrics['deposited'][c]['r_sub'] for c in mq.STRAINS), 3)}–"
             f"{mq.f2(max(metrics['deposited'][c]['r_sub'] for c in mq.STRAINS), 3)}); counts as "
             "the deposited search")
    hits = mq.published_hits()
    m = {"worst_comparisons": ",".join(worst_c), "comparisons": {}}
    for c in mq.COMPARISON_NAMES:
        o, d = metrics["ours"][c], metrics["deposited"][c]
        m["comparisons"][c] = {
            "r_ours": o["r_sub"], "r_deposited": d["r_sub"], "half_split": hs[c][1],
            "slope_ours": o["slope"], "slope_deposited": d["slope"],
            "direction_ours": o["direction"], "direction_deposited": d["direction"],
            "up_ours": o["up"], "down_ours": o["down"], "up_deposited": d["up"],
            "down_deposited": d["down"], "up_published": int((hits[c] == 1).sum()),
            "down_published": int((hits[c] == -1).sum()),
            "recall_ours": o["recall"], "recall_deposited": d["recall"],
            "precision_ours": o["precision"], "precision_deposited": d["precision"],
            "verdict": verdicts[c]}
    return mq.answer(QUESTION, "\n".join(lines), figures, short, verdict, criterion, m)


if __name__ == "__main__":
    print(main())
