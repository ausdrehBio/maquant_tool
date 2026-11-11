"""Q3. How do the identifications compare?

Counts protein groups in each table, matches them between tables (connected components of
shared protein IDs, see mq_common.match_groups) and compares the identification depth per raw
file between our search and the authors' deposited one: groups detected, identified by MS/MS,
transferred by match between runs, quantified by LFQ, and the summed MS/MS counts.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mq_common as mq

QUESTION = "How do the identifications compare?"
PAPER = "6,052 protein groups after the paper's filter (Supplementary Data 5)"
NOTE = "per raw file vs the deposited search; group matching by shared protein IDs"


def totals(layer):
    clean = layer.clean
    pg = layer.pg
    return {
        "protein groups": len(pg),
        "decoys": int(layer.decoy.sum()),
        "contaminants": int(layer.contaminant.sum()),
        "only identified by site": int(layer.site_only.sum()),
        "clean set": int(clean.sum()),
        "paper filter (Intensity detection)": int(mq.paper_filter(layer).sum()),
        "paper filter (LFQ detection)": int(mq.paper_filter(layer, "LFQ intensity").sum()),
    }


def depth(layer):
    """Per-sample identification depth of the clean set."""
    clean = layer.clean.values
    idt = layer.matrix("Identification type", numeric=False)[clean]
    return pd.DataFrame({
        "detected": (layer.matrix("Intensity")[clean] > 0).sum(),
        "by MS/MS": (idt == "By MS/MS").sum(),
        "by matching": (idt == "By matching").sum(),
        "LFQ": (layer.matrix("LFQ intensity")[clean] > 0).sum(),
        "MS/MS count": layer.matrix("MS/MS count")[clean].sum(),
    })


def coverage(components, n_b):
    """Share of B groups (published) that have a partner in A (1:1 or split/merge)."""
    has = components[components["category"] != "only B"]
    covered = sum(len(b) for b in has["b"])
    return covered / n_b


def not_in_database(ids_list):
    known = mq.fasta_ids()
    return sum(1 for s in ids_list if s and not any(t in known for t in s))


def extra_run_only(layer):
    """Clean groups identified by MS/MS only in runs the paper does not use."""
    if not layer.extra.get("Identification type"):
        return 0
    used = (layer.matrix("Identification type", numeric=False) == "By MS/MS").any(axis=1)
    extra = (layer.extra_matrix("Identification type", numeric=False) == "By MS/MS").any(axis=1)
    return int((layer.clean & extra & ~used).sum())


def viral_grouping(layer):
    out = {}
    for acc in mq.pr8_proteins():
        rows = [s for s in layer.ids if acc in s]
        if not rows:
            out[acc] = "not identified"
            continue
        members = sorted(rows[0])
        partners = [mq.entry_name(t) for t in members if t != acc]
        out[acc] = ", ".join(partners) if partners else "alone"
    return out


def plot(ratio, cats):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [2.2, 1]})
    ax = axes[0]
    x = np.arange(len(mq.SAMPLES))
    for k, (name, color) in enumerate([("by MS/MS", mq.SERIES[0]), ("detected", mq.INK2),
                                       ("LFQ", mq.SERIES[2]), ("by matching", mq.SERIES[1])]):
        ax.plot(x, ratio[name].values, "o", ms=3.5, color=color, label=name)
    ax.axhline(1, color=mq.AXIS, lw=0.8)
    ax.set_xticks(x, mq.SAMPLES)
    mq.rotate_xticks(ax, size=6, angle=90)
    ax.set_ylabel("ours / deposited, per raw file")
    ax.set_title("Identification depth per raw file (clean set)")
    ax.legend(ncol=4, loc="lower left")
    mq.value_grid(ax)
    ax = axes[1]
    order = ["identical", "1:1 overlap", "split/merge", "only A", "only B"]
    colors = [mq.SERIES[0], "#86b6ef", mq.SERIES[3], mq.SERIES[1], mq.MUTED]
    labels = list(cats)
    left = np.zeros(len(labels))
    for c, color in zip(order, colors):
        vals = np.array([cats[l][c] for l in labels], dtype=float)
        ax.barh(labels, vals, left=left, color=color, label=c, **mq.bar_kw())
        left += vals
    ax.invert_yaxis()
    ax.set_xlabel("components")
    ax.set_title("Matching of protein groups")
    ax.legend(fontsize=7, loc="lower right")
    mq.value_grid(ax, "x")
    fig.tight_layout()
    return mq.save_fig(fig, "q03_identifications")


def main():
    lay = mq.layers()
    pub = mq.published()
    tot = pd.DataFrame({name: totals(layer) for name, layer in lay.items()})
    tot_rows = [(k, *[mq.n(v) for v in tot.loc[k]]) for k in tot.index]
    tot_rows.append(("published table (after the paper's filter)", *(
        [mq.n(len(pub.matrix))] + ["" for _ in range(len(lay) - 1)])))

    pairs = [("ours", "published"), ("deposited", "published"), ("ours", "deposited")]
    if "ours-67" in lay:
        pairs += [("ours-67", "published"), ("ours-67", "deposited"), ("ours", "ours-67")]
    ids = {name: list(layer.ids) for name, layer in lay.items()}
    ids["published"] = list(pub.ids)
    cats, cov, unmatched = {}, {}, {}
    pub_mean = pub.matrix.where(pub.real).mean(axis=1)
    pub_rank = pub_mean.rank(pct=True)
    for a, b in pairs:
        comp, _ = mq.match_groups(ids[a], ids[b])
        cats[f"{a} vs {b}"] = mq.category_counts(comp)
        if b == "published":
            cov[a] = coverage(comp, len(ids[b]))
            only_b = [x[0] for x in comp.loc[comp["category"] == "only B", "b"]]
            only_a = [x[0] for x in comp.loc[comp["category"] == "only A", "a"]]
            kept = mq.paper_filter(lay[a]).values
            unmatched[a] = {
                "n_pub": len(only_b),
                "rank": float(pub_rank.iloc[only_b].median()),
                "real": float(pub.real.iloc[only_b].mean(axis=1).median()),
                "n_own": len(only_a),
                "own_filtered": float(1 - kept[only_a].mean()) if only_a else float("nan"),
            }
    match_rows = [(k, *[mq.n(v[c]) for c in ["identical", "1:1 overlap", "split/merge",
                                             "only A", "only B"]]) for k, v in cats.items()]
    pub_missing_db = not_in_database(ids["published"])

    d_ours, d_dep = depth(lay["ours"]), depth(lay["deposited"])
    ratio = d_ours / d_dep
    med = ratio.median()
    dev = (ratio - med).abs().max()
    depth_rows = [(k, mq.n(d_dep[k].median()), mq.n(d_ours[k].median()), mq.f2(med[k], 3),
                   f"{mq.f2(ratio[k].min(), 3)}–{mq.f2(ratio[k].max(), 3)}", mq.f2(dev[k], 3))
                  for k in d_ours.columns]
    mbr_ours = (d_ours["by matching"] / d_ours["detected"]).median()
    mbr_dep = (d_dep["by matching"] / d_dep["detected"]).median()
    extra_only = extra_run_only(lay["deposited"])

    viral = pd.DataFrame({name: viral_grouping(layer) for name, layer in lay.items()})
    viral_rows = [(mq.entry_name(acc), *viral.loc[acc]) for acc in viral.index]
    grouped_diff = [acc for acc in viral.index if viral.at[acc, "ours"] != viral.at[acc, "deposited"]]

    figure = plot(ratio, cats)

    msms_med = med["by MS/MS"]
    max_dev = max(dev["by MS/MS"], dev["MS/MS count"])
    cov_gap = cov["deposited"] - cov["ours"]
    if mq.band(msms_med, 0.9, 1.1) and max_dev <= 0.05 and cov_gap <= 0.05:
        verdict = "consistent"
    elif mq.band(msms_med, 0.8, 1.25) and cov_gap <= 0.10:
        verdict = "minor deviation"
    else:
        verdict = "inconsistent"
    criterion = ("Consistent if the groups identified by MS/MS per raw file are within ±10% of the "
                 "deposited search (median ratio), no file deviates from that median ratio by more "
                 "than 0.05 (a version change affects all files alike; an outlier file points to a "
                 "workflow error), and our coverage of the published groups is at most 5 points "
                 "below the deposited search's.")
    lines = [
        mq.md_table(tot_rows, ["Protein groups", *lay.keys()]),
        "",
        f"**Matching.** Groups are matched through shared protein IDs (Ensembl version suffixes "
        f"removed). The published table covers {mq.n(len(pub.matrix))} filtered groups. Our "
        f"search has a partner for **{mq.pct(cov['ours'])}** of them, the deposited search for "
        f"{mq.pct(cov['deposited'])}"
        + (f", ours-67 for {mq.pct(cov['ours-67'])}" if "ours-67" in cov else "")
        + f". {pub_missing_db} published groups have no ID in the authors' FASTA files as "
        "deposited; they probably come from the newer Ensembl release behind the versioned IDs of "
        "the published table and cannot be found by any search with the deposited database.",
        "",
        mq.md_table(match_rows, ["Tables (A vs B)", "Identical", "1:1 overlap", "Split/merge",
                                 "Only A", "Only B"]),
        "",
        mq.md_table([(a, u["n_pub"], mq.pct(u["rank"], 0), mq.pct(u["real"], 0), u["n_own"],
                      mq.pct(u["own_filtered"], 0)) for a, u in unmatched.items()],
                    ["Search", "Published groups without partner",
                     "Their median abundance rank among published", "Their median share of real "
                     "values", "Own groups without published partner", "Of these, removed by the "
                     "paper's filter"]),
        "",
        f"Published groups without a partner in our search sit low in the abundance ranking "
        f"(median rank {mq.pct(unmatched['ours']['rank'], 0)}), although most of their published "
        f"values are real ({mq.pct(unmatched['ours']['real'], 0)}): the authors' search found "
        f"them, ours did not, near the detection limit. Of our groups without a published "
        f"partner, {mq.pct(unmatched['ours']['own_filtered'], 0)} are removed by the paper's "
        "filter anyway. The deposited search shows the same pattern at a slightly smaller scale.",
        "",
        "**Depth per raw file** (clean set; median over the 42 runs; deviation: largest distance "
        "of a file's ratio from the median ratio):",
        "",
        mq.md_table(depth_rows, ["Per raw file", "Deposited", "Ours", "Ours/deposited",
                                 "Range", "Max. deviation"]),
        "",
        f"Our search identifies slightly more groups by MS/MS per file ({mq.f2(msms_med, 3)}×) "
        f"but transfers fewer by match between runs: {mq.pct(mbr_ours)} of the detected groups "
        f"per file vs {mq.pct(mbr_dep)} in the deposited search. Plausible causes are the "
        "stricter match-between-runs algorithm of MaxQuant 2.x and the 42 instead of 67 donor "
        f"runs; Part B separates the two. {extra_only} groups of the deposited table were "
        "identified by MS/MS only in the 25 runs the paper does not use; they cannot appear in a "
        "42-run search except by matching.",
        "",
        "**Viral protein groups:** members grouped with each PR8 protein (entry names):",
        "",
        mq.md_table(viral_rows, ["PR8 protein", *viral.columns]),
        "",
        (f"Grouped differently: {', '.join(mq.entry_name(a) for a in grouped_diff)} (MaxQuant "
         "2.9.0 merges it with its SC35M homolog, 1.5.6.2 keeps them apart). "
         if grouped_diff else "All PR8 proteins are grouped alike in both searches. ")
        + "The questions below pick, for each PR8 protein, the group that contains it.",
    ]
    short = (f"MS/MS identifications per file {mq.f2(msms_med, 3)}× deposited; "
             f"{mq.pct(cov['ours'])} of published groups matched")
    metrics = {"msms_ratio": float(msms_med), "max_file_deviation": float(max_dev),
               "coverage_ours": float(cov["ours"]), "coverage_deposited": float(cov["deposited"]),
               "msms_ours": float(d_ours["by MS/MS"].median()),
               "msms_deposited": float(d_dep["by MS/MS"].median()),
               "msms_count_ours": float(d_ours["MS/MS count"].median()),
               "msms_count_deposited": float(d_dep["MS/MS count"].median())}
    return mq.answer(QUESTION, "\n".join(lines), figure, short, verdict, criterion, metrics)


if __name__ == "__main__":
    print(main())
