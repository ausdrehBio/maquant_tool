"""Q6. Do the paper's proteomics claims hold with our data?

Every statement of the paper about proteins (sections 3.1-3.5 and 4) is split into parts that
can be checked one by one. A claim about a differentially abundant protein is reproduced if the
protein (its Table S3 accession) is a DEP in the stated direction in our data under the paper
pipeline; viral and structural claims have their own tests, stated in the table. Judged on the
run-filtered layer; the as-reported layer is shown alongside.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import lb_common as lb

QUESTION = "Do the paper's proteomics claims hold with our data?"
PAPER = "sections 3.1-3.5 and the discussion"
NOTE = "claims about transcripts (RNA-seq, qPCR) are outside this reanalysis"
VIRAL_MIN_SHARE = 0.5   # "quantified": detected in at least half of a group's samples

# (claim, comparison, gene as in Table S3, expected direction)
DEP_CLAIMS = [
    ("DEPs common to both infections and tissues (3.1)", c, g, +1)
    for g in ("OASL", "MX1", "IFIT3", "ISG15", "DRG2")
    for c in ("UT sw vs ctr", "UT hu vs ctr", "LT sw vs ctr", "LT hu vs ctr")
] + [
    # MX1 and ISG15 in the lower trachea are parts of the claim above already
    ("antiviral ISGs MX2 and IFI6 up in the lower trachea after both infections (3.2)", c, g, +1)
    for g in ("MX2", "IFI6") for c in ("LT sw vs ctr", "LT hu vs ctr")
] + [
    ("apoptosis-related proteins up after swH1N1, lower trachea (3.5)", "LT sw vs ctr", "ANXA5",
     +1),
    ("apoptosis-related proteins up after swH1N1, lower trachea (3.5)", "LT sw vs ctr",
     "LOC100522887", +1),
    ("huH1N1, lower trachea (3.5, 4)", "LT hu vs ctr", "ANXA5", +1),
    ("huH1N1, lower trachea (3.5, 4)", "LT hu vs ctr", "PDCD10", -1),
    ("huH1N1, lower trachea (3.5, 4)", "LT hu vs ctr", "CCL26", -1),
    ("huH1N1, lower trachea (3.5, 4)", "LT hu vs ctr", "IL13", -1),
    ("DFFA down in the upper trachea after both infections (3.5)", "UT sw vs ctr", "DFFA", -1),
    ("DFFA down in the upper trachea after both infections (3.5)", "UT hu vs ctr", "DFFA", -1),
    ("swH1N1 vs huH1N1, lower trachea (3.4, 4)", "LT sw vs hu", "STAT1", -1),
    ("swH1N1 vs huH1N1, lower trachea (3.4, 4)", "LT sw vs hu", "EP300", +1),
    ("swH1N1 vs huH1N1, lower trachea (3.4, 4)", "LT sw vs hu", "KRT2", +1),
    ("swH1N1 vs huH1N1, lower trachea (3.4, 4)", "LT sw vs hu", "PLXNA4", +1),
    ("swH1N1 vs huH1N1, lower trachea (3.4, 4)", "LT sw vs hu", "S100A1", -1),
    ("swH1N1 vs huH1N1, lower trachea (3.4, 4)", "LT sw vs ctr", "PLEC", +1),
    ("serpins up after swH1N1, lower trachea (4)", "LT sw vs ctr", "SERPINB10", +1),
    ("serpins up after swH1N1, lower trachea (4)", "LT sw vs ctr", "LOC106504547", +1),
    ("cytoskeleton proteins after huH1N1, lower trachea (4)", "LT hu vs ctr", "GSN", -1),
    ("cytoskeleton proteins after huH1N1, lower trachea (4)", "LT hu vs ctr", "PLEC", +1),
    ("cytoskeleton proteins after huH1N1, lower trachea (4)", "LT hu vs ctr", "STRIP1", +1),
]
GENE_LABELS = {"LOC100522887": "caspase-13 (LOC100522887)",
               "LOC106504547": "serpin A3-8 (LOC106504547)"}


def dep_part(comparison, gene, direction, key="ours"):
    pub = lb.table_s3()
    row = pub[(pub.comparison == comparison) & (pub.gene == gene)].iloc[0]
    out = {"paper": f"{row.lfc:+.2f}"}
    s = lb.search(key)
    gid = s.group_of_accession(row.acc)
    ok = {}
    for layer in lb.LAYERS:
        fr = lb.comparison_frame(key, layer, comparison)
        if gid is None or gid not in fr.index:
            out[layer] = "not quantified"
            ok[layer] = False
            continue
        r = fr.loc[gid]
        out[layer] = f"{r.lfc:+.2f} (p {r.p:.2g})"
        ok[layer] = int(r.call) == direction
    return out, ok


def group_count(layer, comparison):
    return int((lb.comparison_frame("ours", layer, comparison).call != 0).sum())


def go_share(layer, comparison, go):
    s = lb.search("ours")
    fr = lb.comparison_frame("ours", layer, comparison)
    dep = fr.index[fr.call != 0]
    if len(dep) == 0:
        return np.nan
    hit = [gid for gid in dep if any(lb.accession_token(i) in go
                                     for i in s.pg.loc[gid, "Majority protein IDs"].split(";"))]
    return len(hit) / len(dep)


def main():
    s = lb.search("ours")
    rows = []

    # viral claims (3.1): quantified = detected in >= half of the group's samples
    det = s.detected
    viral = s.viral_groups

    def viral_share(protein, tissue, group):
        gids = [g for g, members in viral.items() if any(p == protein for _, p in members)]
        if not gids:
            return 0.0
        c = lb.samples_of(tissue, group)
        return max(det.loc[g, c].mean() for g in gids)
    for protein in ("NP", "HA"):
        for t in lb.TISSUES:
            sh = viral_share(protein, t, "sw")
            rows.append(("viral proteins (3.1)", f"{protein} quantified after swH1N1, "
                         f"{lb.TISSUE_LABELS[t]}", "quantified", f"detected in {lb.pct(sh, 0)} "
                         "of samples", "", sh >= VIRAL_MIN_SHARE, sh >= VIRAL_MIN_SHARE))
    hu_max = max((det.loc[g, lb.samples_of(t, "hu")].mean() for g in viral for t in lb.TISSUES),
                 default=0.0)
    rows.append(("viral proteins (3.1)", "no viral protein quantified after huH1N1", "none",
                 f"best viral group detected in {lb.pct(hu_max, 0)} of a huH1N1 tissue's samples",
                 "", hu_max < VIRAL_MIN_SHARE, hu_max < VIRAL_MIN_SHARE))
    found = sorted({p for members in viral.values() for _, p in members})
    need = ["M1", "NS1", "PB1", "HA", "NP"]
    missing = [p for p in need if p not in found]
    rows.append(("viral proteins (3.1)", "M1, NS1, PB1, HA and NP identified", "identified",
                 "identified: " + ", ".join(found) + (f"; missing: {', '.join(missing)}"
                                                       if missing else ""),
                 "", not missing, not missing))

    # number of DEPs: upper > lower trachea for each strain (3.1)
    pub = lb.table_s3()
    for g in ("sw", "hu"):
        pu = int((pub.comparison == f"UT {g} vs ctr").sum())
        pl = int((pub.comparison == f"LT {g} vs ctr").sum())
        res = {}
        for layer in lb.LAYERS:
            u, l_ = group_count(layer, f"UT {g} vs ctr"), group_count(layer, f"LT {g} vs ctr")
            res[layer] = (f"{u} vs {l_}", u > l_)
        rows.append(("more DEPs in the upper than the lower trachea (3.1)",
                     f"{g}H1N1: upper > lower", f"{pu} vs {pl}", res["run-filtered"][0],
                     res["as reported"][0], res["run-filtered"][1], res["as reported"][1]))

    # DEP claims
    for claim, comp, gene, direction in DEP_CLAIMS:
        out, ok = dep_part(comp, gene, direction)
        label = f"{GENE_LABELS.get(gene, gene)} {'up' if direction > 0 else 'down'}, {comp}"
        rows.append((claim, label, out["paper"], out["run-filtered"], out["as reported"],
                     ok["run-filtered"], ok["as reported"]))

    # pigs 9 and 10 (lowest viral load) cluster with the infected pigs (3.2, Fig. 2A)
    isg = lb.isg_score("ours")
    cmax = isg[lb.samples_of("LT", "ctr")].max()
    p9, p10 = isg["LT_sw_p09"], isg["LT_sw_p10"]
    ok_p = bool(p9 > cmax and p10 > cmax)
    rows.append(("pigs 9 and 10 cluster with the infected pigs in proteomics (3.2)",
                 "ISG score of pigs 9 and 10 above all controls, lower trachea", "yes (Fig. 2A)",
                 f"{p9:.2f} and {p10:.2f} vs control max {cmax:.2f}", "", ok_p, ok_p))

    # GO:0051607 share larger in the lower trachea (3.1)
    go = lb.go_set()
    if go is not None:
        for g in ("sw", "hu"):
            res = {}
            for layer in lb.LAYERS:
                su = go_share(layer, f"UT {g} vs ctr", go)
                sl = go_share(layer, f"LT {g} vs ctr", go)
                res[layer] = (f"{lb.pct(sl, 0)} vs {lb.pct(su, 0)}", sl > su)
            paper = {"sw": "9% vs 5%", "hu": "9% vs 4%"}[g]
            rows.append(("larger share of GO:0051607 DEPs in the lower trachea (3.1)",
                         f"{g}H1N1: lower > upper", paper, res["run-filtered"][0],
                         res["as reported"][0], res["run-filtered"][1], res["as reported"][1]))

    df = pd.DataFrame(rows, columns=["claim", "part", "paper", "ours run-filtered",
                                     "ours as reported", "ok", "ok_as_reported"])
    lb.save_table(df, "q06_claim_parts", index=False)
    n_ok, n = int(df.ok.sum()), len(df)
    n_ok_rep = int(df.ok_as_reported.sum())
    share = n_ok / n

    claim_rows = []
    for claim, x in df.groupby("claim", sort=False):
        claim_rows.append((claim, f"{int(x.ok.sum())}/{len(x)}",
                           f"{int(x.ok_as_reported.sum())}/{len(x)}"))
    part_rows = [(r.claim, r.part, r.paper, r[3], r[4], "yes" if r.ok else "no")
                 for r in df.itertuples(index=False)]

    # figure: ISG proteins per sample (as Fig. 2, global proteomics)
    gids = lb.isg_groups("ours")
    m = np.log2(s.lfq.loc[list(gids.values())])
    m.index = list(gids)
    fig, axes = plt.subplots(1, 2, figsize=(10, 2.6), sharey=True)
    for ax, t in zip(axes, lb.TISSUES):
        cols = lb.samples_of(t)
        x = m[cols]
        z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
        im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_xticks(range(len(cols)), [str(lb.pig_of(c)) for c in cols], fontsize=7)
        ax.set_yticks(range(len(z)), z.index)
        for i, c in enumerate(cols):
            ax.add_patch(plt.Rectangle((i - 0.5, -1.3), 1, 0.6,
                                       color=lb.GROUP_COLORS[lb.group_of(c)], clip_on=False))
        ax.set_title(lb.TISSUE_LABELS[t], pad=14)
        ax.set_xlabel("pig")
    fig.colorbar(im, ax=axes, shrink=0.8, label="z-score")
    fig.suptitle("Interferon-stimulated proteins of Fig. 2, ours as reported (log2 LFQ; white: "
                 "no value)", fontsize=10, fontweight="bold", y=1.08)
    figure = lb.save_fig(fig, "q06_isg_proteins")

    # numbers for the text
    single = df[df.claim.str.contains("apoptosis|huH1N1, lower|DFFA|swH1N1 vs huH1N1|serpins|"
                                      "cytoskeleton")]
    n_single, n_single_ok = len(single), int(single.ok.sum())
    n_single_presence = int(single.paper.map(lambda v: abs(float(v)) >= lb.PRESENCE_LFC
                                             if v[:1] in "+-" else False).sum())
    not_ident = single[single["ours run-filtered"] == "not quantified"].part.str.split().str[0]
    not_identified_text = (f"{', '.join(sorted(set(not_ident)))} is not identified and "
                           if len(not_ident) else "")
    drg2 = df[df.part.str.startswith("DRG2")]
    drg2_text = "; ".join(f"{r.part.split(', ')[1]}: {r[3]}" for r in drg2.itertuples(index=False))
    isg_parts = df[df.part.str.match(r"(IFIT3|OASL) up")]
    lost = isg_parts[~isg_parts.ok & isg_parts.ok_as_reported]
    isg_layer_note = (
        f"The two layers differ for low-abundance ISGs: {len(lost)} parts of IFIT3 and OASL are "
        "DEPs as reported but not run-filtered. Run filtering removes their values in controls "
        "where no precursor passes; set to 0, these inflate the variance of the control group, "
        "and Tukey's test loses significance although the fold change grows."
        if len(lost) else "")

    verdict = lb.grade(share, 0.8, 0.6)
    criterion = ("Consistent if at least 80% of the claim parts are reproduced on the "
                 "run-filtered layer; minor deviation at 60%; inconsistent below.")
    lines = [
        f"**{n_ok} of {n} claim parts reproduced** ({lb.pct(share, 0)}; as reported: "
        f"{n_ok_rep} of {n}).",
        "",
        lb.md_table(claim_rows, ["Claim", "Reproduced (run-filtered)", "As reported"]),
        "",
        "Parts (log2 FC as group-mean difference with missing values set to 0, Tukey-adjusted "
        "p):",
        "",
        lb.md_table(part_rows, ["Claim", "Part", "Paper", "Ours, run-filtered",
                                "Ours, as reported", "Reproduced"]),
        "",
        "**What holds:** the interferon response. MX1, MX2, ISG15 and IFI6 rise after both "
        "infections in our data, the swH1N1 tissue carries viral protein and the huH1N1 tissue "
        "does not, the upper trachea has more DEPs than the lower, and the GO:0051607 share is "
        f"larger in the lower trachea. {isg_layer_note}",
        "",
        f"**What does not hold:** the single-protein claims about apoptosis, the cytoskeleton, "
        f"STAT1/EP300 and the serpins ({n_single_ok} of {n_single} parts). "
        f"{n_single_presence} of these rest on published log2 FC of 4 or more, i.e. on a "
        "protein with values in one group and none in another in the authors' data; in ours "
        f"{not_identified_text}the others show no significant difference. DRG2, the fifth "
        "protein of the paper's common core, does not change in our data either "
        f"({drg2_text}).",
    ]
    short = f"{n_ok}/{n} claim parts reproduced (as reported: {n_ok_rep}/{n})"
    metrics = {"parts": n, "reproduced": n_ok, "reproduced_as_reported": n_ok_rep, "share": share}
    return lb.answer(QUESTION, "\n".join(lines), figure, short, verdict, criterion, metrics)


if __name__ == "__main__":
    print(main())
