"""Q7. Does the reanalysis reproduce the proteins reported in the paper?

Table 1 of Liljedahl et al. 2017 lists 23 proteins with significant abundance differences
(one-way ANOVA with Tukey; directions from the shotgun data, most verified by SRM). Column 6
gives db/db vehicle vs healthy control (here contrast VEH-CTRL), column 7 the treatment effect
vs db/db vehicle (here LIRA-VEH and INS-VEH). Directions marked as trends (P < 0.1) in the
paper are flagged.

Proteins are looked up by accession in the current Swiss-Prot FASTA, else by gene symbol,
because some accessions of the paper were TrEMBL entries or isoforms, and some genes were
renamed (Slc9a3r1/2 -> Nherf1/2, Car2/3 -> Ca2/3). A direction agrees when the log2FC has the
same sign; it is confirmed when, in addition, the adjusted p-value is below 0.05.
"""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import mq_common as mq

QUESTION = "Does the reanalysis reproduce the proteins reported in the paper?"
PAPER = "23 proteins in Table 1"
NOTE = "same sign of log2FC; confirmed = adj. p < 0.05 in MSstats"

# (accession in the paper, short name in the paper, gene symbols (current first),
#  db/db vs db/+ ("up", "down", "ns"), {contrast: direction}, trend flags)
TABLE1 = [
    ("Q9CPU0", "GLO1", ["Glo1"], "up", {}, ()),
    ("Q61847", "MEP1B", ["Mep1b"], "down", {}, ()),
    ("P16015", "CAR3", ["Ca3", "Car3"], "down", {}, ()),
    ("P15947", "KLK1", ["Klk1"], "down", {}, ()),
    ("P16045", "LGAL1", ["Lgals1"], "down", {}, ()),
    ("Q9WVM8", "AADAT", ["Aadat"], "down", {}, ()),
    ("P05202", "GOT2", ["Got2"], "down", {}, ()),
    ("P26350", "PTMA", ["Ptma"], "up", {}, ()),
    ("P54071", "IDH2", ["Idh2"], "down", {}, ()),
    ("P14206", "LAMR1", ["Rpsa"], "down", {}, ()),
    ("Q61171", "PRDX2", ["Prdx2"], "down", {"LIRA-VEH": "up"}, ()),
    ("Q8CHT0", "ALDH4A1", ["Aldh4a1"], "down", {"LIRA-VEH": "up"}, ()),
    ("P70441-2", "NHERF1", ["Nherf1", "Slc9a3r1"], "up", {"LIRA-VEH": "down"}, ("LIRA-VEH",)),
    ("Q3TWL8", "PSAP", ["Psap"], "up", {"LIRA-VEH": "down"}, ("VEH-CTRL",)),
    ("P83940", "ELOC", ["Eloc", "Tceb1"], "up", {"LIRA-VEH": "down"}, ()),
    ("Q3UI46", "NEDD8", ["Nedd8"], "ns", {"LIRA-VEH": "down"}, ()),
    ("P00920", "CAR2", ["Ca2", "Car2"], "down", {"INS-VEH": "up"}, ("INS-VEH",)),
    ("Q99K41", "EMIL1", ["Emilin1"], "down", {"INS-VEH": "up"}, ()),
    ("Q9CQM9", "GLRX3", ["Glrx3"], "ns", {"INS-VEH": "up"}, ()),
    ("Q9JHL1", "NHERF2", ["Nherf2", "Slc9a3r2"], "down", {"INS-VEH": "up", "LIRA-VEH": "up"}, ()),
    ("Q9EQ21", "HEPC", ["Hamp"], "ns", {"INS-VEH": "down", "LIRA-VEH": "up"}, ()),
    ("P61458", "PCBD1", ["Pcbd1"], "ns", {"INS-VEH": "up", "LIRA-VEH": "down"}, ()),
    ("Q8K385", "FRRS1", ["Frrs1"], "ns", {"INS-VEH": "up", "LIRA-VEH": "down"}, ()),
]
MODEL_CONTRAST = "VEH-CTRL"


def resolve(acc, genes, fasta_genes):
    """Accession in the current FASTA: the paper's (canonical) accession, else by gene."""
    canonical = acc.split("-")[0]
    if canonical in fasta_genes:
        return canonical
    wanted = {g.lower() for g in genes}
    hits = [a for a, g in fasta_genes.items() if g.lower() in wanted]
    return hits[0] if len(hits) == 1 else None


def verdict(direction, fc, adj):
    if direction == "ns" or fc is None or np.isnan(fc):
        return "–"
    same = (fc > 0) == (direction == "up")
    if not same:
        return "opposite"
    return "confirmed" if adj < mq.P_CUTOFF else "same direction"


def main():
    cr = mq.load_msstats("ComparisonResult")
    label_col, fc_col, adj_col = mq.col(cr, "Label"), mq.col(cr, "log2FC"), mq.col(cr, "adj.pvalue")
    prot_col = mq.col(cr, "Protein")
    stats = {(p, lab): (fc, adj) for p, lab, fc, adj in
             zip(cr[prot_col].astype(str), cr[label_col], cr[fc_col].astype(float), cr[adj_col].astype(float))}
    pg = mq.load_protein_groups()
    pg_members = set()
    for ids in pg[mq.col(pg, "Protein IDs")].astype(str):
        pg_members.update(i.split("-")[0] for i in ids.split(";"))
    fasta_genes = mq.fasta_genes()
    labels = [c[0] for c in mq.contrasts()]

    rows, grid, sig_grid, names = [], [], [], []
    counts = {"model": [0, 0, 0], "treat": [0, 0, 0]}  # [tested, same direction, confirmed]
    n_quant = 0
    for acc, short, genes, model_dir, treat, trends in TABLE1:
        acc_now = resolve(acc, genes, fasta_genes)
        in_msstats = acc_now is not None and any((acc_now, lab) in stats for lab in labels)
        names.append(f"{short} ({fasta_genes.get(acc_now, genes[0]) if acc_now else genes[0]})")
        grid.append([stats.get((acc_now, lab), (np.nan, np.nan))[0] for lab in labels])
        sig_grid.append([stats.get((acc_now, lab), (np.nan, np.nan))[1] < mq.P_CUTOFF for lab in labels])
        if not in_msstats:
            status = "identified, not quantified by MSstats" if acc_now in pg_members else "not identified"
            rows.append((short, acc_now or acc, status, model_dir, "", "", fmt_treat(treat, trends), ""))
            continue
        n_quant += 1
        fc, adj = stats.get((acc_now, MODEL_CONTRAST), (np.nan, np.nan))
        v_model = verdict(model_dir, fc, adj)
        if model_dir != "ns" and not np.isnan(fc):
            counts["model"][0] += 1
            counts["model"][1] += v_model in ("same direction", "confirmed")
            counts["model"][2] += v_model == "confirmed"
        ours_treat, v_treat = [], []
        for lab, direction in treat.items():
            tfc, tadj = stats.get((acc_now, lab), (np.nan, np.nan))
            v = verdict(direction, tfc, tadj)
            ours_treat.append(f"{lab} {fmt_fc(tfc, tadj)}")
            v_treat.append(v)
            if not np.isnan(tfc):
                counts["treat"][0] += 1
                counts["treat"][1] += v in ("same direction", "confirmed")
                counts["treat"][2] += v == "confirmed"
        rows.append((short, acc_now, "quantified", model_dir + ("*" if MODEL_CONTRAST in trends else ""),
                     fmt_fc(fc, adj), v_model, fmt_treat(treat, trends),
                     "; ".join(f"{o}: {v}" for o, v in zip(ours_treat, v_treat)) or "–"))

    figure = plot(names, labels, np.array(grid, dtype=float), np.array(sig_grid))
    m, t = counts["model"], counts["treat"]
    text = (
        f"**{n_quant} of {len(TABLE1)}** proteins of the paper's Table 1 were quantified by MSstats. "
        f"db/db vs db/+ (VEH-CTRL): the direction agrees for **{m[1]} of {m[0]}** proteins with "
        f"a reported direction, confirmed (adj. p < 0.05) for {m[2]}. Treatment effects "
        f"(LIRA-VEH, INS-VEH): the direction agrees for **{t[1]} of {t[0]}**, confirmed for "
        f"{t[2]}. Cells show log2FC (adjusted p); * marks directions reported as trends "
        "(P < 0.1) in the paper.\n\n"
        + mq.md_table(rows, ["Protein", "Accession (current)", "Status", "Paper: db/db vs db/+",
                             "VEH-CTRL log2FC (adj. p)", "Agreement", "Paper: treatment effect",
                             "This analysis: treatment effect"])
    )
    short = (f"{n_quant}/{len(TABLE1)} quantified; db/db vs db/+ direction agrees {m[1]}/{m[0]} "
             f"(confirmed {m[2]}); treatment effects agree {t[1]}/{t[0]} (confirmed {t[2]})")
    return mq.answer(QUESTION, text, figure, short)


def fmt_fc(fc, adj):
    if fc is None or np.isnan(fc):
        return "–"
    if np.isinf(fc):
        return f"{'+' if fc > 0 else '-'}inf (missing in one group)"
    return f"{fc:+.2f} ({adj:.2g})"


def fmt_treat(treat, trends):
    if not treat:
        return "ns"
    return ", ".join(f"{lab.split('-')[0]} {d}{'*' if lab in trends else ''}" for lab, d in treat.items())


def plot(names, labels, grid, sig):
    lim = 2.0
    shown = np.clip(np.where(np.isinf(grid), np.sign(grid) * lim, grid), -lim, lim)
    cmap = LinearSegmentedColormap.from_list("div", mq.DIVERGING)
    cmap.set_bad(mq.GRID)
    fig, ax = mq.plt.subplots(figsize=(6.2, 8.2))
    im = ax.imshow(np.ma.masked_invalid(shown), cmap=cmap, vmin=-lim, vmax=lim, aspect="auto")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if not np.isnan(grid[i, j]):
                txt = f"{grid[i, j]:+.1f}" if np.isfinite(grid[i, j]) else "inf"
                ax.text(j, i, txt + ("*" if sig[i, j] else ""), ha="center", va="center",
                        fontsize=7, color=mq.INK)
    ax.set_xticks(range(len(labels)), labels)
    ax.set_yticks(range(len(names)), names)
    ax.tick_params(length=0, labelsize=8)
    ax.xaxis.tick_top()
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03, label="log2 fold change (clipped at ±2)")
    ax.set_title("Paper's Table 1 proteins in this analysis\n(* adj. p < 0.05, grey: not quantified)",
                 pad=28)
    return mq.save_fig(fig, "q07_paper_proteins")


if __name__ == "__main__":
    print(main())
