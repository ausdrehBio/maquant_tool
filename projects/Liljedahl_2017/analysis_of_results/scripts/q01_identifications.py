"""Q1. How many protein groups and peptides were identified?

Protein groups = rows of proteinGroups.txt, split into target groups, contaminant groups with
a mouse protein, contaminant-only groups and decoys; peptides = rows of peptides.txt.
Per run: share of identified MS/MS spectra and identified peptide sequences (summary.txt).
"""
import mq_common as mq

QUESTION = "How many protein groups and peptides were identified?"
PAPER = "1390 proteins (Progenesis QI and Mascot, 5% FDR)"
NOTE = "MaxQuant with 1% PSM and protein FDR; protein groups of the clean set"


def main():
    pg = mq.load_protein_groups()
    ids = pg[mq.col(pg, "Protein IDs")]
    decoy = mq.is_decoy(pg)
    con = mq.is_contaminant(pg) & ~decoy
    target = mq.has_target(ids) & ~decoy
    clean = mq.clean_mask(pg)
    site_only = target & ~clean
    rows_pg = [
        ("Target (mouse FASTA)", mq.n((target & ~con).sum())),
        ("Contaminant with mouse protein", mq.n((con & target).sum())),
        ("Contaminant only (CON__)", mq.n((con & ~target).sum())),
        ("Decoy (REV__)", mq.n(decoy.sum())),
        ("**Total**", f"**{mq.n(len(pg))}**"),
    ]

    pep = mq.load_peptides()
    pep_decoy = mq.is_decoy(pep)
    pep_con = mq.is_contaminant(pep) & ~pep_decoy
    n_pep_target = int((~pep_decoy & ~pep_con).sum())

    s = mq.load_summary()
    rate = s[mq.col(s, "MS/MS identified [%]")].astype(float)
    seqs = s[mq.col(s, "Peptide sequences identified")].astype(float)
    by_cond = (
        s.assign(cond=[mq.condition_of(e) for e in s.index], rate=rate, seqs=seqs)
        .groupby("cond")[["rate", "seqs"]].mean().loc[mq.CONDITIONS]
    )

    fig, axes = mq.plt.subplots(2, 1, figsize=(10, 6.4), sharex=True)
    mq.run_bars(axes[0], rate, "MS/MS identified [%]")
    axes[0].set_title("Identified MS/MS spectra per run")
    mq.run_bars(axes[1], seqs, "peptide sequences")
    axes[1].set_title("Identified peptide sequences per run")
    mq.condition_legend(axes[0], loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    fig.tight_layout()
    figure = mq.save_fig(fig, "q01_identifications")

    cond_rows = [
        (c, mq.CONDITION_LABELS[c], mq.pct(by_cond.loc[c, "rate"]), mq.n(round(by_cond.loc[c, "seqs"])))
        for c in mq.CONDITIONS
    ]
    text = (
        f"**{mq.n(len(pg))} protein groups** and **{mq.n(len(pep))} peptide sequences** were "
        f"found (rows of proteinGroups.txt and peptides.txt; {mq.n(n_pep_target)} peptides "
        "without decoys and contaminants).\n\n"
        + mq.md_table(rows_pg, ["Category", "Protein groups"])
        + f"\n\nThe **clean set** (target member, no decoy, not only identified by site) has "
        f"**{mq.n(clean.sum())} protein groups** ({mq.n(site_only.sum())} removed as only "
        "identified by site). It is used for all quantitative questions. Contaminant groups "
        "with a mouse protein stay in it.\n\n"
        f"Per run, {mq.pct(rate.min())} to {mq.pct(rate.max())} of the MS/MS spectra were "
        f"identified (median {mq.pct(rate.median())}), and {mq.n(seqs.min())} to "
        f"{mq.n(seqs.max())} peptide sequences (median {mq.n(seqs.median())}).\n\n"
        + mq.md_table(cond_rows, ["Group", "Mice", "MS/MS identified (mean)", "Peptide sequences (mean)"])
        + "\n\nThe paper reports 1390 proteins from Progenesis QI and Mascot at 5% FDR against "
        "UniProt mouse 2013.05. Here the FDR is 1% and the database is the current Swiss-Prot "
        "canonical mouse proteome, so the numbers are related but not equal."
    )
    short = f"{mq.n(clean.sum())} protein groups (clean set), {mq.n(len(pg))} in total"
    return mq.answer(QUESTION, text, figure, short)


if __name__ == "__main__":
    print(main())
