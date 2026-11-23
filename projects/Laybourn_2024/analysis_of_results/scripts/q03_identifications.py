"""Q3. How do the identifications compare with the paper?

The paper reports 8,935 protein groups (76,690 peptides) identified, of which 3,830 (17,526
peptides) were quantified "due to the stringent search parameters used for quantification",
and a library of 93,131 peptides and 10,003 proteins (sections 2.3.2 and 3.1). It does not say
what "quantified" means; protein groups detected in every run are the closest analogue in our
data. Also counted: how many of the proteins in Table S3 our search identifies at all.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import lb_common as lb

QUESTION = "How do the identifications compare with the paper?"
PAPER = "8,935 protein groups (76,690 peptides) identified; 3,830 (17,526 peptides) quantified"
NOTE = "Spectronaut groups proteins differently; 'quantified' is not defined in the paper"
PAPER_COUNTS = {"groups": 8935, "peptides": 76690, "quantified_groups": 3830,
                "quantified_peptides": 17526, "library_peptides": 93131, "library_proteins": 10003}


def counts(key):
    s = lb.search(key)
    ev = s.evidence
    clean_ev = ev[~ev["Potential contaminant"]]
    ok = clean_ev[clean_ev["passes"]]
    n = len(lb.samples())
    pep_runs = ok.groupby("Sequence").Experiment.nunique()
    det = s.detected
    return {
        "groups": len(s.pg),
        "peptides": s.peptide_count(),
        "groups_all_runs": int((det.sum(axis=1) == n).sum()),
        "peptides_all_runs": int((pep_runs == n).sum()),
        "groups_per_run_reported": float(s.lfq.notna().sum().median()),
        "groups_per_run_detected": float(det.sum().median()),
        "contaminants": s.n_contaminants,
    }


def pools_counts():
    txt = lb.POOLS_DIR / "results" / "combined" / "txt"
    if not (txt / "proteinGroups.txt").exists():
        return None
    pg = pd.read_csv(txt / "proteinGroups.txt", sep="\t", low_memory=False,
                     usecols=["Decoy", "Potential contaminant", "Only identified by site"])
    clean = ~(pg.fillna("").eq("+").any(axis=1))
    pep = pd.read_csv(txt / "peptides.txt", sep="\t", low_memory=False,
                      usecols=lambda c: c in ("Sequence", "Decoy", "Potential contaminant"))
    pclean = ~(pep.drop(columns="Sequence").fillna("").eq("+").any(axis=1))
    return {"groups": int(clean.sum()), "peptides": int(pclean.sum())}


def coverage(key):
    s = lb.search(key)
    pub = lb.table_s3()
    accs = pub.acc.unique()
    found = [a for a in accs if s.group_of_accession(a) is not None]
    return len(found) / len(accs), len(found), len(accs), sorted(set(accs) - set(found))


def main():
    keys = lb.available_searches()
    c = {k: counts(k) for k in keys}
    cov = {k: coverage(k) for k in keys}
    pools = pools_counts()
    ours = c["ours"]
    ratio_groups = ours["groups"] / PAPER_COUNTS["groups"]
    ratio_peptides = ours["peptides"] / PAPER_COUNTS["peptides"]
    cov_ours = cov["ours"][0]

    def col(k):
        x = c[k]
        return [lb.thousands(x["groups"]), lb.thousands(x["peptides"]),
                lb.thousands(x["groups_all_runs"]), lb.thousands(x["peptides_all_runs"]),
                lb.thousands(x["groups_per_run_reported"]),
                lb.thousands(x["groups_per_run_detected"]),
                f"{cov[k][1]}/{cov[k][2]} ({lb.pct(cov[k][0])})"]
    labels = ["protein groups identified", "peptides identified",
              "protein groups detected in all 44 runs", "peptides detected in all 44 runs",
              "protein groups per run, as reported (median)",
              "protein groups per run, detected (median)", "Table S3 proteins identified"]
    paper_col = [lb.thousands(PAPER_COUNTS["groups"]), lb.thousands(PAPER_COUNTS["peptides"]),
                 f"{lb.thousands(PAPER_COUNTS['quantified_groups'])} ('quantified')",
                 f"{lb.thousands(PAPER_COUNTS['quantified_peptides'])} ('quantified')", "", "",
                 "all"]
    table = [(lab, p, *[col(k)[i] for k in keys])
             for i, (lab, p) in enumerate(zip(labels, paper_col))]

    # missing published proteins: are they presence/absence DEPs?
    pub = lb.table_s3()
    missing = set(cov["ours"][3])
    miss_rows = pub[pub.acc.isin(missing)]
    miss_presence = (miss_rows.lfc.abs() >= lb.PRESENCE_LFC).mean()
    all_presence = (pub.lfc.abs() >= lb.PRESENCE_LFC).mean()

    # figure: per-run detected groups by sample group, as reported vs run-filtered
    s = lb.search("ours")
    fig, ax = plt.subplots(figsize=(8, 3.2))
    x = np.arange(len(lb.samples()))
    ax.bar(x, s.lfq.notna().sum()[lb.samples()], color=lb.AXIS, width=0.8, label="as reported")
    ax.bar(x, s.detected.sum()[lb.samples()], width=0.5,
           color=[lb.GROUP_COLORS[lb.group_of(smp)] for smp in lb.samples()],
           label="detected in the run (q <= 0.01, PEP <= 0.1)")
    ax.axhline(PAPER_COUNTS["quantified_groups"], color=lb.INK, lw=0.8, ls="--")
    ax.text(len(x) - 0.5, PAPER_COUNTS["quantified_groups"] + 80, "paper: 3,830 quantified",
            ha="right", fontsize=7)
    ax.set_xticks(x, [f"{smp[:2]} {lb.pig_of(smp)}" for smp in lb.samples()], rotation=90,
                  fontsize=7)
    ax.set_ylabel("protein groups")
    ax.set_title("Protein groups per run (ours)")
    lb.value_grid(ax)
    ax.legend(loc="lower right", fontsize=7)
    figure = lb.save_fig(fig, "q03_groups_per_run")

    v_groups = lb.grade(min(ratio_groups, ratio_peptides), 0.8, 0.6)
    v_cov = lb.grade(cov_ours, 0.9, 0.75)
    verdict = lb.worst(v_groups, v_cov)
    criterion = ("Consistent if our identified protein groups and peptides are each at least 80% "
                 "of the paper's and we identify at least 90% of the proteins in Table S3; "
                 "minor deviation at 60% and 75%; inconsistent below.")
    lines = [
        lb.md_table(table, ["", "Paper (Spectronaut)", *keys]),
        "",
        f"Our search identifies {lb.pct(ratio_peptides)} of the paper's peptides but only "
        f"{lb.pct(ratio_groups)} of its protein groups. With practically the same peptides, the "
        "difference most likely lies in protein inference: MaxQuant merges proteins whose "
        "peptides are all shared into one group (the pig proteome holds many unreviewed TrEMBL "
        "isoforms of the same gene), while Spectronaut appears to keep more isoform-level groups "
        "apart (Q5 shows published accessions that sit inside larger MaxQuant groups). The two "
        "counts are therefore not a measure of sensitivity.",
        "",
        f"The paper's 3,830 'quantified' protein groups compare best with the "
        f"{lb.thousands(ours['groups_all_runs'])} groups we detect in all 44 runs; the paper's "
        f"17,526 'quantified' peptides do not match our {lb.thousands(ours['peptides_all_runs'])} "
        "peptides detected in all runs, so the authors' definition remains unclear. Per run, "
        f"MaxDIA reports a value for {lb.thousands(ours['groups_per_run_reported'])} groups "
        f"(median) but detects {lb.thousands(ours['groups_per_run_detected'])} (Q1, Q2).",
        "",
        f"Of the {cov['ours'][2]} distinct proteins in Table S3, {cov['ours'][1]} are in our "
        f"protein groups. Of the published DEP entries whose protein we do not identify, "
        f"{lb.pct(miss_presence)} have |log2 FC| >= {lb.PRESENCE_LFC:.0f}, i.e. were present in "
        f"one group only in the authors' data (all entries: {lb.pct(all_presence)}).",
    ]
    if pools:
        lines += ["", f"Library: the authors' Spectronaut library from the pools holds "
                  f"{lb.thousands(PAPER_COUNTS['library_peptides'])} peptides and "
                  f"{lb.thousands(PAPER_COUNTS['library_proteins'])} proteins; our MaxDIA search "
                  f"of the 45 pool runs (../Laybourn_2024_pools, the pool library of Q8) "
                  f"identifies {lb.thousands(pools['peptides'])} peptides in "
                  f"{lb.thousands(pools['groups'])} protein groups."]
    short = (f"protein groups {lb.pct(ratio_groups)} and peptides {lb.pct(ratio_peptides)} of "
             f"the paper; {lb.pct(cov_ours)} of the Table S3 proteins identified")
    metrics = {"ratio_groups": ratio_groups, "ratio_peptides": ratio_peptides,
               "coverage": cov_ours, "groups_all_runs": ours["groups_all_runs"],
               **{f"{k}_groups": c[k]["groups"] for k in keys}}
    return lb.answer(QUESTION, "\n".join(lines), figure, short, verdict, criterion, metrics)


if __name__ == "__main__":
    print(main())
