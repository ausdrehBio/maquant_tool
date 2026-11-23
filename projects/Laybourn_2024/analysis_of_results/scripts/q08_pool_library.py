"""Q8 (Part B). What changes with a library from the pooled samples?

The authors searched the samples against a Spectronaut library built from DIA searches of
pooled samples (93,131 peptides, 10,003 proteins). Laybourn_2024 uses MaxDIA's library predicted
from the FASTA files instead. Laybourn_2024_pools searches the 45 pool runs (3 pools, 10 high-pH
fractions and 5 unfractionated runs each) with the same settings, and Laybourn_2024_poollib
searches the 44 samples against its evidence.txt and msms.txt as library ("MaxQuant DDA
library"), everything else unchanged.
"""
import numpy as np
import pandas as pd

import lb_common as lb
import q03_identifications as q03
import q06_paper_claims as q06

QUESTION = "What changes with a library from the pooled samples?"
PAPER = "Spectronaut library from DIA searches of the pools: 93,131 peptides, 10,003 proteins"
NOTE = "variant of our workflow; no verdict"
VARIANT = "pool library"


def library_size():
    txt = lb.POOLS_DIR / "results" / "combined" / "txt"
    ev = pd.read_csv(txt / "evidence.txt", sep="\t", low_memory=False,
                     usecols=lambda c: c in ("Modified sequence", "Charge", "Sequence", "Decoy",
                                             "Potential contaminant"))
    bad = ev.get("Decoy", pd.Series("", index=ev.index)).fillna("").eq("+")
    bad |= ev["Potential contaminant"].fillna("").eq("+")
    ev = ev[~bad]
    return {"precursors": len(ev.drop_duplicates(["Modified sequence", "Charge"])),
            "peptides": ev.Sequence.nunique()}


def isg_parts(key):
    """The interferon parts of Q6 (core DEPs and MX2/IFI6) reproduced in a search."""
    ok, n = 0, 0
    for claim, comp, gene, direction in q06.DEP_CLAIMS:
        if "ISG" in claim or "common" in claim:
            ok += q06.dep_part(comp, gene, direction, key)[1][lb.PRIMARY_LAYER]
            n += 1
    return ok, n


def main():
    if VARIANT not in lb.available_searches():
        return lb.answer(QUESTION, "The search ../Laybourn_2024_poollib has no results.", None,
                         "not run", None, "")
    keys = ["ours", VARIANT]
    c = {k: q03.counts(k) for k in keys}
    pools = q03.pools_counts()
    lib = library_size()
    summ = {(k, layer): lb.agreement_summary(k, layer) for k in keys for layer in lb.LAYERS}
    vp = {k: lb.viral_precursors(k) for k in keys}
    isg = {k: isg_parts(k) for k in keys}
    overlap = {n: lb.dep_overlap("ours", VARIANT, lb.PRIMARY_LAYER, n) for n in lb.COMPARISON_NAMES}
    rev = {n: lb.dep_overlap(VARIANT, "ours", lb.PRIMARY_LAYER, n) for n in lb.COMPARISON_NAMES}
    r_lfq, n_matched = lb.lfq_correlation("ours", VARIANT)
    present = {k: (lb.search(k).lfq.notna().mean().mean(),
                   lb.search(k).layer("run-filtered").notna().mean().mean()) for k in keys}

    def sw_viral(k):
        sw = [x for x in lb.samples() if lb.group_of(x) == "sw"]
        return int(vp[k].loc[sw, "total"].sum()), int((vp[k].loc[sw, "total"] >= 2).sum())

    rows = [
        ("library", "predicted from the FASTA files",
         f"pool search (without contaminants): {lb.thousands(lib['peptides'])} peptides, "
         f"{lb.thousands(lib['precursors'])} precursors, {lb.thousands(pools['groups'])} "
         "protein groups"),
        ("protein groups", *[lb.thousands(c[k]["groups"]) for k in keys]),
        ("peptides", *[lb.thousands(c[k]["peptides"]) for k in keys]),
        ("protein groups detected in all 44 runs",
         *[lb.thousands(c[k]["groups_all_runs"]) for k in keys]),
        ("protein groups per run, detected (median)",
         *[lb.thousands(c[k]["groups_per_run_detected"]) for k in keys]),
        ("values present: as reported / run-filtered",
         *[f"{lb.pct(present[k][0])} / {lb.pct(present[k][1])}" for k in keys]),
        ("passing viral precursors in the 16 swH1N1 samples (sum)",
         *[sw_viral(k)[0] for k in keys]),
        ("virus-positive swH1N1 samples", *[f"{sw_viral(k)[1]}/16" for k in keys]),
        ("interferon claim parts of Q6 reproduced", *[f"{isg[k][0]}/{isg[k][1]}" for k in keys]),
    ]
    for layer in lb.LAYERS:
        rows.append((f"published DEPs in the same direction ({layer})",
                     *[lb.pct(summ[(k, layer)]["direction"]) for k in keys]))
        rows.append((f"published DEPs found as DEPs ({layer})",
                     *[lb.pct(summ[(k, layer)]["recall"]) for k in keys]))
    dep_rows = [(n, overlap[n]["a"], overlap[n]["b"], lb.pct(overlap[n]["a_in_b"], 0),
                 lb.pct(rev[n]["a_in_b"], 0)) for n in lb.COMPARISON_NAMES]
    med_overlap = np.median([overlap[n]["a_in_b"] for n in lb.COMPARISON_NAMES])
    d_rec = summ[(VARIANT, lb.PRIMARY_LAYER)]["recall"] - summ[("ours", lb.PRIMARY_LAYER)]["recall"]
    d_dir = (summ[(VARIANT, lb.PRIMARY_LAYER)]["direction"]
             - summ[("ours", lb.PRIMARY_LAYER)]["direction"])

    lines = [
        lb.md_table(rows, ["", "ours (predicted library)", "pool library"]),
        "",
        f"The pool library is smaller than the authors' ({lb.thousands(lib['peptides'])} vs "
        "93,131 peptides): MaxDIA identifies fewer peptides in the 45 pool runs than Spectronaut "
        "put into its library. Searched against it, the samples yield "
        f"{lb.thousands(c[VARIANT]['peptides'])} instead of {lb.thousands(c['ours']['peptides'])} "
        "peptides, but finds the library peptides more consistently: "
        f"{lb.thousands(c[VARIANT]['groups_all_runs'])} protein groups are detected in every run "
        f"(predicted library: {lb.thousands(c['ours']['groups_all_runs'])}; paper: 3,830 "
        "'quantified').",
        "",
        f"Agreement with Table S3 ({lb.PRIMARY_LAYER}) changes by {100 * d_dir:+.1f} percentage "
        f"points in direction and {100 * d_rec:+.1f} in recall"
        + (": the library is not what keeps the published DEP lists from being reproduced "
           "(Q5)." if abs(d_rec) < 0.03 and abs(d_dir) < 0.03 else "."),
        "",
        f"**DEPs** ({lb.PRIMARY_LAYER}): share of one search's DEPs that are DEPs of the same "
        "sign in the other:",
        "",
        lb.md_table(dep_rows, ["Comparison", "DEPs, predicted library", "DEPs, pool library",
                               "predicted found in pool", "pool found in predicted"]),
        "",
        f"A median of {lb.pct(med_overlap, 0)} of our DEPs recur with the pool library "
        f"(matched groups: {lb.thousands(n_matched)}, log2 LFQ r = {r_lfq:.3f} per sample): "
        "the library changes which proteins pass the paper pipeline's thresholds far more "
        "than it changes the agreement with Table S3.",
    ]
    short = (f"recall {100 * d_rec:+.1f} and direction {100 * d_dir:+.1f} points vs ours; "
             f"{lb.thousands(c[VARIANT]['groups_all_runs'])} groups in all runs; "
             f"{lb.pct(med_overlap, 0)} of DEPs recur")
    metrics = {"recall_variant": summ[(VARIANT, lb.PRIMARY_LAYER)]["recall"],
               "direction_variant": summ[(VARIANT, lb.PRIMARY_LAYER)]["direction"],
               "groups_all_runs": c[VARIANT]["groups_all_runs"], "dep_overlap": med_overlap,
               "library_peptides": lib["peptides"], "lfq_r": r_lfq}
    return lb.answer(QUESTION, "\n".join(lines), None, short, None, "", metrics)


if __name__ == "__main__":
    print(main())
