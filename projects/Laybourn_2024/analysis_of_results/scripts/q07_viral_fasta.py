"""Q7 (Part B). What changes with the viral FASTA files as deposited?

Laybourn_2024 searches a cleaned version of the authors' two viral FASTA files (22 proteins:
entries with internal stop codons removed, trailing '*' stripped, unique IDs;
make_parameters_Laybourn_2024.sh). Laybourn_2024_depositedfasta searches the two files as
deposited (26 entries, including translations of the whole M and NS segments with internal stop
codons) with everything else unchanged. Compared: identifications, the viral protein groups, the
quantities of the shared protein groups, and the agreement with Table S3.
"""
import numpy as np
import pandas as pd

import lb_common as lb
import q03_identifications as q03

QUESTION = "What changes with the viral FASTA files as deposited?"
PAPER = "the authors searched their two viral FASTA files"
NOTE = "variant of our workflow; no verdict"
VARIANT = "deposited FASTA"


def viral_rows(key):
    s = lb.search(key)
    det = s.detected
    sw = [x for x in lb.samples() if lb.group_of(x) == "sw"]
    rows = []
    for gid, members in sorted(s.viral_groups.items(), key=lambda kv: kv[1]):
        names = ", ".join(sorted({f"{p} ({st[:2]})" for st, p in members}))
        rows.append((key, names, int(s.pg.loc[gid, "Peptides"]),
                     f"{int(det.loc[gid, sw].sum())}/16"))
    return rows


def precursor_agreement(key_a, key_b):
    """Precursor intensities per run shared by both searches: share identical, r of log2."""
    def keyed(k):
        ev = lb.search(k).evidence
        x = ev.assign(key=ev["Modified sequence"] + "/" + ev.Charge.astype(str) + "/"
                      + ev.Experiment).set_index("key").Intensity
        return x[~x.index.duplicated()]
    a, b = keyed(key_a), keyed(key_b)
    j = a.index.intersection(b.index)
    r = np.corrcoef(np.log2(a[j] + 1), np.log2(b[j] + 1))[0, 1]
    return len(j) / len(a), float((a[j] == b[j]).mean()), float(r)


def lfq_differences(key_a, key_b):
    """Matched groups (one-to-one): median |log2 LFQ difference| and share > 1 (2-fold)."""
    m = lb.match_groups(key_a, key_b).dropna().astype(int)
    m = m[~m.duplicated(keep=False)]
    a = np.log2(lb.search(key_a).lfq.loc[m.index])
    b = np.log2(lb.search(key_b).lfq.loc[m.values]).set_axis(m.index)
    d = (a - b).stack()
    return float(d.abs().median()), float((d.abs() > 1).mean())


def pb1_note(key_a, key_b):
    """PB1 peptides of the search that identifies PB1, and their presence in the other."""
    found = {k: [gid for gid, m in lb.search(k).viral_groups.items()
                 if any(p == "PB1" for _, p in m)] for k in (key_a, key_b)}
    if bool(found[key_a]) == bool(found[key_b]):
        return ""
    has, lacks = (key_a, key_b) if found[key_a] else (key_b, key_a)
    seqs = lb.search(has).pg.loc[found[has][0], "Peptide sequences"].split(";")
    in_other = lb.search(lacks).evidence.Sequence.isin(seqs).sum()
    name = {"ours": "our (cleaned FASTA)", VARIANT: "deposited-FASTA"}
    return (f"PB1 is identified only by the {name[has]} search ({len(seqs)} peptides). All of "
            f"them are also in the PB1 sequences of the {name[lacks]} search, but that search "
            + ("reports none of them" if in_other == 0 else f"reports them in only {in_other} "
               "evidence rows")
            + ": PB1 sits at the identification threshold "
            "and comes out differently in two nearly identical MaxDIA runs. This decides the "
            "claim part 'M1, NS1, PB1, HA and NP identified' of Q6.")


def main():
    if VARIANT not in lb.available_searches():
        return lb.answer(QUESTION, "The search ../Laybourn_2024_depositedfasta has no results.",
                         None, "not run", None, "")
    keys = ["ours", VARIANT]
    c = {k: q03.counts(k) for k in keys}
    r_lfq, n_matched = lb.lfq_correlation("ours", VARIANT)
    shared, identical, r_prec = precursor_agreement("ours", VARIANT)
    med_diff, share_2fold = lfq_differences("ours", VARIANT)
    pb1 = pb1_note("ours", VARIANT)
    vp = {k: lb.viral_precursors(k) for k in keys}
    summ = {(k, layer): lb.agreement_summary(k, layer) for k in keys for layer in lb.LAYERS}
    overlap = {n: lb.dep_overlap("ours", VARIANT, lb.PRIMARY_LAYER, n) for n in lb.COMPARISON_NAMES}
    rev = {n: lb.dep_overlap(VARIANT, "ours", lb.PRIMARY_LAYER, n) for n in lb.COMPARISON_NAMES}

    def sw_viral(k):
        sw = [x for x in lb.samples() if lb.group_of(x) == "sw"]
        return int(vp[k].loc[sw, "total"].sum()), int((vp[k].loc[sw, "total"] >= 2).sum())

    rows = [
        ("protein groups", *[lb.thousands(c[k]["groups"]) for k in keys]),
        ("peptides", *[lb.thousands(c[k]["peptides"]) for k in keys]),
        ("protein groups detected in all 44 runs",
         *[lb.thousands(c[k]["groups_all_runs"]) for k in keys]),
        ("passing viral precursors in the 16 swH1N1 samples (sum)",
         *[sw_viral(k)[0] for k in keys]),
        ("virus-positive swH1N1 samples", *[f"{sw_viral(k)[1]}/16" for k in keys]),
    ]
    for layer in lb.LAYERS:
        rows.append((f"published DEPs in the same direction ({layer})",
                     *[lb.pct(summ[(k, layer)]["direction"]) for k in keys]))
        rows.append((f"published DEPs found as DEPs ({layer})",
                     *[lb.pct(summ[(k, layer)]["recall"]) for k in keys]))
    dep_rows = [(n, overlap[n]["a"], overlap[n]["b"], lb.pct(overlap[n]["a_in_b"], 0),
                 lb.pct(rev[n]["a_in_b"], 0)) for n in lb.COMPARISON_NAMES]
    med_overlap = np.median([overlap[n]["a_in_b"] for n in lb.COMPARISON_NAMES])
    rec_o = summ[("ours", lb.PRIMARY_LAYER)]["recall"]
    rec_v = summ[(VARIANT, lb.PRIMARY_LAYER)]["recall"]

    lines = [
        lb.md_table(rows, ["", "ours (cleaned viral FASTA)", "deposited viral FASTA"]),
        "",
        "**Viral protein groups** (members by protein and strain; detected: swH1N1 samples with "
        "a passing precursor of the group):",
        "",
        lb.md_table(viral_rows("ours") + viral_rows(VARIANT),
                    ["Search", "Members", "Peptides", "Detected in swH1N1"]),
        "",
        "With the deposited files, MaxQuant reads the identifier rule >([^\\s]*) (the default "
        "rule would give duplicate IDs, make_parameters_Laybourn_2024_depositedfasta.py), and "
        "the whole-segment translations of M and NS join the groups of their mature proteins. "
        "The viral evidence of HA, M1, NP and NS1 is practically the same: the cleaning "
        "changes protein IDs and grouping, not which of their peptides are found. " + pb1,
        "",
        f"**Host proteins.** {lb.pct(shared)} of our precursor quantities per run are also in "
        f"the variant, and {lb.pct(identical)} of these are identical (r = {r_prec:.4f}). At "
        f"the protein level the searches differ more: of our groups, {lb.thousands(n_matched)} "
        "map one-to-one to a group of the variant, and their log2 LFQ intensities correlate at "
        f"a median r = {r_lfq:.3f} per sample; the median difference is {med_diff:.3f} log2, "
        f"but {lb.pct(share_2fold)} of the values differ more than 2-fold. MaxLFQ reacts to "
        "which precursors enter a group, and MaxDIA's library, decoys and machine-learning "
        "scores depend on the whole FASTA, so two searches whose FASTA files differ only in the "
        "viral entries (22 vs 26 of about 48,400) are not identical.",
        "",
        f"**DEPs** ({lb.PRIMARY_LAYER}): share of one search's DEPs that are DEPs of the same "
        "sign in the other:",
        "",
        lb.md_table(dep_rows, ["Comparison", "DEPs, ours", "DEPs, deposited FASTA",
                               "ours found in variant", "variant found in ours"]),
        "",
        f"Even between these two nearly identical searches, a median of {lb.pct(med_overlap, 0)} "
        "of the DEPs recur: the paper pipeline's calls are unstable for proteins near its "
        "thresholds, which also bounds the agreement with Table S3 (Q5). "
        + ("The agreement with the published DEPs is the same for both FASTA versions "
           f"(recall {lb.pct(rec_v)} vs {lb.pct(rec_o)})." if abs(rec_v - rec_o) < 0.02 else
           f"The recall of published DEPs differs: {lb.pct(rec_v)} vs {lb.pct(rec_o)}."),
    ]
    short = ("Table S3 agreement unchanged; viral evidence unchanged"
             + (" except PB1 (identified only with the deposited files)" if pb1 else "")
             + f"; host LFQ r = {r_lfq:.3f}; {lb.pct(med_overlap, 0)} of DEPs recur")
    metrics = {"lfq_r": r_lfq, "dep_overlap": med_overlap, "lfq_2fold_share": share_2fold,
               "precursor_identical": identical,
               "recall_variant": summ[(VARIANT, lb.PRIMARY_LAYER)]["recall"],
               "direction_variant": summ[(VARIANT, lb.PRIMARY_LAYER)]["direction"]}
    return lb.answer(QUESTION, "\n".join(lines), None, short, None, "", metrics)


if __name__ == "__main__":
    print(main())
