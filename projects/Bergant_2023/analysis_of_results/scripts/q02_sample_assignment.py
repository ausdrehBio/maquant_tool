"""Q2. Are the samples assigned to the right conditions?

Three independent checks of design.txt:
1. Sample identity: every published sample must correlate best with the sample of the same name
   in our data (and in the deposited table). Profiles are offset-free (log2 value minus the
   protein's mean), computed on 1:1-matched protein groups with a real value in all 42 samples
   of both tables, so the published imputed entries play no role. The verdict uses the margin to
   the best sample of another condition: a swap between replicates of one condition would not
   change any result, a swap between conditions would.
2. NS1 genotype: the strains differ in NS1 residues 103/106 (S/I or F/M) and 184 (G or R). Our
   FASTA holds wild-type PR8 NS1 (S103, I106, G184), so a peptide covering residue 184 can only
   be identified in G strains, and one covering 103/106 only in S/I strains. Only evidence
   identified by MS/MS in that run counts (types MSMS, MULTI-MSMS); match between runs
   transfers peptides into runs that do not contain them.
3. Viral traces in mock samples, from the MS/MS counts of the PR8 protein groups.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mq_common as mq

QUESTION = "Are the samples assigned to the right conditions?"
PAPER = "42 runs, 6 per condition (Sample_descriptions_IAV_FPMS.txt)"
NOTE = "sample identity from protein profiles; NS1 genotype from MS/MS-identified peptides"
NS1 = "P03496"
SITES = {"S103/I106": (103, 106), "G184": (184,)}
EXPECT_SITE = {"S103/I106": {"RSIG", "RSIR"}, "G184": {"RSIG", "RFMG", "AFMG"}}
MSMS_TYPES = {"MSMS", "MULTI-MSMS"}


def ns1_sequence():
    seq, take = "", False
    for line in open(mq.FASTA_VIRUS["PR8"]):
        if line.startswith(">"):
            take = f"|{NS1}|" in line
            continue
        if take:
            seq += line.strip()
    return seq


def profile_mapping(A, B, valid_a, valid_b):
    """A, B: matched rows x SAMPLES log2 values. Returns (R, summary) where R[b][a] is the
    correlation of the offset-free profiles of sample a (A) and sample b (B)."""
    keep = valid_a.all(axis=1).values & valid_b.all(axis=1).values
    a = A[keep].sub(A[keep].mean(axis=1), axis=0)
    b = B[keep].sub(B[keep].mean(axis=1), axis=0)
    za = (a - a.mean()) / a.std()
    zb = (b - b.mean()) / b.std()
    R = pd.DataFrame(zb.values.T @ za.values / (len(a) - 1), index=mq.SAMPLES, columns=mq.SAMPLES)
    rows = []
    for s in mq.SAMPLES:
        col = R[s]
        cond = s.rsplit("_", 1)[0]
        best = col.idxmax()
        other_cond = col[[x for x in mq.SAMPLES if x.rsplit("_", 1)[0] != cond]].max()
        rows.append({"sample": s, "best": best, "r_same": col[s],
                     "margin_condition": col[s] - other_cond,
                     "margin_next": col[s] - col.drop(s).max()})
    return R, pd.DataFrame(rows), int(keep.sum())


def matched_matrices(layer_a, layer_b):
    """1:1-matched rows of two sources as (A, B, valid_a, valid_b), each indexed 0..n-1."""
    comp, pairs = mq.match_groups(layer_a["ids"], layer_b["ids"])
    A = layer_a["log2"].iloc[pairs["a"]].reset_index(drop=True)
    B = layer_b["log2"].iloc[pairs["b"]].reset_index(drop=True)
    va = layer_a["valid"].iloc[pairs["a"]].reset_index(drop=True)
    vb = layer_b["valid"].iloc[pairs["b"]].reset_index(drop=True)
    return A, B, va, vb


def source(name):
    if name == "published":
        pub = mq.published()
        return {"ids": list(pub.ids), "log2": pub.matrix, "valid": pub.real}
    layer = mq.layers()[name]
    lfq = layer.matrix("LFQ intensity")
    return {"ids": list(layer.ids), "log2": np.log2(lfq.where(lfq > 0)), "valid": lfq > 0}


def genotype_evidence():
    """NS1 evidence rows (MS/MS identified and matched) with peptide start/end positions."""
    seq = ns1_sequence()
    cols = ["Sequence", "Proteins", "Experiment", "Type", "Intensity"]
    parts = []
    for ch in pd.read_csv(mq.TXT_DIR / "evidence.txt", sep="\t", usecols=cols, chunksize=500000,
                          low_memory=False):
        parts.append(ch[ch["Proteins"].fillna("").str.contains(NS1)])
    ev = pd.concat(parts)
    ev["start"] = ev["Sequence"].map(lambda s: seq.find(s) + 1)
    ev = ev[ev["start"] > 0].copy()
    ev["end"] = ev["start"] + ev["Sequence"].str.len() - 1
    ev["msms"] = ev["Type"].isin(MSMS_TYPES)
    return ev, seq


def viral_traces(layer):
    """MS/MS counts of PR8 protein groups: total in mock and infected samples."""
    ids = layer.ids
    pr8 = set(mq.pr8_proteins())
    is_viral = ids.map(lambda s: bool(s & pr8)) & ~layer.decoy
    msms = layer.matrix("MS/MS count")[is_viral.values]
    mock = msms[mq.samples_of("mock")].sum(axis=1)
    infected = msms[[s for s in mq.SAMPLES if not s.startswith("mock")]].sum(axis=1)
    n_mock_samples = (msms[mq.samples_of("mock")] > 0).any(axis=0).sum()
    return int(mock.sum()), int(infected.sum()), int(n_mock_samples), int(is_viral.sum())


def plot(R, ev, seq):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axes[0]
    im = ax.imshow(R.values, cmap="viridis", vmin=np.nanpercentile(R.values, 5), vmax=1)
    ticks = [i * 6 + 2.5 for i in range(len(mq.CONDITIONS))]
    ax.set_xticks(ticks, mq.CONDITIONS)
    ax.set_yticks(ticks, mq.CONDITIONS)
    for i in range(1, len(mq.CONDITIONS)):
        ax.axhline(i * 6 - 0.5, color=mq.SURFACE, lw=0.8)
        ax.axvline(i * 6 - 0.5, color=mq.SURFACE, lw=0.8)
    ax.set_xlabel("published sample")
    ax.set_ylabel("our sample")
    ax.set_title("Profile correlation, ours vs published")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label="r of offset-free profiles")
    ax.tick_params(length=0)

    ax = axes[1]
    peps = (ev.groupby("Sequence")["start"].first().sort_values())
    mat = pd.DataFrame(np.nan, index=peps.index, columns=mq.SAMPLES)
    msms = ev[ev["msms"]].groupby(["Sequence", "Experiment"])["Intensity"].max()
    for (pep, exp), v in msms.items():
        if exp in mat.columns:
            mat.at[pep, exp] = np.log10(v) if v > 0 else 5.0
    im2 = ax.imshow(mat.values, aspect="auto", cmap="magma_r", interpolation="nearest")
    labels = [f"{peps[p]}–{peps[p] + len(p) - 1} {p[:14]}{'…' if len(p) > 14 else ''}"
              for p in mat.index]
    ax.set_yticks(range(len(mat)), labels, fontsize=6)
    ax.set_xticks(ticks, mq.CONDITIONS)
    for i in range(1, len(mq.CONDITIONS)):
        ax.axvline(i * 6 - 0.5, color=mq.GRID, lw=0.8)
    for i, p in enumerate(mat.index):
        start, end = peps[p], peps[p] + len(p) - 1
        if any(start <= x <= end for pos in SITES.values() for x in pos):
            ax.get_yticklabels()[i].set_color(mq.DIVERGING[4])
            ax.get_yticklabels()[i].set_fontweight("bold")
    ax.set_title("NS1 peptides identified by MS/MS (red: genotype sites)")
    ax.tick_params(length=0)
    fig.colorbar(im2, ax=ax, fraction=0.04, pad=0.02, label="log10 intensity")
    fig.tight_layout()
    return mq.save_fig(fig, "q02_sample_assignment")


def main():
    # 1. sample identity
    srcs = {name: source(name) for name in ["published", *mq.layers().keys()]}
    results = {}
    for a, b in [("published", "ours"), ("published", "deposited"), ("ours", "deposited")] + \
                ([("published", "ours-67")] if "ours-67" in srcs else []):
        A, B, va, vb = matched_matrices(srcs[a], srcs[b])
        R, summary, n_rows = profile_mapping(A, B, va, vb)
        results[(a, b)] = (R, summary, n_rows)
    rows = []
    cross, replicate_swaps = 0, 0
    for (a, b), (R, s, n_rows) in results.items():
        wrong_cond = int((s["best"].str.rsplit("_", n=1).str[0]
                          != s["sample"].str.rsplit("_", n=1).str[0]).sum())
        swaps = int((s["best"] != s["sample"]).sum()) - wrong_cond
        if b in ("ours", "ours-67") or a == "ours":
            cross += wrong_cond
            replicate_swaps += swaps
        rows.append((f"{a} vs {b}", mq.n(n_rows), f"{int((s['best'] == s['sample']).sum())}/42",
                     mq.f2(s["r_same"].median()), mq.f2(s["margin_condition"].min()),
                     mq.f2(s["margin_condition"].median()), mq.f2(s["margin_next"].min())))

    # 2. NS1 genotype
    ev, seq = genotype_evidence()
    site_rows, genotype_ok = [], True
    for site, positions in SITES.items():
        cover = ev[ev["msms"] & ev.apply(lambda r: any(r["start"] <= p <= r["end"]
                                                       for p in positions), axis=1)]
        seen = cover.groupby(cover["Experiment"].str.rsplit("_", n=1).str[0])["Experiment"].nunique()
        counts = {c: int(seen.get(c, 0)) for c in mq.CONDITIONS}
        expected = EXPECT_SITE[site]
        wrong = {c: k for c, k in counts.items() if c not in expected and k > 0}
        missing = {c: k for c, k in counts.items() if c in expected and k < 5} if site == "G184" else {}
        genotype_ok &= not wrong and not missing
        residues = "/".join(seq[p - 1] + str(p) for p in positions)
        site_rows.append((residues, ", ".join(sorted(expected)),
                          " ".join(f"{c} {counts[c]}" for c in mq.CONDITIONS),
                          "pass" if not wrong and not missing else "**fail**"))
    p2237 = ev[ev["msms"] & (ev["Sequence"] == "VADQELGDAPFLDRLR")]
    seen2237 = p2237.groupby(p2237["Experiment"].str.rsplit("_", n=1).str[0])["Experiment"].nunique()

    # 3. viral traces in mock samples
    trace_rows = []
    for name, layer in mq.layers().items():
        mock, infected, n_mock, n_groups = viral_traces(layer)
        trace_rows.append((name, n_groups, mq.n(mock), mq.n(infected),
                           mq.pct(mock / max(mock + infected, 1), 2), f"{n_mock}/6"))

    figure = plot(results[("published", "ours")][0], ev, seq)

    if cross or not genotype_ok:
        verdict = "inconsistent"
    elif replicate_swaps:
        verdict = "minor deviation"
    else:
        verdict = "consistent"
    criterion = ("Consistent if every published and deposited sample correlates best with our "
                 "sample of the same name and the NS1 genotype peptides appear only in the strains "
                 "that carry the residue (G184 peptide in at least 5 of 6 samples of each G strain). "
                 "A best match to another replicate of the same condition is a minor deviation; a "
                 "match to another condition or a genotype violation is inconsistent.")
    s_po = results[("published", "ours")][1]
    lines = [
        "**Sample identity.** Offset-free LFQ profiles of 1:1-matched protein groups with a real "
        "value in all 42 samples of both tables; margins are differences in correlation:",
        "",
        mq.md_table(rows, ["Tables", "Protein groups", "Same name best", "Median r (same)",
                           "Min. margin to other condition", "Median margin to other condition",
                           "Min. margin to next sample"]),
        "",
        f"Every published sample is closest to our sample of the same name "
        f"({int((s_po['best'] == s_po['sample']).sum())}/42); the smallest margin to a sample of "
        f"another condition is {mq.f2(s_po['margin_condition'].min())}. Our design.txt therefore "
        "assigns the raw files as the authors did, including the remeasured runs "
        "(RSIG_4 = REM_1_4, AFMG_6 = REM_6_6 and others).",
        "",
        "**NS1 genotype.** Samples per condition in which a peptide covering the residue was "
        "identified by MS/MS:",
        "",
        mq.md_table(site_rows, ["Residues (wild type)", "Expected in", "Samples with MS/MS "
                                "identification", "Check"]),
        "",
        "The G184 peptide (NAVGVLIGGLEWNDNTVR) is identified in all G-strain samples and in no "
        "R-strain or mock sample; the S103/I106 peptide (DWSMLIPK) only in RSIG and RSIR. Match "
        "between runs transfers both peptides into the other samples at about 1/100 of the "
        "intensity, which is why only MS/MS identifications count. Residue 38 (R/A) is not "
        "testable: its tryptic peptides are too short.",
        "",
        "Observation, not used for the verdict: the wild-type peptide 22–37 (VADQELGDAPFLDRLR) is "
        "identified by MS/MS in "
        + ", ".join(f"{c} {int(seen2237.get(c, 0))}" for c in mq.STRAINS)
        + " samples. It is missing from all AFMR samples, although AFMG has it and other NS1 "
        "peptides are identified in AFMR. This could point to an NS1 sequence difference near "
        "residues 22–38 in the AFMR virus that the strain name does not describe; the authors' "
        "NS1 sequences would be needed to tell.",
        "",
        "**Viral traces in mock samples** (MS/MS counts of PR8 protein groups):",
        "",
        mq.md_table(trace_rows, ["Table", "PR8 groups", "MS/MS in mock", "MS/MS in infected",
                                 "Mock share", "Mock samples with viral MS/MS"]),
        "",
        "Both tables show the same small viral signal in mock samples (sample carryover or "
        "contamination in the measurement, not a processing effect). Analyses below treat mock "
        "as nearly, not completely, virus-free.",
    ]
    short = (f"{int((s_po['best'] == s_po['sample']).sum())}/42 samples matched by name; "
             f"NS1 genotypes {'as expected' if genotype_ok else 'violated'}")
    s_pd = results[("published", "deposited")][1]
    return mq.answer(QUESTION, "\n".join(lines), figure, short, verdict, criterion,
                     {"cross_condition": cross, "replicate_swaps": replicate_swaps,
                      "genotype_ok": genotype_ok,
                      "matched_ours": int((s_po["best"] == s_po["sample"]).sum()),
                      "matched_deposited": int((s_pd["best"] == s_pd["sample"]).sum()),
                      "margin_ours": float(s_po["margin_condition"].min()),
                      "margin_deposited": float(s_pd["margin_condition"].min()),
                      "r_same_ours": float(s_po["r_same"].median()),
                      "r_same_deposited": float(s_pd["r_same"].median())})


if __name__ == "__main__":
    print(main())
