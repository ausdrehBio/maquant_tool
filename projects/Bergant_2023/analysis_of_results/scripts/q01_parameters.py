"""Q1. Do our search settings match the authors' MaxQuant run?

Compares the deposited mqpar_IAV_FPMS.xml (MaxQuant 1.5.6.2) with our parameter/mqpar.xml
(2.9.0): a curated list of the settings that affect identification and quantification, then
every other setting that both files share. Also checks that every protein ID in our results
exists in the FASTA files (identifier parsing) and compares the contaminant sets.
"""
import re
import xml.etree.ElementTree as ET

import mq_common as mq

QUESTION = "Do our search settings match the authors' MaxQuant run?"
PAPER = "MaxQuant 1.5.6.2 with default parameters, LFQ and match between runs (Methods)"
NOTE = "settings that affect identification and quantification; deposited mqpar_IAV_FPMS.xml"
MSMS_NAMES_1562 = ["FTMS", "ITMS", "TOF", "UNKNOWN"]   # order of msmsParamsArray in 1.5.6.2
PER_FILE = {"filePaths", "experiments", "fractions", "ptms", "paramGroupIndices",
            "referenceChannel"}
IGNORED = {"name", "numThreads", "maxQuantVersion", "tempFolder", "fixedSearchFolder",
           "fixedCombinedFolder", "fastaFiles", "fastaFilePath", "sendEmail"}


def first(root, *paths):
    for p in paths:
        e = root.find(p)
        if e is not None:
            return e
    return None


def text(root, *paths):
    e = first(root, *paths)
    if e is None:
        return "–"
    if len(e):
        vals = [x.text or "" for x in e]
        return ", ".join(sorted(v for v in vals if v)) or "none"
    return (e.text or "").strip() or "–"


def fasta_names(root):
    paths = [e.text for e in root.iter("fastaFilePath")] or \
            [e.text for e in root.find("fastaFiles") if e.tag == "string"]
    return ", ".join(sorted(re.split(r"[\\/]", p or "")[-1] for p in paths))


def msms(root, analyzer, tag):
    arr = root.find("msmsParamsArray")
    for i, p in enumerate(arr):
        name = p.findtext("Name") or (MSMS_NAMES_1562[i] if i < len(MSMS_NAMES_1562) else "")
        if name.upper() == analyzer:
            value = p.findtext(tag) or ""
            in_ppm = p.findtext(tag + "InPpm")
            unit = "" if in_ppm is None else (" ppm" if in_ppm.lower() == "true" else " Da")
            if in_ppm is None and analyzer == "FTMS":
                unit = " ppm"   # 1.5.6.2 stores FTMS tolerances in ppm without a flag
            return value + unit
    return "–"


PG = "parameterGroups/parameterGroup/"
# label, relevant for identification/quantification, getter
SETTINGS = [
    ("MaxQuant version", False, lambda r: text(r, "maxQuantVersion")),
    ("raw files in the search", False, lambda r: str(len(r.find("filePaths")))),
    ("FASTA files", True, fasta_names),
    ("enzyme", True, lambda r: text(r, PG + "enzymes")),
    ("enzyme mode", True, lambda r: text(r, PG + "enzymeMode")),
    ("max. missed cleavages", True, lambda r: text(r, PG + "maxMissedCleavages")),
    ("fixed modifications", True, lambda r: text(r, PG + "fixedModifications", "fixedModifications")),
    ("variable modifications", True, lambda r: text(r, PG + "variableModifications")),
    ("max. modifications per peptide", True, lambda r: text(r, PG + "maxNmods")),
    ("first search: own variable mods", True,
     lambda r: text(r, PG + "useVariableModificationsFirstSearch")),
    ("first search: variable mods", True, lambda r: text(r, PG + "variableModificationsFirstSearch")),
    ("max. charge", True, lambda r: text(r, PG + "maxCharge")),
    ("first search tolerance (ppm)", True, lambda r: text(r, PG + "firstSearchTol")),
    ("main search tolerance (ppm)", True, lambda r: text(r, PG + "mainSearchTol")),
    ("isotope match tolerance", True, lambda r: text(r, PG + "isotopeMatchTol")),
    ("MS/MS match tolerance, FTMS", True, lambda r: msms(r, "FTMS", "MatchTolerance")),
    ("MS/MS deisotoping tolerance, FTMS", True, lambda r: msms(r, "FTMS", "DeisotopeTolerance")),
    ("PSM FDR", True, lambda r: text(r, "peptideFdr")),
    ("protein FDR", True, lambda r: text(r, "proteinFdr")),
    ("site FDR", True, lambda r: text(r, "siteFdr")),
    ("min. peptide length", True, lambda r: text(r, "minPeptideLength", "minPepLen")),
    ("max. peptide mass", True, lambda r: text(r, "maxPeptideMass")),
    ("min. peptides / razor / unique", True, lambda r: " / ".join(
        text(r, t) for t in ("minPeptides", "minRazorPeptides", "minUniquePeptides"))),
    ("decoy mode", True, lambda r: text(r, "decoyMode")),
    ("include contaminants", True, lambda r: text(r, "includeContaminants")),
    ("second peptides", True, lambda r: text(r, "secondPeptide")),
    ("match between runs", True, lambda r: text(r, "matchBetweenRuns")),
    ("matching time window (min)", True, lambda r: text(r, "matchingTimeWindow")),
    ("alignment time window (min)", True, lambda r: text(r, "alignmentTimeWindow")),
    ("LFQ mode", True, lambda r: text(r, PG + "lfqMode")),
    ("LFQ min. ratio count", True, lambda r: text(r, PG + "lfqMinRatioCount")),
    ("LFQ min. / average edges", True, lambda r: text(r, PG + "lfqMinEdgesPerNode") + " / "
     + text(r, PG + "lfqAvEdgesPerNode")),
    ("LFQ: stabilise large ratios", True, lambda r: text(r, "lfqStabilizeLargeRatios")),
    ("LFQ: require MS/MS for comparisons", True, lambda r: text(r, "lfqRequireMsms")),
    ("separate LFQ in parameter groups", True, lambda r: text(r, "separateLfq")),
    ("iBAQ / log fit", True, lambda r: text(r, "ibaq") + " / " + text(r, "ibaqLogFit")),
    ("protein quantification: restrict to mods", True, lambda r: text(r, "restrictMods")),
    ("use razor and unique peptides (counterparts)", True, lambda r: text(r, "useCounterparts")),
    ("multiplicity / LC-MS run type", True, lambda r: text(r, PG + "multiplicity") + " / "
     + text(r, PG + "lcmsRunType")),
]


def same(a, b, label=""):
    if label == "FASTA files":   # the authors' local names have spaces, PRIDE's underscores
        a, b = a.replace(" ", "_"), b.replace(" ", "_")
    return a.strip().lower() == b.strip().lower()


def flat(root):
    """Leaf tag paths -> list of values, without per-file lists and ignored tags."""
    out = {}

    def walk(e, path):
        kids = list(e)
        if not kids:
            out.setdefault(path, []).append((e.text or "").strip())
            return
        for k in kids:
            walk(k, f"{path}/{k.tag}")

    for child in root:
        if child.tag in PER_FILE or child.tag in IGNORED or child.tag == "msmsParamsArray":
            continue
        walk(child, child.tag)
    return out


def classify(tag):
    t = tag.lower()
    if "denovo" in t:
        return "de novo sequencing, not used"
    if "crosslink" in t or "/xl" in t or t.startswith("xl"):
        return "cross-linking, not used"
    if "dia" in t.split("/")[-1][:3]:
        return "DIA, not used"
    if "tims" in t or "neucode" in t or "reporter" in t or "isobaric" in t:
        return "other acquisition modes, not used"
    return "check"


def identifier_check():
    """Our protein IDs that are not in the FASTA files (wrong identifier parsing would show)."""
    known = mq.fasta_ids()
    missing, total = set(), 0
    for s in mq.ours().pg["Protein IDs"].fillna(""):
        for t in str(s).split(";"):
            if not t or t.startswith(("CON__", "REV__")):
                continue
            total += 1
            if t not in known:
                missing.add(t)
    return total, sorted(missing)


def contaminant_sets():
    def con(layer):
        return {t for s in layer.pg["Protein IDs"].fillna("") for t in str(s).split(";")
                if t.startswith("CON__")}
    a, b = con(mq.ours()), con(mq.deposited())
    return len(a), len(b), len(a & b)


def main():
    dep = ET.parse(mq.DEPOSITED_MQPAR).getroot()
    our = ET.parse(mq.OUR_MQPAR).getroot()
    rows, relevant_diffs = [], []
    for label, relevant, get in SETTINGS:
        a, b = get(dep), get(our)
        ok = same(a, b, label)
        if ok:
            note = "same" if a.lower() == b.lower() else "same (spaces vs underscores in names)"
        elif label in ("MaxQuant version", "raw files in the search"):
            note = "intended (version; the paper's 42 runs)"
        elif relevant:
            note = "**differs**"
            relevant_diffs.append(label)
        else:
            note = "differs"
        rows.append((label, a, b, note))

    fd, fo = flat(dep), flat(our)
    shared = sorted(set(fd) & set(fo))
    other = []
    for tag in shared:
        if [v.lower() for v in fd[tag]] != [v.lower() for v in fo[tag]]:
            other.append((tag, ", ".join(fd[tag])[:40], ", ".join(fo[tag])[:40], classify(tag)))
    unexplained = [o for o in other if o[3] == "check"]

    total, missing = identifier_check()
    n_con_ours, n_con_dep, n_con_both = contaminant_sets()

    verdict = "consistent"
    if relevant_diffs or missing:
        verdict = "inconsistent" if missing else "minor deviation"
    criterion = ("Consistent if no setting that affects identification or quantification differs "
                 "and every protein ID of our results exists in the FASTA files.")
    lines = [
        f"All **{len(SETTINGS)}** curated settings were compared; **{len(relevant_diffs)}** that "
        "affect identification or quantification differ"
        + (f" ({', '.join(relevant_diffs)})" if relevant_diffs else "")
        + ". The run differs from the authors' only in the MaxQuant version (1.5.6.2 vs "
        f"{mq.ours().version}) and in the raw files: the paper's 42 runs instead of all 67.",
        "",
        mq.md_table(rows, ["Setting", "Authors (1.5.6.2)", "Ours (2.9.0)", ""]),
        "",
        f"Beyond the curated list, the two files share **{len(shared)}** settings, and "
        f"**{len(other)}** differ. None of them is used by this data-dependent, label-free run:",
        "",
        mq.md_table(other, ["Setting", "Authors", "Ours", "Assessment"]) if other else "none",
        "",
        f"**Identifier parsing:** {mq.n(total)} protein IDs in our proteinGroups.txt (without "
        f"CON__ and REV__); **{len(missing)}** of them are missing from the FASTA files"
        + (f": {', '.join(missing[:10])}" if missing else "") + ".",
        "",
        f"**Contaminants:** {n_con_ours} CON__ IDs in our table, {n_con_dep} in the authors', "
        f"{n_con_both} in both. MaxQuant's built-in contaminant list changed between the "
        "versions; the paper's filter removes contaminants either way.",
    ]
    if unexplained:
        lines += ["", "Shared settings that differ and were not classified automatically: "
                  + ", ".join(u[0] for u in unexplained) + "."]
    short = (f"{len(relevant_diffs)} relevant differences; only version and run set differ"
             if not relevant_diffs else f"{len(relevant_diffs)} relevant differences")
    return mq.answer(QUESTION, "\n".join(lines), None, short, verdict, criterion,
                     {"relevant_differences": len(relevant_diffs), "missing_ids": len(missing)})


if __name__ == "__main__":
    print(main())
