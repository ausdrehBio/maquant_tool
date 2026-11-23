"""Q1. Do our search settings match the authors' Spectronaut analysis?

The authors deposited no parameter file; their settings are those stated in the paper (sections
2.3.2 and 2.3.3) and, where the paper says "default settings", Spectronaut's defaults. Ours are
read from parameter/mqpar.xml and results/combined/txt/parameters.txt.
"""
import xml.etree.ElementTree as ET

import lb_common as lb

QUESTION = "Do our search settings match the authors' Spectronaut analysis?"
PAPER = "Spectronaut 17.4, library from DIA searches of pools; settings in sections 2.3.2-2.3.3"
NOTE = "no parameter file deposited; Spectronaut defaults where the paper says default"


def strings(group, tag):
    node = group.find(tag)
    return [] if node is None else [s.text for s in node.findall("string")]


def fasta_count(path):
    with open(path) as fh:
        return sum(line.startswith(">") for line in fh)


def main():
    root = ET.parse(lb.OUR_MQPAR).getroot()
    grp = root.find("parameterGroups/parameterGroup")
    txt = lb.search("ours").txt_dir
    version = lb.maxquant_version(txt)
    fastas = [e.findtext("fastaFilePath") for e in root.find("fastaFiles")]
    n_pig = fasta_count(lb.FASTA_PIG)
    n_vir = fasta_count(lb.FASTA_VIRAL_CLEAN)
    n_dep = [fasta_count(p) for p in lb.FASTA_VIRAL_DEPOSITED]
    g = grp.findtext
    r = root.findtext
    lib_type = {"0": "predicted from the FASTA files", "2": "MaxQuant search results"}
    # (setting, authors, ours, class); class: same, intended, software, default-differs, differs
    rows = [
        ("software", "Spectronaut 17.4", f"MaxQuant {version} ({g('lcmsRunType')})", "intended"),
        ("spectral library", "DIA searches of the pooled samples (93,131 peptides, 10,003 "
         "proteins)", lib_type.get(g("diaLibraryType"), g("diaLibraryType")),
         "intended (Q8 uses a library from the pools)"),
        ("pig FASTA", "UniProt UP000008227, 2022-04, 48,349 entries",
         f"{lb.FASTA_PIG.name} ({n_pig:,} entries)", "same"),
        ("viral FASTA", f"two in-house files ({n_dep[0]} + {n_dep[1]} entries)",
         f"the same proteins, cleaned ({n_vir} entries; parameter/viral_proteins.fasta)",
         "intended (Q7 uses the files as deposited)"),
        ("enzyme", "Trypsin/P and LysC", ", ".join(strings(grp, "enzymes")) +
         f" (mode {g('enzymeMode')}, specific)", "same (LysC sites are a subset of Trypsin/P)"),
        ("missed cleavages", "2 (Spectronaut default)",
         f"{g('maxMissedCleavages')} (library: {g('maxMissedCleavagesDia')})", "same"),
        ("fixed modifications", "Carbamidomethyl (C)",
         ", ".join(strings(grp, "fixedModifications")), "same"),
        ("variable modifications", "Acetyl (Protein N-term), Oxidation (M), Deamidation (N)",
         ", ".join(sorted(strings(grp, "variableModifications"))), "same"),
        ("variable modifications per peptide", "5 (Spectronaut default; not stated)",
         f"{g('diaMaxModifications')} in the library (maxNmods {g('maxNmods')})",
         "default differs"),
        ("peptide length", "7-52 (Spectronaut default)",
         f"min. {r('minPeptideLength')}, max. mass {r('maxPeptideMass')} Da", "software"),
        ("precursor charge", "from the library (DIA data)",
         f"{g('diaMinCharge')}-{g('diaMaxCharge')}", "software"),
        ("precursor / peptide FDR", "q-value 0.01 (default), precursor PEP cutoff 0.1",
         f"peptide FDR {r('peptideFdr')}", "same"),
        ("protein FDR", "q-value 0.01 per experiment (default), 0.05 per run, PEP cutoff 0.5",
         f"protein FDR {r('proteinFdr')}", "same at the experiment level"),
        ("values in runs where a precursor was not identified",
         "removed by the run-wise q-value filter; no imputation",
         f"transferred while the run's q-value is below {g('diaTransferQvalue')} "
         "(diaTransferQvalue)", "software (run-filtered layer, see Q2)"),
        ("protein quantification", "Spectronaut's MaxLFQ, precursor quantification on MS level",
         f"MaxLFQ (lfqMode {g('lfqMode')}, min. ratio count {g('lfqMinRatioCountDia')})",
         "software"),
        ("normalisation", "cross-run normalisation (default)",
         f"MaxLFQ normalisation (lfqNormType {g('lfqNormType')})", "software"),
        ("contaminants", "no contaminant database (not stated)",
         f"MaxQuant contaminants searched ({r('includeContaminants')}), removed before the "
         "analysis", "software"),
    ]
    stated_differs = [x for x in rows if x[3] == "differs"]
    default_differs = [x for x in rows if x[3] == "default differs"]
    if stated_differs:
        verdict = "inconsistent"
    elif default_differs:
        verdict = "minor deviation"
    else:
        verdict = "consistent"
    criterion = ("Consistent if every setting that the paper states and that defines the search "
                 "space (FASTA, enzyme, modifications) is the same and all other differences are "
                 "inherent to the software or intended; minor deviation if a Spectronaut default "
                 "that the paper does not state differs in a way that changes the search space; "
                 "inconsistent if a stated setting differs.")
    lines = [
        f"All {len(rows)} settings that the paper states or that Spectronaut sets by default were "
        f"compared with ours. **{len(stated_differs)}** stated settings differ; "
        f"**{len(default_differs)}** unstated default differs: "
        + (", ".join(x[0] for x in default_differs) if default_differs else "none") + ".",
        "",
        lb.md_table([x[:3] + (x[3],) for x in rows],
                    ["Setting", "Authors (paper; Spectronaut default)", "Ours", ""]),
        "",
        "The search space is the authors': the same pig FASTA file, the same viral proteins, the "
        "same enzyme specificity and the same modifications. MaxDIA limits library peptides to "
        f"{g('diaMaxModifications')} variable modifications, where Spectronaut's default allows "
        "5; peptides with three or more of acetylation, oxidation and deamidation are rare, so "
        "the effect on identifications is small.",
        "",
        "The difference that matters for the comparison is how values enter the protein table. "
        "Spectronaut, as the authors ran it, keeps a precursor's quantity in a run only where it "
        "passes the run-wise q-value filter, and the authors imputed nothing, so a protein not "
        "identified in a run has no value there. MaxDIA transfers each precursor that passes the "
        "experiment-wide FDR into every run where its q-value in that run is below "
        f"{g('diaTransferQvalue')}, so nearly every protein has a value in nearly every run (Q2). "
        "The questions below therefore analyse two layers of our search: **as reported** and "
        "**run-filtered**, the latter keeping a value only where the protein group has a "
        f"precursor with q-value <= {lb.Q_RUN} and PEP <= {lb.PEP_RUN} in that run, as an "
        "approximation of the authors' processing.",
    ]
    short = (f"{len(stated_differs)} stated settings differ; unstated default differs: "
             + (", ".join(x[0] for x in default_differs) if default_differs else "none"))
    metrics = {"stated_differs": len(stated_differs), "default_differs": len(default_differs),
               "pig_fasta_entries": n_pig, "viral_entries": n_vir}
    return lb.answer(QUESTION, "\n".join(lines), None, short, verdict, criterion, metrics)


if __name__ == "__main__":
    print(main())
