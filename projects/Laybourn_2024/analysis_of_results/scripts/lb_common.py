"""Shared code for the Laybourn_2024 comparison suite (q01 to q08).

Question of the suite: does our MaxQuant 2.9.0 (MaxDIA) workflow reproduce the published global
proteomics results of Laybourn et al. (2024, Front Immunol 15:1432743, doi
10.3389/fimmu.2024.1432743)? Pig tracheal tissue 3 days after infection with a swine-adapted
(swH1N1) or human-adapted (huH1N1) H1N1 virus or mock; upper (UT) and lower trachea (LT) of 22
pigs (6 control, 8 swH1N1, 8 huH1N1): 44 DIA runs, PRIDE PXD052229. The authors used Spectronaut
17.4 with a library from DIA searches of pooled samples; no MaxQuant output is deposited.

Searches (read only):
- ours:            results/combined/txt/ (MaxDIA, library predicted from the FASTA files, the
                   cleaned viral proteins of parameter/viral_proteins.fasta)
- deposited FASTA: ../Laybourn_2024_depositedfasta/ (the authors' two viral FASTA files as
                   deposited instead of the cleaned ones; everything else the same)
- pool library:    ../Laybourn_2024_poollib/ (library from a MaxDIA search of the 45 pool runs,
                   ../Laybourn_2024_pools/, as the authors built theirs from pooled samples)
- published:       data/deposited/Laybourn_2024_Table_S3.xlsx (Supplementary Table S3: the
                   differentially expressed proteins, DEPs, of six comparisons with log2 fold
                   change and adjusted p)
Further inputs: parameter/design.txt, parameter/mqpar.xml, parameter/deposited/ (GEO GSE268254
series matrix: pig -> group), the FASTA files, annotation/ (GO:0051607 set from UniProt, written
by make_annotation.py). Outputs: analysis_of_results/answers.md and verdicts.tsv (in git),
figures/ and tables/ (not in git). Paths are printed relative to the project folder.

Definitions used by all scripts:
- samples: the experiments of design.txt, <tissue>_<group>_p<pig>, e.g. LT_sw_p09.
- comparisons: the six sheets of Table S3, named as there: "UT sw vs ctr" = mean(sw) - mean(ctr)
  of the upper trachea, likewise hu vs ctr and sw vs hu, and the same for LT.
- clean protein groups: not decoy, not potential contaminant, not only identified by site.
- run-level detection: MaxDIA reports a quantity for a precursor in every run once it passes
  the experiment-wide FDR, also in runs where it was not found (the viral nucleoprotein has an
  LFQ intensity in all 12 controls). Its evidence.txt holds a q-value and PEP per precursor and
  run. A protein group counts as detected in a run if one of its precursors has q-value <= 0.01
  (Spectronaut's default precursor q-value cutoff) and PEP <= 0.1 (the authors' precursor PEP
  cutoff) in that run.
- layers of a search: "as reported" (LFQ intensity > 0) and "run-filtered" (LFQ intensity where
  the group is detected in that run, missing elsewhere). The run-filtered layer approximates the
  authors' Spectronaut processing (run-wise q-value filtering, "no imputation method was used").
- paper pipeline (section 2.7 of the paper; details reconstructed from Table S3 in q04): per
  tissue, all protein groups with a value in one of the 22 samples; log2; missing values set to
  0; one-way ANOVA per protein followed by Tukey's HSD; DEP: Tukey-adjusted p < 0.05 and
  |log2 fold change| >= 1, the fold change being the difference of the group means. No
  correction across proteins.
- matching: a published accession maps to the clean protein group of ours that contains it
  (preferring one where it is a majority protein); viral IDs of the deposited FASTA
  (ssi|MT666905|A_swine_..._NP) are matched by their accession token.
- verdicts: consistent, minor deviation, inconsistent; every question states its rule. Part B
  (q07, q08) compares variants of our workflow and carries no verdict.
"""
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import openpyxl  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

SCRIPTS_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPTS_DIR.parent
PROJECT_DIR = ANALYSIS_DIR.parent
PROJECTS_DIR = PROJECT_DIR.parent
OUT_DIR = Path(os.environ.get("LB_OUT_DIR", ANALYSIS_DIR))
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"
ANNOTATION_DIR = ANALYSIS_DIR / "annotation"
DESIGN = PROJECT_DIR / "parameter" / "design.txt"
OUR_MQPAR = PROJECT_DIR / "parameter" / "mqpar.xml"
GEO_MATRIX = PROJECT_DIR / "parameter" / "deposited" / "GSE268254_series_matrix.txt.gz"
TABLE_S3 = PROJECT_DIR / "data" / "deposited" / "Laybourn_2024_Table_S3.xlsx"
FASTA_PIG = PROJECT_DIR / "data" / "uniprot-proteom_3AUP000008227_reviewed_3Ano.fasta"
FASTA_VIRAL_CLEAN = PROJECT_DIR / "parameter" / "viral_proteins.fasta"
FASTA_VIRAL_DEPOSITED = [PROJECT_DIR / "data" / "A_Denmark_238_2020_aa.fasta",
                         PROJECT_DIR / "data" / "A_swine_Denmark_10298-4_2017_aa.fasta"]
GO_TABLE = ANNOTATION_DIR / "go0051607_defense_response_to_virus.tsv"
POOLS_DIR = PROJECTS_DIR / "Laybourn_2024_pools"
SEARCH_DIRS = {
    "ours": PROJECT_DIR,
    "deposited FASTA": PROJECTS_DIR / "Laybourn_2024_depositedfasta",
    "pool library": PROJECTS_DIR / "Laybourn_2024_poollib",
}
LAYERS = ["as reported", "run-filtered"]
PRIMARY_LAYER = "run-filtered"

Q_RUN = 0.01            # precursor q-value per run (Spectronaut default)
PEP_RUN = 0.1           # precursor PEP per run (paper, section 2.3.3)
FC_CUTOFF = 1.0         # paper: log2 fold change of +-1
P_CUTOFF = 0.05         # paper: adjusted p below 0.05
PRESENCE_LFC = 4.0      # published |log2 FC| >= 4 (16-fold): driven by missing values set to 0
PAPER_REF = ("Laybourn HA, Hellemann Polhaus C, Kristensen C, et al. (2024). Multi-omics analysis "
             "reveals the impact of influenza A virus host adaptation on immune signatures in pig "
             "tracheal tissue. Front Immunol 15:1432743, doi:10.3389/fimmu.2024.1432743")
TISSUES = ["UT", "LT"]
GROUPS = ["ctr", "sw", "hu"]
GROUP_LABELS = {"ctr": "control (mock)", "sw": "swH1N1 (swine-adapted)",
                "hu": "huH1N1 (human-adapted)"}
TISSUE_LABELS = {"UT": "upper trachea", "LT": "lower trachea"}
COMPARISONS = [(t, a, b) for t in TISSUES for a, b in (("sw", "ctr"), ("hu", "ctr"), ("sw", "hu"))]
COMPARISON_NAMES = [f"{t} {a} vs {b}" for t, a, b in COMPARISONS]
VERDICTS = ["consistent", "minor deviation", "inconsistent"]

# style (the same palette as the Bergant_2023 suite)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514c"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
GROUP_COLORS = {"ctr": MUTED, "sw": "#c23b3a", "hu": "#2a78d6"}
SEARCH_COLORS = {"ours": SERIES[0], "deposited FASTA": SERIES[2], "pool library": SERIES[3],
                 "published": INK}
OUTCOME_COLORS = {"DEP, same direction": "#1baf7a", "same direction, not DEP": "#a8dcc6",
                  "opposite direction": "#c23b3a", "not quantified": AXIS}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "figure.dpi": 110,
    "savefig.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
})


def value_grid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def save_fig(fig, name):
    """Write figures/<name>.png; return its path relative to the answers file."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"figure written: {rel(path)}")
    return path.relative_to(OUT_DIR).as_posix()


def save_table(df, name, **kw):
    """Write tables/<name>.tsv (not in git); return its path relative to the answers file."""
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / f"{name}.tsv"
    df.to_csv(path, sep="\t", **kw)
    print(f"table written: {rel(path)}")
    return path.relative_to(OUT_DIR).as_posix()


def rel(path):
    """Path relative to the project folder for printing; outside it, relative to projects/."""
    p = Path(path).resolve()
    for base, prefix in ((PROJECT_DIR, ""), (PROJECTS_DIR, "../")):
        try:
            return prefix + p.relative_to(base).as_posix()
        except ValueError:
            pass
    return p.name


def md_table(rows, header):
    def cell(c):
        return str(c).replace("|", "\\|")
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(cell(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def pct(x, digits=1):
    return "n/a" if x is None or not np.isfinite(x) else f"{100 * x:.{digits}f}%"


def f2(x, digits=2):
    return "n/a" if x is None or not np.isfinite(x) else f"{x:.{digits}f}"


def thousands(n):
    return f"{int(n):,}"


@dataclass
class Answer:
    """str() gives the Markdown section body. `short` and `verdict` feed the summary tables;
    Part B answers have verdict None."""
    question: str
    text: str
    figure: object
    short: str
    verdict: object
    criterion: str
    metrics: dict = field(default_factory=dict)

    def __str__(self):
        figs = self.figure if isinstance(self.figure, (list, tuple)) else [self.figure]
        head = (f"**Verdict: {self.verdict}.** {self.criterion}\n\n" if self.verdict
                else (f"{self.criterion}\n\n" if self.criterion else ""))
        body = f"{head}{self.text}\n"
        for fig in figs:
            if fig:
                body += f"\n![{self.question}]({fig})\n"
        return body


def answer(question, text, figure, short, verdict, criterion, metrics=None):
    if verdict is not None and verdict not in VERDICTS:
        raise ValueError(f"unknown verdict {verdict!r}")
    return Answer(question, text, figure, short, verdict, criterion, metrics or {})


def grade(value, consistent, minor, higher_is_better=True):
    """Verdict of one metric against two thresholds."""
    if not higher_is_better:
        value, consistent, minor = -value, -consistent, -minor
    if value >= consistent:
        return "consistent"
    return "minor deviation" if value >= minor else "inconsistent"


def worst(*verdicts):
    return max(verdicts, key=VERDICTS.index)


# ---------------------------------------------------------------- design and annotation
@cache
def design():
    d = pd.read_csv(DESIGN, sep="\t", dtype=str)
    d["tissue"] = d["Experiment"].str.split("_").str[0]
    d["group"] = d["Experiment"].str.split("_").str[1]
    d["pig"] = d["Experiment"].str.rsplit("_p", n=1).str[1].astype(int)
    d["run"] = d["Name"].str.extract(r"_MS(\d+)$", expand=False).astype(int)
    order = {(t, g): i for i, (t, g) in enumerate((t, g) for t in TISSUES for g in GROUPS)}
    d["order"] = [order[(t, g)] * 100 + p for t, g, p in zip(d.tissue, d.group, d.pig)]
    return d.sort_values("order").reset_index(drop=True)


@cache
def samples():
    return list(design()["Experiment"])


def samples_of(tissue, group=None):
    d = design()
    sel = d.tissue == tissue if group is None else (d.tissue == tissue) & (d.group == group)
    return list(d.loc[sel, "Experiment"])


def group_of(sample):
    return sample.split("_")[1]


def pig_of(sample):
    return int(sample.rsplit("_p", 1)[1])


@cache
def geo_groups():
    """pig -> group from the GEO series matrix (Sample_title "Pig 7, swH1N1")."""
    import gzip
    with gzip.open(GEO_MATRIX, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
                break
    names = {"control": "ctr", "swH1N1": "sw", "huH1N1": "hu"}
    out = {}
    for t in titles:
        m = re.match(r"Pig (\d+), (\S+)", t)
        out[int(m.group(1))] = names[m.group(2)]
    return out


@cache
def pig_genes():
    """UniProt accession -> gene symbol (GN=) of the pig FASTA; entry count."""
    genes = {}
    with open(FASTA_PIG) as fh:
        for line in fh:
            if line.startswith(">"):
                acc = line.split("|")[1]
                m = re.search(r" GN=(\S+)", line)
                genes[acc] = m.group(1) if m else ""
    return genes


def viral_entries(search_key):
    """Viral protein ID of the search's FASTA -> (strain, protein)."""
    out = {}
    if search_key == "deposited FASTA":
        for path in FASTA_VIRAL_DEPOSITED:
            strain = "swH1N1" if "swine" in path.name else "huH1N1"
            for line in open(path):
                if line.startswith(">"):
                    pid = line[1:].split()[0]
                    out[pid] = (strain, pid.rsplit("|", 1)[-1].rsplit("_", 1)[-1])
    else:
        for line in open(FASTA_VIRAL_CLEAN):
            if line.startswith(">"):
                w = line[1:].split()
                out[w[0]] = (w[1], w[2])
    return out


def accession_token(pid):
    """MT666905 for ssi|MT666905|A_swine_..._NP, the ID itself otherwise."""
    parts = pid.split("|")
    return parts[1] if len(parts) == 3 else pid


@cache
def go_set():
    """UniProt accessions of pig proteins annotated with GO:0051607 (annotation/)."""
    if not GO_TABLE.exists():
        return None
    return frozenset(pd.read_csv(GO_TABLE, sep="\t", dtype=str, comment="#")["accession"])


# ---------------------------------------------------------------- searches
EVIDENCE_COLUMNS = ["Experiment", "Sequence", "Modified sequence", "Charge", "Proteins",
                    "Protein group IDs", "Q-value", "PEP", "Intensity", "Decoy",
                    "Potential contaminant"]


class Search:
    """One MaxDIA search of the 44 sample runs: clean protein groups, LFQ, run-level detection."""

    def __init__(self, key):
        self.key = key
        self.project_dir = SEARCH_DIRS[key]
        self.txt_dir = self.project_dir / "results" / "combined" / "txt"
        pg = pd.read_csv(self.txt_dir / "proteinGroups.txt", sep="\t", low_memory=False)
        self.n_all = len(pg)
        flags = {c: pg[c].fillna("").eq("+") for c in ("Decoy", "Potential contaminant",
                                                          "Only identified by site")}
        self.n_contaminants = int(flags["Potential contaminant"].sum())
        self.n_decoys = int(flags["Decoy"].sum())
        clean = ~(flags["Decoy"] | flags["Potential contaminant"]
                  | flags["Only identified by site"])
        self.pg = pg[clean].set_index("id")
        lfq = self.pg[[f"LFQ intensity {s}" for s in samples()]].astype(float)
        lfq.columns = samples()
        self.lfq = lfq.where(lfq > 0)
        self.viral = viral_entries(key)
        self._evidence = None
        self._detected = None

    # -- protein IDs and labels
    @property
    def protein_ids(self):
        return self.pg["Protein IDs"].str.split(";")

    @property
    def majority_ids(self):
        return self.pg["Majority protein IDs"].str.split(";")

    @property
    @cache
    def labels(self):
        """Gene symbol of the first majority protein; viral: protein (strains)."""
        genes = pig_genes()
        out = {}
        for gid, ids in self.majority_ids.items():
            vir = [self.viral[i] for i in ids if i in self.viral]
            if vir:
                strains = "/".join(sorted({s[:2] for s, _ in vir}))
                out[gid] = f"{vir[0][1]} ({strains}H1N1)"
            else:
                out[gid] = next((genes[i] for i in ids if genes.get(i)), ids[0])
        return pd.Series(out)

    @property
    @cache
    def viral_groups(self):
        """Group ids that contain a viral protein -> sorted (strain, protein) members."""
        out = {}
        for gid, ids in self.protein_ids.items():
            vir = sorted({self.viral[i] for i in ids if i in self.viral})
            if vir:
                out[gid] = vir
        return out

    @property
    @cache
    def accession_index(self):
        """Accession -> group id (majority member first, then the group with most peptides)."""
        cand = {}
        majority = self.majority_ids
        peptides = self.pg["Peptides"].astype(int)
        for gid, ids in self.protein_ids.items():
            maj = set(majority[gid])
            for pid in ids:
                for key in {pid, accession_token(pid)}:
                    cand.setdefault(key, []).append((pid not in maj, -peptides[gid], gid))
        return {k: min(v)[2] for k, v in cand.items()}

    def group_of_accession(self, acc):
        return self.accession_index.get(acc)

    # -- evidence and run-level detection
    @property
    def evidence(self):
        if self._evidence is None:
            ev = pd.read_csv(self.txt_dir / "evidence.txt", sep="\t", low_memory=False,
                             usecols=lambda c: c in EVIDENCE_COLUMNS)
            for c in ("Decoy", "Potential contaminant"):
                ev[c] = ev[c].fillna("").eq("+") if c in ev.columns else False
            ev = ev[~ev["Decoy"]]
            ev["passes"] = (ev["Q-value"] <= Q_RUN) & (ev["PEP"] <= PEP_RUN)
            self._evidence = ev
        return self._evidence

    @property
    def detected(self):
        """Clean groups x samples: True if a precursor of the group passes in that run."""
        if self._detected is None:
            ok = self.evidence[self.evidence["passes"]][["Experiment", "Protein group IDs"]]
            ok = ok.assign(gid=ok["Protein group IDs"].astype(str).str.split(";")).explode("gid")
            ok["gid"] = ok["gid"].astype(int)
            det = ok.groupby(["gid", "Experiment"]).size().unstack(fill_value=0) > 0
            self._detected = det.reindex(index=self.pg.index, columns=samples(), fill_value=False)
        return self._detected

    def layer(self, name):
        if name == "as reported":
            return self.lfq
        if name == "run-filtered":
            return self.lfq.where(self.detected)
        raise ValueError(name)

    def peptide_count(self):
        """Clean peptide sequences of peptides.txt (no decoys, no contaminants)."""
        pep = pd.read_csv(self.txt_dir / "peptides.txt", sep="\t", low_memory=False,
                          usecols=lambda c: c in ("Sequence", "Decoy", "Potential contaminant"))
        bad = pep.get("Decoy", pd.Series("", index=pep.index)).fillna("").eq("+")
        bad |= pep["Potential contaminant"].fillna("").eq("+")
        return int((~bad).sum())


@cache
def search(key):
    print(f"loading search: {key} ({rel(SEARCH_DIRS[key] / 'results' / 'combined' / 'txt')})")
    return Search(key)


def available_searches():
    return [k for k, d in SEARCH_DIRS.items()
            if (d / "results" / "combined" / "txt" / "proteinGroups.txt").exists()]


def maxquant_version(txt_dir):
    p = pd.read_csv(txt_dir / "parameters.txt", sep="\t", dtype=str)
    return p.set_index(p.columns[0]).iloc[:, 0].get("Version", "?")


# ---------------------------------------------------------------- published table
@cache
def table_s3():
    """Supplementary Table S3: comparison, acc, gene, desc, lfc, padj, expr, plus the count in
    each sheet title ("List of 121 differentially expressed proteins")."""
    wb = openpyxl.load_workbook(TABLE_S3, read_only=True)
    recs, titles = [], {}
    for ws in wb.worksheets:
        rows = list(ws.iter_rows(values_only=True))
        titles[ws.title] = int(re.search(r"List of (\d+)", rows[0][0]).group(1))
        for r in rows[2:]:
            if r[0] is None or r[3] is None:
                continue
            recs.append((ws.title, str(r[0]).strip(), r[1], r[2], float(r[3]), float(r[4]), r[5]))
    df = pd.DataFrame(recs, columns=["comparison", "acc", "gene", "desc", "lfc", "padj", "expr"])
    if sorted(df["comparison"].unique()) != sorted(COMPARISON_NAMES):
        raise ValueError(f"unexpected sheets in Table S3: {sorted(df['comparison'].unique())}")
    df.attrs["titles"] = titles
    return df


# ---------------------------------------------------------------- paper pipeline
@cache
def _tukey_grid(k=3, df=19):
    """log10 of the studentized range survival function on a q grid (fast interpolation)."""
    q = np.concatenate([np.arange(0, 10, 0.01), np.arange(10, 60.01, 0.1)])
    p = stats.studentized_range.sf(q, k, df)
    return q, np.log10(np.maximum(p, 1e-300))


def tukey_p(q, k=3, df=19):
    grid_q, grid_lp = _tukey_grid(k, df)
    lp = np.interp(np.minimum(q, grid_q[-1]), grid_q, grid_lp)
    p = np.power(10.0, lp)
    return np.where(q >= grid_q[-1], 0.0, p)


def paper_stats(matrix, tissue):
    """Paper pipeline for one tissue. matrix: groups x samples, NaN = missing. Returns a
    DataFrame (index: analysed groups) with n_<group>, mean_<group>, anova_p and per comparison
    lfc, p and call (+1 up, -1 down, 0)."""
    cols = {g: samples_of(tissue, g) for g in GROUPS}
    sub = matrix[sum(cols.values(), [])]
    sub = sub[sub.notna().any(axis=1)]
    y = np.log2(sub).fillna(0.0)
    out = pd.DataFrame(index=sub.index)
    means, ss_within, n = {}, 0.0, {}
    for g, c in cols.items():
        v = y[c].values
        n[g] = v.shape[1]
        means[g] = v.mean(axis=1)
        ss_within = ss_within + ((v - means[g][:, None]) ** 2).sum(axis=1)
        out[f"n_{g}"] = sub[c].notna().sum(axis=1)
        out[f"mean_{g}"] = means[g]
    k, n_tot = len(cols), sum(n.values())
    dfw = n_tot - k
    mse = ss_within / dfw
    grand = y.values.mean(axis=1)
    ss_between = sum(n[g] * (means[g] - grand) ** 2 for g in cols)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = (ss_between / (k - 1)) / mse
        out["anova_p"] = np.where(mse > 0, stats.f.sf(f, k - 1, dfw), 0.0)
        for t, a, b in COMPARISONS:
            if t != tissue:
                continue
            name = f"{t} {a} vs {b}"
            lfc = means[a] - means[b]
            q = np.abs(lfc) / np.sqrt(mse / 2 * (1 / n[a] + 1 / n[b]))
            p = np.where(mse > 0, tukey_p(np.nan_to_num(q, nan=0.0, posinf=1e9), k, dfw),
                         np.where(lfc == 0, 1.0, 0.0))
            out[f"{name}|lfc"] = lfc
            out[f"{name}|p"] = p
            sig = (p < P_CUTOFF) & (np.abs(lfc) >= FC_CUTOFF)
            out[f"{name}|call"] = np.where(sig, np.sign(lfc), 0).astype(int)
    return out


@cache
def pipeline(search_key, layer):
    """Paper pipeline on both tissues of a search layer: dict tissue -> paper_stats frame."""
    m = search(search_key).layer(layer)
    return {t: paper_stats(m, t) for t in TISSUES}


def comparison_frame(search_key, layer, name):
    """Per analysed group: lfc, p, call of one comparison."""
    t = name.split()[0]
    st = pipeline(search_key, layer)[t]
    return st[[f"{name}|lfc", f"{name}|p", f"{name}|call"]].set_axis(["lfc", "p", "call"], axis=1)


def outcome_of(pub_lfc, ours_row):
    if ours_row is None:
        return "not quantified"
    if np.sign(ours_row["lfc"]) != np.sign(pub_lfc) or ours_row["lfc"] == 0:
        return "opposite direction"
    return "DEP, same direction" if ours_row["call"] != 0 else "same direction, not DEP"


@cache
def agreement(search_key, layer):
    """Every published DEP with its outcome in the search layer (paper pipeline)."""
    s = search(search_key)
    pub = table_s3()
    rows = []
    frames = {n: comparison_frame(search_key, layer, n) for n in COMPARISON_NAMES}
    for r in pub.itertuples():
        gid = s.group_of_accession(r.acc)
        fr = frames[r.comparison]
        ours = fr.loc[gid] if gid is not None and gid in fr.index else None
        rows.append({
            "comparison": r.comparison, "acc": r.acc, "gene": r.gene, "published_lfc": r.lfc,
            "published_padj": r.padj, "presence": abs(r.lfc) >= PRESENCE_LFC,
            "group_id": gid, "identified": gid is not None,
            "our_lfc": np.nan if ours is None else ours["lfc"],
            "our_p": np.nan if ours is None else ours["p"],
            "our_call": 0 if ours is None else int(ours["call"]),
            "outcome": outcome_of(r.lfc, ours),
        })
    return pd.DataFrame(rows)


def agreement_summary(search_key, layer, comparisons=None):
    """Direction agreement and recall of published DEPs (pooled over the comparisons given)."""
    a = agreement(search_key, layer)
    if comparisons is not None:
        a = a[a.comparison.isin(comparisons)]
    quant = a[a.outcome != "not quantified"]
    same = quant.outcome.isin(["DEP, same direction", "same direction, not DEP"])
    return {
        "published": len(a),
        "quantified": len(quant),
        "direction": same.mean() if len(quant) else np.nan,
        "recall": (a.outcome == "DEP, same direction").mean() if len(a) else np.nan,
        "recall_quantitative": (a[~a.presence].outcome == "DEP, same direction").mean(),
        "recall_presence": (a[a.presence].outcome == "DEP, same direction").mean(),
    }


def our_dep_counts(search_key, layer):
    out = {}
    for name in COMPARISON_NAMES:
        fr = comparison_frame(search_key, layer, name)
        out[name] = {"up": int((fr.call > 0).sum()), "down": int((fr.call < 0).sum())}
    return out


# ---------------------------------------------------------------- two searches compared
def match_groups(key_a, key_b):
    """Group id of search a -> group id of search b holding a's leading majority protein."""
    a, b = search(key_a), search(key_b)
    lead = a.majority_ids.str[0]
    return lead.map(lambda pid: b.group_of_accession(pid)
                    if b.group_of_accession(pid) is not None
                    else b.group_of_accession(accession_token(pid)))


def lfq_correlation(key_a, key_b):
    """Median over samples of the Pearson r of log2 LFQ (as reported) of matched groups, and
    the number of groups compared."""
    m = match_groups(key_a, key_b).dropna().astype(int)
    m = m[~m.duplicated(keep=False)]
    xa = np.log2(search(key_a).lfq.loc[m.index])
    xb = np.log2(search(key_b).lfq.loc[m.values]).set_axis(m.index)
    rs = [xa[smp].corr(xb[smp]) for smp in samples()]
    return float(np.median(rs)), len(m)


def dep_overlap(key_a, key_b, layer, comparison):
    """Share of a's DEPs that are DEPs of the same sign in b, and the DEP counts."""
    fa = comparison_frame(key_a, layer, comparison)
    fb = comparison_frame(key_b, layer, comparison)
    m = match_groups(key_a, key_b)
    dep_a = fa[fa.call != 0]
    same = 0
    for gid, r in dep_a.iterrows():
        gb = m.get(gid)
        if gb is not None and not pd.isna(gb) and int(gb) in fb.index:
            same += int(fb.loc[int(gb), "call"] == r.call)
    return {"a": len(dep_a), "b": int((fb.call != 0).sum()),
            "a_in_b": same / len(dep_a) if len(dep_a) else np.nan}


# ---------------------------------------------------------------- viral detection per sample
def viral_precursors(search_key):
    """Passing viral precursors per sample by strain specificity: DataFrame samples x
    (sw_only, hu_only, shared, total, log2_intensity)."""
    s = search(search_key)
    ev = s.evidence
    prots = ev["Proteins"].fillna("").str.split(";")
    vir_mask = prots.map(lambda a: any(x in s.viral for x in a))
    v = ev[vir_mask & ev["passes"]].copy()
    v["strain"] = prots[v.index].map(
        lambda a: "+".join(sorted({s.viral[x][0] for x in a if x in s.viral}))
        + ("+host" if any(x not in s.viral for x in a) else ""))
    v["prec"] = v["Modified sequence"] + "/" + v["Charge"].astype(str)
    rows = {}
    for smp in samples():
        x = v[v.Experiment == smp]
        rows[smp] = {
            "sw_only": x[x.strain == "swH1N1"].prec.nunique(),
            "hu_only": x[x.strain == "huH1N1"].prec.nunique(),
            "shared": x[x.strain == "huH1N1+swH1N1"].prec.nunique(),
            "with_host": x[x.strain.str.contains("host")].prec.nunique(),
            "total": x[~x.strain.str.contains("host")].prec.nunique(),
            "log2_intensity": float(np.log2(x[~x.strain.str.contains("host")].Intensity.sum() + 1)),
        }
    return pd.DataFrame(rows).T


ISG_GENES = ["MX1", "MX2", "ISG15", "OASL", "IFIT3", "IFI6"]  # Fig. 2, global proteomics


def isg_groups(search_key):
    """Gene -> group id for the ISG proteins of Fig. 2 (via their Table S3 accessions)."""
    pub = table_s3()
    s = search(search_key)
    out = {}
    for g in ISG_GENES:
        accs = pub.loc[pub.gene == g, "acc"].unique()
        gid = next((s.group_of_accession(a) for a in accs if s.group_of_accession(a)), None)
        if gid is not None:
            out[g] = gid
    return out


def isg_score(search_key, layer="as reported"):
    """Mean z-score per tissue of the log2 ISG proteins (missing -> tissue minimum)."""
    s = search(search_key)
    m = np.log2(s.layer(layer).loc[list(isg_groups(search_key).values())])
    score = {}
    for t in TISSUES:
        cols = samples_of(t)
        x = m[cols]
        x = x.apply(lambda r: r.fillna(r.min()), axis=1)
        z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
        score.update(z.mean().to_dict())
    return pd.Series(score)[samples()]


# ---------------------------------------------------------------- self test
def selftest():
    """Quick checks of the reconstructed pipeline before the questions run."""
    # Tukey interpolation against the exact survival function
    for q in (0.5, 2.0, 3.6, 5.0, 12.0):
        exact = stats.studentized_range.sf(q, 3, 19)
        assert abs(tukey_p(np.array([q]))[0] - exact) <= 1e-3 * exact + 1e-12, q
    # Tukey on one protein against scipy.stats.tukey_hsd
    rng = np.random.default_rng(1)
    g = [rng.normal(20, 1, 6), rng.normal(21, 1, 8), rng.normal(20.5, 1, 8)]
    ref = stats.tukey_hsd(*g).pvalue
    cols = samples_of("UT")
    m = pd.DataFrame([np.power(2.0, np.concatenate(g))], columns=cols)
    st = paper_stats(m, "UT")
    assert abs(st["UT sw vs ctr|p"].iloc[0] - ref[1, 0]) < 1e-3
    assert abs(st["UT sw vs hu|p"].iloc[0] - ref[1, 2]) < 1e-3
    # design: 44 samples, 6/8/8 per tissue, groups as in GEO
    d = design()
    assert len(d) == 44 and all(len(samples_of(t, "ctr")) == 6 for t in TISSUES)
    assert all(geo_groups()[p] == g for p, g in zip(d.pig, d.group))
    print("self test passed")
