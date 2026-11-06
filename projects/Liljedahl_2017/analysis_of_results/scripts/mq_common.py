"""Shared helpers for the Liljedahl_2017 analysis scripts q01 to q07.

Inputs (read only):
- MaxQuant output: results/combined/txt/{proteinGroups,peptides,evidence,summary,parameters}.txt
- MSstats output of run_msstats.R: analysis_of_results/msstats/
- sample annotation and contrasts: parameter/sample_annotation.tsv, parameter/comparison_matrix.tsv
- FASTA: the one data/*.fasta (gene and protein names)
Outputs: analysis_of_results/figures/*.png and analysis_of_results/answers.md

For a test on other MaxQuant output, inputs and outputs can be redirected with the
environment variables LJ_TXT_DIR, LJ_MSSTATS_DIR and LJ_OUT_DIR.
Paths are printed relative to the project folder, never as complete host paths.

Definitions used by all scripts:
- groups: CTRL = db/+ vehicle (healthy control), VEH = db/db vehicle, LIRA = db/db liraglutide,
  INS = db/db insulin; 5 mice per group, 2 technical replicates (runs) per mouse
- target protein: a member of a protein group that is neither CON__ nor REV__ (the FASTA is
  mouse only); contaminant groups with a target member, such as albumin or keratins, are kept,
  because in kidney tissue they can be genuine mouse proteins
- clean set: protein groups with a target member, without decoys and without groups only
  identified by site
- significant: MSstats adjusted p-value < 0.05; up/down: log2FC beyond +-0.58 (1.5-fold)
"""
from dataclasses import dataclass
from functools import cache
from pathlib import Path
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SCRIPTS_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPTS_DIR.parent
PROJECT_DIR = ANALYSIS_DIR.parent
TXT_DIR = Path(os.environ.get("LJ_TXT_DIR", PROJECT_DIR / "results" / "combined" / "txt"))
MSSTATS_DIR = Path(os.environ.get("LJ_MSSTATS_DIR", ANALYSIS_DIR / "msstats"))
OUT_DIR = Path(os.environ.get("LJ_OUT_DIR", ANALYSIS_DIR))
FIG_DIR = OUT_DIR / "figures"
ANNOTATION = PROJECT_DIR / "parameter" / "sample_annotation.tsv"
COMPARISON_MATRIX = PROJECT_DIR / "parameter" / "comparison_matrix.tsv"

FC_CUTOFF = 0.58  # log2(1.5)
P_CUTOFF = 0.05   # adjusted p-value threshold

CONDITIONS = ["CTRL", "VEH", "LIRA", "INS"]
CONDITION_LABELS = {
    "CTRL": "db/+ vehicle (healthy control)",
    "VEH": "db/db vehicle",
    "LIRA": "db/db liraglutide",
    "INS": "db/db insulin",
}

# ---------------------------------------------------------------- plot style
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
NEUTRAL = INK2
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")  # categorical slots 1-4, fixed order
DIVERGING = ["#184f95", "#86b6ef", "#f0efec", "#f19c9b", "#c23b3a"]  # blue - gray - red

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlecolor": INK,
    "xtick.color": AXIS, "ytick.color": AXIS, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "legend.frameon": False, "legend.fontsize": 9,
})


def value_grid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def rotate_xticks(ax, size=8, angle=60):
    """Rotated tick labels that end under their tick (anchored at the right end)."""
    ax.tick_params(axis="x", length=0, labelsize=size)
    for t in ax.get_xticklabels():
        t.set_rotation(angle)
        t.set_horizontalalignment("right")
        t.set_rotation_mode("anchor")


def bar_kw():
    return {"edgecolor": SURFACE, "linewidth": 1.0}


def label_bars(ax, bars, texts, **kw):
    ax.bar_label(bars, labels=texts, padding=3, color=INK2, fontsize=kw.pop("fontsize", 9), **kw)


def condition_legend(ax, **kw):
    """Legend with one patch per condition."""
    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=f"{k}: {CONDITION_LABELS[k]}") for k, c in condition_colors().items()]
    ax.legend(handles=handles, **kw)


def run_bars(ax, values, ylabel):
    """One bar per run (experiment order), coloured by condition."""
    exps = list(values.index)
    colors = [condition_colors()[condition_of(e)] for e in exps]
    bars = ax.bar(range(len(exps)), values.values, color=colors, **bar_kw())
    ax.set_xticks(range(len(exps)), exps)
    rotate_xticks(ax, size=7, angle=90)
    ax.set_ylabel(ylabel)
    value_grid(ax)
    return bars


def save_fig(fig, name):
    """Write figures/<name>.png; return its path relative to the answers file."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"figure written: {rel(path)}")
    return path.relative_to(OUT_DIR).as_posix()


# ---------------------------------------------------------------- paths, io
def rel(path):
    """Path relative to the project folder for printing; outside it, only the file name."""
    p = Path(path).resolve()
    try:
        return p.relative_to(PROJECT_DIR).as_posix()
    except ValueError:
        return p.name


def col(df, name):
    """Column by header name, case-insensitive."""
    hits = [c for c in df.columns if c.lower() == name.lower()]
    if len(hits) != 1:
        raise KeyError(f"expected one column named {name!r} (any case), found {hits}")
    return hits[0]


@cache
def _read(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def load_protein_groups():
    return _read(TXT_DIR / "proteinGroups.txt").copy()


def load_peptides():
    return _read(TXT_DIR / "peptides.txt").copy()


def load_evidence():
    return _read(TXT_DIR / "evidence.txt").copy()


def load_summary():
    """summary.txt, one row per raw file (the Total row removed), with an Experiment column
    taken from the sample annotation."""
    s = _read(TXT_DIR / "summary.txt").copy()
    raw = s[col(s, "Raw file")].astype(str)
    s = s[raw.str.lower() != "total"].copy()
    lookup = dict(zip(annotation()["Raw.file"], annotation()["Experiment"]))
    s["experiment"] = s[col(s, "Raw file")].map(lookup)
    if s["experiment"].isna().any():
        raise ValueError("summary.txt has raw files that are not in the sample annotation")
    return s.set_index("experiment").loc[experiments()]


def load_msstats(name):
    """MSstats table written by run_msstats.R, e.g. 'ComparisonResult'."""
    return _read(MSSTATS_DIR / f"{name}.tsv").copy()


def msstats_info():
    """Key-value lines of msstats/run_info.txt."""
    info = {}
    for line in (MSSTATS_DIR / "run_info.txt").read_text().splitlines():
        key, _, value = line.partition("\t")
        info[key] = value if value else key
    return info


def maxquant_version():
    """MaxQuant version from the Version row of parameters.txt, e.g. 2.9.0.0."""
    for line in (TXT_DIR / "parameters.txt").read_text().splitlines():
        key, _, value = line.partition("\t")
        if key.strip().lower() == "version":
            return value.strip()
    return "version unknown"


# ---------------------------------------------------------------- row classes
def plus(df, name):
    return df[col(df, name)].fillna("").eq("+")


def is_decoy(df):
    """MaxQuant 2.9.0 has a Decoy column, older versions a Reverse column."""
    name = "Decoy" if any(c.lower() == "decoy" for c in df.columns) else "Reverse"
    return plus(df, name)


def is_contaminant(df):
    return plus(df, "Potential contaminant")


def has_target(series):
    """True where a ';'-separated ID list has a member that is neither CON__ nor REV__."""
    def any_target(ids):
        return any(i and not i.startswith(("CON__", "REV__")) for i in str(ids).split(";"))
    return series.fillna("").map(any_target)


def clean_mask(pg):
    """Clean set: target member, no decoy, not only identified by site."""
    mask = has_target(pg[col(pg, "Protein IDs")]) & ~is_decoy(pg)
    if any(c.lower() == "only identified by site" for c in pg.columns):
        mask &= ~plus(pg, "Only identified by site")
    return mask


# ---------------------------------------------------------------- samples, conditions
@cache
def _annotation():
    a = pd.read_csv(ANNOTATION, sep="\t", dtype=str)
    a["order"] = a["Condition"].map(CONDITIONS.index)
    return a.sort_values(["order", "Experiment"]).drop(columns="order").reset_index(drop=True)


def annotation():
    return _annotation().copy()


def experiments():
    """Experiment names (= MaxQuant sample names), CTRL first, then VEH, LIRA, INS."""
    return list(_annotation()["Experiment"])


def condition_of(experiment):
    return dict(zip(_annotation()["Experiment"], _annotation()["Condition"]))[experiment]


def mouse_of(experiment):
    return dict(zip(_annotation()["Experiment"], _annotation()["BioReplicate"]))[experiment]


def partner_of(experiment):
    """The other technical replicate of the same mouse, according to the SDRF."""
    a = _annotation()
    mouse = mouse_of(experiment)
    other = a[(a["BioReplicate"] == mouse) & (a["Experiment"] != experiment)]["Experiment"]
    return other.iloc[0]


def mice():
    """Mouse IDs (BioReplicate) in condition order."""
    return list(dict.fromkeys(_annotation()["BioReplicate"]))


def condition_of_mouse(mouse):
    return dict(zip(_annotation()["BioReplicate"], _annotation()["Condition"]))[mouse]


def condition_colors():
    return {c: SERIES[i] for i, c in enumerate(CONDITIONS)}


def experiment_columns(pg, prefix):
    """{experiment: '<prefix> <experiment>'}, e.g. prefix 'Intensity' or 'LFQ intensity'."""
    out = {}
    for c in pg.columns:
        m = re.fullmatch(rf"{re.escape(prefix)} (.+)", c, re.IGNORECASE)
        if m:
            out[m.group(1)] = c
    missing = [e for e in experiments() if e not in out]
    if missing:
        raise KeyError(f"no '{prefix} <experiment>' column for {missing}")
    return {e: out[e] for e in experiments()}


def log2_matrix(pg, prefix="LFQ intensity"):
    """log2 values (0 -> NaN) of the clean set, rows = protein groups, columns = experiments."""
    clean = pg[clean_mask(pg)]
    cols = experiment_columns(clean, prefix)
    m = clean[list(cols.values())].astype(float).replace(0, np.nan)
    m.columns = list(cols.keys())
    m.index = clean[col(clean, "Protein IDs")]
    return np.log2(m)


def mouse_matrix(m):
    """Mean of the technical replicates per mouse (NaN only if both are missing)."""
    by_mouse = m.T.groupby(lambda e: mouse_of(e)).mean().T
    return by_mouse[mice()]


def contrasts():
    """[(label, numerator, denominator)] from the comparison matrix, in file order."""
    cm = pd.read_csv(COMPARISON_MATRIX, sep="\t")
    out = []
    for _, row in cm.iterrows():
        up = [c for c in cm.columns[1:] if row[c] > 0]
        down = [c for c in cm.columns[1:] if row[c] < 0]
        out.append((row.iloc[0], up[0], down[0]))
    return out


# ---------------------------------------------------------------- protein names
@cache
def _fasta():
    paths = sorted((PROJECT_DIR / "data").glob("*.fasta"))
    if len(paths) != 1:
        raise FileNotFoundError(f"expected one FASTA in data/, found {len(paths)}")
    table = {}
    for line in paths[0].read_text().splitlines():
        if not line.startswith(">"):
            continue
        parts = line[1:].split("|")
        acc = parts[1] if len(parts) > 2 else parts[0].split()[0]
        gene = re.search(r"\bGN=(\S+)", line)
        name = re.search(r"^>\S+ (.*?) OS=", line)
        table[acc] = (gene.group(1) if gene else "", name.group(1) if name else "")
    return table


def fasta_genes():
    """{accession: gene} of the FASTA."""
    return {acc: g for acc, (g, _) in _fasta().items()}


def accession(protein_id):
    """UniProt accession of the first member: 'P07724', 'sp|P07724|ALBU_MOUSE' or 'CON__P02769'."""
    pid = str(protein_id).split(";")[0]
    if "|" in pid:
        return pid.split("|")[1]
    return pid.removeprefix("CON__").removeprefix("REV__").split("-")[0]


def gene_of(protein_id):
    gene, _ = _fasta().get(accession(protein_id), ("", ""))
    return gene or accession(protein_id)


def name_of(protein_id):
    _, name = _fasta().get(accession(protein_id), ("", ""))
    return name


# ---------------------------------------------------------------- formatting, answers
def n(x):
    return f"{int(x):,}"


def pct(x, digits=1):
    return f"{x:.{digits}f}%"


def md_table(rows, header):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


@dataclass
class Answer:
    """str() gives the Markdown section body; `short` feeds the comparison with the paper."""
    question: str
    text: str
    figure: str
    short: str

    def __str__(self):
        return f"{self.text}\n\n![{self.question}]({self.figure})\n"


def answer(question, text, figure, short):
    return Answer(question, text, figure, short)
