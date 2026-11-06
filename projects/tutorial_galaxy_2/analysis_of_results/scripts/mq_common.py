"""Shared helpers for the tutorial_galaxy_2 analysis scripts q01 to q10.

Inputs (read only):
- MaxQuant output: results/combined/txt/{proteinGroups,evidence}.txt
- MSstats output of run_msstats.R: analysis_of_results/msstats/
- FASTA: data/input_protein_database.fasta (gene and protein names)
Outputs: analysis_of_results/figures/*.png and analysis_of_results/answers.md

For a test on other MaxQuant output, inputs and outputs can be redirected with the
environment variables TG2_TXT_DIR, TG2_MSSTATS_DIR and TG2_OUT_DIR.
Paths are printed relative to the project folder, never as complete host paths.

Definitions used by all scripts:
- condition of a run: from the MSstats annotation (metast or RDEB)
- human protein: "HUMAN" in the protein ID (the tutorial's Select step)
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
TXT_DIR = Path(os.environ.get("TG2_TXT_DIR", PROJECT_DIR / "results" / "combined" / "txt"))
MSSTATS_DIR = Path(os.environ.get("TG2_MSSTATS_DIR", ANALYSIS_DIR / "msstats"))
OUT_DIR = Path(os.environ.get("TG2_OUT_DIR", ANALYSIS_DIR))
FIG_DIR = OUT_DIR / "figures"
FASTA = PROJECT_DIR / "data" / "input_protein_database.fasta"
COMPARISON_MATRIX = PROJECT_DIR / "parameter" / "input_comparison_matrix.tabular"

FC_CUTOFF = 0.58  # log2(1.5), the tutorial's fold-change threshold
P_CUTOFF = 0.05   # adjusted p-value threshold

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


def load_evidence():
    return _read(TXT_DIR / "evidence.txt").copy()


def load_msstats(name):
    """MSstats table written by run_msstats.R, e.g. 'ProteinLevelData'."""
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


def plus(df, name):
    return df[col(df, name)].fillna("").eq("+")


def is_decoy(df):
    """MaxQuant 2.9.0 has a Decoy column, older versions a Reverse column."""
    name = "Decoy" if any(c.lower() == "decoy" for c in df.columns) else "Reverse"
    return plus(df, name)


def is_contaminant(df):
    return plus(df, "Potential contaminant")


def has_human(series):
    return series.fillna("").astype(str).str.contains("HUMAN", regex=False)


# ---------------------------------------------------------------- samples, conditions
def core_name(x):
    """Run name without Galaxy suffixes and PRIDE prefix: 'metast_cSCC1'."""
    x = re.sub(r"(\.raw)?(\.thermo)?(\.raw)?$", "", str(x))
    return re.sub(r"^Experiment[0-9]+_", "", x)


@cache
def _annotation():
    a = _read(MSSTATS_DIR / "annotation_used.tsv").copy()
    a["sample"] = a[col(a, "Raw.file")].map(core_name)
    return a


def samples():
    """Sample names in annotation order (metast first, then RDEB)."""
    return list(_annotation()["sample"])


def conditions():
    """Condition names in annotation order."""
    return list(dict.fromkeys(_annotation()[col(_annotation(), "Condition")]))


def sample_of_run(run):
    """Sample name for an MSstats run label. MSstats strips dots from file names
    ('x.raw.thermo' -> 'xrawthermo'), so map through the annotation, not by string rules."""
    a = _annotation()
    raw = a[col(a, "Raw.file")].astype(str)
    lookup = dict(zip(raw, a["sample"]))
    lookup.update(zip(raw.str.replace(".", "", regex=False), a["sample"]))
    lookup.update(zip(a[col(a, "Run")].astype(str), a["sample"]))
    return lookup[str(run)]


def condition_of(name):
    a = _annotation()
    return dict(zip(a["sample"], a[col(a, "Condition")]))[core_name(name)]


def condition_colors():
    return {c: SERIES[i] for i, c in enumerate(conditions())}


def sample_intensity_columns(pg):
    """{sample: 'Intensity <experiment>'} for the protein groups table."""
    out = {}
    for c in pg.columns:
        m = re.fullmatch(r"intensity (.+)", c, re.IGNORECASE)
        if m:
            out[core_name(m.group(1))] = c
    return {s: out[s] for s in samples()}


def contrast():
    """(name, numerator, denominator) of the tutorial's comparison matrix."""
    cm = pd.read_csv(COMPARISON_MATRIX, sep="\t")
    row = cm.iloc[0]
    up = [c for c in cm.columns[1:] if row[c] > 0]
    down = [c for c in cm.columns[1:] if row[c] < 0]
    return row.iloc[0], up[0], down[0]


# ---------------------------------------------------------------- protein names
@cache
def _fasta():
    table = {}
    for line in FASTA.read_text().splitlines():
        if not line.startswith(">"):
            continue
        parts = line[1:].split("|")
        acc = parts[1] if len(parts) > 2 else parts[0].split()[0]
        gene = re.search(r"\bGN=(\S+)", line)
        name = re.search(r"_HUMAN (.*?) OS=", line)
        table[acc] = (gene.group(1) if gene else "", name.group(1) if name else "")
    return table


def accession(protein_id):
    """UniProt accession from 'sp|P02768|ALBU_HUMAN' or 'CON__P02768-1'."""
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
    """str() gives the Markdown section body; `short` feeds the tutorial comparison."""
    question: str
    text: str
    figure: str
    short: str

    def __str__(self):
        return f"{self.text}\n\n![{self.question}]({self.figure})\n"


def answer(question, text, figure, short):
    return Answer(question, text, figure, short)
