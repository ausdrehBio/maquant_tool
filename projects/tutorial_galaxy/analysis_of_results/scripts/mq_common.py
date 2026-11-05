"""Shared helpers for the tutorial_galaxy analysis scripts q01 to q12.

Reads MaxQuant 2.9.0 output from results/combined/txt/ (read only) and
writes figures to analysis_of_results/figures/. All paths derive from this
file's location and are printed relative to the project folder.

Definitions used by all scripts:
- decoy: MaxQuant's "Decoy" column is "+". A group that holds a REV__
  member next to a target protein is not a decoy.
- contaminant: MaxQuant's "Potential contaminant" column is "+".
- clean set: all protein groups except decoys and contaminant groups that
  contain no protein from the human FASTA (e.g. bovine CON__ entries).
  Human serum proteins on the contaminant list, such as albumin, stay,
  because in serum they are genuine sample proteins.
- abundance: "Intensity <sample>"; share = intensity / sample total * 100.
"""
from dataclasses import dataclass
from functools import cache
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SCRIPTS_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPTS_DIR.parent
PROJECT_DIR = ANALYSIS_DIR.parent
TXT_DIR = PROJECT_DIR / "results" / "combined" / "txt"
FIG_DIR = ANALYSIS_DIR / "figures"

# High-abundance serum proteins removed by common depletion columns (the
# 14 targets of a top-14 depletion), as gene names. IGLC* matches by prefix.
TOP_SERUM_GENES = frozenset({
    "ALB",                                        # albumin
    "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGKC",   # IgG heavy chains, kappa light chain
    "SERPINA1",                                   # alpha-1-antitrypsin
    "IGHA1", "IGHA2",                             # IgA
    "TF",                                         # serotransferrin
    "HP",                                         # haptoglobin
    "FGA", "FGB", "FGG",                          # fibrinogen
    "A2M",                                        # alpha-2-macroglobulin
    "ORM1", "ORM2",                               # alpha-1-acid glycoprotein
    "IGHM",                                       # IgM
    "APOA1", "APOA2",                             # apolipoproteins A-I, A-II
    "C3",                                         # complement C3
    "TTR",                                        # transthyretin
})

# ---------------------------------------------------------------- plot style
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
NEUTRAL = INK2  # marks that are not about one sample
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")  # categorical slots 1-4, fixed order

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "text.color": INK,
    "axes.labelcolor": INK2,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.titlecolor": INK,
    "xtick.color": AXIS,
    "ytick.color": AXIS,
    "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2,
    "legend.frameon": False,
    "legend.fontsize": 9,
})


def sample_colors(names):
    """Fixed colour per sample, in column order, identical in every figure."""
    return {name: SERIES[i] for i, name in enumerate(names)}


def value_grid(ax, axis="y"):
    """Solid hairline grid on the value axis only, behind the marks."""
    ax.grid(axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def bar_kw():
    """Surface-coloured edge = visible gap between adjacent bars."""
    return {"edgecolor": SURFACE, "linewidth": 1.5}


def label_bars(ax, bars, texts):
    ax.bar_label(bars, labels=texts, padding=3, color=INK2, fontsize=9)


def save_fig(fig, name):
    """Write figures/<name>.png; return its path relative to analysis_of_results/."""
    FIG_DIR.mkdir(exist_ok=True)
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"figure written: {rel(path)}")
    return path.relative_to(ANALYSIS_DIR).as_posix()


# ---------------------------------------------------------------- paths, io
def rel(path):
    """Path relative to the project folder, for printing. Never a complete host path."""
    return Path(path).resolve().relative_to(PROJECT_DIR).as_posix()


def col(df, name):
    """Column by header name, case-insensitive (MaxQuant changes case between versions)."""
    hits = [c for c in df.columns if c.lower() == name.lower()]
    if len(hits) != 1:
        raise KeyError(f"expected one column named {name!r} (any case), found {hits}")
    return hits[0]


@cache
def _read(file_name):
    return pd.read_csv(TXT_DIR / file_name, sep="\t", low_memory=False)


def load_protein_groups():
    return _read("proteinGroups.txt").copy()


def load_peptides():
    return _read("peptides.txt").copy()


def load_msms():
    return _read("msms.txt").copy()


def load_summary():
    """summary.txt, one row per raw file (the 'Total' row removed)."""
    s = _read("summary.txt").copy()
    raw = s[col(s, "Raw file")].fillna("")
    return s[raw.ne("") & raw.ne("Total")]


def mqpar_value(tag):
    """Value of a tag in the project's parameter/mqpar.xml (read only)."""
    text = (PROJECT_DIR / "parameter" / "mqpar.xml").read_text()
    m = re.search(rf"<{tag}>([^<]*)</{tag}>", text)
    return m.group(1) if m else None


# ---------------------------------------------------------------- flags, sets
def plus(df, name):
    """True where a MaxQuant flag column holds '+'."""
    return df[col(df, name)].fillna("").eq("+")


def is_decoy(df):
    return plus(df, "Decoy")


def is_contaminant(df):
    return plus(df, "Potential contaminant")


def has_db_protein(pg):
    """True if the group contains at least one protein from the FASTA (not CON__, not REV__)."""
    ids = pg[col(pg, "Protein IDs")].fillna("").str.split(";")
    return ids.map(lambda xs: any(x and not x.startswith(("CON__", "REV__")) for x in xs)).astype(bool)


def clean(pg):
    """Protein groups without decoys and without contaminant groups lacking a human protein."""
    keep = ~is_decoy(pg) & ~(is_contaminant(pg) & ~has_db_protein(pg))
    return pg[keep].copy()


def samples(df):
    """Sample (experiment) names from the 'Intensity <name>' columns, in column order."""
    out = [m.group(1) for c in df.columns if (m := re.fullmatch(r"intensity (.+)", c, re.IGNORECASE))]
    if not out:
        raise KeyError("no 'Intensity <sample>' columns found")
    return out


def intensity(df, sample):
    return df[col(df, f"Intensity {sample}")].fillna(0).astype(float)


def lfq(df, sample):
    return df[col(df, f"LFQ intensity {sample}")].fillna(0).astype(float)


def share(df, sample):
    """Intensity as % of the sample's summed intensity within df."""
    x = intensity(df, sample)
    return x / x.sum() * 100


def genes(pg):
    return pg[col(pg, "Gene names")].fillna("").str.split(";")


def is_albumin(pg):
    return genes(pg).map(lambda g: "ALB" in g).astype(bool)


def is_top_serum(pg):
    return genes(pg).map(lambda g: any(x in TOP_SERUM_GENES or x.startswith("IGLC") for x in g)).astype(bool)


def identified_in(pg, sample):
    """At least one identified peptide in this sample."""
    return pg[col(pg, f"Peptides {sample}")].fillna(0).gt(0)


def labels(pg, maxlen=32):
    """Short display name: first protein name, else first gene, else first protein ID."""
    names = pg[col(pg, "Protein names")].fillna("").str.split(";").str[0]
    gene = genes(pg).str[0]
    ids = pg[col(pg, "Majority protein IDs")].fillna("").str.split(";").str[0]
    out = names.where(names.ne(""), gene.where(gene.ne(""), ids))
    return out.map(lambda s: s if len(s) <= maxlen else s[: maxlen - 1] + "…")


def depleted_sample(pg):
    """(depleted, undepleted): the sample with the lower albumin share is the depleted one."""
    cl = clean(pg)
    alb = is_albumin(cl)
    if alb.sum() != 1:
        raise ValueError(f"expected one albumin group, found {int(alb.sum())}")
    names = samples(cl)
    if len(names) != 2:
        raise ValueError(f"expected two samples, found {names}")
    alb_share = {s: share(cl, s)[alb].iloc[0] for s in names}
    dep = min(alb_share, key=alb_share.get)
    und = next(s for s in names if s != dep)
    return dep, und


# ---------------------------------------------------------------- formatting
def n(x):
    return f"{int(x):,}"


def pct(x, digits=1):
    return f"{x:.{digits}f}%"


def sci(x):
    return f"{x:.2e}"


def md_table(rows, header):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


@dataclass
class Answer:
    """One question's result. str() gives the Markdown section body; `short` is a
    one-line summary for the comparison with the Galaxy tutorial."""
    question: str
    text: str
    figure: str  # relative to analysis_of_results/
    short: str

    def __str__(self):
        return f"{self.text}\n\n![{self.question}]({self.figure})\n"


def answer(question, text, figure, short):
    return Answer(question, text, figure, short)
