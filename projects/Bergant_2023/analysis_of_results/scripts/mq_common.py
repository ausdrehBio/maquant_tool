"""Shared code for the Bergant_2023 comparison suite (q01 to q08).

Question of the suite: does our MaxQuant 2.9.0 workflow reproduce the published results of
Bergant et al. (2023, Nat Commun 14:4906, doi 10.1038/s41467-023-40469-6) for the IAV_FPMS
dataset (PRIDE PXD040456: A549 cells, mock-infected or infected with six PR8 NS1 mutants,
6 replicates each, 42 runs)?

Data layers (read only):
- ours:      results/combined/txt/ (MaxQuant 2.9.0, the paper's 42 runs)
- deposited: data/deposited/proteinGroups_IAV_FPMS.txt (the authors' MaxQuant 1.5.6.2 output;
             67 runs, 25 of which the paper does not use)
- published: data/deposited/Supplementary_Data_5_IAV_FPMS.xlsx (the paper's filtered, imputed
             log2 LFQ matrix of the 42 samples with fold changes, p and adjusted p values; it
             comes from a search that is not deposited)
- ours-67:   ../Bergant_2023_67runs/results/combined/txt/ (optional, Part B: MaxQuant 2.9.0 on
             all 67 runs); questions add it when it exists
Further inputs: parameter/design.txt, parameter/mqpar.xml, parameter/deposited/, the FASTA files
in data/, annotation/ (gene symbols and the GO:0051607 set, written by make_annotation.py).
Outputs: analysis_of_results/answers.md and verdicts.tsv (in git), figures/ and tables/ (not in
git). Paths are printed relative to the project folder, never as complete host paths.
Redirects for tests: BG_TXT_DIR, BG_TXT67_DIR, BG_OUT_DIR; BG_SEEDS sets the number of
imputation seeds (default 100).

Definitions used by all scripts:
- conditions: mock and six NS1 mutants named by NS1 residues 38/103/106/184 (R or A; S/I or F/M;
  G or R). RSIG is the PR8 wild type. CPSF30-binding: RFMG, AFMG. G184R: RSIR, RFMR, AFMR.
- comparisons, with the paper's signs: strain - mock for the six strains; RSIX = RSIG - RSIR,
  RFMX = RFMG - RFMR, AFMX = AFMG - AFMR.
- protein IDs: Ensembl version suffixes are removed (ENSP00000000233.5 -> ENSP00000000233);
  REV__ members are ignored when protein groups are matched.
- matching: groups of two tables that share a majority protein ID form connected components;
  1:1 components are compared quantitatively, others are counted as split/merge.
- paper pipeline (Methods; details reverse-engineered from Supplementary Data 5 in q05): no
  decoys, contaminants or site-only groups; >= 2 peptides; detected (Intensity > 0) in >= 5 of 6
  replicates of at least one condition; log2 LFQ; missing values drawn from a normal distribution
  at mean - 1.8 SD with width 0.3 SD of all observed values; per-protein linear model; strain vs
  mock: t test with 35 df; G vs R: z test; Benjamini-Hochberg per comparison; significant:
  |log2FC| > 0.75 and adjusted p < 0.01.
- real and imputed published entries: MaxQuant writes intensities with 5 significant digits, so a
  real published log2 value v satisfies 2^v = m * 10^k with a whole number m; imputed draws do not.
- verdicts: consistent, minor deviation, inconsistent; every question states its rule.
"""
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
import os
import re
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

SCRIPTS_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPTS_DIR.parent
PROJECT_DIR = ANALYSIS_DIR.parent
TXT_DIR = Path(os.environ.get("BG_TXT_DIR", PROJECT_DIR / "results" / "combined" / "txt"))
TXT67_DIR = Path(os.environ.get(
    "BG_TXT67_DIR", PROJECT_DIR.parent / "Bergant_2023_67runs" / "results" / "combined" / "txt"))
OUT_DIR = Path(os.environ.get("BG_OUT_DIR", ANALYSIS_DIR))
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"
ANNOTATION_DIR = ANALYSIS_DIR / "annotation"
DESIGN = PROJECT_DIR / "parameter" / "design.txt"
OUR_MQPAR = PROJECT_DIR / "parameter" / "mqpar.xml"
DEPOSITED_MQPAR = PROJECT_DIR / "parameter" / "deposited" / "mqpar_IAV_FPMS.xml"
SAMPLE_DESCRIPTIONS = PROJECT_DIR / "parameter" / "deposited" / "Sample_descriptions_IAV_FPMS.txt"
DEPOSITED_PG = PROJECT_DIR / "data" / "deposited" / "proteinGroups_IAV_FPMS.txt"
PUBLISHED_XLSX = PROJECT_DIR / "data" / "deposited" / "Supplementary_Data_5_IAV_FPMS.xlsx"
FASTA_HUMAN = PROJECT_DIR / "data" / "ensembl_human_GRCh37-75.fa"
FASTA_VIRUS = {"PR8": PROJECT_DIR / "data" / "Flu_H1N1_PR8.fasta",
               "SC35M": PROJECT_DIR / "data" / "Flu_H7N7_SC35M.fasta"}
GENE_TABLE = ANNOTATION_DIR / "ensp_gene.tsv"
GO_TABLE = ANNOTATION_DIR / "go0051607_defense_response_to_virus.tsv"

N_SEEDS = int(os.environ.get("BG_SEEDS", "100"))
SEED_ROOT = 20230814          # root of np.random.SeedSequence (publication date of the paper)
FC_CUTOFF = 0.75              # paper: |log2 fold change| > 0.75
P_CUTOFF = 0.01               # paper: FDR-adjusted p < 0.01
PAPER_REF = ("Bergant V, Schnepf D, de Andrade Krätzig N, et al. (2023). mRNA 3'UTR lengthening "
             "by alternative polyadenylation attenuates inflammatory responses and correlates "
             "with virulence of Influenza A virus. Nat Commun 14:4906, "
             "doi:10.1038/s41467-023-40469-6")

CONDITIONS = ["mock", "RSIG", "RSIR", "RFMG", "RFMR", "AFMG", "AFMR"]
STRAINS = CONDITIONS[1:]
REPLICATES = range(1, 7)
SAMPLES = [f"{c}_{r}" for c in CONDITIONS for r in REPLICATES]
CONDITION_LABELS = {
    "mock": "mock-infected",
    "RSIG": "R/SI/G, PR8 wild type",
    "RSIR": "R/SI/R, G184R",
    "RFMG": "R/FM/G, CPSF30-binding",
    "RFMR": "R/FM/R, FM with G184R",
    "AFMG": "A/FM/G, dsRNA-binding mutant, CPSF30-binding",
    "AFMR": "A/FM/R, dsRNA-binding mutant with G184R",
}
CPSF_BINDING = ("RFMG", "AFMG")
G184R = ("RSIR", "RFMR", "AFMR")
# name, condition a, condition b, test (fold change = a - b)
COMPARISONS = [(s, s, "mock", "t") for s in STRAINS] + [
    ("RSIX", "RSIG", "RSIR", "z"), ("RFMX", "RFMG", "RFMR", "z"), ("AFMX", "AFMG", "AFMR", "z")]
COMPARISON_NAMES = [c[0] for c in COMPARISONS]
VERDICTS = ["consistent", "minor deviation", "inconsistent"]


def comparison_label(name):
    _, a, b, _ = next(c for c in COMPARISONS if c[0] == name)
    return f"{a} − {b}"


def samples_of(condition):
    return [f"{condition}_{r}" for r in REPLICATES]


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
# one hue per NS1 background, dark for G184, light for R184; mock in gray
CONDITION_COLORS = {"mock": MUTED, "RSIG": "#184f95", "RSIR": "#86b6ef", "RFMG": "#b34a12",
                    "RFMR": "#f2a77e", "AFMG": "#0f7a54", "AFMR": "#86d6b4"}
LAYER_COLORS = {"published": INK, "deposited": SERIES[1], "ours": SERIES[0], "ours-67": SERIES[2]}

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


def condition_legend(ax, conditions=CONDITIONS, **kw):
    from matplotlib.patches import Patch
    handles = [Patch(color=CONDITION_COLORS[c], label=c) for c in conditions]
    ax.legend(handles=handles, **kw)


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
def _read(path, usecols=None):
    return pd.read_csv(path, sep="\t", low_memory=False,
                       usecols=(lambda c: c in usecols) if usecols else None)


def maxquant_version(txt_dir=TXT_DIR):
    """MaxQuant version from the Version row of parameters.txt, e.g. 2.9.0.0."""
    for line in (Path(txt_dir) / "parameters.txt").read_text().splitlines():
        key, _, value = line.partition("\t")
        if key.strip().lower() == "version":
            return value.strip()
    return "version unknown"


def versions():
    """Software versions of this analysis, for answers.md."""
    import platform
    import matplotlib as mpl
    import openpyxl
    import scipy
    return {"Python": platform.python_version(), "pandas": pd.__version__,
            "numpy": np.__version__, "scipy": scipy.__version__,
            "matplotlib": mpl.__version__, "openpyxl": openpyxl.__version__}


@cache
def design():
    """parameter/design.txt with the raw-file key (e.g. '1_1', 'REM_1_4')."""
    d = pd.read_csv(DESIGN, sep="\t", dtype=str)
    d["key"] = d["Name"].str.extract(r"_VBP1_(.+)$")[0]
    if d["key"].isna().any() or set(d["Experiment"]) != set(SAMPLES):
        raise ValueError("design.txt does not list the 42 expected experiments")
    return d.set_index("Experiment").loc[SAMPLES].reset_index()


@cache
def sample_descriptions():
    """The authors' sample description (67 IAV FP-MS runs) with raw-file key and 'used' flag."""
    s = pd.read_csv(SAMPLE_DESCRIPTIONS, sep="\t", dtype=str, encoding="utf-8-sig")
    s = s[s["raw file"].str.lower().str.endswith(".raw") & (s["dataset"] == "IAV FP-MS")].copy()
    s["key"] = s["raw file"].str.extract(r"_VBP1_(.+)\.raw$", flags=re.IGNORECASE)[0]
    s["used"] = s["Comment"].fillna("").str.strip() != "not used"
    return s.reset_index(drop=True)


# ---------------------------------------------------------------- protein IDs
def norm_id(token):
    """Remove an Ensembl version suffix; keep everything else (isoform suffixes, CON__)."""
    m = re.fullmatch(r"(ENSP\d{11})\.\d+", token)
    return m.group(1) if m else token


def id_set(ids, drop_rev=True):
    """Normalised protein IDs of a ';'-separated list as a frozenset."""
    out = []
    for t in str(ids if ids is not None else "").split(";"):
        t = t.strip()
        if not t or t == "nan" or (drop_rev and t.startswith("REV__")):
            continue
        out.append(norm_id(t))
    return frozenset(out)


@cache
def fasta_ids():
    """All protein IDs in the three FASTA files (human ENSP, PR8 and SC35M accessions)."""
    ids = set()
    for path in [FASTA_HUMAN, *FASTA_VIRUS.values()]:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    head = line[1:].split()[0]
                    parts = head.split("|")
                    ids.add(parts[1] if len(parts) >= 2 else head)
    return frozenset(ids)


@cache
def viral_entries():
    """{accession: (entry name, strain, gene)} from the virus FASTA headers, e.g.
    P03496 -> ('NS1_I34A1', 'PR8', 'NS')."""
    out = {}
    for strain, path in FASTA_VIRUS.items():
        for line in open(path):
            if line.startswith(">"):
                head = line[1:].strip()
                parts = head.split()[0].split("|")
                gene = re.search(r"\bGN=(\S+)", head)
                out[parts[1]] = (parts[2], strain, gene.group(1) if gene else "")
    return out


def pr8_proteins():
    """PR8 accessions in a fixed order (polymerase, surface, matrix, NS, accessory)."""
    order = ["PB2", "RDRP", "PB1F2", "PA", "PAX", "HEMA", "NRAM", "NCAP", "M1", "M2", "NS1", "NEP"]
    ent = {acc: name.split("_")[0] for acc, (name, strain, _) in viral_entries().items()
           if strain == "PR8"}
    return sorted(ent, key=lambda a: order.index(ent[a]) if ent[a] in order else 99)


def entry_name(acc):
    return viral_entries().get(acc, (acc, "", ""))[0]


# ---------------------------------------------------------------- MaxQuant tables (layers)
FAMILIES = ["Intensity", "LFQ intensity", "MS/MS count", "Identification type", "iBAQ",
            "Peptides", "Razor + unique peptides", "Unique peptides"]


@dataclass
class Layer:
    """One MaxQuant proteinGroups table. Matrices have the columns SAMPLES (42 runs); `extra`
    lists further runs (the 25 runs the paper does not use) by raw-file key or experiment."""
    name: str
    pg: pd.DataFrame
    columns: dict                      # family -> {sample: column}
    extra: dict = field(default_factory=dict)  # family -> {key: column}
    version: str = ""

    def matrix(self, family, numeric=True):
        m = self.pg[[self.columns[family][s] for s in SAMPLES]].set_axis(SAMPLES, axis=1)
        return m.astype(float) if numeric else m

    def extra_matrix(self, family, numeric=True):
        cols = self.extra.get(family, {})
        m = self.pg[list(cols.values())].set_axis(list(cols), axis=1)
        return m.astype(float) if numeric else m

    @property
    def ids(self):
        return self.pg["Majority protein IDs"].map(id_set)

    @property
    def decoy(self):
        name = "Decoy" if "Decoy" in self.pg.columns else "Reverse"
        return self.pg[name].fillna("").eq("+")

    @property
    def contaminant(self):
        return self.pg["Potential contaminant"].fillna("").eq("+")

    @property
    def site_only(self):
        return self.pg["Only identified by site"].fillna("").eq("+")

    @property
    def clean(self):
        """No decoy, no contaminant, not only identified by site, has a target member."""
        target = self.pg["Protein IDs"].fillna("").map(
            lambda s: any(t and not t.startswith(("CON__", "REV__")) for t in str(s).split(";")))
        return target & ~self.decoy & ~self.contaminant & ~self.site_only


def _layer_columns(pg, sample_column):
    """family -> {sample: column}; sample_column(family, sample) returns the column or None."""
    out = {}
    for fam in FAMILIES:
        cols = {}
        for s in SAMPLES:
            c = sample_column(fam, s)
            if c is None:
                raise KeyError(f"no column for {fam!r} and {s}")
            cols[s] = c
        out[fam] = cols
    return out


@cache
def ours(txt_dir=None):
    """Our proteinGroups.txt (42 runs, experiment names = SAMPLES)."""
    txt = Path(txt_dir) if txt_dir else TXT_DIR
    pg = _read(txt / "proteinGroups.txt").copy()
    cols = _layer_columns(pg, lambda fam, s: f"{fam} {s}" if f"{fam} {s}" in pg.columns else None)
    known = {f"{fam} {s}" for fam in FAMILIES for s in SAMPLES}
    extra = {}
    for fam in FAMILIES:
        hits = {c[len(fam) + 1:]: c for c in pg.columns
                if c.startswith(fam + " ") and c not in known
                and re.fullmatch(r"[A-Za-z0-9]+_\d+(rem)?", c[len(fam) + 1:])}
        if hits:
            extra[fam] = hits
    name = "ours" if txt == TXT_DIR else "ours-67"
    return Layer(name, pg, cols, extra, maxquant_version(txt))


def has_ours67():
    return (TXT67_DIR / "proteinGroups.txt").exists()


@cache
def deposited():
    """The authors' proteinGroups (67 runs); columns renamed to our sample names via the
    raw-file keys of design.txt (the REM columns have no trailing underscore)."""
    pg = _read(DEPOSITED_PG).copy()
    key_of = dict(zip(design()["Experiment"], design()["key"]))

    def column(fam, s):
        pat = re.compile(rf"{re.escape(fam)} A549Easy12_VBP1_{re.escape(key_of[s])}_?")
        hits = [c for c in pg.columns if pat.fullmatch(c)]
        return hits[0] if len(hits) == 1 else None

    cols = _layer_columns(pg, column)
    extra = {}
    unused = sample_descriptions()[~sample_descriptions()["used"]]["key"]
    for fam in FAMILIES:
        e = {}
        for key in unused:
            pat = re.compile(rf"{re.escape(fam)} A549Easy12_VBP1_{re.escape(key)}_?")
            hits = [c for c in pg.columns if pat.fullmatch(c)]
            if len(hits) == 1:
                e[key] = hits[0]
        extra[fam] = e
    return Layer("deposited", pg, cols, extra, "1.5.6.2")


def layers(include_67=True):
    """The MaxQuant layers available: ours, deposited, and ours-67 when Part B has run."""
    out = {"ours": ours(), "deposited": deposited()}
    if include_67 and has_ours67():
        out["ours-67"] = ours(str(TXT67_DIR))
    return out


# ---------------------------------------------------------------- published table
def is_maxquant_value(v):
    """True where a log2 value comes from a MaxQuant intensity (5 significant digits)."""
    v = np.asarray(v, dtype=float)
    x = np.power(2.0, v)
    k = np.floor(np.log10(x)) - 4
    m = x / np.power(10.0, k)
    return np.abs(m - np.round(m)) < 1e-4


@dataclass
class Published:
    df: pd.DataFrame          # the sheet, columns renamed (MOCK -> mock)
    matrix: pd.DataFrame      # log2 LFQ incl. imputed values, columns SAMPLES
    real: pd.DataFrame        # True for MaxQuant values, False for imputed draws
    ids: pd.Series            # frozensets of normalised majority IDs
    genes: pd.Series          # gene names as given
    fc: pd.DataFrame          # comparison -> log2 fold change
    p: pd.DataFrame
    q: pd.DataFrame


@cache
def published():
    df = pd.read_excel(PUBLISHED_XLSX, header=1)
    df.columns = [re.sub(r"^MOCK", "mock", c) if isinstance(c, str) else c for c in df.columns]
    matrix = df[SAMPLES].astype(float)
    real = pd.DataFrame(is_maxquant_value(matrix.values), index=matrix.index, columns=SAMPLES)
    ids = df["Majority.protein.IDs"].map(id_set)
    genes = df["gene_names"].fillna("").map(str)
    fc = pd.DataFrame({c: df[c].astype(float) for c in COMPARISON_NAMES})
    p = pd.DataFrame({c: df[f"{c}.p"].astype(float) for c in COMPARISON_NAMES})
    q = pd.DataFrame({c: df[f"{c}.p.adj"].astype(float) for c in COMPARISON_NAMES})
    return Published(df, matrix, real, ids, genes, fc, p, q)


def published_hits():
    """Published calls per comparison: +1 up, -1 down, 0 not significant."""
    pub = published()
    sig = (pub.fc.abs() > FC_CUTOFF) & (pub.q < P_CUTOFF)
    return np.sign(pub.fc).where(sig, 0).astype(int)


# ---------------------------------------------------------------- matching of protein groups
def match_groups(ids_a, ids_b):
    """Connected components of groups that share a protein ID.

    ids_a, ids_b: sequences of frozensets (row order of the two tables).
    Returns (components, pairs):
    - components: DataFrame with columns a (list of positions), b, category
      (identical, 1:1 overlap, split/merge, only A, only B)
    - pairs: DataFrame a, b, identical for the 1:1 components (positions)."""
    ids_a, ids_b = list(ids_a), list(ids_b)
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    owner = {}
    for tag, sets in (("a", ids_a), ("b", ids_b)):
        for i, s in enumerate(sets):
            if not s:          # decoy groups (only REV__ members) take no part
                continue
            node = (tag, i)
            find(node)
            for t in s:
                if t in owner:
                    ra, rb = find(node), find(owner[t])
                    if ra != rb:
                        parent[ra] = rb
                else:
                    owner[t] = node
    comp = {}
    for node in list(parent):
        comp.setdefault(find(node), ([], []))[0 if node[0] == "a" else 1].append(node[1])
    rows = []
    for a, b in comp.values():
        a, b = sorted(a), sorted(b)
        if a and b:
            if len(a) == 1 and len(b) == 1:
                cat = "identical" if ids_a[a[0]] == ids_b[b[0]] else "1:1 overlap"
            else:
                cat = "split/merge"
        else:
            cat = "only A" if a else "only B"
        rows.append((a, b, cat))
    components = pd.DataFrame(rows, columns=["a", "b", "category"])
    one = components[components["category"].isin(["identical", "1:1 overlap"])]
    pairs = pd.DataFrame({"a": [x[0] for x in one["a"]], "b": [x[0] for x in one["b"]],
                          "identical": (one["category"] == "identical").values})
    return components, pairs.sort_values("a").reset_index(drop=True)


def category_counts(components):
    cats = ["identical", "1:1 overlap", "split/merge", "only A", "only B"]
    c = components["category"].value_counts()
    out = {k: int(c.get(k, 0)) for k in cats}
    sm = components[components["category"] == "split/merge"]
    out["split/merge groups A"] = int(sum(len(x) for x in sm["a"]))
    out["split/merge groups B"] = int(sum(len(x) for x in sm["b"]))
    return out


# ---------------------------------------------------------------- paper pipeline
def paper_filter(layer, detect="Intensity"):
    """Rows kept by the paper's filter (clean, >= 2 peptides, detected in >= 5 of 6 replicates
    of at least one condition); detect is 'Intensity' (default) or 'LFQ intensity'."""
    det = layer.matrix(detect) > 0
    enough = np.column_stack([det[samples_of(c)].sum(axis=1) >= 5 for c in CONDITIONS]).any(axis=1)
    return layer.clean & (layer.pg["Peptides"].astype(float) >= 2) & enough


def impute(X, rng):
    """Whole-matrix imputation as in the paper: normal draws at mean - 1.8 SD, width 0.3 SD of
    all observed values. X: 2-D array with NaN for missing values."""
    obs = X[~np.isnan(X)]
    mu, sd = obs.mean(), obs.std()
    Y = X.copy()
    miss = np.isnan(Y)
    Y[miss] = rng.normal(mu - 1.8 * sd, 0.3 * sd, int(miss.sum()))
    return Y


def bh(p):
    """Benjamini-Hochberg adjusted p values (1-D)."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def paper_stats(Y):
    """Fold changes, p and adjusted p for the 9 comparisons. Y: n x 42 array, columns SAMPLES.
    Strain vs mock: t test of the linear-model coefficient (35 df); G vs R: z test."""
    n = Y.shape[0]
    Y3 = Y.reshape(n, len(CONDITIONS), len(REPLICATES))
    means = Y3.mean(axis=2)
    df_res = Y.shape[1] - len(CONDITIONS)
    s2 = ((Y3 - means[:, :, None]) ** 2).sum(axis=(1, 2)) / df_res
    se = np.sqrt(s2 * 2 / len(REPLICATES))
    fc, p, q = {}, {}, {}
    for name, a, b, test in COMPARISONS:
        d = means[:, CONDITIONS.index(a)] - means[:, CONDITIONS.index(b)]
        z = np.abs(d / se)
        pv = 2 * (stats.t.sf(z, df_res) if test == "t" else stats.norm.sf(z))
        fc[name], p[name], q[name] = d, pv, bh(pv)
    return fc, p, q


def calls(fc, q):
    """+1 up, -1 down, 0 not significant (paper thresholds); arrays of equal shape."""
    sig = (np.abs(fc) > FC_CUTOFF) & (q < P_CUTOFF)
    return np.where(sig, np.sign(fc), 0).astype(np.int8)


def seed_sequences(n=None):
    return np.random.SeedSequence(SEED_ROOT).spawn(n or N_SEEDS)


@dataclass
class PipelineResult:
    """Paper pipeline over n seeds. rows: positions in the source table (sorted by majority IDs);
    observed: n_rows x 42; fc and calls: seeds x n_rows x 9 (COMPARISON_NAMES order)."""
    name: str
    rows: np.ndarray
    observed: np.ndarray
    log2: np.ndarray
    fc: np.ndarray
    calls: np.ndarray

    def counts(self):
        """seeds x 9 x 2 (up, down)."""
        return np.stack([(self.calls == 1).sum(axis=1), (self.calls == -1).sum(axis=1)], axis=2)

    def imputed_share(self):
        return 1 - self.observed.mean(axis=1)


def _run(name, rows, X, n_seeds):
    observed = ~np.isnan(X)
    fcs, cls = [], []
    for ss in seed_sequences(n_seeds):
        Y = impute(X, np.random.default_rng(ss))
        fc, _, q = paper_stats(Y)
        f = np.column_stack([fc[c] for c in COMPARISON_NAMES])
        fcs.append(f.astype(np.float32))
        cls.append(calls(f, np.column_stack([q[c] for c in COMPARISON_NAMES])))
    return PipelineResult(name, rows, observed, X, np.stack(fcs), np.stack(cls))


@cache
def pipeline(layer_name, detect="Intensity", n_seeds=None):
    """Paper pipeline on a MaxQuant layer ('ours', 'deposited', 'ours-67')."""
    layer = layers()[layer_name]
    keep = paper_filter(layer, detect).values
    rows = np.flatnonzero(keep)
    order = np.argsort(layer.pg["Majority protein IDs"].fillna("").map(str).values[rows],
                       kind="stable")
    rows = rows[order]
    lfq = layer.matrix("LFQ intensity").values[rows]
    X = np.log2(np.where(lfq > 0, lfq, np.nan))
    return _run(layer_name, rows, X, n_seeds or N_SEEDS)


@cache
def published_reimputed(n_seeds=None):
    """The published matrix with its imputed entries removed, re-imputed with our seeds."""
    pub = published()
    rows = np.arange(len(pub.matrix))
    X = np.where(pub.real.values, pub.matrix.values, np.nan)
    return _run("published re-imputed", rows, X, n_seeds or N_SEEDS)


# ---------------------------------------------------------------- statistics helpers
def rma_slope(x, y):
    """Reduced major axis slope of y on x: sign(r) * SD(y) / SD(x)."""
    r = np.corrcoef(x, y)[0, 1]
    return float(np.sign(r) * np.std(y) / np.std(x))


def spearman_brown(r, factor=2):
    return factor * r / (1 + (factor - 1) * r)


def mad(x):
    x = np.asarray(x, dtype=float)
    return float(np.median(np.abs(x - np.median(x))))


def worst(verdicts):
    vs = [v for v in verdicts if v]
    return max(vs, key=VERDICTS.index) if vs else "consistent"


def band(value, lo, hi):
    return lo <= value <= hi


# ---------------------------------------------------------------- annotation (genes, GO)
@cache
def gene_table():
    """ENSP -> gene symbol, HGNC ID, Swiss-Prot accessions (annotation/ensp_gene.tsv)."""
    if not GENE_TABLE.exists():
        raise FileNotFoundError(f"{rel(GENE_TABLE)} missing; run scripts/make_annotation.py")
    t = pd.read_csv(GENE_TABLE, sep="\t", comment="#", dtype=str).fillna("")
    return t.set_index("ensp")


@cache
def go_defense():
    """UniProt accessions and symbols annotated to GO:0051607 or a descendant."""
    if not GO_TABLE.exists():
        raise FileNotFoundError(f"{rel(GO_TABLE)} missing; run scripts/make_annotation.py")
    t = pd.read_csv(GO_TABLE, sep="\t", comment="#", dtype=str).fillna("")
    return frozenset(t["uniprot"]), frozenset(t["symbol"])


def group_genes(ids):
    """Gene symbols of a protein group (frozenset of IDs), sorted and unique."""
    gt = gene_table()
    out = set()
    for t in ids:
        if t in gt.index and gt.at[t, "gene"]:
            out.add(gt.at[t, "gene"])
        elif t in viral_entries():
            out.add(entry_name(t))
    return sorted(out)


def group_in_go(ids):
    """True if a member of the group maps to a UniProt accession of the GO:0051607 set."""
    acc, _ = go_defense()
    gt = gene_table()
    for t in ids:
        if t in gt.index:
            if any(a in acc for a in gt.at[t, "uniprot"].split(";") if a):
                return True
    return False


# ---------------------------------------------------------------- formatting, answers
def n(x):
    return f"{int(x):,}"


def pct(x, digits=1):
    return f"{100 * x:.{digits}f}%"


def f2(x, digits=2):
    return "–" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{digits}f}"


def md_table(rows, header):
    def cell(c):
        return str(c).replace("|", "\\|")
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(cell(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


@dataclass
class Answer:
    """str() gives the Markdown section body. `short` and `verdict` feed the summary tables."""
    question: str
    text: str
    figure: object
    short: str
    verdict: str
    criterion: str
    metrics: dict = field(default_factory=dict)

    def __str__(self):
        figs = self.figure if isinstance(self.figure, (list, tuple)) else [self.figure]
        body = f"**Verdict: {self.verdict}.** {self.criterion}\n\n{self.text}\n"
        for fig in figs:
            if fig:
                body += f"\n![{self.question}]({fig})\n"
        return body


def answer(question, text, figure, short, verdict, criterion, metrics=None):
    if verdict not in VERDICTS:
        raise ValueError(f"unknown verdict {verdict!r}")
    return Answer(question, text, figure, short, verdict, criterion, metrics or {})


# ---------------------------------------------------------------- self test
def selftest():
    """Checks of the building blocks on synthetic data; run_all.py calls it first."""
    # ID normalisation
    assert norm_id("ENSP00000000233.5") == "ENSP00000000233"
    assert norm_id("CON__P01044-1") == "CON__P01044-1"
    assert id_set("ENSP00000000233.5;REV__ENSP00000000412;P03496") == frozenset(
        {"ENSP00000000233", "P03496"})
    # matching: identical, overlap, subset, split, only A, only B
    A = [frozenset("ab"), frozenset("cd"), frozenset("e"), frozenset("fg"), frozenset("x")]
    B = [frozenset("ab"), frozenset("c"), frozenset("ez"), frozenset("f"), frozenset("g"),
         frozenset("y")]
    comp, pairs = match_groups(A, B)
    cats = dict(zip(comp["a"].map(tuple), comp["category"]))
    assert cats[(0,)] == "identical" and cats[(1,)] == "1:1 overlap" and cats[(2,)] == "1:1 overlap"
    assert cats[(3,)] == "split/merge" and cats[(4,)] == "only A"
    assert (comp["category"] == "only B").sum() == 1 and len(pairs) == 3
    # chained overlaps form one component
    comp, _ = match_groups([frozenset("ab"), frozenset("cd")], [frozenset("bc")])
    assert len(comp) == 1 and comp["category"].iloc[0] == "split/merge"
    # Benjamini-Hochberg against a direct computation
    p = np.array([0.01, 0.04, 0.03, 0.2])
    ref = np.array([0.04, 0.0533333, 0.0533333, 0.2])
    assert np.allclose(bh(p), ref, atol=1e-6)
    # MaxQuant rounding rule
    assert is_maxquant_value(np.log2([51139000.0, 2780200000.0, 12345.0])).all()
    rng = np.random.default_rng(1)
    assert is_maxquant_value(rng.normal(24, 0.8, 10000)).mean() < 0.01
    # paper statistics: equal groups give fold change 0; a shifted group gives the shift
    Y = np.tile(np.arange(42, dtype=float) % 6, (2, 1))
    Y[1, 6:12] += 3.0
    fc, p_, q_ = paper_stats(Y)
    assert abs(fc["RSIG"][0]) < 1e-12 and abs(fc["RSIG"][1] - 3.0) < 1e-12
    assert abs(fc["RSIX"][1] - 3.0) < 1e-12
    print("selftest passed")
