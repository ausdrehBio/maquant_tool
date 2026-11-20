r"""Evidence for DPR sequences in the Bergant_2023_adapted search.

Question: which peptides identified by MaxQuant can only come from the DPR entries of
parameter/additional_proteins.fasta (PR8_3frame_uniq_*), and how many of them survive filters of
increasing strictness? Every filter is reported with the peptides it removes, so a hit lost at a
strict level stays visible.

Inputs (read only):
- results/combined/txt/: peptides.txt, evidence.txt, proteinGroups.txt, parameters.txt
- parameter/additional_proteins.fasta, parameter/design.txt, results/run_info.txt
- the canonical FASTA files (human Ensembl GRCh37, PR8, SC35M) from data/ of this project or of
  ../Bergant_2023_67runs or ../Bergant_2023; their SHA-256 must match results/run_info.txt
- MaxQuant's contaminants.fasta from MaxQuant_v.2.9.0/ at the repository root
- optional: ../Bergant_2023/results/combined/txt/proteinGroups.txt (the search without the
  additional file), for the side effects on canonical viral proteins
Outputs: analysis_of_results/answers.md (in git); figures/ and tables/ (not in git), including
tables/dpr_evidence.tsv, which doubles as a cache of the evidence.txt rows needed (evidence.txt
is several GB; the first run reads it in a few minutes).
Paths are printed relative to the project folder. Redirects: DPR_TXT_DIR, DPR_OUT_DIR,
DPR_BASELINE_TXT_DIR.

Definitions:
- DPR entry: an entry PR8_3frame_uniq_* of additional_proteins.fasta. Its VIRAL_* (Mount Sinai
  PR8) and CTRL_* (human) entries count as canonical sequences here.
- DPR-specific peptide: a peptide at MaxQuant's 1% FDR that is a tryptic peptide of at least one
  DPR entry (Trypsin/P, up to 2 missed cleavages, protein N-terminal Met removal) and occurs
  nowhere in a canonical sequence (human, PR8, SC35M, contaminants, VIRAL_*, CTRL_*; I = L).
  The mapping comes from the FASTA files, not from peptides.txt: peptides.txt lists only the
  proteins kept in protein groups, and MaxQuant drops canonical proteins whose peptides all went
  to longer DPR entries (M1 in this search), so their peptides look DPR-specific there.
- decoys: MaxQuant's decoy sequences ("revert": reversed, then each K/R swapped with the residue
  before it). A decoy peptide is DPR-specific if it occurs in a reverted DPR entry and in no
  reverted canonical sequence. Each filter treats these decoys exactly like the targets, so
  decoys / targets estimates the share of false hits left at each level.
- origin (a description, not a filter), from the PR8 proteins (UniProt and VIRAL_* entries):
  junction: the peptide begins in one viral sequence and ends in another, or in a distant part of
  the same one, each side >= 5 residues, together covering the peptide (a deletion junction);
  frame continuation: one end matches >= 5 residues of a viral sequence, the rest is new (a
  junction into another reading frame, or an extension);
  substitution: the peptide differs from a viral sequence in one residue (more likely a sequence
  variant of the virus stock than a DPR);
  no viral anchor: neither end matches >= 5 residues of a viral sequence.
- runs and conditions come from evidence.txt. "MS/MS" means the evidence types MSMS,
  MULTI-MSMS and MULTI-SECPEP; "main search" means MSMS and MULTI-MSMS (MULTI-SECPEP is the
  second-peptide search of chimeric spectra; MULTI-MATCH and MULTI-MATCH-MSMS are transfers by
  match between runs). RSIG_3rem, a second measurement of RSIG_3, counts as an RSIG run.
- filters, applied in the order of FILTERS below; level Lk keeps the peptides that pass the
  first k filters. A peptide's "lost at" is the first filter it fails; "fails only" names the
  filter when it is the only one the peptide fails.

Usage (from any directory):  python3 <path to>/scripts/dpr_evidence.py
"""
import hashlib
import os
import platform
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# ---------------------------------------------------------------- paths
SCRIPTS_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPTS_DIR.parent
PROJECT_DIR = ANALYSIS_DIR.parent
REPO_DIR = PROJECT_DIR.parent.parent
TXT_DIR = Path(os.environ.get("DPR_TXT_DIR", PROJECT_DIR / "results" / "combined" / "txt"))
OUT_DIR = Path(os.environ.get("DPR_OUT_DIR", ANALYSIS_DIR))
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"
BASELINE_TXT = Path(os.environ.get(
    "DPR_BASELINE_TXT_DIR", PROJECT_DIR.parent / "Bergant_2023" / "results" / "combined" / "txt"))
OWN_FASTA = PROJECT_DIR / "parameter" / "additional_proteins.fasta"
DESIGN = PROJECT_DIR / "parameter" / "design.txt"
RUN_INFO = PROJECT_DIR / "results" / "run_info.txt"
DATA_DIRS = [PROJECT_DIR / "data", PROJECT_DIR.parent / "Bergant_2023_67runs" / "data",
             PROJECT_DIR.parent / "Bergant_2023" / "data"]
DATABASES = {"human": "ensembl_human_GRCh37-75.fa", "PR8": "Flu_H1N1_PR8.fasta",
             "SC35M": "Flu_H7N7_SC35M.fasta"}
CONTAMINANTS = REPO_DIR / "MaxQuant_v.2.9.0" / "bin" / "conf" / "contaminants.fasta"

# ---------------------------------------------------------------- settings
PEP_MAX = 0.01          # best main-search PEP
MIN_RUNS = 3            # runs with a main-search identification
ANCHOR = 5              # residues that must match a viral sequence to count as an anchor
LOOKALIKE_PPM = 10      # mass window of the isobaric-twin check
MSMS_TYPES = {"MSMS", "MULTI-MSMS", "MULTI-SECPEP"}
MAIN_TYPES = {"MSMS", "MULTI-MSMS"}
MOCK = {"mock", "mock2"}
SC35M = {"SC35M", "SC35MdNS1"}
FILTERS = [  # key, short name, rule
    ("msms", "MS/MS",
     "identified by an MS/MS spectrum in at least one run, not only transferred by match between "
     "runs"),
    ("main", "main search",
     "at least one identification from the main search, not only from the second-peptide search "
     "of chimeric spectra"),
    ("pep", f"PEP < {PEP_MAX}", f"best main-search PEP below {PEP_MAX}"),
    ("no_mock", "not in mock", "no MS/MS identification in a mock or mock2 run"),
    ("not_sc35m_only", "not SC35M-only",
     "MS/MS identifications not only in SC35M or SC35MdNS1 runs (H7N7-infected cells)"),
    ("no_isobaric", "no isobaric twin",
     f"no canonical sequence of the same mass (within {LOOKALIKE_PPM} ppm) that shares at least "
     "half of the residues as prefix or suffix (any part of a viral protein, or an identified "
     "canonical peptide)"),
    ("no_deamid_twin", "not a deamidated canonical",
     "at least one main-search identification whose modified form is not a canonical sequence "
     "(deamidated N has the composition of D, deamidated Q that of E)"),
    ("runs", f">= {MIN_RUNS} runs", f"main-search identifications in at least {MIN_RUNS} runs"),
    ("no_deamid", "one spectrum without deamidation",
     "at least one main-search identification without deamidation (strictest: the new residues "
     "must not depend on deamidation to fit the spectrum)"),
]
ORIGINS = ["junction", "frame continuation", "no viral anchor", "substitution"]
UNIPROT_NAMES = {"PB2": "PB2", "RDRP": "PB1", "PB1F2": "PB1-F2", "PA": "PA", "PAX": "PA-X",
                 "HEMA": "HA", "NCAP": "NP", "NRAM": "NA", "M1": "M1", "M2": "M2", "NS1": "NS1",
                 "NEP": "NEP"}
MONO = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768,
            C=160.03065, L=113.08406, I=113.08406, N=114.04293, D=115.02694, Q=128.05858,
            K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841, R=156.10111,
            Y=163.06333, W=186.07931, U=150.95364, O=237.14773)
WATER = 18.01056

# ---------------------------------------------------------------- plot style (as Bergant_2023)
SURFACE, INK, INK2, MUTED, GRID, AXIS = ("#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9",
                                         "#c3c2b7")
ORIGIN_COLORS = dict(zip(ORIGINS, ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")))
RUN_RAMP = ["#f0efec", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95",
            "#0d366b"]  # 0 runs (neutral), then blue steps 100 to 700
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "savefig.dpi": 200, "savefig.bbox": "tight", "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlecolor": INK, "xtick.color": AXIS, "ytick.color": AXIS, "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2, "legend.frameon": False, "legend.fontsize": 9,
})


# ---------------------------------------------------------------- helpers
def fail(msg):
    sys.exit(f"ERROR: {msg}")


def rel(path):
    """Path relative to the project folder (or the repository) for printing, never a host path."""
    p = Path(path).resolve()
    for base, prefix in ((PROJECT_DIR, ""), (REPO_DIR, "<repo>/")):
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


def save_fig(fig, name):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"figure written: {rel(path)}")
    return path.relative_to(OUT_DIR).as_posix()


def save_table(df, name):
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / f"{name}.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(f"table written: {rel(path)}")
    return path.relative_to(OUT_DIR).as_posix()


def L(s):
    return s.replace("I", "L")


def decoy_column(df):
    return "Decoy" if "Decoy" in df.columns else "Reverse"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path, rule, prefix=""):
    """{id: sequence}, {id: header}; rule 'pipe' = MaxQuant's >[^|]*\\|(.*?)\\|, 'first' = first word."""
    seqs, headers, key = {}, {}, None
    for line in open(path):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if rule == "pipe":
                m = re.match(r">[^|]*\|(.*?)\|", line)
                key = prefix + m.group(1) if m else None
            else:
                key = prefix + line[1:].split()[0]
            if key in seqs:          # MaxQuant keeps the first entry of a duplicated ID
                key = None
                continue
            if key:
                seqs[key], headers[key] = [], line[1:]
        elif key:
            seqs[key].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}, headers


def revert(seq):
    """MaxQuant's decoy sequence (decoy mode 'revert', special amino acids K and R)."""
    r = list(seq[::-1])
    for i in range(1, len(r)):
        if r[i] in "KR":
            r[i - 1], r[i] = r[i], r[i - 1]
    return "".join(r)


def tryptic_map(seqs, wanted):
    """{peptide (I = L): set of entry ids} for the peptides in `wanted` that are tryptic peptides
    of an entry (Trypsin/P, up to 2 missed cleavages, with and without the N-terminal Met)."""
    hits = defaultdict(set)
    for k, seq in seqs.items():
        s = L(seq)
        cuts = sorted({0, len(s)} | {m.end() for m in re.finditer("[KR]", s)})
        for i in range(len(cuts) - 1):
            starts = [0, 1] if cuts[i] == 0 and s.startswith("M") else [cuts[i]]
            for j in range(i + 1, min(i + 4, len(cuts))):
                for a in starts:
                    p = s[a:cuts[j]]
                    if p in wanted:
                        hits[p].add(k)
    return hits


def mass(seq):
    return sum(MONO.get(c, 0.0) for c in seq) + WATER


def deamidated_form(modified):
    """Sequence of a modified peptide with deamidated N as D and deamidated Q as E (identical
    compositions); other modifications are dropped."""
    s = modified.replace("N(Deamidation (NQ))", "D").replace("Q(Deamidation (NQ))", "E")
    s = re.sub(r"\([^()]*\([^()]*\)\)", "", s)
    s = re.sub(r"\([^()]*\)", "", s)
    return s.strip("_")


def condition_of(experiment):
    return re.sub(r"_(\d+|3rem)$", "", experiment)


# ---------------------------------------------------------------- inputs
def find_database(name):
    for folder in DATA_DIRS:
        if (folder / name).exists():
            return folder / name
    fail(f"{name} not found in " + ", ".join(rel(d) for d in DATA_DIRS)
         + "; run get_data_Bergant_2023_adapted.sh")


def check_checksums(paths):
    """Compare the FASTA files read here with the checksums of the run (results/run_info.txt)."""
    if not RUN_INFO.exists():
        return ["results/run_info.txt is missing, so the FASTA files were not checked against the "
                "run."]
    recorded = {}
    for line in open(RUN_INFO):
        parts = line.split()
        if len(parts) == 2 and re.fullmatch(r"[0-9a-f]{64}", parts[0]):
            recorded[Path(parts[1]).name] = parts[0]
    notes = []
    for p in paths:
        want = recorded.get(p.name)
        if want is None:
            notes.append(f"{p.name}: not listed in results/run_info.txt")
        elif sha256(p) != want:
            fail(f"{rel(p)} differs from the file the run used (results/run_info.txt); the mapping "
                 "would not match the search")
    return notes


class Databases:
    def __init__(self):
        self.paths = {k: find_database(v) for k, v in DATABASES.items()}
        self.notes = check_checksums([*self.paths.values(), OWN_FASTA])
        if not CONTAMINANTS.exists():
            fail(f"{rel(CONTAMINANTS)} missing (MaxQuant_v.2.9.0/ belongs at the repository root)")
        canonical = {}
        for key, path in self.paths.items():
            canonical.update(read_fasta(path, "pipe")[0])
        self.viral = {}                                   # display name -> sequence (I = L)
        pr8, pr8_headers = read_fasta(self.paths["PR8"], "pipe")
        for acc, header in pr8_headers.items():
            mnemonic = header.split("|")[2].split("_")[0] if header.count("|") >= 2 else acc
            self.viral[UNIPROT_NAMES.get(mnemonic, mnemonic)] = L(pr8[acc])
        canonical.update(read_fasta(CONTAMINANTS, "first", "CON__")[0])
        own, self.own_headers = read_fasta(OWN_FASTA, "first")
        self.dpr = {k: v for k, v in own.items() if k.startswith("PR8_3frame")}
        self.other_own = {k: v for k, v in own.items() if not k.startswith("PR8_3frame")}
        for k, v in self.other_own.items():
            if k.startswith("VIRAL_"):
                name = k.split("_")[1]
                if L(v) != self.viral.get(name):
                    self.viral[f"{name} (Mount Sinai)"] = L(v)
        self.core_text = "|".join(L(s) for s in canonical.values())          # without VIRAL/CTRL
        self.canon_text = "|".join([self.core_text] + [L(s) for s in self.other_own.values()])
        self.canon_rev_text = "|".join(L(revert(s)) for s in
                                       list(canonical.values()) + list(self.other_own.values()))
        self.viral_own_text = "|".join(L(v) for k, v in self.other_own.items()
                                       if k.startswith("VIRAL_"))
        self.n_canonical = len(canonical) + len(self.other_own)


def dpr_fields(header):
    return {k: int(v) for k, v in re.findall(r"(len|n_orfs|n_delvg|span_any)=(\d+)", header)}


# ---------------------------------------------------------------- sequence descriptions
class Describer:
    def __init__(self, db, canonical_peptides):
        self.viral = db.viral
        self.index = defaultdict(list)      # rounded mass -> (mass, sequence, source)
        for name, v in self.viral.items():
            for i in range(len(v)):
                m = WATER
                for j in range(i, min(i + 45, len(v))):
                    m += MONO.get(v[j], 0.0)
                    self.index[round(m)].append((m, v[i:j + 1], name))
        for s in canonical_peptides:
            m = mass(s)
            self.index[round(m)].append((m, L(s), "identified peptide"))

    def _longest(self, p, side):
        for n in range(len(p), 0, -1):
            part = p[:n] if side == "start" else p[-n:]
            hits = [(name, v.find(part)) for name, v in self.viral.items() if part in v]
            if hits:
                return n, hits
        return 0, []

    @staticmethod
    def _where(hits, n):
        seen, out = set(), []
        for name, start in hits:
            base = name.replace(" (Mount Sinai)", "")
            if base not in seen:
                seen.add(base)
                out.append(f"{base} {start + 1}-{start + n}")
        return "/".join(out[:3])

    def origin(self, p):
        p = L(p)
        n1, h1 = self._longest(p, "start")
        n2, h2 = self._longest(p, "end")
        if n1 >= ANCHOR and n2 >= ANCHOR and n1 + n2 >= len(p):
            return "junction", f"{self._where(h1, n1)} | {self._where(h2, n2)}"
        for name, v in self.viral.items():
            for i in range(len(v) - len(p) + 1):
                diff = []
                for j in range(len(p)):
                    if v[i + j] != p[j]:
                        diff.append(j)
                        if len(diff) > 1:
                            break
                if len(diff) == 1:
                    j = diff[0]
                    base = name.replace(" (Mount Sinai)", "")
                    return "substitution", f"{base} {v[i + j]}{i + j + 1}{p[j]}"
        if max(n1, n2) >= ANCHOR:
            if n1 >= n2:
                return "frame continuation", f"{self._where(h1, n1)} + {len(p) - n1} new"
            return "frame continuation", f"{len(p) - n2} new + {self._where(h2, n2)}"
        return "no viral anchor", f"longest viral match {max(n1, n2)} of {len(p)} residues"

    def isobaric_twin(self, p):
        p = L(p)
        m = mass(p)
        best = None
        for key in (round(m) - 1, round(m), round(m) + 1):
            for mm, s, src in self.index.get(key, []):
                if s == p or abs(mm - m) / m * 1e6 > LOOKALIKE_PPM:
                    continue
                k = min(len(s), len(p))
                pre = next((i for i in range(k) if s[i] != p[i]), k)
                suf = next((i for i in range(k) if s[-1 - i] != p[-1 - i]), k)
                share = max(pre, suf) / len(p)
                if share >= 0.5 and (best is None or share > best[0]):
                    best = (share, f"{s} ({src}, {(m - mm) / m * 1e6:+.1f} ppm)")
        return best[1] if best else ""


# ---------------------------------------------------------------- evidence
EV_COLS = ["Raw file", "Experiment", "Sequence", "Modifications", "Modified sequence", "Decoy",
           "Charge", "Mass error [ppm]", "Score", "PEP", "Type", "MS/MS count",
           "MS/MS scan number", "Peptide ID", "id"]


def evidence_for(peptide_ids):
    """evidence.txt rows of the given peptides; cached in tables/dpr_evidence.tsv."""
    path = TXT_DIR / "evidence.txt"
    if not path.exists():
        fail(f"{rel(path)} missing")
    st = path.stat()
    ids = sorted(int(i) for i in peptide_ids)
    key = hashlib.sha1(f"{st.st_size}:{st.st_mtime_ns}:{ids}".encode()).hexdigest()
    cache = TABLE_DIR / "dpr_evidence.tsv"
    if cache.exists():
        with open(cache) as fh:
            if fh.readline().strip() == f"# cache key {key}":
                print(f"evidence rows from the cache {rel(cache)}")
                return pd.read_csv(cache, sep="\t", skiprows=1, low_memory=False)
    print(f"reading {rel(path)} ({st.st_size / 1e9:.1f} GB), a few minutes ...")
    wanted = set(ids)
    parts = [chunk[chunk["Peptide ID"].isin(wanted)]
             for chunk in pd.read_csv(path, sep="\t", usecols=lambda c: c in EV_COLS,
                                      chunksize=500_000, low_memory=False)]
    ev = pd.concat(parts, ignore_index=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache, "w") as fh:
        fh.write(f"# cache key {key}\n")
        ev.to_csv(fh, sep="\t", index=False)
    print(f"evidence rows cached: {rel(cache)}")
    return ev


def peptide_filters(cands, ev, db, describer, decoy):
    """One row per candidate peptide with its evidence summary and the filter results."""
    by_peptide = {k: g for k, g in ev.groupby("Peptide ID")}
    empty = ev.iloc[0:0]
    rows = []
    for _, r in cands.iterrows():
        e = by_peptide.get(r["id"], empty)
        ms = e[e["Type"].isin(MSMS_TYPES)]
        main = e[e["Type"].isin(MAIN_TYPES)]
        ms_conditions = {condition_of(x) for x in ms["Experiment"]}
        main_runs = main.groupby("Experiment").size()
        twin_free = any(L(deamidated_form(m)) not in db.canon_text for m in main["Modified sequence"])
        rows.append({
            "id": r["id"], "Sequence": r["Sequence"],
            "best PEP": main["PEP"].min() if len(main) else np.nan,
            "best score": main["Score"].max() if len(main) else np.nan,
            "MS/MS spectra": int(len(ms)), "main-search spectra": int(len(main)),
            "runs MS/MS": int(ms["Experiment"].nunique()), "runs main": int(len(main_runs)),
            "runs any": int(e["Experiment"].nunique()),
            "conditions MS/MS": ",".join(sorted(ms_conditions, key=CONDITION_ORDER.index)),
            "runs per condition": {c: int(sum(n > 0 for x, n in main_runs.items()
                                              if condition_of(x) == c)) for c in CONDITION_ORDER},
            "modified forms": "; ".join(f"{k} ({v})" for k, v in
                                        main["Modifications"].value_counts().items()),
            "isobaric twin": describer.isobaric_twin(r["Sequence"]),
            "msms": len(ms) > 0,
            "main": len(main) > 0,
            "pep": bool(len(main)) and main["PEP"].min() < PEP_MAX,
            "no_mock": not (ms_conditions & MOCK),
            "not_sc35m_only": not (ms_conditions and ms_conditions <= SC35M),
            "no_deamid_twin": twin_free,
            "runs": len(main_runs) >= MIN_RUNS,
            "no_deamid": bool(len(main)) and
            (~main["Modifications"].fillna("").str.contains("Deamidation")).any(),
        })
    df = pd.DataFrame(rows)
    df["no_isobaric"] = df["isobaric twin"].eq("")
    df["decoy"] = decoy
    return df


def apply_cascade(df):
    keep = pd.Series(True, index=df.index)
    df["lost at"] = ""
    for key, short, _ in FILTERS:
        lost = keep & ~df[key]
        df.loc[lost, "lost at"] = short
        keep &= df[key]
    df["level"] = 0
    for i, (key, _, _) in enumerate(FILTERS, 1):
        passed_all = df[[k for k, _, _ in FILTERS[:i]]].all(axis=1)
        df.loc[passed_all, "level"] = i
    failing = df[[k for k, _, _ in FILTERS]].apply(
        lambda row: [s for (k, s, _) in FILTERS if not row[k]], axis=1)
    df["fails"] = failing.map(", ".join)
    df["fails only"] = failing.map(lambda f: f[0] if len(f) == 1 else "")
    return df


# ---------------------------------------------------------------- side effects
def viral_groups(txt_dir, db):
    """Per PR8 protein: is it in a protein group, with how many peptides and MS/MS."""
    path = txt_dir / "proteinGroups.txt"
    if not path.exists():
        return None
    pg = pd.read_csv(path, sep="\t", low_memory=False,
                     usecols=lambda c: c in {"Protein IDs", "Peptides", "Razor + unique peptides",
                                             "MS/MS count", "Number of proteins", "Decoy",
                                             "Reverse"})
    pg = pg[pg[decoy_column(pg)] != "+"]
    pr8, headers = read_fasta(db.paths["PR8"], "pipe")
    out = {}
    for acc, header in headers.items():
        mnemonic = header.split("|")[2].split("_")[0]
        name = UNIPROT_NAMES.get(mnemonic, mnemonic)
        g = pg[pg["Protein IDs"].str.split(";").map(lambda ids: acc in ids)]
        out[name] = None if g.empty else g.iloc[0]
    return out


# ---------------------------------------------------------------- figures
def figure_cascade(levels):
    labels = [f"L{i}  {name}" for i, name in enumerate(levels["filter"])]
    fig, ax = plt.subplots(figsize=(9, 0.42 * len(levels) + 1.6))
    y = np.arange(len(levels))[::-1]
    left = np.zeros(len(levels))
    for origin in ORIGINS:
        vals = levels[origin].to_numpy(dtype=float)
        ax.barh(y, vals, left=left, height=0.62, color=ORIGIN_COLORS[origin], label=origin,
                edgecolor=SURFACE, linewidth=1.0)
        left += vals
    top = levels["targets"].max()
    for yi, (_, row) in zip(y, levels.iterrows()):
        ax.text(row["targets"] + top * 0.01, yi,
                f"{row['targets']}  (decoys {row['decoys']}, entries {row['entries']})",
                va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(y, labels, fontsize=9)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, top * 1.38)
    ax.set_xlabel("DPR-specific target peptides left")
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_title("Filters of increasing strictness: what is left after each step", pad=26)
    ax.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=8.5,
              handlelength=1.2, columnspacing=1.2)
    return save_fig(fig, "dpr_filter_cascade")


def figure_runs(df):
    """Main-search runs per condition for the peptides passing all filters and the near misses."""
    show = pd.concat([df[df["level"] == len(FILTERS)].assign(group="passes all filters"),
                      df[(df["level"] < len(FILTERS)) & df["fails only"].ne("")]
                      .assign(group=lambda d: "fails only: " + d["fails only"])])
    if show.empty:
        return None
    show = show.sort_values(["group", "runs main"], ascending=[True, False])
    show = pd.concat([show[show["group"] == "passes all filters"],
                      show[show["group"] != "passes all filters"]])
    M = np.array([[r["runs per condition"][c] for c in CONDITION_ORDER] for _, r in show.iterrows()])
    cmap = ListedColormap(RUN_RAMP)
    norm = BoundaryNorm(np.arange(-0.5, len(RUN_RAMP)), cmap.N)
    h = 0.24 * len(show) + 1.8
    fig, ax = plt.subplots(figsize=(10.5, h))
    ax.imshow(np.minimum(M, len(RUN_RAMP) - 1), cmap=cmap, norm=norm, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j]:
                ax.text(j, i, M[i, j], ha="center", va="center", fontsize=6.5,
                        color=SURFACE if M[i, j] >= 4 else INK)
    ylab = [f"{s[:30] + '…' if len(s) > 31 else s}  ·  {o}"
            for s, o in zip(show["Sequence"], show["origin"])]
    ax.set_yticks(range(len(show)), ylab, fontsize=7, family="monospace")
    ax.set_xticks(range(len(CONDITION_ORDER)), CONDITION_ORDER, fontsize=8)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    groups = show["group"].tolist()
    for i in range(1, len(groups)):
        if groups[i] != groups[i - 1]:
            ax.axhline(i - 0.5, color=INK2, linewidth=0.8)
    for i, g in enumerate(groups):
        if i == 0 or groups[i] != groups[i - 1]:
            ax.text(len(CONDITION_ORDER) - 0.4, i - 0.1, g, fontsize=7.5, color=INK2,
                    va="top", ha="left")
    ax.set_title("Runs with a main-search identification, per condition", pad=28)
    handles = [Patch(color=RUN_RAMP[k], label=str(k) if k < len(RUN_RAMP) - 1 else f"{k}+")
               for k in range(len(RUN_RAMP))]
    ax.legend(handles=handles, title="runs", ncol=len(handles), loc="upper left",
              bbox_to_anchor=(0, 0), borderaxespad=0.4, fontsize=7.5, title_fontsize=8)
    return save_fig(fig, "dpr_peptides_by_condition")


# ---------------------------------------------------------------- main
def main():
    start = time.time()
    global CONDITION_ORDER
    design = pd.read_csv(DESIGN, sep="\t")
    CONDITION_ORDER = list(dict.fromkeys(condition_of(e) for e in design["Experiment"]))

    print("reading the FASTA files ...")
    db = Databases()
    print(f"DPR entries: {len(db.dpr)}; canonical sequences: {db.n_canonical}")

    pep = pd.read_csv(TXT_DIR / "peptides.txt", sep="\t", low_memory=False,
                      usecols=lambda c: c in {"id", "Sequence", "Decoy", "Reverse", "PEP",
                                              "MS/MS Count"})
    dcol = decoy_column(pep)
    targets = pep[pep[dcol] != "+"]
    decoys = pep[pep[dcol] == "+"]

    # DPR-specific target peptides
    t_keys = {L(s) for s in targets["Sequence"]}
    t_map = tryptic_map(db.dpr, t_keys)
    in_dpr = targets[targets["Sequence"].map(lambda s: L(s) in t_map)]
    canonical_hit = in_dpr["Sequence"].map(lambda s: L(s) in db.canon_text)
    T = in_dpr[~canonical_hit].copy()
    shared_with_canonical = int(canonical_hit.sum())
    # peptides that only Mount Sinai VIRAL entries explain (sequence variants, not DPRs)
    variant = targets[targets["Sequence"].map(
        lambda s: L(s) in db.viral_own_text and L(s) not in db.core_text)]

    # DPR-specific decoy peptides
    d_keys = {L(s) for s in decoys["Sequence"]}
    d_map = tryptic_map({k: revert(v) for k, v in db.dpr.items()}, d_keys)
    D = decoys[decoys["Sequence"].map(lambda s: L(s) in d_map and L(s) not in db.canon_rev_text)].copy()
    print(f"DPR-specific peptides: {len(T)} targets, {len(D)} decoys "
          f"({shared_with_canonical} more peptides of DPR entries also occur in canonical sequences)")

    ev = evidence_for(list(T["id"]) + list(D["id"]))
    canonical_peptides = [s for s in targets["Sequence"] if L(s) not in t_map]
    describer = Describer(db, canonical_peptides)
    T = T.merge(pd.DataFrame([dict(zip(("origin", "origin detail"), describer.origin(s)), id=i)
                              for i, s in zip(T["id"], T["Sequence"])]), on="id")
    T["DPR entries"] = T["Sequence"].map(lambda s: ";".join(sorted(t_map[L(s)])))
    T["n DPR entries"] = T["Sequence"].map(lambda s: len(t_map[L(s)]))
    TF = apply_cascade(peptide_filters(T, ev, db, describer, decoy=False).merge(
        T[["id", "origin", "origin detail", "DPR entries", "n DPR entries", "PEP", "MS/MS Count"]]
        .rename(columns={"PEP": "peptides.txt PEP", "MS/MS Count": "peptides.txt MS/MS"}),
        on="id"))
    DF = apply_cascade(peptide_filters(D, ev, db, describer, decoy=True))

    # cascade table
    rows = []
    for i in range(len(FILTERS) + 1):
        kt, kd = TF[TF["level"] >= i], DF[DF["level"] >= i]
        entries = defaultdict(int)
        for ents in kt["DPR entries"]:
            for e in ents.split(";"):
                entries[e] += 1
        rows.append({"level": f"L{i}", "filter": "MaxQuant 1% FDR" if i == 0 else FILTERS[i - 1][1],
                     "targets": len(kt), "lost": 0 if i == 0 else int((TF["level"] == i - 1).sum()),
                     "decoys": len(kd),
                     "decoys/targets": f"{len(kd) / len(kt):.3f}" if len(kt) else "-",
                     "entries": len(entries), "entries >= 2": sum(v >= 2 for v in entries.values()),
                     **{o: int((kt["origin"] == o).sum()) for o in ORIGINS}})
    levels = pd.DataFrame(rows)

    # each filter on its own
    alone = pd.DataFrame([{"filter": short, "rule": rule,
                           "targets failing": int((~TF[key]).sum()),
                           "fail only this one": int((TF["fails only"] == short).sum()),
                           "lost here in the cascade": int((TF["lost at"] == short).sum()),
                           "decoys failing": int((~DF[key]).sum())}
                          for key, short, rule in FILTERS])

    # entries at the strictest level: smallest set explaining the peptides (greedy)
    final = TF[TF["level"] == len(FILTERS)]
    left = {r["Sequence"]: set(r["DPR entries"].split(";")) for _, r in final.iterrows()}
    cover = []
    while left:
        count = defaultdict(list)
        for s, ents in left.items():
            for e in ents:
                count[e].append(s)
        best = max(count, key=lambda e: (len(count[e]), -len(db.dpr[e]), e))
        cover.append((best, sorted(count[best])))
        for s in count[best]:
            left.pop(s)
    cover_df = pd.DataFrame([{"entry": e, **dpr_fields(db.own_headers[e]),
                              "peptides": len(p), "sequences": ";".join(p)} for e, p in cover])

    # tables
    out_cols = (["Sequence", "origin", "origin detail", "level", "lost at", "fails", "fails only",
                 "best PEP", "best score", "MS/MS spectra", "main-search spectra", "runs MS/MS",
                 "runs main", "runs any", "conditions MS/MS", "modified forms", "isobaric twin",
                 "n DPR entries", "DPR entries", "peptides.txt PEP", "peptides.txt MS/MS"]
                + [k for k, _, _ in FILTERS])
    TF = TF.sort_values(["level", "runs main", "MS/MS spectra"], ascending=False)
    t_peptides = save_table(TF[out_cols], "dpr_peptides")
    save_table(DF.drop(columns=["runs per condition"]), "dpr_decoys")
    t_levels = save_table(levels, "dpr_levels")
    save_table(alone, "dpr_filters")
    save_table(cover_df, "dpr_entries_strict")
    fig_cascade = figure_cascade(levels)
    fig_runs = figure_runs(TF)

    # canonical viral proteins with and without the additional file
    ours = viral_groups(TXT_DIR, db)
    base = viral_groups(BASELINE_TXT, db) if BASELINE_TXT != TXT_DIR else None

    # ------------------------------------------------------------ answers.md
    n0, nk = len(TF), len(final)
    d0, dk = len(DF), int((DF["level"] == len(FILTERS)).sum())
    lvl = levels.set_index("level")
    k_final = f"L{len(FILTERS)}"
    version = "?"
    params = TXT_DIR / "parameters.txt"
    if params.exists():
        for line in open(params):
            if line.startswith("Version\t"):
                version = line.split("\t")[1].strip()

    def pep_rows(df, limit=None):
        df = df if limit is None else df.head(limit)
        return [(f"`{r['Sequence']}`", r["origin"], r["origin detail"],
                 f"{r['best PEP']:.1e}" if pd.notna(r["best PEP"]) else "-",
                 r["main-search spectra"], r["runs main"], r["conditions MS/MS"] or "-",
                 r["n DPR entries"]) for _, r in df.iterrows()]

    pep_header = ["Peptide", "Origin", "Where", "Best PEP", "Main-search spectra", "Runs",
                  "Conditions (MS/MS)", "DPR entries"]
    parts = [
        "# Bergant_2023_adapted: evidence for DPR sequences",
        "",
        "Generated by `analysis_of_results/scripts/dpr_evidence.py`; do not edit by hand. "
        "Definitions are in the script's docstring; per-peptide tables go to `tables/` and figures "
        "to `figures/` (both not in git).",
        "",
        f"Search: MaxQuant {version}, 67 IAV_FPMS runs (PRIDE PXD040456), human Ensembl GRCh37, "
        f"PR8, SC35M and `parameter/additional_proteins.fasta` with {len(db.dpr):,} DPR entries "
        f"(`PR8_3frame_uniq_*`) plus {len(db.other_own)} VIRAL/CTRL entries.",
        "",
        "## Summary",
        "",
        f"- **{n0} peptides** at MaxQuant's 1% FDR occur only in DPR entries and in no canonical "
        f"sequence; {dk if dk else 'none'} of the {d0} decoys of this class "
        f"{'survives' if dk == 1 else 'survive'} all filters.",
        f"- In this class, decoys / targets is {lvl.loc['L0', 'decoys/targets']} before any filter, "
        f"against {len(decoys) / max(len(targets), 1):.4f} for all peptides of the search: MaxQuant's "
        "1% FDR holds for the search as a whole, not for this small class.",
        f"- **{nk} peptides pass all {len(FILTERS)} filters**: "
        + ", ".join(f"{int(lvl.loc[k_final, o])} {o}" for o in ORIGINS)
        + f". Decoys left: {dk}, so about {dk} of the {nk} may be "
        + ("a random match" if dk == 1 else "random matches")
        + f" (decoys / targets {lvl.loc[k_final, 'decoys/targets']}).",
        f"- Explaining them needs {len(cover_df)} DPR entries; "
        f"{int((cover_df['peptides'] >= 2).sum()) if len(cover_df) else 0} of them have two or more "
        "peptides at the strictest level.",
        f"- {len(variant)} peptides are explained only by the Mount Sinai VIRAL entries "
        "(sequence variants of PR8, not DPRs).",
        "",
        "## Filters and what each one removes",
        "",
        "Filters are applied in this order; each level keeps the peptides that pass all filters "
        "so far. Decoys go through the same filters.",
        "",
        md_table([(f"{i}", short, rule) for i, (_, short, rule) in enumerate(FILTERS, 1)],
                 ["#", "Filter", "Rule"]),
        "",
        md_table([(r["level"], r["filter"], r["targets"], r["lost"], r["decoys"],
                   r["decoys/targets"], r["entries"], r["entries >= 2"],
                   *[r[o] for o in ORIGINS]) for _, r in levels.iterrows()],
                 ["Level", "After", "Peptides", "Lost", "Decoys", "Decoys/peptides",
                  "DPR entries with one of them", "DPR entries with >= 2", *ORIGINS]),
        "",
        "Each filter on its own: peptides failing it, peptides failing only this filter (they would "
        "pass if it were dropped), and peptides it removes in the cascade above.",
        "",
        md_table([(r["filter"], r["targets failing"], r["fail only this one"],
                   r["lost here in the cascade"], r["decoys failing"]) for _, r in alone.iterrows()],
                 ["Filter", "Failing", "Failing only this one", "Lost here in the cascade",
                  "Decoys failing"]),
        "",
        f"![Filter cascade]({fig_cascade})",
        "",
        f"## Peptides that pass all filters ({nk})",
        "",
        md_table(pep_rows(final), pep_header) if nk else "None.",
        "",
    ]
    near = TF[(TF["level"] < len(FILTERS)) & TF["fails only"].ne("")]
    parts += [
        f"## Near misses: peptides that fail exactly one filter ({len(near)})",
        "",
        md_table([(f"`{r['Sequence']}`", r["fails only"], r["origin"], r["origin detail"],
                   f"{r['best PEP']:.1e}" if pd.notna(r["best PEP"]) else "-",
                   r["runs main"], r["conditions MS/MS"] or "-",
                   r["isobaric twin"] or r["modified forms"] or "-")
                  for _, r in near.sort_values(["fails only", "runs main"],
                                               ascending=[True, False]).iterrows()],
                 ["Peptide", "Fails only", "Origin", "Where", "Best PEP", "Runs",
                  "Conditions (MS/MS)", "Twin or modified forms"]) if len(near) else "None.",
        "",
    ]
    if fig_runs:
        parts += [f"![Runs per condition]({fig_runs})", ""]
    parts += ["## Peptides lost at each filter", "",
              "The five peptides with the most main-search spectra among those each filter removes "
              f"in the cascade; all of them are in `{t_peptides}`.", ""]
    for _, short, _ in FILTERS:
        lost = TF[TF["lost at"] == short]
        if lost.empty:
            continue
        parts += [f"**{short}** ({len(lost)} lost): " + "; ".join(
            f"`{r['Sequence']}` ({r['origin']}, {r['main-search spectra']} spectra"
            + (f", twin {r['isobaric twin'].split(' (')[0]}" if short == "no isobaric twin" else "")
            + (f", MS/MS in {r['conditions MS/MS']}" if short in ("not in mock", "not SC35M-only")
               else "") + ")"
            for _, r in lost.sort_values("main-search spectra", ascending=False).head(5).iterrows()),
            ""]
    parts += [
        "## DPR entries at the strictest level",
        "",
        "The smallest set of DPR entries that explains the peptides passing all filters (a peptide "
        "often occurs in several overlapping entries; ties go to the shortest entry). `n_delvg`, "
        "`n_orfs`, `span_any` and `len` come from the FASTA header.",
        "",
        md_table([(r["entry"], r.get("len", ""), r.get("n_delvg", ""), r.get("n_orfs", ""),
                   r.get("span_any", ""), r["peptides"], ", ".join(f"`{s}`" for s in
                                                                  r["sequences"].split(";")))
                  for _, r in cover_df.iterrows()],
                 ["Entry", "len", "n_delvg", "n_orfs", "span_any", "Peptides", "Sequences"])
        if len(cover_df) else "None.",
        "",
        "## Side effect: canonical viral proteins",
        "",
        "Long DPR entries contain most of a canonical viral protein plus a few residues more, so "
        "MaxQuant's protein inference can hand a canonical protein's peptides to a DPR entry and "
        "drop the canonical protein (it needs at least one razor peptide or a peptide set that "
        "another protein group contains completely). Viral protein quantities from this search are "
        "therefore not comparable with a search without the additional file.",
        "",
    ]
    if ours is not None:
        def cell(g):
            if g is None:
                return "in no protein group"
            n = int(g["Number of proteins"])
            return (f"group of {n} protein{'s' if n > 1 else ''}: {int(g['Peptides'])} peptides, "
                    f"{int(g['MS/MS count'])} MS/MS")
        parts += ["Per PR8 protein, the protein group that contains it (peptides and MS/MS are those "
                  "of the whole group):", ""]
        hdr = ["Protein", "This search (67 runs)"]
        rows = [(name, cell(g)) for name, g in ours.items()]
        if base is not None:
            hdr.append("Bergant_2023 (42 runs, no additional file)")
            rows = [(name, c, cell(base.get(name))) for name, c in rows]
        parts += [md_table(rows, hdr), ""]
        if base is not None:
            parts += ["The two searches differ in the number of runs as well (67 vs 42), so counts "
                      "should be higher here, not lower.", ""]
    parts += [
        "## Limitations",
        "",
        f"- Decoy counts in this class are small ({d0} before filtering), so the decoy ratio per "
        "level is a rough estimate.",
        "- The isobaric-twin check compares unmodified masses with parts of viral proteins and with "
        "identified canonical peptides; spectra explained by a modified or non-tryptic human "
        "peptide are not excluded. The spectra themselves (msms.txt) have not been inspected.",
        "- Canonical means the databases of the search: a human sequence missing from Ensembl "
        "GRCh37 (a variant or a newer isoform) can look DPR-specific; this mostly concerns peptides "
        "without a viral anchor.",
        "- The origin classes come from comparisons with the PR8 protein sequences only; the "
        "position of each DPR's junction is not in the FASTA header.",
        "",
        "## Inputs and software",
        "",
        "Canonical FASTA files: " + ", ".join(rel(p) for p in db.paths.values())
        + f", {rel(CONTAMINANTS)}; checked against results/run_info.txt"
        + ("" if not db.notes else " (" + "; ".join(db.notes) + ")") + ".",
        "",
        f"Python {platform.python_version()}, pandas {pd.__version__}, numpy {np.__version__}, "
        f"matplotlib {matplotlib.__version__}; runtime {time.time() - start:.0f} s.",
        "",
    ]
    out = OUT_DIR / "answers.md"
    out.write_text("\n".join(parts))
    print(f"answers written: {rel(out)}")
    print(levels[["level", "filter", "targets", "lost", "decoys", "entries"]].to_string(index=False))


if __name__ == "__main__":
    main()
