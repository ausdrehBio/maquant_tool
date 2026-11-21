"""Create parameter/design.txt, parameter/mqpar.xml and parameter/viral_proteins.fasta of
Laybourn_2024_pools.

The authors built their Spectronaut library from DIA searches of pooled samples (paper, section
2.3.2). This project is the MaxQuant counterpart: the search of Laybourn_2024 (MaxDIA, library
predicted from the FASTA files, the same settings) applied to the 45 pool runs instead of the 44
sample runs. Its evidence.txt and msms.txt are the library of ../Laybourn_2024_poollib.

The pools (PRIDE PXD052229; get_data_Laybourn_2024_pools.sh): control (CTR), swH1N1 (SWINE) and
huH1N1 (HUMAN) samples, each fractionated into 10 high-pH fractions (..._pool<GROUP>_frac<n>) and
measured unfractionated 5 times (..._pool<GROUP>_<n>).
design.txt: a fractionated sample per pool with the fraction number of the file name, e.g.
Experiment CTR_frac, Fraction 3; each unfractionated run a sample of its own, e.g. CTR_2.

All three files are derived from ../Laybourn_2024/parameter/ without running MaxQuant:
- viral_proteins.fasta: unchanged
- mqpar.xml: the six per-file lists (filePaths, experiments, fractions, ptms, paramGroupIndices,
  referenceChannel) list the 45 pool runs instead of the 44 sample runs; every other line stays
  byte-identical (checked).
The runner's --dryrun validates mqpar.xml on the server.

Usage:  python3 make_parameters_Laybourn_2024_pools.py           first time
        FORCE=1 python3 make_parameters_Laybourn_2024_pools.py   overwrite existing files
"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
SOURCE = PROJECT.parent / "Laybourn_2024" / "parameter"
POOLS = ["CTR", "SWINE", "HUMAN"]
EXPECTED = 45
NO_FRACTION = 32767                       # MaxQuant's value for "not fractionated"
PER_FILE = {  # list tag -> entry template
    "filePaths": "<string>/work/data/{name}.raw</string>",
    "experiments": "<string>{experiment}</string>",
    "fractions": "<short>{fraction}</short>",
    "ptms": "<boolean>False</boolean>",
    "paramGroupIndices": "<int>0</int>",
    "referenceChannel": "<string></string>",
}


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def design_rows():
    """[(name, fraction or None, experiment)] from the raw files of data_manifest.tsv."""
    manifest = PROJECT / "data_manifest.tsv"
    if not manifest.exists():
        fail("data_manifest.tsv missing, run get_data_Laybourn_2024_pools.sh (MANIFEST_ONLY=1 is enough)")
    rows = []
    for line in manifest.read_text().splitlines()[1:]:
        kind, path = line.split("\t")[:2]
        if kind != "raw":
            continue
        name = path[len("data/"):-len(".raw")]
        m = re.search(r"_pool([A-Z]+)_(frac)?(\d+)$", name)
        if not m or m.group(1) not in POOLS:
            fail(f"unexpected raw file {name}")
        pool, frac, n = m.group(1), m.group(2), int(m.group(3))
        if frac:
            rows.append(((POOLS.index(pool), 0, n), name, n, f"{pool}_frac"))
        else:
            rows.append(((POOLS.index(pool), 1, n), name, None, f"{pool}_{n}"))
    rows.sort()
    if len(rows) != EXPECTED:
        fail(f"{len(rows)} raw files in data_manifest.tsv, expected {EXPECTED}")
    return [r[1:] for r in rows]


def set_list(text, tag, entries):
    pat = re.compile(rf"^( *)<{tag}>\n(.*?)^( *)</{tag}>\n", re.M | re.S)
    hits = pat.findall(text)
    if len(hits) != 1:
        fail(f"list <{tag}> found {len(hits)} times, expected once")
    indent = hits[0][0]
    body = "".join(f"{indent}   {e}\n" for e in entries)
    return pat.sub(lambda m: f"{indent}<{tag}>\n{body}{indent}</{tag}>\n", text)


def without_lists(text):
    for tag in PER_FILE:
        text = re.sub(rf"^( *)<{tag}>\n.*?^ *</{tag}>\n", f"<{tag}/>\n", text, flags=re.M | re.S)
    return text


def main():
    out = PROJECT / "parameter"
    targets = [out / "design.txt", out / "mqpar.xml", out / "viral_proteins.fasta"]
    if any(t.exists() for t in targets) and os.environ.get("FORCE") != "1":
        fail("parameter files exist, rerun with FORCE=1 to overwrite")
    for name in ("mqpar.xml", "viral_proteins.fasta"):
        if not (SOURCE / name).exists():
            fail(f"Laybourn_2024/parameter/{name} missing, run make_parameters_Laybourn_2024.sh first")
    rows = design_rows()
    text = (SOURCE / "mqpar.xml").read_text(encoding="utf-8")

    new = text
    for tag, template in PER_FILE.items():
        entries = [template.format(name=name, experiment=exp,
                                   fraction=NO_FRACTION if frac is None else frac)
                   for name, frac, exp in rows]
        new = set_list(new, tag, entries)

    # checks: only the per-file lists changed, and they list the design
    if without_lists(new) != without_lists(text):
        fail("mqpar.xml changed outside the per-file lists")
    root = ET.fromstring(new)
    got = [(p.text, e.text, int(f.text)) for p, e, f in
           zip(root.find("filePaths"), root.find("experiments"), root.find("fractions"))]
    want = [(f"/work/data/{n}.raw", e, NO_FRACTION if f is None else f) for n, f, e in rows]
    if got != want or any(len(root.find(t)) != EXPECTED for t in PER_FILE):
        fail("per-file lists in mqpar.xml do not match the design")

    out.mkdir(exist_ok=True)
    with open(out / "design.txt", "w") as fh:
        fh.write("Name\tFraction\tExperiment\tPTM\n")
        for name, frac, exp in rows:
            fh.write(f"{name}\t{'' if frac is None else frac}\t{exp}\tFalse\n")
    (out / "viral_proteins.fasta").write_bytes((SOURCE / "viral_proteins.fasta").read_bytes())
    (out / "mqpar.xml").write_text(new, encoding="utf-8")
    n_frac = sum(f is not None for _, f, _ in rows)
    print(f"written: parameter/design.txt ({len(rows)} runs: {n_frac} fractions, "
          f"{len(rows) - n_frac} unfractionated), parameter/mqpar.xml, parameter/viral_proteins.fasta")


if __name__ == "__main__":
    main()
