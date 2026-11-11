"""Create parameter/design.txt and parameter/mqpar.xml of Bergant_2023_67runs.

The 67-run search must differ from the 42-run search of Bergant_2023 only in its raw files, so
both files are derived from ../Bergant_2023/parameter/ without running MaxQuant:
- design.txt: the 42 rows of Bergant_2023, then the 25 runs that the paper does not use, named
  by the rule of make_parameters_Bergant_2023.sh (condition from the description, replicate =
  last number of the raw file name): A149V_1..6, SC35M_1..6, SC35MdNS1_1..6, mock2_1..6, and
  RSIG_3rem for REM_1_3, a second measurement of the RSIG_3 sample.
- mqpar.xml: the 42-run file with 25 entries appended to each of the six per-file lists
  (filePaths, experiments, fractions, ptms, paramGroupIndices, referenceChannel). Every other
  line stays byte-identical, which the script checks.
The runner's --dryrun validates the result on the server.
Usage:  python3 make_parameters_Bergant_2023_67runs.py           first time
        FORCE=1 python3 make_parameters_Bergant_2023_67runs.py   overwrite existing files
"""
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
SOURCE = PROJECT.parent / "Bergant_2023" / "parameter"
DESCRIPTION = SOURCE / "deposited" / "Sample_descriptions_IAV_FPMS.txt"
EXPECTED = 67
PER_FILE = {  # list tag -> entry template
    "filePaths": "<string>/work/data/{name}.raw</string>",
    "experiments": "<string>{experiment}</string>",
    "fractions": "<short>32767</short>",
    "ptms": "<boolean>False</boolean>",
    "paramGroupIndices": "<int>0</int>",
    "referenceChannel": "<string></string>",
}
EXTRA_ORDER = ["A149V", "SC35M", "SC35MdNS1", "mock2", "RSIG"]


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def condition(desc):
    desc = desc.strip()
    m = re.fullmatch(r"IAV-PR8 \((\w+)\)", desc)
    if m:
        return m.group(1)
    if desc in ("SC35M", "SC35MdNS1"):
        return desc
    if desc == "mock 2":
        return "mock2"
    fail(f"unexpected description {desc!r}")


def design_rows():
    with open(SOURCE / "design.txt", newline="") as fh:
        base = list(csv.DictReader(fh, delimiter="\t"))
    used = {r["Name"] for r in base}
    names = {r["Experiment"] for r in base}
    extra = []
    with open(DESCRIPTION, encoding="utf-8-sig", newline="") as fh:
        for cols in csv.reader(fh, delimiter="\t"):
            if len(cols) < 3 or not cols[0].lower().endswith(".raw") or cols[1] != "IAV FP-MS":
                continue
            name = cols[0][:-4]
            if name in used:
                continue
            cond = condition(cols[2])
            rep = int(re.search(r"_(\d+)$", name).group(1))
            exp = f"{cond}_{rep}"
            if exp in names:
                exp = f"{cond}_{rep}rem"
            if exp in names:
                fail(f"experiment name {exp} not unique")
            names.add(exp)
            extra.append((EXTRA_ORDER.index(cond), rep, {"Name": name, "Fraction": "",
                                                          "Experiment": exp, "PTM": "False"}))
    rows = base + [r for _, _, r in sorted(extra, key=lambda x: (x[0], x[1], x[2]["Name"]))]
    if len(rows) != EXPECTED:
        fail(f"{len(rows)} runs, expected {EXPECTED}")
    return rows, len(base)


def extend_mqpar(text, new_rows):
    out = text
    for tag, template in PER_FILE.items():
        m = re.search(rf"\n([ \t]*)</{tag}>", out)
        if not m:
            fail(f"<{tag}> not found in the source mqpar.xml")
        entry = re.search(rf"<{tag}>\s*\n([ \t]*)<", out)
        indent = entry.group(1) if entry else m.group(1) + "   "
        lines = "".join(f"\n{indent}" + template.format(name=r["Name"], experiment=r["Experiment"])
                        for r in new_rows)
        out = out[:m.start()] + lines + out[m.start():]
    return out


def without_lists(text):
    for tag in PER_FILE:
        text = re.sub(rf"<{tag}>.*?</{tag}>", f"<{tag}/>", text, flags=re.S)
    return text


def main():
    out_design = PROJECT / "parameter" / "design.txt"
    out_mqpar = PROJECT / "parameter" / "mqpar.xml"
    if (out_design.exists() or out_mqpar.exists()) and os.environ.get("FORCE") != "1":
        fail("parameter/design.txt or mqpar.xml exists, rerun with FORCE=1 to overwrite")
    rows, n_base = design_rows()
    src = (SOURCE / "mqpar.xml").read_text()
    new = extend_mqpar(src, rows[n_base:])

    # checks: list lengths, names, file paths, everything else unchanged
    root = ET.fromstring(new.encode())
    for tag in PER_FILE:
        n = len(root.find(tag))
        if n != EXPECTED:
            fail(f"<{tag}> has {n} entries, expected {EXPECTED}")
    exps = [e.text for e in root.find("experiments")]
    paths = [e.text for e in root.find("filePaths")]
    if len(set(exps)) != EXPECTED or len(set(paths)) != EXPECTED:
        fail("experiment names or file paths in mqpar.xml are not unique")
    # the 42-run file lists its runs in MaxQuant's order, so compare the mapping, not the order
    design_map = {f"/work/data/{r['Name']}.raw": r["Experiment"] for r in rows}
    if dict(zip(paths, exps)) != design_map:
        fail("file paths and experiments in mqpar.xml do not match design.txt")
    if without_lists(new) != without_lists(src):
        fail("settings outside the per-file lists changed")
    manifest = PROJECT / "data_manifest.tsv"
    if manifest.exists():
        with open(manifest, newline="") as fh:
            local = {r["local_path"] for r in csv.DictReader(fh, delimiter="\t")}
        missing = [r["Name"] for r in rows if f"data/{r['Name']}.raw" not in local]
        if missing:
            fail(f"raw files not in data_manifest.tsv: {', '.join(missing[:5])}")

    out_design.parent.mkdir(exist_ok=True)
    with open(out_design, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["Name", "Fraction", "Experiment", "PTM"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    out_mqpar.write_text(new)
    extra = rows[n_base:]
    print(f"written: parameter/design.txt ({len(rows)} runs: {n_base} of Bergant_2023 + "
          f"{len(extra)} more)")
    print("written: parameter/mqpar.xml (per-file lists extended, all other settings unchanged)")
    print("added experiments: " + ", ".join(r["Experiment"] for r in extra))


if __name__ == "__main__":
    main()
