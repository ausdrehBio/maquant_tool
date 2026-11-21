r"""Create parameter/design.txt and parameter/mqpar.xml of Laybourn_2024_depositedfasta.

Laybourn_2024 searches the viral proteins in parameter/viral_proteins.fasta, a cleaned version
of the authors' two viral FASTA files (make_parameters_Laybourn_2024.sh). This variant searches
the two files as deposited, to see what the cleaning changes:
    data/A_Denmark_238_2020_aa.fasta              (huH1N1, 13 entries)
    data/A_swine_Denmark_10298-4_2017_aa.fasta    (swH1N1, 13 entries)
With MaxQuant's default identifier rule >[^|]*\|(.*?)\| both files would give the IDs M1, M2,
NEP/NS2, NS1 and PAX, and one header (>OQ062647|A_Denmark_238_2020_PB2) no ID at all. Both files
therefore get the identifier rule >([^\s]*): the whole first word of the header is the ID
(ssi|MT666905|A_swine_Denmark_10298-4_2017_NP), unique across both files. Everything else stays
as deposited: sequences ending in '*', and the translations of the whole M and NS segments with
5 stop codons inside each.

Both files are derived from ../Laybourn_2024/parameter/ without running MaxQuant:
- design.txt: unchanged (44 runs)
- mqpar.xml: the <FastaFileInfo> of /work/parameter/viral_proteins.fasta is replaced by one for
  each deposited file, with the rule >([^\s]*) and otherwise the parse rules of the replaced
  entry; every other line stays byte-identical (checked).
The runner's --dryrun validates mqpar.xml on the server.

Usage:  python3 make_parameters_Laybourn_2024_depositedfasta.py           first time
        FORCE=1 python3 make_parameters_Laybourn_2024_depositedfasta.py   overwrite existing files
"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
SOURCE = PROJECT.parent / "Laybourn_2024" / "parameter"
OWN = "/work/parameter/viral_proteins.fasta"
DEPOSITED = ["/work/data/A_Denmark_238_2020_aa.fasta",
             "/work/data/A_swine_Denmark_10298-4_2017_aa.fasta"]
FIRST_WORD = r">([^\s]*)"


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def main():
    out_design, out_mqpar = PROJECT / "parameter" / "design.txt", PROJECT / "parameter" / "mqpar.xml"
    if (out_design.exists() or out_mqpar.exists()) and os.environ.get("FORCE") != "1":
        fail("parameter/design.txt or mqpar.xml exists, rerun with FORCE=1 to overwrite")
    for f in (SOURCE / "design.txt", SOURCE / "mqpar.xml"):
        if not f.exists():
            fail(f"{f.relative_to(PROJECT.parent)} missing, run make_parameters_Laybourn_2024.sh first")
    design = (SOURCE / "design.txt").read_bytes()
    text = (SOURCE / "mqpar.xml").read_text(encoding="utf-8")

    # the <FastaFileInfo> block of the own viral file, replaced by one block per deposited file
    pat = re.compile(r"^( *)<FastaFileInfo>\n((?:(?!</FastaFileInfo>).)*?<fastaFilePath>"
                     + re.escape(OWN) + r"</fastaFilePath>.*?)^ *</FastaFileInfo>\n", re.M | re.S)
    hits = list(pat.finditer(text))
    if len(hits) != 1:
        fail(f"<FastaFileInfo> of {OWN} found {len(hits)} times, expected once")
    block = hits[0].group(0)
    if f"<identifierParseRule>{FIRST_WORD}</identifierParseRule>" not in block:
        fail(f"the {OWN} entry does not have the identifier rule {FIRST_WORD}")
    new_blocks = "".join(block.replace(OWN, path) for path in DEPOSITED)
    new = text[:hits[0].start()] + new_blocks + text[hits[0].end():]

    # check: the only change is the replaced block
    old_lines, new_lines = text.splitlines(True), new.splitlines(True)
    n_block = block.count("\n")
    start = text[:hits[0].start()].count("\n")
    if (old_lines[:start] != new_lines[:start]
            or old_lines[start + n_block:] != new_lines[start + 2 * n_block:]):
        fail("mqpar.xml changed outside the FASTA entries")
    root = ET.fromstring(new)
    paths = [e.findtext("fastaFilePath") for e in root.find("fastaFiles")]
    rules = [e.findtext("identifierParseRule") for e in root.find("fastaFiles")]
    want = ["/work/data/uniprot-proteom_3AUP000008227_reviewed_3Ano.fasta"] + DEPOSITED
    if paths != want or rules[1:] != [FIRST_WORD, FIRST_WORD]:
        fail(f"unexpected FASTA entries: {list(zip(paths, rules))}")

    out_design.parent.mkdir(exist_ok=True)
    out_design.write_bytes(design)
    out_mqpar.write_text(new, encoding="utf-8")
    print("written: parameter/design.txt (unchanged from Laybourn_2024), parameter/mqpar.xml")
    for p, r in zip(paths, rules):
        print(f"  FASTA {p}  identifier rule {r}")


if __name__ == "__main__":
    main()
