"""Create parameter/design.txt, parameter/mqpar.xml and parameter/viral_proteins.fasta of
Laybourn_2024_poollib.

Laybourn_2024 searches the 44 sample runs with a spectral library that MaxQuant predicts from the
FASTA files (library type OTF). The authors used an empirical library instead: Spectronaut built
it from DIA searches ("DIA search archives") of pooled samples, then searched the samples against
it (93,131 peptides, 10,003 proteins; paper, section 2.3.2). This variant does the same in
MaxQuant: ../Laybourn_2024_pools searches the pooled samples with MaxDIA, and its evidence.txt
and msms.txt are the library here (library type "MaxQuant DDA library", which also reads
evidence.txt and msms.txt of a DIA search: "Library generation from DIA samples", changelog
2.3.0.0). get_data_Laybourn_2024_poollib.sh copies them to data/library/.

All three files are derived from ../Laybourn_2024/parameter/ without running MaxQuant:
- design.txt and viral_proteins.fasta: unchanged (44 runs; the cleaned viral proteins)
- mqpar.xml: diaLibraryType 0 (OTF) -> 2 (MaxQuant DDA library), diaEvidencePaths =
  /work/data/library/evidence.txt, diaMsmsPaths = /work/data/library/msms.txt; every other line
  stays byte-identical (checked). (MaxQuantCmd --create --pathLibrary writes both files into both
  lists; the GUI pairs them, one file per list.)
The runner checks that the library files exist, records their SHA-256 in results/run_info.txt
and copies them with --workdir; its --dryrun validates mqpar.xml on the server.

Usage:  python3 make_parameters_Laybourn_2024_poollib.py           first time
        FORCE=1 python3 make_parameters_Laybourn_2024_poollib.py   overwrite existing files
"""
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
SOURCE = PROJECT.parent / "Laybourn_2024" / "parameter"
EVIDENCE = "/work/data/library/evidence.txt"
MSMS = "/work/data/library/msms.txt"


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def replace(text, old, new):
    if text.count(old) != 1:
        fail(f"{old!r} found {text.count(old)} times, expected once")
    return text.replace(old, new)


def main():
    out = PROJECT / "parameter"
    targets = [out / "design.txt", out / "mqpar.xml", out / "viral_proteins.fasta"]
    if any(t.exists() for t in targets) and os.environ.get("FORCE") != "1":
        fail("parameter files exist, rerun with FORCE=1 to overwrite")
    for name in ("design.txt", "mqpar.xml", "viral_proteins.fasta"):
        if not (SOURCE / name).exists():
            fail(f"Laybourn_2024/parameter/{name} missing, run make_parameters_Laybourn_2024.sh first")
    text = (SOURCE / "mqpar.xml").read_text(encoding="utf-8")

    new = replace(text, "<diaLibraryType>0</diaLibraryType>", "<diaLibraryType>2</diaLibraryType>")
    new = replace(new, "         <diaEvidencePaths>\n         </diaEvidencePaths>\n",
                  f"         <diaEvidencePaths>\n            <string>{EVIDENCE}</string>\n"
                  "         </diaEvidencePaths>\n")
    new = replace(new, "         <diaMsmsPaths>\n         </diaMsmsPaths>\n",
                  f"         <diaMsmsPaths>\n            <string>{MSMS}</string>\n"
                  "         </diaMsmsPaths>\n")

    # check: exactly the three intended changes (one changed line, two added lines)
    old_lines, new_lines = text.splitlines(), new.splitlines()
    added = [l.strip() for l in new_lines if l not in old_lines]
    removed = [l.strip() for l in old_lines if l not in new_lines]
    if (len(new_lines) != len(old_lines) + 2 or removed != ["<diaLibraryType>0</diaLibraryType>"]
            or sorted(added) != sorted(["<diaLibraryType>2</diaLibraryType>",
                                        f"<string>{EVIDENCE}</string>", f"<string>{MSMS}</string>"])):
        fail(f"unexpected changes: removed {removed}, added {added}")
    group = ET.fromstring(new).find("parameterGroups/parameterGroup")
    if ([e.text for e in group.find("diaEvidencePaths")] != [EVIDENCE]
            or [e.text for e in group.find("diaMsmsPaths")] != [MSMS]
            or len(group.find("diaLabelIndsForLibraryMatch")) != 1):
        fail("library lists in mqpar.xml are not as intended")

    out.mkdir(exist_ok=True)
    for name in ("design.txt", "viral_proteins.fasta"):
        (out / name).write_bytes((SOURCE / name).read_bytes())
    (out / "mqpar.xml").write_text(new, encoding="utf-8")
    print("written: parameter/design.txt, parameter/viral_proteins.fasta (both unchanged from "
          "Laybourn_2024), parameter/mqpar.xml (library: data/library/evidence.txt, msms.txt)")


if __name__ == "__main__":
    main()
