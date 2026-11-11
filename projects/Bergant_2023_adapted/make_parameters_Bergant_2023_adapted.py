r"""Create parameter/design.txt and parameter/mqpar.xml of Bergant_2023_adapted, and check the own
FASTA file parameter/additional_proteins.fasta.

Bergant_2023_adapted is Bergant_2023_67runs (MaxQuant 2.9.0, all 67 IAV_FPMS runs, same settings)
with one more search database: the protein sequences in parameter/additional_proteins.fasta.
Both files are derived from ../Bergant_2023_67runs/parameter/ without running MaxQuant:
- design.txt: unchanged (67 runs)
- mqpar.xml: one more <FastaFileInfo> for /work/parameter/additional_proteins.fasta, with the
  identifier rule >([^\s]*) (the first word of the header is the protein ID) and otherwise the
  parse rules of the other FASTA files; every other line stays byte-identical (checked).
The runner's --dryrun validates mqpar.xml on the server.

Structure of parameter/additional_proteins.fasta:
    >PR8_3frame_uniq_000001 len=752 n_orfs=1 n_delvg=1 span_any=1
    MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL
    VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV
    ...
- One entry per protein: a header line starting with '>', then the sequence on one or more
  lines (any line length; 60 is customary). No comment lines inside an entry.
- Header: '>ID description'. The identifier rule >([^\s]*) of this file takes the first word
  (everything up to the first space) as the protein ID, PR8_3frame_uniq_000001 above; the rest
  is free text. UniProt-style headers also work, with the whole first word as the ID
  (>sp|P42212|GFP_AEQVI ... gives sp|P42212|GFP_AEQVI). 'OS=', 'OX=' (taxonomy ID) and 'GN='
  (gene name) are read when present.
- IDs: unique in the file and across all FASTA files of the search (human Ensembl IDs are
  ENSP..., the influenza proteins UniProt accessions such as P03496); no CON__ or REV__ prefix,
  no ';' (MaxQuant separates IDs with it).
- Sequence: one-letter amino-acid codes only (ACDEFGHIKLMNPQRSTVWY, and B, J, O, U, X, Z),
  upper case, no '*', digits or gaps.
- MaxQuant adds decoys (reversed sequences) and contaminants itself; do not add them.
After editing the file, check it with --check; mqpar.xml need not be recreated, because only the
path is stored there. The runner records the file's checksum in results/run_info.txt.

Usage:  python3 make_parameters_Bergant_2023_adapted.py           first time
        FORCE=1 python3 make_parameters_Bergant_2023_adapted.py   overwrite existing files
        python3 make_parameters_Bergant_2023_adapted.py --check   only check the own FASTA file
"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
SOURCE = PROJECT.parent / "Bergant_2023_67runs" / "parameter"
OWN_FASTA = "parameter/additional_proteins.fasta"
OWN_PATH = "/work/" + OWN_FASTA
DATABASES = ["ensembl_human_GRCh37-75.fa", "Flu_H1N1_PR8.fasta", "Flu_H7N7_SC35M.fasta"]
HEADER = re.compile(r"^>[^|\s]*\|([^|\s]+)\|")     # rule of the other FASTA files
OWN_RULE = r">([^\s]*)"                               # rule of the own file: first word
AMINO = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBJOUXZ]+$")


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def read_fasta(path):
    """[(header, id, sequence)]; fails on lines outside entries."""
    entries, header, seq = [], None, []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        line = line.rstrip("\r")
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq)))
            header, seq = line, []
        elif not line.strip():
            continue
        elif header is None:
            fail(f"{OWN_FASTA}, line {n}: sequence before the first header")
        else:
            seq.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


def database_ids():
    """IDs of the other FASTA files, from the first data/ folder that holds all three."""
    for folder in (PROJECT / "data", PROJECT.parent / "Bergant_2023_67runs" / "data",
                   PROJECT.parent / "Bergant_2023" / "data"):
        if all((folder / db).exists() for db in DATABASES):
            ids = set()
            for db in DATABASES:
                for line in open(folder / db):
                    if line.startswith(">"):
                        m = HEADER.match(line)
                        if m:
                            ids.add(m.group(1))
            return ids, f"../{folder.parent.name}/data" if folder.parent != PROJECT else "data"
    return None, None


def check_fasta():
    path = PROJECT / OWN_FASTA
    if not path.exists():
        fail(f"{OWN_FASTA} missing")
    entries = read_fasta(path)
    if not entries:
        fail(f"{OWN_FASTA} has no entries")
    ids = []
    for header, seq in entries:
        pid = header[1:].split()[0] if header[1:].split() else ""
        if not pid:
            fail(f"header without an ID (first word): {header[:80]!r}")
        if pid.startswith(("CON__", "REV__")) or ";" in pid:
            fail(f"ID {pid} starts with CON__ or REV__ or contains ';'")
        if not seq:
            fail(f"{pid} has no sequence")
        if not AMINO.match(seq):
            bad = sorted(set(re.sub(r"[ACDEFGHIKLMNPQRSTVWYBJOUXZ]", "", seq)))
            fail(f"{pid}: sequence has characters other than amino-acid codes: {''.join(bad)!r}")
        ids.append(pid)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        fail(f"IDs occur more than once: {', '.join(dup)}")
    known, where = database_ids()
    if known is None:
        print("note: the other FASTA files are not downloaded yet; ID clashes not checked")
    else:
        clash = sorted(set(ids) & known)
        if clash:
            fail(f"IDs also in the other FASTA files ({where}): {', '.join(clash)}")
    lengths = [len(s) for _, s in entries]
    print(f"{OWN_FASTA}: {len(entries)} entries, {min(lengths)}–{max(lengths)} amino acids, "
          f"IDs {', '.join(ids[:5])}{' ...' if len(ids) > 5 else ''}"
          + ("" if known is None else f"; no clash with {len(known)} IDs of the other databases"))


def extend_mqpar(text):
    blocks = list(re.finditer(r"([ \t]*)<FastaFileInfo>.*?</FastaFileInfo>\n", text, flags=re.S))
    if not blocks:
        fail("no <FastaFileInfo> in the source mqpar.xml")
    if OWN_PATH in text:
        fail(f"{OWN_PATH} is already listed in the source mqpar.xml")
    template = blocks[0].group(0)
    new_block = re.sub(r"<fastaFilePath>.*?</fastaFilePath>",
                       f"<fastaFilePath>{OWN_PATH}</fastaFilePath>", template, count=1)
    new_block = re.sub(r"<identifierParseRule>.*?</identifierParseRule>",
                       lambda m: f"<identifierParseRule>{OWN_RULE}</identifierParseRule>",
                       new_block, count=1)
    end = blocks[-1].end()
    return text[:end] + new_block + text[end:], len(blocks), template, new_block


def without_fasta_files(text):
    return re.sub(r"<fastaFiles>.*?</fastaFiles>", "<fastaFiles/>", text, flags=re.S)


def main():
    check_fasta()
    if "--check" in sys.argv[1:]:
        return
    out_design = PROJECT / "parameter" / "design.txt"
    out_mqpar = PROJECT / "parameter" / "mqpar.xml"
    if (out_design.exists() or out_mqpar.exists()) and os.environ.get("FORCE") != "1":
        fail("parameter/design.txt or mqpar.xml exists, rerun with FORCE=1 to overwrite")
    src = (SOURCE / "mqpar.xml").read_text()
    new, n_before, template, new_block = extend_mqpar(src)

    root = ET.fromstring(new.encode())
    paths = [e.text for e in root.iter("fastaFilePath")]
    if len(paths) != n_before + 1 or paths.count(OWN_PATH) != 1:
        fail("the extended mqpar.xml does not list the own FASTA file exactly once")
    if without_fasta_files(new) != without_fasta_files(src):
        fail("settings outside <fastaFiles> changed")
    strip = re.compile(r"<(fastaFilePath|identifierParseRule)>.*?</\1>")
    if strip.sub("", new_block) != strip.sub("", template):
        fail("the new <FastaFileInfo> differs from the others in more than path and ID rule")
    if f"<identifierParseRule>{OWN_RULE}</identifierParseRule>" not in new_block:
        fail("the identifier rule of the own FASTA file was not set")

    out_design.write_text((SOURCE / "design.txt").read_text())
    out_mqpar.write_text(new)
    n_runs = len(root.find("filePaths"))
    print(f"written: parameter/design.txt ({n_runs} runs, as Bergant_2023_67runs)")
    print(f"written: parameter/mqpar.xml (FASTA files: {', '.join(Path(p).name for p in paths)}; "
          "all other settings unchanged)")


if __name__ == "__main__":
    main()
