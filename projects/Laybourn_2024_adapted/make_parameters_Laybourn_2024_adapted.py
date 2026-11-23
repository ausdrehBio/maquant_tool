r"""Create parameter/design.txt, parameter/mqpar.xml and parameter/viral_proteins.fasta of
Laybourn_2024_adapted, and check the own FASTA file parameter/additional_proteins.fasta.

Laybourn_2024_adapted is Laybourn_2024 (MaxQuant 2.9.0 MaxDIA, the paper's 44 runs, library
predicted from the FASTA files, the cleaned viral proteins) with one more search database: the
protein sequences in parameter/additional_proteins.fasta (provided as
swH1N1_3frame_unique_DPRs.fasta). MaxDIA predicts the library from all FASTA files of the search,
so the additional sequences enter the library without further settings.
All three files are derived from ../Laybourn_2024/parameter/ without running MaxQuant:
- design.txt and viral_proteins.fasta: unchanged (44 runs; the cleaned viral proteins)
- mqpar.xml: one more <FastaFileInfo> for /work/parameter/additional_proteins.fasta, a copy of
  the one of viral_proteins.fasta (identifier rule >([^\s]*): the first word of the header is
  the protein ID) with only the path changed; every other line stays byte-identical (checked).
The runner's --dryrun validates mqpar.xml on the server.

Content of parameter/additional_proteins.fasta: 2,946 entries swH1N1_3frame_uniq_<n>, protein
sequences derived from deletion variants of the swH1N1 genome (three-frame translations, unique
sequences), e.g.
    >swH1N1_3frame_uniq_000001 len=749 n_orfs=1 n_delvg=1 span_any=1
Three of them are identical to a canonical protein (M1, NS1, NEP of swH1N1), so their peptides
are never specific to a deletion variant. The file as provided also held 11 VIRAL_<protein>_swH1N1
entries (copies of the canonical swH1N1 proteins of viral_proteins.fasta) and 6 CTRL_<gene>_PIG
entries (Swiss-Prot pig GAPDH, ACTB, TUBA1B, VIM, HSP90AA1, RPLP0); they were removed so that
the search differs from Laybourn_2024 by the deletion-variant proteins only. The VIRAL entries
added no peptide, the CTRL entries none except 6 of vimentin, but both would have added IDs to
existing protein groups.
File rules (checked by --check): one entry per protein, a header '>ID description' and the
sequence on one or more lines; IDs unique in the file and across all FASTA files of the search,
no CON__ or REV__ prefix, no ';'; one-letter amino-acid codes only, no '*'. MaxQuant adds decoys
and contaminants itself. After editing the file, run --check; mqpar.xml need not be recreated,
because only the path is stored there. The runner records the file's checksum in
results/run_info.txt.

Usage:  python3 make_parameters_Laybourn_2024_adapted.py           first time
        FORCE=1 python3 make_parameters_Laybourn_2024_adapted.py   overwrite existing files
        python3 make_parameters_Laybourn_2024_adapted.py --check   only check the own FASTA file
"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
SOURCE = PROJECT.parent / "Laybourn_2024" / "parameter"
OWN_FASTA = "parameter/additional_proteins.fasta"
OWN_PATH = "/work/" + OWN_FASTA
VIRAL_PATH = "/work/parameter/viral_proteins.fasta"
PIG_FASTA = "uniprot-proteom_3AUP000008227_reviewed_3Ano.fasta"
PIG_HEADER = re.compile(r"^>[^|\s]*\|([^|\s]+)\|")   # default rule, used for the pig FASTA
OWN_RULE = r">([^\s]*)"                              # rule of viral_proteins.fasta and the own file
AMINO = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBJOUXZ]+$")
MIN_LENGTH = 7                                       # minPeptideLength of the search


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def read_fasta(path, label):
    """[(header, sequence)]; fails on lines outside entries."""
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
            fail(f"{label}, line {n}: sequence before the first header")
        else:
            seq.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


def other_databases():
    """IDs and sequences of the other FASTA files of the search: viral_proteins.fasta and the pig
    FASTA from the first data/ folder that holds it. Returns (ids, sequence -> ID, note)."""
    ids, seqs = set(), {}
    viral = SOURCE / "viral_proteins.fasta"
    if not viral.exists():
        fail("../Laybourn_2024/parameter/viral_proteins.fasta missing, run "
             "make_parameters_Laybourn_2024.sh first")
    for header, seq in read_fasta(viral, "viral_proteins.fasta"):
        pid = header[1:].split()[0]
        ids.add(pid)
        seqs.setdefault(seq, pid)
    for folder in (PROJECT / "data", PROJECT.parent / "Laybourn_2024" / "data"):
        pig = folder / PIG_FASTA
        if pig.exists():
            for header, seq in read_fasta(pig, PIG_FASTA):
                m = PIG_HEADER.match(header)
                if m:
                    ids.add(m.group(1))
                    seqs.setdefault(seq, m.group(1))
            where = "data" if folder.parent == PROJECT else "../Laybourn_2024/data"
            return ids, seqs, f"viral_proteins.fasta and {where}/{PIG_FASTA}"
    return ids, seqs, "viral_proteins.fasta only (the pig FASTA is not downloaded yet)"


def check_fasta():
    path = PROJECT / OWN_FASTA
    if not path.exists():
        fail(f"{OWN_FASTA} missing")
    entries = read_fasta(path, OWN_FASTA)
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
    seen, dup = set(), set()
    for i in ids:
        (dup if i in seen else seen).add(i)
    if dup:
        fail(f"IDs occur more than once: {', '.join(sorted(dup))}")
    known, known_seqs, where = other_databases()
    clash = sorted(set(ids) & known)
    if clash:
        fail(f"IDs also in the other FASTA files ({where}): {', '.join(clash)}")
    lengths = [len(s) for _, s in entries]
    prefixes = {}
    for i in ids:
        key = re.sub(r"_\d+$", "_<n>", i) if re.search(r"_\d+$", i) else i.split("_")[0] + "_*"
        prefixes[key] = prefixes.get(key, 0) + 1
    print(f"{OWN_FASTA}: {len(entries)} entries ("
          + ", ".join(f"{k}: {v}" for k, v in sorted(prefixes.items(), key=lambda kv: -kv[1]))
          + f"), {min(lengths)}-{max(lengths)} amino acids; no ID clash with {len(known)} IDs of "
          f"{where}")
    # notes, not errors: duplicates of other sequences and entries too short to give a peptide
    same = [(header[1:].split()[0], known_seqs[seq])
            for header, seq in entries if seq in known_seqs]
    if same:
        print(f"note: {len(same)} entries have the sequence of an entry in the other FASTA files "
              f"(their peptides are shared with it): "
              + ", ".join(f"{a} = {b}" for a, b in same[:6]) + (" ..." if len(same) > 6 else ""))
    seqs = {}
    for header, seq in entries:
        seqs.setdefault(seq, []).append(header[1:].split()[0])
    twins = [v for v in seqs.values() if len(v) > 1]
    if twins:
        print(f"note: {len(twins)} sequences occur more than once in the file: "
              + "; ".join(" = ".join(v) for v in twins[:6]))
    short = [header[1:].split()[0] for header, seq in entries if len(seq) < MIN_LENGTH]
    if short:
        print(f"note: {len(short)} entries are shorter than {MIN_LENGTH} amino acids and cannot "
              "give a peptide")


def extend_mqpar(text):
    blocks = list(re.finditer(r"([ \t]*)<FastaFileInfo>.*?</FastaFileInfo>\n", text, flags=re.S))
    if not blocks:
        fail("no <FastaFileInfo> in the source mqpar.xml")
    if OWN_PATH in text:
        fail(f"{OWN_PATH} is already listed in the source mqpar.xml")
    viral = [b for b in blocks if f"<fastaFilePath>{VIRAL_PATH}</fastaFilePath>" in b.group(0)]
    if len(viral) != 1:
        fail(f"<FastaFileInfo> of {VIRAL_PATH} found {len(viral)} times, expected once")
    template = viral[0].group(0)
    if f"<identifierParseRule>{OWN_RULE}</identifierParseRule>" not in template:
        fail(f"the {VIRAL_PATH} entry does not have the identifier rule {OWN_RULE}")
    new_block = template.replace(f"<fastaFilePath>{VIRAL_PATH}</fastaFilePath>",
                                 f"<fastaFilePath>{OWN_PATH}</fastaFilePath>")
    end = blocks[-1].end()
    return text[:end] + new_block + text[end:], len(blocks), template, new_block


def without_fasta_files(text):
    return re.sub(r"<fastaFiles>.*?</fastaFiles>", "<fastaFiles/>", text, flags=re.S)


def main():
    check_fasta()
    if "--check" in sys.argv[1:]:
        return
    out = PROJECT / "parameter"
    targets = [out / "design.txt", out / "mqpar.xml", out / "viral_proteins.fasta"]
    if any(t.exists() for t in targets) and os.environ.get("FORCE") != "1":
        fail("parameter files exist, rerun with FORCE=1 to overwrite")
    for name in ("design.txt", "mqpar.xml", "viral_proteins.fasta"):
        if not (SOURCE / name).exists():
            fail(f"Laybourn_2024/parameter/{name} missing, run make_parameters_Laybourn_2024.sh "
                 "first")
    src = (SOURCE / "mqpar.xml").read_text(encoding="utf-8")
    new, n_before, template, new_block = extend_mqpar(src)

    root = ET.fromstring(new.encode())
    paths = [e.text for e in root.iter("fastaFilePath")]
    if len(paths) != n_before + 1 or paths.count(OWN_PATH) != 1:
        fail("the extended mqpar.xml does not list the own FASTA file exactly once")
    if without_fasta_files(new) != without_fasta_files(src):
        fail("settings outside <fastaFiles> changed")
    if new_block.replace(OWN_PATH, VIRAL_PATH) != template:
        fail("the new <FastaFileInfo> differs from that of viral_proteins.fasta in more than "
             "the path")

    for name in ("design.txt", "viral_proteins.fasta"):
        (out / name).write_bytes((SOURCE / name).read_bytes())
    (out / "mqpar.xml").write_text(new, encoding="utf-8")
    n_runs = len(root.find("filePaths"))
    print(f"written: parameter/design.txt ({n_runs} runs) and parameter/viral_proteins.fasta, both "
          "unchanged from Laybourn_2024")
    print(f"written: parameter/mqpar.xml (FASTA files: {', '.join(Path(p).name for p in paths)}; "
          "all other settings unchanged)")


if __name__ == "__main__":
    main()
