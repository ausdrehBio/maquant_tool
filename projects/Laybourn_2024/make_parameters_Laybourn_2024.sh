#!/usr/bin/env bash
# Create parameter/viral_proteins.fasta, parameter/design.txt and parameter/mqpar.xml of
# Laybourn_2024 for MaxQuant 2.9.0.
#
# The authors analysed their DIA data with Spectronaut 17.4 (no mqpar.xml deposited). This script
# creates a MaxDIA mqpar.xml with --create and transfers the authors' search settings (paper,
# section 2.3.2, and the PRIDE data processing protocol) where MaxQuant has an equivalent:
#   - variable modifications Oxidation (M), Acetyl (Protein N-term), Deamidation (N); fixed
#     Carbamidomethyl (C) (the --create default)
#   - 2 missed cleavages also for the DIA library (maxMissedCleavagesDia; --create default 1),
#     as in Spectronaut's default
#   - Trypsin/P only: the authors' "Trypsin/P and LysC" cleaves at the same sites, since LysC
#     cleaves after K, which Trypsin/P does as well
#   - label-free quantification with MaxLFQ (Spectronaut's MaxLFQ in the paper); no imputation
#     (nothing to set: MaxQuant does not impute)
# Settings without an equivalent stay at the MaxQuant defaults: the library is predicted from the
# FASTA files (on-the-fly prediction, library type OTF) instead of the authors' Spectronaut library
# from DIA search archives of pooled samples (that variant: ../Laybourn_2024_poollib), PSM and
# protein FDR 0.01 (the authors: precursor PEP < 0.1, protein q < 0.05 per run, protein PEP < 0.5),
# min. peptide length 7, max. 2 variable modifications per peptide in the DIA library, precursor
# charges 1-4, MaxQuant's contaminants added.
#
# Run type MaxDIA, not FAIMS-MaxDIA: the data were acquired with FAIMS at a single compensation
# voltage (-45 V). MaxQuant 2.9.0 does not quantify FAIMS-MaxDIA runs label-free (the dry run stops
# with "We do not support FAIMS DIA LFQ"), and it treats data with one CV as data without FAIMS
# (FAIMS is set only for more than one CV, changelog 2.8.1).
#
# FASTA: the three databases of the authors' search (PRIDE PXD052229, in data/, downloaded by
# get_data_Laybourn_2024.sh):
#   - uniprot-proteom_3AUP000008227_reviewed_3Ano.fasta: Sus scrofa, 48,349 UniProt entries
#     (>tr|ID|...), MaxQuant's default identifier rule, used as it is
#   - A_Denmark_238_2020_aa.fasta (huH1N1) and A_swine_Denmark_10298-4_2017_aa.fasta (swH1N1):
#     not usable as they are: with the default identifier rule both files give the IDs M1, M2,
#     NEP/NS2, NS1 and PAX, sequences end in '*', each file holds two translations of a whole
#     segment (M, NS) through 5 stop codons, and one header lacks the first field
#     (>OQ062647|A_Denmark_238_2020_PB2). Step 1 writes parameter/viral_proteins.fasta (in git):
#     the 22 proteins (11 per strain), '*' removed, the two segment translations per strain left
#     out, ID = the GenBank accession of the header where it has one (as in the paper's
#     Supplementary Table S3: MT666905 = swH1N1 NP), otherwise <strain>_<protein> (swH1N1_M1).
#     Identifier rule of this file: >([^\s]*), the first word of the header. The deposited files
#     stay unchanged in data/ (searched as they are in ../Laybourn_2024_depositedfasta).
#
# design.txt: the 44 runs of the paper, one sample per run, Experiment = <tissue>_<group>_p<pig>,
# e.g. LT_sw_p08. Sources of the assignment:
#   - group of each pig: GEO GSE268254 (RNA-seq of the same pigs, sample titles "Pig 7, swH1N1"):
#     control = pigs 1-6, swH1N1 = 7-14, huH1N1 = 19-26; Supplementary Table S4 lists the same pigs
#   - run and tissue: the authors' Spectronaut experiment 20231229_Trachea_DIA_hu+sw_final (PRIDE,
#     .sne, 14 GB) holds exactly the 44 runs ..._MS<n>.raw with n = 1-14, 19-26, 43-56, 61-68, and
#     the condition label of MS1 is "Upper_ctr". So n <= 26 is the upper trachea (UT) of pig n and
#     n >= 43 the lower trachea (LT) of pig n - 42 (the same gaps as the pig numbers). The run list
#     is in the processing log in the first 8 MB of the .sne (UTF-16), the label of MS1 at byte
#     262,069,946 (curl -r gives both without the whole file). The lower-trachea assignment rests
#     on the numbering; the viral proteins, found only in the swH1N1 pigs, check the group
#     assignment of both tissues after the search.
#
# The raw files need not be downloaded: --create and --dryrun do not open them, so both run in a
# temporary folder with empty placeholders named after design.txt. --create fails now and then
# with "Exception has been thrown by the target of an invocation" and a truncated file (seen on
# the Mac, x64 emulated); it is repeated up to 3 times.
#
# Usage:  ./make_parameters_Laybourn_2024.sh           first time
#         FORCE=1 ./make_parameters_Laybourn_2024.sh   overwrite existing files
# Writes: parameter/viral_proteins.fasta, parameter/design.txt, parameter/mqpar.xml,
#         parameter/mqpar_default.xml (--create output before the changes, not in git)
set -euo pipefail

IMAGE="localhost/maxquant:2.9.0"
NUM_THREADS=8
PIG_FASTA="uniprot-proteom_3AUP000008227_reviewed_3Ano.fasta"
HU_FASTA="A_Denmark_238_2020_aa.fasta"
SW_FASTA="A_swine_Denmark_10298-4_2017_aa.fasta"
VIRAL_FASTA="parameter/viral_proteins.fasta"
GEO_MATRIX="parameter/deposited/GSE268254_series_matrix.txt.gz"
EXPECTED_RUNS=44

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
fail() { echo "ERROR: $*" >&2; exit 1; }

# 1. prerequisites
command -v podman > /dev/null || fail "podman not found"
command -v python3 > /dev/null || fail "python3 not found"
podman image exists "$IMAGE" || fail "image $IMAGE missing, see maxquant_intro.md section 1"
for f in "data/$PIG_FASTA" "data/$HU_FASTA" "data/$SW_FASTA" "$GEO_MATRIX" data_manifest.tsv; do
  [[ -f "$f" ]] || fail "$f missing, run get_data_Laybourn_2024.sh (SKIP_RAW=1 is enough)"
done
if [[ -f parameter/mqpar.xml || -f parameter/design.txt || -f "$VIRAL_FASTA" ]]; then
  [[ "${FORCE:-0}" == "1" ]] || fail "parameter files exist, rerun with FORCE=1 to overwrite"
fi

# 2. parameter/viral_proteins.fasta from the two deposited viral files
python3 -I - "data/$HU_FASTA" huH1N1 "A/Denmark/238/2020(H1N1)" \
  "data/$SW_FASTA" swH1N1 "A/swine/Denmark/2017_10298_4_4p1/2017(H1N1)" \
  > "$VIRAL_FASTA.tmp" <<'PY'
import re, sys

AMINO = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
ACCESSION = re.compile(r"^[A-Z]{1,2}\d{5,6}$")        # GenBank nucleotide accession
EXPECTED = {"PB2", "PB1", "PA", "PAX", "HA", "NP", "NA", "M1", "M2", "NEP", "NS1"}
out, ids = [], set()
for path, label, strain in (sys.argv[1:4], sys.argv[4:7]):
    entries, header = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
            entries.append([header, ""])
        elif line:
            entries[-1][1] += line
    proteins = set()
    for header, seq in entries:
        fields = header.split()[0].split("|")
        # last field: <strain>_<protein>, e.g. A_Denmark_238_2020_NP
        protein = fields[-1].rsplit("_", 1)[-1]
        if "*" in seq.rstrip("*"):
            sys.stderr.write(f"left out: {path}: {header} (stop codons inside, segment translation)\n")
            continue
        seq = seq.rstrip("*")
        if not AMINO.match(seq):
            sys.exit(f"ERROR: {path}: {header}: characters other than amino acids")
        acc = [f for f in fields if ACCESSION.match(f)]
        pid = acc[0] if acc else f"{label}_{protein}"
        if pid in ids:
            sys.exit(f"ERROR: duplicate ID {pid}")
        ids.add(pid)
        proteins.add(protein)
        out.append(f">{pid} {label} {protein} OS=Influenza A virus ({strain}) GN={protein}")
        out += [seq[i:i + 60] for i in range(0, len(seq), 60)]
    if proteins != EXPECTED:
        sys.exit(f"ERROR: {path}: proteins {sorted(proteins)}, expected {sorted(EXPECTED)}")
print("\n".join(out))
PY
mv "$VIRAL_FASTA.tmp" "$VIRAL_FASTA"
echo "written: $VIRAL_FASTA ($(grep -c '^>' "$VIRAL_FASTA") proteins)"

# 3. design.txt: raw files from data_manifest.tsv, groups from the GEO series matrix
python3 -I - data_manifest.tsv "$GEO_MATRIX" "$EXPECTED_RUNS" > parameter/design.txt.tmp <<'PY'
import gzip, re, sys
from collections import Counter

manifest, matrix, expected = sys.argv[1], sys.argv[2], int(sys.argv[3])
groups = {"control": "ctr", "swH1N1": "sw", "huH1N1": "hu"}
pig_group = {}
for line in gzip.open(matrix, "rt"):
    if line.startswith("!Sample_title"):
        for title in line.rstrip("\n").split("\t")[1:]:
            m = re.fullmatch(r'"Pig (\d+), (\w+)"', title)
            if not m or m.group(2) not in groups:
                sys.exit(f"ERROR: unexpected GEO sample title {title}")
            pig_group[int(m.group(1))] = groups[m.group(2)]
if Counter(pig_group.values()) != Counter({"ctr": 6, "sw": 8, "hu": 8}):
    sys.exit(f"ERROR: GEO groups {Counter(pig_group.values())}, expected 6/8/8 pigs")
rows = []
for line in list(open(manifest))[1:]:
    kind, path = line.split("\t")[:2]
    if kind != "raw":
        continue
    name = path[len("data/"):-len(".raw")]
    n = int(re.search(r"_MS(\d+)$", name).group(1))
    tissue, pig = ("UT", n) if n <= 26 else ("LT", n - 42)
    if pig not in pig_group:
        sys.exit(f"ERROR: {name}: pig {pig} not in GEO GSE268254")
    group = pig_group[pig]
    rows.append(((tissue != "UT", list(groups.values()).index(group), pig), name,
                 f"{tissue}_{group}_p{pig:02d}"))
count = Counter(exp.rsplit("_", 1)[0] for _, _, exp in rows)
want = {f"{t}_{g}": c for t in ("UT", "LT") for g, c in (("ctr", 6), ("sw", 8), ("hu", 8))}
if len(rows) != expected or len({r[2] for r in rows}) != expected or dict(count) != want:
    sys.exit(f"ERROR: {len(rows)} runs, per tissue and group {dict(count)}, expected {want}")
print("Name\tFraction\tExperiment\tPTM")
for _, name, exp in sorted(rows):
    print(f"{name}\t\t{exp}\tFalse")
PY
mv parameter/design.txt.tmp parameter/design.txt
echo "written: parameter/design.txt ($EXPECTED_RUNS runs)"

# 4. --create in a temporary folder with placeholder raw files and copies of the FASTA files
tmp="$(mktemp -d "$PROJECT_DIR/parameter/.create.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/data" "$tmp/parameter"
cp parameter/design.txt "$VIRAL_FASTA" "$tmp/parameter/"
cp "data/$PIG_FASTA" "$tmp/data/"
tail -n +2 parameter/design.txt | cut -f1 | while IFS= read -r name; do
  touch "$tmp/data/$name.raw"
done
mq_tmp() { podman run --rm --pull never --platform linux/amd64 --security-opt label=disable \
             --ulimit nofile=-1:-1 -v "$tmp":/work "$IMAGE" "$@"; }
for attempt in 1 2 3; do
  rm -f "$tmp/parameter/mqpar.xml"
  if mq_tmp --create --newMqpar /work/parameter/mqpar.xml \
       --LCMSType DIA --instrumentType TO \
       --pathFasta "/work/data/$PIG_FASTA" "/work/$VIRAL_FASTA" \
       --pathRawFileFolder /work/data \
       --experimentalDesign /work/parameter/design.txt \
       --useLFQ --numThreads "$NUM_THREADS" > "$tmp/create.log" 2>&1 \
     && python3 -I -c 'import sys, xml.etree.ElementTree as ET; ET.parse(sys.argv[1])' \
          "$tmp/parameter/mqpar.xml" 2> /dev/null; then
    break
  fi
  (( attempt < 3 )) || { cat "$tmp/create.log" >&2; fail "mqpar.xml creation failed 3 times"; }
  echo "--create failed (attempt $attempt of 3), repeating: $(grep -m1 -i error "$tmp/create.log" || true)"
done
cp "$tmp/parameter/mqpar.xml" parameter/mqpar_default.xml

# 5. transfer the authors' settings; every change must apply exactly once
python3 -I - parameter/mqpar_default.xml "$tmp/parameter/mqpar.xml" "/work/$VIRAL_FASTA" <<'PY'
import re, sys

src, dst, viral = sys.argv[1:4]
xml = open(src, encoding="utf-8").read()

def replace(old, new):
    global xml
    n = xml.count(old)
    if n != 1:
        sys.exit(f"ERROR: {old!r} found {n} times, expected once")
    xml = xml.replace(old, new)

def expect(text):
    """A --create default this script relies on, unchanged."""
    if xml.count(text) != 1:
        sys.exit(f"ERROR: {text!r} found {xml.count(text)} times, expected once")

def set_list(tag, values):
    """Set the <string> entries of the list <tag>...</tag> (the tag name must match exactly)."""
    global xml
    pat = re.compile(rf"^( *)<{tag}>\n(.*?)^ *</{tag}>\n", re.M | re.S)
    hits = pat.findall(xml)
    if len(hits) != 1:
        sys.exit(f"ERROR: list <{tag}> found {len(hits)} times, expected once")
    indent = hits[0][0]
    body = "".join(f"{indent}   <string>{v}</string>\n" for v in values)
    xml = pat.sub(lambda m: f"{indent}<{tag}>\n{body}{indent}</{tag}>\n", xml)

expect("<lcmsRunType>MaxDIA</lcmsRunType>")
expect("<diaLibraryType>0</diaLibraryType>")          # OTF: library predicted from the FASTA files
expect("<maxMissedCleavages>2</maxMissedCleavages>")
expect("<lfqMode>1</lfqMode>")
# identifier rule of the own viral file: the first word of the header
replace(f"<fastaFilePath>{viral}</fastaFilePath>\n"
        r"         <identifierParseRule>>[^|]*\|(.*?)\|</identifierParseRule>",
        f"<fastaFilePath>{viral}</fastaFilePath>\n"
        r"         <identifierParseRule>>([^\s]*)</identifierParseRule>")
set_list("variableModifications", ["Oxidation (M)", "Acetyl (Protein N-term)", "Deamidation (N)"])
replace("<maxMissedCleavagesDia>1</maxMissedCleavagesDia>",
        "<maxMissedCleavagesDia>2</maxMissedCleavagesDia>")
open(dst, "w", encoding="utf-8").write(xml)
PY

# 6. checks: runs and experiments as in design.txt, FASTA paths, dry run
python3 -I - "$tmp/parameter/mqpar.xml" parameter/design.txt <<'PY'
import sys, xml.etree.ElementTree as ET

root = ET.parse(sys.argv[1]).getroot()
files = [s.text for s in root.find("filePaths")]
exps = [s.text or "" for s in root.find("experiments")]
design = {}
for line in list(open(sys.argv[2]))[1:]:
    name, _, exp, _ = line.rstrip("\n").split("\t")
    design[f"/work/data/{name}.raw"] = exp
got = dict(zip(files, exps))
if len(files) != len(exps) or got != design:
    sys.exit("ERROR: filePaths/experiments in mqpar.xml do not match design.txt")
print(f"checked: {len(files)} runs, experiments as in design.txt")
PY
mq_tmp /work/parameter/mqpar.xml --dryrun > "$tmp/dryrun.txt" 2>&1 \
  || { cat "$tmp/dryrun.txt" >&2; fail "dry run failed"; }
echo "dry run ok: $(awk -F'\t' '$1 ~ /^[0-9]+$/' "$tmp/dryrun.txt" | wc -l | tr -d ' ') jobs"

cp "$tmp/parameter/mqpar.xml" parameter/mqpar.xml
echo "written: parameter/mqpar.xml; changes against the --create defaults:"
diff parameter/mqpar_default.xml parameter/mqpar.xml || true
