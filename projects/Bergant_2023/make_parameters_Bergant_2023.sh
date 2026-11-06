#!/usr/bin/env bash
# Create parameter/design.txt and parameter/mqpar.xml of Bergant_2023 for MaxQuant 2.9.0.
#
# The authors' mqpar_IAV_FPMS.xml (parameter/deposited/) was written by MaxQuant 1.5.6.2 and
# does not run with 2.9.0. This script creates a new mqpar.xml with --create and transfers the
# authors' settings that differ from the 2.9.0 defaults:
#   - variable modifications Oxidation (M), Acetyl (Protein N-term), Deamidation (NQ)
#   - first search with Oxidation (M) and Acetyl (Protein N-term) only
#   - LFQ min. ratio count 2 (--create default 1)
#   - iBAQ on, with log fit
#   - match time window 0.7 min (--create default 0.4)
# Settings that are the same in both versions stay at the 2.9.0 defaults: Trypsin/P, 2 missed
# cleavages, Carbamidomethyl (C), 20/4.5 ppm first/main search, 20 ppm FTMS MS/MS, max. charge
# 7, min. peptide length 7, PSM and protein FDR 0.01, min. 1 peptide, LFQ min. 3/6 edges,
# protein quantification on unmodified, Oxidation (M) and Acetyl (Protein N-term) peptides.
# Match between runs and LFQ are on, as in the authors' run.
#
# design.txt: the 42 runs of parameter/deposited/Sample_descriptions_IAV_FPMS.txt without the
# comment "not used", one sample per run. Experiment name = condition_replicate, e.g. RSIG_4:
# condition from the description ("IAV-PR8 (RSIG)" -> RSIG, "mock 1" -> mock), replicate = the
# last number of the raw file name, also for remeasured runs (..._VBP1_REM_1_4 -> RSIG_4).
#
# FASTA: the authors' databases in data/ (downloaded by get_data_Bergant_2023.sh). Their headers
# have the form >sp|ID|..., so MaxQuant's default identifier rule applies.
#
# The raw files need not be downloaded: --create and --dryrun do not open them, so both run in
# a temporary folder with empty placeholders named after design.txt. The mqpar.xml then lists
# exactly the runs of the design.
#
# Usage:  ./make_parameters_Bergant_2023.sh           first time
#         FORCE=1 ./make_parameters_Bergant_2023.sh   overwrite existing design.txt and mqpar.xml
# Writes: parameter/design.txt, parameter/mqpar.xml, parameter/mqpar_default.xml (--create
#         output before the changes, not in git)
set -euo pipefail

IMAGE="localhost/maxquant:2.9.0"
NUM_THREADS=8
FASTA_FILES=(ensembl_human_GRCh37-75.fa Flu_H1N1_PR8.fasta Flu_H7N7_SC35M.fasta)
SAMPLE_DESCRIPTION="parameter/deposited/Sample_descriptions_IAV_FPMS.txt"
EXPECTED_RUNS=42

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
fail() { echo "ERROR: $*" >&2; exit 1; }

# 1. prerequisites
command -v podman > /dev/null || fail "podman not found"
command -v python3 > /dev/null || fail "python3 not found"
podman image exists "$IMAGE" || fail "image $IMAGE missing, see maxquant_intro.md section 1"
[[ -f "$SAMPLE_DESCRIPTION" ]] || fail "$SAMPLE_DESCRIPTION missing, run get_data_Bergant_2023.sh"
for db in "${FASTA_FILES[@]}"; do
  [[ -f "data/$db" ]] || fail "data/$db missing, run get_data_Bergant_2023.sh"
done
if [[ -f parameter/mqpar.xml || -f parameter/design.txt ]]; then
  [[ "${FORCE:-0}" == "1" ]] || fail "parameter/mqpar.xml or design.txt exists, rerun with FORCE=1 to overwrite"
fi

# 2. design.txt from the authors' sample description
python3 -I - "$SAMPLE_DESCRIPTION" "$EXPECTED_RUNS" > parameter/design.txt.tmp <<'PY'
import re, sys

path, expected = sys.argv[1], int(sys.argv[2])
order = ["mock", "RSIG", "RSIR", "RFMG", "RFMR", "AFMG", "AFMR"]
rows = []
for line in open(path, encoding="utf-8-sig"):
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) < 3 or not cols[0].lower().endswith(".raw") or cols[1] != "IAV FP-MS":
        continue
    if len(cols) > 3 and cols[3].strip() == "not used":
        continue
    name, desc = cols[0][:-4], cols[2].strip()
    m = re.fullmatch(r"IAV-PR8 \((\w+)\)", desc)
    cond = m.group(1) if m else ("mock" if desc == "mock 1" else None)
    if cond not in order:
        sys.exit(f"ERROR: unexpected description {desc!r} for {name}")
    rep = int(re.search(r"_(\d+)$", name).group(1))
    rows.append((order.index(cond), rep, name, f"{cond}_{rep}"))
experiments = [r[3] for r in rows]
if len(rows) != expected or len(set(experiments)) != len(rows):
    sys.exit(f"ERROR: {len(rows)} runs, {len(set(experiments))} distinct experiments, expected {expected}")
print("Name\tFraction\tExperiment\tPTM")
for _, _, name, exp in sorted(rows):
    print(f"{name}\t\t{exp}\tFalse")
PY
mv parameter/design.txt.tmp parameter/design.txt
echo "written: parameter/design.txt ($EXPECTED_RUNS runs)"

# 3. --create in a temporary folder with placeholder raw files and copies of the FASTA files
tmp="$(mktemp -d "$PROJECT_DIR/parameter/.create.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/data" "$tmp/parameter"
cp parameter/design.txt "$tmp/parameter/"
tail -n +2 parameter/design.txt | cut -f1 | while IFS= read -r name; do
  touch "$tmp/data/$name.raw"
done
fasta_args=()
for db in "${FASTA_FILES[@]}"; do
  cp "data/$db" "$tmp/data/"
  fasta_args+=("/work/data/$db")
done
mq_tmp() { podman run --rm --pull never --platform linux/amd64 --security-opt label=disable \
             -v "$tmp":/work "$IMAGE" "$@"; }
mq_tmp --create --newMqpar /work/parameter/mqpar.xml \
  --LCMSType ST --instrumentType TO \
  --pathFasta "${fasta_args[@]}" \
  --pathRawFileFolder /work/data \
  --experimentalDesign /work/parameter/design.txt \
  --useLFQ --useMBR --numThreads "$NUM_THREADS" > "$tmp/create.log" 2>&1 \
  || { cat "$tmp/create.log" >&2; fail "mqpar.xml creation failed"; }
cp "$tmp/parameter/mqpar.xml" parameter/mqpar_default.xml

# 4. transfer the authors' settings; every change must apply exactly once
python3 -I - parameter/mqpar_default.xml "$tmp/parameter/mqpar.xml" <<'PY'
import re, sys

src, dst = sys.argv[1], sys.argv[2]
xml = open(src, encoding="utf-8").read()

def replace(old, new):
    global xml
    n = xml.count(old)
    if n != 1:
        sys.exit(f"ERROR: {old!r} found {n} times, expected once")
    xml = xml.replace(old, new)

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

set_list("variableModifications", ["Oxidation (M)", "Acetyl (Protein N-term)", "Deamidation (NQ)"])
replace("<useVariableModificationsFirstSearch>False</useVariableModificationsFirstSearch>",
        "<useVariableModificationsFirstSearch>True</useVariableModificationsFirstSearch>")
set_list("variableModificationsFirstSearch", ["Oxidation (M)", "Acetyl (Protein N-term)"])
replace("<lfqMinRatioCount>1</lfqMinRatioCount>", "<lfqMinRatioCount>2</lfqMinRatioCount>")
replace("<ibaq>False</ibaq>", "<ibaq>True</ibaq>")
replace("<ibaqLogFit>False</ibaqLogFit>", "<ibaqLogFit>True</ibaqLogFit>")
replace("<matchingTimeWindow>0.4</matchingTimeWindow>", "<matchingTimeWindow>0.7</matchingTimeWindow>")
open(dst, "w", encoding="utf-8").write(xml)
PY

# 5. checks: runs and experiments as in design.txt, FASTA paths, dry run
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
echo "dry run ok: $(( $(wc -l < "$tmp/dryrun.txt") - 1 )) jobs"

cp "$tmp/parameter/mqpar.xml" parameter/mqpar.xml
echo "written: parameter/mqpar.xml; changes against the --create defaults:"
diff parameter/mqpar_default.xml parameter/mqpar.xml || true
