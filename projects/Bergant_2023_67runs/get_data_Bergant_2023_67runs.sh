#!/usr/bin/env bash
# Download the input data of Bergant_2023_67runs into data/ and parameter/deposited/.
#   Control project for the comparison in projects/Bergant_2023/analysis_of_results/: the same
#   search as Bergant_2023 (MaxQuant 2.9.0, same settings), but over all 67 IAV_FPMS runs like
#   the authors' deposited search, so that version and run-set effects can be separated.
#   raw files: all 67 IAV_FPMS runs of Sample_descriptions_IAV_FPMS.txt (PRIDE PXD040456): the
#              paper's 42 plus the 25 marked "not used" (SC35M, SC35MdNS1, PR8 A149V, mock 2 and
#              two remeasured runs). Raw files already present in ../Bergant_2023/data/ with the
#              size of the manifest are copied instead of downloaded, then checked like downloads.
#   FASTA:     the databases of the authors' IAV_FPMS search, saved in data/:
#              ensembl_human_GRCh37-75.fa, Flu_H1N1_PR8.fasta, Flu_H7N7_SC35M.fasta
#   deposited: the authors' mqpar_IAV_FPMS.xml and Sample_descriptions_IAV_FPMS.txt in
#              parameter/deposited/
# Usage:  ./get_data_Bergant_2023_67runs.sh                 download (or copy) what is missing
#         VERIFY=1 ./get_data_Bergant_2023_67runs.sh        also recheck checksums of present files
#         MANIFEST_ONLY=1 ./get_data_Bergant_2023_67runs.sh create data_manifest.tsv and stop
# Every file is checked against data_manifest.tsv (size and SHA-1 from PRIDE). The manifest is
# created on the first run and kept in git, so later runs verify against the same record.
# Details: maxquant_intro.md, section 6
set -euo pipefail

ACCESSION="PXD040456"
PRIDE_API="https://www.ebi.ac.uk/pride/ws/archive/v3/projects/$ACCESSION/files/all"
PRIDE_HTTPS="https://ftp.pride.ebi.ac.uk"
SAMPLE_DESCRIPTION="Sample_descriptions_IAV_FPMS.txt"
DATASET="IAV FP-MS"
EXPECTED_RAW=67
PARALLEL=4
MARGIN_GB=20        # free space kept after the download
MQ_FACTOR=2.2       # MaxQuant working files + combined/, per byte of raw data (tutorial_galaxy_2)

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
MANIFEST="data_manifest.tsv"
# raw files of the 42-run project are copied from there when present
SIBLING_DATA="$(cd "$PROJECT_DIR/../Bergant_2023/data" 2> /dev/null && pwd || true)"
export SIBLING_DATA
fail() { echo "ERROR: $*" >&2; exit 1; }

size_of() { wc -c < "$1" | tr -d ' '; }

# verify FILE BYTES sha1:HEX; returns 1 for a wrong size, 2 for a wrong checksum
verify() {
  local f="$1" bytes="$2" checksum="$3"
  local algo="${checksum%%:*}" want="${checksum#*:}" size got
  size="$(size_of "$f")"
  [[ "$size" == "$bytes" ]] || { echo "ERROR: $f has $size bytes, expected $bytes" >&2; return 1; }
  got="$(shasum -a "${algo#sha}" "$f" | cut -d' ' -f1)"
  [[ "$got" == "$want" ]] || { echo "ERROR: $f fails the $algo check" >&2; return 2; }
}

# download_one LOCAL_PATH BYTES CHECKSUM URL: resumable download into LOCAL_PATH.part,
# verified, then renamed. A verified file is never downloaded again. A raw file of the same
# name and size in ../Bergant_2023/data/ is copied into the .part file instead. curl's --retry
# skips dropped connections (curl: (56) ...), so each of up to 5 attempts is a new curl call
# that continues the .part file.
download_one() {
  local dest="$1" bytes="$2" checksum="$3" url="$4" part="$1.part"
  mkdir -p "$(dirname "$dest")"
  if [[ -f "$dest" ]]; then
    if [[ "${VERIFY:-0}" == "1" ]]; then
      verify "$dest" "$bytes" "$checksum" || return 1
    else
      [[ "$(size_of "$dest")" == "$bytes" ]] \
        || { echo "ERROR: $dest has the wrong size; delete it to download it again" >&2; return 1; }
    fi
    echo "present: $dest"
    return 0
  fi
  if [[ -f "$part" && "$(size_of "$part")" -gt "$bytes" ]]; then rm -f "$part"; fi
  local sibling="${SIBLING_DATA:-}/${dest#data/}"
  if [[ -n "${SIBLING_DATA:-}" && "$dest" == data/*.raw && -f "$sibling" \
        && "$(size_of "$sibling")" == "$bytes" ]]; then
    cp "$sibling" "$part" || { echo "ERROR: copying ../Bergant_2023/$dest failed" >&2; return 1; }
    echo "copied from ../Bergant_2023: $dest"
  fi
  local tries=0
  while [[ ! -f "$part" || "$(size_of "$part")" -lt "$bytes" ]]; do
    tries=$(( tries + 1 ))
    curl -fsSL --retry 5 --retry-delay 10 -C - -o "$part" "$url" && break
    (( tries < 5 )) || { echo "ERROR: download failed after $tries attempts: $url" >&2; return 1; }
    echo "failed, continuing in 60 s (attempt $(( tries + 1 )) of 5): $dest" >&2
    sleep 60
  done
  local status=0
  verify "$part" "$bytes" "$checksum" || status=$?
  if (( status == 1 )); then
    rm -f "$part"       # wrong size: start again on the next run
    return 1
  elif (( status == 2 )); then
    mv "$part" "$dest.badchecksum"   # complete but different: kept for inspection
    echo "ERROR: kept as $dest.badchecksum" >&2
    return 1
  fi
  mv "$part" "$dest"
  echo "downloaded: $dest"
}
export -f size_of verify download_one

# 1. prerequisites
for tool in curl shasum python3; do
  command -v "$tool" > /dev/null || fail "$tool not found"
done
mkdir -p data parameter/deposited

# 2. manifest, created once from the PRIDE API and the authors' sample description
if [[ ! -f "$MANIFEST" ]]; then
  echo "creating $MANIFEST from the PRIDE API"
  tmp_json="$(mktemp)"
  tmp_sd="$(mktemp)"
  trap 'rm -f "$tmp_json" "$tmp_sd"' EXIT
  curl -fsS --retry 3 "$PRIDE_API" -o "$tmp_json" || fail "PRIDE API not reachable: $PRIDE_API"
  sd_url="$(python3 -I - "$tmp_json" "$SAMPLE_DESCRIPTION" "$PRIDE_HTTPS" <<'PY'
import json, re, sys
for f in json.load(open(sys.argv[1])):
    if f["fileName"] == sys.argv[2]:
        ftp = [l["value"] for l in f["publicFileLocations"] if l["value"].startswith("ftp://")]
        print(re.sub(r"^ftp://[^/]+", sys.argv[3], ftp[0]))
PY
)"
  [[ -n "$sd_url" ]] || fail "$SAMPLE_DESCRIPTION not found in $ACCESSION"
  curl -fsS --retry 3 "$sd_url" -o "$tmp_sd" || fail "download failed: $sd_url"
  python3 -I - "$tmp_json" "$tmp_sd" "$PRIDE_HTTPS" "$ACCESSION" "$DATASET" "$EXPECTED_RAW" \
    > "$MANIFEST.tmp" <<'PY'
import json, re, sys

path, sd_path, https, accession, dataset, expected = sys.argv[1:7]
expected = int(expected)
# all raw files of the dataset, used in the paper or not
raw = set()
for line in open(sd_path, encoding="utf-8-sig"):
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) < 3 or not cols[0].lower().endswith(".raw") or cols[1] != dataset:
        continue
    raw.add(cols[0])
if len(raw) != expected:
    sys.exit(f"ERROR: the sample description lists {len(raw)} raw files, expected {expected}")

# local path for every PRIDE file this project needs
wanted = {
    "mqpar_IAV_FPMS.xml": ("deposited", "parameter/deposited/mqpar_IAV_FPMS.xml"),
    "Sample_descriptions_IAV_FPMS.txt":
        ("deposited", "parameter/deposited/Sample_descriptions_IAV_FPMS.txt"),
    "ensembl_human_GRCh37-75.fa": ("fasta", "data/ensembl_human_GRCh37-75.fa"),
    "Flu_H1N1_PR8.fasta": ("fasta", "data/Flu_H1N1_PR8.fasta"),
    "Flu_H7N7_SC35M.fasta": ("fasta", "data/Flu_H7N7_SC35M.fasta"),
}
for name in raw:
    wanted[name] = ("raw", "data/" + re.sub(r"\.raw$", "", name, flags=re.IGNORECASE) + ".raw")

order = {"deposited": 0, "fasta": 1, "raw": 2}
rows = []
for f in json.load(open(path)):
    name = f["fileName"]
    if name not in wanted:
        continue
    kind, local = wanted.pop(name)
    ftp = [l["value"] for l in f.get("publicFileLocations", []) if l["value"].startswith("ftp://")]
    if not ftp:
        sys.exit(f"ERROR: no FTP location for {name}")
    sha1 = f.get("checksum") or ""
    if not re.fullmatch(r"[0-9a-f]{40}", sha1):
        sys.exit(f"ERROR: no SHA-1 for {name}")
    url = re.sub(r"^ftp://[^/]+", https, ftp[0])
    rows.append(((order[kind], local), kind, local, f["fileSizeBytes"], "sha1:" + sha1, url,
                 f"PRIDE {accession}, {name}"))
if wanted:
    sys.exit("ERROR: not in PRIDE: " + ", ".join(sorted(wanted)))
print("kind\tlocal_path\tbytes\tchecksum\tsource\tnote")
for row in sorted(rows):
    print("\t".join(str(c) for c in row[1:]))
PY
  mv "$MANIFEST.tmp" "$MANIFEST"
fi
if [[ "${MANIFEST_ONLY:-0}" == "1" ]]; then
  echo "manifest: $MANIFEST ($(( $(wc -l < "$MANIFEST") - 1 )) files); MANIFEST_ONLY=1, stopping"
  exit 0
fi

# 3. disk space: the download plus a margin must fit; MaxQuant needs more later
missing_bytes=0
raw_bytes=0
while IFS=$'\t' read -r kind path bytes _; do
  [[ "$kind" == "kind" ]] && continue
  if [[ ! -f "$path" ]]; then
    part_bytes=0
    [[ -f "$path.part" ]] && part_bytes="$(size_of "$path.part")"
    missing_bytes=$(( missing_bytes + bytes - part_bytes ))
  fi
  if [[ "$kind" == "raw" ]]; then
    raw_bytes=$(( raw_bytes + bytes ))
  fi
done < "$MANIFEST"
free_bytes="$(( $(df -k . | awk 'NR == 2 { print $4 }') * 1024 ))"
echo "to download or copy: $(( missing_bytes / 1000000000 )) GB, free: $(( free_bytes / 1000000000 )) GB"
(( free_bytes > missing_bytes + MARGIN_GB * 1000000000 )) \
  || fail "not enough disk space: need $(( missing_bytes / 1000000000 + MARGIN_GB )) GB (download plus $MARGIN_GB GB margin)"
mq_bytes="$(awk -v r="$raw_bytes" -v f="$MQ_FACTOR" 'BEGIN { printf "%.0f", r * f }')"
if (( free_bytes - missing_bytes < mq_bytes )); then
  echo "WARNING: after the download, $(( (free_bytes - missing_bytes) / 1000000000 )) GB will be free;" \
    "the MaxQuant run needs about $(( mq_bytes / 1000000000 )) GB more (working files and combined/)," \
    "unless run_Bergant_2023_67runs.sh runs with --workdir on another disk."
fi

# 4. all files, $PARALLEL downloads at a time (small files first)
awk -F'\t' 'NR > 1 { print $2; print $3; print $4; print $5 }' "$MANIFEST" \
  | xargs -n 4 -P "$PARALLEL" bash -c 'download_one "$@"' _ \
  || fail "some downloads failed (messages above); rerun the script to resume"

# 5. FASTA checks (maxquant_intro.md, section 6)
for db in data/*.fa data/*.fasta; do
  [[ -f "$db" ]] || continue
  entries="$(grep -c '^>' "$db")"
  decoys="$(grep -ciE '^>(rev_|decoy|##)' "$db" || true)"
  [[ "$decoys" == "0" ]] || fail "$db contains $decoys decoy entries"
  # the default identifier rule >[^|]*\|(.*?)\| needs a header of the form >xx|ID|
  bad="$(grep '^>' "$db" | grep -c -v -E '^>[^|]*\|[^|]+\|' || true)"
  [[ "$bad" == "0" ]] || fail "$db has $bad headers without the form >xx|ID|"
  echo "FASTA ok: $db ($entries entries)"
done

# 6. summary
n_raw="$(find data -maxdepth 1 -name '*.raw' | wc -l | tr -d ' ')"
(( n_raw == EXPECTED_RAW )) || fail "$n_raw raw files in data/, expected $EXPECTED_RAW"
echo "done: $n_raw raw files in data/, FASTA in data/, manifest $MANIFEST"
