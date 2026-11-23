#!/usr/bin/env bash
# Download the input data of Laybourn_2024_adapted into data/ and parameter/deposited/.
#   Laybourn_2024 with one more search database, parameter/additional_proteins.fasta (in git:
#   proteins from deletion variants of the swH1N1 genome; make_parameters_Laybourn_2024_adapted.py).
#   Inputs and their sources as in ../Laybourn_2024
#   (get_data_Laybourn_2024.sh, PRIDE PXD052229): the 44 raw files of the paper, the three FASTA
#   files of the authors' search, the GEO GSE268254 series matrix and Supplementary Table S3.
#   Raw files already present with the size of the manifest in ../Laybourn_2024/data/ are
#   hard-linked (same disk, no extra space) or copied instead of downloaded, then checked.
# Usage:  ./get_data_Laybourn_2024_adapted.sh   download (or copy) what is missing
#   with VERIFY=1          also recheck checksums of present files
#   with SKIP_RAW=1        everything except the raw files
#   with MANIFEST_ONLY=1   create data_manifest.tsv and stop
# Every file is checked against data_manifest.tsv: raw files by size (PRIDE gives no checksums
# for PXD052229; the runner records their SHA-256 in results/run_info.txt), the small files by
# size and SHA-256. The manifest is created on the first run and kept in git.
# Details: maxquant_intro.md, section 6
set -euo pipefail

ACCESSION="PXD052229"
PRIDE_API="https://www.ebi.ac.uk/pride/ws/archive/v3/projects/$ACCESSION/files/all"
PRIDE_HTTPS="https://ftp.pride.ebi.ac.uk"
GEO_MATRIX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268254/matrix/GSE268254_series_matrix.txt.gz"
TABLE_S3="https://www.frontiersin.org/api/v4/articles/1432743/file/Table3.xlsx/1432743_supplementary-materials_tables_3_xlsx/1"
# run numbers n of the paper's raw files ..._MS<n>.raw (see above)
RUNS="$(seq 1 14; seq 19 26; seq 43 56; seq 61 68)"
EXPECTED_RAW=44
PARALLEL=4
MARGIN_GB=20        # free space kept after the download
MQ_FACTOR=2.2       # MaxQuant working files + combined/, per byte of raw data (tutorial_galaxy_2)

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
MANIFEST="data_manifest.tsv"
# raw files of the baseline project are linked or copied from there (':'-separated folders)
SIBLING_DATA="$(for d in ../Laybourn_2024/data; do
  (cd "$PROJECT_DIR/$d" 2> /dev/null && pwd) || true; done | paste -sd: -)"
export SIBLING_DATA
fail() { echo "ERROR: $*" >&2; exit 1; }

size_of() { wc -c < "$1" | tr -d ' '; }

# sibling_file DEST BYTES: the raw file of the same name and size in a related project's data/,
# empty if there is none
sibling_file() {
  local dest="$1" bytes="$2" dir dirs=()
  [[ "$dest" == data/*.raw && -n "${SIBLING_DATA:-}" ]] || return 0
  IFS=: read -r -a dirs <<< "$SIBLING_DATA"
  for dir in "${dirs[@]}"; do
    if [[ -f "$dir/${dest#data/}" && "$(size_of "$dir/${dest#data/}")" == "$bytes" ]]; then
      echo "$dir/${dest#data/}"
      return 0
    fi
  done
}

# verify FILE BYTES CHECKSUM, CHECKSUM sha256:HEX or "none" (size only); returns 1 for a wrong
# size, 2 for a wrong checksum
verify() {
  local f="$1" bytes="$2" checksum="$3"
  local algo="${checksum%%:*}" want="${checksum#*:}" size got
  size="$(size_of "$f")"
  [[ "$size" == "$bytes" ]] || { echo "ERROR: $f has $size bytes, expected $bytes" >&2; return 1; }
  [[ "$checksum" == "none" ]] && return 0
  got="$(shasum -a "${algo#sha}" "$f" | cut -d' ' -f1)"
  [[ "$got" == "$want" ]] || { echo "ERROR: $f fails the $algo check" >&2; return 2; }
}

# download_one KIND LOCAL_PATH BYTES CHECKSUM URL: resumable download into LOCAL_PATH.part,
# verified, then renamed. A verified file is never downloaded again. A raw file of the same
# name and size in ../Laybourn_2024/data/ is hard-linked (or copied) instead. With
# SKIP_RAW=1, raw files are skipped. curl's --retry skips dropped connections (curl: (56) ...),
# so each of up to 5 attempts is a new curl call that continues the .part file.
download_one() {
  local kind="$1" dest="$2" bytes="$3" checksum="$4" url="$5" part="$2.part"
  if [[ "$kind" == "raw" && "${SKIP_RAW:-0}" == "1" ]]; then
    return 0
  fi
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
  local sibling
  sibling="$(sibling_file "$dest" "$bytes")"
  if [[ -n "$sibling" ]]; then
    # a hard link where possible (same disk; raw files are never modified), otherwise a copy
    local from="../$(basename "$(dirname "$(dirname "$sibling")")")"
    if ln "$sibling" "$part" 2> /dev/null; then
      echo "linked from $from: $dest"
    else
      cp "$sibling" "$part" || { echo "ERROR: copying $from/$dest failed" >&2; return 1; }
      echo "copied from $from: $dest"
    fi
  fi
  # each attempt resumes the .part file; only attempts without progress count against the limit
  # of 5 in a row, so a large file survives many broken connections
  local stalls=0 before now
  while [[ ! -f "$part" || "$(size_of "$part")" -lt "$bytes" ]]; do
    before=0
    if [[ -f "$part" ]]; then before="$(size_of "$part")"; fi
    if curl -fsSL --retry 5 --retry-delay 10 -C - -o "$part" "$url"; then
      now=0
      if [[ -f "$part" ]]; then now="$(size_of "$part")"; fi
      if (( now >= bytes )); then break; fi
      # curl can end without an error before the end of the file: it takes a resumed request
      # that the server answers with 416 as complete. Without progress, start the file again.
      echo "curl ended at $now of $bytes bytes without an error: $dest" >&2
      if (( now <= before )); then rm -f "$part"; fi
    fi
    now=0
    if [[ -f "$part" ]]; then now="$(size_of "$part")"; fi
    if (( now > before )); then stalls=0; else stalls=$(( stalls + 1 )); fi
    (( stalls < 5 )) || { echo "ERROR: download failed, 5 attempts in a row without progress: $url" >&2; return 1; }
    echo "failed at $now of $bytes bytes, continuing in 60 s ($stalls of 5 attempts in a row" \
      "without progress): $dest" >&2
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
export -f size_of sibling_file verify download_one

# 1. prerequisites
for tool in curl shasum python3; do
  command -v "$tool" > /dev/null || fail "$tool not found"
done
mkdir -p data parameter/deposited

# 2. manifest, created once: raw files and FASTA from the PRIDE API; the small files (FASTA, GEO
# matrix, Table S3) are downloaded once to record their size and SHA-256
if [[ ! -f "$MANIFEST" ]]; then
  echo "creating $MANIFEST from the PRIDE API"
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' EXIT
  curl -fsS --retry 3 "$PRIDE_API" -o "$tmp/files.json" || fail "PRIDE API not reachable: $PRIDE_API"
  python3 -I - "$tmp/files.json" "$PRIDE_HTTPS" "$ACCESSION" "$EXPECTED_RAW" $RUNS \
    > "$tmp/rows.tsv" <<'PY'
import json, re, sys

path, https, accession, expected = sys.argv[1:5]
runs = {int(n) for n in sys.argv[5:]}
if len(runs) != int(expected):
    sys.exit(f"ERROR: {len(runs)} run numbers, expected {expected}")
fasta = {
    "uniprot-proteom_3AUP000008227_reviewed_3Ano.fasta":
        "data/uniprot-proteom_3AUP000008227_reviewed_3Ano.fasta",
    "A_Denmark_238_2020_aa.fasta": "data/A_Denmark_238_2020_aa.fasta",
    "A_swine_Denmark_10298-4_2017_aa.fasta": "data/A_swine_Denmark_10298-4_2017_aa.fasta",
}
found = {}
for f in json.load(open(path)):
    name = f["fileName"]
    m = re.search(r"_MS(\d+)\.raw$", name, flags=re.IGNORECASE)
    if m and int(m.group(1)) in runs:
        kind, local = "raw", "data/" + re.sub(r"\.raw$", "", name, flags=re.IGNORECASE) + ".raw"
        key = int(m.group(1))
    elif name in fasta:
        kind, local, key = "fasta", fasta[name], name
    else:
        continue
    if key in found:
        sys.exit(f"ERROR: two PRIDE files for {key}: {found[key][2]} and {local}")
    ftp = [l["value"] for l in f.get("publicFileLocations", []) if l["value"].startswith("ftp://")]
    if not ftp:
        sys.exit(f"ERROR: no FTP location for {name}")
    url = re.sub(r"^ftp://[^/]+", https, ftp[0])
    found[key] = (kind, local, f["fileSizeBytes"], url, f"PRIDE {accession}, {name}")
missing = sorted(runs - set(found)) + sorted(set(fasta) - set(found))
if missing:
    sys.exit("ERROR: not in PRIDE: " + ", ".join(str(m) for m in missing))
for kind, local, size, url, note in found.values():
    print("\t".join([kind, local, str(size), url, note]))
PY
  printf '%s\t%s\t%s\t%s\t%s\n' \
    deposited parameter/deposited/GSE268254_series_matrix.txt.gz 0 "$GEO_MATRIX" \
      "GEO GSE268254, series matrix (pig number and group of the RNA-seq samples)" \
    published data/deposited/Laybourn_2024_Table_S3.xlsx 0 "$TABLE_S3" \
      "Laybourn et al. 2024, Front Immunol 15:1432743, Supplementary Table S3" >> "$tmp/rows.tsv"
  {
    printf 'kind\tlocal_path\tbytes\tchecksum\tsource\tnote\n'
    while IFS=$'\t' read -r kind local bytes url note; do
      if [[ "$kind" == "raw" ]]; then
        printf '%s\t%s\t%s\tnone\t%s\t%s\n' "$kind" "$local" "$bytes" "$url" "$note"
        continue
      fi
      curl -fsSL --retry 3 -o "$tmp/file" "$url" || fail "download failed: $url"
      got="$(size_of "$tmp/file")"
      [[ "$bytes" == "0" || "$bytes" == "$got" ]] || fail "$url: $got bytes, PRIDE lists $bytes"
      printf '%s\t%s\t%s\tsha256:%s\t%s\t%s\n' "$kind" "$local" "$got" \
        "$(shasum -a 256 "$tmp/file" | cut -d' ' -f1)" "$url" "$note"
    done < "$tmp/rows.tsv"
  } | python3 -I -c '
import sys
lines = sys.stdin.read().splitlines()
order = {"deposited": 0, "fasta": 1, "published": 2, "raw": 3}
print(lines[0])
for line in sorted(lines[1:], key=lambda l: (order[l.split("\t")[0]], l.split("\t")[1])):
    print(line)' > "$MANIFEST.tmp"
  mv "$MANIFEST.tmp" "$MANIFEST"
fi
if [[ "${MANIFEST_ONLY:-0}" == "1" ]]; then
  echo "manifest: $MANIFEST ($(( $(wc -l < "$MANIFEST") - 1 )) files); MANIFEST_ONLY=1, stopping"
  exit 0
fi

# 3. disk space: the download plus a margin must fit; MaxQuant needs more later. Raw files that
# can be hard-linked from ../Laybourn_2024/data/ need no space.
link_ok=""   # whether a hard link from the related project works, tested once
missing_bytes=0
raw_bytes=0
while IFS=$'\t' read -r kind path bytes _; do
  [[ "$kind" == "kind" ]] && continue
  if [[ "$kind" == "raw" && "${SKIP_RAW:-0}" == "1" ]]; then continue; fi
  if [[ ! -f "$path" && "$kind" == "raw" ]]; then
    sibling="$(sibling_file "$path" "$bytes")"
    if [[ -n "$sibling" && -z "$link_ok" ]]; then
      link_ok=no
      if ln "$sibling" "data/.link_test.$$" 2> /dev/null; then link_ok=yes; fi
      rm -f "data/.link_test.$$"
    fi
    if [[ -n "$sibling" && "$link_ok" == "yes" ]]; then
      raw_bytes=$(( raw_bytes + bytes ))
      continue
    fi
  fi
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
    "unless run_Laybourn_2024_adapted.sh runs with --workdir on another disk."
fi

# 4. all files, $PARALLEL downloads at a time (small files first)
awk -F'\t' 'NR > 1 { print $1; print $2; print $3; print $4; print $5 }' "$MANIFEST" \
  | xargs -n 5 -P "$PARALLEL" bash -c 'download_one "$@"' _ \
  || fail "some downloads failed (messages above); rerun the script to resume"

# 5. FASTA checks (maxquant_intro.md, section 6)
for db in data/*.fasta parameter/*.fasta; do
  [[ -f "$db" ]] || continue
  entries="$(grep -c '^>' "$db")"
  decoys="$(grep -ciE '^>(rev_|decoy|##)' "$db" || true)"
  [[ "$decoys" == "0" ]] || fail "$db contains $decoys decoy entries"
  case "$db" in
    parameter/*)
      # own file, identifier rule >([^\s]*): the first word of the header is the ID
      bad="$(grep '^>' "$db" | grep -c -E '^>([[:space:]]|$)|^>[^[:space:]]*;' || true)"
      [[ "$bad" == "0" ]] || fail "$db has $bad headers without an ID as first word, or with ';' in it"
      bad="$(grep -v '^>' "$db" | grep -c -v -E '^[ACDEFGHIKLMNPQRSTVWYBJOUXZ]*$' || true)"
      [[ "$bad" == "0" ]] || fail "$db has $bad sequence lines with other characters than amino acids"
      ;;
    data/A_*)
      # the deposited viral files: ssi|ACCESSION|STRAIN_PROTEIN headers, sequences ending in '*';
      # make_parameters_Laybourn_2024.sh turns them into parameter/viral_proteins.fasta
      ;;
    *)
      # the default identifier rule >[^|]*\|(.*?)\| needs a header of the form >xx|ID|
      bad="$(grep '^>' "$db" | grep -c -v -E '^>[^|]*\|[^|]+\|' || true)"
      [[ "$bad" == "0" ]] || fail "$db has $bad headers without the form >xx|ID|"
      ;;
  esac
  echo "FASTA ok: $db ($entries entries)"
done

# 6. summary
n_raw="$(find data -maxdepth 1 -name '*.raw' | wc -l | tr -d ' ')"
if [[ "${SKIP_RAW:-0}" == "1" ]]; then
  echo "done (SKIP_RAW=1): $n_raw of $EXPECTED_RAW raw files in data/, FASTA in data/, manifest $MANIFEST"
  exit 0
fi
(( n_raw == EXPECTED_RAW )) || fail "$n_raw raw files in data/, expected $EXPECTED_RAW"
echo "done: $n_raw raw files in data/, FASTA in data/, manifest $MANIFEST"
