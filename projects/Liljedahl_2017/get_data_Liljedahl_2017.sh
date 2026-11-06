#!/usr/bin/env bash
# Download the input data of Liljedahl_2017 into data/ and parameter/.
#   raw files: PRIDE PXD002988 (40 Thermo .RAW files, saved as .raw)
#   design:    PRIDE community SDRF, saved as parameter/PXD002988.sdrf.tsv
#   FASTA:     UniProt Swiss-Prot mouse reference proteome UP000000589, canonical sequences
# Usage:  ./get_data_Liljedahl_2017.sh            download what is missing
#         VERIFY=1 ./get_data_Liljedahl_2017.sh   also recheck checksums of files already present
# Every file is checked against data_manifest.tsv (size and SHA-1 from PRIDE, SHA-256 of the
# FASTA). The manifest is created on the first run and kept in git, so later runs verify
# against the same record. Details: maxquant_intro.md, section 6
set -euo pipefail

ACCESSION="PXD002988"
PRIDE_API="https://www.ebi.ac.uk/pride/ws/archive/v3/projects/$ACCESSION/files/all"
PRIDE_HTTPS="https://ftp.pride.ebi.ac.uk"
EXPECTED_RAW=40
# UniProt query "(proteome:UP000000589) AND (reviewed:true)", URL-encoded
UNIPROT_QUERY="%28proteome%3AUP000000589%29%20AND%20%28reviewed%3Atrue%29"
UNIPROT_STREAM="https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=$UNIPROT_QUERY"
UNIPROT_COUNT="https://rest.uniprot.org/uniprotkb/search?format=fasta&size=0&query=$UNIPROT_QUERY"
FASTA_PREFIX="UP000000589_mouse_swissprot"
PARALLEL=4

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
MANIFEST="data_manifest.tsv"
fail() { echo "ERROR: $*" >&2; exit 1; }

size_of() { wc -c < "$1" | tr -d ' '; }

# verify FILE BYTES ALGO:HEX  (ALGO sha1 or sha256, or "none" for size only); returns 1 for a
# wrong size, 2 for a wrong checksum
verify() {
  local f="$1" bytes="$2" checksum="$3"
  local algo="${checksum%%:*}" want="${checksum#*:}" size got
  size="$(size_of "$f")"
  [[ "$size" == "$bytes" ]] || { echo "ERROR: $f has $size bytes, expected $bytes" >&2; return 1; }
  [[ "$checksum" == "none" ]] && return 0
  got="$(shasum -a "${algo#sha}" "$f" | cut -d' ' -f1)"
  [[ "$got" == "$want" ]] || { echo "ERROR: $f fails the $algo check" >&2; return 2; }
}

# download_one LOCAL_PATH BYTES CHECKSUM URL: resumable download into LOCAL_PATH.part,
# verified, then renamed. A verified file is never downloaded again.
download_one() {
  local dest="$1" bytes="$2" checksum="$3" url="$4" part="$1.part"
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
  if [[ ! -f "$part" || "$(size_of "$part")" -lt "$bytes" ]]; then
    curl -fsSL --retry 5 --retry-delay 10 -C - -o "$part" "$url" \
      || { echo "ERROR: download failed: $url" >&2; return 1; }
  fi
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
mkdir -p data parameter

# 2. manifest of the PRIDE files (raw files and SDRF), created once
if [[ ! -f "$MANIFEST" ]]; then
  echo "creating $MANIFEST from the PRIDE API"
  tmp_json="$(mktemp)"
  trap 'rm -f "$tmp_json"' EXIT
  curl -fsS --retry 3 "$PRIDE_API" -o "$tmp_json" || fail "PRIDE API not reachable: $PRIDE_API"
  python3 -I - "$tmp_json" "$PRIDE_HTTPS" "$ACCESSION" "$EXPECTED_RAW" > "$MANIFEST.tmp" <<'PY'
import json, re, sys

path, https, accession, expected = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
files = json.load(open(path))
rows = []
for f in files:
    kind = f.get("fileCategory", {}).get("value")
    if kind not in ("RAW", "EXPERIMENTAL DESIGN"):
        continue
    name = f["fileName"]
    ftp = [l["value"] for l in f.get("publicFileLocations", []) if l["value"].startswith("ftp://")]
    if not ftp:
        sys.exit(f"ERROR: no FTP location for {name}")
    url = re.sub(r"^ftp://[^/]+", https, ftp[0])
    sha1 = f.get("checksum") or ""
    if re.fullmatch(r"[0-9a-f]{40}", sha1):
        checksum = "sha1:" + sha1
    elif kind == "RAW":
        sys.exit(f"ERROR: no SHA-1 for raw file {name}")
    else:
        checksum = "none"   # PRIDE gives no checksum for the SDRF; size is checked, git keeps it
    if kind == "RAW":
        local = "data/" + re.sub(r"\.raw$", "", name, flags=re.IGNORECASE) + ".raw"
        m = re.search(r"_NG([A-D])(\d+)_", name)
        key = (0, m.group(1), int(m.group(2))) if m else (0, name, 0)
        rows.append((key, "raw", local, f["fileSizeBytes"], checksum, url,
                     f"PRIDE {accession}, {name}"))
    else:
        rows.append(((1, name, 0), "sdrf", f"parameter/{accession}.sdrf.tsv", f["fileSizeBytes"],
                     checksum, url, f"PRIDE {accession}, {name}"))
n_raw = sum(r[1] == "raw" for r in rows)
if n_raw != expected or not any(r[1] == "sdrf" for r in rows):
    sys.exit(f"ERROR: PRIDE lists {n_raw} raw files (expected {expected}) or no SDRF")
print("kind\tlocal_path\tbytes\tchecksum\tsource\tnote")
for row in sorted(rows):
    print("\t".join(str(c) for c in row[1:]))
PY
  mv "$MANIFEST.tmp" "$MANIFEST"
fi

# 3. disk space: about 3x the raw data is needed (raw files, MaxQuant working files, results)
missing_bytes=0
while IFS=$'\t' read -r kind path bytes _; do
  if [[ "$kind" == "raw" && ! -f "$path" ]]; then
    missing_bytes=$(( missing_bytes + bytes ))
  fi
done < "$MANIFEST"
free_bytes="$(( $(df -k . | awk 'NR == 2 { print $4 }') * 1024 ))"
need_bytes="$(( missing_bytes * 3 ))"
echo "raw data to download: $(( missing_bytes / 1000000 )) MB, free: $(( free_bytes / 1000000000 )) GB"
(( free_bytes > need_bytes )) || fail "not enough disk space: need about $(( need_bytes / 1000000000 )) GB"

# 4. raw files and SDRF, $PARALLEL downloads at a time
awk -F'\t' '$1 == "raw" || $1 == "sdrf" { print $2; print $3; print $4; print $5 }' "$MANIFEST" \
  | xargs -n 4 -P "$PARALLEL" bash -c 'download_one "$@"' _ \
  || fail "some downloads failed (messages above); rerun the script to resume"

# 5. FASTA: downloaded once; UniProt serves only the current release
fasta_row="$(awk -F'\t' '$1 == "fasta"' "$MANIFEST")"
if [[ -n "$fasta_row" ]]; then
  IFS=$'\t' read -r _ fasta bytes checksum _ _ <<< "$fasta_row"
  [[ -f "$fasta" ]] || fail "$fasta (listed in $MANIFEST) is missing. UniProt serves only the current release; to use it, remove the fasta row from $MANIFEST and rerun (results will differ)"
  if [[ "${VERIFY:-0}" == "1" ]]; then verify "$fasta" "$bytes" "$checksum"; fi
  echo "present: $fasta"
else
  echo "downloading FASTA from UniProt"
  part="data/$FASTA_PREFIX.fasta.part"
  headers="$(mktemp)"
  count_headers="$(mktemp)"
  curl -fsS --retry 3 -D "$count_headers" -o /dev/null "$UNIPROT_COUNT" || fail "UniProt not reachable"
  curl -fsS --retry 3 -D "$headers" -o "$part" "$UNIPROT_STREAM" || fail "FASTA download failed"
  # value of one HTTP response header, case-insensitive name
  header() { tr -d '\r' < "$1" | awk -v k="$2" 'tolower($1) == tolower(k ":") { print $2 }'; }
  release="$(header "$headers" x-uniprot-release)"
  release_date="$(header "$headers" x-uniprot-release-date)"
  expected="$(header "$count_headers" x-total-results)"
  rm -f "$headers" "$count_headers"
  [[ -n "$release" && -n "$expected" ]] || fail "UniProt release or entry count missing in the response"
  entries="$(grep -c '^>' "$part")"
  [[ "$entries" == "$expected" ]] || fail "FASTA has $entries entries, UniProt reports $expected"
  not_sp="$(grep '^>' "$part" | grep -c -v '^>sp|' || true)"
  [[ "$not_sp" == "0" ]] || fail "FASTA has $not_sp headers not in UniProt Swiss-Prot format"
  decoys="$(grep -ciE '^>(rev_|decoy|##)' "$part" || true)"
  [[ "$decoys" == "0" ]] || fail "FASTA contains $decoys decoy entries"
  fasta="data/${FASTA_PREFIX}_${release}.fasta"
  mv "$part" "$fasta"
  printf 'fasta\t%s\t%s\tsha256:%s\t%s\t%s\n' "$fasta" "$(size_of "$fasta")" \
    "$(shasum -a 256 "$fasta" | cut -d' ' -f1)" "$UNIPROT_STREAM" \
    "UniProt release $release ($release_date), $entries entries, retrieved $(date -u +%Y-%m-%d)" \
    >> "$MANIFEST"
  echo "downloaded: $fasta ($entries entries, UniProt release $release)"
fi

# 6. summary
n_raw="$(find data -maxdepth 1 -name '*.raw' | wc -l | tr -d ' ')"
echo "done: $n_raw raw files in data/, FASTA $fasta, manifest $MANIFEST"
