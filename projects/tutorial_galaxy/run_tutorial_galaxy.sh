#!/usr/bin/env bash
# Reproducible MaxQuant run for one project.
# Usage:  ./run_<project>.sh          first run
#         CLEAN=1 ./run_<project>.sh  rerun, deletes previous results
# Details: maxquant_intro.md, section 5
set -euo pipefail

IMAGE="maxquant:2.9.0"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="$(basename "$PROJECT_DIR")"
CONTAINER="mq_${PROJECT_NAME}"
MQPAR_HOST="$PROJECT_DIR/parameter/mqpar.xml"
MQPAR="/work/parameter/mqpar.xml"
RESULTS="$PROJECT_DIR/results"
LOG="$RESULTS/maxquant_run.log"

cd "$PROJECT_DIR"
mq() { docker run --rm --platform linux/amd64 -v "$PROJECT_DIR":/work "$IMAGE" "$@"; }
fail() { echo "ERROR: $*" >&2; exit 1; }

# 1. prerequisites
command -v docker > /dev/null || fail "docker not found"
docker image inspect "$IMAGE" > /dev/null 2>&1 || fail "image $IMAGE missing, see maxquant_intro.md section 1"
[[ -f parameter/mqpar.xml ]] || fail "parameter/mqpar.xml missing"
[[ -f parameter/design.txt ]] || fail "parameter/design.txt missing"
compgen -G "data/*.raw" > /dev/null || fail "no .raw files in data/"
# raw and FASTA paths in mqpar.xml must be container paths to existing files
mq_paths="$(sed -n -e 's|.*<fastaFilePath>\(.*\)</fastaFilePath>.*|\1|p' \
  -e '/<filePaths>/,/<\/filePaths>/s|.*<string>\(.*\)</string>.*|\1|p' "$MQPAR_HOST")"
[[ -n "$mq_paths" ]] || fail "no raw or FASTA paths found in mqpar.xml"
while IFS= read -r p; do
  [[ "$p" == /work/* ]] || fail "mqpar.xml path is not a /work/... path: $p"
  [[ -e "$PROJECT_DIR/${p#/work/}" ]] || fail "mqpar.xml path does not exist in the project: $p"
done <<< "$mq_paths"

# 2. previous results
if [[ -d results/combined || -d data/combined ]]; then
  [[ "${CLEAN:-0}" == "1" ]] || fail "previous results exist, rerun with CLEAN=1 to delete them"
  echo "CLEAN=1: removing previous results and MaxQuant working files"
  rm -rf results/combined data/combined
  for raw in data/*.raw; do
    rm -rf "${raw%.raw}"
    rm -f "${raw%.raw}.index"
  done
  rm -f data/mqpar.xml
fi
mkdir -p "$RESULTS"

# 3. parameter check
echo "dry run ..."
mq "$MQPAR" --dryrun > "$RESULTS/dryrun.txt" 2>&1 || fail "dry run failed, see results/dryrun.txt"

# 4. run info
{
  echo "project:       $PROJECT_NAME"
  echo "date (UTC):    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "image:         $IMAGE"
  echo "image id:      $(docker image inspect --format '{{.Id}}' "$IMAGE")"
  echo "maxquant:      $(mq --version 2>&1 | head -1)"
  echo "host:          $(uname -sm)"
  echo
  echo "sha256 of inputs:"
  shasum -a 256 parameter/mqpar.xml parameter/design.txt data/*.raw
  for db in data/*.fasta data/*.fa data/*.faa; do
    [[ -f "$db" ]] && shasum -a 256 "$db"
  done
} > "$RESULTS/run_info.txt" || true

# 5. MaxQuant
echo "MaxQuant running, log: $LOG"
docker run --rm --name "$CONTAINER" --platform linux/amd64 \
  -v "$PROJECT_DIR":/work "$IMAGE" "$MQPAR" 2>&1 | tee "$LOG"

# 6. collect output
if [[ -d data/combined && ! -d results/combined ]]; then
  mv data/combined results/combined
fi
[[ -d results/combined/txt ]] || fail "no results/combined/txt, check $LOG"
echo "done: $RESULTS/combined/txt"