#!/usr/bin/env bash
# Reproducible MaxQuant run for one project, in podman.
# Usage:  ./run_<project>.sh [--clean] [--threads N] [--skip-disk-check]
#   no option          first run; continues an unfinished run (crash, stop, SLURM time limit)
#                      if parameter/mqpar.xml and design.txt are unchanged
#   --clean            deletes previous results and MaxQuant working files, starts at job 1
#   --threads N        MaxQuant threads (default: $SLURM_CPUS_PER_TASK, else numThreads in mqpar.xml)
#   --skip-disk-check  starts despite too little free disk space
#   CLEAN=1 and SKIP_DISK_CHECK=1 work like --clean and --skip-disk-check.
# Details: maxquant_intro.md, section 5
set -euo pipefail

IMAGE="localhost/maxquant:2.9.0"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="$(basename "$PROJECT_DIR")"
CONTAINER="mq_${PROJECT_NAME}"
MQPAR_HOST="$PROJECT_DIR/parameter/mqpar.xml"
MQPAR_RUN="/work/data/mqpar_run.xml"            # mqpar.xml with this run's thread count
RESULTS="$PROJECT_DIR/results"
LOG="$RESULTS/maxquant_run.log"
FINISHED="data/combined/proc/#runningTimes.txt"  # MaxQuant's list of finished jobs

cd "$PROJECT_DIR"
fail() { echo "ERROR: $*" >&2; exit 1; }
mq() { podman run --rm --pull never --platform linux/amd64 --security-opt label=disable \
         -v "$PROJECT_DIR":/work "$IMAGE" "$@"; }
sha256() { if command -v sha256sum > /dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
# id of a job in results/dryrun.txt; names are compared without spaces and punctuation, since
# MaxQuant's proc files write "MS/MS main search" as MSMS_main_search
job_id() {
  awk -F'\t' -v n="$1" 'function norm(s) { gsub(/[^A-Za-z0-9]/, "", s); return s }
    NR > 1 && norm($2) == norm(n) { hits++; id = $1 } END { if (hits == 1) print id }' results/dryrun.txt
}

clean="${CLEAN:-0}"
skip_disk="${SKIP_DISK_CHECK:-0}"
threads="${SLURM_CPUS_PER_TASK:-}"
while (( $# > 0 )); do
  case "$1" in
    --clean) clean=1 ;;
    --skip-disk-check) skip_disk=1 ;;
    --threads)
      [[ "${2:-}" =~ ^[1-9][0-9]*$ ]] || fail "--threads needs a positive number"
      threads="$2"
      shift ;;
    *) fail "unknown option: $1 (usage at the top of run_$PROJECT_NAME.sh)" ;;
  esac
  shift
done

# 1. prerequisites
command -v podman > /dev/null || fail "podman not found"
podman image exists "$IMAGE" || fail "image $IMAGE missing, see maxquant_intro.md section 1"
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
mqpar_threads="$(sed -n 's|.*<numThreads>\([0-9]*\)</numThreads>.*|\1|p' "$MQPAR_HOST")"
[[ "$mqpar_threads" =~ ^[0-9]+$ ]] || fail "mqpar.xml needs exactly one <numThreads>"
threads="${threads:-$mqpar_threads}"

# 2. previous results: --clean deletes them, an unfinished run with unchanged parameters continues
# a manual --dryrun leaves an empty data/combined/proc; rmdir removes only empty folders
rmdir data/combined/proc data/combined 2>/dev/null || true
resume=""   # job id to continue at; empty for a new run
if [[ "$clean" == "1" ]]; then
  if [[ -d results/combined || -d data/combined ]]; then
    echo "--clean: removing previous results and MaxQuant working files"
    rm -rf results/combined data/combined
    for raw in data/*.raw; do
      rm -rf "${raw%.raw}"
      rm -f "${raw%.raw}.index"
    done
    rm -f data/mqpar.xml data/mqpar_run.xml
  fi
elif [[ -d results/combined ]]; then
  fail "a finished run exists in results/combined; rerun with --clean to delete it and start again"
elif [[ -d data/combined ]]; then
  for f in parameter/mqpar.xml parameter/design.txt; do
    recorded="$(awk -v f="$f" '$2 == f { print $1 }' results/run_info.txt 2>/dev/null || true)"
    [[ -n "$recorded" && "$(sha256 "$f" | cut -d' ' -f1)" == "$recorded" ]] \
      || fail "data/combined holds an unfinished run, but $f differs from results/run_info.txt; rerun with --clean"
  done
  [[ -f results/dryrun.txt ]] || fail "data/combined holds an unfinished run without results/dryrun.txt; rerun with --clean"
  last=""
  if [[ -f "$FINISHED" ]]; then
    last="$(tail -n +2 "$FINISHED" | tr -d '\r' \
      | awk -F'\t' 'NF { name = $1 } END { sub(/ +$/, "", name); print name }')"
  fi
  # a failed job leaves "<job> <n>.error.txt" in proc and is listed as finished all the same
  failed=""
  if [[ -d data/combined/proc ]]; then
    failed="$(find data/combined/proc -maxdepth 1 -name '*.error.txt' \
      | sed -e 's|.*/||' -e 's/ [0-9]*\.error\.txt$//' | head -1 || true)"
  fi
  if [[ -n "$failed" ]]; then
    resume="$(job_id "$failed")"
    [[ -n "$resume" ]] || fail "failed job \"$failed\" not found in results/dryrun.txt; rerun with --clean"
    echo "unfinished run: job $resume ($failed) had failed and is repeated, error in results/maxquant_run.log"
  elif [[ -z "$last" ]]; then
    resume=1
    echo "unfinished run: no job had finished, starting at job 1"
  else
    id="$(job_id "$last")"
    [[ -n "$id" ]] || fail "last finished job \"$last\" not found in results/dryrun.txt; rerun with --clean"
    resume=$(( id + 1 ))
    echo "unfinished run: continuing at job $resume (last finished: $last)"
  fi
  # markers and error files of the interrupted job; MaxQuant writes new ones
  if [[ -d data/combined/proc ]]; then
    find data/combined/proc -maxdepth 1 -type f ! -name '#runningTimes.txt' -delete
  fi
fi
mkdir -p "$RESULTS"

# 3. disk space: MaxQuant's working files and combined/ need up to 2.2 times the raw data
raw_kb="$(du -ck data/*.raw | tail -1 | cut -f1)"
used_kb=0
if [[ -n "$resume" ]]; then
  used_kb="$(du -sk data/combined | cut -f1)"
  for raw in data/*.raw; do
    for w in "${raw%.raw}" "${raw%.raw}.index"; do
      if [[ -e "$w" ]]; then used_kb=$(( used_kb + $(du -sk "$w" | cut -f1) )); fi
    done
  done
fi
need_kb=$(( raw_kb * 22 / 10 - used_kb ))
if (( need_kb < 0 )); then need_kb=0; fi
free_kb="$(df -Pk data | awk 'NR == 2 { print $4 }')"
disk="$(awk -v n="$need_kb" -v r="$raw_kb" -v u="$used_kb" -v f="$free_kb" 'BEGIN {
  printf "about %.1f GB needed (2.2 x %.1f GB raw data", n / 1e6, r / 1e6
  if (u > 0) printf " minus %.1f GB already written", u / 1e6
  printf "), %.1f GB free", f / 1e6 }')"
if (( free_kb < need_kb )); then
  [[ "$skip_disk" == "1" ]] \
    || fail "not enough disk space: $disk; free space, or rerun with --skip-disk-check"
  echo "WARNING: not enough disk space: $disk; --skip-disk-check, starting anyway"
else
  echo "disk space: $disk"
fi

# 4. this run's copy of mqpar.xml, and the parameter check (new runs only)
sed "s|<numThreads>[0-9]*</numThreads>|<numThreads>$threads</numThreads>|" parameter/mqpar.xml \
  > data/mqpar_run.xml
echo "threads: $threads"
if [[ -z "$resume" ]]; then
  echo "dry run ..."
  mq "$MQPAR_RUN" --dryrun > "$RESULTS/dryrun.txt" 2>&1 || fail "dry run failed, see results/dryrun.txt"
fi
n_jobs="$(awk -F'\t' 'NR > 1 && $1 ~ /^[0-9]+$/' "$RESULTS/dryrun.txt" | wc -l | tr -d ' ')"

# 5. run info: written by a new run, extended by a continued one
if [[ -z "$resume" ]]; then
  {
    echo "project:       $PROJECT_NAME"
    echo "date (UTC):    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "image:         $IMAGE"
    echo "image id:      $(podman image inspect --format '{{.Id}}' "$IMAGE")"
    echo "maxquant:      $(mq --version 2>&1 | head -1)"
    echo "host:          $(uname -sm)"
    echo "threads:       $threads"
    echo
    echo "sha256 of inputs:"
    sha256 parameter/mqpar.xml parameter/design.txt data/*.raw
    for db in data/*.fasta data/*.fa data/*.faa; do
      if [[ -f "$db" ]]; then sha256 "$db"; fi
    done
  } > "$RESULTS/run_info.txt" || true
else
  {
    echo
    echo "continued (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ) at job $resume of $n_jobs," \
      "threads $threads, host $(uname -sm), image id $(podman image inspect --format '{{.Id}}' "$IMAGE")"
  } >> "$RESULTS/run_info.txt" || true
fi

# 6. MaxQuant; a container left over from an interrupted run would block the name
if [[ -n "$resume" ]] && (( resume > n_jobs )); then
  echo "all $n_jobs MaxQuant jobs had finished"
else
  if podman container exists "$CONTAINER"; then
    [[ "$(podman container inspect --format '{{.State.Running}}' "$CONTAINER")" != "true" ]] \
      || fail "container $CONTAINER is running already, see podman ps"
    podman rm -f "$CONTAINER" > /dev/null
  fi
  if [[ -z "$resume" ]]; then
    : > "$LOG"
    echo "MaxQuant running, log: results/maxquant_run.log"
  else
    echo "=== continued at job $resume ($(date -u +%Y-%m-%dT%H:%M:%SZ))" >> "$LOG"
    echo "MaxQuant continuing at job $resume of $n_jobs, log: results/maxquant_run.log"
  fi
  podman run --rm --name "$CONTAINER" --pull never --platform linux/amd64 --security-opt label=disable \
    -v "$PROJECT_DIR":/work "$IMAGE" "$MQPAR_RUN" --partial-processing="${resume:-1}" 2>&1 | tee -a "$LOG" \
    || fail "MaxQuant stopped before the end, see results/maxquant_run.log; run this script again to continue"
fi

# 7. collect output
if [[ -d data/combined && ! -d results/combined ]]; then
  mv data/combined results/combined
fi
[[ -d results/combined/txt ]] || fail "no results/combined/txt, check results/maxquant_run.log"
echo "done: results/combined/txt"
