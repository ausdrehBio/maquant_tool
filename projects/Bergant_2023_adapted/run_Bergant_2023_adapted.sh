#!/usr/bin/env bash
# Reproducible MaxQuant run for one project, in podman.
# Usage:  ./run_<project>.sh [--clean] [--threads N] [--skip-disk-check] [--workdir DIR]
#   no option          first run; continues an unfinished run (crash, stop, SLURM time limit)
#                      if parameter/mqpar.xml and design.txt are unchanged
#   --clean            deletes previous results and MaxQuant working files, starts at job 1
#   --threads N        MaxQuant threads (default: $SLURM_CPUS_PER_TASK, else numThreads in mqpar.xml)
#   --skip-disk-check  starts despite too little free disk space
#   --workdir DIR      runs MaxQuant on copies of the input files in DIR, for example on a compute
#                      node's local disk when the project folder is on slow storage; combined/ is
#                      copied to results/, and the copies are deleted after a finished run
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

cd "$PROJECT_DIR"
fail() { echo "ERROR: $*" >&2; exit 1; }
# $work is the folder mounted as /work: the project folder, or its copy with --workdir.
# --ulimit nofile=-1:-1: open-files limit as high as this process may set it; a fixed
# default_ulimits in a server's containers.conf can exceed that, and the container fails to start
mq() { podman run --rm --pull never --platform linux/amd64 --security-opt label=disable \
         --ulimit nofile=-1:-1 -v "$work":/work "$IMAGE" "$@"; }
sha256() { if command -v sha256sum > /dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
# true if file $2 exists and has the size of file $1
same_size() { [[ -f "$2" && "$(wc -c < "$1" | tr -d ' ')" == "$(wc -c < "$2" | tr -d ' ')" ]]; }
# MaxQuant's working files in a data folder: combined/, a folder and an .index file per raw
# file, the mqpar copies. Raw and FASTA files stay.
remove_working_files() {
  local d="$1" raw
  rm -rf "$d/combined"
  for raw in "$d"/*.raw; do
    [[ -e "$raw" ]] || continue
    rm -rf "${raw%.raw}"
    rm -f "${raw%.raw}.index"
  done
  rm -f "$d/mqpar.xml" "$d/mqpar_run.xml"
}
# id of a job in results/dryrun.txt; names are compared without spaces and punctuation, since
# MaxQuant's proc files write "MS/MS main search" as MSMS_main_search. Only lines with a numeric
# id count: the header and podman warnings (servers without a login session) are skipped.
job_id() {
  awk -F'\t' -v n="$1" 'function norm(s) { gsub(/[^A-Za-z0-9]/, "", s); return s }
    $1 ~ /^[0-9]+$/ && norm($2) == norm(n) { hits++; id = $1 } END { if (hits == 1) print id }' results/dryrun.txt
}

clean="${CLEAN:-0}"
skip_disk="${SKIP_DISK_CHECK:-0}"
threads="${SLURM_CPUS_PER_TASK:-}"
workdir=""
while (( $# > 0 )); do
  case "$1" in
    --clean) clean=1 ;;
    --skip-disk-check) skip_disk=1 ;;
    --threads)
      [[ "${2:-}" =~ ^[1-9][0-9]*$ ]] || fail "--threads needs a positive number"
      threads="$2"
      shift ;;
    --workdir)
      [[ -n "${2:-}" ]] || fail "--workdir needs a folder"
      workdir="$2"
      shift ;;
    *) fail "unknown option: $1 (usage at the top of run_$PROJECT_NAME.sh)" ;;
  esac
  shift
done

# the folder mounted as /work. With --workdir, a folder of its own in DIR holds copies of the
# raw and FASTA files and all of MaxQuant's working files; results/ stays in the project folder.
work="$PROJECT_DIR"
copy_mode=0
if [[ -n "$workdir" ]]; then
  [[ -d "$workdir" && -w "$workdir" ]] || fail "--workdir $workdir is not a writable folder"
  work_name="maxquant_$(id -u)_$PROJECT_NAME"
  work="$(cd "$workdir" && pwd)/$work_name"
  copy_mode=1
  echo "work folder: $work_name in $workdir on $(uname -n)"
fi
FINISHED="$work/data/combined/proc/#runningTimes.txt"  # MaxQuant's list of finished jobs
continue_hint="run this script again to continue"
if (( copy_mode )); then
  continue_hint="run this script again with --workdir $workdir on $(uname -n) to continue"
fi

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
  [[ -f "$PROJECT_DIR/${p#/work/}" ]] || fail "mqpar.xml path does not exist in the project: $p"
done <<< "$mq_paths"
# the FASTA files of the search, relative to the project folder (data/ or elsewhere, for
# example an own database in parameter/); their checksums go into results/run_info.txt
fasta_files="$(sed -n 's|.*<fastaFilePath>/work/\(.*\)</fastaFilePath>.*|\1|p' "$MQPAR_HOST")"
mqpar_threads="$(sed -n 's|.*<numThreads>\([0-9]*\)</numThreads>.*|\1|p' "$MQPAR_HOST")"
[[ "$mqpar_threads" =~ ^[0-9]+$ ]] || fail "mqpar.xml needs exactly one <numThreads>"
threads="${threads:-$mqpar_threads}"

# 2. previous results: --clean deletes them, an unfinished run with unchanged parameters continues
# a manual --dryrun leaves an empty data/combined/proc; rmdir removes only empty folders
for d in "$PROJECT_DIR" "$work"; do
  rmdir "$d/data/combined/proc" "$d/data/combined" 2>/dev/null || true
done
resume=""   # job id to continue at; empty for a new run
if [[ "$clean" == "1" ]]; then
  if [[ -d results/combined || -d data/combined || -d "$work/data/combined" ]]; then
    echo "--clean: removing previous results and MaxQuant working files"
    rm -rf results/combined results/combined.part
    remove_working_files data
    if (( copy_mode )); then remove_working_files "$work/data"; fi
  fi
elif [[ -d results/combined ]]; then
  fail "a finished run exists in results/combined; rerun with --clean to delete it and start again"
elif (( copy_mode )) && [[ -d data/combined ]]; then
  fail "data/combined in the project folder holds an unfinished run; continue it without --workdir, or rerun with --clean to delete it"
elif [[ -d "$work/data/combined" ]]; then
  while IFS= read -r f; do
    [[ -n "$f" ]] || continue
    recorded="$(awk -v f="$f" '$2 == f { print $1 }' results/run_info.txt 2>/dev/null || true)"
    [[ -n "$recorded" && "$(sha256 "$f" | cut -d' ' -f1)" == "$recorded" ]] \
      || fail "data/combined holds an unfinished run, but $f differs from results/run_info.txt; rerun with --clean"
  done <<< "$(printf '%s\n' parameter/mqpar.xml parameter/design.txt "$fasta_files")"
  [[ -f results/dryrun.txt ]] || fail "data/combined holds an unfinished run without results/dryrun.txt; rerun with --clean"
  last=""
  if [[ -f "$FINISHED" ]]; then
    last="$(tail -n +2 "$FINISHED" | tr -d '\r' \
      | awk -F'\t' 'NF { name = $1 } END { sub(/ +$/, "", name); print name }')"
  fi
  # a failed job leaves "<job> <n>.error.txt" in proc and is listed as finished all the same
  failed=""
  if [[ -d "$work/data/combined/proc" ]]; then
    failed="$(find "$work/data/combined/proc" -maxdepth 1 -name '*.error.txt' \
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
  if [[ -d "$work/data/combined/proc" ]]; then
    find "$work/data/combined/proc" -maxdepth 1 -type f ! -name '#runningTimes.txt' -delete
  fi
elif (( copy_mode )) && [[ -f results/run_info.txt ]]; then
  # a run interrupted on another node left its work folder there
  echo "no unfinished run in the work folder on $(uname -n); starting a new run"
fi
mkdir -p "$RESULTS"

# 3. disk space: MaxQuant's working files and combined/ need up to 2.2 times the raw data; with
# --workdir, the work folder also needs room for the input files not yet copied there
raw_kb="$(du -ck data/*.raw | tail -1 | cut -f1)"
used_kb=0
if [[ -n "$resume" ]]; then
  used_kb="$(du -sk "$work/data/combined" | cut -f1)"
  for raw in "$work"/data/*.raw; do
    for w in "${raw%.raw}" "${raw%.raw}.index"; do
      if [[ -e "$w" ]]; then used_kb=$(( used_kb + $(du -sk "$w" | cut -f1) )); fi
    done
  done
fi
need_kb=$(( raw_kb * 22 / 10 - used_kb ))
if (( need_kb < 0 )); then need_kb=0; fi
copy_kb=0
if (( copy_mode )); then
  while IFS= read -r p; do
    if ! same_size "$PROJECT_DIR/${p#/work/}" "$work/${p#/work/}"; then
      copy_kb=$(( copy_kb + $(du -k "$PROJECT_DIR/${p#/work/}" | cut -f1) ))
    fi
  done <<< "$mq_paths"
  need_kb=$(( need_kb + copy_kb ))
  mkdir -p "$work/data"
fi
free_kb="$(df -Pk "$work/data" | awk 'NR == 2 { print $4 }')"
disk="$(awk -v n="$need_kb" -v r="$raw_kb" -v u="$used_kb" -v c="$copy_kb" -v f="$free_kb" \
  -v w="$copy_mode" 'BEGIN {
  printf "about %.1f GB needed (2.2 x %.1f GB raw data", n / 1e6, r / 1e6
  if (u > 0) printf " minus %.1f GB already written", u / 1e6
  if (c > 0) printf " plus %.1f GB of input files to copy", c / 1e6
  printf "), %.1f GB free%s", f / 1e6, (w ? " in the work folder" : "") }')"
if (( free_kb < need_kb )); then
  [[ "$skip_disk" == "1" ]] \
    || fail "not enough disk space: $disk; free space, or rerun with --skip-disk-check"
  echo "WARNING: not enough disk space: $disk; --skip-disk-check, starting anyway"
else
  echo "disk space: $disk"
fi
# copies of the raw and FASTA files listed in mqpar.xml; a copy with the size of the original
# counts as done, an interrupted one (.part) is replaced
if (( copy_kb > 0 )); then
  echo "copying $(awk -v c="$copy_kb" 'BEGIN { printf "%.1f", c / 1e6 }') GB of input files to the work folder ..."
  start=$SECONDS
  while IFS= read -r p; do
    src="$PROJECT_DIR/${p#/work/}"
    dst="$work/${p#/work/}"
    if same_size "$src" "$dst"; then continue; fi
    mkdir -p "$(dirname "$dst")"
    cp "$src" "$dst.part" || fail "copying ${p#/work/} to the work folder failed"
    mv "$dst.part" "$dst"
  done <<< "$mq_paths"
  echo "copied in $(( (SECONDS - start) / 60 )) min"
fi

# 4. this run's copy of mqpar.xml, and the parameter check (new runs only)
sed "s|<numThreads>[0-9]*</numThreads>|<numThreads>$threads</numThreads>|" parameter/mqpar.xml \
  > "$work/data/mqpar_run.xml"
echo "threads: $threads"
if [[ -z "$resume" ]]; then
  echo "dry run ..."
  mq "$MQPAR_RUN" --dryrun > "$RESULTS/dryrun.txt" 2>&1 || fail "dry run failed, see results/dryrun.txt"
fi
n_jobs="$(awk -F'\t' 'NR > 1 && $1 ~ /^[0-9]+$/' "$RESULTS/dryrun.txt" | wc -l | tr -d ' ')"

# 5. run info: written by a new run, extended by a continued one; with --workdir, the
# checksums are those of the copies MaxQuant read
if [[ -z "$resume" ]]; then
  {
    echo "project:       $PROJECT_NAME"
    echo "date (UTC):    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "image:         $IMAGE"
    echo "image id:      $(podman image inspect --format '{{.Id}}' "$IMAGE")"
    echo "maxquant:      $(mq --version 2>&1 | grep -m1 '^MaxQuantCmd')"
    echo "host:          $(uname -sm)"
    echo "threads:       $threads"
    if (( copy_mode )); then echo "inputs:        copies in a work folder (--workdir)"; fi
    echo
    echo "sha256 of inputs:"
    sha256 parameter/mqpar.xml parameter/design.txt
    (
      cd "$work"
      sha256 data/*.raw
      while IFS= read -r db; do
        if [[ -n "$db" && -f "$db" ]]; then sha256 "$db"; fi
      done <<< "$fasta_files"
    )
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
    --ulimit nofile=-1:-1 -v "$work":/work "$IMAGE" \
    "$MQPAR_RUN" --partial-processing="${resume:-1}" 2>&1 | tee -a "$LOG" \
    || fail "MaxQuant stopped before the end, see results/maxquant_run.log; $continue_hint"
fi

# 7. collect output: combined/ to results/; with --workdir it is copied, then the work folder deleted
if [[ -d "$work/data/combined" && ! -d results/combined ]]; then
  if (( copy_mode )); then
    echo "copying combined/ to results/ ..."
    rm -rf results/combined.part
    cp -R "$work/data/combined" results/combined.part \
      || fail "copying combined/ to results/ failed; the work folder is kept, $continue_hint"
    mv results/combined.part results/combined
  else
    mv data/combined results/combined
  fi
fi
[[ -d results/combined/txt ]] || fail "no results/combined/txt, check results/maxquant_run.log"
if (( copy_mode )); then
  rm -rf "$work"
  echo "work folder deleted"
fi
echo "done: results/combined/txt"
