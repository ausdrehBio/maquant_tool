#!/usr/bin/env bash
# Server entry point: environment, MaxQuant image, data, then the MaxQuant run of this project.
# On a fresh clone, with MaxQuant_v.2.9.0/ copied to the repository root, one command does all:
#   bash projects/<project>/run_server.sh [--dry-run | --download-only] [runner options]
#     --dry-run         checks only: environment, image, podman; no download, no MaxQuant
#     --download-only   everything except MaxQuant
#     runner options    passed on to run_<project>.sh: --clean, --threads N, --skip-disk-check
# Rerunning is safe: every step skips what exists, and run_<project>.sh continues an unfinished
# MaxQuant run (for example after a SLURM time limit) unless called with --clean.
# Needs: bash, curl, podman (rootless), internet. MaxQuant runs in the podman image.
# Details: maxquant_intro.md, section 8
set -euo pipefail

IMAGE="localhost/maxquant:2.9.0"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="$(basename "$PROJECT_DIR")"
REPO_DIR="$(cd "$PROJECT_DIR/../.." && pwd)"
ENV_FILE="$PROJECT_DIR/environment.yml"
MICROMAMBA_LOCAL="$HOME/.local/bin/micromamba"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

log() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "[$(date +%H:%M:%S)] ERROR: $*" >&2; exit 1; }

# package manager and environment lookup, same as in setup.sh
find_manager() {
  if command -v micromamba > /dev/null 2>&1; then command -v micromamba
  elif [[ -x "$MICROMAMBA_LOCAL" ]]; then echo "$MICROMAMBA_LOCAL"
  elif command -v conda > /dev/null 2>&1; then command -v conda
  elif [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then echo "$CONDA_EXE"
  fi
}
env_prefix() {
  "$manager" env list 2>/dev/null \
    | awk -v n="$env_name" '{ p = $NF } p ~ /^\// { k = split(p, a, "/"); if (a[k] == n) { print p; exit } }' \
    || true
}

# 1. arguments: --dry-run and --download-only are ours, everything else goes to the runner
mode=run
runner_args=()
passed=""
for arg in "$@"; do
  case "$arg" in
    --dry-run|--download-only)
      [[ "$mode" == run ]] || die "use only one of --dry-run and --download-only"
      mode="${arg#--}" ;;
    *)
      runner_args+=("$arg")
      passed="$passed $arg" ;;
  esac
done
log "projects/$PROJECT_NAME/run_server.sh, mode: $mode"

# 2. paths, before anything expensive
for f in setup.sh environment.yml "get_data_$PROJECT_NAME.sh" "run_$PROJECT_NAME.sh" \
         parameter/mqpar.xml parameter/design.txt; do
  [[ -f "$PROJECT_DIR/$f" ]] || die "projects/$PROJECT_NAME/$f missing"
done
for f in Dockerfile .dockerignore; do
  [[ -f "$REPO_DIR/$f" ]] || die "$f missing at the repository root"
done
command -v podman > /dev/null || die "podman not found; it must be installed (rootless), see maxquant_intro.md section 8"
log "paths ok"
# rootless podman keeps its runtime files in $XDG_RUNTIME_DIR. A SLURM job inherits the login
# session's /run/user/<uid>, which does not exist on a compute node without a login session.
# Unset, podman uses a folder in $TMPDIR or /tmp instead.
if [[ -n "${XDG_RUNTIME_DIR:-}" && ! -w "$XDG_RUNTIME_DIR" ]]; then
  log "XDG_RUNTIME_DIR $XDG_RUNTIME_DIR is missing or not writable (no login session): unset for podman"
  unset XDG_RUNTIME_DIR
fi

# 3. micromamba or conda, conda environment, MaxQuant image
log "setup.sh"
bash "$PROJECT_DIR/setup.sh" || die "setup.sh failed, see the messages above"

# 4. activate the environment by its prefix (conda activate is not available in batch jobs)
env_name="$(awk '/^name:/ { print $2; exit }' "$ENV_FILE")"
manager="$(find_manager)"
[[ -n "$manager" ]] || die "neither micromamba nor conda found after setup.sh"
prefix="$(env_prefix)"
[[ -n "$prefix" && -d "$prefix/bin" ]] || die "conda environment $env_name not found after setup.sh"
export CONDA_PREFIX="$prefix"
export PATH="$prefix/bin:$PATH"
log "conda environment $env_name active"

# 5. tools, and podman: the image starts and can write files that belong to the user
for tool in podman curl python3 shasum; do
  command -v "$tool" > /dev/null || die "$tool not found after activating $env_name"
done
# MaxQuantCmd prints its version on stderr and exits with 1; podman warnings can come before it
out="$(podman run --rm --pull never --platform linux/amd64 "$IMAGE" --version 2>&1 || true)"
version="$(printf '%s\n' "$out" | grep -m1 '^MaxQuantCmd' || true)"
[[ "$version" == "MaxQuantCmd 2.9.0"* ]] \
  || die "the MaxQuant image does not start; podman's output:"$'\n'"$(printf '%s\n' "$out" | tail -20)"
probe=".podman_write_test"
rm -f "$PROJECT_DIR/$probe"
podman run --rm --pull never --platform linux/amd64 --security-opt label=disable \
  -v "$PROJECT_DIR":/work --entrypoint sh "$IMAGE" -c "touch /work/$probe" \
  || die "podman cannot write to the mounted project folder"
if [[ ! -O "$PROJECT_DIR/$probe" ]]; then
  rm -f "$PROJECT_DIR/$probe"
  die "files written in the container do not belong to you; rootless podman is needed"
fi
rm -f "$PROJECT_DIR/$probe"
log "podman ok: $version, writes to the project folder"
if [[ "$mode" == dry-run ]]; then
  log "dry run done: ready to download and run"
  exit 0
fi

# 6. data, skips verified files
log "get_data_$PROJECT_NAME.sh"
bash "$PROJECT_DIR/get_data_$PROJECT_NAME.sh" || die "download failed, see the messages above; rerun to resume"
if [[ "$mode" == download-only ]]; then
  log "download done; MaxQuant not started (--download-only)"
  exit 0
fi

# 7. MaxQuant: the runner replaces this process, so signals reach it directly
log "run_$PROJECT_NAME.sh$passed"
exec bash "$PROJECT_DIR/run_$PROJECT_NAME.sh" ${runner_args[@]+"${runner_args[@]}"}
