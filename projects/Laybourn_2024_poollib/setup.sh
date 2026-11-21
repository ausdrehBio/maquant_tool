#!/usr/bin/env bash
# Install what run_server.sh needs. Every step is skipped when its result exists:
#   1. micromamba, only if neither micromamba nor conda is found
#      (to ~/.local/bin/micromamba; environments under $MAMBA_ROOT_PREFIX, default ~/micromamba)
#   2. the conda environment from environment.yml, recreated when environment.yml changed
#   3. the MaxQuant image localhost/maxquant:2.9.0 in podman, built from the repository's
#      Dockerfile and MaxQuant_v.2.9.0/ (the unzipped release, copied to the repository root by hand),
#      rebuilt when Dockerfile or .dockerignore changed; then a check that MaxQuant's native
#      libraries find every library they need in the image
# Usage:  bash projects/<project>/setup.sh   (run_server.sh calls it)
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

# package manager: micromamba (preferred), else conda; same lookup as in run_server.sh
find_manager() {
  if command -v micromamba > /dev/null 2>&1; then command -v micromamba
  elif [[ -x "$MICROMAMBA_LOCAL" ]]; then echo "$MICROMAMBA_LOCAL"
  elif command -v conda > /dev/null 2>&1; then command -v conda
  elif [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then echo "$CONDA_EXE"
  fi
}
# prefix of the environment $env_name from the manager's list, empty if it does not exist
env_prefix() {
  "$manager" env list 2>/dev/null \
    | awk -v n="$env_name" '{ p = $NF } p ~ /^\// { k = split(p, a, "/"); if (a[k] == n) { print p; exit } }' \
    || true
}

env_name="$(awk '/^name:/ { print $2; exit }' "$ENV_FILE")"
[[ -n "$env_name" ]] || die "no name: in projects/$PROJECT_NAME/environment.yml"

# 1. micromamba or conda
manager="$(find_manager)"
if [[ -n "$manager" ]]; then
  log "package manager: $(basename "$manager")"
else
  case "$(uname -s)-$(uname -m)" in
    Linux-x86_64)  platform=linux-64 ;;
    Linux-aarch64) platform=linux-aarch64 ;;
    Darwin-arm64)  platform=osx-arm64 ;;
    Darwin-x86_64) platform=osx-64 ;;
    *) die "no micromamba build for $(uname -sm); install micromamba or conda by hand" ;;
  esac
  for tool in curl tar bzip2; do
    command -v "$tool" > /dev/null || die "$tool not found (needed to install micromamba)"
  done
  log "installing micromamba to ~/.local/bin"
  tmp="$(mktemp -d)"
  curl -fsSL --retry 3 "https://micro.mamba.pm/api/micromamba/$platform/latest" \
    | tar -xj -C "$tmp" bin/micromamba || die "micromamba download failed"
  mkdir -p "$(dirname "$MICROMAMBA_LOCAL")"
  mv "$tmp/bin/micromamba" "$MICROMAMBA_LOCAL"
  rm -rf "$tmp"
  manager="$MICROMAMBA_LOCAL"
  log "micromamba $("$manager" --version) installed"
fi

# 2. conda environment; a checksum of environment.yml in the environment marks it as current
yml_sum="$(cksum < "$ENV_FILE" | cut -d' ' -f1)"
prefix="$(env_prefix)"
if [[ -n "$prefix" && "$(cat "$prefix/.environment_yml.cksum" 2>/dev/null || true)" == "$yml_sum" ]]; then
  log "conda environment $env_name present"
else
  if [[ -n "$prefix" ]]; then
    log "environment.yml changed: recreating conda environment $env_name"
    "$manager" env remove -y -n "$env_name" > /dev/null || die "removing $env_name failed"
  else
    log "creating conda environment $env_name"
  fi
  case "$(basename "$manager")" in
    micromamba) "$manager" create -q -y -n "$env_name" -f "$ENV_FILE" ;;
    *)          "$manager" env create -q -n "$env_name" -f "$ENV_FILE" ;;
  esac || die "creating conda environment $env_name failed"
  prefix="$(env_prefix)"
  [[ -n "$prefix" ]] || die "conda environment $env_name not found after creating it"
  echo "$yml_sum" > "$prefix/.environment_yml.cksum"
  log "conda environment $env_name created"
fi

# 3. MaxQuant image; a checksum of Dockerfile and .dockerignore in an image label marks it as current
command -v podman > /dev/null || die "podman not found; it must be installed (rootless), see maxquant_intro.md section 8"
build_sum="$(cat "$REPO_DIR/Dockerfile" "$REPO_DIR/.dockerignore" | cksum | cut -d' ' -f1)"
image_sum=""
if podman image exists "$IMAGE"; then
  image_sum="$(podman image inspect --format '{{ index .Labels "maxquant_tool.build_cksum" }}' "$IMAGE")" \
    || die "podman image inspect $IMAGE failed"
fi
if [[ "$image_sum" == "$build_sum" ]]; then
  log "MaxQuant image $IMAGE present"
else
  [[ -f "$REPO_DIR/MaxQuant_v.2.9.0/bin/MaxQuantCmd.dll" ]] \
    || die "MaxQuant_v.2.9.0/ missing at the repository root: copy the unzipped MaxQuant 2.9.0 release there"
  if podman image exists "$IMAGE"; then
    log "Dockerfile or .dockerignore changed since the MaxQuant image was built: rebuilding $IMAGE"
  else
    log "building the MaxQuant image $IMAGE (pulls the .NET 8 base image)"
  fi
  podman build --platform linux/amd64 --label "maxquant_tool.build_cksum=$build_sum" \
    -t "$IMAGE" -f "$REPO_DIR/Dockerfile" "$REPO_DIR" \
    || die "building the MaxQuant image failed"
  log "MaxQuant image $IMAGE built"
fi

# 4. MaxQuant's native libraries (*.so in bin/) need libraries of the image, for example
# libxgboost.so (MaxDIA library prediction) needs libgomp.so.1. A missing one only shows when a
# job loads the library, possibly hours into a run. ldd lists it as "<library> => not found".
# --ulimit nofile=-1:-1 as in run_<project>.sh
out="$(podman run --rm --pull never --platform linux/amd64 --ulimit nofile=-1:-1 \
  --entrypoint sh "$IMAGE" -c \
  'for f in /opt/maxquant/bin/*.so; do
     ldd "$f" | sed -n "s|^[[:space:]]*\([^ ]*\) => not found.*|MISSING: ${f##*/} needs \1|p"
   done' 2>&1)" \
  || die "checking the native libraries in $IMAGE failed; podman's output:"$'\n'"$(printf '%s\n' "$out" | tail -20)"
missing="$(printf '%s\n' "$out" | sed -n 's/^MISSING: //p')"
[[ -z "$missing" ]] \
  || die "libraries missing in $IMAGE (add their packages to the Dockerfile):"$'\n'"$missing"
log "MaxQuant's native libraries find their dependencies"
log "setup done"
