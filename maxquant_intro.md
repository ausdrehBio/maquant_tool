# MaxQuant in this repository: workflow reference

This file describes how MaxQuant is run in this repository. Read it fully before processing a dataset.

## Overview

- MaxQuant 2.9.0 runs inside a local podman image, `localhost/maxquant:2.9.0`.
- Each dataset lives in its own folder under `projects/`.
- A MaxQuant run needs two input files besides the data: `mqpar.xml` (all analysis parameters) and `design.txt` (sample assignment of the raw files). Both live in the project's `parameter/` folder.
- Every project has a runner script, `run_<project_name>.sh`, which performs the MaxQuant run. **Full runs are only started through this script**, so every result in the repository is reproducible with one command.
- On a Linux server, `run_server.sh` in the project folder sets up everything on a fresh clone and then calls the runner script (section 8).
- If parameter files or the runner script are missing, create them as described below.
- **Never start a full MaxQuant run without the user's confirmation of the design and the key parameters.** A run takes from minutes to many hours, and a wrong parameter means repeating it.

## Repository layout

```
repo/
├── maxquant_intro.md
├── Dockerfile                     image definition (section 1)
├── .dockerignore
├── MaxQuant_v.2.9.0/              MaxQuant release (local only, not in git)
└── projects/
    └── <project_name>/
        ├── run_<project_name>.sh  runner script: reproduces the MaxQuant run
        ├── get_data_<project_name>.sh  optional: downloads public data (section 6)
        ├── data_manifest.tsv      optional: sizes and checksums of the downloaded files
        ├── make_parameters_<project_name>.sh  optional: creates parameter/ from a deposited mqpar.xml
        ├── run_server.sh          optional: server entry point (section 8)
        ├── setup.sh               optional: micromamba, conda environment, image (section 8)
        ├── environment.yml        optional: host tools for run_server.sh (section 8)
        ├── data/                  .raw files and FASTA file(s), input only
        ├── parameter/             mqpar.xml, design.txt
        ├── results/               MaxQuant output, run log, run info
        └── analysis_of_results/   downstream scripts, plots, tables
```

Roles of the folders:

| Folder | Content | Written by |
|---|---|---|
| `data/` | raw files, FASTA | the user. Raw and FASTA files are never modified. MaxQuant adds working files (per-raw-file folders, `.index` files, a copy of mqpar.xml, `combined/` during the run), the runner script adds `mqpar_run.xml` (section 5) |
| `parameter/` | `mqpar.xml`, `design.txt` | generated once, then edited only on purpose |
| `results/` | `combined/`, `maxquant_run.log`, `run_info.txt` | only the runner script |
| `analysis_of_results/` | QC, plots, statistics | analysis scripts, which read `results/` and never modify it (one exception: PTXQC, section 7) |

Reference project: `projects/tutorial_galaxy/` (Galaxy Training Network tutorial "Label-free data analysis using MaxQuant", two human serum samples).

---

## 1. The podman image

### Why a container

MaxQuant's GUI runs only on Windows. The command-line version runs on .NET 8, but on native macOS with Apple Silicon, Thermo's raw file reader (`ThermoFisher.CommonCore.RawFileReader.dll`) fails to load. Do not attempt native runs on the Mac. Inside a Linux x64 container, which is MaxQuant's officially supported platform, it works. .NET exists only inside the image, so the host needs podman but no .NET.

The container engine is podman, on the Mac and on servers. On a Linux x86_64 server, rootless podman runs the image natively. On the Mac, podman runs containers in a Linux VM, the podman machine, and emulates x64 there. Create the machine once with enough CPUs and memory, since MaxQuant needs about 3 GB per thread:

```bash
brew install podman
podman machine init --cpus 8 --memory 24576    # memory in MB
podman machine start
```

The project folders must be inside the home folder, which the machine mounts by default.

### Check that the image exists

```bash
podman image exists localhost/maxquant:2.9.0 && echo "image present"
```

Podman stores locally built images under the `localhost/` prefix. Scripts use the full name together with `--pull never`, so podman never looks for the image in a registry.

### Rebuild only if it is missing

Requires the unzipped release in `MaxQuant_v.2.9.0/` at the repository root. On a server, `setup.sh` builds the image (section 8). If `Dockerfile` or `.dockerignore` are missing at the root, create them with this content.

`Dockerfile`:

```dockerfile
FROM --platform=linux/amd64 mcr.microsoft.com/dotnet/runtime:8.0
COPY MaxQuant_v.2.9.0 /opt/maxquant
WORKDIR /work
ENTRYPOINT ["dotnet", "/opt/maxquant/bin/MaxQuantCmd.dll"]
```

`.dockerignore` (keeps raw data out of the build; podman reads it like Docker):

```
projects/
.git/
```

```bash
podman build --platform linux/amd64 -t localhost/maxquant:2.9.0 .
podman run --rm --pull never --platform linux/amd64 localhost/maxquant:2.9.0 --version   # expect: MaxQuantCmd 2.9.0.0 (on stderr)
```

The build pulls the .NET 8 base image from `mcr.microsoft.com` once, so it needs internet access.

---

## 2. How to run the container

### Mount convention

The **project folder** is mounted at `/work`. Every path MaxQuant sees must be a container path:

| On the host | In the container |
|---|---|
| `projects/<name>/data/` | `/work/data/` |
| `projects/<name>/parameter/` | `/work/parameter/` |
| `projects/<name>/results/` | `/work/results/` |

Host paths (`/Users/...`) inside `mqpar.xml` or in MaxQuant arguments will fail. Because all paths are relative to the project mount, a project folder stays valid when moved or copied.

### Base command for short tasks

For generating, editing and checking parameter files (sections 3 and 4). Run from the repository root; the later examples use this `mq` function:

```bash
P=projects/<project_name>
mq() { podman run --rm --pull never --platform linux/amd64 --security-opt label=disable \
         -v "$(pwd)/$P":/work localhost/maxquant:2.9.0 "$@"; }
mq <MaxQuant arguments>
```

| Option | Why |
|---|---|
| `--pull never` | use the local image, never contact a registry |
| `--platform linux/amd64` | MaxQuant's platform; emulated on the Mac, native on x86_64 servers |
| `--security-opt label=disable` | on servers with SELinux, the mounted project folder is otherwise not accessible; no effect elsewhere |

`ENTRYPOINT` already calls `MaxQuantCmd.dll`, so only MaxQuant's arguments follow the image name. All MaxQuant options: `mq --help`.

Rootless podman maps the container's root user to the user who starts it, so files MaxQuant writes belong to that user.

Full runs go through the runner script (section 5), not through ad hoc `podman run` calls.

---

## 3. mqpar.xml

### Purpose

`mqpar.xml` is MaxQuant's parameter file. It holds every setting of an analysis in one XML file: input files, FASTA databases, sample assignment, enzyme, modifications, mass tolerances, FDR thresholds, quantification method and output options. The command-line version takes its settings only from this file. The Windows GUI and the Galaxy MaxQuant tool also write an mqpar.xml before each run, so it is the common exchange format. MaxQuant copies it into its output, which makes it the record of how the data were processed.

### Properties

- Plain XML. Each tag is one parameter, for example `<maxMissedCleavages>2</maxMissedCleavages>`. Lists use `<string>` children.
- Contains absolute paths (FASTA and raw files). In this repository these must be `/work/...` paths.
- **Version-bound.** An mqpar.xml runs only with the MaxQuant version that created it. Check `<maxQuantVersion>`.
- Parameters inside `<parameterGroups>` apply per group of raw files. Normally all files share one group.
- The same modification names appear in several blocks (searched modifications, modifications used for quantification). Edit only the block you intend to change.

### Generate a new one

`design.txt` must exist first (section 4).

```bash
mq --create \
  --newMqpar /work/parameter/mqpar.xml \
  --LCMSType ST --instrumentType TO \
  --pathFasta /work/data/<database>.fasta \
  --pathRawFileFolder /work/data \
  --experimentalDesign /work/parameter/design.txt \
  --useLFQ \
  --numThreads 4
```

| Option | Meaning | Values |
|---|---|---|
| `--LCMSType` | acquisition and quantification type | `ST` standard DDA, `R2`/`R3` TMT reporter ions MS2/MS3, `DIA` MaxDIA, `TD` TIMS-DDA, more in `--help` |
| `--instrumentType` | sets instrument-specific defaults (mass tolerances, peak detection) | `TO` Orbitrap, `TA` Astral, `BT` Bruker TIMS, `BQ` Bruker Q-TOF, `SC` Sciex, `AQ` Agilent, `WQ` Waters |
| `--pathFasta` | one or several FASTA files, space-separated | |
| `--pathFastaFolder` | add all FASTA files in a folder | |
| `--pathRawFileFolder` | all raw files in this folder are included | |
| `--useLFQ` | label-free quantification (MaxLFQ) | |
| `--useMBR` | match between runs | |
| `--numThreads` | CPU threads, at or below `podman info --format '{{.Host.CPUs}}'`; the runner script can override it per run (`--threads`) | |

Only standard DDA on Thermo Orbitrap data has been tested with this image.

`--create` takes the raw files from `--pathRawFileFolder`, not from `design.txt`: a raw file missing from `design.txt` gets an empty experiment name without a warning, and a `design.txt` row without a raw file is dropped. `--create` and `--dryrun` do not open raw files (the FASTA files must exist), so parameters can be prepared against empty placeholder files named after `design.txt` before a download has finished. Example: `projects/Bergant_2023/make_parameters_Bergant_2023.sh`.

### Output location

The mqpar.xml of MaxQuant 2.9.0 has no parameter for the location of `combined/`. MaxQuant writes `combined/` into the raw file folder, `data/`, together with per-raw-file working folders and `.index` files next to the raw files. After a successful run, the runner script moves `data/combined/` to `results/combined/` (section 5).

### Default search parameters

These defaults suit standard bottom-up experiments. Keep them unless the user specifies otherwise.

| Parameter | Tag | Default |
|---|---|---|
| Enzyme | `<enzymes>` | Trypsin/P (cleaves after K/R, also before P) |
| Max. missed cleavages | `<maxMissedCleavages>` | 2 |
| Fixed modifications | `<fixedModifications>` | Carbamidomethyl (C) |
| Variable modifications | `<variableModifications>` | Oxidation (M), Acetyl (Protein N-term) |
| Min. peptide length | `<minPeptideLength>` | 7 |
| Min. unique peptides | `<minUniquePeptides>` | 0 |
| PSM and protein FDR | | 0.01 |
| Identifier parse rule | `<identifierParseRule>` | UniProt accession, regex below |
| Description parse rule | `<descriptionParseRule>` | `>(.*)` |

Default identifier parse rule, as it appears in mqpar.xml:

```xml
<identifierParseRule>>[^|]*\|(.*?)\|</identifierParseRule>
```

It takes the accession between the first two `|` of a UniProt header (`>sp|P02768|ALBU_HUMAN ...` gives `P02768`). If FASTA headers are not in UniProt format, the rule must be adapted, otherwise protein IDs are wrong.

The `tutorial_galaxy` project deviates from these defaults on purpose to reproduce the tutorial (1 missed cleavage, no variable modifications, min. 1 unique peptide, description rule `>(.*) OS`, LFQ min. ratio count 2 instead of the `--create` default 1). Do not copy these settings to other projects.

### Editing

Single values:

```bash
mq --changeParameter /work/parameter/mqpar.xml numThreads 8
```

Use the full parameter name. `--changeParameter` matches partial names, and if a name matches several parameters it asks interactively, which a plain `podman run` (without `-it`) cannot answer.

Direct edits with `sed` (macOS syntax `-i ''`). Always back up first and diff afterwards:

```bash
cp $P/parameter/mqpar.xml $P/parameter/mqpar_default.xml
sed -i '' 's|<maxMissedCleavages>2</maxMissedCleavages>|<maxMissedCleavages>1</maxMissedCleavages>|' $P/parameter/mqpar.xml
diff $P/parameter/mqpar_default.xml $P/parameter/mqpar.xml
```

Removing entries from one list block only:

```bash
sed -i '' -e '/<variableModifications>/,/<\/variableModifications>/{' -e '/<string>/d' -e '}' $P/parameter/mqpar.xml
```

### Reusing an existing mqpar.xml

1. Check `<maxQuantVersion>`. It must be 2.9.0.x for this image. Otherwise generate a new file and transfer the non-default settings. Example: `projects/Bergant_2023/` (authors' mqpar.xml of MaxQuant 1.5.6.2, transferred by `make_parameters_Bergant_2023.sh`).
2. Check that all paths are `/work/...` paths and point to existing files.
3. Check that `<filePaths>` lists exactly the raw files of this dataset and `<experiments>` the intended sample names, in the same order.
4. To update FASTA and raw paths in one step (arguments: new file, FASTA folder, raw folder):

```bash
mq /work/parameter/mqpar_template.xml --changeFolder /work/parameter/mqpar.xml /work/data /work/data
```

If the number or names of raw files changed, regenerating the file with a new `design.txt` is less error-prone than editing `<filePaths>` and `<experiments>` by hand.

### Verify

```bash
grep -n -A3 -E "<maxQuantVersion>|<fastaFilePath>|<filePaths>|<experiments>" $P/parameter/mqpar.xml
mq /work/parameter/mqpar.xml --dryrun
```

`--dryrun` parses the file and prints the job list without processing. A malformed file fails here within seconds. The runner script repeats this check before every run, and before it checks that every raw and FASTA path in mqpar.xml exists.

---

## 4. design.txt

### Purpose

The experimental design file assigns raw files to samples ("experiments"). MaxQuant treats raw files and samples separately because one sample can span several raw files, for example after fractionation. Each experiment name becomes its own intensity column in the output tables (`Intensity <name>`, `LFQ intensity <name>`). Without a correct assignment, raw files can end up pooled into one column, and per-sample quantification is lost.

### Format

Tab-separated, with a header row and one row per raw file:

| Column | Content |
|---|---|
| `Name` | raw file name without `.raw`, exactly as on disk |
| `Fraction` | fraction number, empty if not fractionated |
| `Experiment` | sample name, becomes the output column name |
| `PTM` | `True` for PTM-enriched runs (for example phosphopeptide enrichment), otherwise `False` |

Unfractionated example:

```
Name	Fraction	Experiment	PTM
sample1		control_1	False
sample2		treated_1	False
```

Fractionated example (one sample measured as three raw files):

```
Name	Fraction	Experiment	PTM
s1_f01	1	sample1	False
s1_f02	2	sample1	False
s1_f03	3	sample1	False
```

### Properties and pitfalls

- **Real tabs are required.** Write the file with `printf` and check with `cat -t` (tabs appear as `^I`):

```bash
printf 'Name\tFraction\tExperiment\tPTM\nsample1\t\tsample1\tFalse\n' > $P/parameter/design.txt
cat -t $P/parameter/design.txt
```

- Contains no paths, so it needs no changes for the container.
- Experiment names should be short and without spaces. They become column headers.
- Every raw file in `data/` must appear in the file.

---

## 5. The runner script

### Purpose

`run_<project_name>.sh` is the single, documented way to run MaxQuant for a project. It fixes everything that is not stored in `mqpar.xml`: which image is used, how the project is mounted, where output goes, how earlier results are handled, and what is logged. Together with `data/` and `parameter/`, it makes a result reproducible with one command: anyone with the image can rerun the project and get the same output.

### What it does

1. Checks prerequisites: podman available, image present, `mqpar.xml` and `design.txt` present, every raw and FASTA path in `mqpar.xml` is a `/work/...` path to an existing file.
2. Handles earlier results:
   - A finished run (`results/combined/`) is never overwritten without `--clean`.
   - An unfinished run (`data/combined/` left by a failure, a stop or a SLURM time limit) is continued if `mqpar.xml` and `design.txt` have the checksums recorded in `results/run_info.txt`. If the parameters changed, it stops and asks for `--clean`. MaxQuant restarts with `--partial-processing` (job numbers from `results/dryrun.txt`):
     - after a stop, at the job after the last finished one (`data/combined/proc/#runningTimes.txt`);
     - after a failure, at the failed job, which `proc/` marks with `<job> <n>.error.txt`. MaxQuant lists a failed job in `#runningTimes.txt` like a finished one, so this check comes first. Proc file names differ from the job names (`MSMS_main_search` for "MS/MS main search"), so names are compared without spaces and punctuation.

     Markers and error files of the interrupted job are removed before the restart. The error itself stays in `results/maxquant_run.log`.
   - `--clean` deletes `results/combined/` and MaxQuant's working files in `data/` (per-raw-file folders, `.index` files, the mqpar copies). It never touches raw or FASTA files.
3. Checks free disk space: MaxQuant's working files and `combined/` need up to 2.2 times the size of the raw files, minus what a continued run has already written. With less free space it stops, unless called with `--skip-disk-check`. A run that fills the disk fails hours in.
4. Writes `data/mqpar_run.xml`, a copy of `mqpar.xml` with this run's thread count: `--threads N`, otherwise the CPUs of the SLURM job (`SLURM_CPUS_PER_TASK`), otherwise `numThreads` from `mqpar.xml`. The thread count changes only the speed, not the results. A new run then checks the parameters with `--dryrun`.
5. Writes `results/run_info.txt`: date, image ID, MaxQuant version, threads, SHA-256 checksums of `mqpar.xml`, `design.txt`, raw files and FASTA files. The checksums prove which exact inputs produced the results. A continued run appends one line.
6. Runs MaxQuant in the container and writes the full log to `results/maxquant_run.log` (appended when continued).
7. Moves `combined/` from `data/` to `results/` (MaxQuant 2.9.0 always writes it into `data/`).

It does not run downstream analysis. QC and statistics belong in `analysis_of_results/`.

### Rules

- One runner script per project, at the project root.
- The script determines its own location, so it works from any working directory.
- No project-specific settings in the script beyond what is listed in its header. Analysis settings belong in `mqpar.xml`.
- Change the script only when the run procedure changes, not for parameter changes.
- Never edit or overwrite a runner script while it runs. Bash reads a script while executing it and, after MaxQuant finishes, continues at the old byte position in the changed file: the final step (moving `combined/` to `results/`) is then lost or replaced by other commands. Before changing a runner, check that no run of that project is active: `pgrep -fl run_<project_name>.sh` and `podman ps --filter name=mq_<project_name>`.

### Usage

```bash
projects/<project_name>/run_<project_name>.sh                     # first run, or continue an unfinished one
projects/<project_name>/run_<project_name>.sh --clean             # rerun from job 1, deletes previous results
projects/<project_name>/run_<project_name>.sh --threads 16        # thread count for this run
projects/<project_name>/run_<project_name>.sh --skip-disk-check   # start despite the disk space check
```

`CLEAN=1` and `SKIP_DISK_CHECK=1` work like `--clean` and `--skip-disk-check`.

Long runs: start in the background and follow the log. On the Mac, `caffeinate -i` keeps it awake while the runner runs:

```bash
P=projects/<project_name>
mkdir -p $P/results
nohup caffeinate -i $P/run_<project_name>.sh > $P/results/nohup.out 2>&1 &
tail -f $P/results/maxquant_run.log
podman ps --filter name=mq_<project_name>                    # still running?
```

`nohup.out` holds the runner's own messages. If a check fails before MaxQuant starts (steps 1 to 4), the error appears only there, and `maxquant_run.log` is not written.

A finished run ends with `Finish writing tables` in the log. On failure, the log names a job. Details are in `combined/proc/<job>.error.txt`: under `data/` if the script stopped before moving `combined/`, otherwise under `results/`.

To stop a run: `podman stop mq_<project_name>`. To continue it later, start the runner again without `--clean`. Closing the Mac's lid pauses a run, which then continues after waking.

### Template

For a new project, copy this file to `projects/<project_name>/run_<project_name>.sh` and make it executable with `chmod +x`. Nothing needs to be edited unless the image tag changes.

```bash
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
```

Notes on the template:

- `set -euo pipefail` stops the script at the first error. A failing MaxQuant run inside the `tee` pipe ends with a message that running the script again continues the run.
- The script is bash and also runs with macOS's bash 3.2. Run it as `./run_<project>.sh` or with `bash`, not by `source` in zsh.
- MaxQuantCmd prints `--version` and `--help` to stderr, hence `2>&1` in the run info.
- Messages name paths relative to the project folder, never complete host paths (section 10).
- Checksums use `sha256sum` (Linux) or `shasum` (macOS). Checksumming large raw files takes a few seconds per GB.
- MaxQuant's `.index` files must be deleted on reruns with `--clean`, since a broken index from a failed run can be reused silently.
- `--dryrun` creates an empty `data/combined/proc`. The `rmdir` line removes it, so a manual dry run, or a runner that stopped before MaxQuant wrote anything, does not count as an unfinished run.
- A continued run skips the dry run and uses `results/dryrun.txt` of the original run for the job numbers; the unchanged checksums guarantee that the job list is the same. It also works when every job had finished and only step 7 was missing.
- A container left over from an interrupted run (for example after a SLURM job was killed) would block the container name. The script removes it, unless it is still running, in which case it stops: two MaxQuant processes on the same folder would corrupt the run.
- The disk space check runs after the `--clean` cleanup, so space freed by deleting previous results counts. The factor 2.2 is the larger of the two values measured in the tutorial projects (1.6 and 2.2 times the raw data). It is an estimate: the true need can be lower, which `--skip-disk-check` allows for, at the risk of a run that fails when the disk is full.

---

## 6. Agent workflow for a new dataset

1. Check that the podman image exists (section 1).
2. Inventory `data/`: raw files, FASTA files, sizes. A raw file of a few KB is a failed download. For a public dataset (for example PRIDE), download it with a script `get_data_<project_name>.sh` that checks each file's size and checksum and records them in a committed `data_manifest.tsv`, so the inputs can be fetched and verified again. Example: `projects/Liljedahl_2017/`. Check disk space first: MaxQuant's working files in `data/` and `combined/` need another 1.6 to 2.2 times the raw data size (measured in the tutorial projects), on top of the raw files.
3. Check every FASTA file:

```bash
grep -c ">" <db>.fasta                                 # entry count
head -2 <db>.fasta                                     # header format
grep -ciE "^>(rev_|decoy|##)" <db>.fasta               # must be 0
```

MaxQuant generates its own decoys (`REV__`) and adds common contaminants (`CON__`). A FASTA with decoys gives a wrong FDR.

4. Look for `parameter/mqpar.xml` and `parameter/design.txt`. Reuse them after the checks in section 3, or prepare new ones.
5. **Ask the user** about everything that cannot be derived from the data (list below). Present the proposed `design.txt` and the key parameters in one message and wait for confirmation.
6. Generate or update the parameter files, apply edits, verify with `--dryrun`.
7. Look for `run_<project_name>.sh`. If it is missing, create it from the template in section 5. For runs on a server, add `run_server.sh`, `setup.sh` and `environment.yml` (section 8).
8. Start the run through the runner script, in the background for long runs, and monitor the log.
9. Report: where results are, protein group and peptide counts from `summary.txt` or `proteinGroups.txt`, anything unusual in the log.
10. Downstream analysis goes into `analysis_of_results/`, reading from `results/combined/txt/`.

### Ask the user about

- **Instrument and acquisition type** (Orbitrap, Astral, timsTOF, DDA, DIA, TMT). Not derivable from file names.
- **Sample names and grouping**: which raw files belong to which sample or condition, whether replicates exist.
- **Fractionation**: whether several raw files form one sample.
- **Database**: whether the FASTA is complete for the organisms in the sample (for example host plus pathogen proteome), and whether additional databases should be added.
- **Modifications**: any deviation from the defaults, PTM enrichment, labels (SILAC, TMT).
- **Match between runs**: useful for several runs of similar samples, misleading for samples of very different composition.
- **Any conflict** between an existing mqpar.xml and the data (different file names, version, paths).
- **Before `--clean`**: deleting existing results, or discarding an unfinished run.

### May decide without asking

- Default search parameters from section 3.
- `numThreads` based on the CPUs podman can use (`podman info --format '{{.Host.CPUs}}'`).
- Experiment names derived from file names, if the user gave no names and the proposal is shown before the run.
- Creating a missing runner script from the template.

---

## 7. Output

`results/` after a run:

```
results/
├── combined/
│   ├── txt/            result tables
│   └── proc/           job logs, *.error.txt on failure
├── maxquant_run.log    full MaxQuant log
├── nohup.out           runner messages of background runs, not in git
├── dryrun.txt          job list
└── run_info.txt        image, version, input checksums
```

Main tables in `combined/txt/`:

| File | Content |
|---|---|
| `proteinGroups.txt` | protein level: identification, intensities, LFQ |
| `peptides.txt` | peptide level |
| `evidence.txt` | every peptide feature per raw file |
| `msms.txt` | every identified MS/MS spectrum |
| `summary.txt` | identification statistics per raw file |
| `parameters.txt` | parameter summary |

### MaxQuant 2.9.0 specifics

Relevant when adapting older tutorials or scripts:

- `proteinGroups.txt` has **no `Reverse` column**. Decoys are marked `+` in the `Decoy` column. Do not filter on `REV__` anywhere in `Protein IDs`: a target group can contain a `REV__` member (in `tutorial_galaxy`, gelsolin: `P06396;CON__Q3SX14;REV__Q6TDU7`).
- `Potential contaminant` is a column, `+` marks contaminants.
- `evidence.txt` protein columns are `Proteins`, `Leading proteins`, `Leading razor proteins` (plural).
- **LFQ intensity is 0 for proteins without at least one peptide quantified in two samples**, since MaxLFQ is built from pairwise peptide ratios. Older versions reported LFQ values for proteins seen in only one sample. For within-sample comparisons, use `Intensity`, or iBAQ. iBAQ is off in files generated with `--create` (`<ibaq>False</ibaq>`) and must be switched on before the run.
- Column positions change with version and sample count. Select columns by header name, case-insensitively (`Protein names` and `Protein Names` both occur).
- In pandas 3, `astype(str)` keeps missing values as `NaN`. Convert with `map(str, ...)` before joining text.
- MaxQuant's contaminant list contains human proteins (albumin, keratins). In human samples these can be genuine sample proteins. Filter only non-human contaminants in that case.

### Quality control

PTXQC runs natively in R on the host, not in the container. Place the call in `analysis_of_results/`:

```bash
Rscript -e "library(PTXQC); createReport('$(pwd)/projects/<project_name>/results/combined/txt')"
```

Use an absolute path. With a relative path, report rendering fails. Install once with `install.packages("PTXQC")`. PTXQC writes its report files into the txt folder it reads. This is the one accepted exception to the rule that analysis does not write into `results/`: it only adds files and leaves MaxQuant's tables unchanged. A rerun with `--clean` deletes the report together with `combined/`. PTXQC thresholds are tuned for cell lysate runs, so low-complexity samples (plasma, serum, purified complexes) fail some metrics without a quality problem.

---

## 8. Server runs

### Purpose

`run_server.sh` turns a fresh clone on a Linux server into a MaxQuant run with one command: it installs the host tools, builds the image, downloads the data and hands over to the runner script. It only prepares; the run itself is the runner's job (section 5). In a batch system such as SLURM, the batch script only needs to call `bash projects/<project_name>/run_server.sh`.

### Files

Three files in the project folder. Example: `projects/Bergant_2023/`.

| File | Content |
|---|---|
| `run_server.sh` | the entry point. Generic: it takes the project name from its folder. |
| `setup.sh` | installs micromamba if neither micromamba nor conda exists (to `~/.local/bin`, environments under `$MAMBA_ROOT_PREFIX`, default `~/micromamba`); creates the conda environment from `environment.yml` and recreates it when the file changed; builds the image if it is missing |
| `environment.yml` | host tools for the download and the checksums: python, curl, perl (`shasum`). MaxQuant needs nothing from it. `name:` must differ between projects, for example `maxquant_<project_name>`. |

For a new project, copy the three files and change `name:` in `environment.yml`. The project also needs `get_data_<project_name>.sh`, `run_<project_name>.sh`, `parameter/mqpar.xml` and `parameter/design.txt`.

### Requirements

- Linux x86_64 with bash, curl, tar and bzip2 (the last three for installing micromamba).
- Rootless podman (tested with 4.9.3). No .NET, Java or MaxQuant installation on the host.
- Internet access: micromamba and conda-forge, `mcr.microsoft.com` for the .NET base image, and the data source (PRIDE).
- `MaxQuant_v.2.9.0/`, copied by hand to the repository root (not in git, section 10).
- Disk space in the project folder for the raw data plus 2.2 times its size.

### Usage

```bash
git clone <repository> && cd maquant_tool                     # then copy MaxQuant_v.2.9.0/ here
bash projects/<project_name>/run_server.sh --dry-run          # checks only: no download, no MaxQuant
bash projects/<project_name>/run_server.sh --download-only    # environment, image and data
bash projects/<project_name>/run_server.sh                    # everything, then the MaxQuant run
bash projects/<project_name>/run_server.sh --threads 32       # runner options are passed through
```

Run `--dry-run` first, once on the login node and once inside a batch job: podman can behave differently on the two.

### What it does

1. Separates its own options (`--dry-run`, `--download-only`) from the runner's, which it passes through.
2. Checks that the project's scripts and parameter files, `Dockerfile`, `.dockerignore` and podman exist, before anything expensive.
3. Runs `setup.sh`.
4. Activates the conda environment by its prefix: it finds the environment in `micromamba env list` or `conda env list`, puts its `bin/` first in `PATH` and sets `CONDA_PREFIX`. `conda activate` needs shell initialisation, which batch jobs lack, and a name alone is ambiguous when a server has several conda installations.
5. Checks the tools (`podman`, `curl`, `python3`, `shasum`) and podman itself: the image must start (`MaxQuantCmd 2.9.0`), and a file written in the container into the project folder must belong to the user.
6. Downloads the data with `get_data_<project_name>.sh`, which skips verified files. `--dry-run` stops before this step, `--download-only` after it.
7. Hands over to `run_<project_name>.sh` with `exec`: the runner replaces the process and receives signals directly, for example from SLURM at the time limit.

Every message has the form `[HH:MM:SS] step` or `[HH:MM:SS] ERROR: ...`, so a hanging run can be located in the batch log. Every step skips what exists, so a rerun repeats only what is missing.

### Batch jobs (SLURM)

- MaxQuant uses CPUs only, no GPU.
- Threads: the runner uses `SLURM_CPUS_PER_TASK` unless `--threads` is given. Memory: plan about 3 GB per thread (a run with 8 threads peaked at 20 GB).
- Time limit: when a job is killed, submit the same batch script again. The runner continues at the interrupted job (section 5). `--clean` starts from job 1 instead.
- If compute nodes have no internet, run `--download-only` on the login node first.
- Podman stores images per user, in the folder that `podman info --format '{{.Store.GraphRoot}}'` prints. If it differs between the login node and the compute nodes, each node needs its own image, which `setup.sh` then builds there.

### Checking podman on a server

```bash
podman --version
podman info --format 'rootless={{.Host.Security.Rootless}} storage={{.Store.GraphRoot}}'
grep "^$USER:" /etc/subuid /etc/subgid
```

| Check | Expected | Otherwise |
|---|---|---|
| `rootless=` | `true` | podman runs as root; files written by MaxQuant would belong to root |
| `/etc/subuid`, `/etc/subgid` | both list the user | rootless podman needs these entries; ask the administrators |
| `storage=` | the same folder on login and compute nodes | each node builds its own image (see above) |

`run_server.sh --dry-run` then tests the rest: the image starts, and podman writes files into the project folder that belong to the user.

---

## 9. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `Could not load file or assembly 'ThermoFisher.CommonCore.RawFileReader'` | MaxQuant was run natively on macOS. Use the container. |
| `... does not exist` for a path | Host path instead of `/work/...`, or the file is outside the mounted project folder. |
| `no .raw files in data/` from the runner, although raw files are there | Upper-case extension `.RAW`, as in many PRIDE datasets. The runner's `data/*.raw` pattern is case-sensitive. Save or rename the files as `.raw`. |
| `a finished run exists in results/combined` from the runner | Results from an earlier run. Confirm with the user, then rerun with `--clean`. |
| `data/combined holds an unfinished run, but ... differs from results/run_info.txt` | `mqpar.xml` or `design.txt` changed after the unfinished run started, so it cannot be continued. Confirm with the user, then rerun with `--clean`. |
| `MaxQuant stopped before the end` | MaxQuant was stopped (`podman stop`, a SLURM time limit) or failed. After a stop, run the script again: it continues. After a failure, read the log and `data/combined/proc/<job>.error.txt` first; continuing repeats the failed job. |
| `container mq_<project_name> is running already` | Another run of this project is active, possibly from an earlier batch job. Check `podman ps --filter name=mq_<project_name>`; stop it with `podman stop mq_<project_name>` if it is not wanted. |
| `not enough disk space` from the runner | Free space on the disk that holds the project, for example MaxQuant's working files in `data/` of finished projects (their results stay; a rerun with `--clean` recreates them). Ask the user before deleting anything. `--skip-disk-check` starts anyway. |
| Deleting files frees less disk space than expected (Mac) | The Linux VM that runs the containers can keep files open that a container has used, even after they are deleted on the Mac (`lsof +L1 \| grep <project_name>` lists them). Restart the VM: `podman machine stop && podman machine start`. This can also hide the space freed by `--clean` from the runner's disk space check. |
| `the MaxQuant image does not start` from `run_server.sh` | Run `podman run --rm --pull never --platform linux/amd64 localhost/maxquant:2.9.0 --version` by hand to see podman's error. A broken build: `podman rmi localhost/maxquant:2.9.0`, then `setup.sh` builds it again. |
| `files written in the container do not belong to you` | Podman does not run rootless, or the user has no entries in `/etc/subuid` and `/etc/subgid` (section 8). |
| `podman cannot write to the mounted project folder` | Missing write permission on the project folder, or a security module other than SELinux blocks the mount. Run the failing `podman run` command by hand to see the error. |
| Podman errors in batch jobs about `XDG_RUNTIME_DIR` or `/run/user` | Rootless podman needs a runtime folder, which batch jobs may lack. Ask the administrators for the recommended setting; `--dry-run` inside a batch job shows whether it works. |
| `zsh: no matches found` | zsh aborts a command when a wildcard matches nothing. Name files explicitly, or use the bash runner script. |
| A job fails | Read `combined/proc/<job>.error.txt`. |
| Run is slow (Mac) | x64 emulation in the podman machine costs speed. Check the machine's CPUs and memory with `podman machine inspect`. On a Linux x86_64 server, the image runs without emulation. |
| Multi-line command breaks | A `\` must be the last character of the line. |

## 10. Repository notes

- **No complete host paths in any committed file**: no `/Users/<name>/...` or other home folder paths in scripts, logs, docs or parameter files. Scripts print paths relative to the project folder. The `/work/...` container paths in mqpar.xml are required by MaxQuant and reveal nothing about the host. Check before every commit, the command must print nothing:

```bash
git grep -n -I -E "/Users/[A-Za-z]|/home/[a-z]" ; git ls-files -o --exclude-standard | xargs grep -n -I -E "/Users/[A-Za-z]|/home/[a-z]"
```

- `MaxQuant_v.2.9.0/` contains MaxQuant's binaries and Thermo's RawFileReader DLLs. Both are distributed under their own licenses, which restrict redistribution. Keep this folder out of git.
- Raw files and MaxQuant output are large. Keep them out of git. Keep `run_info.txt`, `dryrun.txt` and the log, since they document the run. Commit them only after the run has finished. `nohup.out` is ignored: it duplicates the log plus a few runner messages.
- Papers (PDF) may be kept in a project folder for reference, but stay out of git.

`.gitignore` of this repository:

```
.DS_Store
/MaxQuant_v.2.9.0/
/projects/**/data/
/projects/**/results/combined/
/projects/**/results/nohup.out
/projects/**/parameter/mqpar_default.xml
/projects/**/analysis_of_results/figures/
/projects/**/analysis_of_results/msstats/FeatureLevelData.tsv
/projects/*/*.pdf
```
