# MaxQuant in this repository: workflow reference

This file describes how MaxQuant is run in this repository. Read it fully before processing a dataset.

## Overview

- MaxQuant 2.9.0 runs inside a local Docker image, `maxquant:2.9.0`.
- Each dataset lives in its own folder under `projects/`.
- A MaxQuant run needs two input files besides the data: `mqpar.xml` (all analysis parameters) and `design.txt` (sample assignment of the raw files). Both live in the project's `parameter/` folder.
- Every project has a runner script, `run_<project_name>.sh`, which performs the MaxQuant run. **Full runs are only started through this script**, so every result in the repository is reproducible with one command.
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
        ├── data/                  .raw files and FASTA file(s), input only
        ├── parameter/             mqpar.xml, design.txt
        ├── results/               MaxQuant output, run log, run info
        └── analysis_of_results/   downstream scripts, plots, tables
```

Roles of the folders:

| Folder | Content | Written by |
|---|---|---|
| `data/` | raw files, FASTA | the user. Never modified, except MaxQuant's per-file working folders |
| `parameter/` | `mqpar.xml`, `design.txt` | generated once, then edited only on purpose |
| `results/` | `combined/`, `maxquant_run.log`, `run_info.txt` | only the runner script |
| `analysis_of_results/` | QC, plots, statistics | analysis scripts, which read `results/` and never modify it |

Reference project: `projects/tutorial_galaxy/` (Galaxy Training Network tutorial "Label-free data analysis using MaxQuant", two human serum samples).

---

## 1. The Docker image

### Why Docker

MaxQuant's GUI runs only on Windows. The command-line version runs on .NET 8, but on native macOS with Apple Silicon, Thermo's raw file reader (`ThermoFisher.CommonCore.RawFileReader.dll`) fails to load. Do not attempt native runs on the Mac. Inside a Linux x64 container, which is MaxQuant's officially supported platform, it works. Docker Desktop runs the container through Rosetta x64 emulation.

### Check that the image exists

```bash
docker image inspect maxquant:2.9.0 > /dev/null && echo "image present"
```

### Rebuild only if it is missing

Requires the unzipped release in `MaxQuant_v.2.9.0/` at the repository root. If `Dockerfile` or `.dockerignore` are missing at the root, create them with this content.

`Dockerfile`:

```dockerfile
FROM --platform=linux/amd64 mcr.microsoft.com/dotnet/runtime:8.0
COPY MaxQuant_v.2.9.0 /opt/maxquant
WORKDIR /work
ENTRYPOINT ["dotnet", "/opt/maxquant/bin/MaxQuantCmd.dll"]
```

`.dockerignore` (keeps raw data out of the build):

```
projects/
.git/
```

```bash
docker build --platform linux/amd64 -t maxquant:2.9.0 .
docker run --rm --platform linux/amd64 maxquant:2.9.0 --version   # expect: MaxQuantCmd 2.9.0.0
```

The linter warning `FromPlatformFlagConstDisallowed` is harmless.

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

For generating, editing and checking parameter files (sections 3 and 4). Run from the repository root:

```bash
P=projects/<project_name>
docker run --rm --platform linux/amd64 -v "$(pwd)/$P":/work maxquant:2.9.0 <MaxQuant arguments>
```

`ENTRYPOINT` already calls `MaxQuantCmd.dll`, so only MaxQuant's arguments follow the image name. All MaxQuant options:

```bash
docker run --rm --platform linux/amd64 maxquant:2.9.0 --help
```

Full runs go through the runner script (section 5), not through ad hoc `docker run` calls.

---

## 3. mqpar.xml

### Purpose

`mqpar.xml` is MaxQuant's parameter file. It holds every setting of an analysis in one XML file: input files, FASTA databases, sample assignment, enzyme, modifications, mass tolerances, FDR thresholds, quantification method and output options. The command-line version takes its settings only from this file. The Windows GUI and the Galaxy MaxQuant tool also write an mqpar.xml before each run, so it is the common exchange format. MaxQuant copies it into its output, which makes it the record of how the data were processed.

### Properties

- Plain XML. Each tag is one parameter, for example `<maxMissedCleavages>2</maxMissedCleavages>`. Lists use `<string>` children.
- Contains absolute paths (FASTA, raw files, output folder). In this repository these must be `/work/...` paths.
- **Version-bound.** An mqpar.xml runs only with the MaxQuant version that created it. Check `<maxQuantVersion>`.
- Parameters inside `<parameterGroups>` apply per group of raw files. Normally all files share one group.
- The same modification names appear in several blocks (searched modifications, modifications used for quantification). Edit only the block you intend to change.

### Generate a new one

`design.txt` must exist first (section 4).

```bash
docker run --rm --platform linux/amd64 -v "$(pwd)/$P":/work maxquant:2.9.0 \
  --create \
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
| `--numThreads` | CPU threads, at or below `docker info --format '{{.NCPU}}'` | |

Only standard DDA on Thermo Orbitrap data has been tested with this image.

### Output location

By default, MaxQuant writes `combined/` into the raw file folder (`data/`). To direct it to `results/`, check whether the file has a fixed combined folder parameter and set it:

```bash
grep -n "fixedCombinedFolder" $P/parameter/mqpar.xml
sed -i '' 's|<fixedCombinedFolder>.*</fixedCombinedFolder>|<fixedCombinedFolder>/work/results</fixedCombinedFolder>|' $P/parameter/mqpar.xml
```

If the tag does not exist, MaxQuant writes to `data/combined/`. The runner script then moves it to `results/` after the run (section 5). After the first run of a project, verify where `combined/` was actually written. Per-raw-file working folders appear next to the raw files in `data/`.

### Default search parameters

These defaults suit standard bottom-up experiments. Keep them unless the user specifies otherwise.

| Parameter | Tag | Default |
|---|---|---|
| Enzyme | `<enzymes>` | Trypsin/P (cleaves after K/R, also before P) |
| Max. missed cleavages | `<maxMissedCleavages>` | 2 |
| Fixed modifications | `<fixedModifications>` | Carbamidomethyl (C) |
| Variable modifications | `<variableModifications>` | Oxidation (M), Acetyl (Protein N-term) |
| Min. peptide length | `<minPepLen>` | 7 |
| Min. unique peptides | `<minUniquePeptides>` | 0 |
| PSM and protein FDR | | 0.01 |
| Identifier parse rule | `<identifierParseRule>` | `>[^|]*\|(.*?)\|` (UniProt headers) |
| Description parse rule | `<descriptionParseRule>` | `>(.*)` |

If FASTA headers are not in UniProt format (`>sp|P02768|ALBU_HUMAN ...`), the identifier parse rule must be adapted, otherwise protein IDs are wrong.

The `tutorial_galaxy` project deviates from these defaults on purpose to reproduce the tutorial (1 missed cleavage, no variable modifications, min. 1 unique peptide, description rule `>(.*) OS`). Do not copy these settings to other projects.

### Editing

Single values:

```bash
docker run --rm --platform linux/amd64 -v "$(pwd)/$P":/work maxquant:2.9.0 \
  --changeParameter /work/parameter/mqpar.xml numThreads 8
```

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

1. Check `<maxQuantVersion>`. It must be 2.9.0.x for this image. Otherwise generate a new file and transfer the non-default settings.
2. Check that all paths are `/work/...` paths and point to existing files.
3. Check that `<filePaths>` lists exactly the raw files of this dataset and `<experiments>` the intended sample names, in the same order.
4. To update FASTA and raw paths in one step (arguments: new file, FASTA folder, raw folder):

```bash
docker run --rm --platform linux/amd64 -v "$(pwd)/$P":/work maxquant:2.9.0 \
  /work/parameter/mqpar_template.xml --changeFolder /work/parameter/mqpar.xml /work/data /work/data
```

If the number or names of raw files changed, regenerating the file with a new `design.txt` is less error-prone than editing `<filePaths>` and `<experiments>` by hand.

### Verify

```bash
grep -n -A3 -E "<maxQuantVersion>|<fastaFilePath>|<filePaths>|<experiments>|<fixedCombinedFolder>" $P/parameter/mqpar.xml
docker run --rm --platform linux/amd64 -v "$(pwd)/$P":/work maxquant:2.9.0 /work/parameter/mqpar.xml --dryrun
```

`--dryrun` parses the file and prints the job list without processing. A broken file or path fails here within seconds. The runner script repeats this check before every run.

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

1. Checks prerequisites: Docker available, image present, `mqpar.xml` and `design.txt` present, no host paths in `mqpar.xml`.
2. Refuses to overwrite existing results unless called with `CLEAN=1`. With `CLEAN=1`, it deletes `results/combined/` and MaxQuant's working files in `data/` (per-raw-file folders, `.index` files, MaxQuant's mqpar copy). It never touches raw or FASTA files.
3. Runs `--dryrun` as a fast parameter check.
4. Writes `results/run_info.txt`: date, image ID, MaxQuant version, SHA-256 checksums of `mqpar.xml`, `design.txt`, raw files and FASTA files. The checksums prove which exact inputs produced the results.
5. Runs MaxQuant in the container and writes the full log to `results/maxquant_run.log`.
6. Moves `combined/` from `data/` to `results/` if MaxQuant wrote it there.

It does not run downstream analysis. QC and statistics belong in `analysis_of_results/`.

### Rules

- One runner script per project, at the project root.
- The script determines its own location, so it works from any working directory.
- No project-specific settings in the script beyond what is listed in its header. Analysis settings belong in `mqpar.xml`.
- Change the script only when the run procedure changes, not for parameter changes.

### Usage

```bash
projects/<project_name>/run_<project_name>.sh                # first run
CLEAN=1 projects/<project_name>/run_<project_name>.sh        # rerun, deletes previous results
```

Long runs: start in the background and follow the log:

```bash
nohup projects/<project_name>/run_<project_name>.sh > /dev/null 2>&1 &
tail -f projects/<project_name>/results/maxquant_run.log
docker ps --filter name=mq_<project_name>                    # still running?
```

A finished run ends with `Finish writing tables` in the log. On failure, the log names a job. Details are in `results/combined/proc/<job>.error.txt` (or `data/combined/proc/`). The Mac must stay awake during a run, since sleep pauses the container.

### Template

For a new project, copy this file to `projects/<project_name>/run_<project_name>.sh` and make it executable with `chmod +x`. Nothing needs to be edited unless the image tag changes.

```bash
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
! grep -q "/Users/" "$MQPAR_HOST" || fail "mqpar.xml contains host paths, use /work/..."
compgen -G "data/*.raw" > /dev/null || fail "no .raw files in data/"

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
  echo "maxquant:      $(mq --version | head -1)"
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
```

Notes on the template:

- `set -euo pipefail` stops the script at the first error, including a failing MaxQuant run inside the `tee` pipe.
- The script is bash. Run it as `./run_<project>.sh` or with `bash`, not by `source` in zsh.
- Checksumming large raw files takes a few seconds per GB.
- MaxQuant's `.index` files must be deleted on reruns, since a broken index from a failed run can be reused silently.

---

## 6. Agent workflow for a new dataset

1. Check that the Docker image exists (section 1).
2. Inventory `data/`: raw files, FASTA files, sizes. A raw file of a few KB is a failed download.
3. Check every FASTA file:

```bash
grep -c ">" <db>.fasta                                 # entry count
head -2 <db>.fasta                                     # header format
grep -ciE "^>(rev_|decoy|##)" <db>.fasta               # must be 0
```

MaxQuant generates its own decoys (`REV__`) and adds common contaminants (`CON__`). A FASTA with decoys gives a wrong FDR.

4. Look for `parameter/mqpar.xml` and `parameter/design.txt`. Reuse them after the checks in section 3, or prepare new ones.
5. **Ask the user** about everything that cannot be derived from the data (list below). Present the proposed `design.txt` and the key parameters in one message and wait for confirmation.
6. Generate or update the parameter files, apply edits, set the output folder to `/work/results`, verify with `--dryrun`.
7. Look for `run_<project_name>.sh`. If it is missing, create it from the template in section 5.
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
- **Before `CLEAN=1`**: deleting existing results.

### May decide without asking

- Default search parameters from section 3.
- `numThreads` based on the Docker CPU limit.
- Output location in `results/`.
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

- `proteinGroups.txt` has **no `Reverse` column**. Identify decoys by the `REV__` prefix in `Protein IDs`.
- `Potential contaminant` is a column, `+` marks contaminants.
- `evidence.txt` protein columns are `Proteins`, `Leading proteins`, `Leading razor proteins` (plural).
- **LFQ intensity is 0 for proteins without at least one peptide quantified in two samples**, since MaxLFQ is built from pairwise peptide ratios. Older versions reported LFQ values for proteins seen in only one sample. For within-sample comparisons, use `Intensity` or iBAQ.
- Column positions change with version and sample count. Select columns by header name, case-insensitively (`Protein names` and `Protein Names` both occur).
- In pandas 3, `astype(str)` keeps missing values as `NaN`. Convert with `map(str, ...)` before joining text.
- MaxQuant's contaminant list contains human proteins (albumin, keratins). In human samples these can be genuine sample proteins. Filter only non-human contaminants in that case.

### Quality control

PTXQC runs natively in R on the host, not in the container. Place the call in `analysis_of_results/`:

```bash
Rscript -e "library(PTXQC); createReport('$(pwd)/projects/<project_name>/results/combined/txt')"
```

Use an absolute path. With a relative path, report rendering fails. Install once with `install.packages("PTXQC")`. PTXQC writes its report into the txt folder. PTXQC thresholds are tuned for cell lysate runs, so low-complexity samples (plasma, serum, purified complexes) fail some metrics without a quality problem.

---

## 8. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `Could not load file or assembly 'ThermoFisher.CommonCore.RawFileReader'` | MaxQuant was run natively on macOS. Use the container. |
| `... does not exist` for a path | Host path instead of `/work/...`, or the file is outside the mounted project folder. |
| `previous results exist` from the runner | Results from an earlier run. Confirm with the user, then rerun with `CLEAN=1`. |
| `zsh: no matches found` | zsh aborts a command when a wildcard matches nothing. Name files explicitly, or use the bash runner script. |
| A job fails | Read `combined/proc/<job>.error.txt`. |
| Run is slow | x64 emulation costs speed. Check CPU and memory in Docker Desktop → Settings → Resources. |
| Container name already in use | A previous run is still running or was not removed. Check `docker ps -a --filter name=mq_<project_name>`. |
| Multi-line command breaks | A `\` must be the last character of the line. |

## 9. Repository notes

- `MaxQuant_v.2.9.0/` contains MaxQuant's binaries and Thermo's RawFileReader DLLs. Both are distributed under their own licenses, which restrict redistribution. Keep this folder out of git.
- Raw files and MaxQuant output are large. Keep them out of git. Keep `run_info.txt`, `dryrun.txt` and the log, since they document the run.

Suggested `.gitignore`:

```
MaxQuant_v*/
projects/*/data/
projects/*/results/combined/
projects/*/parameter/mqpar_default.xml
```
