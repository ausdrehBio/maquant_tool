#!/usr/bin/env Rscript
# MSstats analysis for tutorial_galaxy_2, following the Galaxy tutorial
# "MaxQuant and MSstats for the analysis of label-free data" (MSstats 4.0.0 there).
#
# Usage (from any directory):
#   Rscript run_msstats.R                      MaxQuant output of this project
#   Rscript run_msstats.R EVIDENCE PROTEINGROUPS OUT_DIR   other MaxQuant output (e.g. a test)
#
# Reads (read only): results/combined/txt/{evidence,proteinGroups}.txt,
#   parameter/input_annotation_file.tabular, parameter/input_comparison_matrix.tabular
# Writes: analysis_of_results/msstats/ (tables, MSstats log, QC and volcano plots)
#
# Steps, as in the tutorial:
#   1. keep only lines with a human protein ("HUMAN" in the protein IDs); this removes
#      non-human contaminants and keeps human ones such as albumin and keratins
#   2. clear the "Potential contaminant" flag, so MSstats keeps those human proteins
#   3. MaxQtoMSstatsFormat with Leading razor protein as protein ID and proteins with
#      only one peptide removed
#   4. dataProcess with MSstats defaults (log2, equalizeMedians, TMP summarisation,
#      censored values imputed); featureSubset = "all", the default of MSstats 4.0.0
#   5. groupComparison with the tutorial's comparison matrix (metast vs RDEB)

suppressPackageStartupMessages({
  library(MSstats)
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
project_dir <- normalizePath(file.path(dirname(normalizePath(script_file)), "..", ".."))
rel <- function(path) {
  # path relative to the project folder for printing; never a complete host path
  p <- normalizePath(path, mustWork = FALSE)
  if (startsWith(p, paste0(project_dir, "/"))) substring(p, nchar(project_dir) + 2) else basename(p)
}

args <- commandArgs(trailingOnly = TRUE)
txt_dir <- file.path(project_dir, "results", "combined", "txt")
evidence_path <- if (length(args) >= 1) args[1] else file.path(txt_dir, "evidence.txt")
pg_path <- if (length(args) >= 2) args[2] else file.path(txt_dir, "proteinGroups.txt")
out_dir <- if (length(args) >= 3) args[3] else file.path(project_dir, "analysis_of_results", "msstats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out_dir, "MSstats_log.txt")
if (file.exists(log_file)) file.remove(log_file)

read_mq <- function(path) {
  # check.names gives the read.table-style names MSstats expects ("Raw.file", ...);
  # integer64 = "double": intensities above 2^31 would otherwise become integer64,
  # which MSstats' aggregation cannot handle
  fread(path, sep = "\t", quote = "", data.table = FALSE, check.names = TRUE,
        integer64 = "double", showProgress = FALSE)
}

harmonise <- function(df) {
  # MaxQuant 2.9.0 marks decoys in "Decoy"; MSstats looks for "Reverse"
  if (!"Reverse" %in% names(df) && "Decoy" %in% names(df)) df$Reverse <- df$Decoy
  # MaxQuant 2.9.0 names the evidence column "Leading razor proteins" (plural)
  if (!"Leading.razor.protein" %in% names(df) && "Leading.razor.proteins" %in% names(df)) {
    df$Leading.razor.protein <- df$Leading.razor.proteins
  }
  df
}

keep_human <- function(df, id_col) {
  # tutorial: Select lines matching "HUMAN", then replace "+" in Potential contaminant
  keep <- grepl("HUMAN", df[[id_col]], fixed = TRUE)
  out <- df[keep, , drop = FALSE]
  out$Potential.contaminant <- ""
  list(data = out, removed = sum(!keep), before = nrow(df))
}

cat("reading MaxQuant output:", rel(evidence_path), "and", rel(pg_path), "\n")
ev_raw <- harmonise(read_mq(evidence_path))
pg_raw <- harmonise(read_mq(pg_path))
pg_h <- keep_human(pg_raw, "Protein.IDs")
ev_h <- keep_human(ev_raw, "Proteins")
cat(sprintf("protein groups: %d, without a human protein (removed): %d\n", pg_h$before, pg_h$removed))
cat(sprintf("evidence rows:  %d, without a human protein (removed): %d\n", ev_h$before, ev_h$removed))

# annotation: match the tutorial's Galaxy file names (metast_cSCC1.raw.thermo) to the
# raw file names in evidence.txt (metast_cSCC1 here; Experiment105_metast_cSCC1.raw.thermo
# in the tutorial's own MaxQuant output)
core_name <- function(x) sub("^Experiment[0-9]+_", "", sub("(\\.raw)?(\\.thermo)?(\\.raw)?$", "", x))
annot <- fread(file.path(project_dir, "parameter", "input_annotation_file.tabular"), data.table = FALSE)
names(annot)[tolower(names(annot)) == "run"] <- "Run"
raw_files <- unique(ev_h$data$Raw.file)
idx <- match(core_name(annot$Raw.file), core_name(raw_files))
if (anyNA(idx) || anyDuplicated(idx) || length(raw_files) != nrow(annot)) {
  stop("annotation and evidence raw files do not match one to one:\n  annotation: ",
       paste(core_name(annot$Raw.file), collapse = ", "), "\n  evidence:   ",
       paste(sort(core_name(raw_files)), collapse = ", "))
}
annot$Raw.file <- raw_files[idx]
fwrite(annot, file.path(out_dir, "annotation_used.tsv"), sep = "\t")
cat(sprintf("annotation: %d runs, conditions: %s\n", nrow(annot),
            paste(names(table(annot$Condition)), table(annot$Condition), sep = " x", collapse = ", ")))

quant <- MaxQtoMSstatsFormat(
  evidence = ev_h$data, annotation = annot, proteinGroups = pg_h$data,
  proteinID = "Leading.razor.protein",
  removeProtein_with1Peptide = TRUE,
  use_log_file = TRUE, append = FALSE, verbose = FALSE, log_file_path = log_file
)

processed <- dataProcess(
  quant, featureSubset = "all",
  use_log_file = TRUE, append = TRUE, verbose = FALSE, log_file_path = log_file
)
fwrite(processed$FeatureLevelData, file.path(out_dir, "FeatureLevelData.tsv"), sep = "\t")
fwrite(processed$ProteinLevelData, file.path(out_dir, "ProteinLevelData.tsv"), sep = "\t")

sample_q <- quantification(processed, type = "Sample", format = "matrix",
                           use_log_file = TRUE, append = TRUE, verbose = FALSE, log_file_path = log_file)
group_q <- quantification(processed, type = "Group", format = "matrix",
                          use_log_file = TRUE, append = TRUE, verbose = FALSE, log_file_path = log_file)
fwrite(sample_q, file.path(out_dir, "SampleQuantificationMatrix.tsv"), sep = "\t")
fwrite(group_q, file.path(out_dir, "GroupQuantificationMatrix.tsv"), sep = "\t")

dataProcessPlots(processed, type = "QCPlot", which.Protein = "allonly", text.angle = 0,
                 address = file.path(out_dir, ""))

cmp <- fread(file.path(project_dir, "parameter", "input_comparison_matrix.tabular"), data.table = FALSE)
contrast <- as.matrix(cmp[, -1, drop = FALSE])
rownames(contrast) <- cmp[[1]]
grp <- processed$ProteinLevelData$GROUP
groups <- if (is.factor(grp)) levels(grp) else sort(unique(as.character(grp)))
if (!setequal(colnames(contrast), groups)) {
  stop("comparison matrix columns (", paste(colnames(contrast), collapse = ", "),
       ") do not match the conditions (", paste(groups, collapse = ", "), ")")
}
contrast <- contrast[, groups, drop = FALSE]  # MSstats expects the condition order

comparison <- groupComparison(
  contrast.matrix = contrast, data = processed, save_fitted_models = FALSE,
  use_log_file = TRUE, append = TRUE, verbose = FALSE, log_file_path = log_file
)
fwrite(comparison$ComparisonResult, file.path(out_dir, "ComparisonResult.tsv"), sep = "\t")
groupComparisonPlots(comparison$ComparisonResult, type = "VolcanoPlot", FCcutoff = 1.5,
                     logBase.pvalue = 2, ProteinName = FALSE, address = file.path(out_dir, ""))

writeLines(c(
  sprintf("R\t%s.%s", R.version$major, R.version$minor),
  sprintf("MSstats\t%s", packageVersion("MSstats")),
  sprintf("MSstatsConvert\t%s", packageVersion("MSstatsConvert")),
  sprintf("protein_groups_total\t%d", pg_h$before),
  sprintf("protein_groups_without_human_protein\t%d", pg_h$removed),
  sprintf("evidence_rows_total\t%d", ev_h$before),
  sprintf("evidence_rows_without_human_protein\t%d", ev_h$removed)
), file.path(out_dir, "run_info.txt"))

n_prot <- length(unique(processed$ProteinLevelData$Protein))
cat(sprintf("proteins after dataProcess: %d | compared: %d\n", n_prot,
            length(unique(comparison$ComparisonResult$Protein))))
cat("written to:", rel(out_dir), "\n")
