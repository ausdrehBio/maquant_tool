#!/usr/bin/env Rscript
# MSstats analysis for Liljedahl_2017: four groups of mice, 5 mice per group, two technical
# replicates (runs) per mouse.
#
# Usage (from any directory):
#   Rscript run_msstats.R                                  MaxQuant output of this project
#   Rscript run_msstats.R EVIDENCE PROTEINGROUPS OUT_DIR   other MaxQuant output (e.g. a test)
#
# Reads (read only): results/combined/txt/{evidence,proteinGroups}.txt,
#   parameter/sample_annotation.tsv, parameter/comparison_matrix.tsv
# Writes: analysis_of_results/msstats/ (tables, MSstats log, QC and volcano plots)
#
# Steps:
#   1. keep only rows with a target protein (a member that is neither CON__ nor REV__; the
#      FASTA is mouse only); this removes contaminant-only groups and keeps contaminant groups
#      with a mouse protein, such as albumin, whose contaminant flag is cleared so that
#      MSstats keeps them
#   2. MaxQtoMSstatsFormat with Leading razor protein as protein ID and proteins with only
#      one peptide removed
#   3. dataProcess with MSstats defaults (log2, equalizeMedians, TMP summarisation, censored
#      values imputed, featureSubset = "all")
#   4. groupComparison with the contrasts of parameter/comparison_matrix.tsv; BioReplicate is
#      the mouse, so the two runs of a mouse are treated as technical replicates

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

has_target <- function(ids) {
  # TRUE where a ";"-separated ID list has a member that is neither CON__ nor REV__
  vapply(strsplit(as.character(ids), ";", fixed = TRUE),
         function(x) any(nzchar(x) & !startsWith(x, "CON__") & !startsWith(x, "REV__")),
         logical(1))
}

keep_target <- function(df, id_col) {
  keep <- has_target(df[[id_col]])
  out <- df[keep, , drop = FALSE]
  out$Potential.contaminant <- ""
  list(data = out, removed = sum(!keep), before = nrow(df))
}

cat("reading MaxQuant output:", rel(evidence_path), "and", rel(pg_path), "\n")
ev_raw <- harmonise(read_mq(evidence_path))
pg_raw <- harmonise(read_mq(pg_path))
pg_t <- keep_target(pg_raw, "Protein.IDs")
ev_t <- keep_target(ev_raw, "Proteins")
cat(sprintf("protein groups: %d, without a target protein (removed): %d\n", pg_t$before, pg_t$removed))
cat(sprintf("evidence rows:  %d, without a target protein (removed): %d\n", ev_t$before, ev_t$removed))

# annotation: only the columns MSstats uses; MSstats renames Raw.file to Run itself
annot_all <- fread(file.path(project_dir, "parameter", "sample_annotation.tsv"), data.table = FALSE)
annot <- annot_all[, c("Raw.file", "Condition", "BioReplicate", "IsotopeLabelType")]
raw_files <- unique(ev_t$data$Raw.file)
if (!setequal(annot$Raw.file, raw_files) || anyDuplicated(annot$Raw.file)) {
  stop("annotation and evidence raw files do not match one to one:\n  only in annotation: ",
       paste(setdiff(annot$Raw.file, raw_files), collapse = ", "), "\n  only in evidence:   ",
       paste(setdiff(raw_files, annot$Raw.file), collapse = ", "))
}
fwrite(annot, file.path(out_dir, "annotation_used.tsv"), sep = "\t")
cat(sprintf("annotation: %d runs, %d mice, conditions: %s\n", nrow(annot),
            length(unique(annot$BioReplicate)),
            paste(names(table(annot$Condition)), table(annot$Condition), sep = " x", collapse = ", ")))

quant <- MaxQtoMSstatsFormat(
  evidence = ev_t$data, annotation = annot, proteinGroups = pg_t$data,
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

dataProcessPlots(processed, type = "QCPlot", which.Protein = "allonly", text.angle = 90,
                 address = file.path(out_dir, ""))

cmp <- fread(file.path(project_dir, "parameter", "comparison_matrix.tsv"), data.table = FALSE)
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
  sprintf("protein_groups_total\t%d", pg_t$before),
  sprintf("protein_groups_without_target_protein\t%d", pg_t$removed),
  sprintf("evidence_rows_total\t%d", ev_t$before),
  sprintf("evidence_rows_without_target_protein\t%d", ev_t$removed)
), file.path(out_dir, "run_info.txt"))

n_prot <- length(unique(processed$ProteinLevelData$Protein))
cat(sprintf("proteins after dataProcess: %d | compared: %d | contrasts: %s\n", n_prot,
            length(unique(comparison$ComparisonResult$Protein)),
            paste(unique(as.character(comparison$ComparisonResult$Label)), collapse = ", ")))
cat("written to:", rel(out_dir), "\n")
