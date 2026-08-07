#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[i + 1L]]
}
root <- normalizePath(value_after("--project-root", "."), winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(value_after("--output-dir", file.path(root, "figures", "Figure_05", "output")), winslash = "/", mustWork = TRUE)
src_dir <- normalizePath(value_after("--source-dir", file.path(root, "figures", "Figure_05", "data")), winslash = "/", mustWork = TRUE)

read_src <- function(name) fread(file.path(src_dir, name), sep = "\t", na.strings = c("", "NA"))
atlas <- read_src("Figure_5_cell_context_atlas.tsv")
exact <- read_src("Figure_5_exact_event_cross_context.tsv")
snuc <- read_src("Figure_5_single_nucleus_eqtl.tsv")
disease <- read_src("Figure_5_disease_state_rna.tsv")
structure <- read_src("Figure_5_ec2_structure.tsv")
events <- read_src("Figure_5_transcript_events.tsv")

stem <- file.path(out_dir, "Figure_5_cell_context_structure_v9")
exports <- paste0(stem, c(".png", ".pdf", ".svg", ".tiff"))
svg_text <- paste(readLines(paste0(stem, ".svg"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")

checks <- data.table(
  check = c(
    "atlas contains four evidence classes", "atlas separates four biological contexts",
    "one lead variant across five exact-event contexts", "single-nucleus eQTL has four estimates",
    "single-nucleus eQTL effects are positive", "disease-state display has four estimates",
    "all displayed disease-state results are FDR non-significant",
    "SEA-AD intervals cross the null", "EC2 residue track is complete",
    "exact coding boundary is AA150/151", "boundary pLDDT is high confidence",
    "png export exists", "pdf export exists", "svg export exists", "tiff export exists",
    "uppercase panel labels", "no local paths in SVG"
  ),
  passed = c(
    uniqueN(atlas$evidence_class) == 4L,
    uniqueN(atlas$context) == 4L,
    uniqueN(exact$snp) == 1L && nrow(exact) == 5L,
    nrow(snuc) == 4L,
    all(snuc$beta > 0),
    nrow(disease) == 4L,
    all(disease$fdr > 0.05),
    all(disease[interval_available == TRUE, lo < 0 & hi > 0]),
    nrow(structure[residue >= 114 & residue <= 232]) == 119L,
    events[event == "project_exon5_exon6", protein_transition] == "AA150/151",
    mean(structure[residue %in% c(150, 151), pLDDT]) > 90,
    file.exists(exports[1]) && file.info(exports[1])$size > 10000,
    file.exists(exports[2]) && file.info(exports[2])$size > 10000,
    file.exists(exports[3]) && file.info(exports[3])$size > 10000,
    file.exists(exports[4]) && file.info(exports[4])$size > 10000,
    all(vapply(paste0(">", LETTERS[1:5], "<"), grepl, logical(1), x = svg_text, fixed = TRUE)),
    !grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", svg_text, perl = TRUE)
  ),
  detail = c(
    uniqueN(atlas$evidence_class), uniqueN(atlas$context), paste(unique(exact$snp), nrow(exact), sep = "; "),
    nrow(snuc), paste(range(snuc$beta), collapse = " to "), nrow(disease),
    paste(range(disease$fdr), collapse = " to "),
    paste(disease[interval_available == TRUE, cell_label], collapse = ", "),
    nrow(structure[residue >= 114 & residue <= 232]),
    events[event == "project_exon5_exon6", protein_transition],
    round(mean(structure[residue %in% c(150, 151), pLDDT]), 2),
    file.info(exports[1])$size, file.info(exports[2])$size, file.info(exports[3])$size, file.info(exports[4])$size,
    "A-E", "portable SVG"
  )
)

fwrite(checks, file.path(out_dir, "Figure_5_QA.tsv"), sep = "\t")
if (any(!checks$passed)) {
  print(checks[passed == FALSE])
  stop("Figure 5 QA failed", call. = FALSE)
}
message("Figure 5 QA passed: ", nrow(checks), " checks")
