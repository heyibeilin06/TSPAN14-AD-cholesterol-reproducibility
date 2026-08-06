#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", winslash = "/", mustWork = TRUE)
source <- file.path(root, "audit", "reviewer_revision")
dest <- file.path(root, "figures", "Figure_03", "data")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

donor <- fread(file.path(source, "count_level_acceptor_choice_donors.tsv"))
models <- fread(file.path(source, "count_level_acceptor_choice_models.tsv"))
depth <- fread(file.path(source, "count_level_acceptor_choice_depth_sensitivity.tsv"))

stopifnot(nrow(donor) == 147L, all(c(0L, 1L, 2L) %in% donor$genotype),
          nrow(models[grepl("Firth|Beta-binomial", analysis)]) == 3L, nrow(depth) == 5L)

fwrite(donor, file.path(dest, "Figure_3_ba24_canonical_cryptic_counts.tsv"), sep = "\t")
fwrite(models, file.path(dest, "Figure_3_canonical_cryptic_models.tsv"), sep = "\t")
fwrite(depth, file.path(dest, "Figure_3_detection_depth_sensitivity.tsv"), sep = "\t")
message("Prepared revised Figure 3B source tables")
