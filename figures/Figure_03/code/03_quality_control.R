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
out_dir <- normalizePath(value_after("--output-dir", file.path(root, "figures", "Figure_03", "output")), winslash = "/", mustWork = TRUE)
source_dir <- normalizePath(value_after("--source-dir", file.path(root, "figures", "Figure_03", "data")), winslash = "/", mustWork = TRUE)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, passed = isTRUE(passed), detail = as.character(detail))
}

replication <- fread(file.path(source_dir, "Figure_3_exact_replication_matrix.tsv"))
donor <- fread(file.path(source_dir, "Figure_3_ba24_donor_ier.tsv"))
depth <- fread(file.path(source_dir, "Figure_3_depth_sensitivity_delta_ier.tsv"))
counts <- fread(file.path(source_dir, "Figure_3_brain_junction_counts.tsv"))
cousage <- fread(file.path(source_dir, "Figure_3_brain_cousage_summary.tsv"))
events <- fread(file.path(source_dir, "Figure_3_splice_events.tsv"))

add_check("six variants by four tissues", nrow(replication) == 24 && uniqueN(replication$snp) == 6 && uniqueN(replication$tissue_label) == 4, nrow(replication))
add_check("all exact-event NES positive", all(replication$nes > 0), range(replication$nes))
add_check("147 BA24 donors", nrow(donor) == 147 && uniqueN(donor$donor_id) == 147, nrow(donor))
add_check("three genotype groups", setequal(unique(donor$genotype), 0:2), paste(sort(unique(donor$genotype)), collapse = ","))
add_check("raw IER bounded", all(donor$raw_ier >= 0 & donor$raw_ier <= 1), range(donor$raw_ier))
add_check("primary DeltaIER exceeds one point", depth[minimum_cluster_reads == 1, per_alt_allele_delta_ier_percentage_points] > 1, depth[minimum_cluster_reads == 1, per_alt_allele_delta_ier_percentage_points])
add_check("primary DeltaIER CI excludes one point", depth[minimum_cluster_reads == 1, bootstrap_95ci_low] > 1, depth[minimum_cluster_reads == 1, bootstrap_95ci_low])
add_check("five depth thresholds", nrow(depth) == 5, nrow(depth))
add_check("2642 junction samples", nrow(counts) == 2642, nrow(counts))
add_check("junction co-usage strong", cor(counts$log_ex5_6, counts$log_ex6_7, method = "spearman") > 0.92, cor(counts$log_ex5_6, counts$log_ex6_7, method = "spearman"))
add_check("13 brain-region summaries", nrow(cousage) == 13, nrow(cousage))
add_check("three splice event classes", setequal(unique(events$event_class), c("exact", "competing", "adjacent")), paste(events$event_class, collapse = ","))

stem <- file.path(out_dir, "Figure_3_exact_splicing_consistency_v10")
for (extension in c("png", "pdf", "svg", "tiff")) {
  path <- paste0(stem, ".", extension)
  add_check(paste(extension, "export exists"), file.exists(path) && file.info(path)$size > 10000, if (file.exists(path)) file.info(path)$size else 0)
}

svg_path <- paste0(stem, ".svg")
if (file.exists(svg_path)) {
  svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  add_check("uppercase panel labels", all(vapply(c("A", "B", "C", "D"), function(x) grepl(paste0(">", x, "<"), svg_text, fixed = TRUE), logical(1))), "A-D")
  add_check("no local paths in SVG", !grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", svg_text, perl = TRUE), "portable SVG")
}

qa <- rbindlist(checks)
fwrite(qa, file.path(out_dir, "Figure_3_QA.tsv"), sep = "\t")
if (any(!qa$passed)) {
  print(qa[passed == FALSE])
  stop("Figure 3 QA failed", call. = FALSE)
}
message("Figure 3 QA passed: ", nrow(qa), " checks")
